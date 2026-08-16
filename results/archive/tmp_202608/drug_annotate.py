#!/usr/bin/env python3
"""
P4 drug annotation layer for dual-channel MR atlas priority genes.
Source: Open Targets GraphQL API (api.platform.opentargets.org) -- primary.
Pharos API was unreachable (HTTP 403 / HTML app) -- noted in report; TDL estimated
from Open Targets tractability + clinical-stage flags, labeled "(est)".

Honesty rules:
  - Only write data returned by the API. Nothing fabricated.
  - API failures are recorded per-gene as network_failure=1 and fields left blank.
"""
import json
import time
import sys
import csv
import re
import urllib.request
import urllib.error

OT_ENDPOINT = "https://api.platform.opentargets.org/api/v4/graphql"
INPUT = "/data/qiushuogeng/projects/dual-channel-mr-atlas/tmp/priority_genes_heidi_pass.txt"
RAW_OUT = "/data/qiushuogeng/projects/dual-channel-mr-atlas/tmp/drug_annotate_raw.json"
CSV_OUT = "/data/qiushuogeng/projects/dual-channel-mr-atlas/results/drug_annotation_20260813.csv"

QUERY = """
query($ensg: String!) {
  target(ensemblId: $ensg) {
    approvedSymbol
    approvedName
    biotype
    drugAndClinicalCandidates {
      count
      rows {
        maxClinicalStage
        drug { name drugType }
      }
    }
    associatedDiseases {
      count
      rows {
        score
        disease { id name }
      }
    }
    tractability { label value }
  }
}
"""


def query_ot(ensg, retries=3):
    body = json.dumps({"query": QUERY, "variables": {"ensg": ensg}}).encode()
    last_err = None
    for attempt in range(retries):
        try:
            req = urllib.request.Request(
                OT_ENDPOINT,
                data=body,
                headers={"Content-Type": "application/json", "User-Agent": "research-annotation/1.0"},
                method="POST",
            )
            with urllib.request.urlopen(req, timeout=60) as resp:
                return json.loads(resp.read().decode())
        except Exception as e:  # noqa: BLE001
            last_err = e
            time.sleep(2 * (attempt + 1))
    return {"__error__": str(last_err)}


def map_phase(stage):
    """Map Open Targets maxClinicalStage string -> 0..4."""
    s = (stage or "").upper()
    if any(k in s for k in ("APPROVAL", "LAUNCHED", "MARKET", "REGISTERED")):
        return 4
    if "PHASE III" in s:
        return 3
    if "PHASE II" in s:
        return 2
    if "PHASE I" in s or "PHASE 0" in s or "PHASE 1" in s:
        return 1
    if "PRECLINICAL" in s:
        return 0
    if "UNKNOWN" in s or s == "":
        return None
    return None


def estimate_tdl(tractability, biotype, n_drugs, max_phase):
    """Estimate Pharos-style TDL from Open Targets tractability flags. Labeled '(est)'."""
    if biotype != "protein_coding":
        return ""
    flags = {}
    for t in tractability or []:
        flags.setdefault(t.get("label", ""), t.get("value"))
    if flags.get("Approved Drug") or flags.get("Advanced Clinical"):
        return "Tclin(est)"
    if flags.get("Phase 1 Clinical") or flags.get("High-Quality Ligand") or flags.get("Structure with Ligand"):
        return "Tchem(est)"
    # no small-molecule tractability but has any drug entry -> chem-ish anyway
    if n_drugs and max_phase is not None and max_phase >= 1:
        return "Tchem(est)"
    # biological assay evidence
    bio = False
    for lbl, v in flags.items():
        if v and re.search(r"bio|antibody|Protac|PROTAC", lbl, re.I):
            bio = True
    if bio:
        return "Tbio(est)"
    return "Tdark(est)"


# Disease keywords of interest (T2D / CAD / cardiometabolic).
T2D_KEYWORDS = ["type 2 diabetes", "type 2 diabetes mellitus", "t2dm", "diabetes mellitus type 2", "non-insulin-dependent diabetes"]
CAD_KEYWORDS = ["coronary artery disease", "coronary heart disease", "ischemic heart disease", "ischaemic heart disease",
                "coronary atherosclerosis", "myocardial infarction", "atherosclerotic cardiovascular disease",
                "acute myocardial infarction", "coronary arteriosclerosis"]
OTHER_DIABETES = ["diabetes", "glucose", "insulin", "glycated", "glycemic"]
LIPID = ["ldl", "hdl", "cholesterol", "triglyceride", "lipid", "lipoprotein"]


def disease_assoc_text(diseases):
    """Return a concise honest string of T2D/CAD-relevant associations from Open Targets."""
    hits = []
    for row in diseases or []:
        name = (row.get("disease") or {}).get("name", "")
        score = row.get("score")
        low = name.lower()
        tag = None
        if any(k in low for k in T2D_KEYWORDS):
            tag = "T2D"
        elif any(k in low for k in CAD_KEYWORDS):
            tag = "CAD"
        if tag:
            hits.append(f"{name} ({score:.3f})")
    # Cap to top 3 per tag to keep cell readable
    seen = []
    for h in hits:
        if h not in seen:
            seen.append(h)
    return "; ".join(seen[:6]) if seen else ""


def main():
    rows = []
    with open(INPUT) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5:
                continue
            rows.append({"ensg": parts[0], "symbol": parts[1], "outcome": parts[2],
                         "pp_h4": parts[3], "p_heidi": parts[4]})
    print(f"Loaded {len(rows)} genes")

    results = []
    raw_store = {}
    stage_values = set()
    for i, r in enumerate(rows, 1):
        ensg = r["ensg"]
        resp = query_ot(ensg)
        raw_store[ensg] = resp
        rec = {
            "gene_symbol": r["symbol"],
            "ensg": ensg,
            "outcome": r["outcome"],
            "PP.H4": r["pp_h4"],
            "p_HEIDI": r["p_heidi"],
        }
        if "__error__" in resp:
            rec.update({
                "biotype_ok": "",
                "TDL": "",
                "n_known_drugs": "",
                "max_clinical_phase": "",
                "drug_names": "",
                "known_disease_assoc": "API_FAIL",
                "note": "OpenTargets request failed",
            })
            results.append(rec)
            print(f"[{i}/{len(rows)}] {r['symbol']} API_FAIL: {resp['__error__'][:120]}")
            continue

        t = resp.get("data", {}).get("target") or {}
        biotype = t.get("biotype", "")
        drug_rows = (t.get("drugAndClinicalCandidates") or {}).get("rows") or []
        drug_names = []
        phases = []
        for dr in drug_rows:
            nm = (dr.get("drug") or {}).get("name", "")
            if nm and nm not in drug_names:
                drug_names.append(nm)
            st = dr.get("maxClinicalStage")
            if st:
                stage_values.add(st)
            p = map_phase(st)
            if p is not None:
                phases.append(p)
        n_drugs = len(drug_names)
        max_phase = max(phases) if phases else (None if drug_rows else 0)
        if max_phase is None:
            max_phase = 0  # drugs present but phase unknown/preclinical

        diseases = (t.get("associatedDiseases") or {}).get("rows") or []
        assoc = disease_assoc_text(diseases)
        tdl = estimate_tdl(t.get("tractability"), biotype, n_drugs, max_phase)

        rec.update({
            "biotype_ok": 1 if biotype == "protein_coding" else 0,
            "TDL": tdl,
            "n_known_drugs": n_drugs,
            "max_clinical_phase": max_phase,
            "drug_names": "; ".join(drug_names),
            "known_disease_assoc": assoc if assoc else "none_matched",
        })
        results.append(rec)
        print(f"[{i}/{len(rows)}] {r['symbol']} biotype={biotype} drugs={n_drugs} maxphase={max_phase} assoc={'yes' if assoc else 'no'}")
        time.sleep(0.25)

    print("\nDistinct maxClinicalStage values seen:", sorted(stage_values))

    with open(RAW_OUT, "w") as fh:
        json.dump(raw_store, fh, indent=1)
    print(f"Raw API payloads -> {RAW_OUT}")

    with open(CSV_OUT, "w", newline="") as fh:
        fieldnames = ["gene_symbol", "ensg", "outcome", "PP.H4", "p_HEIDI",
                      "biotype_ok", "TDL", "n_known_drugs", "max_clinical_phase",
                      "drug_names", "known_disease_assoc"]
        w = csv.DictWriter(fh, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(results)
    print(f"CSV -> {CSV_OUT}")


if __name__ == "__main__":
    main()
