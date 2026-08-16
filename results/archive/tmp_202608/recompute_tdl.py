#!/usr/bin/env python3
"""Recompute TDL (estimated from Open Targets) using approved-drug phase first,
then rewrite the CSV in place. Does not re-query the API."""
import csv
import json
import re

RAW = "/data/qiushuogeng/projects/dual-channel-mr-atlas/tmp/drug_annotate_raw.json"
CSV = "/data/qiushuogeng/projects/dual-channel-mr-atlas/results/drug_annotation_20260813.csv"


def map_phase(stage):
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
    return None


def estimate_tdl(biotype, tractability, n_drugs, max_phase):
    if biotype != "protein_coding":
        return ""
    if n_drugs and max_phase is not None:
        if max_phase >= 4:
            return "Tclin(est)"
        if max_phase >= 1:
            return "Tchem(est)"
    flags = {}
    for t in tractability or []:
        flags.setdefault(t.get("label", ""), t.get("value"))
    if flags.get("Approved Drug") or flags.get("Advanced Clinical"):
        return "Tclin(est)"
    if (flags.get("Phase 1 Clinical") or flags.get("High-Quality Ligand")
            or flags.get("Structure with Ligand")):
        return "Tchem(est)"
    bio = any(v and re.search(r"bio|antibody|Protac", lbl, re.I) for lbl, v in flags.items())
    if bio:
        return "Tbio(est)"
    return "Tdark(est)"


def main():
    raw = json.load(open(RAW))
    rows = list(csv.DictReader(open(CSV)))
    for r in rows:
        resp = raw.get(r["ensg"])
        if not resp or "__error__" in resp:
            continue
        t = resp.get("data", {}).get("target") or {}
        stages = [(dr.get("drug") or {}).get("name") and dr.get("maxClinicalStage")
                  for dr in (t.get("drugAndClinicalCandidates") or {}).get("rows") or []]
        phases = [map_phase(s) for s in stages if s]
        phases = [p for p in phases if p is not None]
        n_drugs = int(r["n_known_drugs"] or 0)
        max_phase = max(phases) if phases else (None if n_drugs else 0)
        r["TDL"] = estimate_tdl(t.get("biotype", ""), t.get("tractability"), n_drugs, max_phase)
    with open(CSV, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print("TDL recomputed. Done.")
    for r in rows:
        if int(r["n_known_drugs"] or 0) > 0:
            print(r["gene_symbol"], r["TDL"], r["max_clinical_phase"], r["drug_names"][:60])


if __name__ == "__main__":
    main()
