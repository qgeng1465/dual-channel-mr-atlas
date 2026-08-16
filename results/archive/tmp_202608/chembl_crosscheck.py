#!/usr/bin/env python3
"""
ChEMBL cross-check for protein-coding genes that had 0 drugs from Open Targets.
Only clinical-stage (max_phase>=1) mechanisms-of-action are kept.
Pure API data; failures recorded per gene.
"""
import json
import time
import urllib.request
import urllib.parse

CHEMBL = "https://www.ebi.ac.uk/chembl/api/data"
RAW = "/data/qiushuogeng/projects/dual-channel-mr-atlas/tmp/drug_annotate_raw.json"
OUT = "/data/qiushuogeng/projects/dual-channel-mr-atlas/tmp/chembl_crosscheck.json"


def get(url):
    req = urllib.request.Request(url, headers={"User-Agent": "research-annotation/1.0"})
    with urllib.request.urlopen(req, timeout=40) as resp:
        return json.loads(resp.read().decode())


def main():
    raw = json.load(open(RAW))
    targets = []
    for ensg, resp in raw.items():
        if "__error__" in resp:
            continue
        t = resp.get("data", {}).get("target") or {}
        if t.get("biotype") != "protein_coding":
            continue
        if (t.get("drugAndClinicalCandidates") or {}).get("rows"):
            continue  # already has OT drugs
        targets.append((ensg, t.get("approvedSymbol") or ensg))

    print(f"Cross-checking {len(targets)} protein-coding genes with 0 OT drugs")
    out = {}
    for i, (ensg, sym) in enumerate(targets, 1):
        rec = {"ensg": ensg, "symbol": sym, "status": "ok", "mechanisms": []}
        try:
            # 1. find ChEMBL target by gene symbol
            q = urllib.parse.urlencode({"target_synonym": sym, "limit": 50})
            tr = get(f"{CHEMBL}/target.json?{q}")
            hs = [x for x in tr.get("targets", []) if x.get("organism") == "Homo sapiens"]
            # prefer single-protein target
            hs.sort(key=lambda x: 0 if x.get("target_type") == "SINGLE PROTEIN" else 1)
            if not hs:
                rec["status"] = "no_hs_target"
                out[ensg] = rec
                print(f"[{i}/{len(targets)}] {sym}: no human ChEMBL target")
                continue
            tid = hs[0]["target_chembl_id"]
            rec["chembl_target_id"] = tid
            # 2. mechanisms
            q = urllib.parse.urlencode({"target_chembl_id": tid, "limit": 200})
            mech = get(f"{CHEMBL}/mechanism.json?{q}")
            for m in mech.get("mechanisms", []):
                mp = m.get("max_phase")
                try:
                    mp = float(mp) if mp is not None else 0.0
                except (TypeError, ValueError):
                    mp = 0.0
                if mp >= 1:
                    mid = m.get("molecule_chembl_id")
                    # 3. drug name
                    name = None
                    try:
                        q2 = urllib.parse.urlencode({"molecule_chembl_id": mid})
                        mol = get(f"{CHEMBL}/molecule.json?{q2}")
                        if mol.get("molecules"):
                            name = mol["molecules"][0].get("pref_name")
                    except Exception:
                        name = None
                    rec["mechanisms"].append({
                        "molecule": mid,
                        "name": name,
                        "max_phase": int(mp),
                        "moea": m.get("mechanism_of_action"),
                    })
                    time.sleep(0.1)
            out[ensg] = rec
            flag = f" CLINICAL x{len(rec['mechanisms'])}" if rec["mechanisms"] else ""
            print(f"[{i}/{len(targets)}] {sym}: {len(rec['mechanisms'])} clinical mech{flag}")
        except Exception as e:  # noqa: BLE001
            rec["status"] = f"error: {e}"
            out[ensg] = rec
            print(f"[{i}/{len(targets)}] {sym}: ERROR {e}")
        time.sleep(0.15)

    json.dump(out, open(OUT, "w"), indent=1)
    print(f"\nSaved -> {OUT}")
    n_clin = sum(1 for v in out.values() if v.get("mechanisms"))
    print(f"Genes with clinical ChEMBL mechanisms: {n_clin}")


if __name__ == "__main__":
    main()
