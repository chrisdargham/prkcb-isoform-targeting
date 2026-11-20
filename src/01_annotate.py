#!/usr/bin/env python3
"""
01_annotate.py
--------------
Domain annotation of PKC-beta from validated sources.

This replaces v1's homemade regular expressions. The reason is not that regex is
crude but that it is circular: v1's patterns were written by reading PKC-beta's
sequence and were then used to search PKC-beta. A pattern derived from a
sequence necessarily matches that sequence, so the exercise confirmed nothing,
produced no measure of statistical significance, and reported the extent of the
author's pattern rather than the extent of the domain.

Sources used here, in order of authority:

  1. UniProt P05771 (reviewed)  - expert-curated feature table.
  2. InterProScan                - 13 member databases run over the sequence,
                                   each an independent prediction with its own
                                   scoring model. Profile HMMs (Pfam, SMART)
                                   carry E-values: the number of hits this good
                                   expected by chance in a database this size.

Agreement across independent databases is the evidence. Where Pfam, PROSITE,
SMART and the UniProt curators converge on a boundary, that boundary is real in
a way no single method establishes.

Inputs : data/raw/interproscan_PKCB1.tsv
         data/raw/uniprot_P05771.txt
         data/raw/PKCB1_human_NP_997700.fasta
Outputs: out/annotation.json
         out/annotation_report.txt
         out/v1_regex_comparison.tsv
"""

import json
import re
import sys
from collections import defaultdict
from pathlib import Path

# python -I does not put the script's own directory on sys.path, so a sibling
# module has to be found explicitly. Done here rather than by dropping -I,
# because -I also keeps stray user site-packages out of the run.
sys.path.insert(0, str(Path(__file__).resolve().parent))
import provenance  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
OUT = ROOT / "out"
OUT.mkdir(exist_ok=True)

IPR_TSV = RAW / "interproscan_PKCB1.tsv"
UNIPROT = RAW / "uniprot_P05771.txt"
FASTA = RAW / "PKCB1_human_NP_997700.fasta"

# Databases whose hits are profile-based and carry interpretable E-values.
# Ordered by how much weight we give them when forming a consensus.
PROFILE_DBS = ["PROSITE profiles", "SMART", "Pfam", "CDD", "SUPERFAMILY"]

# ---------------------------------------------------------------------------
# The coordinates v1 reported, reproduced here verbatim so the comparison is
# reproducible rather than asserted. Source: v1 out/protein_analysis.json.
# 'None' means v1 had no pattern for that feature and therefore missed it.
# ---------------------------------------------------------------------------
V1_REGEX = {
    "C1A zinc finger": (37, 71),
    "C1B zinc finger": (115, 135),
    "C2 domain": None,
    "Protein kinase domain": (348, 511),   # v1 called this "catalytic core"
    "AGC-kinase C-terminal": None,
    "Active site (proton acceptor)": (466, 466),
    "PDK1 phosphosite": (500, 500),
}

# Which curated feature each v1 entry should be compared against.
V1_TO_CURATED = {
    "C1A zinc finger": ("ZN_FING", "Phorbol-ester/DAG-type 1"),
    "C1B zinc finger": ("ZN_FING", "Phorbol-ester/DAG-type 2"),
    "C2 domain": ("DOMAIN", "C2"),
    "Protein kinase domain": ("DOMAIN", "Protein kinase"),
    "AGC-kinase C-terminal": ("DOMAIN", "AGC-kinase C-terminal"),
    "Active site (proton acceptor)": ("ACT_SITE", "Proton acceptor"),
    "PDK1 phosphosite": ("MOD_RES", "Phosphothreonine; by PDPK1"),
}


# ---------------------------------------------------------------------------
# parsers
# ---------------------------------------------------------------------------
def parse_interproscan(path):
    """InterProScan TSV -> list of hit dicts.

    The format is 11-15 tab-separated columns; trailing columns (InterPro
    accession, description, GO, pathways) are present only for signatures
    integrated into InterPro. Rows without an InterPro accession are kept but
    flagged, since they are predictions no integrated entry vouches for.
    """
    hits = []
    for line in path.read_text().splitlines():
        if not line.strip():
            continue
        f = line.split("\t")
        if len(f) < 9:
            continue
        try:
            start, stop = int(f[6]), int(f[7])
        except ValueError:
            continue
        score = f[8].strip()
        # Column 9 is an E-value for profile methods. Some rows carry a
        # different statistic close to 1.0; those are not E-values and must not
        # be read as such.
        try:
            evalue = float(score)
        except ValueError:
            evalue = None
        hits.append(dict(
            db=f[3].strip(),
            signature=f[4].strip(),
            description=f[5].strip(),
            start=start,
            stop=stop,
            evalue=evalue,
            interpro=f[11].strip() if len(f) > 11 and f[11].strip() != "-" else None,
            interpro_desc=f[12].strip() if len(f) > 12 and f[12].strip() != "-" else None,
        ))
    return hits


def parse_uniprot_features(path):
    """UniProt flat file -> list of curated features.

    FT lines come as a key line carrying the position, optionally followed by
    continuation lines holding /note="...". Isoform-scoped positions such as
    'P05771-2:641' are recorded but marked, since they index a different
    sequence and must not be mixed with canonical coordinates.
    """
    feats = []
    current = None
    for line in path.read_text().splitlines():
        if not line.startswith("FT"):
            continue
        body = line[2:].strip()
        m = re.match(r"^([A-Z_]+)\s+(\S+)", body)
        if m:
            if current:
                feats.append(current)
            key, pos = m.group(1), m.group(2)
            isoform = None
            if ":" in pos:
                isoform, pos = pos.split(":", 1)
            rng = pos.split("..")
            try:
                start = int(rng[0].lstrip("<>?"))
                stop = int(rng[-1].lstrip("<>?"))
            except ValueError:
                current = None
                continue
            current = dict(key=key, start=start, stop=stop,
                           note="", isoform=isoform)
        elif current is not None:
            n = re.search(r'/note="([^"]*)"?', body)
            if n:
                current["note"] = (current["note"] + " " + n.group(1)).strip()
    if current:
        feats.append(current)
    return feats


def load_sequence(path):
    lines = [l.strip() for l in path.read_text().splitlines() if l.strip()]
    return "".join(l for l in lines if not l.startswith(">"))


# ---------------------------------------------------------------------------
# consensus
# ---------------------------------------------------------------------------
def build_consensus(ipr_hits, uniprot_feats):
    """Group independent predictions by InterPro entry and report the spread.

    The point of this function is that it does NOT pick a winner. It reports
    every database's boundary call for each domain, so disagreement between
    methods is visible rather than hidden behind a single number.
    """
    groups = defaultdict(list)
    for h in ipr_hits:
        if h["interpro"] and h["db"] in PROFILE_DBS:
            groups[(h["interpro"], h["interpro_desc"])].append(h)

    consensus = []
    for (ipr, desc), hits in sorted(groups.items()):
        # Several C1 fingers share one InterPro entry, so cluster by position
        # before averaging; otherwise two domains collapse into one.
        hits = sorted(hits, key=lambda h: h["start"])
        clusters = []
        for h in hits:
            placed = False
            for c in clusters:
                if abs(h["start"] - c[0]["start"]) < 40:
                    c.append(h)
                    placed = True
                    break
            if not placed:
                clusters.append([h])

        for c in clusters:
            starts = [h["start"] for h in c]
            stops = [h["stop"] for h in c]
            evalues = [h["evalue"] for h in c
                       if h["evalue"] is not None and h["evalue"] < 0.01]
            consensus.append(dict(
                interpro=ipr,
                description=desc,
                n_databases=len({h["db"] for h in c}),
                databases=sorted({h["db"] for h in c}),
                start_min=min(starts), start_max=max(starts),
                stop_min=min(stops), stop_max=max(stops),
                boundary_spread=max(starts) - min(starts) + max(stops) - min(stops),
                best_evalue=min(evalues) if evalues else None,
                calls=[dict(db=h["db"], signature=h["signature"],
                            start=h["start"], stop=h["stop"],
                            evalue=h["evalue"]) for h in c],
            ))

    # attach the curated UniProt call where one matches
    for c in consensus:
        match = None
        for f in uniprot_feats:
            if f["isoform"]:
                continue
            if f["key"] not in ("DOMAIN", "ZN_FING"):
                continue
            if abs(f["start"] - c["start_min"]) < 40:
                match = f
                break
        c["uniprot"] = (dict(key=match["key"], start=match["start"],
                             stop=match["stop"], note=match["note"])
                        if match else None)
    return sorted(consensus, key=lambda c: c["start_min"])


def compare_v1(uniprot_feats):
    """v1's regex coordinates against the curated truth."""
    rows = []
    for name, v1 in V1_REGEX.items():
        key, note_frag = V1_TO_CURATED[name]
        curated = None
        for f in uniprot_feats:
            if f["isoform"]:
                continue
            if f["key"] == key and note_frag.lower() in f["note"].lower():
                curated = f
                break
        row = dict(feature=name,
                   v1_start=v1[0] if v1 else None,
                   v1_end=v1[1] if v1 else None,
                   curated_start=curated["start"] if curated else None,
                   curated_end=curated["stop"] if curated else None)
        if v1 is None:
            row["verdict"] = "MISSED ENTIRELY by v1"
            row["start_error"] = row["end_error"] = None
            row["residues_missed"] = (curated["stop"] - curated["start"] + 1
                                      if curated else None)
        elif curated is None:
            row["verdict"] = "no curated counterpart found"
            row["start_error"] = row["end_error"] = None
            row["residues_missed"] = None
        else:
            se = v1[0] - curated["start"]
            ee = v1[1] - curated["stop"]
            row["start_error"] = se
            row["end_error"] = ee
            covered = min(v1[1], curated["stop"]) - max(v1[0], curated["start"]) + 1
            total = curated["stop"] - curated["start"] + 1
            row["residues_missed"] = total - max(0, covered)
            row["verdict"] = ("exact" if se == 0 and ee == 0
                              else f"off by {se:+d} / {ee:+d}")
        rows.append(row)
    return rows


# ---------------------------------------------------------------------------
def main():
    seq = load_sequence(FASTA)
    ipr = parse_interproscan(IPR_TSV)
    up = parse_uniprot_features(UNIPROT)
    consensus = build_consensus(ipr, up)
    comparison = compare_v1(up)

    prov = provenance.stamp(__file__, [IPR_TSV, UNIPROT, FASTA], root=ROOT)

    payload = dict(
        provenance=prov,
        sequence_length=len(seq),
        sources=dict(
            interproscan=dict(file=IPR_TSV.name, n_hits=len(ipr),
                              databases=sorted({h["db"] for h in ipr})),
            uniprot=dict(file=UNIPROT.name, accession="P05771",
                         n_features=len(up)),
        ),
        consensus_domains=consensus,
        curated_features=[f for f in up if not f["isoform"]],
        v1_regex_comparison=comparison,
    )
    (OUT / "annotation.json").write_text(json.dumps(payload, indent=2))

    # ---- readable report ---------------------------------------------------
    L = []
    L.append("PKC-beta domain annotation from validated sources")
    L.append("=" * 76)
    L.append(provenance.summary_line(prov))
    L.append("")
    L.append(f"Sequence: {len(seq)} aa")
    L.append(f"InterProScan: {len(ipr)} hits across "
             f"{len({h['db'] for h in ipr})} member databases")
    L.append(f"UniProt P05771: {len([f for f in up if not f['isoform']])} "
             f"canonical features")
    L.append("")
    L.append("CONSENSUS DOMAINS (independent databases agreeing)")
    L.append("")
    for c in consensus:
        ev = f"{c['best_evalue']:.1e}" if c["best_evalue"] else "n/a"
        L.append(f"  {c['description'] or c['interpro']}")
        L.append(f"    InterPro {c['interpro']}   best E-value {ev}")
        L.append(f"    {c['n_databases']} databases: {', '.join(c['databases'])}")
        if c["start_min"] == c["start_max"] and c["stop_min"] == c["stop_max"]:
            L.append(f"    boundaries: {c['start_min']}-{c['stop_min']} "
                     f"(unanimous)")
        else:
            L.append(f"    boundaries: {c['start_min']}-{c['start_max']} to "
                     f"{c['stop_min']}-{c['stop_max']} "
                     f"(spread {c['boundary_spread']} aa)")
        if c["uniprot"]:
            L.append(f"    UniProt curated: {c['uniprot']['start']}-"
                     f"{c['uniprot']['stop']}  \"{c['uniprot']['note']}\"")
        L.append("")

    L.append("")
    L.append("WHAT v1's REGEX GOT WRONG")
    L.append("")
    L.append(f"  {'feature':<32}{'v1':<14}{'curated':<14}{'verdict'}")
    for r in comparison:
        v1s = f"{r['v1_start']}-{r['v1_end']}" if r["v1_start"] else "not found"
        cur = (f"{r['curated_start']}-{r['curated_end']}"
               if r["curated_start"] else "?")
        L.append(f"  {r['feature']:<32}{v1s:<14}{cur:<14}{r['verdict']}")
    L.append("")
    missed = sum(r["residues_missed"] or 0 for r in comparison)
    exact = sum(1 for r in comparison if r["verdict"] == "exact")
    L.append(f"  Exactly correct: {exact} of {len(comparison)} features")
    L.append(f"  Residues inside real domains that v1 failed to cover: {missed}")

    (OUT / "annotation_report.txt").write_text("\n".join(L))

    with (OUT / "v1_regex_comparison.tsv").open("w") as fh:
        cols = ["feature", "v1_start", "v1_end", "curated_start", "curated_end",
                "start_error", "end_error", "residues_missed", "verdict"]
        fh.write("\t".join(cols) + "\n")
        for r in comparison:
            fh.write("\t".join("" if r.get(c) is None else str(r.get(c))
                               for c in cols) + "\n")

    print("\n".join(L))


if __name__ == "__main__":
    main()
