#!/usr/bin/env python3
"""
Investigate why sequences were dropped in the protein (CDS) alignment step.
Run from your ebolaseq output directory (e.g. myanalyses) or pass its path.

Usage:
  cd myanalyses && python scripts/investigate_dropped.py
  python scripts/investigate_dropped.py /path/to/myanalyses
  python scripts/investigate_dropped.py --ids id1 id2 ...   # from current dir

Prints sequence length for each dropped ID. Complete Ebola genomes are ~18 800–19 000 nt;
shorter sequences are likely partial and missing the 5' region (NP/VP35).
"""

from __future__ import print_function
import argparse
import os
import sys

# Complete Ebola genome ref lengths (nt)
REF_LENGTHS = {"Sudan": 18875, "Zaire": 18959, "Bundibugyo": 18940, "Tai Forest": 18935, "Reston": 18891}
COMPLETE_MIN = 18000  # below this consider partial


def main():
    ap = argparse.ArgumentParser(description="Investigate dropped sequences: show length from combined FASTA.")
    ap.add_argument("run_dir", nargs="?", default=os.getcwd(), help="Ebolaseq output directory (default: cwd)")
    ap.add_argument("--ids", nargs="*", help="Dropped IDs to look up (default: parse ebolaseq_run.log)")
    args = ap.parse_args()

    run_dir = os.path.abspath(args.run_dir)
    combined_fasta = os.path.join(run_dir, "Alignment", "FASTA", "Ebola_Combined.fasta")
    log_path = os.path.join(run_dir, "ebolaseq_run.log")

    if not os.path.isfile(combined_fasta):
        print("Not found: %s" % combined_fasta, file=sys.stderr)
        print("Run this script from your ebolaseq output dir (e.g. myanalyses) or pass it as first argument.", file=sys.stderr)
        sys.exit(1)

    id_to_len = {}
    try:
        from Bio import SeqIO
        for rec in SeqIO.parse(combined_fasta, "fasta"):
            id_to_len[rec.id] = len(rec.seq)
    except Exception as e:
        print("Error reading FASTA: %s" % e, file=sys.stderr)
        sys.exit(1)

    if args.ids:
        dropped_ids = args.ids
    else:
        dropped_ids = []
        if os.path.isfile(log_path):
            in_dropped = False
            with open(log_path) as f:
                for line in f:
                    line = line.rstrip()
                    if "Dropped (protein coverage not met)" in line:
                        in_dropped = True
                        continue
                    if in_dropped and line.startswith("    "):
                        sid = line.strip()
                        if sid and " dropped" not in sid:
                            dropped_ids.append(sid)
                    elif in_dropped and line.startswith("  ") and not line.startswith("    "):
                        continue
                    else:
                        in_dropped = False
        if not dropped_ids:
            print("No dropped IDs in ebolaseq_run.log. Pass IDs with --ids.", file=sys.stderr)
            sys.exit(1)

    # Deduplicate preserving order
    seen = set()
    unique = []
    for i in dropped_ids:
        if i not in seen:
            seen.add(i)
            unique.append(i)

    print("Dropped sequence lengths (complete Ebola genome ~18 800–19 000 nt)")
    print("Short sequences likely lack 5' region → NP/VP35 coverage < 50%.\n")
    for did in unique:
        L = id_to_len.get(did)
        if L is None:
            print("  %s  L=?" % did)
        else:
            note = " (partial)" if L < COMPLETE_MIN else ""
            print("  %s  L=%d nt%s" % (did, L, note))
    print("\nRef lengths: %s" % ", ".join("%s %d nt" % (k, v) for k, v in sorted(REF_LENGTHS.items())))


if __name__ == "__main__":
    main()
