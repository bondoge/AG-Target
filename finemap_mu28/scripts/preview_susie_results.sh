#!/usr/bin/env bash
set -euo pipefail

OUT_DIR="${1:-${OUT:-/home/ec2-user/susie/data/finemap_mu28}}"
NFILES="${2:-5}"
NLINES="${3:-10}"

RESULTS="$OUT_DIR/results"

echo "Results directory: $RESULTS"
echo

echo "Result file counts:"
echo -n "  PIP files:      "
find "$RESULTS" -maxdepth 1 -name "*.susie_pip.tsv" | wc -l
echo -n "  Summary files:  "
find "$RESULTS" -maxdepth 1 -name "*.susie_summary.txt" | wc -l
echo -n "  RDS fit files:  "
find "$RESULTS" -maxdepth 1 -name "*.susie_fit.rds" | wc -l
echo

echo "All result files:"
ls -lh "$RESULTS"
echo

echo "Preview first $NFILES PIP files, first $NLINES lines each:"
echo "============================================================"

find "$RESULTS" -maxdepth 1 -name "*.susie_pip.tsv" | sort | head -n "$NFILES" | while read -r f
do
  echo
  echo "FILE: $f"
  echo "SIZE: $(du -h "$f" | cut -f1)"
  echo "LINES: $(wc -l < "$f")"
  echo "---- head ----"
  head -n "$NLINES" "$f"
done

echo
echo "Preview first $NFILES summary files:"
echo "============================================================"

find "$RESULTS" -maxdepth 1 -name "*.susie_summary.txt" | sort | head -n "$NFILES" | while read -r f
do
  echo
  echo "FILE: $f"
  echo "---- content ----"
  cat "$f"
done
