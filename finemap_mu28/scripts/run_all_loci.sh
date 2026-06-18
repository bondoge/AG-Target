#!/usr/bin/env bash
set -euo pipefail

tail -n +2 "$OUT/loci.tsv" | while IFS=$'\t' read -r locus_id lead_snpid chr lead_pos start end lead_ea lead_nea lead_p N
do
  echo "Processing $locus_id chr=$chr start=$start end=$end"

  LDIR="$OUT/loci/$locus_id"
  mkdir -p "$LDIR"

  if [ -f "$OUT/results/$locus_id.susie_pip.tsv" ]; then
    echo "Skipping $locus_id: SuSiE result already exists"
    continue
  fi

  awk -v c="$chr" -v s="$start" -v e="$end" 'BEGIN{OFS="\t"} ($1==c && $4>=s && $4<=e){print NR,$0}' \
    ${PLINK_PREFIX}.bim > "$LDIR/bim_window.tsv"

  python $SCRIPTS/harmonize_locus.py \
    --gwas "$GWAS" \
    --bim-window "$LDIR/bim_window.tsv" \
    --chr "$chr" \
    --start "$start" \
    --end "$end" \
    --out-prefix "$LDIR/$locus_id"

  nvar=$(wc -l < "$LDIR/$locus_id.extract_ids.txt")
  if [ "$nvar" -lt 2 ]; then
    echo "Skipping $locus_id: fewer than 2 harmonized variants"
    continue
  fi

  plink \
    --bfile "$PLINK_PREFIX" \
    --extract "$LDIR/$locus_id.extract_ids.txt" \
    --make-bed \
    --out "$LDIR/$locus_id.plink"

  python $SCRIPTS/check_order.py \
    --harmonized "$LDIR/$locus_id.harmonized.tsv" \
    --extracted-bim "$LDIR/$locus_id.plink.bim"

  plink \
    --bfile "$LDIR/$locus_id.plink" \
    --r square gz \
    --out "$LDIR/$locus_id"

  Rscript $SCRIPTS/run_susie_rss.R \
    --sumstats "$LDIR/$locus_id.harmonized.tsv" \
    --ld "$LDIR/$locus_id.ld.gz" \
    --out-prefix "$OUT/results/$locus_id" \
    --L 10 \
    --ref-n 10000

  echo "Done $locus_id"
done 2>&1 | tee "$OUT/logs/run_all_loci.log"

