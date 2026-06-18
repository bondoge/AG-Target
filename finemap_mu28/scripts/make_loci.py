import argparse
import pandas as pd
import re

def safe_id(x):
    return re.sub(r'[^A-Za-z0-9_.-]+', '_', str(x))

parser = argparse.ArgumentParser()
parser.add_argument("--leads", required=True)
parser.add_argument("--window", type=int, default=1_000_000)
parser.add_argument("--out", required=True)
args = parser.parse_args()

leads = pd.read_csv(args.leads, sep="\t", compression="infer")

required = ["SNPID", "CHR", "POS", "EA", "NEA", "P", "N"]
missing = [c for c in required if c not in leads.columns]
if missing:
    raise ValueError(f"Lead file is missing columns: {missing}")

rows = []
for i, r in leads.iterrows():
    chrom = str(r["CHR"]).replace("chr", "")
    pos = int(r["POS"])
    start = max(1, pos - args.window)
    end = pos + args.window
    locus_id = f"chr{chrom}_{pos}_{safe_id(r['SNPID'])}"
    rows.append({
        "locus_id": locus_id,
        "lead_snpid": r["SNPID"],
        "chr": chrom,
        "lead_pos": pos,
        "start": start,
        "end": end,
        "lead_ea": r["EA"],
        "lead_nea": r["NEA"],
        "lead_p": r["P"],
        "N": r["N"],
    })

out = pd.DataFrame(rows)
out.to_csv(args.out, sep="\t", index=False)
print(out)
