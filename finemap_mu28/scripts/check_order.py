import argparse
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("--harmonized", required=True)
parser.add_argument("--extracted-bim", required=True)
args = parser.parse_args()

h = pd.read_csv(args.harmonized, sep="\t")
b = pd.read_csv(
    args.extracted_bim,
    sep="\t",
    header=None,
    names=["CHR", "ID", "CM", "POS", "A1", "A2"]
)

if len(h) != len(b):
    raise SystemExit(f"Length mismatch: harmonized={len(h)} extracted_bim={len(b)}")

if list(h["BIM_ID"].astype(str)) != list(b["ID"].astype(str)):
    bad = [(i, h.loc[i, "BIM_ID"], b.loc[i, "ID"]) for i in range(min(len(h), len(b))) if str(h.loc[i, "BIM_ID"]) != str(b.loc[i, "ID"])]
    raise SystemExit(f"Order mismatch. First mismatches: {bad[:10]}")

print(f"Order OK: {len(h)} variants")
