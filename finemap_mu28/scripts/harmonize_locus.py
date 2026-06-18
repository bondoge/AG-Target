import argparse
import pandas as pd
import numpy as np

def clean_chr(x):
    return str(x).replace("chr", "").replace("CHR", "")

def clean_allele(x):
    return str(x).upper().strip()

def allele_key(a, b):
    a = clean_allele(a)
    b = clean_allele(b)
    return "|".join(sorted([a, b]))

parser = argparse.ArgumentParser()
parser.add_argument("--gwas", required=True)
parser.add_argument("--bim-window", required=True)
parser.add_argument("--chr", required=True)
parser.add_argument("--start", type=int, required=True)
parser.add_argument("--end", type=int, required=True)
parser.add_argument("--out-prefix", required=True)
args = parser.parse_args()

# Read GWAS parquet.
g = pd.read_parquet(args.gwas)

rename = {
    "rsid": "GWAS_ID",
    "chromosome": "CHR",
    "base_pair_location": "POS",
    "effect_allele": "EA",
    "other_allele": "NEA",
    "effect_allele_frequency": "EAF",
    "beta": "BETA",
    "standard_error": "SE",
    "p_value": "P",
    "n": "N",
}
missing = [c for c in rename if c not in g.columns]
if missing:
    raise ValueError(f"GWAS file missing columns: {missing}")

g = g.rename(columns=rename)
g["CHR"] = g["CHR"].map(clean_chr)
g["POS"] = g["POS"].astype(int)
g["EA"] = g["EA"].map(clean_allele)
g["NEA"] = g["NEA"].map(clean_allele)

# Filter GWAS to locus.
chrom = clean_chr(args.chr)
g = g[(g["CHR"] == chrom) & (g["POS"] >= args.start) & (g["POS"] <= args.end)].copy()

if g.empty:
    print("No GWAS variants in locus window.")
    pd.DataFrame().to_csv(args.out_prefix + ".harmonized.tsv", sep="\t", index=False)
    open(args.out_prefix + ".extract_ids.txt", "w").close()
    raise SystemExit(0)

g["ALLELE_KEY"] = [allele_key(a, b) for a, b in zip(g["EA"], g["NEA"])]
g["Z_ORIG"] = g["BETA"] / g["SE"]

# Read BIM window. We include BIM_ORDER from awk NR.
b = pd.read_csv(
    args.bim_window,
    sep="\t",
    header=None,
    names=["BIM_ORDER", "BIM_CHR", "BIM_ID", "CM", "BIM_POS", "BIM_A1", "BIM_A2"],
    dtype={"BIM_CHR": str, "BIM_ID": str}
)
b["BIM_CHR"] = b["BIM_CHR"].map(clean_chr)
b["BIM_POS"] = b["BIM_POS"].astype(int)
b["BIM_A1"] = b["BIM_A1"].map(clean_allele)
b["BIM_A2"] = b["BIM_A2"].map(clean_allele)
b["ALLELE_KEY"] = [allele_key(a, b2) for a, b2 in zip(b["BIM_A1"], b["BIM_A2"])]

# Merge by chr, pos, unordered allele pair.
m = g.merge(
    b,
    left_on=["CHR", "POS", "ALLELE_KEY"],
    right_on=["BIM_CHR", "BIM_POS", "ALLELE_KEY"],
    how="inner"
)

if m.empty:
    print("No harmonized variants after CHR/POS/allele matching.")
    pd.DataFrame().to_csv(args.out_prefix + ".harmonized.tsv", sep="\t", index=False)
    open(args.out_prefix + ".extract_ids.txt", "w").close()
    raise SystemExit(0)

# Remove exact duplicate BIM IDs if any.
m = m.drop_duplicates(subset=["BIM_ID"])

# Align Z and beta to BIM_A1.
same = (m["EA"] == m["BIM_A1"]) & (m["NEA"] == m["BIM_A2"])
flip = (m["EA"] == m["BIM_A2"]) & (m["NEA"] == m["BIM_A1"])

bad = ~(same | flip)
if bad.any():
    raise ValueError("Unexpected allele mismatch after allele-key merge.")

m["ORIENTATION"] = np.where(same, "same_EA_is_BIM_A1", "flip_EA_is_BIM_A2")
m["BETA_ALIGNED_TO_BIM_A1"] = np.where(same, m["BETA"], -m["BETA"])
m["Z_ALIGNED_TO_BIM_A1"] = np.where(same, m["Z_ORIG"], -m["Z_ORIG"])
m["EAF_ALIGNED_TO_BIM_A1"] = np.where(same, m["EAF"], 1 - m["EAF"])

# Sort by original BIM order so it should match PLINK output order.
m = m.sort_values("BIM_ORDER").reset_index(drop=True)

outcols = [
    "BIM_ID", "GWAS_ID",
    "CHR", "POS",
    "EA", "NEA",
    "BIM_A1", "BIM_A2",
    "ORIENTATION",
    "EAF", "EAF_ALIGNED_TO_BIM_A1",
    "BETA", "SE", "Z_ORIG",
    "BETA_ALIGNED_TO_BIM_A1", "Z_ALIGNED_TO_BIM_A1",
    "P", "N",
    "BIM_ORDER"
]

m[outcols].to_csv(args.out_prefix + ".harmonized.tsv", sep="\t", index=False)

with open(args.out_prefix + ".extract_ids.txt", "w") as f:
    for vid in m["BIM_ID"]:
        f.write(str(vid) + "\n")

print(f"GWAS variants in window: {len(g)}")
print(f"BIM variants in window: {len(b)}")
print(f"Harmonized variants: {len(m)}")
print(f"Same orientation: {same.sum()}")
print(f"Flipped orientation: {flip.sum()}")
