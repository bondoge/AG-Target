import glob
import pandas as pd

files = glob.glob("./data/finemap_mu28/results/*.susie_pip.tsv")
dfs = []
for f in files:
    d = pd.read_csv(f, sep="\t")
    d["source_file"] = f
    dfs.append(d)

allres = pd.concat(dfs, ignore_index=True)
allres["locus_id"] = allres["source_file"].str.split("/").str[-1].str.replace(".susie_pip.tsv", "", regex=False)
allres["CS"] = allres["CS"].fillna("").astype(str)

allres["credible_set_id"] = ""
mask = allres["CS"] != ""
allres.loc[mask, "credible_set_id"] = allres.loc[mask, "locus_id"] + "_" + allres.loc[mask, "CS"]

# Example input sets.
credible = allres[allres["CS"].fillna("") != ""].copy()
pip01 = allres[allres["PIP"] >= 0.01].copy()
pip05 = allres[allres["PIP"] >= 0.05].copy()

allres.to_csv("./data/finemap_mu28/results/alphagenome_input_all_finemapped_variants.tsv", sep="\t", index=False)
credible.to_csv("./data/finemap_mu28/results/alphagenome_input_susie_credible_sets.tsv", sep="\t", index=False)
pip01.to_csv("./data/finemap_mu28/results/alphagenome_input_susie_pip_ge_0.01.tsv", sep="\t", index=False)
pip05.to_csv("./data/finemap_mu28/results/alphagenome_input_susie_pip_ge_0.05.tsv", sep="\t", index=False)

print("credible-set variants:", credible.shape)
print("PIP >= 0.01 variants:", pip01.shape)
print("PIP >= 0.05 variants:", pip05.shape)

print("credible-set variants:", credible.shape)
print("credible sets:", credible["credible_set_id"].nunique())
print("PIP >= 0.01 variants:", pip01.shape)
print("PIP >= 0.05 variants:", pip05.shape)

print("\nTop credible sets by number of variants:")
print(
    credible["credible_set_id"]
    .value_counts()
    .head(20)
)
