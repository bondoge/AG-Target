#!/usr/bin/env python3

import argparse
import gzip
from bisect import bisect_right
from pathlib import Path

import pandas as pd
import pysam


def open_text(path):
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def clean_chrom(x):
    s = str(x).strip()

    if s.endswith(".0"):
        s = s[:-2]

    if s.lower().startswith("chr"):
        chrom = "chr" + s[3:]
    else:
        chrom = "chr" + s

    if chrom == "chr23":
        return "chrX"
    if chrom == "chr24":
        return "chrY"
    if chrom in {"chr25", "chrMT"}:
        return "chrM"

    return chrom


def parse_pos(x):
    s = str(x).strip()
    if s.endswith(".0"):
        s = s[:-2]
    return int(s)


def clean_allele(x):
    if pd.isna(x):
        return None

    a = str(x).strip().upper()

    if a in {"", ".", "NA", "NAN"}:
        return None

    return a


def is_dna_allele(a):
    if a is None:
        return False
    return all(base in {"A", "C", "G", "T", "N"} for base in a)


def fasta_chrom_name(fasta, chrom):
    refs = set(fasta.references)

    if chrom in refs:
        return chrom

    no_chr = chrom[3:] if chrom.startswith("chr") else chrom
    if no_chr in refs:
        return no_chr

    if chrom == "chrM":
        for candidate in ["chrM", "MT", "M"]:
            if candidate in refs:
                return candidate

    raise ValueError(f"Chromosome {chrom} not found in FASTA.")


def fetch_ref(fasta, chrom, pos_1based, length):
    ref_chrom = fasta_chrom_name(fasta, chrom)
    start0 = pos_1based - 1
    end0 = start0 + length
    return fasta.fetch(ref_chrom, start0, end0).upper()


def parse_allele_pairs(s):
    """
    Example:
      EA,NEA;BIM_A1,BIM_A2
    """
    pairs = []

    for part in s.split(";"):
        part = part.strip()
        if not part:
            continue

        a, b = part.split(",")
        pairs.append((a.strip(), b.strip()))

    return pairs


def resolve_ref_alt_from_pos_only(row, fasta, chrom_col, pos_col, allele_pairs):
    """
    Uses ONLY:
      - CHR/POS from table
      - allele columns from table
      - hg38 FASTA

    Does NOT parse BIM_ID or GWAS_ID.
    """

    chrom = clean_chrom(row[chrom_col])
    pos = parse_pos(row[pos_col])

    checked = []

    for a_col, b_col in allele_pairs:
        if a_col not in row.index or b_col not in row.index:
            checked.append(f"missing allele columns: {a_col}/{b_col}")
            continue

        a1 = clean_allele(row[a_col])
        a2 = clean_allele(row[b_col])

        if not is_dna_allele(a1) or not is_dna_allele(a2):
            checked.append(f"{a_col}/{b_col}: invalid allele {a1}/{a2}")
            continue

        if a1 == a2:
            checked.append(f"{a_col}/{b_col}: identical alleles {a1}/{a2}")
            continue

        try:
            ref_for_a1_len = fetch_ref(fasta, chrom, pos, len(a1))
            ref_for_a2_len = fetch_ref(fasta, chrom, pos, len(a2))
        except Exception as e:
            checked.append(f"{a_col}/{b_col}: FASTA fetch failed: {e}")
            continue

        if ref_for_a1_len == a1:
            return {
                "ok": True,
                "CHROM": chrom,
                "POS": pos,
                "REF": a1,
                "ALT": a2,
                "allele_source": f"{a_col}/{b_col}",
                "refalt_status": "matched_first_allele_as_REF",
            }

        if ref_for_a2_len == a2:
            return {
                "ok": True,
                "CHROM": chrom,
                "POS": pos,
                "REF": a2,
                "ALT": a1,
                "allele_source": f"{a_col}/{b_col}",
                "refalt_status": "matched_second_allele_as_REF",
            }

        checked.append(
            f"{a_col}/{b_col} {a1}/{a2}: "
            f"hg38_ref_len_{len(a1)}={ref_for_a1_len}, "
            f"hg38_ref_len_{len(a2)}={ref_for_a2_len}"
        )

    return {
        "ok": False,
        "failure_reason": "no table allele matched hg38 reference at CHR/POS",
        "checked": " | ".join(checked),
    }


def parse_gtf_attrs(attr_text):
    attrs = {}

    for part in attr_text.strip().split(";"):
        part = part.strip()

        if not part:
            continue

        if " " in part:
            key, value = part.split(" ", 1)
            attrs[key] = value.strip().strip('"')
        elif "=" in part:
            key, value = part.split("=", 1)
            attrs[key] = value.strip().strip('"')

    return attrs


def load_gene_index(gtf_path, protein_coding_only=False):
    genes = []

    with open_text(gtf_path) as f:
        for line in f:
            if not line or line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")

            if len(fields) < 9:
                continue

            chrom, source, feature, start, end, score, strand, frame, attrs_text = fields

            if feature != "gene":
                continue

            attrs = parse_gtf_attrs(attrs_text)

            gene_id = attrs.get("gene_id", "")
            gene_name = attrs.get("gene_name", gene_id)
            gene_type = (
                attrs.get("gene_type")
                or attrs.get("gene_biotype")
                or attrs.get("biotype")
                or ""
            )

            if protein_coding_only and gene_type != "protein_coding":
                continue

            start = int(start)
            end = int(end)
            chrom = clean_chrom(chrom)

            tss = start if strand == "+" else end

            genes.append(
                {
                    "chrom": chrom,
                    "start": start,
                    "end": end,
                    "tss": tss,
                    "strand": strand,
                    "gene_id": gene_id,
                    "gene_name": gene_name,
                    "gene_type": gene_type,
                }
            )

    if not genes:
        raise ValueError(f"No genes loaded from GTF: {gtf_path}")

    by_chrom = {}

    for chrom, gdf in pd.DataFrame(genes).groupby("chrom"):
        gdf = gdf.sort_values("start").reset_index(drop=True)

        starts = gdf["start"].tolist()
        ends = gdf["end"].tolist()

        prefix_max_end = []
        prefix_max_idx = []

        best_end = -1
        best_idx = -1

        for i, end in enumerate(ends):
            if end > best_end:
                best_end = end
                best_idx = i

            prefix_max_end.append(best_end)
            prefix_max_idx.append(best_idx)

        by_chrom[chrom] = {
            "df": gdf,
            "starts": starts,
            "prefix_max_end": prefix_max_end,
            "prefix_max_idx": prefix_max_idx,
        }

    return by_chrom


def nearest_gene(gene_index, chrom, pos):
    if chrom not in gene_index:
        return {
            "nearest_gene": "",
            "nearest_gene_id": "",
            "nearest_gene_type": "",
            "nearest_gene_distance_bp": "",
            "nearest_gene_relation": "no_genes_on_chrom",
            "nearest_gene_start": "",
            "nearest_gene_end": "",
            "nearest_gene_strand": "",
        }

    idx = gene_index[chrom]
    gdf = idx["df"]
    starts = idx["starts"]

    insert_i = bisect_right(starts, pos)

    candidates = []

    if insert_i > 0:
        prev_idx = idx["prefix_max_idx"][insert_i - 1]
        gene = gdf.iloc[prev_idx]

        if int(gene["start"]) <= pos <= int(gene["end"]):
            dist = 0
            relation = "inside_gene"
        else:
            dist = pos - int(gene["end"])
            relation = "nearest_by_gene_body"

        tss_dist = abs(pos - int(gene["tss"]))
        candidates.append((dist, tss_dist, relation, gene))

    if insert_i < len(gdf):
        gene = gdf.iloc[insert_i]

        if int(gene["start"]) <= pos <= int(gene["end"]):
            dist = 0
            relation = "inside_gene"
        else:
            dist = int(gene["start"]) - pos
            relation = "nearest_by_gene_body"

        tss_dist = abs(pos - int(gene["tss"]))
        candidates.append((dist, tss_dist, relation, gene))

    if not candidates:
        return {
            "nearest_gene": "",
            "nearest_gene_id": "",
            "nearest_gene_type": "",
            "nearest_gene_distance_bp": "",
            "nearest_gene_relation": "no_candidate_gene",
            "nearest_gene_start": "",
            "nearest_gene_end": "",
            "nearest_gene_strand": "",
        }

    candidates.sort(key=lambda x: (x[0], x[1]))
    dist, tss_dist, relation, gene = candidates[0]

    return {
        "nearest_gene": gene["gene_name"],
        "nearest_gene_id": gene["gene_id"],
        "nearest_gene_type": gene["gene_type"],
        "nearest_gene_distance_bp": int(dist),
        "nearest_gene_relation": relation,
        "nearest_gene_start": int(gene["start"]),
        "nearest_gene_end": int(gene["end"]),
        "nearest_gene_strand": gene["strand"],
    }


def make_variant_ids(df):
    return (
        df["CHROM"].astype(str)
        + "_"
        + df["POS"].astype(str)
        + "_"
        + df["REF"].astype(str)
        + "_"
        + df["ALT"].astype(str)
        + "_b38"
    )


def process_file(
    input_path,
    output_path,
    failures_path,
    fasta,
    gene_index,
    chrom_col,
    pos_col,
    allele_pairs,
    on_fail,
):
    df = pd.read_csv(input_path, sep="\t", dtype=str)

    if chrom_col not in df.columns:
        raise ValueError(f"Missing chromosome column: {chrom_col}")

    if pos_col not in df.columns:
        raise ValueError(f"Missing position column: {pos_col}")

    ok_rows = []
    bad_rows = []

    for i, row in df.iterrows():
        resolved = resolve_ref_alt_from_pos_only(
            row=row,
            fasta=fasta,
            chrom_col=chrom_col,
            pos_col=pos_col,
            allele_pairs=allele_pairs,
        )

        original = row.to_dict()

        # Preserve original CHR/POS under explicit names.
        original["input_CHR"] = original.get(chrom_col, "")
        original["input_POS"] = original.get(pos_col, "")

        # Remove original CHR/POS so AG CHROM/POS are unambiguous.
        original.pop(chrom_col, None)
        original.pop(pos_col, None)

        if not resolved["ok"]:
            bad = dict(original)
            bad["failure_reason"] = resolved["failure_reason"]
            bad["checked"] = resolved["checked"]
            bad_rows.append(bad)

            if on_fail == "fail":
                raise ValueError(
                    f"Failed row {i} in {input_path}:\n"
                    f"{resolved['failure_reason']}\n"
                    f"{resolved['checked']}"
                )

            if on_fail == "drop":
                continue

            if on_fail == "keep":
                out = {
                    "CHROM": clean_chrom(row[chrom_col]),
                    "POS": parse_pos(row[pos_col]),
                    "REF": "",
                    "ALT": "",
                    "canonical_variant": "",
                    "nearest_gene": "",
                    "nearest_gene_id": "",
                    "nearest_gene_type": "",
                    "nearest_gene_distance_bp": "",
                    "nearest_gene_relation": "",
                    "allele_source": "",
                    "refalt_status": "FAILED_REFALT_RESOLUTION",
                    **original,
                }
                ok_rows.append(out)

            continue

        gene = nearest_gene(
            gene_index=gene_index,
            chrom=resolved["CHROM"],
            pos=resolved["POS"],
        )

        out = {
            "CHROM": resolved["CHROM"],
            "POS": resolved["POS"],
            "REF": resolved["REF"],
            "ALT": resolved["ALT"],
            "canonical_variant": (
                f"{resolved['CHROM']}:{resolved['POS']}:"
                f"{resolved['REF']}>{resolved['ALT']}"
            ),
            **gene,
            "allele_source": resolved["allele_source"],
            "refalt_status": resolved["refalt_status"],
            **original,
        }

        ok_rows.append(out)

    out_df = pd.DataFrame(ok_rows)
    bad_df = pd.DataFrame(bad_rows)

    if len(out_df) > 0:
        out_df.insert(0, "variant_id", make_variant_ids(out_df))

        # If duplicates exist across loci/credible sets, keep rows but make variant_id unique.
        duplicated = out_df["variant_id"].duplicated(keep=False)

        if duplicated.any():
            out_df.loc[duplicated, "variant_id"] = (
                out_df.loc[duplicated, "variant_id"]
                + "_"
                + out_df.loc[duplicated].groupby("variant_id").cumcount().add(1).astype(str)
            )

        first_cols = [
            "variant_id",
            "CHROM",
            "POS",
            "REF",
            "ALT",
            "canonical_variant",
            "nearest_gene",
            "nearest_gene_id",
            "nearest_gene_type",
            "nearest_gene_distance_bp",
            "nearest_gene_relation",
            "allele_source",
            "refalt_status",
        ]

        remaining = [c for c in out_df.columns if c not in first_cols]
        out_df = out_df[first_cols + remaining]

    out_df.to_csv(output_path, sep="\t", index=False)
    bad_df.to_csv(failures_path, sep="\t", index=False)

    print()
    print(f"Input:          {input_path}")
    print(f"Output:         {output_path}")
    print(f"Failures:       {failures_path}")
    print(f"Valid variants: {len(out_df)}")
    print(f"Failed rows:    {len(bad_df)}")


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Convert SuSiE candidate TSV files to AlphaGenome hg38 input. "
            "Uses only CHR/POS and allele columns. Never parses BIM_ID/GWAS_ID."
        )
    )

    parser.add_argument(
        "inputs",
        nargs="+",
        help="Input TSV files.",
    )

    parser.add_argument(
        "--fasta",
        required=True,
        help="hg38 FASTA path.",
    )

    parser.add_argument(
        "--gtf",
        required=True,
        help="hg38 GTF file, for example GENCODE GTF.",
    )

    parser.add_argument(
        "--outdir",
        default=None,
        help="Output directory. Default: same directory as input.",
    )

    parser.add_argument(
        "--chrom-col",
        default="CHR",
        help="Input chromosome column. Default: CHR.",
    )

    parser.add_argument(
        "--pos-col",
        default="POS",
        help="Input position column. Default: POS.",
    )

    parser.add_argument(
        "--allele-pairs",
        default="EA,NEA;BIM_A1,BIM_A2",
        help=(
            "Allele column pairs to test against hg38 reference. "
            "Default: EA,NEA;BIM_A1,BIM_A2"
        ),
    )

    parser.add_argument(
        "--on-fail",
        choices=["drop", "fail", "keep"],
        default="drop",
        help="What to do if REF/ALT cannot be resolved. Default: drop.",
    )

    parser.add_argument(
        "--protein-coding-only",
        action="store_true",
        help="Use only protein-coding genes for nearest-gene annotation.",
    )

    args = parser.parse_args()

    fasta_path = Path(args.fasta)

    if not fasta_path.exists():
        raise FileNotFoundError(f"FASTA not found: {fasta_path}")

    fai_path = Path(str(fasta_path) + ".fai")

    if not fai_path.exists():
        print(f"Creating FASTA index: {fai_path}")
        pysam.faidx(str(fasta_path))

    fasta = pysam.FastaFile(str(fasta_path))

    print(f"Loading gene annotation: {args.gtf}")
    gene_index = load_gene_index(
        args.gtf,
        protein_coding_only=args.protein_coding_only,
    )

    allele_pairs = parse_allele_pairs(args.allele_pairs)

    for input_file in args.inputs:
        input_path = Path(input_file)

        if args.outdir is None:
            outdir = input_path.parent
        else:
            outdir = Path(args.outdir)
            outdir.mkdir(parents=True, exist_ok=True)

        stem = input_path.name
        if stem.endswith(".tsv"):
            stem = stem[:-4]

        output_path = outdir / f"{stem}.ag_hg38.tsv"
        failures_path = outdir / f"{stem}.ag_hg38.failures.tsv"

        process_file(
            input_path=input_path,
            output_path=output_path,
            failures_path=failures_path,
            fasta=fasta,
            gene_index=gene_index,
            chrom_col=args.chrom_col,
            pos_col=args.pos_col,
            allele_pairs=allele_pairs,
            on_fail=args.on_fail,
        )


if __name__ == "__main__":
    main()
