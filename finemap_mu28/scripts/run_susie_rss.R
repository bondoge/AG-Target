suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(susieR)
})

option_list <- list(
  make_option("--sumstats", type="character"),
  make_option("--ld", type="character"),
  make_option("--L", type="integer", default=10),
  make_option("--out-prefix", type="character", dest="out_prefix"),
  make_option("--ref-n", type="integer", default=10000, dest="ref_n")
)

opt <- parse_args(OptionParser(option_list=option_list))

ss <- fread(opt$sumstats)
if (nrow(ss) < 2) {
  stop("Need at least 2 variants for LD fine-mapping.")
}

R <- as.matrix(fread(opt$ld, header=FALSE))

if (nrow(R) != nrow(ss) || ncol(R) != nrow(ss)) {
  stop(sprintf("LD dimension mismatch: LD is %d x %d, sumstats has %d variants",
               nrow(R), ncol(R), nrow(ss)))
}

if (any(!is.finite(R))) stop("LD matrix contains non-finite values.")
if (any(!is.finite(ss$Z_ALIGNED_TO_BIM_A1))) stop("Z contains non-finite values.")

# Basic LD sanity checks.
diag_dev <- max(abs(diag(R) - 1), na.rm=TRUE)
sym_dev <- max(abs(R - t(R)), na.rm=TRUE)

cat("Variants:", nrow(ss), "\n")
cat("Max diag deviation:", diag_dev, "\n")
cat("Max symmetry deviation:", sym_dev, "\n")

if (diag_dev > 1e-4) warning("LD diagonal is not very close to 1.")
if (sym_dev > 1e-4) warning("LD matrix is not very symmetric.")

z <- ss$Z_ALIGNED_TO_BIM_A1
n <- round(median(ss$N, na.rm=TRUE))

# Use external/reference LD, so keep estimate_residual_variance = FALSE.
# susieR docs recommend estimate_residual_variance=TRUE only for in-sample R.
args <- list(
  z = z,
  R = R,
  n = n,
  L = opt$L,
  coverage = 0.95,
  estimate_residual_variance = FALSE
)

# Newer susieR versions support finite-reference LD correction.
# Your LD panel has 10,000 samples, so pass ref-n if supported.
if ("R_finite" %in% names(formals(susie_rss))) {
  args$R_finite <- opt$ref_n
}

cat("Starting susie_rss()\n")
flush.console()

fit <- do.call(susie_rss, args)

cat("Finished susie_rss()\n")
cat("fit class:", class(fit), "\n")
cat("PIP length:", length(fit$pip), "\n")
flush.console()

pip <- fit$pip

cs_membership <- rep("", length(pip))
if (!is.null(fit$sets) && !is.null(fit$sets$cs)) {
  for (i in seq_along(fit$sets$cs)) {
    idx <- fit$sets$cs[[i]]
    cs_membership[idx] <- ifelse(
      cs_membership[idx] == "",
      paste0("CS", i),
      paste0(cs_membership[idx], ";CS", i)
    )
  }
}

out <- copy(ss)
out[, PIP := pip]
out[, CS := cs_membership]

setorder(out, -PIP)

cat("About to write files to:", opt$out_prefix, "\n")
flush.console()
fwrite(out, paste0(opt$out_prefix, ".susie_pip.tsv"), sep="\t")

saveRDS(fit, paste0(opt$out_prefix, ".susie_fit.rds"))

sink(paste0(opt$out_prefix, ".susie_summary.txt"))
cat("Input variants:", nrow(ss), "\n")
cat("GWAS N:", n, "\n")
cat("L:", opt$L, "\n")
cat("Reference LD N:", opt$ref_n, "\n")
cat("Max LD diag deviation:", diag_dev, "\n")
cat("Max LD symmetry deviation:", sym_dev, "\n")
cat("\nCredible sets:\n")
print(fit$sets)
cat("\nTop PIPs:\n")
print(head(out[, .(BIM_ID, GWAS_ID, CHR, POS, EA, NEA, BIM_A1, BIM_A2, P, Z_ALIGNED_TO_BIM_A1, PIP, CS)], 30))
sink()
