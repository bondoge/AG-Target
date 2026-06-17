# AG-Target

## Data organization

```text
susie/
├── README.md
└── data/
    ├── ImputedGenotypes_merged_british_hg38.10k.nodup.bed              # PLINK genotype matrix (~91 GB)
    ├── ImputedGenotypes_merged_british_hg38.10k.nodup.bim              # Variant metadata (~1.1 GB)
    ├── ImputedGenotypes_merged_british_hg38.10k.nodup.fam              # Sample metadata (~245 KB)
    └── gwas_sumstat_lead/
        ├── clumped_hg19/
        ├── clumped_hg38/
        ├── lead_hg19/
        ├── lead_hg38/
        └── sumstats_hg38/
