#!/bin/bash
# GRCh38 version of 20260616_keep_phased_only.sh.
# Keep only phased sites (GT uses '|') for each per-chromosome PGandRTGphasetransfer
# VCF (combined, snps, indels). Uses `bcftools view -p` (--phased), which keeps
# records where the sample genotype is phased and drops unphased ('/') ones.
# Writes new *.phased.vcf.gz files (originals untouched) + tabix indices.
set -euo pipefail

BCFTOOLS=/u/home/w/weijian/.conda/envs/scRET/bin/bcftools
TABIX=/u/home/w/weijian/.conda/envs/scIsoSim/bin/tabix

OUT_DIR=/u/project/gxxiao/weijian/project/scEditSim/Genomes/GIAB_NA12878_HG001/NISTv3.3.2_GRCh38/by_chrom
PREFIX=HG001_GRCh38_GIAB_highconf_v.3.3.2_highconf_PGandRTGphasetransfer

CHROMS=( $(seq 1 22) X )

for c in "${CHROMS[@]}"; do
    echo "=== chr${c} ==="
    for tag in "" "_snps" "_indels"; do
        in="$OUT_DIR/${PREFIX}_chr${c}${tag}.vcf.gz"
        out="$OUT_DIR/${PREFIX}_chr${c}${tag}.phased.vcf.gz"
        "$BCFTOOLS" view -p -O z -o "$out" "$in"
        "$TABIX" -f -p vcf "$out"
        tot=$("$BCFTOOLS" view --no-header "$in"  | wc -l)
        kept=$("$BCFTOOLS" view --no-header "$out" | wc -l)
        printf '  %-8s %d -> %d phased\n' "chr${c}${tag}:" "$tot" "$kept"
    done
done

echo "Done. Phased-only VCFs written with .phased.vcf.gz suffix in $OUT_DIR"
