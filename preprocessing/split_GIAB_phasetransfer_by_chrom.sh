#!/bin/bash
# GRCh38 version of 20260616_split_GIAB_phasetransfer_by_chrom.sh.
# Split the GIAB HG001 PGandRTGphasetransfer truth-set VCF (GRCh38 build):
# per-chromosome -> separate snps/indels -> add AF tag -> bgzip + tabix index.
#
# NOTE: unlike the GRCh37 source, the GRCh38 VCF already uses chr-prefixed
# contig names (chr1..chrX), so the --rename-chrs step is dropped and regions
# are extracted with the chr prefix.
set -euo pipefail

BCFTOOLS=/u/home/w/weijian/.conda/envs/scRET/bin/bcftools
TABIX=/u/home/w/weijian/.conda/envs/scIsoSim/bin/tabix

GIAB_DIR=/u/project/gxxiao/weijian/project/scEditSim/Genomes/GIAB_NA12878_HG001/NISTv3.3.2_GRCh38
SRC=$GIAB_DIR/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.vcf.gz
OUT_DIR=$GIAB_DIR/by_chrom
PREFIX=HG001_GRCh38_GIAB_highconf_v.3.3.2_highconf_PGandRTGphasetransfer

CHROMS=( $(seq 1 22) X )

mkdir -p "$OUT_DIR"

for c in "${CHROMS[@]}"; do
    echo "=== chr${c} ==="
    base="$OUT_DIR/${PREFIX}_chr${c}"

    # 1. Extract this chromosome (chr-prefixed region in GRCh38 source)
    "$BCFTOOLS" view -r "chr${c}" -O z -o "${base}.vcf.gz" "$SRC"
    "$TABIX" -f -p vcf "${base}.vcf.gz"

    for vtype in snps indels; do
        # 2. Split by variant type
        "$BCFTOOLS" view -v "$vtype" -O z -o "${base}_${vtype}.vcf.gz" "${base}.vcf.gz"
        # 3. Add AF tag -> 4. recompress, then index
        "$BCFTOOLS" +fill-tags "${base}_${vtype}.vcf.gz" -O u -- -t AF \
            | "$BCFTOOLS" view -O z -o "${base}_${vtype}.tmp.vcf.gz"
        mv -f "${base}_${vtype}.tmp.vcf.gz" "${base}_${vtype}.vcf.gz"
        "$TABIX" -f -p vcf "${base}_${vtype}.vcf.gz"
        n=$("$BCFTOOLS" view --no-header "${base}_${vtype}.vcf.gz" | wc -l)
        echo "  ${vtype}: ${n} records -> ${base}_${vtype}.vcf.gz"
    done
done

echo "Done. Outputs in $OUT_DIR (prefix ${PREFIX})"
