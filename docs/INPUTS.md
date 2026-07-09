# Input files

Every path is set in the YAML config: copy `examples/config.template.yaml` (all
placeholder paths) and fill it in, or see `examples/config.chr19_GM12878.yaml` for a
concrete filled-in example.
Large inputs are **not** shipped with this repo; the `preprocessing/` scripts below
build them. All coordinates are 1-based genomic; the example uses GRCh38 with
`chr`-prefixed contig names.

| Config key | File | Format | Produced by |
|:-----------|:-----|:-------|:------------|
| `inputs.bamfile` | real scRNA-seq BAM (per-chrom) | BAM with `CB`/`UB` tags, indexed | your aligner (e.g. STARsolo) + split by chromosome |
| `inputs.cells_barcode_file` | cell barcodes | one barcode per line (`barcodes.tsv`) | STARsolo `filtered/` output |
| `inputs.reference_genome_file` | reference genome | FASTA + `.fai` (`samtools faidx`) | genome provider (GENCODE/UCSC) |
| `inputs.genome_size_file` | chromosome sizes | TSV `name<TAB>length` | `samtools faidx … ; cut -f1,2 *.fai` |
| `inputs.genome_annotation` | gene annotation | GFF3 (e.g. GENCODE basic) | genome provider |
| `inputs.red_site` | RNA-editing sites + levels | TSV, DARNED-style (see below) | `preprocessing/rna_editing_snps_profile_from_real_data.sh` |
| `inputs.snp_site` | phased SNPs (per-chrom, cov≥3) | bgzipped VCF + `.tbi`, phased `GT`, `AF` | `preprocessing/split_GIAB_phasetransfer_by_chrom.sh` → `keep_phased_only.sh` → `filter_GIAB_by_scRNAseq_coverage.py` |
| `inputs.indel_site` | phased indels (per-chrom, cov≥3) | bgzipped VCF + `.tbi`, phased `GT`, `AF` | same GIAB chain as SNPs |
| `inputs.allelic_ratio_rds` | allelic-ratio pool | R `.rds`, flat numeric vector in [0,1] | `preprocessing/build_allelic_ratio_pool.R` |
| `inputs.cosmic_coding_vcf` | COSMIC coding SNVs | bgzipped VCF, numeric contigs (`1`) | COSMIC download + `bcftools norm` |
| `inputs.cosmic_noncoding_vcf` | COSMIC non-coding SNVs | bgzipped VCF, numeric contigs | COSMIC download + `bcftools norm` |

---

## BAM, barcodes, reference, annotation

The BAM is a real scRNA-seq alignment carrying **`CB`** (cell barcode) and **`UB`**
(UMI) tags, split per chromosome (the pipeline runs one chromosome at a time via the
`chr` config key). scVarSim learns expression heterogeneity, coverage, splice
structure, and the sequencing-error profile from this BAM. Confirm the tags:

```bash
samtools view input.bam | head -1 | tr '\t' '\n' | grep -E "^(CB|UB):"
```

Build the chromosome-sizes file from the FASTA index:

```bash
samtools faidx GRCh38.p13.genome.fa
cut -f1,2 GRCh38.p13.genome.fa.fai > GRCh38.p13.genoms.sizes
```

## RNA-editing sites (`red_site`)

Tab-separated, read with a header. The **essential** columns are `Region`,
`Position`, `Ref`, `Ed`, `Strand`, and a **tissue** column holding the editing level
(default tissue name `Lung`). Extended DARNED-style annotation columns
(`Accession`, `type`, `Func.wgEncodeGencodeBasicV45`, `Gene.wgEncodeGencodeBasicV45`)
are used when present. Editing level is filtered to `[min_edit, max_edit]`
(the runner passes `min_edit=0.01, max_edit=0.99`).

```
Region   Position   Ref   Ed   Strand   Lung
chr19    1000123    A     G    +        0.42
chr19    1000456    A     G    -        0.31
```

Build the empirical mismatch profile from the real BAM with
`preprocessing/rna_editing_snps_profile_from_real_data.sh` (samtools mpileup → base
counts). That produces `mismatch_counts.tsv`; annotate/threshold it into the
DARNED-style `red_site` table with your editing-site catalog.

## Phased GIAB SNPs / indels (`snp_site`, `indel_site`)

Per-chromosome, **phased** (`GT` uses `|`), coverage-filtered VCFs with an `AF` INFO
tag. The haplotype label comes from the genotype: `0|1 → '1'`, `1|0 → '0'`,
`1|1 → homo`; multiallelic sites (`len(ALT)>1` or a `GT` referencing allele ≥ 2 like
`1|2`) are dropped. Build them in three steps:

```bash
# 1. Split the GIAB PGandRTGphasetransfer truth set per chromosome, split
#    snps/indels, add AF, bgzip + tabix:
bash preprocessing/split_GIAB_phasetransfer_by_chrom.sh

# 2. Keep only phased sites (GT with '|'):
bash preprocessing/keep_phased_only.sh

# 3. Keep only positions covered by >= 3 reads in the real scRNA-seq BAM
#    (writes the by_chrom_cov3/*.phased_cov3.vcf.gz used by the config):
python preprocessing/filter_GIAB_by_scRNAseq_coverage.py
```

Edit the paths at the top of each script for your layout. Source truth set:
[GIAB HG001/NA12878 NISTv3.3.2 (GRCh38)](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/NA12878_HG001/).

## Allelic-ratio pool (`allelic_ratio_rds`)

A flat numeric R vector of heterozygous allelic ratios `alt/(ref+alt)` in `(0,1)`,
pooled across samples. Each gene/intergene feature draws one ratio from this pool as
its haplotype-0 read fraction. Build (and plot the fit) with:

```bash
Rscript preprocessing/build_allelic_ratio_pool.R allelic_ratio_all_samples.pdf
# -> writes allelic_ratio_all_pooled.rds (cached), used by inputs.allelic_ratio_rds
```

On first use the engine converts the `.rds` to a memory-mapped `.npy` cache next to
it (`<rds>.npy`) via `Rscript`; keep an `Rscript` on `$PATH`.

## COSMIC catalogs (`cosmic_*`, optional)

Used only when `somatic_mutation: true`. Download the GRCh38
**Genome Screens Mutant** (coding) and **Non-Coding Variants** catalogs from
[COSMIC](https://cancer.sanger.ac.uk/cosmic/download), then normalize:

```bash
bcftools norm -f GRCh38.p13.genome.fa -m -any cosmic_raw.vcf.gz -Oz -o cosmic_norm.vcf.gz
```

The loader keeps single-base `ACGT` SNVs on the target chromosome (COSMIC uses
numeric contig names such as `1`; the pipeline maps them to `chr1`), and only sites
whose REF matches the reference genome and that have pseudobulk depth ≥
`min_pseudobulk_reads` (10).
