# Output files

All outputs are written to `outdirectory` (from the config). In the example,
`filename: GM12878`, so the combined read prefix is
`GM12878.syntheticBAM.combined` — written `<combined>` below.

## Primary outputs

| File | Description |
|:-----|:------------|
| `<combined>.EditingIncluded_ErrorIncluded.read1.bed2fa.sorted.fq` | **Final Read 1** — cell barcode (16 bp) + UMI (10 bp). |
| `<combined>.EditingIncluded_ErrorIncluded.read2.bed2fa.sorted.fq` | **Final Read 2** — genomic sequence carrying germline SNPs/indels, RNA edits, somatic SNVs, and sequencing errors. |
| `variants.vcf.gz` (+ `.tbi`) | **Ground-truth VCF** — combined germline SNVs + indels actually introduced (bgzipped, tabix-indexed). |
| `somatic_mutations.vcf` | Ground-truth somatic SNVs with `CLASS`/`FEATURE`/`CELLTYPE`/`VAF`/`NCELLS` INFO (only if `somatic_mutation: true`). |
| `somatic_mutations_ground_truth.tsv` | Per-somatic-SNV detail: feature, allele, cell type, target/realized carrier-cell counts, pseudobulk depth, carrier barcodes. |

## Ground-truth editing tables

| File | Description |
|:-----|:------------|
| `ErrorIncluded_Synthetic_RNA_editing_events.csv` | RNA-editing events realized in the error-included reads (the ground truth to score editing callers against). |
| `ErrorIncluded_edited_events_by_cell_matrix.csv` | Per-cell × per-site **edited**-read counts. |
| `ErrorIncluded_unedited_events_by_cell_matrix.csv` | Per-cell × per-site **unedited**-read counts. |
| `Synthetic_RNA_editing_events.csv` | Editing events before the error pass (intermediate). |

## Variant / allele intermediate tables (also useful as ground truth)

| File | Description |
|:-----|:------------|
| `SNPs.txt` | Germline SNPs with haplotype label `allelic` (`homo`/`0`/`1`), duplicated per strand. |
| `indels.txt` | Germline indels with haplotype label. |
| `RNA_editing_sites_levels.txt` | RNA-editing sites kept (SNP/indel-overlapping sites removed), with editing level. |
| `feature_allelic_ratios.txt` | Per-feature sampled allelic ratio + read counts (**ASE ground truth**). |
| `<combined>.read_alleles.txt` | Per-read haplotype assignment (`read_name`, `allele` 0/1). |
| `read_editing_positions.csv` | Per-read within-read positions of every germline/editing/somatic change (post indel-shift). |
| `<combined>.introduced_indel.csv` | Indels actually introduced, per read. |

## scReadSim intermediates

| File | Description |
|:-----|:------------|
| `scReadSim.Gene.bed`, `scReadSim.InterGene.bed` | Foreground (gene) and background (intergene) feature sets. |
| `<filename>.gene.countmatrix`, `<filename>.intergene.countmatrix` | Real UMI count matrices from the input BAM. |
| `*.LouvainClusterResults.txt` | Louvain clustering of the real cells (cell-type labels for training). |
| `*.scDesign2Simulated.txt` | Synthetic count matrices from scDesign2. |
| `*.scDesign2Simulated.CellTypeLabel.txt` | Synthetic-cell → cell-type labels. |
| `synthetic_cell_barcode.txt` | Synthetic cell barcodes. |
| `synthetic_cell_barcode.txt.withSynthCluster` | Synthetic barcodes + assigned cluster/cell type (used by the somatic step). |
| `<gene\|intergene>.syntheticBAM.*.read.bed` / `.read.bed12` | Synthetic read coordinates (BED6 and splice-aware BED12). |
| `<combined>.ErrorIncluded.read{1,2}.bed2fa.sorted.fq` | Reads with sequencing errors, before editing/variants are applied (intermediate). |

## Other

- A sequence dictionary `*.dict` is created next to the reference FASTA (needed by
  fgbio/htsjdk). It is written once, atomically.
- On first use, an allelic-ratio cache `<allelic_ratio_rds>.npy` is written next to
  the input `.rds`.

> **Note.** The pipeline only produces output when `embed_seq_error: true` (the
> error-free path was removed upstream). The scored ground truth for reads is the
> `EditingIncluded_ErrorIncluded` FASTQ pair + `variants.vcf.gz` + the
> `ErrorIncluded_*` editing CSVs + the somatic VCF/TSV.
