# Tutorial: chr19 GM12878, end to end

This walks through one chromosome (chr19) of the GM12878 example. The same steps
apply to any chromosome — change `chr` and the per-chromosome paths in the config.

## 0. Install

Follow [INSTALL.md](INSTALL.md):

```bash
conda env create -f environment.yml
conda activate scvarsim
pip install -r requirements.txt
pip install -e .
```

## 1. Gather / build the inputs

See [INPUTS.md](INPUTS.md) for formats. You need, for `chr19`:

- a real scRNA-seq BAM for chr19 with `CB`/`UB` tags (+ `.bai`),
- the cell `barcodes.tsv`,
- the reference FASTA (+ `.fai`), a chrom-sizes file, and a GENCODE GFF3,
- an RNA-editing site list with levels,
- **phased** GIAB SNP + indel VCFs for chr19 filtered to coverage ≥ 3,
- the allelic-ratio `.rds`,
- (optional) COSMIC coding + non-coding VCFs.

Build the phased GIAB VCFs and the allelic-ratio pool with the shipped scripts:

```bash
bash   preprocessing/split_GIAB_phasetransfer_by_chrom.sh
bash   preprocessing/keep_phased_only.sh
python preprocessing/filter_GIAB_by_scRNAseq_coverage.py
Rscript preprocessing/build_allelic_ratio_pool.R
# RNA-editing profile from the real BAM:
bash   preprocessing/rna_editing_snps_profile_from_real_data.sh
```

(Edit the paths at the top of each script first.)

## 2. Configure

Copy the example config and point it at your files:

```bash
cp examples/config.chr19_GM12878.yaml my_config.yaml
$EDITOR my_config.yaml         # set inputs.* , tools.* , outdirectory
```

Keep `chr: chr19`, `read_len: 96`, and both toggles `embed_seq_error: true` /
`somatic_mutation: true` (set `somatic_mutation: false` to skip COSMIC SNVs).

## 3. Run

Locally (uses `$NSLOTS` cores if set, else 16):

```bash
python examples/run_simulation.py --config my_config.yaml
```

On a UGE/SGE cluster, edit `CONFIG`/`CONDA_ENV` in `examples/submit_hoffman2.sh`
and submit:

```bash
qsub examples/submit_hoffman2.sh
```

Progress prints per stage: feature sets → count matrices → phased variant lists →
scDesign2 training → synthetic read coordinates → ASE allele assignment → RNA-editing
placement → (somatic SNVs) → FASTQ extraction → indels → sequencing errors → final
mutate → ground-truth VCF.

## 4. Inspect the outputs

In `outdirectory` (see [OUTPUTS.md](OUTPUTS.md) for the full list):

```bash
# Final synthetic reads:
ls *EditingIncluded_ErrorIncluded.read1.bed2fa.sorted.fq \
   *EditingIncluded_ErrorIncluded.read2.bed2fa.sorted.fq

# Ground-truth variants:
zcat variants.vcf.gz | grep -v '^##' | head
cat  somatic_mutations.vcf                | grep -v '^##' | head   # if enabled

# Ground-truth editing + ASE:
head Synthetic_RNA_editing_events.csv feature_allelic_ratios.txt
```

## 5. Downstream benchmarking

The two FASTQs are a standard paired-end 10x-style library (Read 1 = 16 bp barcode +
10 bp UMI; Read 2 = cDNA). Typical benchmark loop:

1. Align Read 2 (e.g. STARsolo with the barcodes) to the reference.
2. Run your variant caller / RNA-editing detector / somatic caller.
3. Compare calls to the ground truth:
   - SNPs + indels → `variants.vcf.gz`
   - RNA editing → `ErrorIncluded_Synthetic_RNA_editing_events.csv` and the per-cell
     `ErrorIncluded_{edited,unedited}_events_by_cell_matrix.csv`
   - somatic SNVs → `somatic_mutations.vcf` / `somatic_mutations_ground_truth.tsv`
   - allele-specific expression → `feature_allelic_ratios.txt` and
     `<combined>.read_alleles.txt`

Always score against **that run's** ground-truth files (see
[REPRODUCIBILITY.md](REPRODUCIBILITY.md) — runs are statistically equivalent but not
byte-identical).
