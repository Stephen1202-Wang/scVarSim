# Quickstart — run the simulation

This assumes you **already have** all input files ready (a per-chromosome BAM with
`CB`/`UB` tags, reference FASTA, barcodes, phased GIAB SNP/indel VCFs, an RNA-editing
list, the allelic-ratio `.rds`, and — if using somatic SNVs — the COSMIC VCFs) and
that scVarSim is installed. If not, see [INSTALL.md](INSTALL.md) and
[INPUTS.md](INPUTS.md) first.

## 1. Activate the environment

```bash
conda activate scvarsim        # the env from environment.yml
```

## 2. Make a config pointing at your files

```bash
cp examples/config.template.yaml my_config.yaml
```

Open `my_config.yaml` and replace every `/path/to/...` with your file's path. The
keys you must set:

```yaml
chr: chr19                      # chromosome to simulate (must match your per-chrom files)
filename: MySample              # output-file prefix
read_len: 96                    # Read 2 length

tools:                          # point at your installed tools (also need them on $PATH)
  samtools_directory: "/path/to/env/bin"
  bedtools_directory: "/path/to/env/bin"
  seqtk_directory:    "/path/to/env/bin/"
  gffread_dir:        "/path/to/env/bin/"
  fgbio_jarfile:      "/path/to/fgbio.jar"
  picard_jarfile:     "/path/to/picard.jar"

inputs:                         # your input files
  bamfile:               "/path/to/mysample.chr19.bam"
  cells_barcode_file:    "/path/to/barcodes.tsv"
  reference_genome_file: "/path/to/genome.fa"
  genome_size_file:      "/path/to/genome.sizes"
  genome_annotation:     "/path/to/annotation.gff3"
  red_site:              "/path/to/rna_editing_sites.chr19.txt"
  snp_site:              "/path/to/giab.chr19_snps.phased_cov3.vcf.gz"
  indel_site:            "/path/to/giab.chr19_indels.phased_cov3.vcf.gz"
  allelic_ratio_rds:     "/path/to/allelic_ratio_all_pooled.rds"
  cosmic_coding_vcf:     "/path/to/Cosmic_GenomeScreensMutant.vcf.gz"
  cosmic_noncoding_vcf:  "/path/to/Cosmic_NonCodingVariants.vcf.gz"

outdirectory: "/path/to/output_dir"
```

To skip somatic SNVs, set `somatic_mutation: false` (the `cosmic_*` paths are then
ignored). Keep `embed_seq_error: true` — it is required to produce output.

## 3. Run

```bash
python examples/run_simulation.py --config my_config.yaml
```

It uses `$NSLOTS` cores under a scheduler, otherwise 16. On a UGE/SGE cluster, set
`CONFIG`/`CONDA_ENV` in `examples/submit_hoffman2.sh` and run `qsub examples/submit_hoffman2.sh`.

## 4. Results

Everything lands in `outdirectory`. The main files (full list in
[OUTPUTS.md](OUTPUTS.md)):

```bash
cd /path/to/output_dir
ls *EditingIncluded_ErrorIncluded.read1.bed2fa.sorted.fq   # final Read 1 (barcode+UMI)
ls *EditingIncluded_ErrorIncluded.read2.bed2fa.sorted.fq   # final Read 2 (cDNA + variants)
zcat variants.vcf.gz | grep -v '^##' | head                # ground-truth SNPs + indels
head Synthetic_RNA_editing_events.csv                      # ground-truth RNA edits
cat  somatic_mutations.vcf | grep -v '^##' | head          # ground-truth somatic SNVs
```

That's it — align Read 2, run your caller, and score against these ground-truth files.
