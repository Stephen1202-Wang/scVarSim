# Installation

scVarSim has three dependency layers:

1. **CLI tools + R stack** (conda): samtools, bedtools, seqtk, gffread, fgbio,
   picard, gatk4, and R (>=4.2) with the packages `scReadSim` / `scDesign2` need.
2. **Python packages** (pip, pinned): numpy, pandas, pysam, biopython, vcfpy,
   gffpandas, joblib, tqdm, plus `scReadSim` and PyYAML.
3. **The `scvarsim` package itself** (this repo).

The versions below are those used to produce the manuscript results (see
[REPRODUCIBILITY.md](REPRODUCIBILITY.md)). `numpy==1.26.4` is the one pin that
matters for matching random draws.

---

## 1. Create the conda environment

```bash
conda env create -f environment.yml
conda activate scvarsim
```

This installs the CLI tools (samtools 1.21, bedtools 2.29.1, seqtk, gffread,
fgbio, picard 3.4.0, gatk4) and the R 4.2.3 stack + rpy2.

> **Make sure the tools are on `$PATH`.** The engine locates `samtools`,
> `bedtools`, `seqtk`, `Rscript`, and `java` via `shutil.which(...)`. Activating
> the conda env puts them on `$PATH`; if you wrap the run in a batch script, keep
> the env active (see `examples/submit_hoffman2.sh`).

## 2. Install the pinned python dependencies

```bash
pip install -r requirements.txt
```

This installs the pinned python packages **and** `scReadSim==1.4.1`.

## 3. Install scVarSim

```bash
pip install -e .        # from the repo root
```

Verify:

```bash
python -c "import scvarsim; print(scvarsim.__version__)"
python -c "import scvarsim; scvarsim.prepare_SNP_RED_list_phased"   # no error = OK
```

---

## R / scDesign2 (required by scReadSim)

`scReadSim.GenerateSyntheticCount` trains an [scDesign2](https://github.com/JSB-UCLA/scDesign2)
model in R via rpy2. Install scDesign2 once into the env's R:

```bash
R --vanilla -e 'if(!require("scDesign2")) { \
  install.packages("remotes", repos="https://cloud.r-project.org"); \
  remotes::install_github("JSB-UCLA/scDesign2", upgrade="never") }'
R --vanilla -e 'library(scDesign2); packageVersion("scDesign2")'
```

The allelic-ratio `.rds` loader (`_load_allelic_ratio_pool`) shells out to
`Rscript` once to export the vector to a cached `.npy`; no rpy2 is needed for that
step, only an `Rscript` on `$PATH`.

---

## Troubleshooting

### rpy2 import error / segfault
Symptoms: `libicuuc.so.*: cannot open shared object file`,
`cannot import name 'SexpVectorCCompatibleAbstract'`, or a segfault on import.
Install rpy2 from conda (not pip) so it links the env's R:

```bash
pip uninstall -y rpy2
conda install -c conda-forge rpy2 --force-reinstall -y
```

If you must use pip, set `R_HOME` first: `export R_HOME=$(R RHOME)`.

### rpy2 can't find R at runtime (batch jobs)
Export R locations before launching (mirrors `examples/submit_hoffman2.sh`):

```bash
export R_HOME="${CONDA_PREFIX}/lib/R"
export LD_LIBRARY_PATH="${R_HOME}/lib:${CONDA_PREFIX}/lib:${LD_LIBRARY_PATH:-}"
```

### scDesign2 install hits GitHub API rate limit (HTTP 403)
Install from a downloaded archive instead:

```bash
cd /tmp && curl -L -o scDesign2.zip \
  https://github.com/JSB-UCLA/scDesign2/archive/refs/heads/master.zip
unzip scDesign2.zip && R CMD INSTALL scDesign2-master
```

### R prompts "Update all/some/none?" in a non-interactive run
Disable update prompts: `remotes::install_github(..., upgrade="never")` and set
`options(repos = c(CRAN = "https://cloud.r-project.org"))`.

### R package compilation errors (gert / Rsubread / systemfonts)
`undefined symbol: getentropy`, `conflicting types for 'msgqu_init'`, or a missing
`systemfonts` are glibc/compiler mismatches — install those packages via conda
rather than compiling from source:

```bash
conda install -c conda-forge r-devtools r-gert r-systemfonts -y
conda install -c bioconda bioconductor-rsubread -y
```

### `zlib.h: No such file or directory` when R compiles a package
```bash
conda install -c conda-forge zlib -y
```

### Java heap space (picard / fgbio / gatk)
```bash
export JAVA_OPTS="-Xmx16g"
```

### Missing CB/UB tags in the BAM
The input BAM must carry cell-barcode (`CB`) and UMI (`UB`) tags:

```bash
samtools view input.bam | head -1 | tr '\t' '\n' | grep -E "^(CB|UB):"
```

### Seurat v5 vs v4
scReadSim expects Seurat v4-style assays. If you hit v5 assay errors, force v3
behavior in R (`options(Seurat.object.assay.version = "v3")`) or install Seurat 4:
`remotes::install_version("Seurat", version = "4.4.0")`.
