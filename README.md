# scVarSim

[![Python](https://img.shields.io/badge/python-3.9-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

**scVarSim** (Single-Cell Variant Simulator) generates realistic single-cell RNA-seq
reads with **explicit ground truth** for genomic variants (phased SNPs and short
indels), **allele-specific expression (ASE)**, post-transcriptional **RNA-editing**
events, and optional **COSMIC-guided somatic SNVs**. Trained directly on a real
scRNA-seq BAM, it learns empirical gene-expression heterogeneity, genome-wide read
coverage, splice structure, and position-specific sequencing-error profiles, then
embeds haplotype-resolved variants and probabilistic RNA editing into synthetic
reads. This enables controlled, reproducible benchmarking of variant-calling and
RNA-editing-detection methods at both single-cell and pseudobulk resolution.

Key features:

- **Phased, haplotype-aware variants** — SNP/indel haplotype comes from the real
  GIAB genotype (`GT`), not a random draw, so phase is preserved across a read.
- **Allele-specific expression** — each gene/intergene feature draws an allelic
  ratio from an empirical pool; reads are assigned to a haplotype accordingly.
- **Splice-aware** — read/exon-block structure (BED12) is taken from real reads, so
  editing sites in introns are excluded and within-read positions are exact.
- **RNA editing** — DARNED-style sites embedded at empirical editing levels.
- **Somatic SNVs (optional)** — COSMIC-guided, cell-type-specific, VAF-controlled.
- **Sequencing errors** — position-specific error model via fgbio.
- **Ground truth for everything** — a combined `variants.vcf.gz`, per-site editing
  CSVs, per-cell editing matrices, and a somatic-mutation VCF/TSV.

> **Reproducibility first.** The simulation engine (`scvarsim/engine.py`) is a
> **frozen, byte-for-byte** copy of the research code — it is *not* refactored, so
> results do not change. See [docs/REPRODUCIBILITY.md](docs/REPRODUCIBILITY.md).

---

## Repository layout

```
scVarSim/
├── scvarsim/
│   ├── engine.py            # frozen simulation engine (do not edit)
│   └── __init__.py          # re-exports the engine's public functions
├── examples/
│   ├── run_simulation.py           # config-driven runner (portable entry point)
│   ├── config.chr19_GM12878.yaml   # example config (chr19 GM12878)
│   ├── reproduce_chr19_GM12878.py  # verbatim manuscript driver (record)
│   └── submit_hoffman2.sh          # UGE/SGE (qsub) wrapper
├── preprocessing/           # scripts that build the input files (see docs/INPUTS.md)
├── docs/                    # INSTALL, TUTORIAL, INPUTS, OUTPUTS, REPRODUCIBILITY
├── environment.yml          # conda env (tools + R stack)
├── requirements.txt         # pinned python deps + scReadSim
└── pyproject.toml
```

---

## Installation (summary)

scVarSim depends on external CLI tools (samtools, bedtools, seqtk, fgbio, picard, R)
and the [`scReadSim`](https://pypi.org/project/scReadSim/) package (which uses rpy2 +
R/scDesign2). The tested path is a conda env plus pinned pip deps:

```bash
conda env create -f environment.yml     # tools + R stack -> env "scvarsim"
conda activate scvarsim
pip install -r requirements.txt          # pinned python deps + scReadSim
pip install -e .                         # install the scvarsim package (editable)
```

Full step-by-step instructions and troubleshooting (rpy2, scDesign2, Seurat, R
compilation issues) are in **[docs/INSTALL.md](docs/INSTALL.md)**.

---

## Quick start

```bash
# 1. Edit examples/config.chr19_GM12878.yaml to point at your inputs
#    (BAM with CB/UB tags, reference FASTA, barcodes, phased GIAB VCFs,
#     RNA-editing list, allelic-ratio .rds, COSMIC VCFs). See docs/INPUTS.md.

# 2. Run the pipeline (uses $NSLOTS cores under a scheduler, else 16)
python examples/run_simulation.py --config examples/config.chr19_GM12878.yaml
```

On a UGE/SGE cluster: `qsub examples/submit_hoffman2.sh`.

The main outputs are the paired FASTQs
`*.EditingIncluded_ErrorIncluded.read{1,2}.bed2fa.sorted.fq` and the ground-truth
`variants.vcf.gz`. Every output file is described in
**[docs/OUTPUTS.md](docs/OUTPUTS.md)**.

---

## Python API

The engine functions are re-exported at the package top level (equivalently
`from scvarsim import engine as scIsoSim`):

```python
import scvarsim as scIsoSim

scIsoSim.prepare_SNP_RED_list_phased(snp_vcf, indel_vcf, red_site, outdir,
                                     min_edit=0.01, max_edit=0.99)
scIsoSim.assign_read_alleles_ASE(outdir, combined_prename, allelic_ratio_rds,
                                 per_cell=False, seed=2023)
# ... see examples/run_simulation.py for the full ordered pipeline.
```

The end-to-end pipeline also calls `scReadSim` (`Utility`, `GenerateSyntheticCount`,
`scRNA_GenerateBAM`) for feature sets, count matrices, and synthetic-count training —
`examples/run_simulation.py` shows the exact sequence.

---

## Documentation

| Doc | Contents |
|:----|:---------|
| [docs/INSTALL.md](docs/INSTALL.md) | conda env, scReadSim/rpy2/R/scDesign2 setup, troubleshooting |
| [docs/TUTORIAL.md](docs/TUTORIAL.md) | end-to-end chr19 GM12878 walkthrough |
| [docs/INPUTS.md](docs/INPUTS.md) | every input file: format + how to produce it (`preprocessing/`) |
| [docs/OUTPUTS.md](docs/OUTPUTS.md) | every output file explained |
| [docs/REPRODUCIBILITY.md](docs/REPRODUCIBILITY.md) | frozen engine, pinned versions, seeds, determinism caveats |

---

## Citation

If you use scVarSim, please cite the manuscript (in preparation). A BibTeX entry
will be added here on publication.

## License

MIT — see [LICENSE](LICENSE).

## Contact

Please open an issue at <https://github.com/gxiaolab/scVarSim/issues>.
