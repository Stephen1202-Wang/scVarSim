"""
scVarSim -- Single-Cell Variant Simulator.

A data-driven simulator for single-cell RNA-seq that generates synthetic reads
with explicit ground truth for genomic variants (phased SNPs / short indels),
allele-specific expression (ASE), RNA-editing events, and optional COSMIC-guided
somatic SNVs.

The simulation engine lives in :mod:`scvarsim.engine`, which is a *frozen,
byte-for-byte* copy of the research code (``scEditSim_20240708.py``). It is NOT
refactored: the algorithm is preserved exactly so results stay reproducible.

Two equivalent ways to reach the engine functions::

    import scvarsim as scIsoSim          # scIsoSim.prepare_SNP_RED_list_phased(...)
    from scvarsim import engine as scIsoSim

See ``examples/run_simulation.py`` for a config-driven end-to-end pipeline and
``examples/reproduce_chr19_GM12878.py`` for the verbatim manuscript driver.
"""

from . import engine

# Public pipeline API -- the functions the end-to-end simulation calls, re-exported
# at the top level so `scvarsim.<fn>` works. These names live in engine.py unchanged.
from .engine import (
    prepare_SNP_RED_list_phased,
    scRNA_GenerateBAMCoord_spliced,
    combine_bed12,
    assign_read_alleles_ASE,
    prepare_editing_parallel_spliced_ASE,
    introduce_somatic_mutations_cosmic,
    scRNA_BED2FASTQ_spliced,
    process_indels_in_reads_spliced,
    scRNA_ErrorBase,
    mutate_fastq_sequences_parallel,
    output_vcf_from_csvs,
)

__version__ = "0.1.0"

__all__ = [
    "engine",
    "prepare_SNP_RED_list_phased",
    "scRNA_GenerateBAMCoord_spliced",
    "combine_bed12",
    "assign_read_alleles_ASE",
    "prepare_editing_parallel_spliced_ASE",
    "introduce_somatic_mutations_cosmic",
    "scRNA_BED2FASTQ_spliced",
    "process_indels_in_reads_spliced",
    "scRNA_ErrorBase",
    "mutate_fastq_sequences_parallel",
    "output_vcf_from_csvs",
    "__version__",
]
