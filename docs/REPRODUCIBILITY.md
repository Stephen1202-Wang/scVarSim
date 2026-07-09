# Reproducibility

scVarSim separates two ideas that are easy to conflate:

1. **Packaging fidelity** — does this package compute the *same thing* as the
   original research code? **Guaranteed.**
2. **Run-to-run determinism** — does re-running produce a *byte-identical* FASTQ?
   **Partial, by design of the original code** (some steps are intentionally
   stochastic). Not changed by packaging.

---

## 1. Packaging fidelity (guaranteed)

- **`scvarsim/engine.py` is a byte-for-byte copy** of the research engine
  (`scEditSim_20240708.py`). Verify at any time:

  ```bash
  sha256sum scvarsim/engine.py
  # b9c32495553c700ac72f42420803945eb7e887ce84f488f40e9838252aff75a3
  ```

- **`examples/run_simulation.py` calls the same functions, with the same
  arguments, in the same order, with the same seeds and toggles** as the verbatim
  manuscript driver `examples/reproduce_chr19_GM12878.py`. The differences are that
  paths/toggles come from a YAML config, the engine is imported as a package, and the
  engine alias is `scVarSim` (the record uses `scIsoSim`) — the same object. Confirm
  the ordered call sequence matches (normalizing the alias):

  ```bash
  norm() { sed -E 's/\b(scIsoSim|scVarSim)\./ENGINE./g' "$1" \
      | grep -oE '(ENGINE|Utility|GenerateSyntheticCount|scRNA_GenerateBAM)\.[A-Za-z_]+'; }
  diff <(norm examples/reproduce_chr19_GM12878.py) \
       <(norm examples/run_simulation.py) && echo "IDENTICAL call sequence"
  ```

No engine logic was modified. We deliberately did **not** add seeds or otherwise
"fix" determinism, because that would change results.

## 2. Pinned software

`numpy==1.26.4` is the reproducibility-critical pin: the seeded steps draw from
numpy's RNG, whose stream can differ across numpy versions. Other pins
(`requirements.txt`) and the conda tool versions (`environment.yml`) match the
environment used for the manuscript:

| Package | Version | | Tool | Version |
|:--------|:--------|-|:-----|:--------|
| numpy | 1.26.4 | | samtools | 1.21 |
| pandas | 2.2.0 | | bedtools | 2.29.1 |
| pysam | 0.23.3 | | seqtk | 1.4 |
| biopython | 1.85 | | gffread | 0.12.7 |
| vcfpy | 0.13.8 | | picard | 3.4.0 |
| gffpandas | 1.2.0 | | fgbio | 3.0.0 |
| joblib | 1.5.3 | | gatk4 | 4.6.2.0 |
| tqdm | 4.66.1 | | R | 4.2.3 |
| scReadSim | 1.4.1 | | python | 3.9.19 |

## 3. What is deterministic vs. stochastic

Given the **same inputs** and the pinned software:

**Deterministic** (fixed seed or no RNG):

| Step | Why |
|:-----|:----|
| `prepare_SNP_RED_list_phased` | no RNG — haplotype comes from the phased `GT`. |
| `assign_read_alleles_ASE` | `np.random.seed(2023)` + `random.seed(2023)`. |
| `introduce_somatic_mutations_cosmic` | `np.random.RandomState(2024)`. |

**Stochastic run-to-run** (unseeded global numpy RNG or per-worker `random`):

| Step | Source |
|:-----|:-------|
| synthetic count training (scReadSim `scRNA_GenerateSyntheticCount` → scDesign2) | external R model, not seeded here. |
| synthetic read coordinates / UMI sampling / ±jitter (`scRNA_GenerateBAMCoord_spliced`) | `np.random.choice/randint` on the global RNG. |
| RNA-editing realization (`random.random() < editing_level` in loky workers) | per-worker Python `random`, not seeded. |
| sequencing-error injection (`scRNA_ErrorBase`/`ErroneousRead`) | `np.random.binomial/choice` on the global RNG. |

**Consequence.** Re-running produces reads that are **statistically equivalent** but
not byte-identical. This is expected and inherited from the original code — it is not
introduced by packaging.

## 4. The ground truth is always self-consistent

Whatever a given run produced, its ground-truth files describe *that run exactly*:
`variants.vcf.gz`, `ErrorIncluded_Synthetic_RNA_editing_events.csv`, the per-cell
editing matrices, `feature_allelic_ratios.txt`, `<combined>.read_alleles.txt`, and
`somatic_mutations_ground_truth.tsv` are written from the same objects that generated
the reads. Always score a run against **its own** ground-truth files, not a previous
run's.

## 5. If you need a fully repeatable run

That is a deliberate, results-changing change and is intentionally **out of scope**
here (it would violate the "don't change results" requirement this package was built
under). It would involve seeding the global numpy RNG before read-coordinate
generation, seeding the loky workers, and fixing the scDesign2 seed inside scReadSim.
Do it as a separate, clearly-versioned modification if required.
