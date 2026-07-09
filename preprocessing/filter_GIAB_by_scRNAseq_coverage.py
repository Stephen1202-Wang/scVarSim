"""
Filter GIAB phased SNP and indel VCF files by scRNA-seq coverage.
Keeps only variants at positions covered by >= MIN_DEPTH reads
in the real scRNA-seq BAM (forward + reverse strand depth combined).

Input:  *_highconf_PGandRTGphasetransfer_{chrom}_{snps,indels}.phased.vcf.gz
Coverage: GM12878 NGS depth files, available for chr1, chr2, chr3, chr19
Output: GM12878_by_chrom_cov3/ subdirectory with .phased_cov3.vcf.gz files + tabix indices

conda env: scIsoSim
"""

import gzip
import os
import subprocess

GIAB_DIR   = "/u/project/gxxiao/weijian/project/scEditSim/Genomes/GIAB_NA12878_HG001/NISTv3.3.2_GRCh38/by_chrom"
COV_DIR    = "/u/project/gxxiao/weijian/project/scEditSim/data/coverage/GM12878"
OUT_DIR    = "/u/project/gxxiao/weijian/project/scEditSim/Genomes/GIAB_NA12878_HG001/NISTv3.3.2_GRCh38/by_chrom_cov3"
BGZIP      = "/u/home/w/weijian/.conda/envs/scIsoSim/bin/bgzip"
TABIX      = "/u/home/w/weijian/.conda/envs/scIsoSim/bin/tabix"
MIN_DEPTH  = 3
CHROMS     = ["chr1", "chr2", "chr3", "chr19"]

os.makedirs(OUT_DIR, exist_ok=True)


def load_depth(path):
    """Return {pos(int): depth(int)} from a 3-column depth file (chr, pos, depth, no header)."""
    depth = {}
    with open(path) as fh:
        for line in fh:
            parts = line.split()
            depth[int(parts[1])] = int(parts[2])
    return depth


def filter_vcf(in_vcf_gz, out_vcf_gz, fwd, rev, min_depth):
    """
    Stream through a gzipped VCF, keep header lines and variant lines
    where fwd[pos] + rev[pos] >= min_depth. Write bgzipped output.
    Returns (total_variants, kept_variants).
    """
    tmp_vcf = out_vcf_gz.replace(".vcf.gz", ".vcf")
    total = kept = 0

    with gzip.open(in_vcf_gz, "rt") as inp, open(tmp_vcf, "w") as out:
        for line in inp:
            if line.startswith("#"):
                out.write(line)
                continue
            total += 1
            pos = int(line.split("\t")[1])
            if fwd.get(pos, 0) + rev.get(pos, 0) >= min_depth:
                out.write(line)
                kept += 1

    subprocess.run([BGZIP, "-f", tmp_vcf], check=True)
    subprocess.run([TABIX, "-p", "vcf", out_vcf_gz], check=True)
    return total, kept


for chrom in CHROMS:
    print(f"\n=== {chrom} ===")

    fwd_path = os.path.join(COV_DIR, f"GRCh38_NGS_GM12878_read_{chrom}_forward_depth.txt")
    rev_path = os.path.join(COV_DIR, f"GRCh38_NGS_GM12878_read_{chrom}_reverse_depth.txt")

    print(f"  Loading coverage files...")
    fwd = load_depth(fwd_path)
    rev = load_depth(rev_path)
    print(f"  Forward positions: {len(fwd):,}  |  Reverse positions: {len(rev):,}")

    prefix = f"HG001_GRCh38_GIAB_highconf_v.3.3.2_highconf_PGandRTGphasetransfer_{chrom}"

    for vtype in ("snps", "indels"):
        in_vcf  = os.path.join(GIAB_DIR, f"{prefix}_{vtype}.phased.vcf.gz")
        out_vcf = os.path.join(OUT_DIR,  f"{prefix}_{vtype}.phased_cov{MIN_DEPTH}.vcf.gz")

        total, kept = filter_vcf(in_vcf, out_vcf, fwd, rev, MIN_DEPTH)
        pct = 100 * kept / total if total else 0
        print(f"  {vtype}: {total:>7,} -> {kept:>7,} kept ({pct:.1f}%)")

print("\nDone. Filtered VCFs written to:", OUT_DIR)
