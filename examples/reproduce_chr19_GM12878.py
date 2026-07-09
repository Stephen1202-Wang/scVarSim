# code for scEditSim: scRNA-Seq Editing Simulator for simulate only one celltype with GRCh38
# Splicing-aware + HAPLOTYPE / ALLELE-SPECIFIC-EXPRESSION (ASE) version.
#   - SNPs/indels come from the PHASED GIAB VCFs (GM12878_by_chrom_cov3); each variant's
#     haplotype is taken from the real genotype (GT) instead of being random.
#   - Each gene/intergene gets an allelic ratio sampled from allelic_ratio_all_pooled.rds
#     (e.g. 0.8 -> 80% of that feature's reads from one haplotype, 20% from the other).
#   - Reads are assigned to a haplotype per that ratio, and variants are introduced so each
#     read carries exactly one haplotype's variants (plus homozygous), preserving phase.
# Differs from 20260602_..._Splicing.py by: phased variant inputs, three swapped calls
#     prepare_SNP_RED_list            -> prepare_SNP_RED_list_phased
#     assign_read_alleles             -> assign_read_alleles_ASE
#     prepare_editing_parallel_spliced-> prepare_editing_parallel_spliced_ASE
# plus an optional COSMIC somatic-mutation block (introduce_somatic_mutations_cosmic), and
# removal of the error-free output path (only error-included reads are produced).
# Conda ENV: scIsoSim
import scReadSim.Utility as Utility
import scReadSim.GenerateSyntheticCount as GenerateSyntheticCount
import scReadSim.scRNA_GenerateBAM as scRNA_GenerateBAM
import pkg_resources
import os
import pandas as pd
import numpy as np
import sys
os.chdir("/u/home/w/weijian/project-gxxiao/project/scEditSim/")
sys.path.append('/u/home/w/weijian/project-gxxiao/project/scEditSim/')
import scEditSim_20240708 as scIsoSim

# Parallel cores: under qsub this auto-matches the granted slots ($NSLOTS from
# `-pe shared N`); defaults to 16 when run outside a scheduler job.
NCORES = int(os.environ.get("NSLOTS", "16"))

#####################################################################
############################## User Input #########################
#####################################################################
chr = 'chr19'
samtools_directory="~/.conda/envs/scIsoSim/bin"
bedtools_directory="~/.conda/envs/scIsoSim/bin"
seqtk_directory = "~/.conda/envs/scIsoSim/bin/"
fgbio_jarfile = "/u/home/w/weijian/.conda/envs/scIsoSim/share/fgbio/fgbio.jar"
gffread_dir = "~/.conda/envs/scIsoSim/bin/"
picard_jarfile = "~/.conda/envs/scIsoSim/share/picard-3.4.0-0/picard.jar"

INPUT_cells_barcode_file = "/u/project/gxxiao/weijian/project/scEditSim/data/NGS_GM12878/GRCh38/STARsolo_alignment/GENCODE_Genome/GM12878.Solo.out/Gene/filtered/barcodes.tsv"
filename = "GM12878"
INPUT_bamfile = "/u/project/gxxiao/weijian/project/scEditSim/data/NGS_GM12878/GRCh38/STARsolo_alignment/GENCODE_Genome/split_by_chromosomes/GM12878_" + chr + ".CBattached.filtered.bam"
INPUT_genome_size_file = "/u/home/w/weijian/project-gxxiao/project/scEditSim/Genomes/GRCh38/GRCh38.p13.genoms.sizes"
INPUT_genome_annotation = "/u/home/w/weijian/project-gxxiao/project/scEditSim/Genomes/GRCh38/gencode.v42.primary_assembly.basic.annotation.gff3" # Put these into server FileGator for tutorial downloading

referenceGenome_file = "/u/home/w/weijian/project-gxxiao/project/scEditSim/Genomes/GRCh38/GRCh38.p13.genome.fa"


# SNPs from GIAB (phased) and RNA editing sites from DARNED
read_len = 96
INPUT_RED_site = "/u/project/gxxiao/weijian/project/scEditSim/result/Real_data_GRCh38/GM12878/" + chr +"/DARNED_user_list_20260618_" + chr + "_3_rc_3.txt"
outdirectory = "/u/home/w/weijian/project-gxxiao/project/scEditSim/result/Synthetic_data/GRCh38_20260624_User_Input_scReadSim_based_Estimate_From_Real_Data_One_Celltype_read_len_" + str(read_len)   + "_" + chr + '_sequence_depth_one_times_GM12878_Haplotype_ASE'

# Phased GIAB variants (GM12878 NGS coverage >= 3), genotype-resolved haplotypes
INPUT_SNP_site = "/u/project/gxxiao/weijian/project/scEditSim/Genomes/GIAB_NA12878_HG001/NISTv3.3.2_GRCh38/by_chrom_cov3/" + "HG001_GRCh38_GIAB_highconf_v.3.3.2_highconf_PGandRTGphasetransfer_"+ chr + "_snps.phased_cov3.vcf.gz"
INPUT_INDEL_site = "/u/project/gxxiao/weijian/project/scEditSim/Genomes/GIAB_NA12878_HG001/NISTv3.3.2_GRCh38/by_chrom_cov3/" + "HG001_GRCh38_GIAB_highconf_v.3.3.2_highconf_PGandRTGphasetransfer_"+ chr + "_indels.phased_cov3.vcf.gz"

# Empirical per-gene allelic ratio pool (R .rds, flat numeric vector in [0,1])
ALLELIC_RATIO_RDS = "/u/project/gxxiao/weijian/project/scEditSim/scripts/allelic_ratio_all_pooled.rds"

# COSMIC somatic mutation catalogs (GRCh38, normalized). Used only when somatic_mutation=True.
COSMIC_CODING_VCF    = "/u/project/gxxiao/weijian/project/scEditSim/Genomes/COSMIC/Cosmic_GenomeScreensMutant_Normal_v104_GRCh38.vcf.gz"
COSMIC_NONCODING_VCF = "/u/project/gxxiao/weijian/project/scEditSim/Genomes/COSMIC/Cosmic_NonCodingVariants_Normal_v104_GRCh38.vcf.gz"

#####################################################################
############################ Main function #########################
#####################################################################


Embed_seq_error = True
somatic_mutation = True   # toggle: COSMIC-guided, cell-type-specific somatic SNVs
if not os.path.isdir(outdirectory):
    os.mkdir(outdirectory)
else:
    # Defensive cleanup of any stale BED files from a previous run. (Not strictly required:
    # scRNA_GenerateBAMCoord_spliced truncates its per-prename *.read.bed/*.read.bed12 with
    # open(...,'w') before writing, and the combined files are rewritten with `cat >`, so a
    # re-run is already safe. This is belt-and-suspenders.)
    import glob
    for _stale in glob.glob(outdirectory + "/*.read.bed") + glob.glob(outdirectory + "/*.read.bed12"):
        os.remove(_stale)

######################## Generate Feature Set ########################
gene_bedfile = outdirectory + "/" + "scReadSim.Gene.bed"
intergene_bedfile = outdirectory + "/" + "scReadSim.InterGene.bed"
Utility.scRNA_CreateFeatureSets(INPUT_bamfile, samtools_directory, bedtools_directory, outdirectory, INPUT_genome_annotation, INPUT_genome_size_file)

######################## Generate UMI Count matrix ########################
UMI_gene_count_mat_filename = "%s.gene.countmatrix" % filename
UMI_intergene_count_mat_filename = "%s.intergene.countmatrix" % filename
# Construct count matrix for foregroud features
Utility.scRNA_bam2countmat_paral(cells_barcode_file=INPUT_cells_barcode_file, bed_file=gene_bedfile, INPUT_bamfile=INPUT_bamfile, outdirectory=outdirectory, count_mat_filename=UMI_gene_count_mat_filename, UMI_modeling=True, UMI_tag = "UB:Z", n_cores=NCORES)
Utility.scRNA_bam2countmat_paral(cells_barcode_file=INPUT_cells_barcode_file, bed_file=intergene_bedfile, INPUT_bamfile=INPUT_bamfile, outdirectory=outdirectory, count_mat_filename=UMI_intergene_count_mat_filename, UMI_modeling=True, UMI_tag = "UB:Z", n_cores=NCORES)

######################## Prepare for RNA editing sites and PHASED SNPs/indels ########################
scIsoSim.prepare_SNP_RED_list_phased(INPUT_SNP_site,INPUT_INDEL_site,INPUT_RED_site, outdirectory,min_edit = 0.01,max_edit=0.99)

######################## Synthetic Matrix Training ########################
GenerateSyntheticCount.scRNA_GenerateSyntheticCount(count_mat_filename=UMI_gene_count_mat_filename, directory=outdirectory, outdirectory=outdirectory ,n_cores = NCORES)
celllabel_file = outdirectory + "/" + UMI_gene_count_mat_filename + ".LouvainClusterResults.txt"
GenerateSyntheticCount.scRNA_GenerateSyntheticCount(count_mat_filename=UMI_intergene_count_mat_filename, directory=outdirectory, outdirectory=outdirectory, celllabel_file=celllabel_file, n_cores = NCORES)

# Specify the names of synthetic count matrices (generated by GenerateSyntheticCount.scRNA_GenerateSyntheticCount)
synthetic_countmat_gene_file = UMI_gene_count_mat_filename + ".scDesign2Simulated.txt"
synthetic_countmat_intergene_file = UMI_intergene_count_mat_filename + ".scDesign2Simulated.txt"
# Specify the base name of bed files containing synthetic reads
OUTPUT_cells_barcode_file = "synthetic_cell_barcode.txt"
gene_read_bedfile_prename = "%s.syntheticBAM.gene" % filename
intergene_read_bedfile_prename = "%s.syntheticBAM.intergene" % filename
BED_filename_combined_pre = "%s.syntheticBAM.combined" % filename
synthetic_cell_label_file = UMI_gene_count_mat_filename + ".scDesign2Simulated.CellTypeLabel.txt"

# Create synthetic read coordinates for genes
# scRNA_GenerateBAMCoord_spliced records CIGAR exon-block structure from each real read and writes:
#   *.read.bed   - BED6 (for downstream RNA editing / allele assignment)
#   *.read.bed12 - BED12 with splice blocks (for splice-aware FASTQ generation)
scIsoSim.scRNA_GenerateBAMCoord_spliced(
        bed_file=gene_bedfile, UMI_count_mat_file=outdirectory + "/" + synthetic_countmat_gene_file, synthetic_cell_label_file=outdirectory + "/" + synthetic_cell_label_file, read_bedfile_prename=gene_read_bedfile_prename, INPUT_bamfile=INPUT_bamfile, outdirectory=outdirectory, OUTPUT_cells_barcode_file=OUTPUT_cells_barcode_file, jitter_size=5, read_len=read_len)
# Create synthetic read coordinates for intergenes
scIsoSim.scRNA_GenerateBAMCoord_spliced(
        bed_file=intergene_bedfile, UMI_count_mat_file=outdirectory + "/" + synthetic_countmat_intergene_file, synthetic_cell_label_file=outdirectory + "/" + synthetic_cell_label_file, read_bedfile_prename=intergene_read_bedfile_prename, INPUT_bamfile=INPUT_bamfile, outdirectory=outdirectory, OUTPUT_cells_barcode_file=OUTPUT_cells_barcode_file, jitter_size=5, read_len=read_len)

# Combine BED6 files. NOTE: vestigial in this splice-aware pipeline -- all downstream steps
# (assign_read_alleles_ASE, prepare_editing_parallel_spliced_ASE, process_indels, somatic) read
# the BED12 below; the combined BED6 is currently unused. Kept for compatibility.
scRNA_GenerateBAM.scRNA_CombineBED(outdirectory=outdirectory, gene_read_bedfile_prename=gene_read_bedfile_prename, intergene_read_bedfile_prename=intergene_read_bedfile_prename, BED_filename_combined_pre=BED_filename_combined_pre)
# Combine BED12 files (used by scRNA_BED2FASTQ_spliced for splice-aware sequence extraction)
scIsoSim.combine_bed12(outdirectory=outdirectory, gene_read_bedfile_prename=gene_read_bedfile_prename, intergene_read_bedfile_prename=intergene_read_bedfile_prename, BED_filename_combined_pre=BED_filename_combined_pre)

synthetic_fastq_prename = BED_filename_combined_pre

# Assign each synthetic read a haplotype via per-gene/intergene allele-specific expression.
# Sampled per-feature ratios are written to feature_allelic_ratios.txt (ground truth).
# This writes <combined>.read_alleles.txt, which the editing and indel steps then only READ
# (prepare_editing_parallel_spliced_ASE does NOT re-randomize it).
scIsoSim.assign_read_alleles_ASE(outdirectory, BED_filename_combined_pre, ALLELIC_RATIO_RDS, per_cell=False, seed=2023)
#
# prepare_editing_parallel_spliced_ASE reads the BED12 file so that:
#   - editing sites falling in introns of spliced reads are excluded
#   - within-read positions are calculated correctly from the exon block structure
#   - read alleles come from assign_read_alleles_ASE (haplotype-preserving), not random
scIsoSim.prepare_editing_parallel_spliced_ASE(outdirectory,BED_filename_combined_pre,RNA_editing_list = outdirectory + "/" + 'RNA_editing_sites_levels.txt',SNPs = outdirectory + "/" + 'SNPs.txt', n_jobs = NCORES,mode = "10x")
# Optional: COSMIC-guided, cell-type-specific somatic SNVs. Injected into
# read_editing_positions.csv NOW (before process_indels) so indel shifting + final application
# are handled by the existing machinery. Each somatic SNV is on allele 0, VAF = its gene's
# sampled allelic ratio, in exactly one cell type, in 1..200 cells, with >=10x pseudobulk coverage.
if somatic_mutation == True:
    scIsoSim.introduce_somatic_mutations_cosmic(outdirectory, BED_filename_combined_pre, referenceGenome_file, COSMIC_CODING_VCF, COSMIC_NONCODING_VCF, chr, n_coding=200, n_noncoding=200, max_cells=200, min_pseudobulk_reads=10)
# Convert combined BED12 into FASTQ using bedtools getfasta -split (splice-aware: exon blocks only)
scIsoSim.scRNA_BED2FASTQ_spliced(bedtools_directory=bedtools_directory, seqtk_directory=seqtk_directory, referenceGenome_file=referenceGenome_file, outdirectory=outdirectory, BED_filename_combined=BED_filename_combined_pre, synthetic_fastq_prename=synthetic_fastq_prename)

scIsoSim.process_indels_in_reads_spliced(outdirectory,synthetic_fastq_prename,referenceGenome_file,n_jobs=NCORES)


if Embed_seq_error == True:
        ################## Generate reads with errors in FASTQs
        # Sequence dictionary for fgbio/htsjdk. htsjdk expects the dict named by REPLACING
        # the fasta extension (GRCh38.p13.genome.dict), NOT appending it (.fa.dict). Build it
        # from referenceGenome_file so the path can't drift from the genome actually used.
        # Created once and atomically (write to a per-process temp, then os.replace) so that:
        #   - re-runs don't fail on Picard's refuse-to-overwrite, and
        #   - parallel per-chromosome jobs sharing this genome dir don't race or read a
        #     half-written dict.
        genome_dict_file = os.path.splitext(referenceGenome_file)[0] + ".dict"
        if not os.path.exists(genome_dict_file):
                tmp_dict_file = "%s.tmp.%d.dict" % (os.path.splitext(referenceGenome_file)[0], os.getpid())
                command_error = f"""
                java -jar {picard_jarfile} CreateSequenceDictionary \
                        -R {referenceGenome_file} \
                        -O {tmp_dict_file}
                """
                ret = os.system(command_error)
                if ret != 0:
                        raise RuntimeError("Picard CreateSequenceDictionary failed (exit %d): %s" % (ret, command_error))
                # Publish atomically; if another job created it meanwhile, keep theirs.
                if not os.path.exists(genome_dict_file):
                        os.replace(tmp_dict_file, genome_dict_file)
                else:
                        os.remove(tmp_dict_file)
        scIsoSim.scRNA_ErrorBase(fgbio_jarfile=fgbio_jarfile, INPUT_bamfile=INPUT_bamfile, referenceGenome_file=referenceGenome_file, outdirectory=outdirectory, synthetic_fastq_prename=BED_filename_combined_pre)
        # Mutate error-included reads. Germline SNP/RNA-editing AND the injected somatic SNVs
        # all live in read_editing_positions.csv (somatic rows added before process_indels, so
        # their positions were indel-shifted) and are applied together here.
        scIsoSim.mutate_fastq_sequences_parallel(read_editing_pos = outdirectory + "/" + "read_editing_positions.csv", input_read1 = outdirectory + "/" + BED_filename_combined_pre + ".ErrorIncluded.read1.bed2fa.sorted.fq", input_read2 = outdirectory + "/" + BED_filename_combined_pre + ".ErrorIncluded.read2.bed2fa.sorted.fq", output_read1 = outdirectory + "/" + BED_filename_combined_pre + ".EditingIncluded_ErrorIncluded.read1.bed2fa.sorted.fq", output_read2 = outdirectory + "/" + BED_filename_combined_pre + ".EditingIncluded_ErrorIncluded.read2.bed2fa.sorted.fq", synthetic_editing = outdirectory + "/" +'Synthetic_RNA_editing_events.csv',outdirectory = outdirectory, ErrorInclude = Embed_seq_error)
        # NOTE (20260624): quality model is now "error-only". ErroneousRead() already stamps a
        # lowered, position-specific Q on each injected sequencing-error base and leaves every
        # other base at the bed2fa default F (Q37); the editing step preserves that quality line.
        # So we DELIBERATELY no longer overwrite the whole quality string per position -- that
        # call collapsed every read to an identical low-Q string (Q24-26 body, Q60 ends), which
        # pushed ~21% of bases under REDItools -bq25 and tripped SPRINT's cutoff detector.
        # Final read2 quality is therefore: F everywhere EXCEPT sequencing-error bases.
        # scIsoSim.replace_fastq_quality_with_sequence(fastq_file=outdirectory + "/" + BED_filename_combined_pre + ".EditingIncluded_ErrorIncluded.read2.bed2fa.sorted.fq", error_rate_file = outdirectory + "/" + 'Real.error_rate_by_read_position.txt')
        scIsoSim.output_vcf_from_csvs(snv_csv = outdirectory + "/ErrorIncluded_Synthetic_RNA_editing_events.csv", indel_csv = outdirectory + "/" + BED_filename_combined_pre + ".introduced_indel.csv", out_vcf = outdirectory + "/variants.vcf")
# [Error-free simulation removed] Only the error-included path (inside the Embed_seq_error
# block above) is produced:
#   .EditingIncluded_ErrorIncluded.read{1,2}.bed2fa.sorted.fq   (synthetic reads)
#   ErrorIncluded_Synthetic_RNA_editing_events.csv             (ground-truth edits)
#   ErrorIncluded_{edited,unedited}_events_by_cell_matrix.csv  (per-cell matrices; the
#       error-included mutate call uses output_per_cell_editing_matrix=True by default)
# NOTE: with this removed, the script requires Embed_seq_error=True to produce any output.
