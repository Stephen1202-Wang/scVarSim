#!/bin/bash
#$ -cwd
#$ -j y
#$ -o /u/home/w/weijian/project-gxxiao/project/scEditSim/log/20260618_GM12878_chr1_Estimate_from_real.log
#  Resources requested
#  PLEASE CHANGE THE RESOURCES REQUESTED AS NEEDED:
#$ -l h_data=60G,h_rt=10:00:00,highp
#  PLEASE CHANGE THE NUMBER OF CORES REQUESTED AS NEEDED:
#$ -pe shared 1

#
# Set up job environment:
#
. /u/local/Modules/default/init/modules.sh
module load mamba
eval "$(conda shell.bash hook)"

conda activate scRET        

chr='chr1'
outdirectory=/u/project/gxxiao/weijian/project/scEditSim/result/Real_data_GRCh38/GM12878/${chr}

mkdir -p $outdirectory

reference_genome=~/project-gxxiao/project/scEditSim/Genomes/GRCh38/GRCh38.p13.genome.fa
VCF_FILE=~/project-gxxiao/project/scEditSim/Genomes/00-common_all_${chr}_chr_attached.vcf.gz
input_bam=/u/project/gxxiao/weijian/project/scEditSim/data/NGS_GM12878/GRCh38/STARsolo_alignment/GENCODE_Genome/split_by_chromosomes/GM12878_${chr}.CBattached.filtered.bam

#############
############# Detecting SNPs
#############
bcftools mpileup -Ou \
    -f ${reference_genome} \
    -Q 20 \
    -q 20 \
    --max-depth 1000000 \
    ${input_bam} | \
bcftools call -mv -Oz -o ${outdirectory}/variants_bcftools.vcf.gz

# Step 2: Index VCF
bcftools index ${outdirectory}/variants_bcftools.vcf.gz

# Step 3: Add AF (allele frequency) annotation
bcftools +fill-tags ${outdirectory}/variants_bcftools.vcf.gz \
    -- -t AF | \
bcftools view -Oz -o ${outdirectory}/variants_bcftools_with_AF.vcf.gz

# Step 4: Index the AF-annotated VCF
bcftools index ${outdirectory}/variants_bcftools_with_AF.vcf.gz

# Step 5: Intersect (only on position, ignoring alleles)
bcftools isec -p ${outdirectory} -n =2 -c none \
    ${outdirectory}/variants_bcftools_with_AF.vcf.gz \
    ${VCF_FILE}

bcftools view -v indels ${outdirectory}/variants_bcftools_with_AF.vcf.gz | \
    bcftools norm -f ${reference_genome} -m -any | \
    bcftools norm -f ${reference_genome} | \
    awk 'BEGIN{OFS="\t"}
         /^#/{print; next}
         {
             key=$1"\t"$2
             qual=($6=="."?-1:$6+0)
             if(!(key in best) || qual > best[key]){best[key]=qual; line[key]=$0}
         }
         END{for(k in line) print line[k]}' | \
    bcftools sort -o ${outdirectory}/indels_no_intersection.vcf
bcftools view -v snps ${outdirectory}/0000.vcf > ${outdirectory}/snps_intersection.vcf

############
############ Detecting RNA Editing Sites
############
samtools mpileup -f ${reference_genome} -Q 20 -q 20 -d 1000000 ${input_bam} | \
awk -v OFS='\t' '
BEGIN {
    print "CHROM","POS","REF","ALT","REF_COUNT","ALT_COUNT"
} 
{
    if($4 > 0) {
        chrom = $1
        pos = $2
        ref = toupper($3)
        depth = $4
        bases = toupper($5)
        
        # Initialize counters
        ref_count = 0
        delete alt_counts
        
        # Parse pileup bases
        i = 1
        while(i <= length(bases)) {
            base = substr(bases, i, 1)
            
            if(base == "." || base == ",") {
                # Reference base
                ref_count++
            }
            else if(base == "^") {
                # Skip mapping quality after ^
                i++
            }
            else if(base == "$") {
                # End of read marker, skip
            }
            else if(base == "*") {
                # Deletion, skip for now
            }
            else if(base == "+" || base == "-") {
                # Insertion/deletion, skip the indel sequence
                i++
                indel_len = ""
                while(i <= length(bases) && substr(bases,i,1) ~ /[0-9]/) {
                    indel_len = indel_len substr(bases,i,1)
                    i++
                }
                if(indel_len != "") {
                    i += int(indel_len) - 1
                }
            }
            else if(base ~ /[ATCGN]/) {
                # Alternative base
                if(base != ref && base != "N") {
                    alt_counts[base]++
                }
            }
            i++
        }
        
        # Output results for each alternative base found
        for(alt in alt_counts) {
            if(alt_counts[alt] > 0) {
                print chrom, pos, ref, alt, ref_count, alt_counts[alt]
            }
        }
    }
}' > ${outdirectory}/mismatch_counts.tsv