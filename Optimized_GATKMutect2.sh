#!/bin/bash
#SBATCH --job-name=mutect2_chr
#SBATCH --output=/cluster/home/t922316uhn/PLO/GATK/mutect2_chr_%A_%a.out
#SBATCH --error=/cluster/home/t922316uhn/PLO/GATK/mutect2_chr_%A_%a.err
#SBATCH --time=7-00:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --partition=pughlab
#SBATCH --array=1-24%4   # 22 autosomes + X + Y (adjust if needed)

module load gatk
module load samtools

# List of chromosomes (modify as needed)
CHRS=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY)

CHR=${CHRS[$SLURM_ARRAY_TASK_ID-1]}

gatk Mutect2 \
  -R /cluster/tools/data/genomes/human/GRCh38/iGenomes/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta \
  -I /cluster/projects/pughlab/myeloma/projects/MM_cell_drugs/WGS_Pipeline/Bam/OICRM4CA-07-01-P.bam \
  -tumor OICRM4CA-07-01-P \
  -L $CHR \
  --panel-of-normals /cluster/tools/data/genomes/human/GRCh38/iGenomes/Annotation/GATKBundle/1000g_pon.hg38.vcf.gz \
  --germline-resource /cluster/tools/data/genomes/human/GRCh38/iGenomes/Annotation/GATKBundle/af-only-gnomad.hg38.vcf.gz  \
  --f1r2-tar-gz /cluster/projects/pughlab/myeloma/projects/MM_cell_drugs/WGS_Pipeline/GATK/${CHR}.f1r2.tar.gz \
  --native-pair-hmm-threads 8 \
  -O /cluster/projects/pughlab/myeloma/projects/MM_cell_drugs/WGS_Pipeline/GATK/OICRM4CA-07-01-P_${CHR}.vcf.gz




  merge upon completion
  gatk MergeVcfs \
  -I /output/path/mutect2_chr1.vcf.gz \
  -I /output/path/mutect2_chr2.vcf.gz \
  ... (all chromosomes) ... \
  -O /output/path/mutect2_merged.vcf.gz