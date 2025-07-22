#!/bin/bash
#SBATCH --job-name=mutect2_chr
#SBATCH --output=mutect2_chr_%A_%a.out
#SBATCH --error=mutect2_chr_%A_%a.err
#SBATCH --time=2-00:00:00
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
  -I /path/to/your.bam \
  -tumor YOUR_SAMPLE_NAME \
  -L $CHR \
  --panel-of-normals /path/to/1000g_pon.hg38.vcf.gz \
  --germline-resource /path/to/af-only-gnomad.hg38.vcf.gz \
  --f1r2-tar-gz ${CHR}.f1r2.tar.gz \
  --native-pair-hmm-threads 8 \
  --max-reads-per-alignment-start 20 \
  -O /output/path/mutect2_${CHR}.vcf.gz




  merge upon completion
  gatk MergeVcfs \
  -I /output/path/mutect2_chr1.vcf.gz \
  -I /output/path/mutect2_chr2.vcf.gz \
  ... (all chromosomes) ... \
  -O /output/path/mutect2_merged.vcf.gz