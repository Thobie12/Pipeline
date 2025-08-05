#!/bin/bash
#SBATCH --job-name=gridss_sv_calling
#SBATCH --output=/cluster/home/t922316uhn/PLO/GRIDSS/gridss_CA-08%j.out
#SBATCH --error=/cluster/home/t922316uhn/PLO/GRIDSS/gridss_CA-08%j.err
#SBATCH --time=96:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --partition=pughlab


module load R
module load bwa

# Variables
CRAM="/cluster/projects/pughlab/myeloma/external_data/ultimagen-oicr/Crams/OICRM4CA-08-R-P.cram"
REF="/cluster/projects/pughlab/myeloma/projects/MM_cell_drugs/WGS_Pipeline/ref/WholeGenomeFasta/Homo_sapiens_assembly38.fasta"
OUT_VCF="/cluster/projects/pughlab/myeloma/projects/MM_cell_drugs/WGS_Pipeline/GRIDSS/OICRM4CA-08-R-P/OICRM4CA-08-R-P.sv.vcf"
OUT_ASSEMBLY="OICRM4CA-07-01-P.gridss.assembly.bam"
THREADS=8

echo "Starting GRIDSS run on ${CRAM}"
date

/cluster/home/t922316uhn/gridss/gridss \
  --jar /cluster/home/t922316uhn/gridss/gridss-2.13.2-gridss-jar-with-dependencies.jar \
  --reference "$REF" \
  --reference "$REF" \
  --output "$OUT_VCF" \
  --assembly "$OUT_ASSEMBLY" \
  --threads "$THREADS" \
  "$CRAM"

echo "GRIDSS finished"
date

---
