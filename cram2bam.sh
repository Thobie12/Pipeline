#!/bin/bash
#SBATCH --job-name=cram2bam
#SBATCH --output=/cluster/home/t922316uhn/PLO/cram2bam.out
#SBATCH --error=/cluster/home/t922316uhn/PLO/cram2bam.err
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=3-00:00:00

# Load samtools (if needed on your cluster)
module load samtools

# Input CRAM and output paths
CRAM="/cluster/projects/pughlab/myeloma/external_data/ultimagen-oicr/Crams/OICRM4CA-07-01-P.cram"
REF="/cluster/tools/data/genomes/human/GRCh38/iGenomes/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta"
BAM="/cluster/projects/pughlab/myeloma/projects/MM_cell_drugs/WGS_Pipeline/Bam/OICRM4CA-07-01-P.bam"

# Convert CRAM to BAM
samtools view -@ 8 -m 2G -T $REF -b -o "$BAM" "$CRAM"

# Index BAM
samtools index -@ 8 -m 2G "$BAM"

echo "Conversion and indexing complete for $CRAM"