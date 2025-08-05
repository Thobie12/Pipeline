#Control FREEC using docker
module load singularity
cd /cluster/home/t922316uhn/UGBIO
singularity pull ultimagenomics/ugbio_freec:1.5.5
singularity pull docker://ultimagenomics/ugbio_cnv:1.5.5
 singularity exec -B /cluster:/cluster /cluster/home/t922316uhn/UGBIO/ugbio_freec_1.5.5.sif /bin/bash

# Build the mpileup file

#Create mpileup file for tumor

#!/bin/bash
#SBATCH --job-name=mpileup_freec
#SBATCH --output=mpileup_freec_%j.log
#SBATCH --error=mpileup_freec_%j.err
#SBATCH --time=72:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G

# Load Singularity if your cluster requires module loading
module load singularity

# Paths
SIF_PATH=/cluster/home/t922316uhn/UGBIO/ugbio_freec_1.5.5.sif
REF_FASTA=/cluster/tools/data/genomes/human/GRCh38/iGenomes/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta
INPUT_CRAM=/cluster/projects/pughlab/myeloma/external_data/ultimagen-oicr/Crams/OICRM4CA-07-01-P.cram
OUTPUT_PILEUP=/cluster/projects/pughlab/myeloma/projects/MM_cell_drugs/WGS_Pipeline/FREEC/OICRM4CA-07-01-P/OICRM4CA-07-01-P_minipileup.pileup

# Run mpileup inside Singularity container
singularity exec -B /cluster:/cluster $SIF_PATH \
  samtools mpileup \
    -f $REF_FASTA \
    -d 8000 \
    -Q 0 \
    -q 1 \
    $INPUT_CRAM \
    > $OUTPUT_PILEUP

# Collect coverage for tumor and normal separately from inside freec
#!/bin/bash
#SBATCH --job-name=depth_freec
#SBATCH --output=depth_freec_%j.log
#SBATCH --error=depth_freec_%j.err
#SBATCH --time=72:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G

# Load Singularity if required
module load singularity

# Paths
SIF_PATH=/cluster/home/t922316uhn/UGBIO/ugbio_freec_1.5.5.sif
REF_FASTA=/cluster/tools/data/genomes/human/GRCh38/iGenomes/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta
INPUT_CRAM=/cluster/projects/pughlab/myeloma/external_data/ultimagen-oicr/Crams/OICRM4CA-07-01-P.cram
OUTPUT_BEDGRAPH=/cluster/projects/pughlab/myeloma/projects/MM_cell_drugs/WGS_Pipeline/FREEC/OICRM4CA-07-01-P/OICRM4CA-07-01-P.bedgraph

# Run samtools depth inside Singularity
singularity exec -B /cluster:/cluster $SIF_PATH \
  samtools depth \
    -J \
    -Q 1 \
    --reference $REF_FASTA \
    $INPUT_CRAM | \
  awk '{print $1"\t"($2-1)"\t"$2"\t"$3}' > $OUTPUT_BEDGRAPH



# Bedgraph to CPN
#unzip input bedgraph file if needed:
#!/bin/bash
#SBATCH --job-name=bedtools_map_freec
#SBATCH --output=bedtools_map_freec_%j.log
#SBATCH --error=bedtools_map_freec_%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G

# Load Singularity if needed
module load singularity

# Paths
SIF_PATH=/cluster/home/t922316uhn/UGBIO/ugbio_freec_1.5.5.sif
FAI_PATH=/cluster/tools/data/genomes/human/GRCh38/iGenomes/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta.fai
BINS_BED=/cluster/home/t922316uhn/bed/Homo_sapiens_assembly38.w1000.bed
BEDGRAPH_GZ=/cluster/projects/pughlab/myeloma/projects/MM_cell_drugs/WGS_Pipeline/FREEC/OICRM4CA-07-01-P/OICRM4CA-07-01-P.bedgraph.gz
BEDGRAPH=/cluster/projects/pughlab/myeloma/projects/MM_cell_drugs/WGS_Pipeline/FREEC/OICRM4CA-07-01-P/OICRM4CA-07-01-P.bedgraph
OUTPUT_CPN=/cluster/projects/pughlab/myeloma/projects/MM_cell_drugs/WGS_Pipeline/FREEC/OICRM4CA-07-01-P/OICRM4CA-07-01-P.cpn

# Decompress if .gz exists
if [[ $BEDGRAPH_GZ =~ \.gz$ && -f "$BEDGRAPH_GZ" ]]; then
    gzip -d -c "$BEDGRAPH_GZ" > "$BEDGRAPH"
fi

# Run bedtools map inside Singularity
singularity exec -B /cluster:/cluster $SIF_PATH \
  bedtools map -g "$FAI_PATH" \
    -a "$BINS_BED" \
    -b "$BEDGRAPH" \
    -c 4 -o mean | \
  awk '{if($4=="."){print $1"\t"$2"\t"0}else{print $1"\t"$2"\t"$4}}' | \
  grep -v "chrY" | \
  sed 's/^chr//' > "$OUTPUT_CPN"
