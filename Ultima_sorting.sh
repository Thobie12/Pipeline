#!/bin/bash
Ultima_sorting 

samtools sort -T /cluster/projects/pughlab/myeloma/projects/MM_cell_drugs/WGS_Pipeline/Dedup/temp \
  -o /cluster/projects/pughlab/myeloma/projects/MM_cell_drugs/WGS_Pipeline/Dedup/OICRM4CA-07-01-P.sorted.cram \
  # o = output name
  -O CRAM \
  -@ 8 \  
  # dont change threads and mem
  -m 2G \
  --reference /cluster/tools/data/genomes/human/GRCh38/iGenomes/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta \
  /cluster/projects/pughlab/myeloma/external_data/ultimagen-oicr/Crams/OICRM4CA-07-01-P.cram 
  # input is the last entry

---
#Automated version for mutliple samples
!/bin/bash
# Automated version for multiple samples
INPUT="/cluster/projects/pughlab/myeloma/external_data/ultimagen-oicr/Crams"
OUTPUT="/cluster/projects/pughlab/myeloma/projects/MM_cell_drugs/WGS_Pipeline/Dedup"
REF="/cluster/tools/data/genomes/human/GRCh38/iGenomes/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta"
TEMP="/cluster/projects/pughlab/myeloma/projects/MM_cell_drugs/WGS_Pipeline/Dedup/temp"

for cram in ${INPUT}/OICRM4*.cram; do
    sample=$(basename "$cram" .cram)
    output="${OUTPUT}/${sample}.sorted.cram"
    
    echo "Sorting $sample"

    samtools sort -T "${TEMP}" \
        -o "${output}" \
        -O CRAM \
        -@ 8 \
        -m 2G \
        --reference "${REF}" \
        "$cram"
done