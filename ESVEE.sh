#!bin/bash
#ESVEE pre and post

#Pre GRIDSS code
java -cp esvee.jar com.hartwig.hmftools.esvee.prep.PrepApplication \
  -sample OICRM4CA-07-01-P \
  -bam_file /cluster/projects/pughlab/myeloma/projects/MM_cell_drugs/WGS_Pipeline/Bam/OICRM4CA-07-01-P.bam \
  -ref_genome /cluster/tools/data/genomes/human/GRCh38/iGenomes/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta \
  -ref_genome_version 38 \
  -output_dir /cluster/projects/pughlab/myeloma/projects/MM_cell_drugs/WGS_Pipeline/ESVEE/ 
  -known_fusion_bed #include if I can find it
  -blacklist_bed #include if I can find it

#GRIDSS Code
module load R
module load bwa
BAM="??"
REF="/cluster/projects/pughlab/myeloma/projects/MM_cell_drugs/WGS_Pipeline/ref/WholeGenomeFasta/Homo_sapiens_assembly38.fasta"
OUT_VCF="/cluster/projects/pughlab/myeloma/projects/MM_cell_drugs/WGS_Pipeline/VCF/OICRM4CA-07-01-P.sv.vcf"
OUT_ASSEMBLY="OICRM4CA-07-01-P.gridss.assembly.bam"
THREADS=8
echo "Starting GRIDSS run on ${CRAM}"
date

/cluster/home/t922316uhn/gridss/gridss \
  --jar /cluster/home/t922316uhn/gridss/gridss-2.13.2-gridss-jar-with-dependencies.jar \
  --reference "$REF" \
  --output "$OUT_VCF" \
  --assembly "$OUT_ASSEMBLY" \
  --threads "$THREADS" \
  "$BAM"

echo "GRIDSS finished"
date

#Post GRIDSS Code
java -cp esvee.jar com.hartwig.hmftools.esvee.depth.DepthAnnotator \
  -input_vcf ${gridss_vcf} \
  -output_vcf ${final_vcf} \
  -samples OICRM4CA-07-01-P \
  -bam_files "${tumor_bam}" \
  -ref_genome /cluster/tools/data/genomes/human/GRCh38/iGenomes/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta \
  -ref_genome_version 38 \
  -threads 8 \