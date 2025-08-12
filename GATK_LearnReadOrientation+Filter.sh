#!/bin/bash
#SBATCH --job-name=LearnReadOrientationModel
#SBATCH --output=/cluster/home/t922316uhn/PLO/GATK/Learn_CA-07_%j.out
#SBATCH --error=/cluster/home/t922316uhn/PLO/GRIDSS/Learn_CA-07_%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=256G
#SBATCH --partition=all

module load gatk

# Directory containing all chr*.f1r2.tar.gz files
F1R2_DIR="/cluster/projects/pughlab/myeloma/projects/MM_cell_drugs/WGS_Pipeline/GATK/OICRM4CA-07-01-P"

# Build input args for each file
INPUTS=""
for f in "${F1R2_DIR}"/chr*.f1r2.tar.gz; do
  INPUTS+="-I $f "
done

# Run LearnReadOrientationModel with all inputs
gatk LearnReadOrientationModel \
  $INPUTS \
  -O /cluster/projects/pughlab/myeloma/projects/MM_cell_drugs/WGS_Pipeline/GATK/OICRM4CA-07-01-P/read-orientation-model.tar.gz

echo "Done!"



------
#for filtering calls
#!/bin/bash
#SBATCH --job-name=LearnReadOrientationModel
#SBATCH --output=/cluster/home/t922316uhn/PLO/GATK/filtercalls_CA-07_%j.out
#SBATCH --error=/cluster/home/t922316uhn/PLO/GRIDSS/filtercalls_CA-07_%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --partition=all

module load gatk

 gatk FilterMutectCalls \
   -R /cluster/projects/pughlab/myeloma/projects/MM_cell_drugs/WGS_Pipeline/ref/WholeGenomeFasta/Homo_sapiens_assembly38.fasta \
   -V /cluster/projects/pughlab/myeloma/projects/MM_cell_drugs/WGS_Pipeline/GATK/OICRM4CA-07-01-P/OICRM4CA-07-01-P_all.vcf.gz \
   -ob-priors /cluster/projects/pughlab/myeloma/projects/MM_cell_drugs/WGS_Pipeline/GATK/OICRM4CA-07-01-P/read-orientation-model.tar.gz \
   -O  /cluster/projects/pughlab/myeloma/projects/MM_cell_drugs/WGS_Pipeline/GATK/OICRM4CA-07-01-P/filtered.vcf.gz

   echo "Done!"



-----------
#filtering for just pass
bcftools view -f PASS \
  -r chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY \
  /cluster/projects/pughlab/myeloma/projects/MM_cell_drugs/WGS_Pipeline/GATK/OICRM4CA-07-01-P/filtered.vcf.gz \
  -Oz -o /cluster/projects/pughlab/myeloma/projects/MM_cell_drugs/WGS_Pipeline/GATK/OICRM4CA-07-01-P/filtered_final.vcf.gz