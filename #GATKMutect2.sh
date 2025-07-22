#GATKMutect2 Steps
#!/bin/bash
#SBATCH --job-name=GATK_Mutect2
#SBATCH --output=/cluster/home/t922316uhn/PLO/GATK_%j.out
#SBATCH --error=/cluster/home/t922316uhn/PLO/GATK_%j.err
#SBATCH --time=5-00:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --partition=pughlab

module load gatk

gatk Mutect2 \
  -R /cluster/tools/data/genomes/human/GRCh38/iGenomes/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta \
  -I /cluster/projects/pughlab/myeloma/external_data/ultimagen-oicr/Crams/OICRM4CA-07-01-P.cram \
  -tumor OICRM4CA-07-01-P \
  --panel-of-normals /cluster/tools/data/genomes/human/GRCh38/iGenomes/Annotation/GATKBundle/1000g_pon.hg38.vcf.gz  \
  --germline-resource /cluster/tools/data/genomes/human/GRCh38/iGenomes/Annotation/GATKBundle/af-only-gnomad.hg38.vcf.gz  \
  --f1r2-tar-gz OICRM4CA-07-01-P.f1r2.tar.gz \
  --native-pair-hmm-threads 4 \
  -O /cluster/projects/pughlab/myeloma/projects/MM_cell_drugs/WGS_Pipeline/GATK/OICRM4CA-07-01-P.vcf.gz

