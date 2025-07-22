#Parascopy Setup + Code
mkdir Parascopy
module load python3
conda create --name paras_env parascopy
conda activate paras_env
conda config --add channels bioconda
conda config --add channels conda-forge
conda install -c bioconda parascopy
wget -c https://dl.dropboxusercontent.com/s/okzeedb6gze6zzs/homology_table_hg38.tar #download homology file
tar xf homology_table_hg38.tar
# Calculate background read depth. takes about 2hrs
parascopy depth -I /cluster/projects/pughlab/myeloma/projects/MM_cell_drugs/WGS_Pipeline/Parascopy/input-list.txt -g hg38 -f /cluster/tools/data/genomes/human/hg38/iGenomes/Sequence/WholeGenomeFasta/genome.fa -o /cluster/projects/pughlab/myeloma/projects/MM_cell_drugs/WGS_Pipeline/Parascopy/CA-07/depth
# input-list.txt contains /cluster/projects/pughlab/myeloma/external_data/ultimagen-oicr/Crams/OICRM4CA-07-01-P.cram
# Estimate agCN and psCN for multiple samples, I removed region
parascopy cn -I /cluster/projects/pughlab/myeloma/projects/MM_cell_drugs/WGS_Pipeline/Parascopy/input-list.txt -t /cluster/home/t922316uhn/parascopy/homology/homology_table/hg38.bed.gz -r /cluster/projects/pughlab/myeloma/projects/MM_cell_drugs/WGS_Pipeline/Parascopy/mm_targets.bed -f /cluster/tools/data/genomes/human/hg38/iGenomes/Sequence/WholeGenomeFasta/genome.fa  -d /cluster/projects/pughlab/myeloma/projects/MM_cell_drugs/WGS_Pipeline/Parascopy/CA-07/depth -o /cluster/projects/pughlab/myeloma/projects/MM_cell_drugs/WGS_Pipeline/Parascopy/CA-07/analysis1
# variant calling
parascopy call -p /cluster/projects/pughlab/myeloma/projects/MM_cell_drugs/WGS_Pipeline/Parascopy/CA-07/analysis1 -f genome.fa -t /cluster/home/t922316uhn/parascopy/homology/homology_table/hg38.bed.gz -f /cluster/tools/data/genomes/human/hg38/iGenomes/Sequence/WholeGenomeFasta/genome.fa

