FREEC Setup + Code
mkdir FREEC
cd FREEC
wget https://github.com/BoevaLab/FREEC/archive/refs/tags/v11.6b.tar.gz
tar -zxvf  v11.6b.tar.gz
cd FREEC-11.6b
cd src
make
mkdir mappability
cd mappability
wget http://xfer.curie.fr/get/vyIi4w8EONl/out100m2_hg38.zip
wget --content-disposition "https://cloud.inf.ethz.ch/s/idTaGpZdnS9To5c/download"
#unzip zip files after
#can download test file after for trial run and unzip
/cluster/home/t922316uhn/FREEC/FREEC-11.6b/src/freec    -conf config_chr19.txt #running test file
cat makeGraph_Chromosome.R | R --slave --args 19 2 chr_19.noDup0.pileup.gz_ratio.txt chr_19.noDup0.pileup.gz_BAF.txt 


#for batch submission
#!/bin/bash

#SBATCH --job-name=FREEC
#SBATCH --output=/cluster/home/t922316uhn/PLO_%j.out
#SBATCH --error=/cluster/home/t922316uhn/PLO_%j.err
#SBATCH --time=7-00:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --partition=pughlab

module load samtools
module load R
module load sambamba

/cluster/home/t922316uhn/FREEC/FREEC-11.6b/src/freec    -conf /cluster/home/t922316uhn/FREEC/Configs/CA-07.txt #running test file
