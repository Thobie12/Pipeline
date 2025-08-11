#!/bin/bash
#SBATCH --job-name=mosdepth
#SBATCH --output=/cluster/home/t922316uhn/PLO/mosdepth%j.out
#SBATCH --error=/cluster/home/t922316uhn/PLO/mosdepth%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --partition=pughlab


/cluster/home/t922316uhn/mosdepth//mosdepth --fast-mode --by 1000000 --threads 8 \
  --fasta /cluster/tools/data/genomes/human/GRCh38/iGenomes/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta \
  /cluster/projects/pughlab/myeloma/projects/MM_cell_drugs/WGS_Pipeline/depth/OICRM4CA-07-01-P \
  /cluster/projects/pughlab/myeloma/external_data/ultimagen-oicr/Crams/OICRM4CA-07-01-P.cram


 zcat OICRM4CA-07-01-P.regions.bed.gz | awk '
BEGIN { OFS="\t" }
$1 ~ /^chr([1-9]|1[0-9]|2[0-2]|X)$/ {
    chr = $1; start = $2; end = $3; cov = $4;
    if (NR == 1 || chr != prev_chr || start != prev_end) {
        print "fixedStep chrom=" chr " start=" (start+1) " step=" (end - start) " span=" (end - start);
    }
    print cov;
    prev_chr = chr; prev_end = end;
}
' > OICRM4CA-07-01-P.mosdepth2.wig



snakemake -s /cluster/home/t922316uhn/ichorCNA/scripts/snakemake/ichorCNA.snakefile -p \
    --cluster "sbatch --mem=16G --cpus-per-task=1 --partition=pughlab --time=12:0:0" \
    --jobs 100



module load R && Rscript /cluster/home/t922316uhn/ichorCNA/scripts/runIchorCNA.R \
  --id OICRM4CA-07-01-P \
  --libdir None \
  --WIG /cluster/projects/pughlab/myeloma/projects/MM_cell_drugs/WGS_Pipeline/ichorCNA/readDepth/OICRM4CA-07-01-P.mosdepth2.wig \
  --gcWig /cluster/home/t922316uhn/ichorCNA/inst/extdata/gc_hg38_1000kb.wig \
  --mapWig /cluster/home/t922316uhn/ichorCNA/inst/extdata/map_hg38_1000kb.wig \
  --normalPanel /cluster/home/t922316uhn/ichorCNA/inst/extdata/HD_ULP_PoN_hg38_1Mb_median_normAutosome_median.rds \
  --ploidy "c(2,3)" \
  --normal "c(0.5,0.6,0.7,0.8,0.9,0.95)" \
  --maxCN 5 \
  --includeHOMD False \
  --chrs "paste0('chr', c(1:22, 'X', 'Y'))" \
  --chrTrain "paste0('chr', c(1:22))" \
  --genomeStyle UCSC \
  --genomeBuild hg38 \
  --estimateNormal True \
  --estimatePloidy True \
  --estimateScPrevalence True \
  --scStates "c(1,3)" \
  --centromere /cluster/home/t922316uhn/ichorCNA/inst/extdata/GRCh38.GCA_000001405.2_centromere_acen.txt \
  --exons.bed None \
  --txnE 0.9999 \
  --txnStrength 10000 \
  --minMapScore 0.75 \
  --fracReadsInChrYForMale 0.002 \
  --maxFracGenomeSubclone 0.5 \
  --maxFracCNASubclone 0.7 \
  --plotFileType pdf \
  --plotYLim "c(-1,1)" \
  --outDir /cluster/projects/pughlab/myeloma/projects/MM_cell_drugs/WGS_Pipeline/ichorCNA/test/ \
  > /cluster/projects/pughlab/myeloma/projects/MM_cell_drugs/WGS_Pipeline/ichorCNA/logs/ichorCNA/OICRM4CA-07-01-P.test.log 2>&1