ssh t922316uhn@h4huhnlogin1.uhnresearch.ca # ssh into server
git clone https://github.com/pughlab/pipeline-suite.git #this is for cloning, when cloning repo, do it without salloc session initiated
salloc -c 1 -p short -t 1:0:0 --mem 1G # request a short job if needed, remove -p if not a short session
module load perl # load perl module or any other module you need
#RNASEQ PIPELINE
perl /cluster/home/t922316uhn/pipeline-suite/pughlab_rnaseq_pipeline.pl \ 
-t /cluster/home/t922316uhn/pipeline-suite/configs/rna_pipeline_config.yaml \ 
-d /cluster/home/t922316uhn/pipeline-suite/configs/rna_fastq_config.yaml \
--summarize \
-c slurm \
--remove \
--dry-run # --dry-run is optional, it will not run the pipeline, just show what will be done, remove if you want to run the pipeline


perl /cluster/home/t922316uhn/pipeline-suite/pughlab_rnaseq_pipeline.pl \
-t /cluster/home/t922316uhn/pipeline-suite/configs/rna_pipeline_config.yaml \ 
-d /cluster/projects/pughlab/myeloma/projects/MM_cell_drugs/Combined/rna_fastq_config.yaml \
--summarize \
-c slurm \
--remove \
--dry-run # --dry-run is optional, it will not run the pipeline, just show what will be done, remove if you want to run the pipeline
