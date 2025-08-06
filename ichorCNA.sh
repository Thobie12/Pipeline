git clone https://github.com/broadinstitute/ichorCNA
# cd into snake makefile
# tweak the snakemake file as needed and combine config files into 1 (combined hg38 and slurm into one config and changed the config file directory from the snakemake file)
then submit the code to run 
snakemake -s ichorCNA.snakefile -np --cluster "sbatch --mem={resources.mem}G" --jobs 10 # dry run , remove n for real run
snakemake -s ichorCNA.snakefile -p --cluster "sbatch --mem={resources.mem}G" --jobs 10
# if file is locked, run code below to unlock
snakemake -s ichorCNA.snakefile --unlock