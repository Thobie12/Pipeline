#!/bin/bash
#SBATCH --job-name=LearnReadOrientationModel
#SBATCH --output=/cluster/home/t922316uhn/PLO/GATK/filtercalls_CA-07_%j.out
#SBATCH --error=/cluster/home/t922316uhn/PLO/GRIDSS/filtercalls_CA-07_%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --partition=all

#sig profiler to do
module load python3

pip install SigProfilerAssignment
python
# from SigProfilerMatrixGenerator import install as genInstall
# genInstall.install('GRCh38')

#above is for installing the matrix generator for GRCh38
# Now we can use the SigProfilerAssignment package to analyze the VCF file

from SigProfilerAssignment import Analyzer as Analyze

# Path to your filtered Mutect2 VCF file
samples = "/cluster/projects/pughlab/myeloma/projects/MM_cell_drugs/WGS_Pipeline/GATK/OICRM4CA-07-01-P/filtered_final.vcf.gz"

# Output directory where results will be saved
output = "/cluster/projects/pughlab/myeloma/projects/MM_cell_drugs/WGS_Pipeline/GATK/OICRM4CA-07-01-P"
# Yes, that's likely the reason. The SigProfilerAssignment package uses SigProfilerPlotting to generate plots.
# You should install it before running the analysis:


Analyze.cosmic_fit(
    samples=samples,
    output=output,
    input_type="vcf",            # since you have VCF input
    context_type="96",           # standard SBS96 context
    collapse_to_SBS96=True,
    cosmic_version="3.4",        # latest COSMIC v3.4
    exome=False,                 # set True if WES data
    genome_build="GRCh38",       # Use GRCh38 as requested
    signature_database=None,     # default COSMIC signatures
    exclude_signature_subgroups=None,
    export_probabilities=True,
    export_probabilities_per_mutation=True,
    make_plots=True,             # Set this to True to get plots
    sample_reconstruction_plots=True,
    verbose=True                 # To get detailed logging
)
