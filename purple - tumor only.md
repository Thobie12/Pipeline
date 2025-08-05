# Sample code for purple - tumor only
java -jar purple.jar \
   -tumor OICRM4CA-07-01-P \ 
   # use same tumor name in amber and cobalt
   -amber /cluster/projects/pughlab/myeloma/projects/MM_cell_drugs/WGS_Pipeline/AMBER \
   -cobalt /path/cobalt/OICRM4CA-07-01-P \
   -gc_profile /path/GC_profile.1000bp.37.cnp \
   # download 38 cnp file from (https://resources.hartwigmedicalfoundation.nl/) and use above
   -ref_genome_version 38 \
   -ref_genome /cluster/tools/data/genomes/human/GRCh38/iGenomes/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta \
   -ensembl_data_dir /path_to_ensembl_data_cache/ \
   -somatic_vcf /path/COLO829/COLO829.somatic.vcf.gz \
   # use Delly?? confirm
   -structural_vcf /path/COLO829/COLO829.sv.high_confidence.vcf.gz \
   # use GRIDSS/delly for SV vcf
   # -sv_recovery_vcf /path/COLO829/COLO829.sv.low_confidence.vcf.gz \
   # skipping sv recovery
   -run_drivers \
   # -driver_gene_panel /path/DriverGenePanel.37.tsv \ 
   #skip the driver_gene_panel
   -circos /path/circos-0.69-6/bin/circos \
   # path to circos binary??
   -output_dir /cluster/projects/pughlab/myeloma/projects/MM_cell_drugs/WGS_Pipeline/purple/ \
----

   # AMBER First
   java -Xmx32G -cp amber.jar com.hartwig.hmftools.amber.AmberApplication \
   -tumor OICRM4CA-07-01-P -tumor_bam /cluster/projects/pughlab/myeloma/external_data/ultimagen-oicr/Crams/OICRM4CA-08-R-P.cram \ 
   -output_dir /cluster/projects/pughlab/myeloma/projects/MM_cell_drugs/WGS_Pipeline/AMBER/ \
   -threads 16 \
   -ref_genome_version 38 \
   -ref_genome /cluster/tools/data/genomes/human/GRCh38/iGenomes/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta \
   -loci /path/to/GermlineHetPon.37.vcf.gz 
----

   # COBALT 2nd
   java -jar -Xmx8G cobalt.jar \
   -ref_genome_version 38 \
   -ref_genome /cluster/tools/data/genomes/human/GRCh38/iGenomes/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta \
   -tumor OICRM4CA-07-01-P \
   -tumor_bam /cluster/projects/pughlab/myeloma/external_data/ultimagen-oicr/Crams/OICRM4CA-08-R-P.cram\ 
   -output_dir /cluster/projects/pughlab/myeloma/projects/MM_cell_drugs/WGS_Pipeline/COBALT/ \ 
   -threads 8 \ 
   -gc_profile /ref_data/GC_profile.1000bp.37.cnp