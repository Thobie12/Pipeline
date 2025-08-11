# ichorCNA_Snakefile
configfile: "/cluster/home/t922316uhn/ichorCNA/scripts/snakemake/config/config_final.yaml"
sample_names = list(config["samples"].keys())
rule all:
    input:
        expand(
            "/cluster/projects/pughlab/myeloma/projects/MM_cell_drugs/WGS_Pipeline/ichorCNA/{sample}/{sample}.bin{binSize}.cna.seg",
            sample=sample_names, binSize=config["binSize"]
        ),
        expand(
            "/cluster/projects/pughlab/myeloma/projects/MM_cell_drugs/WGS_Pipeline/ichorCNA/readDepth/{sample}.bin{binSize}.regions.wig",
            sample=sample_names, binSize=config["binSize"]
        )

#rule all:
#    input:
       # expand("/cluster/projects/pughlab/myeloma/projects/MM_cell_drugs/WGS_Pipeline/ichorCNA/{tumor}/{tumor}.cna.seg", tumor=config["samples"]),
        #expand("/cluster/projects/pughlab/myeloma/projects/MM_cell_drugs/WGS_Pipeline/ichorCNA/{sample}/{sample}.bin{binSize}.cna.seg", sample=config["samples"], binSize=str(config["binSize"])),
        #expand("/cluster/projects/pughlab/myeloma/projects/MM_cell_drugs/WGS_Pipeline/ichorCNA/readDepth/{sample}.bin{binSize}.regions.wig", sample=config["samples"], binSize=str(config["binSize"]))

rule mosdepth:
    input:
        cram=lambda wildcards: config["samples"][wildcards.sample],
        ref=config["referenceFasta"]
    output:
        bed="/cluster/projects/pughlab/myeloma/projects/MM_cell_drugs/WGS_Pipeline/ichorCNA/readDepth/{sample}.bin{binSize}.regions.bed.gz",
        wig="/cluster/projects/pughlab/myeloma/projects/MM_cell_drugs/WGS_Pipeline/ichorCNA/readDepth/{sample}.bin{binSize}.regions.wig",
        global_dist="/cluster/projects/pughlab/myeloma/projects/MM_cell_drugs/WGS_Pipeline/ichorCNA/readDepth/{sample}.bin{binSize}.mosdepth.global.dist.txt",
        summary="/cluster/projects/pughlab/myeloma/projects/MM_cell_drugs/WGS_Pipeline/ichorCNA/readDepth/{sample}.bin{binSize}.mosdepth.summary.txt",
        per_base="/cluster/projects/pughlab/myeloma/projects/MM_cell_drugs/WGS_Pipeline/ichorCNA/readDepth/{sample}.bin{binSize}.per-base.bed.gz"
    params:
        mosdepth=config["mosdepth"],
        binSize=config["binSize"],
        outPrefix="/cluster/projects/pughlab/myeloma/projects/MM_cell_drugs/WGS_Pipeline/ichorCNA/readDepth/{sample}.bin{binSize}"
    threads: 8
    resources:
        mem=4
    log:
        "/cluster/projects/pughlab/myeloma/projects/MM_cell_drugs/WGS_Pipeline/ichorCNA/logs/mosdepth/{sample}.bin{binSize}.log"
    shell:
        """
        {params.mosdepth} --fast-mode --by {wildcards.binSize} --threads {threads} --fasta {input.ref} {params.outPrefix} {input.cram} 2> {log}

        zcat {output.bed} | awk '
        BEGIN {{ OFS="\\t" }}
        $1 ~ /^chr([1-9]|1[0-9]|2[0-2]|X)$/ {{
            chr = $1; start = $2; end = $3; cov = $4;
            if (NR == 1 || chr != prev_chr || start != prev_end) {{
                print "fixedStep chrom=" chr " start=" (start+1) " step=" (end - start) " span=" (end - start);
            }}
            print cov;
            prev_chr = chr; prev_end = end;
        }}
        ' > {output.wig}
        """

rule ichorCNA:
    input:
        tum="/cluster/projects/pughlab/myeloma/projects/MM_cell_drugs/WGS_Pipeline/ichorCNA/readDepth/{sample}.bin{binSize}.regions.wig"
    output:
        corrDepth="/cluster/projects/pughlab/myeloma/projects/MM_cell_drugs/WGS_Pipeline/ichorCNA/{sample}/{sample}.bin{binSize}.correctedDepth.txt",
        param="/cluster/projects/pughlab/myeloma/projects/MM_cell_drugs/WGS_Pipeline/ichorCNA/{sample}/{sample}.bin{binSize}.params.txt",
        cna="/cluster/projects/pughlab/myeloma/projects/MM_cell_drugs/WGS_Pipeline/ichorCNA/{sample}/{sample}.bin{binSize}.cna.seg",
        segTxt="/cluster/projects/pughlab/myeloma/projects/MM_cell_drugs/WGS_Pipeline/ichorCNA/{sample}/{sample}.bin{binSize}.seg.txt",
        seg="/cluster/projects/pughlab/myeloma/projects/MM_cell_drugs/WGS_Pipeline/ichorCNA/{sample}/{sample}.bin{binSize}.seg",
        rdata="/cluster/projects/pughlab/myeloma/projects/MM_cell_drugs/WGS_Pipeline/ichorCNA/{sample}/{sample}.bin{binSize}.RData"
    params:
        outDir="/cluster/projects/pughlab/myeloma/projects/MM_cell_drugs/WGS_Pipeline/ichorCNA/{sample}/",
        rscript=config["ichorCNA_rscript"],
        id="{sample}",
        ploidy=config["ichorCNA_ploidy"],
        normal=config["ichorCNA_normal"],
        gcwig=config["ichorCNA_gcWig"],
        mapwig=config["ichorCNA_mapWig"],
        normalpanel=config["ichorCNA_normalPanel"],
        estimateNormal=config["ichorCNA_estimateNormal"],
        estimatePloidy=config["ichorCNA_estimatePloidy"],
        estimateClonality=config["ichorCNA_estimateClonality"],
        scStates=config["ichorCNA_scStates"],
        maxCN=config["ichorCNA_maxCN"],
        includeHOMD=config["ichorCNA_includeHOMD"],
        chrs=config["ichorCNA_chrs"],
        chrTrain=config["ichorCNA_chrTrain"],
        genomeBuild=config["ichorCNA_genomeBuild"],
        genomeStyle=config["ichorCNA_genomeStyle"],
        centromere=config["ichorCNA_centromere"],
        fracReadsChrYMale=config["ichorCNA_fracReadsInChrYForMale"],
        minMapScore=config["ichorCNA_minMapScore"],
        maxFracGenomeSubclone=config["ichorCNA_maxFracGenomeSubclone"],
        maxFracCNASubclone=config["ichorCNA_maxFracCNASubclone"],
        exons=config["ichorCNA_exons"],
        txnE=config["ichorCNA_txnE"],
        txnStrength=config["ichorCNA_txnStrength"],
        plotFileType=config["ichorCNA_plotFileType"],
        plotYlim=config["ichorCNA_plotYlim"],
        libdir=config["ichorCNA_libdir"]
    resources:
        mem=4
    log:
        "/cluster/projects/pughlab/myeloma/projects/MM_cell_drugs/WGS_Pipeline/ichorCNA/logs/ichorCNA/{sample}.bin{binSize}.log"      
    shell:
        """
        module load R && Rscript {params.rscript} \
          --id {params.id} \
          --libdir {params.libdir} \
          --WIG {input.tum} \
          --gcWig {params.gcwig} \
          --mapWig {params.mapwig} \
          --normalPanel {params.normalpanel} \
          --ploidy "{params.ploidy}" \
          --normal "{params.normal}" \
          --maxCN {params.maxCN} \
          --includeHOMD {params.includeHOMD} \
          --chrs "{params.chrs}" \
          --chrTrain "{params.chrTrain}" \
          --genomeStyle {params.genomeStyle} \
          --genomeBuild {params.genomeBuild} \
          --estimateNormal {params.estimateNormal} \
          --estimatePloidy {params.estimatePloidy} \
          --estimateScPrevalence {params.estimateClonality} \
          --scStates "{params.scStates}" \
          --centromere {params.centromere} \
          --exons.bed {params.exons} \
          --txnE {params.txnE} \
          --txnStrength {params.txnStrength} \
          --minMapScore {params.minMapScore} \
