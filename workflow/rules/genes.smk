# Parameters
INSTALL_DIR = config["INSTALL_DIR"]
PROCESS_DIR = config["PROCESS_DIR"]
GENOME_TAGS = config["GENOME_TAGS"]
FEATURE_LIST = config["FEATURE_LIST"]
PROMOTER_WINDOW = config["PROMOTER_WINDOW"]


rule all:
    input:
        expand(
            PROCESS_DIR + "/{build}/gencode.{build}.{feature}.protein_coding.bed",
            build=GENOME_TAGS,
            feature=FEATURE_LIST,
        ),
    default_target: True


rule subset_protein_coding:
    input:
        PROCESS_DIR + "/{build}/gencode.{build}.basic.annotation.gtf.parsed.stripped.bed",
    output:
        temp(PROCESS_DIR + "/{build}/gencode.{build}_protein_coding.bed"),
    conda:
        "../envs/gencode.yaml"
    shell:
        """
        vawk '{{ if ($7=="protein_coding") print $0}}' {input} > {output}
        """

rule filter_mitochondrial:
    input:
        rules.subset_protein_coding.output,
    output:
        temp(PROCESS_DIR + "/{build}/gencode.{build}_protein_coding.no_chrM.bed"),
    conda:
        "../envs/gencode.yaml"
    shell:
        """
        vawk '{{ if ($1!="chrM") print $0}}' {input} > {output}
        """

rule make_genes:
    input:
        rules.filter_mitochondrial.output,
    output:
        PROCESS_DIR + "/{build}/gencode.{build}.genes.protein_coding.bed",
    conda:
        "../envs/gencode.yaml"
    shell:
        """
        vawk '{{ if ($4=="gene") print $1,$2,$3,$8,".",$5,$6}}' {input} > {output}
        """


rule make_exons:
    input:
        rules.subset_protein_coding.output,
    output:
        PROCESS_DIR + "/{build}/gencode.{build}.exons.protein_coding.bed",
    conda:
        "../envs/gencode.yaml"
    shell:
        """
        vawk '{{ if ($4=="exon") print $1,$2,$3,$8,".",$5,$6}}' {input} > {output}
        """


rule make_proms:
    input:
        genes=PROCESS_DIR + "/{build}/gencode.{build}.genes.protein_coding.bed",
        sizes=INSTALL_DIR + "/{build}.fa.sizes",
    output:
        PROCESS_DIR + "/{build}/gencode.{build}.proms.protein_coding.bed",
    params:
        window=PROMOTER_WINDOW,
    conda:
        "../envs/gencode.yaml"
    shell:
        """
        bedtools flank -i {input.genes} -g {input.sizes} -l {params.window} -r 0 -s > {output}
        """


rule make_tss:
    input:
        genes=rules.make_genes.output,
        sizes=INSTALL_DIR + "/{build}.fa.sizes",
    output:
        PROCESS_DIR + "/{build}/gencode.{build}.tss.protein_coding.bed",
    conda:
        "../envs/gencode.yaml"
    shell:
        """
        bedtools flank -i {input.genes} -g {input.sizes} -l 1 -r 0 -s > {output}
        """