# Parameters
BUILDS = config["GENCODE"]["BUILDS"]
PROMOTER_WINDOW = config["GENCODE"]["PROMOTER_WINDOW"]

rule all:
    input:
        expand(
            "results/gencode/{ASSEMBLY}/gencode.{ASSEMBLY}.genes.protein_coding.bed",
            ASSEMBLY=BUILDS,
        ),
        expand(
            "results/gencode/{ASSEMBLY}/gencode.{ASSEMBLY}.exons.protein_coding.bed",
            ASSEMBLY=BUILDS,
        ),
        expand(
            "results/gencode/{ASSEMBLY}/gencode.{ASSEMBLY}.proms.protein_coding.bed",
            ASSEMBLY=BUILDS,
        ),
        expand(
            "results/gencode/{ASSEMBLY}/gencode.{ASSEMBLY}.tss.protein_coding.bed",
            ASSEMBLY=BUILDS,
        ),
        expand(
            "results/gencode/{ASSEMBLY}/gencode.{ASSEMBLY}.transcripts.protein_coding.bed",
            ASSEMBLY=BUILDS,
        ),


rule parse_gencode:
    message:
        """
        Parses gtf formatted raw gencode download. Notes:
        - Formats to BED
        - Strips gene versions
        """
    input:
        "resources/gencode/{ASSEMBLY}/gencode.{ASSEMBLY}.basic.annotation.gtf.gz",
    output:
        temp("results/gencode/{ASSEMBLY}/gencode.{ASSEMBLY}.basic.annotation.gtf.parsed.tsv"),
    log:
        stdout="workflow/logs/parse_gencode-{ASSEMBLY}.stdout",
        stderr="workflow/logs/parse_gencode-{ASSEMBLY}.stderr",
    conda:
        "../envs/gencode.yaml"
    threads: 8
    script:
        "../scripts/parse.py"


rule filter_gencode:
    message:
        """
        Subsets Gencode to protein-coding features and filters mitochondrial
        """
    input:
        rules.parse_gencode.output,
    output:
        temp("results/gencode/{ASSEMBLY}/gencode.{ASSEMBLY}_protein_coding.no_chrM.bed"),
    conda:
        "../envs/gencode.yaml"
    log:
        stdout="workflow/logs/filter_gencode-{ASSEMBLY}.stdout",
        stderr="workflow/logs/filter_gencode-{ASSEMBLY}.stderr",
    shell:
        """
        vawk '{{ if ($7=="protein_coding" && $1!="chrM") print $0}}' {input} > {output}
        """


rule make_genes:
    message:
        """
        Makes genes
        """
    input:
        rules.filter_gencode.output,
    output:
        "results/gencode/{ASSEMBLY}/gencode.{ASSEMBLY}.genes.protein_coding.bed",
    conda:
        "../envs/gencode.yaml"
    log:
        stdout="workflow/logs/make_genes-{ASSEMBLY}.stdout",
        stderr="workflow/logs/make_genes-{ASSEMBLY}.stderr",
    shell:
        """
        vawk '{{ if ($9=="gene") print $1,$2,$3,$4,".",$5, $6}}' {input} > {output}
        """


rule make_exons:
    message:
        """
        Makes exons
        """
    input:
        rules.filter_gencode.output,
    output:
        "results/gencode/{ASSEMBLY}/gencode.{ASSEMBLY}.exons.protein_coding.bed",
    conda:
        "../envs/gencode.yaml"
    log:
        stdout="workflow/logs/make_exons-{ASSEMBLY}.stdout",
        stderr="workflow/logs/make_exons-{ASSEMBLY}.stderr",
    shell:
        """
        vawk '{{ if ($9=="exon") print $1,$2,$3,$4,".",$5, $6}}' {input} > {output}
        """


rule make_proms:
    message:
        """
        Makes promoters
        """
    input:
        genes=rules.make_genes.output,
        sizes="resources/gencode/{ASSEMBLY}/{ASSEMBLY}.fa.sizes",
    output:
        "results/gencode/{ASSEMBLY}/gencode.{ASSEMBLY}.proms.protein_coding.bed",
    params:
        window=PROMOTER_WINDOW,
    conda:
        "../envs/gencode.yaml"
    log:
        stdout="workflow/logs/make_proms-{ASSEMBLY}.stdout",
        stderr="workflow/logs/make_proms-{ASSEMBLY}.stderr",
    shell:
        """
        bedtools flank -i {input.genes} -g {input.sizes} -l {params.window} -r 0 -s > {output}
        """


rule make_tss:
    message:
        """
        Makes TSS
        """
    input:
        genes=rules.make_genes.output,
        sizes="resources/gencode/{ASSEMBLY}/{ASSEMBLY}.fa.sizes",
    output:
        "results/gencode/{ASSEMBLY}/gencode.{ASSEMBLY}.tss.protein_coding.bed",
    conda:
        "../envs/gencode.yaml"
    log:
        stdout="workflow/logs/make_tss-{ASSEMBLY}.stdout",
        stderr="workflow/logs/make_tss-{ASSEMBLY}.stderr",
    shell:
        """
        bedtools flank -i {input.genes} -g {input.sizes} -l 1 -r 0 -s > {output}
        """


rule make_transcripts:
    message:
        """
        Makes transcripts fileset. Prints trasncript name instead of gene in name field.
        """
    input:
        genes=rules.filter_gencode.output,
    output:
        "results/gencode/{ASSEMBLY}/gencode.{ASSEMBLY}.transcripts.protein_coding.bed",
    conda:
        "../envs/gencode.yaml"
    log:
        stdout="workflow/logs/make_transcripts-{ASSEMBLY}.stdout",
        stderr="workflow/logs/make_transcripts-{ASSEMBLY}.stderr",
    shell:
        """
        vawk '{{ if ($9=="gene") print $1,$2,$3,$8,".",$5, $6}}' {input} > {output}
        """
