from snakemake.utils import min_version

# Parameters
BUILDS = config["builds"]
UCSC_URL_HG19 = config["ucsc_urls"]["hg19"]
UCSC_URL_HG38 = config["ucsc_urls"]["hg38"]
PROMOTER_WINDOW = config["promoter_window"]
GENCODE_URL_HG19 = config["gencode_urls"]["hg19"]
GENCODE_URL_HG38 = config["gencode_urls"]["hg38"]

# Settings
min_version("7.32.4")


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


rule download_gencode:
    message:
        """
        Downloads gencode basic build based off of input assembly.
        """
    output:
        "resources/gencode/{BUILD}/gencode.{BUILD}.basic.annotation.gtf.gz",
    params:
        build=lambda wc: wc.BUILD,
        hg19_url=GENCODE_URL_HG19,
        hg38_url=GENCODE_URL_HG38,
    conda:
        "../envs/gencode.yaml"
    log:
        stdout="workflow/logs/download_gencode-{BUILD}.stdout",
        stderr="workflow/logs/download_gencode-{BUILD}.stderr",
    threads: 1
    shell:
        """
        if [ {params.build} == "hg19" ]
        then
            wget -O {output} {params.hg19_url}
        else
            wget -O {output} {params.hg38_url}
        fi
        """


rule download_sizes:
    message:
        """
        Downloads chromosome sizes for BedTools commands.
        """
    output:
        "resources/gencode/{BUILD}/{BUILD}.fa.sizes",
    params:
        build=lambda wc: wc.BUILD,
        hg19_url=UCSC_URL_HG19,
        hg38_url=UCSC_URL_HG38,
    conda:
        "../envs/gencode.yaml"
    log:
        stdout="workflow/logs/download_sizes-{BUILD}.stdout",
        stderr="workflow/logs/download_sizes-{BUILD}.stderr",
    threads: 1
    shell:
        """
        if [ {params.build} == "hg19" ]
        then
            wget -O {output} {params.hg19_url}
        else
            wget -O {output} {params.hg38_url}
        fi
        """


rule parse_gencode:
    message:
        """
        Parses gtf formatted raw gencode download. Notes:
        - Formats to BED
        - Strips gene versions
        """
    input:
        rules.download_gencode.output,
    output:
        temp("results/gencode/{BUILD}/gencode.{BUILD}.basic.annotation.gtf.parsed.tsv"),
    conda:
        "../envs/gencode.yaml"
    log:
        stdout="workflow/logs/parse_gencode-{BUILD}.stdout",
        stderr="workflow/logs/parse_gencode-{BUILD}.stderr",
    threads: 4
    script:
        "../scripts/parse.py"


rule filter_gencode:
    message:
        """
        Subsets Gencode to protein-coding features and filters mitochondrial.
        """
    input:
        rules.parse_gencode.output,
    output:
        temp("results/gencode/{BUILD}/gencode.{BUILD}_protein_coding.no_chrM.bed"),
    conda:
        "../envs/gencode.yaml"
    log:
        stdout="workflow/logs/filter_gencode-{BUILD}.stdout",
        stderr="workflow/logs/filter_gencode-{BUILD}.stderr",
    threads: 1
    shell:
        """
        vawk '{{ if ($7=="protein_coding" && $1!="chrM") print $0 }}' {input} > {output}
        """


rule make_genes:
    message:
        """
        Makes protein-coding gene BED file from filtered GENCODE file.
        """
    input:
        rules.filter_gencode.output,
    output:
        "results/gencode/{BUILD}/gencode.{BUILD}.genes.protein_coding.bed",
    conda:
        "../envs/gencode.yaml"
    log:
        stdout="workflow/logs/make_genes-{BUILD}.stdout",
        stderr="workflow/logs/make_genes-{BUILD}.stderr",
    threads: 1
    shell:
        """
        vawk '{{ if ($9=="gene") print $1, $2, $3, $4, ".", $5, $6 }}' {input} > {output}
        """


rule make_exons:
    message:
        """
        Makes protein-coding exon BED file from filtered GENCODE file.
        """
    input:
        rules.filter_gencode.output,
    output:
        "results/gencode/{BUILD}/gencode.{BUILD}.exons.protein_coding.bed",
    conda:
        "../envs/gencode.yaml"
    log:
        stdout="workflow/logs/make_exons-{BUILD}.stdout",
        stderr="workflow/logs/make_exons-{BUILD}.stderr",
    threads: 1
    shell:
        """
        vawk '{{ if ($9=="exon") print $1, $2, $3, $4, ".", $5, $6 }}' {input} > {output}
        """


rule make_promoters:
    message:
        """
        Makes protein-coding exon BED file from gene file.
        """
    input:
        genes=rules.make_genes.output,
        sizes=rules.download_sizes.output,
    output:
        "results/gencode/{BUILD}/gencode.{BUILD}.proms.protein_coding.bed",
    params:
        window=PROMOTER_WINDOW,
    conda:
        "../envs/gencode.yaml"
    log:
        stdout="workflow/logs/make_proms-{BUILD}.stdout",
        stderr="workflow/logs/make_proms-{BUILD}.stderr",
    threads: 1
    shell:
        """
        bedtools flank -i {input.genes} -g {input.sizes} -l {params.window} -r 0 -s > {output}
        """


rule make_TSSs:
    message:
        """
        Makes transcription start sites for protein-coding genes.
        """
    input:
        genes=rules.make_genes.output,
        sizes=rules.download_sizes.output,
    output:
        "results/gencode/{BUILD}/gencode.{BUILD}.tss.protein_coding.bed",
    conda:
        "../envs/gencode.yaml"
    log:
        stdout="workflow/logs/make_tss-{BUILD}.stdout",
        stderr="workflow/logs/make_tss-{BUILD}.stderr",
    threads: 1
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
        rules.filter_gencode.output,
    output:
        "results/gencode/{BUILD}/gencode.{BUILD}.transcripts.protein_coding.bed",
    conda:
        "../envs/gencode.yaml"
    log:
        stdout="workflow/logs/make_transcripts-{BUILD}.stdout",
        stderr="workflow/logs/make_transcripts-{BUILD}.stderr",
    threads: 1
    shell:
        """
        vawk '{{ if ($9=="gene") print $1, $2, $3, $8, ".", $5, $6 }}' {input} > {output}
        """
