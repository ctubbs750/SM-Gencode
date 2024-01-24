from os import listdir, path
from snakemake.utils import min_version


# Settings
min_version("7.32.4")


# ------------- #
# Config        #
# ------------- #

BUILDS = config["builds"]
INSTALL_DIR = config["install_dir"]
PROCESS_DIR = config["process_dir"]
UCSC_URL_HG19 = config["ucsc_urls"]["hg19"]
UCSC_URL_HG38 = config["ucsc_urls"]["hg38"]
PROMOTER_WINDOW = config["promoter_window"]
GENCODE_URL_HG19 = config["gencode_urls"]["hg19"]
GENCODE_URL_HG38 = config["gencode_urls"]["hg38"]

# ------------- #
# I/O           #
# ------------- #

# Raw gencode download
GENCODE_DOWNLOAD = path.join(
    INSTALL_DIR, "{build}", "gencode.{build}.basic.annotation.gtf.gz"
)

# Chromosome sizes from UCSC
CHROMOSOME_SIZES = path.join(INSTALL_DIR, "{build}", "{build}.fa.size")

# Parsed gencode
GENCODE_PARSED = path.join(
    PROCESS_DIR, "{build}", "gencode.{build}.basic.annotation.gtf.parsed.tsv"
)

# Filtered gencode
GENCODE_FILTERED = path.join(
    PROCESS_DIR, "{build}", "gencode.{build}_protein_coding.no_chrM.bed"
)

# Final target features
GENES = path.join(PROCESS_DIR, "{build}", "gencode.{build}.genes.protein_coding.bed")
EXONS = path.join(PROCESS_DIR, "{build}", "gencode.{build}.exons.protein_coding.bed")
PROMS = path.join(PROCESS_DIR, "{build}", "gencode.{build}.proms.protein_coding.bed")
TSS = path.join(PROCESS_DIR, "{build}", "gencode.{build}.tss.protein_coding.bed")
TRANSCRIPTS = path.join(
    PROCESS_DIR, "{build}", "gencode.{build}.transcripts.protein_coding.bed"
)

# ------------- #
# Rules         #
# ------------- #


rule all:
    input:
        expand(
            GENES,
            ASSEMBLY=BUILDS,
        ),
        expand(
            EXONS,
            ASSEMBLY=BUILDS,
        ),
        expand(
            PROMS,
            ASSEMBLY=BUILDS,
        ),
        expand(
            TSS,
            ASSEMBLY=BUILDS,
        ),
        expand(
            TRANSCRIPTS,
            ASSEMBLY=BUILDS,
        ),


rule download_gencode:
    message:
        """
        Downloads gencode basic build based off of input assembly.
        """
    output:
        GENCODE_DOWNLOAD,
    params:
        build=lambda wc: wc.build,
        hg19_url=GENCODE_URL_HG19,
        hg38_url=GENCODE_URL_HG38,
    conda:
        "../envs/gencode.yaml"
    log:
        stdout="workflow/logs/download_gencode-{build}.stdout",
        stderr="workflow/logs/download_gencode-{build}.stderr",
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
        CHROMOSOME_SIZES,
    params:
        build=lambda wc: wc.BUILD,
        hg19_url=UCSC_URL_HG19,
        hg38_url=UCSC_URL_HG38,
    conda:
        "../envs/gencode.yaml"
    log:
        stdout="workflow/logs/download_sizes-{build}.stdout",
        stderr="workflow/logs/download_sizes-{build}.stderr",
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
        temp(GENCODE_PARSED),
    conda:
        "../envs/gencode.yaml"
    log:
        stdout="workflow/logs/parse_gencode-{build}.stdout",
        stderr="workflow/logs/parse_gencode-{build}.stderr",
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
        temp(GENCODE_FILTERED),
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
        GENES,
    conda:
        "../envs/gencode.yaml"
    log:
        stdout="workflow/logs/make_genes-{build}.stdout",
        stderr="workflow/logs/make_genes-{build}.stderr",
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
        EXONS,
    conda:
        "../envs/gencode.yaml"
    log:
        stdout="workflow/logs/make_exons-{build}.stdout",
        stderr="workflow/logs/make_exons-{build}.stderr",
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
        PROMS,
    params:
        window=PROMOTER_WINDOW,
    conda:
        "../envs/gencode.yaml"
    log:
        stdout="workflow/logs/make_proms-{build}.stdout",
        stderr="workflow/logs/make_proms-{build}.stderr",
    threads: 1
    shell:
        """
        bedtools flank -i {input.genes} -g {input.sizes} -l {params.window} -r 0 -s > {output}
        """


rule make_TSS:
    message:
        """
        Makes transcription start sites for protein-coding genes.
        """
    input:
        genes=rules.make_genes.output,
        sizes=rules.download_sizes.output,
    output:
        TSS,
    conda:
        "../envs/gencode.yaml"
    log:
        stdout="workflow/logs/make_tss-{build}.stdout",
        stderr="workflow/logs/make_tss-{build}.stderr",
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
        TRANSCRIPTS,
    conda:
        "../envs/gencode.yaml"
    log:
        stdout="workflow/logs/make_transcripts-{build}.stdout",
        stderr="workflow/logs/make_transcripts-{build}.stderr",
    threads: 1
    shell:
        """
        vawk '{{ if ($9=="gene") print $1, $2, $3, $8, ".", $5, $6 }}' {input} > {output}
        """
