"""
Ruleset for parsing gtf formats.
"""

from gtfparse import read_gtf

configfile: "../../config/config.yaml"

# Parameters
INSTALL_DIR = config["INSTALL_DIR"]
PROCESS_DIR = config["PROCESS_DIR"]
GENOME_TAGS = config["GENOME_TAGS"]


rule all:
    input:
        expand(PROCESS_DIR + "/{build}/gencode.{build}.basic.annotation.gtf.parsed.stripped.bed", build=GENOME_TAGS)
    default_target:
        True


rule parse_gencode:
    message:
        """
        Parses gtf format
        """
    input:
        INSTALL_DIR + "/gencode.{build}.basic.annotation.gtf.gz",
    output:
        temp(PROCESS_DIR + "/{build}/gencode.{build}.basic.annotation.gtf.parsed.tsv"),
    run:
        # Read gtf with gtfparse and write out to tsv
        gtf = read_gtf(input[0])
        gtf.to_csv(output[0], sep="\t", index=False, header=None)


rule format_2bed:
    message:
        """
        Formats gtf to UCSC bed, retaining minimal info
        """
    input:
        rules.parse_gencode.output,
    output:
        temp(PROCESS_DIR + "/{build}/gencode.{build}.basic.annotation.gtf.parsed.bed"),
    conda:
        "../envs/gencode.yaml"
    shell:
        """
        vawk '{{print $1, $4, $5, $3, $7, $9, $10, $11}}' {input} > {output}
        """

rule strip_gene_versions:
    message:
        """
        Strips "." off of ENSEMBL gene ids
        """
    input:
        rules.format_2bed.output,
    output:
        temp(PROCESS_DIR + "/{build}/gencode.{build}.basic.annotation.gtf.parsed.stripped.bed"),
    conda:
        "../envs/gencode.yaml"
    shell:
        """
        vawk '{{split($6, gene, "."); print $1, $2, $3, $4, $5, gene[1], $7, $8}}' {input} > {output}
        """