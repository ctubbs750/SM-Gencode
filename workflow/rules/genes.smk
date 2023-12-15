# Parameters
INSTALL_DIR = config["GENCODE"]["INSTALL_DIR"]
PROCESS_DIR = config["GENCODE"]["PROCESS_DIR"]
ASSEMBLY = config["GENCODE"]["ASSEMBLY"]
PROMOTER_WINDOW = config["GENCODE"]["PROMOTER_WINDOW"]


rule all:
    input:
        PROCESS_DIR + f"/{ASSEMBLY}/gencode.{ASSEMBLY}.genes.protein_coding.bed",
        PROCESS_DIR + f"/{ASSEMBLY}/gencode.{ASSEMBLY}.exons.protein_coding.bed",
        PROCESS_DIR + f"/{ASSEMBLY}/gencode.{ASSEMBLY}.proms.protein_coding.bed",
        PROCESS_DIR + f"/{ASSEMBLY}/gencode.{ASSEMBLY}.tss.protein_coding.bed",


rule filter_gencode:
    message:
        """
        Subsets Gencode to protein-coding features and filters mitochondrial
        """
    input:
        PROCESS_DIR + f"/{ASSEMBLY}/gencode.{ASSEMBLY}.basic.annotation.gtf.parsed.tsv",
    output:
        temp(PROCESS_DIR + f"/{ASSEMBLY}/gencode.{ASSEMBLY}_protein_coding.no_chrM.bed"),
    conda:
        "../envs/gencode.yaml"
    log:
        stdout=f"workflow/logs/filter_gencode-{ASSEMBLY}.stdout",
        stderr=f"workflow/logs/filter_gencode-{ASSEMBLY}.stderr",
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
        PROCESS_DIR + f"/{ASSEMBLY}/gencode.{ASSEMBLY}.genes.protein_coding.bed",
    conda:
        "../envs/gencode.yaml"
    log:
        stdout=f"workflow/logs/make_genes-{ASSEMBLY}.stdout",
        stderr=f"workflow/logs/make_genes-{ASSEMBLY}.stderr",
    shell:
        """
        vawk '{{ if ($4=="gene") print $1,$2,$3,$8,".",$5,$6}}' {input} > {output}
        """


rule make_exons:
    message:
        """
        Makes exons
        """
    input:
        rules.filter_gencode.output,
    output:
        PROCESS_DIR + f"/{ASSEMBLY}/gencode.{ASSEMBLY}.exons.protein_coding.bed",
    conda:
        "../envs/gencode.yaml"
    log:
        stdout=f"workflow/logs/make_exons-{ASSEMBLY}.stdout",
        stderr=f"workflow/logs/make_exons-{ASSEMBLY}.stderr",
    shell:
        """
        vawk '{{ if ($4=="exon") print $1,$2,$3,$8,".",$5,$6}}' {input} > {output}
        """


rule make_proms:
    message:
        """
        Makes promoters
        """
    input:
        genes=PROCESS_DIR + f"/{ASSEMBLY}/gencode.{ASSEMBLY}.genes.protein_coding.bed",
        sizes=INSTALL_DIR + f"/{ASSEMBLY}.fa.sizes",
    output:
        PROCESS_DIR + f"/{ASSEMBLY}/gencode.{ASSEMBLY}.proms.protein_coding.bed",
    params:
        window=PROMOTER_WINDOW,
    conda:
        "../envs/gencode.yaml"
    log:
        stdout=f"workflow/logs/make_proms-{ASSEMBLY}.stdout",
        stderr=f"workflow/logs/make_proms-{ASSEMBLY}.stderr",
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
        sizes=INSTALL_DIR + f"/{ASSEMBLY}.fa.sizes",
    output:
        PROCESS_DIR + f"/{ASSEMBLY}/gencode.{ASSEMBLY}.tss.protein_coding.bed",
    conda:
        "../envs/gencode.yaml"
    log:
        stdout=f"workflow/logs/make_tss-{ASSEMBLY}.stdout",
        stderr=f"workflow/logs/make_tss-{ASSEMBLY}.stderr",
    shell:
        """
        bedtools flank -i {input.genes} -g {input.sizes} -l 1 -r 0 -s > {output}
        """

rule make_transcripts:
    **could be helpful for cmc.