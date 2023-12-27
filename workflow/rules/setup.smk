# Parameters
UCSC_URL_HG19 = config["GENCODE"]["UCSC_URLS"]["hg19"]
UCSC_URL_HG38 = config["GENCODE"]["UCSC_URLS"]["hg38"]
GENCODE_URL_HG19 = config["GENCODE"]["GENCODE_URLS"]["hg19"]
GENCODE_URL_HG38 = config["GENCODE"]["GENCODE_URLS"]["hg38"]


rule download_gencode:
    message:
        """
        Downloads gencode basic build based off of input assembly
        """
    output:
        "resources/gencode/{ASSEMBLY}/gencode.{ASSEMBLY}.basic.annotation.gtf.gz",
    params:
        assembly=lambda wc: wc.ASSEMBLY,
        hg19_url=GENCODE_URL_HG19,
        hg38_url=GENCODE_URL_HG38,
    log:
        stdout="workflow/logs/download_gencode-{ASSEMBLY}.stdout",
        stderr="workflow/logs/download_gencode-{ASSEMBLY}.stderr",
    threads: 1
    shell:
        """
        if [ {params.assembly} == "hg19" ]; then
            curl {params.hg19_url} -o {output}
        else
            curl {params.hg38_url} -o {output}
        fi
        """


rule download_chromesome_sizes:
    message:
        """
        Downloads chromosome sizes for BedTools commands
        """
    output:
        "resources/gencode/{ASSEMBLY}/{ASSEMBLY}.fa.sizes",
    params:
        assembly=lambda wc: wc.ASSEMBLY,
        hg19_url=UCSC_URL_HG19,
        hg38_url=UCSC_URL_HG38,
    log:
        stdout="workflow/logs/download_chromesome_sizes-{ASSEMBLY}.stdout",
        stderr="workflow/logs/download_chromesome_sizes-{ASSEMBLY}.stderr",
    threads: 1
    shell:
        """
        if [ {params.assembly} == "hg19" ]; then
            curl {params.hg19_url} -o {output}
        else
            curl {params.hg38_url} -o {output}
        fi
        """
