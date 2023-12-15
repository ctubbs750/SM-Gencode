# Parameters
INSTALL_DIR = config["GENCODE"]["INSTALL_DIR"]
PROCESS_DIR = config["GENCODE"]["PROCESS_DIR"]
ASSEMBLY = config["GENCODE"]["ASSEMBLY"]


rule download_gencode:
    message:
        """
        Downloads gencode basic build based off of input assembly
        """
    output:
        INSTALL_DIR + f"/gencode.{ASSEMBLY}.basic.annotation.gtf.gz",
    params:
        assembly=ASSEMBLY,
        url=config["GENCODE"][f"GENCODE_{ASSEMBLY}"],
    log:
        stdout=f"workflow/logs/download_gencode-{ASSEMBLY}.stdout",
        stderr=f"workflow/logs/download_gencode-{ASSEMBLY}.stderr",
    shell:
        """
        curl {params.url} -o {output}
        """


rule download_chromesome_sizes:
    message:
        """
        Downloads chromosome sizes for BedTools commands
        """
    output:
        INSTALL_DIR + f"/{ASSEMBLY}.fa.sizes",
    params:
        assembly=ASSEMBLY,
        url=config["GENCODE"][f"CHROM_SIZES_{ASSEMBLY}"],
    log:
        stdout=f"workflow/logs/download_chromesome_sizes-{ASSEMBLY}.stdout",
        stderr=f"workflow/logs/download_chromesome_sizes-{ASSEMBLY}.stderr",
    shell:
        """
        curl {params.url} -o {output}
        """
