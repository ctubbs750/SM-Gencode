# Parameters
INSTALL_DIR = config["GENCODE"]["INSTALL_DIR"]
PROCESS_DIR = config["GENCODE"]["PROCESS_DIR"]
ASSEMBLY = config["GENCODE"]["ASSEMBLY"]


rule parse_gencode:
    message:
        """
        Parses gtf formatted raw gencode download. Notes:
        - Formats to BED
        - Strips gene versions
        """
    input:
        INSTALL_DIR + f"/gencode.{ASSEMBLY}.basic.annotation.gtf.gz",
    output:
        temp(
            PROCESS_DIR
            + f"/{ASSEMBLY}/gencode.{ASSEMBLY}.basic.annotation.gtf.parsed.tsv"
        ),
    log:
        stdout=f"workflow/logs/parse_gencode-{ASSEMBLY}.stdout",
        stderr=f"workflow/logs/parse_gencode-{ASSEMBLY}.stderr",
    conda:
        "../envs/gencode.yaml"
    threads: 8
    script:
        "../scripts/parse.py"
