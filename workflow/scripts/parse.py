"""D"""

from pandas import read_csv, DataFrame
from gtfparse import read_gtf

# Snakemake
GENCODE_GTF = snakemake.input[0]  # type: ignore
OUTPUT = snakemake.output[0]  # type: ignore


def read_gencode(filepath: str) -> DataFrame:
    """Returns DF of metadata"""
    return read_gtf(filepath)


def main() -> None:
    """Main"""
    # Read gtf with gtfparse and write out to tsv
    gtf = read_gtf(GENCODE_GTF)
    
    print(gtf.head())
    
    # Format to bed
    fields = [0, 3, 5, 2, 6, 8, 9, 10]
    fields = ["seqname","start","end","gene_id","gene_id","gene_id","gene_id","gene_id",]
    gtf = gtf[fields]
    
    # Strip "." of gene name
    gtf[5] = [i[0] for i in gtf[5].str.split(".")]
    
    # Save results
    gtf.to_csv(OUTPUT, sep="\t", header=None)




# ------------- #
# Main          #
# ------------- #

if __name__ == "__main__":
    main()
