"""Script to parse input GTF file using gtfparse"""

from gtfparse import read_gtf

# Snakemake
GTF = snakemake.input[0]  # type: ignore
OUT = snakemake.output[0]  # type: ignore

def main() -> None:
    """Reads GTF, subsets to fields, and cleans gene/transcript names"""
    # Read gtf with gtfparse
    gtf = read_gtf(GTF)
    
    # Format to bed
    fields = ["seqname","start","end","gene_id", "strand","gene_name","gene_type", "transcript_id", "feature"]
    gtf = gtf[fields]
    
    # Strip "." of gene and transcript IDs
    gtf["gene_id"] = [i[0] for i in gtf["gene_id"].str.split(".")]
    gtf["transcript_id"] = [i[0] for i in gtf["transcript_id"].str.split(".")]
    
    # Save results
    gtf.to_csv(OUT, sep="\t", header=False, index=False)


# ------------- #
# Main          #
# ------------- #

if __name__ == "__main__":
    main()
