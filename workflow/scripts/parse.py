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
    gtf = gtf[
        [
            "seqname",
            "start",
            "end",
            "gene_id",
            "strand",
            "gene_name",
            "gene_type",
            "transcript_id",
            "feature",
        ]
    ]

    # Strip "." of gene and transcript IDs
    gtf["gene_id"] = gtf["gene_id"].apply(lambda x: x.split(".")[0])
    gtf["transcript_id"] = gtf["transcript_id"].apply(lambda x: x.split(".")[0])

    # Save results
    with open(OUT, "w") as f:
        gtf.to_csv(f, sep="\t", header=False, index=False)


# ------------- #
# Main          #
# ------------- #

if __name__ == "__main__":
    main()
