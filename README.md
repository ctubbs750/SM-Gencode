# Gencode

A Snakemake module for installing gene and related annotations.

# Notes

This module is designed to be incorporated into an existing Snakemake pipeline for integrative analysis. To use it, you can point towards its Snakefile raw url during module import:

`https://raw.githubusercontent.com/ctubbs750/gencode/main/workflow/Snakefile`

 By default, it will install annotation files into the existing project directory under:

 `/resources/data/gencode/{build}`

 For this module to function, it will look for parameters under the `GENCODE` tag in config file.  Copy and paste the config parameters found in `/config/config.yaml` into the desired project's config file and edit as needed.