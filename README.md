# VSMC_MPRA
This repo contains the code use to analyze and visualize the data collected from MPRAs performed on VSMCs for BioRxiv 10.1101/2024.10.07.617126

As a brief summary: massively parallel reporter assays (MPRAs) are a type of molecular assay used to simultaneously screen thousands of candidate 
sequences for their requlatory activity. In this particular case, the analyses here were performed on an MPRA library testing whether specific 
coronary artery disease-associated SNPs showed an allelic imbalance in their enhancer activity. (For specifics of library design and experimental
design, see the above BioRxiv link.)

The input files for the initial analyses are the raw fastq files generated from NGS. These are processed with MPRAflow to generate a count table
of DNA and RNA barcode counts. These count tables are used for the subsequent analyses presented here.

# Repo Contents

## Codes
The R codes contained herein are either analysis codes (using the MPRAnalyze framework) or visualization codes. The visualization codes can be 
used to replicate the figures from the above BioRxiv preprint. Specifically, figures 2 and 3 which correspond to this assay.

## Data
Some of the codes reference RData files, which are also included.
