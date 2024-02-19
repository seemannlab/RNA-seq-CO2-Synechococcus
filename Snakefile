import pandas as pd

from snakemake.utils import validate, min_version
min_version("7.9.0")

################################################################################
# Copy files from pipelines for easier processing

rule setup:
    output: 
        directory('raw-data')
    shell:
        """
# Shorter paths
P=/home/projects/rth/co2capture/subprojects/RNA-seq_air-to-30-pct-CO2
P2=$P/Public-RNA-Seq/RNA

mkdir raw-data
cd raw-data

# Annotation
cp $P/PCC7002-genome.gff.gz                         PCC7002-genome.gff.gz
cp $P/PCC7002-genome.fna.gz                         PCC7002-genome.fna.gz

# Growth rate etc
cp $P/experimental-data/{{RNA_JAN23.xlsx,RNA_July22.xlsx}} .

# Expression matrices project data
cp $P/RNA-browser/analysis/41_counts.tsv            1-15pct-co2-counts.tsv
cp $P/RNA-browser-dataset2/analysis/41_counts.tsv   co2-adaptation-counts.tsv

# Info on batches
cp $P/RNA-Schlange/samples.csv                      schlange-data1.csv

# MultiQC stat
cp $P/RNA-browser/multiqc/multiqc_data/multiqc_data.json                            browser-data1-multiqc.json
cp $P/RNA-browser-dataset2/multiqc/multiqc_data/multiqc_data.json                   browser-data2-multiqc.json
cp $P/RNA-Schlange/analysis/53_main_multiqc/multiqc_data/multiqc_data.json          schlange-data1-multiqc.json
cp $P/RNA-Schlange-Dataset2/analysis/53_main_multiqc/multiqc_data/multiqc_data.json schlange-data2-multiqc.json

# Public datasets info
mkdir public-datasets
cd public-datasets

cp $P2-browser/analysis/41_counts.tsv                                            public-raw-counts.tsv
cp $P2-Schlange-Illumina/analysis/53_main_multiqc/multiqc_data/multiqc_data.json multiqc-schlange-illumina.json
cp $P2-Schlange-SOLiD/analysis/53_main_multiqc/multiqc_data/multiqc_data.json    multiqc-schlange-solid.json
cp $P2-browser/multiqc/multiqc_data/multiqc_data.json                            multiqc-public-browser.json
        """

################################################################################
# stats on filtering and mapping

rule A_general:
    input:
        script = 'scripts/A_pre-stat.R',
        xs = [
            'raw-data/schlange-data1-multiqc.json',
            'raw-data/schlange-data2-multiqc.json',
            'raw-data/browser-data1-multiqc.json',
            'raw-data/browser-data2-multiqc.json',
            'raw-data/public-datasets/multiqc-schlange-illumina.json',
            'raw-data/public-datasets/multiqc-schlange-solid.json',
            'raw-data/public-datasets/multiqc-public-browser.json'
        ]
    output:
        'data/A_general-stat.tsv',
        'data/A_public-stat.tsv'
    shell:
        "{input.script}"

################################################################################
# process wet lab measurements

rule B_characteristics:
    input:
        script = 'scripts/B_characteristics.R',
        xs = [
          'raw-data/RNA_July22.xlsx',
          'raw-data/RNA_JAN23.xlsx'
        ]
    output:
        'data/B_characteristics.tsv',
        'data/B_growth.jpeg',
    shell:
        "{input.script}"

################################################################################
# Prepare raw-input

rule C_prep:
    threads: 4
    input:
        script = 'scripts/C_prep.R',
        xs = [
            'raw-data/schlange-data1.csv',
            'raw-data/1-15pct-co2-counts.tsv',
            'raw-data/co2-adaptation-counts.tsv',
            'data/B_characteristics.tsv',
            'raw-data/public-datasets/public-raw-counts.tsv'
        ]
    output:
        'data/C_meta.tsv',
        'data/C_raw-counts.tsv',
        'data/C_annotation.tsv',
        'data/C_public-meta.tsv',
        'data/C_public-raw-counts.tsv'
    shell:
        "{input.script}"
        
        
################################################################################
# Differential expression analysis per dataset for comparison

rule D_high_low:
    input:
        script = 'scripts/D_high-low.R',
        xs = [
            'scripts/helper_deg.R',
            'data/C_meta.tsv',
            'data/C_annotation.tsv',
            'data/C_raw-counts.tsv'
        ]
    output:
        'analysis/D_dge-per-dataset.tsv',
        'analysis/D_PCA-per-dataset.jpeg',
        'analysis/D_correlation-logFCs.jpeg',
        'analysis/D_venn-per-dataset.jpeg',
    shell:
        "{input.script}"


################################################################################
################################################################################

rule all:
    input:
        'data/A_general-stat.tsv',
        'data/B_characteristics.tsv',
        'data/C_raw-counts.tsv',
        'analysis/D_dge-per-dataset.tsv'
