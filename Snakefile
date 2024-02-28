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

mkdir raw-data
cd raw-data

# Annotation
cp $P/PCC7002-genome.gff.gz                         PCC7002-genome.gff.gz
cp $P/PCC7002-genome.fna.gz                         PCC7002-genome.fna.gz

# Growth rate etc
cp $P/experimental-data/RNA_JAN23.xlsx .

# Expression matrices project data
cp $P/RNA-browser-dataset2/analysis/41_counts.tsv   co2-adaptation-counts.tsv

# MultiQC stat
cp $P/RNA-browser-dataset2/multiqc/multiqc_data/multiqc_data.json                   browser-data-multiqc.json
cp $P/RNA-Schlange-Dataset2/analysis/53_main_multiqc/multiqc_data/multiqc_data.json schlange-data-multiqc.json
        """

################################################################################
# stats on filtering and mapping

rule A_general:
    input:
        script = 'scripts/A_pre-stat.R',
        xs = [
            'raw-data/schlange-data-multiqc.json',
            'raw-data/browser-data-multiqc.json'
        ]
    output:
        'data/A_general-stat.tsv'
    log: 'logs/A_general.txt'
    shell:
        "Rscript {input.script} > {log}"

################################################################################
# process wet lab measurements

rule B_characteristics:
    input:
        script = 'scripts/B_characteristics.R',
        xs = [
          'raw-data/RNA_JAN23.xlsx'
        ]
    output:
        'data/B_characteristics.tsv',
        'data/B_growth.jpeg',
    log: 'logs/B_characteristics.txt'
    shell:
        "Rscript {input.script} > {log}"

################################################################################
# Prepare raw-input

rule C_prep:
    threads: 4
    input:
        script = 'scripts/C_prep.R',
        xs = [
            'raw-data/co2-adaptation-counts.tsv',
            'raw-data/PCC7002-genome.gff.gz'
        ]
    output:
        'data/C_meta.tsv',
        'data/C_raw-counts.tsv',
        'data/C_annotation.tsv',
    log: 'logs/C_prep.txt'
    shell:
        "Rscript {input.script} > {log}"
        
        
################################################################################
# Differential expression analysis

rule D_deg:
    input:
        script = 'scripts/D_deg.R',
        xs = [
            'scripts/helper_deg.R',
            'data/C_meta.tsv',
            'data/C_annotation.tsv',
            'data/C_raw-counts.tsv'
        ]
    output:
        'analysis/D_PCA.jpeg',
        'analysis/D_stagewise-adjusted-DEGs.tsv',
        'analysis/D_volcano.jpeg',
        'analysis/D_overview-by-type.tsv',
        'analysis/D_overview-by-logFC.tsv',
        'analysis/D_vst-expression.tsv',
        'analysis/D_normalized-counts.tsv',
        'analysis/D_heatmap.jpeg',
        'analysis/D_logFC-cor.jpeg'
    log: 'logs/D_deg.txt'
    shell:
        "Rscript {input.script} > {log}"
        
################################################################################
# export analysis results for visualization in stringApp

rule E_string:
    input:
        script = 'scripts/E_string.R',
        xs = [
            'data/C_annotation.tsv',
            'data/C_meta.tsv',
            'analysis/D_stagewise-adjusted-DEGs.tsv',
            'analysis/D_vst-expression.tsv'
        ]
    output:
        'analysis/E_string-cutoff.jpeg',
        'analysis/E_string-loci.txt',
        'analysis/E_string-z-expression.tsv'
    log: 'logs/E_string.txt'
    shell:
        "Rscript {input.script} > {log}"


################################################################################
################################################################################

rule all:
    input:
        'data/A_general-stat.tsv',
        'data/B_characteristics.tsv',
        'data/C_raw-counts.tsv',
        'analysis/D_stagewise-adjusted-DEGs.tsv',
        'analysis/E_string-loci.txt',
