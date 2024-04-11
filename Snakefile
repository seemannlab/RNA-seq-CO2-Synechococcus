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
# NOTE: Run Cytoscape + stringApp and MCL cluster for the output of this rule

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
# extract clustering results from cytoscape
## NOTE: Cytoscape with the network session loaded needs to be running for this script

rule F_clusters:
    input:
        script = 'scripts/F_cluster-overview.R',
        xs = [
            'data/C_annotation.tsv',
            # implicit dependency, don't trigger rule rerun for minor changes
            # eg of coloration or layout adjustments
            # string.cys
        ]
    output:
        'analysis/F_mcl-clustering.tsv',
        'analysis/F_cluster-enrichment-export-commands.txt',
        directory('analysis/F_cluster-enrichment/'),
        'analysis/F_string-enrichment.tsv',
        'analysis/F_string-enrichment-genes.tsv',
    log: 'logs/F_clusters.txt'
    shell:
        "Rscript {input.script} > {log}"

        
################################################################################
################################################################################
# Investigate amino acid composition

rule G_peptides:
    input:
        script = 'scripts/G_peptides.R',
        xs = [
            'data/C_annotation.tsv',
            'raw-data/co2-adaptation-counts.tsv',
            'raw-data/PCC7002-genome.fna.gz'
        ]
    output:
        'analysis/G_peptides.faa',
        'analysis/G_frequencies.tsv',
    log: 'logs/G_peptides.txt'
    shell:
        "Rscript {input.script} > {log}"


################################################################################
################################################################################
# SignalP
# Prediction of signal peptides for all coding genes in the genome
# https://services.healthtech.dtu.dk/services/SignalP-6.0/
# Challange: The software is only available for academic use, and may not be
# freely shared, therefore this rule needs some additional care.
#
# Suggestion:
# - Create a conda/mamba environment with the signalp6 placeholder
#   mamba create -n signalp6 -c predector signalp6
# - Activate environment and keep stacked to still use snakemake from a base env
#   mamba activate --stack signalp6
# - Manually download the software from the SignalP website
# - "Incooperate" the software into the environment
#   signalp6-register signalp-6.0*.fast.tar.gz
# - Run snakemake for this rule
#   snakemake --cores all H_signalp

rule H_signalp:
    input:  
        'analysis/G_peptides.faa',
    output:
        directory('analysis/H_signalp')
    threads: 8
    shell:
        """
        signalp6 --organism other         \
            --mode fast --format all      \
            --fastafile {input}           \
            --output_dir {output}
        # somehow this parameter seems to not work, but 8 is supposedly default
        # --torch_num_threads
        """

################################################################################
# process the SignalP predictions

rule I_process:
    input:
        script = 'scripts/I_signalp.R',
        xs = [
            'data/C_annotation.tsv',
            'analysis/D_stagewise-adjusted-DEGs.tsv',
            'data/C_meta.tsv',
            'analysis/D_vst-expression.tsv',
            'analysis/H_signalp'
            #'analysis/H_signalp/prediction_results.txt'
        ]
    output:
        'analysis/I_signal-probs.jpeg',
        'analysis/I_signals.tsv',
        'analysis/I_overview.tsv',
        'analysis/I_heatmap.jpeg',
        'analysis/I_signal-enrichment.tsv',
        'analysis/I_signal-enrichment.jpeg',
        'analysis/I_gene2pathway.tsv',
        'analysis/I_enriched-genes.tsv'
    log: 'logs/I_signalp.txt'
    shell:
        "Rscript {input.script} > {log}"
        
        
################################################################################
# inspect potential expression/AA correlation

rule J_AA:
    input:
        script = 'scripts/J_aa-freqs.R',
        xs = [
            'data/C_annotation.tsv',
            'analysis/D_stagewise-adjusted-DEGs.tsv',
            'analysis/D_vst-expression.tsv',
            'data/C_meta.tsv',
            'analysis/G_frequencies.tsv',
        ]
    output:
        'analysis/J_freqs-overall.jpeg',
        'analysis/J_AA-expr-cor.jpeg',
    log: 'logs/J_aa.txt'
    shell:
        "Rscript {input.script} > {log}"
        
################################################################################
# Complement with GSEA plots

rule K_gsea:
    input:
        script = 'scripts/K_gsea.R',
        xs = [
            'data/C_annotation.tsv',
            'data/C_meta.tsv',
            'analysis/D_stagewise-adjusted-DEGs.tsv',
            'analysis/D_normalized-counts.tsv',
            'analysis/I_gene2pathway.tsv'
        ]
    output:
        'analysis/K_gsea-0.04-.30..CO2.jpeg',
        'analysis/K_gsea-0.04-.4..CO2.jpeg',
        'analysis/K_gsea-0.04-.8..CO2.jpeg',
        'analysis/K_gsea-4-.30..CO2.jpeg',
        'analysis/K_gsea-4-.8..CO2.jpeg',
        'analysis/K_gsea-8-.30..CO2.jpeg',
        'analysis/K_gsea.tsv',
        'analysis/K_enrichment-overview.jpeg',
        'analysis/K_expression-overview.jpeg',
    log: 'logs/K_gsea.txt'
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
        'analysis/F_mcl-clustering.tsv',
        'analysis/G_peptides.faa',
        'analysis/J_freqs-overall.jpeg',
        'analysis/K_gsea.tsv',
