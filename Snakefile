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
#P=/home/projects/rth/co2capture/subprojects/RNA-seq_air-to-30-pct-CO2
P=~/remote-server/projects/rth/co2capture/subprojects/RNA-seq_air-to-30-pct-CO2/

mkdir raw-data
cd raw-data

# Annotation
cp $P/PCC7002-genome.gff.gz                         PCC7002-genome.gff.gz
cp $P/PCC7002-genome.fna.gz                         PCC7002-genome.fna.gz

# Growth rate etc
cp $P/experimental-data/RNA_JAN23.xlsx .

# Expression matrices project data
cp $P/RNA-browser-Dataset2/analysis/41_counts.tsv   co2-adaptation-counts.tsv

# MultiQC stat
cp $P/RNA-browser-Dataset2/multiqc/multiqc_data/multiqc_data.json                   browser-data-multiqc.json
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
        'analysis/D_logFC-cor.jpeg',
        'analysis/D_overview.tsv',
        'analysis/D_overview.jpeg'
    log: 'logs/D_deg.txt'
    shell:
        "Rscript {input.script} > {log}"
        
################################################################################
# Create a beautiful treemap of the expressed genes

rule E_treemap:
    input:
        script = 'scripts/E_treemap.R',
        xs = [
            'data/C_annotation.tsv',
            'data/C_meta.tsv',
            'analysis/D_stagewise-adjusted-DEGs.tsv',
            'analysis/D_vst-expression.tsv'
        ]
    output:
        'analysis/E_zvst.tsv',
        'analysis/E_brite-hierarchy.tsv',
        'analysis/E_treemap.jpeg',
    log: 'logs/E_treemap.txt'
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

rule F_peptides:
    input:
        script = 'scripts/F_peptides.R',
        xs = [
            'data/C_annotation.tsv',
            'raw-data/co2-adaptation-counts.tsv',
            'raw-data/PCC7002-genome.fna.gz'
        ]
    output:
        'analysis/F_peptides.faa',
        'analysis/F_amino-frequencies.tsv',
        'analysis/F_amino-counts.tsv'
    log: 'logs/F_peptides.txt'
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
#   snakemake --cores all G_signalp

rule G_signalp:
    input:  
        'analysis/F_peptides.faa',
    output:
        directory('analysis/G_signalp')
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

rule H_process:
    input:
        script = 'scripts/H_signalp.R',
        xs = [
            'data/C_annotation.tsv',
            'analysis/D_stagewise-adjusted-DEGs.tsv',
            'data/C_meta.tsv',
            'analysis/D_vst-expression.tsv',
            'analysis/G_signalp'
            #'analysis/G_signalp/prediction_results.txt'
        ]
    output:
        'analysis/H_signal-probs.jpeg',
        'analysis/H_signals.tsv',
        'analysis/H_overview.tsv',
        'analysis/H_heatmap.jpeg',
        'analysis/H_signal-enrichment.tsv',
        'analysis/H_signal-enrichment.jpeg',
        'analysis/H_gene2pathway.tsv',
        'analysis/H_enriched-genes.tsv'
    log: 'logs/H_signalp.txt'
    shell:
        "Rscript {input.script} > {log}"
        


# add expression patterns for the custom carotenoid figure
rule H_carot:
    input:
        script = 'fig-carotenoid-pathway/add-expression.R',
        xs = [
          'data/C_meta.tsv',
          'analysis/D_vst-expression.tsv',
          'data/C_annotation.tsv',
          'fig-carotenoid-pathway/carotenoid.jpg',
        ]
    output:
        'fig-carotenoid-pathway/carotenoid-expression.jpg'
    log: 'logs/H_carotenoid.txt'
    shell:
        "Rscript {input.script} > {log}"
        
################################################################################
# inspect potential expression/AA correlation

rule I_aa:
    input:
        script = 'scripts/I_aa-composition.R',
        xs = [
            'data/C_annotation.tsv',
            'data/C_meta.tsv',
            'analysis/D_stagewise-adjusted-DEGs.tsv',
            'analysis/D_vst-expression.tsv',
            'analysis/F_amino-counts.tsv',
            'analysis/H_signals.tsv',
        ]
    output:
        'analysis/I_correlation-aa-length.jpeg',
        'analysis/I_length-log-normal.jpeg',
        'analysis/I_binomial-model.jpeg',
        'analysis/I_negative-binomial-model.jpeg',
        'analysis/I_aa-analysis.tsv',
        'analysis/I_volcano-like.jpeg',
        'analysis/I_extreme-AA-deg.jpeg',
        'analysis/I_extreme-AA.jpeg',
    log: 'logs/I_aa.txt'
    shell:
        "Rscript {input.script} > {log}"
        
        
        
################################################################################
# Differential expression analysis
# Focus: 30% CO2 vs all others

rule K_30focused:
    input:
        script = 'scripts/K_30-focused.R',
        xs = [
            'data/C_meta.tsv',
            'data/C_annotation.tsv',
            'data/C_raw-counts.tsv',
            'analysis/D_vst-expression.tsv',
            'analysis/D_stagewise-adjusted-DEGs.tsv'
        ]
    output:
        'analysis/K_logFC-vs-30.tsv',
        'analysis/K_logFC-comparisons.jpeg',
        'analysis/K_heatmap-30-focused.jpeg',
        'analysis/K_gsea.tsv',
        'analysis/K_string-30-focused.txt',
        'analysis/K_overview-30-focused.jpeg'
    log: 'logs/K_30-focused.txt'
    shell:
        "Rscript {input.script} > {log}"
        
################################################################################
# similar to F_clusters, but for the output for the M_30 network
## NOTE: Cytoscape with the network session loaded needs to be running for this script

rule M2_clusters:
    input:
        script = 'scripts/M2_cluster-overview.R',
        xs = [
            'data/C_annotation.tsv',
            # implicit dependency, don't trigger rule rerun for minor changes
            # eg of coloration or layout adjustments
            # string.cys
        ]
    output:
        'analysis/M2_mcl-clustering.tsv',
        'analysis/M2_cluster-enrichment-export-commands.txt',
        directory('analysis/M2_cluster-enrichment/'),
        'analysis/M2_string-enrichment.tsv',
        'analysis/M2_string-enrichment-genes.tsv',
    log: 'logs/M2_clusters.txt'
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
        'analysis/E_treemap.jpeg',
        'analysis/F_amino-counts.tsv'
        'analysis/H_signals.tsv',
        'fig-carotenoid-pathway/carotenoid-expression.jpg',
        'analysis/I_extreme-AA.jpeg',
        'analysis/K_overview-30-focused.jpeg'
        
        'analysis/M_logFC-vs-30.tsv',
        'analysis/M2_mcl-clustering.tsv',
