# Endocardium
Nkx2-5 Lineage-Derived Endocardium Contributes to Mesenchymal Cap Formation and Atrioventricular Septation

This repository contains the scripts used for upstream raw-data processing and downstream bioinformatics analysis in this study. The workflow is organized as follows:

1. Raw-Data Processing

RNAseq_rawdata_pre_processing.sh: Processes RNA-seq raw FASTQ files to generate the final gene expression matrix. It orchestrates the pipeline by calling the following two scripts:

Quantification_featurecounts.R: An R wrapper around featureCounts (Rsubread) that quantifies gene-level counts for a single sorted, deduplicated BAM. Called once per sample within the main script.

Merge_sample_quantification.py: Merges the per-sample count files into a single gene × sample expression matrix (with gene symbol annotation from the GTF) using the sample metadata. Executed once after all samples finish processing.

2. Bioinformatics_analysis.R performs a complete bulk RNA-seq differential expression and visualization workflow comparing GpTn and GpTp samples, including differential expression analysis with DESeq2, Gene Ontology biological process enrichment analysis of the DEGs using clusterProfiler and org.Mm.eg.db, customized volcano plot visualization of the DEGs, and a heatmap of selected genes.
