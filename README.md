# hga_rnaseq
Scripts and reference files used for the RNAseq read processing, analysis, and graphing in [link]
Raw sequencing reads are available in the SRA PRJNA943215

1. rnaseq_pipe.sh - Command line bash script for raw sequencing read quality control and transcript mapping to the reference genome
2. Lp_genbank.fasta - Reference file used for Kallisto transcript mapping [Genbank accession GCA_000008485.1]
3. hga_rnaseq.Rmd - RStudio Markdown file for differential gene expression analysis and visualization
4. lqsr_rnaseq_samples.txt - Sample key to convert from kallisto abundance output to replicate condition
5. tx2gene_Lp.csv - Comma delimited file for custom timxport DESeq2 annotation object
