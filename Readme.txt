This code can be used to reproduce all figures and results from the manuscript "The Nuclear Receptor ERRg Cooperates with the Cardiogenic Factor GATA4 to
Orchestrate Cardiomyocyte Differentiation" by Sakamoto et al.

1. System requirements
    The software was tested in Ubuntu 20.04 with the following dependencies installed:
    
    Python 3.8.8 with the following pip packages:
        jupyterlab                         3.0.14
        tqdm                               4.59.0
        openpyxl                           3.0.7
        xerox                              0.4.1
        scipy                              1.6.2
        pandas                             1.2.4
        numpy                              1.20.1
        logomaker                          0.8
    R 4.1.1 with the following CRAN packages:
        openxlsx 4.2.4 
        reshape2 1.4.4 
        ggfittext 0.9.1 
        gridExtra 2.3 
        VennDiagram 1.6.20 
        pheatmap 1.0.12 
        tidyverse 1.3.1 
        rvest 1.0.1 
        xlsx 0.6.5 
        gsubfn 0.7 
        devtools 2.4.2 
        ashr 2.2-47 
        gplots 3.1.1 
    and the following Bioconductor packages:
        edgeR 3.34.1
        DESeq2 1.32.0
    Homer v4.11.1 with hg38 genome data
    salmon v1.3.0
    Bowtie 2 2.3.5.1
        
2. Installation guide
    There is nothing to install after the prerequisites are installed. The time to install all prerequisites is about 3 hours.
    
3. Demo
    There is no demo mode, the code will only run on the full data.
    
4. Instructions for use
    All the folders below are relative to the src/ folder.
    
    The reference data files should be placed in the reference/ folder:
        reference/Homo_sapiens.GRCh38.dna.toplevel.fa - reference genome
        reference/Homo_sapiens.GRCh38.cdna.all.fa.gz - reference transcriptome
        reference/Homo_sapiens.GRCh38.99.gtf - reference gene annotations
        
    The raw data should be placed in the following places:
        data/ATAC-Seq/hiPSCCM-ERRagKO/FGC2130 - ATAC-Seq .fastq.gz files
        data/ChIP-Seq/ERRg ChIP-seq in iPSC-CMs_Tomoya_2017/raw - ERRg ChIP-Seq .fastq.gz files, each in their own subfolder
        data/ChIP-Seq/H3K27Ac ChIP-seq in OE ERRg and KO ERRg/raw - H3K27Ac ChIP-Seq .fastq.gz files, each in their own subfolder
        data/RNA-Seq/ERRag_KO_hiPSC-CM/raw - RNA-Seq of ERRag KO, .fastq.gz files
        data/RNA-Seq/RNA-seq in OE ERRg iPSC-CM/raw - RNA-Seq of ERRaf overexpression experiment
        
    The code is split into 7 Jupyter notebooks, which should be opened in Jupyter lab:
        ERR_rnaseq_chipseq_atacseq.ipynb - processes the raw RNA-Seq, ChIP-Seq datasets, creates most figures
        ERRag_ATAC.ipynb - processes the ATAC-Seq dataset, makes some ATAC-related figures
        ERRag_ATAC_r.ipynb - makes ATAC fragment size distribution plots
        de_ko.ipynb - ERRag KO RNA-Seq processing with DESeq2
        de_oe.ipynb - ERRag overexpression RNA-Seq processing with DESeq2
        grab_data.ipynb - downloads published GATA4 data and does GATA4, GATA6, and superenhancer data processing
        QC.ipynb - quality control tables and plots for RNA-Seq, ChIP-Seq, and ATAC-Seq
        