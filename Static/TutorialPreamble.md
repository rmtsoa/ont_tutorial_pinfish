# Tutorial objectives

The **Pinfish tutorial** is intended as a functional guide to demonstrate how the Oxford Nanopore Technologies Pinfish software may be used to annotate genes, and their isoforms, against a reference genome sequence using long sequence read cDNA or direct-RNA data. A human cDNA dataset included with the tutorial demonstrates the key aspects of the workflow, and introduces **`gffcompare`** to enable the comparison of identified genes and gene isoforms with known gene annotations. This workflow thus enables the discovery and characterisation of genes and putatively novel isoforms. 

The tutorial is provided with example starting sequence data so that the workflow can be tested to address experimental questions such as 

* What is the quality of the starting sequence collection?
* which fraction of the sequence reads appear full length?
* how well do these full-length sequence reads map against the reference genome?
* visualisation of annotated genes and corresponding mapped sequence reads
* selection of novel and rare gene isoforms

Editing of the workflow's configuration file, **`config.yaml`** will allow the workflow to be run with different starting cDNA sequence libraries, reference genomes, and gene annotations.

## Methods used include 

* **`R`** for statistical analysis and reporting
* **`conda`** for management of bioinformatics software installations
* **`snakemake`** for managing the bioinformatics workflow
* **`pychopper`** for the filtering of candidate full-length transcript sequences
* **`minimap2`** for mapping sequence reads to reference genome
* **`samtools`** for SAM/BAM handling and mapping statistics
* **`pinfish`** for the clustering, polishing, and reduction of spliced sequence reads
* **`racon`** is used by **`pinfish`** for the polishing of isoform consensus sequences
* **`gffcompare`** for comparing gene annotations for novel isoform discovery
* **`igv`** for the visualisation of mapped reads and gene isoforms

## The computational requirements include 

* Computer running Linux (Centos7, Ubuntu 18_10, Fedora 29) 
* At least 16Gb RAM (when using the human genome example data) 
* At least 15Gb spare disk space for analysis and indices
* Runtime with provided example data - approximately 1 hour


\pagebreak

# Software installation

1. Most software dependecies are managed though **`conda`**, install as described at  <br> [https://conda.io/docs/install/quick.html](https://conda.io/docs/install/quick.html). Accept the license agreement during installation, and it is recommended to allow the Conda installer to prepend its path to your `.bashrc` file when asked.
```
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
    bash Miniconda3-latest-Linux-x86_64.sh
    bash
```
2. Download Nanopore tutorial & example files into folder named `Pinfish`. This tutorial requires the **`git-lfs`** large file support capabilities; this should be installed first through **`Conda`**
```
    conda install -c conda-forge git-lfs
    git lfs install
    git clone https://github.com/nanoporetech/ont_tutorial_pinfish.git Pinfish
```
3. Change into the `Pinfish` sub-directory created
```
    cd Pinfish
```
4. Install conda software dependencies with
```
    conda env create --name Pinfish --file environment.yaml
```
5. Initialise conda environment with 
```
    source activate Pinfish
```


\pagebreak




# Introduction

cDNA and direct-RNA sequencing enables the preparation of full-length sequence transcripts from a biological sample. Collections of sequences may be used for both differential gene expression analyses and for the characterisation and annotation of genes and their isoforms. In this tutorial we will look at a workflow that aims to collapse full-length cDNA sequence reads into clusters that can be polished and used to discover and annotate genes.

[Oxford Nanopore Technologies](https://nanoporetech.com) provides a number of [sequencing solutions](https://nanoporetech.com/rna) to allow users to generate a "snapshot"" of gene expression. This can be achieved by both sequencing the mRNA [directly](https://store.nanoporetech.com/catalog/product/view/id/167/s/direct-rna-sequencing-kit/category/28/), or via a complementary DNA ([cDNA](https://store.nanoporetech.com/catalog/product/view/id/177/s/cdna-pcr/category/28/)) proxy. In contrast to short read sequencing technologies, entire mRNA transcripts can be captured as single reads. 

Once sequencing data has been produced from the biological sample(s), the sequence reads can be mapped to the host's reference genome. Using a gapped alignment strategy, it is anticipated that individual sequence reads will map only to a gene's exons. Gene regions can then be identified by filtering the sequence alignments for depth-of-coverage and clustering can be used to define the distinct isoforms.

This tutorial demonstrates a workflow for the de novo annotation and characterisation of genes and their isoforms from collections of cDNA sequence. This is acheived using the Pinfish software. Additional steps are included using the **`IGV`** genome browser to visualise the gene annotations and the **`gffcompare`** software is used to select for known and potentially novel transcripts (when a reference annotation is available)

A pre-prepared dataset is distributed with the tutorial. This contains a collection of cDNA sequences sampled from the publicly available GM12878 cell-line transcriptome (https://github.com/nanopore-wgs-consortium/NA12878/blob/master/RNA.md). Using the provided transcriptome mapping BAM file, cDNA reads from Chromosome 20 have been selected and are analysed within the described workflow.

The goals included for this tutorial include:

* To introduce a literate framework for analysing Oxford Nanopore cDNA data prepared using MinION or PromethION flowcells
* To utilise best data management practices
* To provide basic cDNA sequence QC metrics such that review and consideration of the starting experimental data can be performed
* To use the **`Pychopper`** software to select for the subset of cDNA sequence reads that are most likely to be full-length
* To use the **`Pinfish`** software to annotate genes and their isoforms against the provided reference genome
* To visualise patterns of genes and their isoforms using **`IGV`** and to super-impose these data with raw sequence reads
* To use the **`gffcompare`** software to select for known and novel isoforms within the gene annotation


# Getting started and best practices

This tutorial requires a computer workstation running a Linux operating system. The workflow described has been tested using **`Fedora 29`**, **`Centos 7`**, and **`Ubuntu 18_04`**. This tutorial has been prepared in the **`Rmarkdown`** file format. This utilises *markdown* (an easy-to-write plain text format as used in many Wiki systems) - see @R-rmarkdown for more information about **`rmarkdown`**. The document template contains chunks of embedded **`R code`** that are dynamically executed during the report preparation. 

The described analytical workflow makes extensive use of the **`conda`** package management and the **`snakemake`** workflow software. These software packages and the functionality of **`Rmarkdown`** provide the source for a rich, reproducible, and extensible tutorial document.

The workflow contained within this Tutorial performs a real bioinformatics analysis and uses the whole human genome as an example. There are some considerations in terms of memory and processor requirement. Indexing the whole human genome for sequence read mapping using **`minimap2`** (performed within **`Pinfish`**) will use at least **`18Gb`** of memory. The minimal recommended hardware setup for this tutorial is therefore an 4 threaded computer with at least 16Gb of RAM and 15Gb of storage space. 

There are few dependencies that need to be installed at the system level prior to running the tutorial. The **`conda`** package management software will coordinate the installation of the required bioinformatics software and their dependencies in user space - this is dependent on a robust internet connection.

As a best practice this tutorial will separate primary cDNA sequence data (the base-called fastq files) from the **`Rmarkdown`** source, and the genome reference data. The analysis results and figures will again be placed in a separate working directory. The required layout for the primary data is shown in figure \ref{fig:FolderLayout}. This minimal structure will be prepared over the next few sections of this tutorial. The cDNA sequences to be analysed must be placed within the folder called **`RawData`**, the reference genome and annotation files must be placed in the folder named **`ReferenceData`**.

![](Static/Images/FolderLayout.png) 

The figure above shows an example of the folder layout after the tutorial has been cloned from `github` and reference genome data has been downloaded. 

* `config.yaml` is the file that defines the tutorial configuration - this includes sequence files, reference genomes and parameters to be used for the analyses
* `environment.yaml` defines the files that should be downloaded and installed the Conda environment required for the data analysis
* `Nanopore_Pinfish_Tutorial.Rmd` is the **`Rmarkdown template`** that will be run to prepare the analysis report
* `RawData` is a folder and your sequence files should be placed only in this folder. In the figure, a bzip2 compressed fastq file is described; this is the dataset provided with the tutorial
* `ReferenceData` is a folder that contains genome references; this may be public or proprietary data but these data are used as a reference and are considered immutable within the analysis The folder contains a gzipped fasta format genome reference and gzipped GFF3 genome annotation
* `Snakefile` is the workflow that will be used with the **`snakemake`** command to perform the data analysis
* `Static` is a folder containing various static files; this includes figures used within the tutorial, citations and the text used for the report


# Experimental setup

![](Static/Images/ConfigYamlParameters.png) 

The experimental design for the tutorial is described in a separate file called **`config.yaml`** - an example is provided with the tutorial and the head of the file is shown in the figure above.

* **`genome_fasta`** defines the genome reference to which the sequence reads will be mapped. This can be either a URL to the file at e.g. ENSEMBL/NCBI or can be the name of the file already placed in the `ReferenceData` folder.
* **`genome_annot`** defines the prior gene annotations for the genome sequence described by `genome_fasta`. These annotations will be used for classification of known and potentially novel transcripts. This can be either a URL or a file name for a file already placed in the `ReferenceData` folder.
* **`raw_fastq`** defines the input sequences for the Pinfish workflow. This should be a fastq file, placed in the `RawData` folder that can optionally be gzip or bzip2 compressed.
* **`pychopper`** is a boolean option (`TRUE` or `FALSE`) that defines if **`pychopper`** should be run.
* **`minimum_cluster_size`** defines the number of sequence reads that must appear within a sequence cluster for it to be considered for further annotation and analysis
* **`minimum_mapping_quality`** defines a threshold for the consideration of mapped sequence reads.

Other parameters that can modified include parameters for controlling read mapping to the reference genome, expectations for intron and exon boundaries and rules for collapsing annotations. If you run the analysis without the **`pychopper`** sequence filtering it is recommended to modify the **`minimap2`** parameters which as provided force the method to consider only the forward transcript strand.


## Example dataset

This tutorial is provided with an example dataset that demonstrates the objectives and capabilities of the workflow. The provided dataset is derived from the public domain. The `Nanopore WGS consortium` [citation] released human transcriptome data prepared from a cell line (GM12878) corresponding to the human reference genome (NA12878) https://github.com/nanopore-wgs-consortium/NA12878/blob/master/RNA.md. These publicly provided data were released in November 2017 and were produced from 12 MinION flowcells of cDNA sequencing. The provided BAM file of cDNA read mappings to the reference genome was used as a data source and the sequence reads mapping to Chromosome 20 have been extracted for the subset used in this tutorial. 



# Snakemake

**`Snakemake`** is a workflow management system [citation].


![](Static/Images/dag.png) 


The directed acyclic graph shown above describes the **`Snakemake`** defined workflow associated with this tutorial. The core tasks involved with mapping sequence reads and the preparation of transcript clusters and their annotations have been utilised verbatim from the **`pipeline-pinfish-analysis`** https://github.com/nanoporetech/pipeline-pinfish-analysis workflow.

The **`Snakemake`** workflow imports all of the configuration parameters from the **`config.yaml`** file that was described in the previous section. The `Snakefile` describes the analytical workflow and the locations of the required input and output files. To modify and adapt the workflow, it is recommended to first modify the `config.yaml` file.


The **`Snakemake`** tasks can be summarised as

1. Download (and decompress) the reference genome sequence against which sequence reads will be mapped.
2. Download (and decompress) the reference genome annotations to be used for the analysis.
3. Unpack (or decompress) the `fastq` format sequence files to be used for the analysis
4. Use the **`Pychopper`** software to select for likely full-length transcript sequences and orientate sequence reads so that they are directionally uniform.
5. Create a genome index for read mapping using the **`Minimap2`** software
6. Map the Pychopper oriented sequence reads against the reference genome using Minimap2
7. Prepare a GFF annotation from the mapped read BAM file using Pinfish spliced_bam2gff
8. Correct sequence reads (using **`Racon`**) within transcript clusters with the Pinfish polish_clusters tool
9. Remap the polished and transcript clustered sequence reads against the reference genome using Minimap2 
10. Prepare a GFF annotation from the polished and clustered mapped transcript BAM file using Pinfish spliced_bam2gff
11. Collapse redundant and partial transcript annotations using Pinfish collapse_partials  
12. Prepare a fasta file of non-redundant, polished and collapsed transcripts using **`gffread`**
13. Compare the de novo produced transcript annotations with the reference genome annotations using **`GFFcompare`**. This is used to assess sensitivity and selectivity and can be used to identify apparently novel transcripts and gene isoforms.


To run the **`snakemake`**  workflow

```
# just type snakemake to run the workflow
# don't type <NPROC> but specify the number of processor cores available (e.g. 2 or 4)

snakemake -j <NPROC>
```

The **`-j <nprocs>`** flag instructs the **`snakemake`** process to use a specific number of threads. This value should be adjusted to best correspond to your computer system.



## Render the tutorial report


The analysis of the sequences specified within the  **`config.yaml`** file and rendered through the **`Rmarkdown`** file will be performed as part of the **`knit`** process. This will load the results from the **`snakemake`** process, will prepare a sequence analysis, render figures and prepare the report. To start the analysis, it is only necessary to click the  **`knit`** button in the **`Rstudio`** software - please see figure \ref{fig:KnitIt}.


![](Static/Images/KnitIt.png)


It is also possible to perform the analysis using **`knit`** from the command line; the commands below will knit the document to produce respectively either an html document or a pdf report. 

```
R --slave -e 'rmarkdown::render("Nanopore_Pinfish_Tutorial.Rmd", "html_document")'
```

\pagebreak
