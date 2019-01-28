# Tutorial Objectives

The **Pinfish tutorial** is intended as a functional guide to demonstrate how the Oxford Nanopore Technologies Pinfish software may be used to annotate genes, and their isoforms, against a reference genome sequence using long sequence read cDNA or direct-RNA data. A human cDNA dataset included with the tutorial demonstrates the key aspects of the workflow, and introduces **`gffcompare`** to enable the comparison of called datasets with already known genome annotations. This enables discovery and characterisation of novel isoforms. 

Sufficient information is provided in the tutorial such that the workflow can be tested, validated, and replicated. 

* What is the quality of the starting sequence collection?
* which fraction of the sequence reads appear full length?
* how do these full-length sequence reads map against the reference genome?
* visualisation of annotated genes and corresponding mapped sequence reads
* selection of novel and rare gene isoforms

**Methods used** in this tutorial include 

* **`R`** for statistical analysis and reporting
* **`conda`** for management of bioinformatics software installations
* **`minimap2`** for mapping sequence reads to reference genome
* **`pychopper`** for the filtering of candidate full-length transcript sequences
* **`pinfish`** for the clustering, polishing, and reduction of spliced sequence reads
* **`gffcompare`** for comparing gene annotations for novel isoform discovery
* **`igv`** for the visualisation of mapped reads and gene isoforms

**Computational requirements** for this tutorial include 

* Computer running Linux (Centos7, Ubuntu 18_10, Fedora 29) 
* At least 16Gb RAM (when using the human genome example data) 
* At least 15Gb spare disk space for analysis and indices
* Runtime with provided example data - approximately 1 hour


\pagebreak

# Quick start 

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
3. Change into the created `pinfish` sub-directory of the NanoporeTutorials
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
6. *optional* edit the provided **`config.yaml`** file to match your own study design
7. Run the Snakefile workflow (the command assumes 4 available threads; adjust to match your computer's capabilities)
```
    snakemake -j 4 all
```
8. Use the **`gffcompare`** tool to compare the annotated gene features to the previously downloaded **`gff`** file
```
gffcompare -r ReferenceData/Homo_sapiens.GRCh38.94.chromosome.20.gff3 -R -M Analysis/pinfish/polished_transcripts_collapsed.gff
```
9. . Render the report using results from the analysis above
```
    R --slave -e 'rmarkdown::render("Nanopore_Pinfish_Tutorial.Rmd", "html_document")'
```
\pagebreak




# Introduction

This tutorial aims to demonstrate a workflow for the de novo annotation and characterisation of genes and their isoforms on the basis of transcriptome sequencing. This tutorial is working with Oxford Nanopore cDNA sequences. 

A pre-prepared dataset is distributed with the tutorial. This contains a collection of cDNA sequences prepared from the GM12878 cell-line. Sequence reads were mapped against the human reference genome (https://github.com/nanopore-wgs-consortium/NA12878/blob/master/RNA.md) and reads from Chromosome 20 have been selected and are analysed within the described workflow.

Steps included within the tutorial include

* Selection of full-length transcripts using **`pychopper`**

**`Pinfish`** has its own pipeline automation package - [pipeline-pinfish-analysis](https://github.com/nanoporetech/pipeline-pinfish-analysis). This is a **`snakemake`** and **`conda`** based pipeline on which this tutorial is based. 



# Setup a computational environment for your QC analysis

This tutorial is intended to be simple to install, run, and customise. The analysis is performed using the **`R`** statistical software and further functionality is provided by a number of **`R packages`** and the **`RStudio`** graphical user interface. The following steps describe a simple approach to installing the tutorial and its dependencies.

## Conda package management software

**`Conda`** provides simple software package management. A wide variety of bioinformatics and scientific computing software has been deployed within Conda and it provides a streamlined way to install both software packages and their required dependencies without the requirement for administrative rights. Conda is available for Linux, macOS, and Windows - not all software packages are available for all platforms though. These installation instructions assume that you are using the BASH shell interface to your computer.

Install the Conda software on your computer using the instructions provided at [https://conda.io/docs/install/quick.html](https://conda.io/docs/install/quick.html)

A recommended **`Conda`** installation could be installed on Linux with the following commands

```
# download Python3 version of the Miniconda installer
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

bash Miniconda3-latest-Linux-x86_64.sh

bash
```

The commands here (1) download the Conda installer, (2) Install conda, and (3) reload your bash shell, and make the new Conda installation available.

It is recommended that the Conda installer be allowed to prepend the Conda path to your **`.bashrc`** file. This makes Conda available each time  you log in to the system. If Conda does not add this path, you should either edit your own `.bashrc` file or initialise your Conda environment manually with the instructions provided.

Check that Conda has been installed with the commands

```
echo $PATH
conda --version
```


## Other dependencies

The majority of functionality required for this tutorial is provided by the Conda environment. The github and subversion software are required to download the tutorial files. This tutorial requires an installation from source for the **`pychopper`**, **`pinfish`**, and **`racon`** software. These installations require 

* **`gcc`**
* **`cmake`**
* **`git-lfs`**


## Download the tutorial files

The tutorial documents and example data are contained on the [Oxford Nanopore Technologies Github site](www.github.com/nanoporetech). All tutorials are contained within a project called **`bioinformatics-tutorials`**. It is necessary to have **`git-lfs`** installed on your system to download the larger tutorial sequence and metadata files that are stored in Git Large File Storage. Please see [accompanying Git Large File Storage note](https://github.com/nanoporetech/bioinformatics-tutorials/git-lfs.md) for further information on how this may be best installed.

```
git clone --recursive https://github.com/nanoporetech/bioinformatics-tutorials.git NanoporeTutorials

# change into the new Pinfish working directory
cd NanoporeTutorials/pinfish
```

This command will download a collection of files into a new folder called `NanoporeTutorials`. The Pinfish tutorial is included in a subfolder called `pinfish`. The files downloaded into this sub-folder include

1. The Rmarkdown file, **`Nanopore_Pinfish_Tutorial.Rmd`**, includes this documentation and performs the analysis.
1. **`RawData/NA12878-cDNA-1D.chr20.filt.fastq`** is a fastq sequence file containing 282,507 sequence reads that map to human chromosome 20.
1. **`Static/Images`** is a folder containing some of the descriptive images used in the tutorial.
1. **`environment.yaml`** is a text file (in yaml format) that describes the software packages and computational environment that we would like to build using Conda.
1. **`config.yaml`** is a text file (in yaml format) that describes the Pinfish analysis that we hope to perform. This is the main file to be edited to describe your own analysis.
1. **`Snakefile`** is the workflow that describes how the Pinfish analysis can be performed in a largely automated manner.


## Build a conda environment

In the previous section we downloaded the project files. This download includes the file, `environment.yaml` that can be used to initialise a **`Conda`** working environment. To create a Conda environment arbitrarily named `BasicQC` we should use the commands

```
conda env create --name pinfish --file environment.yaml
source activate pinfish
```

## Download human reference genome and gene annotations for Chromsome 20.

The tutorial is provided with a fastq file containing sequence reads from the human genome that have been filtered for reads that best map to chromsome 20. For the further analysis of these sequence reads, we require both the human reference genome for re-mapping and we will use the gene annotation for chromosome 20 to evaluate the genes and isoform structures that are computed using the Pinfish pipeline. The commands below are used to download the required files and to uncompress them.

```
wget --directory-prefix ReferenceData ftp://ftp.ensembl.org/pub/release-94/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
  
wget --directory-prefix ReferenceData ftp://ftp.ensembl.org/pub/release-94/gff3/homo_sapiens/Homo_sapiens.GRCh38.94.chromosome.20.gff3.gz

gunzip -d ReferenceData/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz 

gunzip -d ReferenceData/Homo_sapiens.GRCh38.94.chromosome.20.gff3.gz

```

## Install **`pychopper`** software

[Pychopper](https://github.com/nanoporetech/pychopper) is a tool that can be used to identify full length cDNA sequence reads. Using the known adapter sequences that are used in the preparation of a sequencing library, the tool selects reads that contain both the 5' and 3' adapters and are thus more likely to represent full length transcripts. Since the primers used in library preparation are directional, the primers are also used to orientate the sequence reads so that they are all placed into the forward direction for subsequent sequence analysis.

```
pip install git+https://github.com/nanoporetech/pychopper.git
```

**`Pychopper`** is run on the provided sequence fastq file with the following command.

```
mkdir -p Analysis/PyChopper

cdna_classifier.py -b ReferenceData/cdna_barcodes.fas -r Analysis/PyChopper/NA12878-cDNA-1D.chr20.filt.pychopper.pdf ./RawData/NA12878-cDNA-1D.chr20.filt.fastq ./RawData/NA12878-cDNA-1D.chr20.filt.pychopper.fastq
```

This command runs the **`Pychopper`** method using the provided **`cdna_barcodes.fas`** fasta file. A summary pdf file is placed in the subfolder **`Analysis/PyChopper`** and the selected, and directionally oriented, sequence reads will be placed in a sequence file called **`RawData/NA12878-cDNA-1D.chr20.filt.pychopper.fastq`**. This file will form the input for the subsequent **`pinfish`** analysis


## Install **`Racon`**

[Racon](https://github.com/isovic/racon) is a fast consensus sequence based method that can be used to polish long sequence reads. The methods are computationally optimised to make the best use of various microprocessor hardware extensions. This means that binary versions of **`racon`** as provided through e.g. .deb / .rpm / bioconda may not work for all computer processors. For this reason we will install **``** racon from source.

```
git clone --recursive https://github.com/isovic/racon.git racon
mkdir racon/build && cd racon/build
cmake -DCMAKE_BUILD_TYPE=Release .. && make
cd ../..
export PATH=./racon/build/bin:$PATH
```

## Install **`Pinfish`**

[Pinfish](https://github.com/nanoporetech/pinfish) is the Oxford Nanopore Technologies tool for the annotation of genomes using long read transcriptomics data. To install **`pinfish`** for this tutorial it is sufficient to run the following command. 

```
git clone https://github.com/nanoporetech/pinfish.git
```


## Run the **`snakemake`** analysis of Snakefile

The provided **`config.yaml`** file provides a template for the sequence analysis. In the **`config.yaml`** file, there are fields for the specification of starting cDNA sequence fastq files, fields for the specification of the reference genome and a variety of other fields that can be used to optimise the mapping of sequence reads to the reference genome using **`minimap2`**, the logical requirements for the definition of isoform clusters. 

As you progress from running this tutorial with provided human sequence data to using your own cDNA sequences and own reference genome information, please edit the **`config.yaml`** file to match your experimental design. It is recommended to use [pipeline-pinfish-analysis](https://github.com/nanoporetech/pipeline-pinfish-analysis) for further automation and exploration of your libraries beyond this tutorial.

![](Static/Images/dag1.png) 

To run the **`snakemake`** 

```
# just type snakemake to run the workflow
# don't type <NPROC> but specify the number of processor cores available (e.g. 2 or 4)

snakemake -j <NPROC>
```

The **`-j <nprocs>`** flag instructs the **`snakemake`** process to use 8 threads. If you are running this tutorial on a laptop computer or an a large compute server, adjust this value to best correspond to your computer system.

## gffcompare 

**`gffcompare`** is a tool used to compare the content within and between genome annotation **`gff`** files. gffcompare can be run on the results of the **`pinfish`** workflow, and contrasted against the previously downloaded genome annotation using the command


```
gffcompare -r ReferenceData/Homo_sapiens.GRCh38.94.chromosome.20.gff3 -R -M Analysis/pinfish/polished_transcripts_collapsed.gff
```


## Render the tutorial report


The analysis of the sequences specified within the  **`config.yaml`** file and rendered through the **`Rmarkdown`** file will be performed as part of the **`knit`** process. This will load the results from the **`snakemake`** process, will prepare a sequence analysis, render figures and prepare the report. To start the analysis, it is only necessary to click the  **`knit`** button in the **`Rstudio`** software - please see figure \ref{fig:KnitIt}.


![](Static/Images/KnitIt.png)


It is also possible to perform the analysis using **`knit`** from the command line; the commands below will knit the document to produce respectively either an html document or a pdf report. 

```
R --slave -e 'rmarkdown::render("Nanopore_Pinfish_Tutorial.Rmd", "html_document")'
```

\pagebreak
