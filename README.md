![.](Static/Images/ONT_logo.png "Oxford Nanopore Technologies")

******************

# 1. Introduction:


### Overview:

The **Pinfish tutorial** is intended as a functional guide to demonstrate how the Oxford Nanopore Technologies Pinfish software may be used to annotate genes, and their isoforms, against a reference genome sequence using long sequence read cDNA or direct-RNA data. A human cDNA dataset included with the tutorial demonstrates the key aspects of the workflow, and introduces **`gffcompare`** to enable the comparison of called datasets with already known genome annotations. This enables discovery and characterisation of novel isoforms. 

### Features:

Sufficient information is provided in the tutorial such that the workflow can be tested, validated, and replicated. 

* What is the quality of the starting sequence collection?
* which fraction of the sequence reads appear full length?
* how do these full-length sequence reads map against the reference genome?
* visualisation of annotated genes and corresponding mapped sequence reads
* selection of novel and rare gene isoforms

******************

# 2. Getting Started:

This tutorial relies on **`Conda`** for the installion of the **`R`**, **`Rstudio`**, **`minimap2`**, and **`samtools`** software. It is necessary to have **`git-lfs`** installed on your system to download the larger tutorial sequence and metadata files that are stored in Git Large File Storage. Please see [accompanying Git Large File Storage note](https://github.com/nanoporetech/bioinformatics-tutorials/blob/master/git-lfs.md) for further information on how this may be best installed.

### Input and Output: 

This tutorial uses the code contained within the Github repository and an experimental design file (config.yaml) that processes a provided cDNA sequence file (in fastq format) with the downstream pinfish analytical workflow. 

### Dependencies:

This tutorial requires a computer running Linux (Centos7, Ubuntu 18_10, Fedora 29). >16Gb of memory would be recommended. The tutorial has been tested on minimal server installs of these operating systems.

Other dependencies include

* **`git`** packages for downloading the tutorial from Github repository
* **`git-lfs`**

### Installation:

1. Most software dependecies are managed though **`conda`**, install as described at  <br> [https://conda.io/docs/install/quick.html](https://conda.io/docs/install/quick.html).
```
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
    bash Miniconda3-latest-Linux-x86_64.sh
    bash
```
2. Download Nanopore tutorials & example files into folder named `NanoporeTutorials`. The tutorials are often using distributed with large sequence or metadata files and the [git-lfs](https://git-lfs.github.com/) extensions to git are used. Please ensure that **`git-lfs`** is first installed on your system.
```
    git clone --recursive https://github.com/nanoporetech/bioinformatics-tutorials.git NanoporeTutorials
```
3. Change into the created `pinfish` sub-directory of the NanoporeTutorials
```
    cd NanoporeTutorials/Pinfish
```
4. Install conda software dependencies with
```
    conda env create --name pinfish --file environment.yaml
```
5. Initialise conda environment with 
```
    source activate pinfish
```
6. Install **`pychopper`** software
```
    pip install git+https://github.com/nanoporetech/pychopper.git
```
7. Install **`Racon`** software, and add to the path. Although **`Racon`** is installed using **`Bioconda`**, the Conda provided version may utilise SIMD instructions that are unavailable to the processor on the computer system being used. Installation from source is strongly advised.
```
   git clone --recursive https://github.com/isovic/racon.git racon
   mkdir racon/build && cd racon/build
   cmake -DCMAKE_BUILD_TYPE=Release .. && make
   cd ../..
   export PATH=./racon/build/bin:$PATH
```
8. Install the **`pinfish`** software
```
   git clone https://github.com/nanoporetech/pinfish.git
```

#### Compilation From Source

This tutorial requires **`Racon`** and **`Pinfish`** - instructions for their installation is included in the instructions in the **`installation`** guide above



### Usage: 

In your Conda environment, and in the tutorial working directory,

1. Download the reference genome (and unzip it)
```
    wget --directory-prefix ReferenceData ftp://ftp.ensembl.org/pub/release-94/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
    gunzip -d ReferenceData/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
    wget --directory-prefix ReferenceData ftp://ftp.ensembl.org/pub/release-94/gff3/homo_sapiens/Homo_sapiens.GRCh38.94.chromosome.20.gff3.gz
    gunzip -d ReferenceData/Homo_sapiens.GRCh38.94.chromosome.20.gff3.gz
```
2. Select full length transcripts and standardise read orientation
```
    mkdir -p Analysis/PyChopper
    cdna_classifier.py -b ReferenceData/cdna_barcodes.fas -r Analysis/PyChopper/NA12878-cDNA-1D.chr20.filt.pychopper.pdf ./RawData/NA12878-cDNA-1D.chr20.filt.fastq ./RawData/NA12878-cDNA-1D.chr20.filt.pychopper.fastq
```
3. *optional* edit the provided **`config.yaml`** file to match your own study design
4. Run the Snakefile workflow (the command assumes 8 available threads; adjust to match your computer's capabilities)
```
    snakemake -j 8 all
```
5. Use the **`gffcompare`** tool to compare the annotated gene features to the previously downloaded **`gff`** file
```
gffcompare -r ReferenceData/Homo_sapiens.GRCh38.94.chromosome.20.gff3 -R -M Analysis/pinfish/polished_transcripts_collapsed.gff
```
6. Render the report using results from the analysis above
```
    R --slave -e 'rmarkdown::render("Nanopore_Pinfish_Tutorial.Rmd", "html_document")'
```

The provided Rmarkdown tutorial script can also be opened directly in Rstudio

```
rstudio Nanopore_cDNA_Tutorial.Rmd
```

The report can be prepared by "knit" from the GUI as shown in the figure

![.](Static/Images/KnitIt.png "Prepare a report using Knit")


******************

# 3. Results

This tutorial workflow will produce a rich description of your sequence library characteristics and the results from the differential expression analysis. Please visit the tutorial page at [https://community.nanoporetech.com/knowledge/bioinformatics]( https://community.nanoporetech.com/knowledge/bioinformatics) for more information

******************

# 4. Help:

### Licence and Copyright:

© 2018 Oxford Nanopore Technologies Ltd.

Bioinformatics-Tutorials is distributed by Oxford Nanopore Technologies under the terms of the MPL-2.0 license.

### FAQs:



### Abbreviations:


* __knit__ is the command to render an Rmarkdown file. The knitr package is used to embed code, the results of R analyses and their figures within the typeset text from the document. 

* __L50__  the number of sequences (or contigs etc) that are longer than, or equal to, the N50 length and therefore include half the bases of the assembly

* __N50__  length such that sequences (or contigs etc) of this length or longer include half the bases of the sequence collection

* __Rmarkdown__ is an extension to markdown. Functional R code can be embedded in a plain-text document and subsequently rendered to other formats including the PDF format of this report.

* __QV__  the quality value - -log10(p) that any given base is incorrect. QV may be either at the individual base level, or may be averaged across whole sequences

* __sequencing_summary.txt__ a summary file describing sequence characteristics following base calling with the Guppy / Albacore software.


### References and Supporting Information:

*  https://community.nanoporetech.com/knowledge/bioinformatics
*  https://www.r-project.org/
*  https://snakemake.readthedocs.io/en/stable/
*  https://bioconda.github.io/

