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

* **`bash`** shell  for interacting with the computer

### Installation:

1. Most software dependecies are managed though **`conda`**, install as described at  <br> [https://conda.io/docs/install/quick.html](https://conda.io/docs/install/quick.html).
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


#### Compilation From Source

This tutorial requires **`Racon`** and **`Pinfish`** - instructions for their installation is included in the instructions in the **`installation`** guide above



### Usage: 

In your Conda environment, and in the tutorial working directory,

1. *optional* edit the provided **`config.yaml`** file to match your own study design
2. Run the Snakefile workflow (the command assumes 4 available threads; adjust to match your computer's capabilities)
```
    snakemake -j 4 all
```
3. Use the **`gffcompare`** tool to compare the annotated gene features to the previously downloaded **`gff`** file
```
gffcompare -r ReferenceData/Homo_sapiens.GRCh38.94.chromosome.20.gff3 -R -M Analysis/pinfish/polished_transcripts_collapsed.gff
```
4. . Render the report using results from the analysis above
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

This tutorial workflow will produce a rich description of your sequence library characteristics and the results from the Pinfish isoform analysis. Please visit the tutorial page at [https://community.nanoporetech.com/knowledge/bioinformatics]( https://community.nanoporetech.com/knowledge/bioinformatics) for more information

******************

# 4. Help:

### Licence and Copyright:

© 2019 Oxford Nanopore Technologies Ltd.

ont_tutorial_pinfish and other Nanopore tutorials are distributed by Oxford Nanopore Technologies under the terms of the MPL-2.0 license.

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

