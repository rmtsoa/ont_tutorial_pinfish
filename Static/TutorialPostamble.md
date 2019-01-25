

# Customize the tutorial and explore your own data

Final thoughts; behind this **`Rmarkdown`** file (and its glossy pdf) is a modest amount of **`R code`** - please explore the **`Rmarkdown`** template; modify it, and run with your own samples.

The **`Nanopore_Pinfish_Tutorial.Rmd`** script contains the R code to perform the data presentation for the **`Pinfish`** analysis. This Rmarkdown script will first import the **`config.yaml`** configuration file to load the appropriate sequence files and references, and to identify where the **`snakemake`** process will have written result files. Please edit the **`config.yaml`** file to point to your own files and re-run the whole **`snakemake`** process. As a recommended best practise, place the raw fastq cDNA sequences files into the **`RawData`** folder within your project directory. Appropriate genome reference sequences should be placed in the **`ReferenceData`** folder.

The **`config.yaml`** and **`Nanopore_Pinfish_Tutorial.Rmd`** files can both be edited directly in Rstudio. The code below shows the command to open the configuration file with the **`Rstudio`** software.

```
rstudio config.yaml
```

The figure below shows a screenshot of the **`config.yaml`** file. Change the values in the fields **`inputFile`**, and **`flowcellId`** to reflect the details of your own sequencing run. Change the value of **`tutorialText`** to `FALSE` and the texts relating to the installation, configuration, and customisation of the tutorial will be hidden.

![](Static/Images/sumstatEditParams.png)<!--
```{r EditMethodParameters, echo=FALSE, include=TRUE, fig.margin=FALSE, fig.fullwidth = FALSE, fig.cap="Screenshot from top of the tutorial Rmarkdown script; the key lines to edit and modify are shown", cache=FALSE}
knitr::include_graphics("Static/Images/sumstatEditParams.png")
```
 -->


To extract the whole set of **`R code`** from the **`Rmarkdown`**, use the **`purl`** command - this will extract the R code into its own file.

```
knitr::purl("Nanopore_Pinfish_Tutorial.Rmd", quiet=TRUE)
```


\pagebreak



# References and Citations

