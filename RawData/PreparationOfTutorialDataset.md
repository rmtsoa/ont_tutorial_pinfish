## Preparation of training set for pinfish evaluation

This min-vignette describes the preparation of a restricted starting dataset for the Pinfish tutorial. The public [NA12878 Nanopore Consortium RNA datasets](https://github.com/nanopore-wgs-consortium/NA12878/blob/master/RNA.md) are used as primary datasource. The RNA dataset includes both direct-RNA and cDNA sequencing from the GM12878 cell-line. For purposes of isoform detection and a demonstration of Pinfish workflow concepts, this dataset appears appropriate

#### Software requirements

* curl
* samtools
* picard

#### diskspace requirements

The downloaded files, indexes, and derived fastq require 20Gb of diskspace

#### time requirements

Using 100Mbit download connection and an aged Xeon server, total time required = 1 hour

### Recommended workflow

Download the reads mapped against the human reference as BAM and corresponding index file

```
curl -O http://s3.amazonaws.com/nanopore-human-wgs/rna/bamFiles/NA12878-cDNA-1D.pass.dedup.fastq.hg38.minimap2.sorted.bam
curl -O http://s3.amazonaws.com/nanopore-human-wgs/rna/bamFiles/NA12878-cDNA-1D.pass.dedup.fastq.hg38.minimap2.sorted.bam.bai
```

Filter the reads that have been mapped to chromosome 20. It may be worth having a quick look at the BAM headers `samtools view -H  NA12878-cDNA-1D.pass.dedup.fastq.hg38.minimap2.sorted.bam` to check the names of the chromosomes used in the reference that reads were mapped against

```
samtools view -b -h NA12878-cDNA-1D.pass.dedup.fastq.hg38.minimap2.sorted.bam "chr20" > NA12878-cDNA-1D.chr20.filt.bam
```

Extract the reads into a fastq file for subsequent analysis

```
picard SamToFastq I=NA12878-cDNA-1D.chr20.filt.bam FASTQ=NA12878-cDNA-1D.chr20.filt.fastq INCLUDE_NON_PRIMARY_ALIGNMENTS=false
```

The final file; `NA12878-cDNA-1D.chr20.filt.fastq`, is the required input for the Pinfish tutorial