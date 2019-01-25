# Pinfish tutorial Snakefile

import os
import re
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider

HTTP = HTTPRemoteProvider()

configfile: "config.yaml"

ProvidedGenomeFastaLink = config["genome_fasta"]
ReferenceData = config["ReferenceData"]
RawData = config["RawData"]

# deprecated #SNAKEDIR = os.path.dirname(workflow.snakefile)


# Validate that the fasta link is on filesystem; if URL then download
# parse on leading https / http / ftp ...
webLink = ""
if (re.search("ftp://", ProvidedGenomeFastaLink) or re.search("http://", ProvidedGenomeFastaLink) or re.search("https://", ProvidedGenomeFastaLink)):
  webLink = re.sub("^[^:]+://", "", ProvidedGenomeFastaLink)
  ProvidedGenomeFasta = os.path.join(ReferenceData, os.path.basename(webLink)) 
else:
  ProvidedGenomeFasta = os.path.join(ReferenceData, ProvidedGenomeFastaLink)
  webLink = "dummy.domain.com/"+ProvidedGenomeFasta

# some methods require a gunzipped fasta file ... 
UnpackedReferenceFasta = re.sub("\.gz$","",ProvidedGenomeFasta)

if (UnpackedReferenceFasta == ProvidedGenomeFasta):
  # workflow has been called with an uncompressed starting fasta file ...
  # This means that cannot unpack ... - and there is a circular dependency ...
  ProvidedGenomeFasta = os.path.join(os.path.dirname(workflow.snakefile), "Snakefile")
  # this is a nonsense entry; but the judicious usage of ancient means should not be called


# Definition for the Minimap2 index
ReferenceIndex = os.path.join("Analysis/Minimap2", os.path.basename(UnpackedReferenceFasta) + ".mmi")


RawFastq = os.path.join(RawData, config["raw_fastq"])
UnpackedRawFastq = RawFastq
if (re.search("\.gz$",RawFastq)):
  UnpackedRawFastq = re.sub("\.gz$", "", RawFastq)
elif (re.search("\.bz2$", RawFastq)):
  UnpackedRawFastq = re.sub("\.bz2$","",RawFastq)


PychopperPDF = os.path.join("Analysis/Pychopper", "pychopper.pdf") 
PychopperFastq = os.path.join("Analysis/Pychopper", re.sub(os.path.splitext(UnpackedRawFastq)[1],"",os.path.basename(UnpackedRawFastq))+".pychopper.fastq") 


print(PychopperPDF)
print(PychopperFastq)

Minimap2BAM = os.path.join("Analysis/Minimap2", re.sub(os.path.splitext(PychopperFastq)[1],"",os.path.basename(PychopperFastq))+".bam")

PinfishRawTranscripts = os.path.join("Analysis/Pinfish" , "raw_transcripts.gff")
PinfishClusteredTranscripts = os.path.join("Analysis/Pinfish" ,"clustered_transcripts.gff")
PinfishClusterMemberships = os.path.join("Analysis/Pinfish" ,"cluster_memberships.tsv")
PinfishClusteredCollapsedTranscripts = os.path.join("Analysis/Pinfish" ,"clustered_transcripts_collapsed.gff")
PinfishPolishedTranscriptFasta = os.path.join("Analysis/Pinfish" ,"polished_transcripts.fas")
MinimapPolishedBam = os.path.join("Analysis/Minimap2", "polished_reads_aln_sorted.bam")
PinfishPolishedTranscripts = os.path.join("Analysis/Pinfish" ,"polished_transcripts.gff")
PinfishPolishedCollapsedTranscripts = os.path.join("Analysis/Pinfish" ,"polished_transcripts_collapsed.gff")
PinfishCorrectedTranscriptome = os.path.join("Analysis/Pinfish" ,"corrected_transcriptome_polished_collapsed.fas")

rule all:
  input:
    PinfishCorrectedTranscriptome,
    PinfishClusteredCollapsedTranscripts
    

# download the reference genome sequence
rule DownloadReferenceGenome:
  input:
    ancient(HTTP.remote(webLink, keep_local=True))
  output:
    ancient(ProvidedGenomeFasta)
  run:
    shell("mv {input} {output}")
  

# gunzip the provided genome reference if it is gzipped ...
rule UnpackReferenceGenome:
  input:
    ancient(ProvidedGenomeFasta)
  output:
    ancient(UnpackedReferenceFasta)
  shell:
    "gunzip --keep -d {input}"
    
    
# gunzip the provided rawSequence if it is gzipped or bzip2 compressed ...
rule UnpackFastqData:
  input:
    RawFastq
  output:
    ancient(UnpackedRawFastq)
  run:
    if (re.search("\.gz$", RawFastq)):
      shell("gunzip --keep -d {input}")
    elif (re.search("\.bz2$", RawFastq)):
      shell("bunzip2 --keep -d {input}")


# use pychopper to score for the full length sequence reads
rule Pychopper:
  input:
    UnpackedRawFastq
  output:
    pdf = ancient(PychopperPDF),
    fastq = ancient(PychopperFastq)
  run:
    shell("cdna_classifier.py -b ReferenceData/cdna_barcodes.fas -r {output.pdf} {input} {output.fastq}")




# build minimap2 index
rule Minimap2Index: 
    input:
        genome = UnpackedReferenceFasta
    output:
        index = ReferenceIndex
    params:
        opts = config["minimap_index_opts"]
    threads: config["threads"]
    shell:"""
        minimap2 -t {threads} {params.opts} -I 1000G -d {output.index} {input.genome}
    """


rule Minimap2: ## map reads using minimap2
    input:
       index = rules.Minimap2Index.output.index,
       fastq = PychopperFastq if (config["pychopper"]==True) else UnpackedRawFastq
    output:
       bam = Minimap2BAM
    params:
        opts = config["minimap2_opts"],
        min_mq = config["minimum_mapping_quality"],
    threads: config["threads"]
    shell:"""
    minimap2 -t {threads} -ax splice {params.opts} {input.index} {input.fastq}\
    | samtools view -q {params.min_mq} -F 2304 -Sb | samtools sort -@ {threads} - -o {output.bam};
    samtools index {output.bam}
    """


rule PinfishRawBAM2GFF: ## convert BAM to GFF
    input:
        bam = rules.Minimap2.output.bam
    output:
        raw_gff = PinfishRawTranscripts
    params:
        opts = config["spliced_bam2gff_opts"]
    threads: config["threads"]
    shell:
        "spliced_bam2gff {params.opts} -t {threads} -M {input.bam} > {output.raw_gff}"


rule PinfishClusterGFF: ## cluster transcripts in GFF
    input:
        raw_gff = rules.PinfishRawBAM2GFF.output.raw_gff
    output:
        cls_gff = PinfishClusteredTranscripts,
        cls_tab = PinfishClusterMemberships,
    params:
        c = config["minimum_cluster_size"],
        d = config["exon_boundary_tolerance"],
        e = config["terminal_exon_boundary_tolerance"],
        min_iso_frac = config["minimum_isoform_percent"],
    threads: config["threads"]
    shell:
        "cluster_gff -p {params.min_iso_frac} -t {threads} -c {params.c} -d {params.d} -e {params.e} -a {output.cls_tab} {input.raw_gff} > {output.cls_gff}"
    
    
    
rule PinfishCollapseRawPartials: ## collapse clustered read artifacts
    input:
        cls_gff = rules.PinfishClusterGFF.output.cls_gff
    output:
        cls_gff_col = PinfishClusteredCollapsedTranscripts
    params:
        d = config["collapse_internal_tol"],
        e = config["collapse_three_tol"],
        f = config["collapse_five_tol"],
    shell:
       "collapse_partials -d {params.d} -e {params.e} -f {params.f} {input.cls_gff} > {output.cls_gff_col}"
    
    
    
    
rule PinfishPolishClusters: ## polish read clusters
    input:
        cls_gff = rules.PinfishClusterGFF.output.cls_gff,
        cls_tab = rules.PinfishClusterGFF.output.cls_tab,
        bam = rules.Minimap2.output.bam
    output:
        pol_trs = PinfishPolishedTranscriptFasta
    params:
        c = config["minimum_cluster_size"]
    threads: config["threads"]
    shell:
        "polish_clusters -t {threads} -a {input.cls_tab} -c {params.c} -o {output.pol_trs} {input.bam}"
    
    
rule MinimapPolishedClusters: ## map polished transcripts to genome
    input:
       index = rules.Minimap2Index.output.index,
       fasta = rules.PinfishPolishClusters.output.pol_trs,
    output:
       pol_bam = MinimapPolishedBam
    params:
        extra = config["minimap2_opts_polished"]
    threads: config["threads"]
    shell:"""
    minimap2 -t {threads} {params.extra} -ax splice {input.index} {input.fasta}\
    | samtools view -Sb -F 2304 | samtools sort -@ {threads} - -o {output.pol_bam};
    samtools index {output.pol_bam}
    """
    
    
rule PinfishPolishedBAM2GFF: ## convert BAM of polished transcripts to GFF
    input:
        bam = rules.MinimapPolishedClusters.output.pol_bam
    output:
        pol_gff = PinfishPolishedTranscripts
    params:
        extra = config["spliced_bam2gff_opts_pol"]
    threads: config["threads"]
    shell:
        "spliced_bam2gff {params.extra} -t {threads} -M {input.bam} > {output.pol_gff}"


rule PinfishCollapsePolishedPartials: ## collapse polished read artifacts
    input:
        pol_gff = rules.PinfishPolishedBAM2GFF.output.pol_gff
    output:
        pol_gff_col = PinfishPolishedCollapsedTranscripts
    params:
        d = config["collapse_internal_tol"],
        e = config["collapse_three_tol"],
        f = config["collapse_five_tol"],
    shell:
        "collapse_partials -d {params.d} -e {params.e} -f {params.f} {input.pol_gff} > {output.pol_gff_col}"    


rule PrepareCorrectedTranscriptomeFasta: ## Generate corrected transcriptome.
    input:
        genome = UnpackedReferenceFasta,
        gff = rules.PinfishCollapsePolishedPartials.output.pol_gff_col,
    output:
        fasta = PinfishCorrectedTranscriptome
    shell:"""
    gffread -g {input.genome} -w {output.fasta} {input.gff}
    """
    
