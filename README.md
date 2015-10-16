Awesome-Bioinformatics
======================

A curated list of awesome Bioinformatics software and libraries. Mostly command line based, and free or open-source.



## Data Processing

### Command Line Utilities

| Program                                                                 | Description
|:----------------------------------------------------------------------- | :------------
| [datamash](http://www.gnu.org/software/datamash/)                       | Data transformations and statistics. 
| [Bioinformatics One Liners](https://github.com/stephenturner/oneliners) | Git repo of useful single line commands.
| [CSVKit](https://github.com/onyxfish/csvkit) | Utilities for working with CSV/Tab-delimited files.
| [Bedtools2](https://github.com/arq5x/bedtools2)                         | a swiss army knife for genome arithmetic

## Next Generation Sequencing

### Pipelines

| Program                                                                 | Description
|:----------------------------------------------------------------------- | :------------
| [bcbio-nextgen](https://github.com/chapmanb/bcbio-nextgen)              | Batteries included genomic analysis pipeline for variant and RNA-Seq analysis, structural variant calling, annotation, and prediction.

### Sequence Processing 

Sequence Processing includes tasks such as demultiplexing raw read data, and trimming low quality bases.

| Program                                                                 | Description
|:----------------------------------------------------------------------- | :------------
| [Fastqp](https://github.com/mdshw5/fastqp) - Fastq and Sam quality control using python.            | Fastq and Sam quality control using python.
| [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)     |  A quality control tool for high throughput sequence data.
| [Fastx Tookit](http://hannonlab.cshl.edu/fastx_toolkit/) | FASTQ/A short-reads pre-processing tools: Demultiplexing, trimming, clipping, quality filtering, and masking utilities.
| [Seqtk](https://github.com/lh3/seqtk) | Toolkit for processing sequences in FASTA/Q formats. |

### Sequence Alignment

__De Novo Alignment__

__DNA Resequencing__

| Program                                                                 | Description
|:----------------------------------------------------------------------- | :------------
| [BWA](https://github.com/lh3/bwa) | Burrow-Wheeler Aligner for pairwise alignment between DNA sequences. 

### Variant Calling

| Program                                                                 | Description
|:----------------------------------------------------------------------- | :------------
| [samtools/bcftools/htslib](https://github.com/samtools/samtools) | A suite of tools for manipulating next-generation sequencing data.
| [freebayes](https://github.com/ekg/freebayes) | Bayesian haplotype-based polymorphism discovery and genotyping.

### VCF File Utilities

| Program                                                                 | Description
|:----------------------------------------------------------------------- | :------------
| [vcflib](https://github.com/ekg/vcflib) | A C++ library for parsing and manipulating VCF files.
| [bcftools](https://github.com/samtools/bcftools) | Set of tools for manipulating vcf files.
| [vcftools](http://vcftools.sourceforge.net/) | VCF manipulation and statistics (e.g. linkage disequilibrium, allele frequency, Fst)


#### Genomic Traits

__Genomic Traits__ are differences in terms of DNA structure or content observed among populations that may be regulated by genetic variation. For example, telomere length or rDNA copy number.

| Program | Description
|:----------|:---------------------
| [Telseq](https://github.com/zd1/telseq) | Telseq is a tool for estimating telomere length from whole genome sequence data.

### Variant Simulation

| Program                                                                 | Description
|:----------------------------------------------------------------------- | :------------
| [wgsim](https://github.com/lh3/wgsim) | __Comes with samtools!__ - Reads simulator
| [Bam Surgeon](https://github.com/adamewing/bamsurgeon) | Tools for adding mutations to existing .bam files, used for testing mutation callers.

### Variant Filtering / Quality Control

### Variant Prediction/Annotation

| Program                                                                 | Description
|:----------------------------------------------------------------------- | :------------
| [SNPeff](http://snpeff.sourceforge.net/) | Genetic variant annotation and effect prediction toolbox. 


### Python Modules

| Program                                                                 | Description
|:----------------------------------------------------------------------- | :------------
| [pyBedTools](https://github.com/daler/pybedtools)                       | Python wrapper for [bedtools](https://github.com/arq5x/bedtools). 
| [pysam](https://github.com/pysam-developers/pysam)                      | Python wrapper for [samtools](https://github.com/samtools/samtools).
| [pyVCF](https://github.com/jamescasbon/PyVCF)                           | A VCF Parser for python.
| [cyVCF](https://github.com/arq5x/cyvcf)                                 | A fast Python library for VCF files leveraging Cython for speed. Based off of pyVCF.
| [Biopython](https://github.com/biopython/biopython) | Large suite of bioinformatic utilities.

## Visualization

### Gene Diagrams

_Tools for drawing gene diagrams / alignments_

* [scribl](https://github.com/chmille4/Scribl) - javascript library for drawing canvas-based gene diagrams. The [Homepage](http://chmille4.github.io/Scribl/) has examples.
* [pileup.js](https://github.com/hammerlab/pileup.js) - javascript library that can be used to generate interactive and highly customizable web-based genome browsers.

### Genome Browsers

* [IGV](https://www.broadinstitute.org/igv/) - Java based browser.  Fast, efficient, scalable visualization tool for genomics data and annotations. Handles a [variety of formats](http://www.broadinstitute.org/igv/FileFormats).
* [pileup.js]

> Integrated Genome Browser - Java based browser |


