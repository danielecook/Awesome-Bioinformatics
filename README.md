Awesome-Bioinformatics [![Awesome](https://cdn.rawgit.com/sindresorhus/awesome/d7305f38d29fed78fa85652e3a63e154dd8e8829/media/badge.svg)](https://github.com/sindresorhus/awesome)
======================

A curated list of awesome Bioinformatics software, resources, and libraries. Mostly command line based, and free or open-source. Please feel free to [contribute](CONTRIBUTING.md)!

**Table of Contents**

  - [Data Processing](#data-processing)
    - [Command Line Utilities](#command-line-utilities)
  - [Next Generation Sequencing](#next-generation-sequencing)
    - [Pipelines](#pipelines)
    - [Sequence Processing](#sequence-processing)
    - [Sequence Alignment](#sequence-alignment)
    - [Variant Calling](#variant-calling)
    - [BAM File Utilities](#bam-file-utilities)
    - [VCF File Utilities](#vcf-file-utilities)
      - [Genomic Traits](#genomic-traits)
    - [Variant Simulation](#variant-simulation)
    - [Variant Filtering / Quality Control](#variant-filtering--quality-control)
    - [Variant Prediction/Annotation](#variant-predictionannotation)
    - [Python Modules](#python-modules)
  - [Visualization](#visualization)
    - [Genome Browsers / Gene diagrams](#genome-browsers--gene-diagrams)
  - [Database Access](#database-access)
  - [Resources](#resources)
    - [Sequencing](#sequencing)
    - [RNA-seq](#rna-seq)
    - [ChIP-seq](#chip-seq)
    - [Miscellaneous](#miscellaneous)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

----

## Data Processing

### Command Line Utilities

* __[datamash](http://www.gnu.org/software/datamash/)__ - Data transformations and statistics.
* __[Bioinformatics One Liners](https://github.com/stephenturner/oneliners)__ - Git repo of useful single line commands.
* __[CSVKit](https://github.com/onyxfish/csvkit)__ - Utilities for working with CSV/Tab-delimited files.
* __[Bedtools2](https://github.com/arq5x/bedtools2)__ - A Swiss Army knife for genome arithmetic.
* __[GNU `parallel`](http://www.gnu.org/software/parallel/)__ - General parallelizer that runs jobs in parallel on a single multi-core machine. [Here](https://www.biostars.org/p/63816/) are some example scripts using GNU `parallel`.

## Next Generation Sequencing

### Pipelines

* __[bcbio-nextgen](https://github.com/chapmanb/bcbio-nextgen)__ - Batteries included genomic analysis pipeline for variant and RNA-Seq analysis, structural variant calling, annotation, and prediction.

### Sequence Processing

Sequence Processing includes tasks such as demultiplexing raw read data, and trimming low quality bases.

* __[Fastqp](https://github.com/mdshw5/fastqp)__ - Fastq and Sam quality control using python.
* __[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)__ - A quality control tool for high throughput sequence data.
* __[Fastx Tookit](http://hannonlab.cshl.edu/fastx_toolkit/)__ - FASTQ/A short-reads pre-processing tools: Demultiplexing, trimming, clipping, quality filtering, and masking utilities.
* __[Seqtk](https://github.com/lh3/seqtk)__ - Toolkit for processing sequences in FASTA/Q formats.

### Sequence Alignment

__De Novo Alignment__

__DNA Resequencing__

* __[BWA](https://github.com/lh3/bwa)__ Burrow-Wheeler Aligner for pairwise alignment between DNA sequences.

### Variant Calling

* __[samtools/bcftools/htslib](https://github.com/samtools/samtools)__ - A suite of tools for manipulating next-generation sequencing data.
* __[freebayes](https://github.com/ekg/freebayes)__ - Bayesian haplotype-based polymorphism discovery and genotyping.

### BAM File Utilities

* [Bamtools](https://github.com/pezmaster31/bamtools) - Collection of tools for working with BAM files.

### VCF File Utilities

* __[vcflib](https://github.com/ekg/vcflib)__ - A C++ library for parsing and manipulating VCF files.
* __[bcftools](https://github.com/samtools/bcftools)__ - Set of tools for manipulating VCF files.
* __[vcftools](http://vcftools.sourceforge.net/)__ - VCF manipulation and statistics (e.g. linkage disequilibrium, allele frequency, Fst).

#### Genomic Traits

__Genomic Traits__ are differences in terms of DNA structure or content observed among populations that may be regulated by genetic variation. For example, telomere length or rDNA copy number.

* __[Telseq](https://github.com/zd1/telseq)__ - Telseq is a tool for estimating telomere length from whole genome sequence data.
* __[bam toolbox](https://github.com/AndersenLab/bam-toolbox)__ MtDNA:Nuclear Coverage; Bam Toolbox can output the ratio of MtDNA:nuclear coverage, a proxy for mitochondrial content.

### Variant Simulation

* __[wgsim](https://github.com/lh3/wgsim)__ - __Comes with samtools!__ - Reads simulator.
* __[Bam Surgeon](https://github.com/adamewing/bamsurgeon)__ - Tools for adding mutations to existing .bam files, used for testing mutation callers.

### Variant Filtering / Quality Control

### Variant Prediction/Annotation

* __[SIFT](http://sift.jcvi.org/)__ - Predicts whether an amino acid substitution affects protein function.
* __[SNPeff](http://snpeff.sourceforge.net/)__ - Genetic variant annotation and effect prediction toolbox.

### Python Modules

* __[cruzdb](https://github.com/brentp/cruzdb)__ - Pythonic access to the UCSC Genome database.
* __[pyfaidx](https://github.com/mdshw5/pyfaidx)__ - Pythonic access to FASTA files.
* __[pyBedTools](https://github.com/daler/pybedtools)__ - Python wrapper for [bedtools](https://github.com/arq5x/bedtools).
* __[pysam](https://github.com/pysam-developers/pysam)__ - Python wrapper for [samtools](https://github.com/samtools/samtools).
* __[pyVCF](https://github.com/jamescasbon/PyVCF)__ - A VCF Parser for python.
* __[cyvcf](https://github.com/arq5x/cyvcf)__ - A port of [pyVCF](https://github.com/jamescasbon/PyVCF) using Cython for speed.
* __[cyvcf2](https://github.com/brentp/cyvcf2)__ - cython + htslib == fast VCF parsing; Even faster parsing than pyVCF.

## Visualization

### Genome Browsers / Gene diagrams

The following tools can be used to visualize genomic data or for constructing customized visualizations of genomic data including sequence data from DNA-Seq, RNA-Seq, and ChIP-Seq, variants, and more.

* __[biodalliance](http://www.biodalliance.org/)__ - Embeddable genome viewer. Integration data from a wide variety of sources, and can load data directly from popular genomics file formats including bigWig, BAM, and VCF.
* __[IGV](https://www.broadinstitute.org/igv/)__ - Java based browser. Fast, efficient, scalable visualization tool for genomics data and annotations. Handles a large [variety of formats](http://www.broadinstitute.org/igv/FileFormats).
* __[Island Plot](https://github.com/lairdm/islandplot)__ - D3 JavaScript based genome viewer. Constructs SVGs.
* __[pileup.js](https://github.com/hammerlab/pileup.js)__ - JavaScript library that can be used to generate interactive and highly customizable web-based genome browsers.
* __[scribl](https://github.com/chmille4/Scribl)__ - JavaScript library for drawing canvas-based gene diagrams. The [Homepage](http://chmille4.github.io/Scribl/) has examples.

## Database Access

* [Entrez Direct: E-utilities on the UNIX command line](http://www.ncbi.nlm.nih.gov/books/NBK179288/) - UNIX command line tools to access NCBI's databases programmatically. Instructions to install and examples are found in the link.

## Resources

### Sequencing

* [Next-Generation Sequencing Technologies - Elaine Mardis (2014)](https://youtu.be/6Is3W7JkFp8) [1:34:35] - Excellent (technical) overview of next-generation and third-generation sequencing technologies, along with some applications in cancer research.
* [Annotated bibliography of *Seq assays](https://liorpachter.wordpress.com/seq/) - List of ~100 papers on various sequencing technologies and assays ranging from transcription to transposable element discovery.
* [For all you seq... (PDF)](http://www.illumina.com/content/dam/illumina-marketing/documents/applications/ngs-library-prep/ForAllYouSeqMethods.pdf) (3456x5471) - Massive infographic by Illumina on illustrating how many sequencing techniques work. Techniques cover protein-protein interactions, RNA transcription, RNA-protein interactions, RNA low-level detection, RNA modifications, RNA structure, DNA rearrangements and markers, DNA low-level detection, epigenetics, and DNA-protein interactions. References included.

### RNA-Seq

* [Review papers on RNA-seq (Biostars)](https://www.biostars.org/p/52152/) - Includes lots of seminal papers on RNA-seq and analysis methods.
* [Informatics for RNA-seq: A web resource for analysis on the cloud](https://github.com/griffithlab/rnaseq_tutorial) - Educational resource on performing RNA-seq analysis in the cloud using Amazon AWS cloud services. Topics include preparing the data, preprocessing, differential expression, isoform discovery, data visualization, and interpretation.
* [RNA-seqlopedia](http://rnaseq.uoregon.edu/) - RNA-seqlopedia provides an awesome overview of RNA-seq and of the choices necessary to carry out a successful RNA-seq experiment.

### ChIP-Seq

* [ChIP-seq analysis notes from Tommy Tang](https://github.com/crazyhottommy/ChIP-seq-analysis) - Resources on ChIP-seq data which include papers, methods, links to software, and analysis.

### Miscellaneous

* [The Leek group guide to genomics papers](https://github.com/jtleek/genomicspapers/) - Expertly curated genomics papers to get up to speed on genomics, RNA-seq, statistics (used in genomics), software development, and more.
