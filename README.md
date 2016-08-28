Awesome Bioinformatics [![Awesome](https://cdn.rawgit.com/sindresorhus/awesome/d7305f38d29fed78fa85652e3a63e154dd8e8829/media/badge.svg)](https://github.com/sindresorhus/awesome)
======================

> Bioinformatics is an interdisciplinary field that develops methods and software tools for understanding biological data. — [Wikipedia](https://en.wikipedia.org/wiki/Bioinformatics)

A curated list of awesome Bioinformatics software, resources, and libraries. Mostly command line based, and free or open-source. Please feel free to [contribute](CONTRIBUTING.md)!

**Table of Contents**

<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->


- [Data Processing](#data-processing)
  - [Command Line Utilities](#command-line-utilities)
- [Next Generation Sequencing](#next-generation-sequencing)
  - [Pipelines/Pipeline frameworks](#pipelinespipeline-frameworks)
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
    - [Data](#data)
    - [Tools](#tools)
- [Visualization](#visualization)
  - [Genome Browsers / Gene diagrams](#genome-browsers--gene-diagrams)
- [Database Access](#database-access)
- [Resources](#resources)
  - [Becoming a Bioinformatician](#becoming-a-bioinformatician)
  - [Sequencing](#sequencing)
  - [RNA-Seq](#rna-seq)
  - [ChIP-Seq](#chip-seq)
  - [YouTube Channels and Playlists](#youtube-channels-and-playlists)
  - [Blogs](#blogs)
  - [Miscellaneous](#miscellaneous)
- [License](#license)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

----

## Data Processing

### Command Line Utilities

* __[datamash](http://www.gnu.org/software/datamash/)__ - Data transformations and statistics.
* __[Bioinformatics One Liners](https://github.com/stephenturner/oneliners)__ - Git repo of useful single line commands.
* __[Bedtools2](https://github.com/arq5x/bedtools2)__ - A Swiss Army knife for genome arithmetic.
* __[CSVKit](https://github.com/onyxfish/csvkit)__ - Utilities for working with CSV/Tab-delimited files.
* __[csvtk](https://github.com/shenwei356/csvtk)__ - Another cross-platform, efficient, practical and pretty CSV/TSV toolkit.
* __[easy_qsub](https://github.com/shenwei356/easy_qsub)__ - Easily submitting PBS jobs with script template. Multiple input files supported.
* __[GNU `parallel`](http://www.gnu.org/software/parallel/)__ - General parallelizer that runs jobs in parallel on a single multi-core machine. [Here](https://www.biostars.org/p/63816/) are some example scripts using GNU `parallel`.



## Next Generation Sequencing

### Pipelines/Pipeline frameworks
* __[Ruffus](http://www.ruffus.org.uk)__ - Computation Pipeline library for python widely used in science and bioinformatics.
* __[Snakemake](https://bitbucket.org/snakemake/snakemake/wiki/Home)__ - A workflow management system in Python that aims to reduce the complexity of creating workflows by providing a fast and comfortable execution environment.
* __[Nextflow](https://www.nextflow.io)__ - A fluent DSL modelled around the UNIX pipe concept, that simplifies writing parallel and scalable pipelines in a portable manner.
* __[BigDataScript](https://pcingola.github.io/BigDataScript)__ - A cross-system scripting language for working with big data pipelines in computer systems of different sizes and capabilities.
* __[Bpipe](http://docs.bpipe.org)__ - A small language for defining pipeline stages and linking them together to make pipelines.
* __[GATK Queue](http://gatkforums.broadinstitute.org/firecloud/discussion/1306/overview-of-queue)__ - A pipelining system built to work natively with GATK as well as other high-throughput sequence analysis software.
* __[SeqWare](https://seqware.github.io/)__ - Hadoop Oozie-based workflow system focused on genomics data analysis in cloud environments.
* __[bcbio-nextgen](https://github.com/chapmanb/bcbio-nextgen)__ - Batteries included genomic analysis pipeline for variant and RNA-Seq analysis, structural variant calling, annotation, and prediction.

### Sequence Processing

Sequence Processing includes tasks such as demultiplexing raw read data, and trimming low quality bases.

* __[fakit](https://github.com/shenwei356/fakit)__ - A cross-platform and efficient toolkit for FASTA/Q file manipulation.
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
* __[SnpEff](http://snpeff.sourceforge.net/)__ - Genetic variant annotation and effect prediction toolbox.

### Python Modules

#### Data

* __[cruzdb](https://github.com/brentp/cruzdb)__ - Pythonic access to the UCSC Genome database.
* __[pyensembl](https://github.com/hammerlab/pyensembl)__ - Pythonic Access to the Ensembl database.

#### Tools

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
* __[IGV js](https://www.broadinstitute.org/igv/)__ - Java based browser. Fast, efficient, scalable visualization tool for genomics data and annotations. Handles a large [variety of formats](http://www.broadinstitute.org/igv/FileFormats). 
* __[Island Plot](https://github.com/lairdm/islandplot)__ - D3 JavaScript based genome viewer. Constructs SVGs.
* __[pileup.js](https://github.com/hammerlab/pileup.js)__ - JavaScript library that can be used to generate interactive and highly customizable web-based genome browsers.
* __[scribl](https://github.com/chmille4/Scribl)__ - JavaScript library for drawing canvas-based gene diagrams. The [Homepage](http://chmille4.github.io/Scribl/) has examples.
* __[DNAism](https://github.com/drio/dnaism)__ - Horizon chart d3-based js library for DNA data.

## Database Access

* [Entrez Direct: E-utilities on the UNIX command line](http://www.ncbi.nlm.nih.gov/books/NBK179288/) - UNIX command line tools to access NCBI's databases programmatically. Instructions to install and examples are found in the link.

## Resources

### Becoming a Bioinformatician

* [What is a bioinformatician](http://blog.fejes.ca/?p=2418)
* [Bioinformatics Curriculum Guidelines: Toward a Definition of Core Competencies](http://www.ploscompbiol.org/article/info:doi%2F10.1371%2Fjournal.pcbi.1003496)
* [Top N Reasons To Do A Ph.D. or Post-Doc in Bioinformatics/Computational Biology](http://caseybergman.wordpress.com/2012/07/31/top-n-reasons-to-do-a-ph-d-or-post-doc-in-bioinformaticscomputational-biology/)
* [A 10-Step Guide to Party Conversation For Bioinformaticians](http://genomebiology.biomedcentral.com/articles/10.1186/gb-2013-14-1-104) - Here is a step-by-step guide on how to convey concepts to people not involved in the field when asked the question: 'So, what do you do?'
* [A History Of Bioinformatics (In The Year 2039)](http://video.open-bio.org/video/1/a-history-of-bioinformatics-in-the-year-2039) - A talk by C. Titus Brown on his take of looking back at bioinformatics from the year 2039. His notes for this talk can be found [here](http://ivory.idyll.org/blog/2014-bosc-keynote.html).
* [A farewell to bioinformatics](http://madhadron.com/posts/2012-03-26-a-farewell-to-bioinformatics.html) - A critical view of the state of bioinformatics.
* [A Series of Interviews with Notable Bioinformaticians](http://www.acgt.me/blog/2014/3/25/101-questions-a-new-series-of-interviews-with-notable-bioinformaticians) - Dr. Keith Bradnam "thought it might be instructive to ask a simple series of questions to a bunch of notable bioinformaticians to assess their feelings on the current state of bioinformatics research, and maybe get any tips they have about what has been useful to their bioinformatics careers."
* [Learning Resources Index](https://gitlab.com/genomic/learning-resources) - [Adrián E. Salatino's](https://gitlab.com/u/Adrianzo) attempt at consolidating useful links and resources he has found helpful in his graduate career, ranging from (but not limited to) programming help, bioinformatics software, and even blogs to follow.

### Sequencing

* [Next-Generation Sequencing Technologies - Elaine Mardis (2014)](https://youtu.be/6Is3W7JkFp8) [1:34:35] - Excellent (technical) overview of next-generation and third-generation sequencing technologies, along with some applications in cancer research.
* [Annotated bibliography of *Seq assays](https://liorpachter.wordpress.com/seq/) - List of ~100 papers on various sequencing technologies and assays ranging from transcription to transposable element discovery.
* [For all you seq... (PDF)](http://www.illumina.com/content/dam/illumina-marketing/documents/applications/ngs-library-prep/ForAllYouSeqMethods.pdf) (3456x5471) - Massive infographic by Illumina on illustrating how many sequencing techniques work. Techniques cover protein-protein interactions, RNA transcription, RNA-protein interactions, RNA low-level detection, RNA modifications, RNA structure, DNA rearrangements and markers, DNA low-level detection, epigenetics, and DNA-protein interactions. References included.

### RNA-Seq

* [Review papers on RNA-seq (Biostars)](https://www.biostars.org/p/52152/) - Includes lots of seminal papers on RNA-seq and analysis methods.
* [Informatics for RNA-seq: A web resource for analysis on the cloud](https://github.com/griffithlab/rnaseq_tutorial) - Educational resource on performing RNA-seq analysis in the cloud using Amazon AWS cloud services. Topics include preparing the data, preprocessing, differential expression, isoform discovery, data visualization, and interpretation.
* [RNA-seqlopedia](http://rnaseq.uoregon.edu/) - RNA-seqlopedia provides an awesome overview of RNA-seq and of the choices necessary to carry out a successful RNA-seq experiment.
* [A survey of best practices for RNA-seq data analysis](http://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0881-8) - Gives awesome roadmap for RNA-seq computational analyses, including challenges/obstacles and things to look out for, but also how you might integrate RNA-seq data with other data types.
* [Stories from the Supplement](https://www.youtube.com/watch?v=5NiFibnbE8o) [46:39] - Dr. Lior Pachter shares his stories from the supplement for well-known RNA-seq analysis software CuffDiff and [Cufflinks](http://cole-trapnell-lab.github.io/cufflinks/) and explains some of their methodologies.
* [List of RNA-seq Bioinformatics Tools](https://en.wikipedia.org/wiki/List_of_RNA-Seq_bioinformatics_tools) - Extensive list on Wikipedia of RNA-seq bioinformatics tools needed in analysis, ranging from all parts of an analysis pipeline from quality control, alignment, splice analysis, and visualizations.
* [RNA-seq Analysis](https://github.com/crazyhottommy/RNA-seq-analysis) - [@crazyhottommy](https://github.com/crazyhottommy)'s notes on various steps and considerations when doing RNA-seq analysis.

### ChIP-Seq

* [ChIP-seq analysis notes from Tommy Tang](https://github.com/crazyhottommy/ChIP-seq-analysis) - Resources on ChIP-seq data which include papers, methods, links to software, and analysis.

### YouTube Channels and Playlists

* [Current Topics in Genome Analysis 2016](https://www.genome.gov/12514288/current-topics-in-genome-analysis-2016-course-syllabus-handouts-and-videos/) - Excellent series of fourteen lectures given at NIH about current topics in genomics ranging from sequence analysis, to sequencing technologies, and even more translational topics such as genomic medicine.
* [GenomeTV](https://www.youtube.com/user/GenomeTV) - "GenomeTV is NHGRI's collection of official video resources from lectures, to news documentaries, to full video collections of meetings that tackle the research, issues and clinical applications of genomic research."
* [Leading Strand](https://www.youtube.com/user/LeadingStrand) - Keynote lectures from Cold Spring Harbor Laboratory (CSHL) Meetings. More on [The Leading Strand](http://theleadingstrand.cshl.edu/).
* [Genomics, Big Data and Medicine Seminar Series](https://www.youtube.com/playlist?list=PLqLDR0CTP9_pboZCk6gR9Zn4kW7h9XWJI) - "Our seminars are dedicated to the critical intersection of GBM, delving into 'bleeding edge' technology and approaches that will deeply shape the future."
* [Rafael Irizarry's Channel](https://www.youtube.com/user/RafalabChannel/videos) - Dr. Rafael Irizarry's lectures and academic talks on statistics for genomics.
* [NIH VideoCasting and Podcasting](https://www.youtube.com/user/nihvcast) - "NIH VideoCast broadcasts seminars, conferences and meetings live to a world-wide audience over the Internet as a real-time streaming video." Not exclusively genomics and bioinformatics video but many great talks on domain specific use of bioinformatics and genomics.

### Blogs

* [ACGT](http://www.acgt.me/) - Dr. Keith Bradnam writes about this "thoughts on biology, genomics, and the ongoing threat to humanity from the bogus use of bioinformatics acroynums."
* [Opiniomics](http://www.opiniomics.org/) - Dr. Mick Watson write on bioinformatics, genomes, and biology.
* [Bits of DNA](https://liorpachter.wordpress.com/) - Dr. Lior Pachter writes review and commentary on computational biology.
* [it is NOT junk](http://www.michaeleisen.org/blog/) - Dr. Michael Eisen writes "a blog about genomes, DNA, evolution, open science, baseball and other important things"

### Miscellaneous

* [The Leek group guide to genomics papers](https://github.com/jtleek/genomicspapers/) - Expertly curated genomics papers to get up to speed on genomics, RNA-seq, statistics (used in genomics), software development, and more.
* [A New Online Computational Biology Curriculum](http://dx.doi.org/10.1371/journal.pcbi.1003662) - "This article introduces a catalog of several hundred free video courses of potential interest to those wishing to expand their knowledge of bioinformatics and computational biology. The courses are organized into eleven subject areas modeled on university departments and are accompanied by commentary and career advice."
* [How Perl Saved the Human Genome Project](http://www.foo.be/docs/tpj/issues/vol1_2/tpj0102-0001.html) - An anecdote by Lincoln D. Stein on the importance of the Perl programming language in the Human Genome Project.
* [Educational Papers from Nature Biotechnology and PLoS Computational Biology](https://liacs.leidenuniv.nl/~hoogeboomhj/mcb/nature_primer.html) - Page of links to primers and short educational articles on various methods used in computational biology and bioinformatics.

## License

[![CC0](http://mirrors.creativecommons.org/presskit/buttons/88x31/svg/cc-zero.svg)](https://creativecommons.org/publicdomain/zero/1.0/)
