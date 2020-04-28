Awesome Bioinformatics [![Awesome](https://cdn.rawgit.com/sindresorhus/awesome/d7305f38d29fed78fa85652e3a63e154dd8e8829/media/badge.svg)](https://github.com/sindresorhus/awesome)
![URL Check](https://github.com/danielecook/Awesome-Bioinformatics/workflows/URL%20Check/badge.svg) ![TOC](https://github.com/danielecook/Awesome-Bioinformatics/workflows/TOC/badge.svg)
======================

> Bioinformatics is an interdisciplinary field that develops methods and software tools for understanding biological data. — [Wikipedia](https://en.wikipedia.org/wiki/Bioinformatics)

A curated list of awesome Bioinformatics software, resources, and libraries. Mostly command line based, and free or open-source. Please feel free to [contribute](CONTRIBUTING.md)!

<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
**Table of Contents**

- [Package suites](#package-suites)
  - [Bioconductor](#bioconductor)
  - [Biopython](#biopython)
  - [Bioconda](#bioconda)
- [Data Tools](#data-tools)
- [Data Processing](#data-processing)
  - [Command Line Utilities](#command-line-utilities)
- [Next Generation Sequencing](#next-generation-sequencing)
  - [Workflow Managers](#workflow-managers)
  - [Pipelines](#pipelines)
  - [Sequence Processing](#sequence-processing)
  - [Data Analysis](#data-analysis)
  - [Sequence Alignment](#sequence-alignment)
  - [Variant Calling](#variant-calling)
  - [BAM File Utilities](#bam-file-utilities)
  - [VCF File Utilities](#vcf-file-utilities)
  - [GFF BED File Utilities](#gff-bed-file-utilities)
  - [Variant Simulation](#variant-simulation)
  - [Variant Prediction/Annotation](#variant-predictionannotation)
  - [Python Modules](#python-modules)
    - [Data](#data)
    - [Tools](#tools)
- [Visualization](#visualization)
  - [Genome Browsers / Gene Diagrams](#genome-browsers--gene-diagrams)
  - [Circos Related](#circos-related)
- [Database Access](#database-access)
- [Resources](#resources)
  - [Becoming a Bioinformatician](#becoming-a-bioinformatician)
  - [Bioinformatics on GitHub](#bioinformatics-on-github)
  - [Sequencing](#sequencing)
  - [RNA-Seq](#rna-seq)
  - [ChIP-Seq](#chip-seq)
  - [YouTube Channels and Playlists](#youtube-channels-and-playlists)
  - [Blogs](#blogs)
  - [Miscellaneous](#miscellaneous)
- [License](#license)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

----

## Package suites

Package suites gather software packages and installation tools for specific languages or platforms. We have some for bioinformatics software.

### Bioconductor

* __[Bioconductor](https://www.bioconductor.org/)__ - A plethora of tools for analysis and comprehension of high-throughput genomic data, including 1500+ software packages.

### Biopython

* __[Biopython](https://biopython.org/)__ - Freely available tools for biological computing in Python, with included cookbook, packaging and thorough documentation. Part of the [Open Bioinformatics Foundation](http://open-bio.org/). Contains the very useful [Entrez](https://biopython.org/DIST/docs/api/Bio.Entrez-module.html) package for API access to the NCBI databases.

### Bioconda

* __[Bioconda](https://bioconda.github.io/)__ - A channel for the [conda package manager](http://conda.pydata.org/docs/intro.html) specializing in bioinformatics software. Includes a repository with 3000+ ready-to-install (with `conda install`) bioinformatics packages.

## Data Tools

* __[GGD](https://github.com/gogetdata/ggd-cli)__ - Go Get Data; A command line interface for obtaining genomic data
* __[SRA-Explorer](https://sra-explorer.info/)__ - Easily get SRA download links and other information.

## Data Processing

### Command Line Utilities

* __[Bioinformatics One Liners](https://github.com/stephenturner/oneliners)__ - Git repo of useful single line commands.
* __[BioNode](https://www.bionode.io/)__ - Modular and universal bioinformatics, Bionode provides pipeable UNIX command line tools and JavaScript APIs for bioinformatics analysis workflows.
* __[bioSyntax](http://www.bioSyntax.org/)__ - Syntax Highlighting for Computational Biology file formats (SAM, VCF, GTF, FASTA, PDB, etc...) in vim/less/gedit/sublime.
* __[CSVKit](https://github.com/wireservice/csvkit)__ - Utilities for working with CSV/Tab-delimited files.
* __[csvtk](https://github.com/shenwei356/csvtk)__ - Another cross-platform, efficient, practical and pretty CSV/TSV toolkit.
* __[datamash](http://www.gnu.org/software/datamash/)__ - Data transformations and statistics.
* __[easy_qsub](https://github.com/shenwei356/easy_qsub)__ - Easily submitting PBS jobs with script template. Multiple input files supported.
* __[GNU `parallel`](http://www.gnu.org/software/parallel/)__ - General parallelizer that runs jobs in parallel on a single multi-core machine. [Here](https://www.biostars.org/p/63816/) are some example scripts using GNU `parallel`.
* __[grabix](https://github.com/arq5x/grabix)__ - A wee tool for random access into BGZF files.
* __[gsort](https://github.com/brentp/gsort)__ - Sort genomic files according to a specified order.
* __[tabix](https://github.com/samtools/tabix)__ - Table file index.
* __[wormtable](https://github.com/wormtable/wormtable)__ - Write-once-read-many table for large datasets.
* __[zindex](https://github.com/mattgodbolt/zindex)__ - Create an index on a compressed text file.

## Next Generation Sequencing

### Workflow Managers

* __[BigDataScript](https://pcingola.github.io/BigDataScript/)__ - A cross-system scripting language for working with big data pipelines in computer systems of different sizes and capabilities.
* __[Bpipe](http://docs.bpipe.org)__ - A small language for defining pipeline stages and linking them together to make pipelines.
* __[Common Workflow Language](http://www.commonwl.org/)__ - a specification for describing analysis workflows and tools that are portable and scalable across a variety of software and hardware environments, from workstations to cluster, cloud, and high performance computing (HPC) environments.
* __[Cromwell](https://github.com/broadinstitute/cromwell)__ - A Workflow Management System geared towards scientific workflows.
* __[Galaxy](https://usegalaxy.org/)__ - a popular open-source, web-based platform for data intensive biomedical research. Has several features, from data analysis to workflow management to visualization tools.
* __[GATK Queue](https://gatkforums.broadinstitute.org/gatk/discussion/1288/howto-run-queue-for-the-first-time)__ - A pipelining system built to work natively with GATK as well as other high-throughput sequence analysis software.
* __[Nextflow](https://www.nextflow.io) (recommended)__ - A fluent DSL modelled around the UNIX pipe concept, that simplifies writing parallel and scalable pipelines in a portable manner.
* __[Ruffus](http://www.ruffus.org.uk)__ - Computation Pipeline library for python widely used in science and bioinformatics.
* __[SeqWare](https://seqware.github.io/)__ - Hadoop Oozie-based workflow system focused on genomics data analysis in cloud environments.
* __[Snakemake](https://bitbucket.org/snakemake/snakemake/wiki/Home)__ - A workflow management system in Python that aims to reduce the complexity of creating workflows by providing a fast and comfortable execution environment.
* __[Workflow Descriptor Language](https://github.com/broadinstitute/wdl)__ - Workflow standard developed by the Broad.

### Pipelines

* __[Awesome-Pipeline](https://github.com/pditommaso/awesome-pipeline)__ - A list of pipeline resources.
* __[bcbio-nextgen](https://github.com/chapmanb/bcbio-nextgen)__ - Batteries included genomic analysis pipeline for variant and RNA-Seq analysis, structural variant calling, annotation, and prediction.
* __[R-Peridot](http://www.bioinformatics-brazil.org/r-peridot/)__ - Customizable pipeline for differential expression analysis with an intuitive GUI.

### Sequence Processing

Sequence Processing includes tasks such as demultiplexing raw read data, and trimming low quality bases.

* __[AfterQC](https://github.com/OpenGene/AfterQC)__ - Automatic Filtering, Trimming, Error Removing and Quality Control for fastq data
* __[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)__ - A quality control tool for high throughput sequence data.
* __[Fastqp](https://github.com/mdshw5/fastqp)__ - FASTQ and SAM quality control using Python.
* __[Fastx Tookit](http://hannonlab.cshl.edu/fastx_toolkit/)__ - FASTQ/A short-reads pre-processing tools: Demultiplexing, trimming, clipping, quality filtering, and masking utilities.
* __[MultiQC](http://multiqc.info/)__ - Aggregate results from bioinformatics analyses across many samples into a single report.
* __[SeqKit](https://github.com/shenwei356/seqkit)__ - A cross-platform and ultrafast toolkit for FASTA/Q file manipulation in Golang.
* __[seqmagick](http://seqmagick.readthedocs.io/en/latest/)__ - file format conversion in Biopython in a convenient way
* __[Seqtk](https://github.com/lh3/seqtk)__ - Toolkit for processing sequences in FASTA/Q formats.
* __[smof](https://github.com/incertae-sedis/smof)__ - UNIX-style FASTA manipulation tools.

### Data Analysis

The following items allow for scalable genomic analysis by introducing specialized databases.

* __[Hail](https://github.com/hail-is/hail)__ - Scalable genomic analysis
* __[GLNexus](https://github.com/dnanexus-rnd/GLnexus)__ - Scalable gVCF merging and joint variant calling for population sequencing projects.

### Sequence Alignment

__De Novo Alignment__

__DNA Resequencing__

* __[Bowtie 2](https://github.com/BenLangmead/bowtie2)__ - An ultrafast and memory-efficient tool for aligning sequencing reads to long reference sequences.
* __[BWA](https://github.com/lh3/bwa)__ - Burrow-Wheeler Aligner for pairwise alignment between DNA sequences.

### Variant Calling

* __[freebayes](https://github.com/ekg/freebayes)__ - Bayesian haplotype-based polymorphism discovery and genotyping.
* __[GATK](https://software.broadinstitute.org/gatk/)__ - Variant Discovery in High-Throughput Sequencing Data
* __[samtools/bcftools/htslib](https://github.com/samtools/samtools)__ - A suite of tools for manipulating next-generation sequencing data.

__Structural variant callers__

* __[Delly](https://github.com/dellytools/delly)__ - Structural variant discovery by integrated paired-end and split-read analysis.
* __[lumpy](https://github.com/arq5x/lumpy-sv)__ - lumpy: a general probabilistic framework for structural variant discovery.
* __[manta](https://github.com/Illumina/manta)__ - Structural variant and indel caller for mapped sequencing data.
* __[gridss](https://github.com/PapenfussLab/gridss)__ - GRIDSS: the Genomic Rearrangement IDentification Software Suite.
* __[smoove](https://github.com/brentp/smoove)__ - structural variant calling and genotyping with existing tools, but,smoothly.

### BAM File Utilities

* __[Bamtools](https://github.com/pezmaster31/bamtools)__ - Collection of tools for working with BAM files.
* __[bam toolbox](https://github.com/AndersenLab/bam-toolbox)__ MtDNA:Nuclear Coverage; BAM Toolbox can output the ratio of MtDNA:nuclear coverage, a proxy for mitochondrial content.
* __[mergesam](https://github.com/DarwinAwardWinner/mergesam)__ - Automate common SAM & BAM conversions.
* __[mosdepth](https://github.com/brentp/mosdepth)__ - fast BAM/CRAM depth calculation for WGS, exome, or targeted sequencing
* __[SAMstat](https://github.com/TimoLassmann/samstat)__ - Displaying sequence statistics for next-generation sequencing.
* __[Somalier](https://github.com/brentp/mosdepth)__ - Fast sample-swap and relatedness checks on BAMs/CRAMs/VCFs/GVCFs.
* __[Telseq](https://github.com/zd1/telseq)__ - Telseq is a tool for estimating telomere length from whole genome sequence data.

### VCF File Utilities

* __[bcftools](https://github.com/samtools/bcftools)__ - Set of tools for manipulating VCF files.
* __[vcfanno](https://github.com/brentp/vcfanno)__ - Annotate a VCF with other VCFs/BEDs/tabixed files.
* __[vcflib](https://github.com/vcflib/vcflib)__ - A C++ library for parsing and manipulating VCF files.
* __[vcftools](https://github.com/vcftools/vcftools)__ - VCF manipulation and statistics (e.g. linkage disequilibrium, allele frequency, Fst).

### GFF BED File Utilities

* __[gffutils](https://github.com/daler/gffutils)__ - GFF and GTF file manipulation and interconversion.
* __[BEDOPS](https://bedops.readthedocs.io/en/latest/index.html)__ - The fast, highly scalable and easily-parallelizable genome analysis toolkit.
* __[Bedtools2](https://github.com/arq5x/bedtools2)__ - A Swiss Army knife for genome arithmetic.

### Variant Simulation

* __[Bam Surgeon](https://github.com/adamewing/bamsurgeon)__ - Tools for adding mutations to existing `.bam` files, used for testing mutation callers.
* __[wgsim](https://github.com/lh3/wgsim)__ - __Comes with samtools!__ - Reads simulator.

### Variant Prediction/Annotation

* __[SIFT](http://sift.jcvi.org/)__ - Predicts whether an amino acid substitution affects protein function.
* __[SnpEff](https://github.com/pcingola/SnpEff)__ - Genetic variant annotation and effect prediction toolbox.

### Python Modules

#### Data

* __[cruzdb](https://github.com/brentp/cruzdb)__ - Pythonic access to the UCSC Genome database.
* __[pyensembl](https://github.com/openvax/pyensembl)__ - Pythonic Access to the Ensembl database.
* __[bioservices](https://github.com/cokelaer/bioservices)__ - Access to Biological Web Services from Python.

#### Tools

* __[cyvcf](https://github.com/arq5x/cyvcf)__ - A port of [pyVCF](https://github.com/jamescasbon/PyVCF) using Cython for speed.
* __[cyvcf2](https://github.com/brentp/cyvcf2)__ - Cython + HTSlib == fast VCF parsing; even faster parsing than pyVCF.
* __[pyBedTools](https://github.com/daler/pybedtools)__ - Python wrapper for [bedtools](https://github.com/arq5x/bedtools).
* __[pyfaidx](https://github.com/mdshw5/pyfaidx)__ - Pythonic access to FASTA files.
* __[pysam](https://github.com/pysam-developers/pysam)__ - Python wrapper for [samtools](https://github.com/samtools/samtools).
* __[pyVCF](https://github.com/jamescasbon/PyVCF)__ - A VCF Parser for Python.

## Visualization

### Genome Browsers / Gene Diagrams

The following tools can be used to visualize genomic data or for constructing customized visualizations of genomic data including sequence data from DNA-Seq, RNA-Seq, and ChIP-Seq, variants, and more.

* __[Squiggle](https://github.com/Lab41/squiggle)__ - Easy-to-use DNA sequence visualization tool that turns FASTA files into browser-based visualizations.
* __[biodalliance](http://www.biodalliance.org/)__ - Embeddable genome viewer. Integration data from a wide variety of sources, and can load data directly from popular genomics file formats including bigWig, BAM, and VCF.
* __[BioJS](http://biojs.net/)__ - BioJS is a library of over hundred JavaScript components enabling you to visualize and process data using current web technologies.
* __[Circleator](https://github.com/jonathancrabtree/Circleator)__ - Flexible circular visualization of genome-associated data with BioPerl and SVG.
* __[DNAism](https://github.com/drio/dnaism)__ - Horizon chart D3-based JavaScript library for DNA data.
* __[IGV js](https://www.broadinstitute.org/igv)__ - Java-based browser. Fast, efficient, scalable visualization tool for genomics data and annotations. Handles a large [variety of formats](http://software.broadinstitute.org/software/igv/fileformats).
* __[Island Plot](https://github.com/lairdm/islandplot)__ - D3 JavaScript based genome viewer. Constructs SVGs.
* __[JBrowse](https://jbrowse.org)__ - JavaScript genome browser that is highly customizable via plugins and track customizations
* __[PHAT](https://github.com/chgibb/PHAT)__ - Point and click, cross platform suite for analysing and visualizing next-generation sequencing datasets.
* __[pileup.js](https://github.com/hammerlab/pileup.js)__ - JavaScript library that can be used to generate interactive and highly customizable web-based genome browsers.
* __[scribl](https://github.com/chmille4/Scribl)__ - JavaScript library for drawing canvas-based gene diagrams. The [Homepage](http://chmille4.github.io/Scribl/) has examples.
* __[Lucid Align](https://lucidalign.com/)__ - A modern sequence alignment viewer

### Circos Related

* __[Circos](http://circos.ca/)__ - Perl package for circular plots, which are well suited for genomic rearrangements.
* __[ClicO FS](https://academic.oup.com/bioinformatics/article/31/22/3685/241292)__ - An interactive web-based service of Circos.
* __[OmicCircos](http://www.bioconductor.org/packages/release/bioc/html/OmicCircos.html)__ -  R package for circular plots for omics data.
* __[J-Circos](http://www.australianprostatecentre.org/research/software/jcircos)__ - A Java application for doing interactive work with circos plots.
* __[rCircos](https://cran.r-project.org/web/packages/RCircos/index.html)__ - R package for circular plots.

## Database Access

* [Entrez Direct: E-utilities on the UNIX command line](http://www.ncbi.nlm.nih.gov/books/NBK179288/) - UNIX command line tools to access NCBI's databases programmatically. Instructions to install and examples are found in the link.

## Resources

### Becoming a Bioinformatician

* [What is a bioinformatician](http://blog.fejes.ca/?p=2418)
* [Bioinformatics Curriculum Guidelines: Toward a Definition of Core Competencies](http://www.ploscompbiol.org/article/info:doi%2F10.1371%2Fjournal.pcbi.1003496)
* [Top N Reasons To Do A Ph.D. or Post-Doc in Bioinformatics/Computational Biology](http://caseybergman.wordpress.com/2012/07/31/top-n-reasons-to-do-a-ph-d-or-post-doc-in-bioinformaticscomputational-biology/)
* [A 10-Step Guide to Party Conversation For Bioinformaticians](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2013-14-1-104) - Here is a step-by-step guide on how to convey concepts to people not involved in the field when asked the question: 'So, what do you do?'
* [A History Of Bioinformatics (In The Year 2039)](https://www.youtube.com/watch?v=uwsjwMO-TEA) - A talk by C. Titus Brown on his take of looking back at bioinformatics from the year 2039. His notes for this talk can be found [here](http://ivory.idyll.org/blog/2014-bosc-keynote.html).
* [A farewell to bioinformatics](http://madhadron.com/posts/2012-03-26-a-farewell-to-bioinformatics.html) - A critical view of the state of bioinformatics.
* [A Series of Interviews with Notable Bioinformaticians](http://www.acgt.me/blog/2014/3/25/101-questions-a-new-series-of-interviews-with-notable-bioinformaticians) - Dr. Keith Bradnam "thought it might be instructive to ask a simple series of questions to a bunch of notable bioinformaticians to assess their feelings on the current state of bioinformatics research, and maybe get any tips they have about what has been useful to their bioinformatics careers."
* [Open Source Society University on Bioinformatics](https://github.com/ossu/bioinformatics) - Solid path for those of you who want to complete a Bioinformatics course on your own time, for free, with courses from the best universities in the World.
* [Rosalind](http://rosalind.info/) - Rosalind is a platform for learning bioinformatics through problem solving.
* [A guide for the lonely bioinformatician](http://www.opiniomics.org/a-guide-for-the-lonely-bioinformatician/) - This guide is aimed at bioinformaticians, and is meant to guide them towards better career development.
* [A brief history of bioinformatics](https://doi.org/10.1093/bib/bby063)

### Bioinformatics on GitHub

* [Awesome-alternative-splicing](https://github.com/HussainAther/awesome-alternative-splicing) - List of resources on alternative splicing including software, databases, and other tools..

### Sequencing

* [Next-Generation Sequencing Technologies - Elaine Mardis (2014)](https://youtu.be/6Is3W7JkFp8) [1:34:35] - Excellent (technical) overview of next-generation and third-generation sequencing technologies, along with some applications in cancer research.
* [Annotated bibliography of \*Seq assays](https://liorpachter.wordpress.com/seq/) - List of ~100 papers on various sequencing technologies and assays ranging from transcription to transposable element discovery.
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
* [A New Online Computational Biology Curriculum](https://doi.org/10.1371/journal.pcbi.1003662) - "This article introduces a catalog of several hundred free video courses of potential interest to those wishing to expand their knowledge of bioinformatics and computational biology. The courses are organized into eleven subject areas modeled on university departments and are accompanied by commentary and career advice."
* [How Perl Saved the Human Genome Project](http://www.foo.be/docs/tpj/issues/vol1_2/tpj0102-0001.html) - An anecdote by Lincoln D. Stein on the importance of the Perl programming language in the Human Genome Project.
* [Educational Papers from Nature Biotechnology and PLoS Computational Biology](https://liacs.leidenuniv.nl/~hoogeboomhj/mcb/nature_primer.html) - Page of links to primers and short educational articles on various methods used in computational biology and bioinformatics.
* [The PeerJ Bioinformatics Software Tools Collection](https://peerj.com/collections/45-bioinformatics-software/) - Collection of tools curated by Keith Crandall and Claus White, aimed at collating the most interesting, innovative, and relevant bioinformatics tools articles in PeerJ.

## License

[![CC0](http://mirrors.creativecommons.org/presskit/buttons/88x31/svg/cc-zero.svg)](https://creativecommons.org/publicdomain/zero/1.0/)
