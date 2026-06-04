# Bioinformatics Tools

![Bioinformatics](https://img.shields.io/badge/Bioinformatics-Tools-blue)
![PRs Welcome](https://img.shields.io/badge/PRs-welcome-brightgreen)

> Collection of bioinformatic tools — a curated list of popular bioinformatics software organized by category.

[**中文版请见 → README.zh-CN.md**](./README.zh-CN.md)

## Table of Contents

- [1. Read Preprocessing / Quality Control](#1-read-preprocessing--quality-control)
- [2. Sequence Alignment](#2-sequence-alignment)
- [3. Sequence Similarity Search](#3-sequence-similarity-search)
- [4. Genome Assembly](#4-genome-assembly)
- [5. Third Generation Sequencing Analysis](#5-third-generation-sequencing-analysis)
- [6. RNA-seq Analysis](#6-rna-seq-analysis)
- [7. Variant Calling & Annotation](#7-variant-calling--annotation)
- [8. CRISPR Targeted Deep Sequencing Analysis](#8-crispr-targeted-deep-sequencing-analysis)
- [9. GWAS (Genome-Wide Association Studies)](#9-gwas-genome-wide-association-studies)
- [10. Population Genetics](#10-population-genetics)
- [11. Structural Variation Detection](#11-structural-variation-detection)
- [12. Metagenomics](#12-metagenomics)
- [13. Epigenomics](#13-epigenomics)
- [14. Proteomics](#14-proteomics)
- [15. Single-Cell Analysis](#15-single-cell-analysis)
- [16. Phylogenetics](#16-phylogenetics)
- [17. Protein Structure & Molecular Dynamics](#17-protein-structure--molecular-dynamics)
- [18. Hi-C & Chromatin Conformation](#18-hi-c--chromatin-conformation)
- [19. General-Purpose & Workflow Management](#19-general-purpose--workflow-management)

---

## 1. Read Preprocessing / Quality Control

Quality control, adapter trimming, read merging, error correction, and filtering tools for raw sequencing data.

- [**FastQC**](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) - [commandLine] - FastQC provides quality control checks on raw sequencing data. It generates modular HTML reports covering per-base quality, GC content, N content, sequence length distribution, overrepresented sequences, and adapter content, helping identify data issues before downstream analysis.
- [**MultiQC**](https://github.com/ewels/MultiQC) - [commandLine] - MultiQC aggregates results from bioinformatics analyses across many samples into a single interactive HTML report. It automatically searches directories for analysis logs and compiles them into a summary with interactive plots, supporting FastQC, STAR, HISAT2, Salmon, and many other tools.
- [**fastp**](https://github.com/OpenGene/fastp) - [commandLine] - An all-in-one FASTQ preprocessing tool that performs QC, adapter trimming, quality filtering, per-read quality pruning, poly-G tail trimming, and generates comprehensive HTML/JSON reports in a single pass, replacing multiple conventional tools.
- [**Trimmomatic**](https://github.com/usadellab/Trimmomatic) - [commandLine] - A flexible read trimming tool for Illumina NGS data. It performs adapter clipping, sliding-window quality trimming, minimum length filtering, and leading/trailing low-quality base removal, supporting both single-end and paired-end data.
- [**Cutadapt**](https://cutadapt.readthedocs.io/) - [commandLine] - Finds and removes adapter sequences, primers, poly-A tails, and other unwanted sequences from high-throughput sequencing reads. Supports both single-end and paired-end data with error-tolerant adapter matching.
- [**FLASH**](https://ccb.jhu.edu/software/FLASH/) - [commandLine] - FLASH (Fast Length Adjustment of SHort reads) merges overlapping paired-end reads when the original DNA fragments are shorter than twice the read length. The resulting longer reads significantly improve genome assemblies and transcriptome assembly.
- [**BBDuk**](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbduk-guide/) - [commandLine] - Part of the BBMap suite, BBDuk performs quality trimming, adapter removal, k-mer-based contaminant filtering, and read normalization, filtering reads containing Illumina artifacts and known contaminants.
- [**PEAR**](https://github.com/tseemann/PEAR) - [commandLine] - PEAR (Paired-End reAd mergeR) is an ultrafast, memory-efficient tool that merges overlapping paired-end reads from Illumina platforms using dynamic programming.
- [**NGmerge**](https://github.com/noahp/NGmerge) - [commandLine] - A modern replacement for PEAR that merges overlapping paired-end reads even with substantial mismatches, offering improved accuracy for high-error-rate reads.
- [**sickle**](https://github.com/najoshi/sickle) - [commandLine] - A windowed adaptive trimming tool for FASTQ files, using a sliding window approach to trim low-quality bases from the 3' end based on average quality within the window.
- [**Skewer**](https://github.com/relipmoc/skewer) - [commandLine] - A fast and accurate adapter trimming tool for paired-end reads using dynamic programming and the Needleman-Wunsch algorithm for optimal adapter alignment.
- [**SOAPnuke**](https://github.com/BGI-flexlab/SOAPnuke) - [commandLine] - BGI-developed QC and preprocessing tool supporting adapter trimming, quality filtering, complexity filtering, and duplication removal for large-scale NGS data with multi-threading.
- [**HTStream**](https://github.com/s4hts/HTStream) - [commandLine] - A stream-based FASTQ preprocessing suite that chains adapter trimming, quality filtering, and read merging without intermediate files while maintaining pairing information.
- [**SeqFu**](https://github.com/telatin/seqfu2) - [commandLine] - A Rust-based toolkit for FASTQ/FASTA manipulation and statistics offering low memory usage and high speed for sequence counting, filtering, and quality metrics.

---

## 2. Sequence Alignment

Tools for mapping sequencing reads to a reference genome or performing multiple sequence alignment.

- [**BWA**](https://github.com/lh3/bwa) - [commandLine] - BWA (Burrows-Wheeler Aligner) maps low-divergent sequences against a large reference genome. It includes three algorithms: BWA-backtrack (for Illumina reads up to 100bp), BWA-SW (70bp-1Mbp), and BWA-MEM (70bp-10Mbp, the most widely used and recommended).
- [**BWA-MEM2**](https://github.com/bwa-mem2/bwa-mem2) - [commandLine] - An optimized version of BWA-MEM providing approximately 2x speed improvement using SSE/AVX2 vectorized instructions while producing identical output. A drop-in replacement for BWA-MEM.
- [**Bowtie 2**](https://bowtie-bio.sourceforge.net/bowtie2/) - [commandLine] - An ultrafast and memory-efficient short-read aligner using FM-index based on the Burrows-Wheeler transform. Supports both gapped and ungapped alignment for reads from 50bp to thousands of bases.
- [**STAR**](https://github.com/alexdobin/STAR) - [commandLine] - STAR (Spliced Transcripts Alignment to a Reference) is an ultrafast RNA-seq aligner using sequential maximum mappable seed search. It detects splice junctions directly from reads, and can identify chimeric and circular RNA without requiring annotation files.
- [**HISAT2**](https://daehwankimlab.github.io/hisat2/) - [commandLine] - Uses hierarchical FM-index and graph-based index for fast and accurate spliced alignment of RNA-seq reads with low memory usage.
- [**Subread**](https://subread.sourceforge.net/) - [commandLine] - A general-purpose read aligner using a novel subread indexing strategy for both DNA-seq and RNA-seq. The same package includes featureCounts, a widely used read quantification tool.
- [**BBMap**](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbmap-guide/) - [commandLine] - A short-read aligner combining k-mer and banded alignment, handling large genomes with high sensitivity. Supports PacBio, ONT, and Illumina reads, known for accuracy with many mismatches or indels.
- [**LASTZ**](https://github.com/lastz/lastz) - [commandLine] - A pairwise sequence aligner designed for large genomic sequences, primarily used for whole-genome alignment between closely related species, supporting both global and local alignment.
- [**BLAT**](https://genome.ucsc.edu/FAQ/FAQblat.html) - [webtools, commandLine] - BLAT (BLAST-Like Alignment Tool) rapidly aligns cDNA, mRNA, and protein sequences to a genome using indexed k-mer approach, typically orders of magnitude faster than BLAST for similar alignments.
- [**MAFFT**](https://mafft.cbrc.jp/alignment/software/) - [commandLine] - MAFFT (Multiple Alignment using Fast Fourier Transform) provides multiple sequence alignment with multiple refinement methods, suitable for aligning large numbers of DNA, RNA, or protein sequences.
- [**Clustal Omega**](https://www.ebi.ac.uk/Tools/msa/clustalo/) - [commandLine, webtools] - A multiple sequence alignment program using seeded guide trees and HMM profile-profile techniques, scalable to thousands of sequences for high-quality protein and DNA alignments.
- [**MUSCLE**](https://www.drive5.com/muscle/) - [commandLine] - MUSCLE (MUltiple Sequence Comparison by Log-Expectation) offers fast and accurate multiple sequence alignment using iterative refinement with an excellent speed/accuracy tradeoff.

---

## 3. Sequence Similarity Search

Tools for searching query sequences against large databases, including BLAST and its faster alternatives.

- [**NCBI BLAST+**](https://blast.ncbi.nlm.nih.gov/) - [commandLine, webtools] - The gold-standard sequence similarity search tool supporting blastn, blastp, blastx, tblastn, and tblastx. Finds local similarity regions between sequences and calculates statistical significance.
- [**DIAMOND**](https://github.com/bbuchfink/diamond) - [commandLine] - An ultra-fast protein alignment tool 100x to 10,000x faster than BLASTX while maintaining comparable sensitivity. Uses double indexing and seed extension optimization, widely used in metagenomics.
- [**MMseqs2**](https://github.com/soedinglab/MMseqs2) - [commandLine] - Ultra-fast and sensitive protein sequence search and clustering tool. 10,000x faster than BLAST for clustering and 400x faster than UBLAST for searching. Includes modules for taxonomic classification and profile-profile search.
- [**VSEARCH**](https://github.com/torognes/vsearch) - [commandLine] - A free and open-source 64-bit tool for metagenomic sequence clustering, searching, chimera detection, dereplication, and shuffling. A widely used open-source alternative to USEARCH.
- [**USEARCH**](https://www.drive5.com/usearch/) - [commandLine] - A high-performance sequence analysis tool offering sequence search, UPARSE clustering, OTU picking, and chimera detection, widely used in 16S rRNA amplicon analysis.
- [**LAST**](https://gitlab.com/mcfrith/last) - [commandLine] - A large-scale sequence alignment tool for efficiently aligning sequenced genomes, particularly good at handling frameshifts in nucleotide alignments for comparing related genomes.
- [**FASTA**](https://fasta.bioch.virginia.edu/fasta_www2/fasta_list2.shtml) - [commandLine, webtools] - A classic sequence similarity search suite predating BLAST, including FASTA (protein-protein), FASTN (nucleotide-nucleotide), and FASTX/FASTY (translated searches) using k-mer heuristic methods.

---

## 4. Genome Assembly

Tools for de novo genome assembly supporting short reads, long reads, or hybrid approaches.

- [**SPAdes**](https://github.com/ablab/spades) - [commandLine] - A de Bruijn graph genome assembler for bacterial, fungal, and small eukaryotic genomes supporting Illumina, IonTorrent, PacBio, and ONT data. Includes metaSPAdes (metagenomes), plasmidSPAdes (plasmids), and rnaSPAdes (transcriptomes).
- [**Flye**](https://github.com/fenderglass/Flye) - [commandLine] - A de novo long-read assembler for PacBio and ONT reads using an extend-repeat-merge approach without error correction. Produces high-quality assemblies and supports metagenome assembly.
- [**Canu**](https://github.com/marbl/canu) - [commandLine] - A long-read assembler for PacBio and ONT following the Celera Assembler methodology with three stages: correction (self-correction), trimming, and assembly (overlap-layout-consensus). Handles high error rates typical of long reads.
- [**Hifiasm**](https://github.com/chhylp123/hifiasm) - [commandLine] - A haplotype-resolved de novo assembler for PacBio HiFi reads. Produces highly accurate phased assemblies of diploid genomes using a graph-based approach.
- [**MEGAHIT**](https://github.com/voutcn/megahit) - [commandLine] - An ultra-fast metagenome assembler for short reads using succinct de Bruijn graphs (SdBG) for low memory usage and high speed, suitable for complex metagenomic datasets with uneven coverage.
- [**ABySS**](https://github.com/bcgsc/abyss) - [commandLine] - A de novo short-read assembler for genomes of any size using de Bruijn graphs with MPI parallel execution, designed for paired-end and mate-pair libraries.
- [**Velvet**](https://www.ebi.ac.uk/~zerbino/velvet/) - [commandLine] - A classic de novo genomic assembler using de Bruijn graphs for short reads, one of the pioneering tools that brought de Bruijn graph assembly to eukaryotic genomes.
- [**SOAPdenovo2**](https://github.com/aquaskyline/SOAPdenovo2) - [commandLine] - A short-read de novo assembler developed by BGI for large genomes with limited memory, using de Bruijn graph approach with gap closure and scaffolding modules.
- [**Trinity**](https://github.com/trinityrnaseq/trinityrnaseq) - [commandLine] - A de novo transcriptome assembler for RNA-seq data that reconstructs full-length transcripts and isoforms through three modules: Inchworm, Chrysalis, and Butterfly.
- [**Unicycler**](https://github.com/rrwick/Unicycler) - [commandLine] - A hybrid assembler for bacterial genomes combining SPAdes short-read assembly with long-read repeat resolution to produce complete circular chromosomes and plasmids.
- [**MaSuRCA**](https://github.com/alekseyzimin/masurca) - [commandLine] - A hybrid assembler generating super-reads from short reads and combining them with long reads for assembly using the CABOG assembler, suitable for both small and large genomes.
- [**Shasta**](https://github.com/paoloshasta/shasta) - [commandLine] - A fast de novo long-read assembler for ONT data using novel run-length encoding and marker representation, designed to assemble human-scale genomes with reasonable compute resources.
- [**Raven**](https://github.com/lbcb-sci/raven) - [commandLine] - A de novo long-read assembler for ONT and PacBio CLR reads using overlap-layout-consensus without error correction, designed to be simple and fast.
- [**IDBA-UD**](https://github.com/loneknightpy/idba) - [commandLine] - An iterative de Bruijn graph assembler for deep metagenomic short-read data with uneven sequencing depths, using multiple k-mer sizes iteratively to recover low-abundance species.

---

## 5. Third Generation Sequencing Analysis

Specialized tools for PacBio and Oxford Nanopore sequencing data.

- [**minimap2**](https://github.com/lh3/minimap2) - [commandLine] - A versatile sequence aligner supporting PacBio/ONT read mapping, long-read overlap detection (up to ~15% error rate), splice-aware alignment, Illumina read alignment, assembly-to-assembly alignment, and cross-species whole-genome alignment.
- [**longshot**](https://github.com/pjedge/longshot) - [commandLine] - A variant caller for diploid genomes using long error-prone reads (PacBio/ONT). Takes aligned BAM as input and outputs a phased VCF with haplotype information. Currently calls SNVs and genotypes indels from input VCF.
- [**Medaka**](https://github.com/nanoporetech/medaka) - [commandLine] - A consensus sequence generator for ONT data using neural networks to correct basecalling errors, also capable of variant calling.
- [**Guppy**](https://community.nanoporetech.com/) - [commandLine] - Oxford Nanopore's official basecaller converting raw electrical signals (FAST5) to nucleotide sequences (FASTQ), supporting GPU acceleration with models for different accuracy/speed tradeoffs.
- [**Dorado**](https://github.com/nanoporetech/dorado) - [commandLine] - ONT's next-generation basecaller replacing Guppy with improved accuracy, POD5 format support, simplex/duplex/base modification detection models.
- [**pbmm2**](https://github.com/PacificBiosciences/pbmm2) - [commandLine] - PacBio's official minimap2 wrapper for SMRT read alignment, providing consistent output formats and PacBio-optimized performance.
- [**CCS**](https://github.com/PacificBiosciences/ccs) - [commandLine] - CCS (Circular Consensus Sequencing) generates highly accurate (≥Q20, typically Q30) consensus sequences from PacBio subreads by computing consensus from multiple polymerase passes.
- [**IsoSeq**](https://github.com/PacificBiosciences/IsoSeq) - [commandLine] - PacBio's full-length transcriptome analysis pipeline processing Iso-Seq data to generate high-quality transcript sequences by clustering and polishing.
- [**Clair3**](https://github.com/HKU-BAL/Clair3) - [commandLine] - A high-accuracy deep learning variant caller for ONT reads using a two-stage pileup and full-alignment strategy, achieving accuracy comparable to short-read pipelines.
- [**NanoPlot**](https://github.com/wdecoster/NanoPlot) - [commandLine] - A plotting tool for long-read sequencing data generating visualizations of read length distributions, quality scores, yield, and time-series analysis.
- [**Pychopper**](https://github.com/nanoporetech/pychopper) - [commandLine] - Identifies, orients, and trims full-length Nanopore cDNA reads, classifying reads as full-length, truncated, or chimeric for transcriptomics preprocessing.

---

## 6. RNA-seq Analysis

Tools for transcriptome analysis including alignment, quantification, differential expression, and splicing.

- [**STAR**](https://github.com/alexdobin/STAR) - [commandLine] - An ultrafast RNA-seq aligner using sequential maximum mappable seed search in suffix arrays. Detects splice junctions, chimeric RNA, and circular RNA directly from reads.
- [**Salmon**](https://github.com/COMBINE-lab/salmon) - [commandLine] - Quantifies transcript-level expression using quasi-mapping (bypassing full alignment) with corrections for GC bias, fragment length bias, and sequence-specific bias.
- [**Kallisto**](https://github.com/pachterlab/kallisto) - [commandLine] - Uses pseudoalignment via k-mer matching for extremely fast transcript quantification, processing 30 million reads in minutes on a standard laptop with bootstrap uncertainty quantification.
- [**RSEM**](https://github.com/deweylab/RSEM) - [commandLine] - Quantifies gene and transcript abundances using EM algorithm, supporting paired-end, strand-specific, and variable-length reads with posterior probability calculations.
- [**featureCounts**](https://subread.sourceforge.net/) - [commandLine] - A highly efficient read quantification program counting reads mapped to genomic features, known for its speed and low memory usage for large-scale RNA-seq.
- [**StringTie**](https://ccb.jhu.edu/software/stringtie/) - [commandLine] - Assembles transcripts from RNA-seq alignments and estimates their abundances, identifying novel isoforms with optional reference annotation integration.
- [**DESeq2**](https://bioconductor.org/packages/DESeq2/) - [R] - R/Bioconductor package for differential expression using negative binomial GLMs with shrinkage estimation for dispersions and fold changes, supporting complex experimental designs.
- [**edgeR**](https://bioconductor.org/packages/edgeR/) - [R] - R/Bioconductor package using empirical Bayes estimation and exact tests based on the negative binomial model, particularly effective for small sample sizes.
- [**limma**](https://bioconductor.org/packages/limma/) - [R] - R package for differential expression using linear models, extended to RNA-seq via voom transformation, supporting complex designs with empirical Bayes moderation.
- [**Cufflinks**](https://github.com/cole-trapnell-lab/cufflinks) - [commandLine] - A pioneering suite for transcript assembly, quantification, and differential expression analysis including Cufflinks, Cuffcompare, and Cuffdiff.
- [**Arriba**](https://github.com/suhrig/arriba) - [commandLine] - Detects gene fusions from RNA-seq using STAR aligner output, with extensive filtering to remove artifacts and visualization for validation.
- [**STAR-Fusion**](https://github.com/STAR-Fusion/STAR-Fusion) - [commandLine] - A comprehensive pipeline for detecting gene fusions from RNA-seq using STAR aligner, widely used in cancer genomics.
- [**rMATS**](https://rnaseq-mats.sourceforge.net/) - [commandLine, webtools] - Detects differential alternative splicing events (exon skipping, intron retention, alternative splice sites, mutually exclusive exons) using a statistical model.
- [**SUPPA2**](https://github.com/comprna/SUPPA2) - [commandLine] - A fast alternative splicing analysis tool calculating PSI values and performing differential splicing analysis, integrating with multiple quantification tools.

---

## 7. Variant Calling & Annotation

Tools for detecting SNVs, indels, and annotating their functional effects.

- [**samtools**](https://github.com/samtools/samtools) - [commandLine] - A suite for working with SAM/BAM/CRAM files including sorting, indexing, merging, filtering, and the mpileup command for generating genotype likelihoods.
- [**bcftools**](https://github.com/samtools/bcftools) - [commandLine] - A set of utilities for variant calling and VCF/BCF manipulation, including SNP/indel calling, filtering, comparison, and format conversion.
- [**GATK**](https://gatk.broadinstitute.org/) - [commandLine] - The industry-standard variant discovery toolkit including HaplotypeCaller (germline), Mutect2 (somatic), joint genotyping, VQSR filtering, and CNV analysis.
- [**DeepVariant**](https://github.com/google/deepvariant) - [commandLine] - A deep learning variant caller using convolutional neural networks on pileup images to achieve high accuracy without manual quality filtering.
- [**freebayes**](https://github.com/freebayes/freebayes) - [commandLine] - A Bayesian haplotype-based variant caller detecting SNPs, indels, MNPs, and complex events, suitable for high-ploidy or pooled samples.
- [**Strelka2**](https://github.com/Illumina/strelka) - [commandLine] - A fast and accurate germline and somatic variant caller from Illumina using haplotype-based small variant detection with empirical Bayes methods.
- [**Mutect2**](https://gatk.broadinstitute.org/) - [commandLine] - GATK's somatic variant detection tool for tumor-normal pairs, using a probabilistic model accounting for purity, errors, and contamination.
- [**VarScan2**](https://github.com/dkoboldt/varscan) - [commandLine] - Detects somatic mutations, CNVs, and LOH from tumor-normal paired sequencing data using heuristic and statistical methods.
- [**VEP**](https://www.ensembl.org/info/docs/tools/vep/) - [commandLine, webtools] - Ensembl's Variant Effect Predictor determines functional effects of variants on genes, transcripts, and proteins using Ensembl annotation.
- [**ANNOVAR**](https://annovar.openbioinformatics.org/) - [commandLine] - Annotates genetic variants with gene-based, region-based, and filter-based annotation supporting RefSeq, UCSC, and ENCODE databases.
- [**SnpEff**](https://pcingola.github.io/SnpEff/) - [commandLine] - A fast variant annotation tool classifying effects as high, moderate, low, or modifier impact, tightly integrated with Galaxy.
- [**Pindel**](https://github.com/genome/pindel) - [commandLine] - Detects breakpoints of large deletions and medium insertions using pattern growth from paired-end reads, effective for mobile element insertions.

---

## 8. CRISPR Targeted Deep Sequencing Analysis

Tools for analyzing CRISPR-Cas9 genome editing outcomes from deep and Sanger sequencing.

- [**CRISPResso2**](https://github.com/pinellolab/CRISPResso2) - [commandLine, webtools] - A comprehensive suite for analyzing genome editing deep sequencing data. Aligns reads to reference, quantifies indels/mutations, classifies modified vs. unmodified reads, and generates intuitive plots. Supports single/paired-end reads, multiple gRNAs, and pooled screens.
- [**CrispRVariants**](https://bioconductor.org/packages/CrispRVariants/) - [R] - R/Bioconductor package for analyzing CRISPR-Cas9 mutagenesis experiments. Localizes variant alleles relative to the Cas9 cut site, plots allele combinations, and calculates mutation rates with flexible filtering.
- [**CRISP-ID**](http://crispid.gbiomed.kuleuven.be) - [webtools] - A web application detecting exact indel size and location from Sanger sequencing, identifying mono-allelic and bi-allelic editing events.
- [**Microhomology-Predictor**](http://www.rgenome.net/mich-calculator/) - [webtools] - A web tool calculating microhomology-associated scores for ZFNs, TALENs, and Cas9 to predict MMEJ repair outcomes.
- [**TIDE**](https://tide.nki.nl/) - [webtools] - Quantifies CRISPR editing efficiency from Sanger sequencing by decomposing chromatograms to identify indel frequencies without needing deep sequencing.
- [**ICE**](https://ice.synthego.com/) - [webtools] - Synthego's web tool for analyzing CRISPR editing outcomes from Sanger data, providing editing efficiency, common outcomes, and knockout scores.
- [**CHOPCHOP**](https://chopchop.cbu.uib.no/) - [webtools] - Designs guide RNAs for CRISPR/Cas9/Cas12a, predicting on-target efficiency and off-target specificity with genomic visualization.
- [**CRISPick**](https://portals.broadinstitute.org/gppx/crispick) - [webtools] - Broad Institute's gRNA design tool combining efficiency and specificity scores, integrated with genome-wide CRISPR knockout libraries.

---

## 9. GWAS (Genome-Wide Association Studies)

Statistical tools for GWAS, fine-mapping, and polygenic risk scoring.

- [**PLINK 1.9**](https://www.cog-genomics.org/plink/) - [commandLine] - A free open-source whole-genome association analysis toolset for data management, QC, population stratification, association testing, and LD calculations.
- [**PLINK 2.0**](https://www.cog-genomics.org/plink/2.0/) - [commandLine] - Next-generation PLINK with multi-threading, PGEN format support, and improved performance for biobank-scale datasets.
- [**BOLT-LMM**](https://alkesgroup.broadinstitute.org/BOLT-LMM/) - [commandLine] - A Bayesian mixed-model association method for complex traits, using priors on effect sizes to increase power while controlling for population structure.
- [**SAIGE**](https://github.com/weizhouUMICH/SAIGE) - [commandLine] - A scalable mixed-model method for binary and quantitative traits using saddlepoint approximation for accurate p-values even with extreme case-control imbalance.
- [**GEMMA**](https://github.com/genetics-statistics/GEMMA) - [commandLine] - Performs genome-wide association analysis using standard linear mixed models, estimating variance components and running BSLMM analysis.
- [**GCTA**](https://yanglab.westlake.edu.cn/software/gcta/) - [commandLine] - Estimates SNP heritability, genetic correlation, and performs conditional/joint GWAS analysis (COJO) and multi-trait analysis.
- [**MAGMA**](https://cncr.nl/research/magma/) - [commandLine, webtools] - Aggregates SNP associations to the gene level for gene-based analysis, accounting for LD and incorporating functional annotations.
- [**FINEMAP**](http://www.christianbenner.com/) - [commandLine] - Bayesian fine-mapping identifying causal variants from GWAS summary statistics using stochastic search with posterior probability computation.
- [**SuSiE**](https://github.com/stephenslab/susieR) - [R] - R package for fine-mapping using the "Sum of Single Effects" regression model to identify credible sets of causal variants.
- [**FUMA**](https://fuma.ctglab.nl/) - [webtools] - A web platform for functional annotation of GWAS results, identifying independent significant SNPs, mapping to genes, and performing enrichment analysis.
- [**Michigan Imputation Server**](https://imputationserver.sph.umich.edu/) - [webtools] - A free cloud-based genotype imputation service using Minimac4 with HRC, TOPMed, and 1000 Genomes reference panels.
- [**Minimac4**](https://github.com/statgen/Minimac4) - [commandLine] - A fast genotype imputation tool using hidden Markov models, highly optimized for large-scale imputation with multi-threading.
- [**Beagle 5**](https://faculty.washington.edu/browning/beagle/beagle.html) - [commandLine] - Java application for phasing, imputation, and association analysis using graphical models of haplotype clusters.
- [**PRSice-2**](https://github.com/choishingwan/PRSice) - [commandLine] - Computes polygenic risk scores by weighting risk alleles by GWAS effect sizes, with automatic p-value threshold optimization and plotting.
- [**LDpred2**](https://github.com/privefl/bigsnpr) - [R] - R package for Bayesian PRS construction accounting for LD between SNPs, with multiple models for improved prediction accuracy.

---

## 10. Population Genetics

Software for coalescent simulation, demographic inference, selection scans, and population structure analysis.

- [**msprime**](https://github.com/tskit-dev/msprime) - [commandLine] - A scalable coalescent simulator using tree sequence data structures for efficient simulation of genome-wide variation with recombination.
- [**SLiM**](https://github.com/MesserLab/SLiM) - [commandLine] - A forward-in-time evolutionary simulator modeling selection, mutation, recombination, and population structure changes with detailed genetic architecture specification.
- [**Fastsimcoal2**](http://cmpg.unibe.ch/software/fastsimcoal2/) - [commandLine] - Simulates genetic diversity under complex demographic scenarios and infers parameters from the SFS using maximum-likelihood estimation.
- [**Dadi**](https://github.com/rajanil/dadi) - [commandLine] - Uses diffusion approximation to the allele frequency spectrum for inferring demographic history including divergence, migration, and population size changes.
- [**PSMC**](https://github.com/lh3/psmc) - [commandLine] - Infers historical effective population size from a single diploid genome using a hidden Markov model, widely used for ancient population dynamics.
- [**MSMC2**](https://github.com/stschiff/msmc2) - [commandLine] - Extends PSMC to multiple phased genomes for higher-resolution inference of population size changes and split times.
- [**SMC++**](https://github.com/terhorst/smcpp) - [commandLine] - Estimates population size history using combined SFS and LD information with composite likelihood, effective for recent history from large samples.
- [**ADMIXTURE**](https://dalexander.github.io/admixture/) - [commandLine] - Maximum-likelihood estimation of individual ancestry proportions from K ancestral populations, computationally efficient for large-scale data.
- [**STRUCTURE**](https://web.stanford.edu/group/pritchardlab/structure.html) - [commandLine] - A foundational Bayesian clustering method for inferring population structure and detecting admixture from genotypic data.
- [**EIGENSOFT**](https://github.com/DReichLab/EIG) - [commandLine] - Includes smartpca for PCA-based population structure detection, Tracy-Widom test, and EIGENSTRAT for association test stratification correction.
- [**TREEMIX**](https://github.com/joepickrell/treemix) - [commandLine] - Infers population splits and admixture events from allele frequency data, building trees with migration edges.
- [**AdmixTools**](https://github.com/DReichLab/AdmixTools) - [commandLine] - Analyzes admixture using f-statistics (f2, f3, f4), D-statistics, and qpWave/qpAdm for quantifying ancestry proportions.
- [**VCFtools**](https://vcftools.github.io/) - [commandLine] - A program package for VCF filtering, comparison, and popgen calculations including Fst, Tajima's D, nucleotide diversity, and LD estimation.
- [**ANGSD**](https://github.com/angsd/angsd) - [commandLine] - Analyzes NGS data for population genetics using genotype likelihoods directly, estimating allele frequencies, Fst, and IBS matrices.
- [**SweepFinder2**](https://github.com/3hunter/SweepFinder2) - [commandLine] - Detects selective sweeps using composite likelihood ratio tests without requiring knowledge of the selected allele.
- [**Selscan**](https://github.com/szpiech/selscan) - [commandLine] - Computes iHS, XP-EHH, and XP-CLR statistics for detecting recent positive selection from phased genotype data.
- [**rehh**](https://cran.r-project.org/web/packages/rehh/) - [R] - R package for detecting selection from EHH statistics, providing iHS and Rsb calculation with genome-wide scan visualization.
- [**pixy**](https://github.com/ksamuk/pixy) - [commandLine] - Calculates nucleotide diversity (pi) and divergence (dxy) from VCF files in sliding windows using a genotype likelihood approach.

---

## 11. Structural Variation Detection

Specialized tools for detecting structural variants including deletions, inversions, and translocations.

- [**Snowman**](https://github.com/broadinstitute/SnowmanSV) - [commandLine] - Detects SVs using genome-wide local assembly around regions of interest for high-accuracy breakpoint identification.
- [**seeksv**](https://github.com/qiukunlong/seeksv) - [commandLine] - Accurate SV and virus integration detection using discordant read pairs and split-reads for both germline and somatic variants.
- [**Genomon SV**](https://github.com/Genomon-Project/GenomonSV) - [commandLine] - Detects somatic SVs from cancer genomes including deletions, inversions, and translocations using paired-end mapping patterns.
- [**BreakDancer**](https://github.com/genome/breakdancer) - [commandLine] - Genome-wide SV detection using read pairs with unexpected separation distances or orientation, predicting five types of SVs.
- [**Manta**](https://github.com/Illumina/manta) - [commandLine] - Calls SVs and medium indels from paired-end sequencing, optimized for rapid analysis of large cohorts.
- [**GRIDSS**](https://github.com/PapenfussLab/gridss) - [commandLine] - Uses local assembly and machine learning to detect SV breakpoints from short reads, distinguishing real SVs from artifacts.
- [**LUMPY**](https://github.com/arq5x/lumpy-sv) - [commandLine] - A probabilistic framework integrating discordant read pairs, split-reads, read depth, and prior knowledge for SV discovery.
- [**Wham**](https://github.com/zeeev/wham) - [commandLine] - Clusters read mapping signatures to identify SVs without requiring a pre-specified SV model.
- [**SVision**](https://github.com/xjtu-omics/SVision) - [commandLine] - Deep learning-based SV detection from long-read alignment patterns with intuitive visualization.
- [**CuteSV**](https://github.com/tjiangHIT/cuteSV) - [commandLine] - A sensitive and fast long-read SV caller using clustering and ensemble learning, supporting PacBio and ONT data.

---

## 12. Metagenomics

Tools for taxonomic and functional analysis of microbial communities.

- [**MetaPhlAn 4**](https://huttenhower.sph.harvard.edu/metaphlan/) - [commandLine] - Strain-level microbial profiling using clade-specific marker genes to estimate species and strain abundance from metagenomes.
- [**Kraken 2**](https://ccb.jhu.edu/software/kraken2/) - [commandLine] - Fast taxonomic classification by matching k-mers against a compressed microbial genome database using FM-index.
- [**Bracken**](https://ccb.jhu.edu/software/bracken/) - [commandLine] - Bayesian species-level abundance estimation from Kraken results, correcting for misclassification between similar genomes.
- [**HUMAnN 3**](https://huttenhower.sph.harvard.edu/humann/) - [commandLine] - Functional profiling estimating gene family and metabolic pathway abundances from metagenomes or metatranscriptomes.
- [**MEGAN6**](https://ab.inf.uni-tuebingen.de/software/megan6/) - [commandLine] - Interactive taxonomic and functional analysis using LCA algorithm with extensive visualization capabilities.
- [**Centrifuge**](https://ccb.jhu.edu/software/centrifuge/) - [commandLine] - Rapid FM-index-based taxonomic classifier for large metagenomic datasets with low memory usage.
- [**MetaBAT 2**](https://bitbucket.org/berkeleylab/metabat) - [commandLine] - Bins metagenomic contigs into draft genomes using tetranucleotide frequency and abundance profiles.
- [**CONCOCT**](https://github.com/BinPro/CONCOCT) - [commandLine] - Unsupervised binning of metagenomic contigs using Gaussian mixture models on k-mer composition and coverage.
- [**MaxBin 2**](https://sourceforge.net/projects/maxbin2/) - [commandLine] - Automated EM-based binning of metagenomic sequences into population genomes.
- [**CheckM**](https://github.com/Ecogenomics/CheckM) - [commandLine] - Assesses MAG quality (completeness and contamination) using lineage-specific marker gene sets.
- [**GTDB-Tk**](https://ecogenomics.github.io/GTDBTk/) - [commandLine] - Assigns taxonomic classifications to bacterial and archaeal genomes using the GTDB reference database.
- [**Mash**](https://mash.readthedocs.io/) - [commandLine] - Fast genome distance estimation using MinHash sketching for genome clustering and database searching.
- [**Pavian**](https://github.com/fbreitwieser/pavian) - [webtools] - Interactive web-based visualization of metagenomic classification results with sankey diagrams and heatmaps.
- [**QIIME 2**](https://qiime2.org/) - [commandLine] - A next-generation microbiome bioinformatics platform supporting amplicon and shotgun metagenomics with end-to-end analysis.

---

## 13. Epigenomics

Tools for analyzing DNA methylation, histone modifications, and chromatin accessibility.

- [**MACS2/3**](https://github.com/macs3-project/MACS) - [commandLine] - Identifies enriched peaks from ChIP-seq, ATAC-seq, and CUT&Tag data using a Poisson-based model relative to controls.
- [**Bismark**](https://www.bioinformatics.babraham.ac.uk/projects/bismark/) - [commandLine] - Maps bisulfite-treated reads and performs single-base resolution DNA methylation calling for WGBS and RRBS.
- [**MethylKit**](https://github.com/al2na/methylKit) - [R] - R/Bioconductor package for differential DNA methylation analysis, annotation, and visualization from bisulfite sequencing.
- [**MethylDackel**](https://github.com/dpryan79/MethylDackel) - [commandLine] - Rapid methylation extraction from Bismark-aligned reads, providing per-CpG methylation metrics.
- [**deepTools**](https://deeptools.readthedocs.io/) - [commandLine] - A Python suite for normalization, correlation, and visualization of ChIP-seq, ATAC-seq, and bisulfite-seq data.
- [**HOMER**](https://homer.ucsd.edu/homer/) - [commandLine] - A suite for peak finding, de novo motif discovery, and annotation for ChIP-seq, ATAC-seq, and Hi-C data.
- [**ChIPseeker**](https://github.com/YuLab-SMU/ChIPseeker) - [R] - R/Bioconductor package for annotating ChIP-seq peaks to genomic features and visualizing coverage profiles.
- [**bwa-meth**](https://github.com/brentp/bwa-meth) - [commandLine] - A specialized aligner for WGBS using modified BWA-MEM to handle C-to-T converted reads.
- [**GemBS**](https://github.com/heathsc/gemBS) - [commandLine] - A comprehensive bisulfite sequencing analysis platform integrating bwa-meth and MethylDackel with QC reporting.
- [**Picard Tools**](https://broadinstitute.github.io/picard/) - [commandLine] - Java tools for marking duplicates, alignment metrics, and quality score recalibration essential in epigenomics workflows.

---

## 14. Proteomics

Software for mass spectrometry-based proteomics.

- [**MaxQuant**](https://www.maxquant.org/) - [commandLine] - Quantitative proteomics software supporting LFQ, SILAC, TMT, and iTRAQ with Andromeda search engine and posterior error probabilities.
- [**FragPipe**](https://fragpipe.nesvilab.org/) - [commandLine] - A complete proteomics pipeline integrating MSFragger identification with quantification and validation, supporting open search and glycomics.
- [**MSFragger**](https://msfragger.nesvilab.org/) - [commandLine] - Ultrafast peptide identification using fragment-ion indexing with open-mass searching for comprehensive PTM identification.
- [**OpenMS**](https://www.openms.de/) - [commandLine] - An open-source C++ framework providing tools for LC-MS data analysis, identification, quantification, and custom pipeline construction.
- [**Skyline**](https://skyline.ms/) - [commandLine] - A Windows application for targeted proteomics method creation and analysis of SRM, PRM, and DIA data.
- [**DIA-NN**](https://github.com/vdemichev/DiaNN) - [commandLine] - Deep neural network-based DIA proteomics analysis with high speed and accuracy, supporting library-based and library-free approaches.
- [**ProteoWizard**](https://proteowizard.sourceforge.io/) - [commandLine] - Converts vendor-specific MS formats to open standards (mzML, mzXML) and processes mass spectrometry data.
- [**Comet**](https://uwpr.github.io/Comet/) - [commandLine] - An open-source tandem MS search engine supporting multiple fragmentation types with cross-correlation scoring.
- [**X!Tandem**](https://www.thegpm.org/tandem/) - [commandLine] - A popular open-source MS/MS search engine with parallel processing, PTM support, and post-analysis tools.
- [**Sage**](https://github.com/lazear/sage) - [commandLine] - A fast open-source peptide search engine with GPU acceleration, written in Rust for performance and memory safety.
- [**pFind**](http://pfind.ict.ac.cn/software/pFind/) - [commandLine] - A high-performance database search tool using open-search strategy for comprehensive peptide identification.

---

## 15. Single-Cell Analysis

Tools for single-cell RNA-seq, ATAC-seq, and spatial transcriptomics analysis.

- [**Seurat**](https://satijalab.org/seurat/) - [R] - An R package for scRNA-seq analysis including QC, normalization, clustering, differential expression, multi-sample integration, and spatial transcriptomics.
- [**Scanpy**](https://scanpy.readthedocs.io/) - [commandLine] - A scalable Python toolkit for single-cell gene expression analysis designed for large-scale datasets with preprocessing, clustering, and trajectory inference.
- [**Cell Ranger**](https://support.10xgenomics.com/) - [commandLine] - 10x Genomics' official pipeline for processing single-cell RNA-seq and ATAC-seq data from demultiplexing to gene expression matrices.
- [**Monocle 3**](https://cole-trapnell-lab.github.io/monocle3/) - [R] - R package for single-cell trajectory inference using reversed graph embedding, clustering, and differential expression.
- [**scvi-tools**](https://scvi-tools.org/) - [commandLine] - A Python framework for deep generative modeling of single-cell omics including scVI, scANVI, and TotalVI.
- [**Harmony**](https://portals.broadinstitute.org/harmony/) - [commandLine] - Fast single-cell dataset integration by iterative PCA-space alignment, correcting batch effects while preserving biological variation.
- [**scVelo**](https://scvelo.readthedocs.io/) - [commandLine] - RNA velocity estimation using a dynamical model of splicing kinetics without requiring steady-state assumptions.
- [**CellChat**](https://github.com/sqjin/CellChat) - [R] - R package for inferring cell-cell communication networks from single-cell data using curated ligand-receptor interactions.
- [**Squidpy**](https://squidpy.readthedocs.io/) - [commandLine] - A Python framework for spatial single-cell analysis integrating transcriptomics with spatial imaging data.

---

## 16. Phylogenetics

Tools for constructing and analyzing phylogenetic trees from molecular sequences.

- [**IQ-TREE 2**](http://www.iqtree.org/) - [commandLine] - Maximum-likelihood phylogenetics with ultrafast bootstrap, automatic ModelFinder, and partition model testing for both nucleotide and protein data.
- [**RAxML-NG**](https://github.com/amkozlov/raxml-ng) - [commandLine] - The next-generation RAxML for maximum-likelihood phylogenetics, optimized for large datasets with rapid bootstrap analysis.
- [**FastTree 2**](http://www.microbesonline.org/fasttree/) - [commandLine] - Builds approximate maximum-likelihood trees from very large alignments (hundreds of thousands of sequences) using heuristics.
- [**MrBayes**](https://nbisweden.github.io/MrBayes/) - [commandLine] - Bayesian phylogenetic inference using MCMC, supporting mixed models and posterior probability estimation.
- [**BEAST 2**](https://www.beast2.org/) - [commandLine] - Bayesian evolutionary analysis for time-calibrated phylogenies, divergence time estimation, and population dynamics.
- [**MEGA 11**](https://megasoftware.net/) - [commandLine] - An integrated desktop software for sequence alignment, phylogenetic tree construction, and molecular evolution analysis with a user-friendly GUI.
- [**APE**](https://cran.r-project.org/package=ape) - [R] - A foundational R package providing comprehensive phylogenetic computation including tree manipulation and comparative methods.
- [**ggtree**](https://bioconductor.org/packages/ggtree/) - [R] - R/Bioconductor package for phylogenetic tree visualization using ggplot2 grammar with annotation support.
- [**TrimAl**](https://github.com/scapella/trimal) - [commandLine] - Automatically removes poorly aligned regions from multiple sequence alignments before phylogenetic analysis.

---

## 17. Protein Structure & Molecular Dynamics

Tools for protein structure prediction and molecular simulation.

- [**AlphaFold2**](https://github.com/google-deepmind/alphafold) - [commandLine] - DeepMind's revolutionary deep learning model predicting protein 3D structures from sequences with near-experimental accuracy using a transformer architecture.
- [**AlphaFold3**](https://github.com/google-deepmind/alphafold3) - [commandLine] - Extends structure prediction to protein-ligand, protein-DNA, and protein-RNA complexes using a diffusion-based architecture.
- [**ColabFold**](https://github.com/sokrypton/ColabFold) - [webtools] - Fast AlphaFold2 implementation via Google Colab using MMseqs2 for rapid MSA generation without extensive computational resources.
- [**ESMFold**](https://github.com/facebookresearch/esm) - [commandLine] - Meta's protein structure prediction model ~60x faster than AlphaFold2 using a language model approach trained on millions of sequences.
- [**GROMACS**](https://www.gromacs.org/) - [commandLine] - A high-performance molecular dynamics engine for simulating proteins, lipids, and nucleic acids with GPU acceleration and parallel computing.
- [**AMBER**](https://ambermd.org/) - [commandLine] - A molecular dynamics simulation suite with force fields for proteins, nucleic acids, and small molecules.
- [**OpenMM**](https://openmm.org/) - [commandLine] - A GPU-accelerated molecular simulation toolkit with Python API for building and running custom simulations.
- [**AutoDock Vina**](https://vina.scripps.edu/) - [commandLine] - An open-source molecular docking tool for virtual screening, predicting how small molecules bind to protein receptors.
- [**PyMOL**](https://pymol.org/) - [commandLine] - A molecular visualization system producing high-quality 3D images of biomolecular structures with ray tracing and scripting.
- [**Rosetta**](https://www.rosettacommons.org/) - [commandLine] - A comprehensive suite for protein structure prediction, design, protein-protein docking, and loop modeling.

---

## 18. Hi-C & Chromatin Conformation

Tools for analyzing 3D chromatin architecture from Hi-C data.

- [**Juicer**](https://github.com/aidenlab/juicer) - [commandLine] - A comprehensive Hi-C data processing pipeline from alignment to contact matrix generation, producing .hic files for Juicebox visualization.
- [**Juicebox**](https://github.com/aidenlab/Juicebox) - [commandLine] - Interactive visualization and analysis of Hi-C contact maps as 2D heatmaps with annotation tracks and condition comparison.
- [**HiC-Pro**](https://github.com/nservant/HiC-Pro) - [commandLine] - An optimized Hi-C processing pipeline with alignment, filtering, contact map generation, and comprehensive QC metrics.
- [**Cooler**](https://github.com/open2c/cooler) - [commandLine] - Sparse compressed storage format and tools for Hi-C contact matrices, supporting multi-resolution analysis and efficient operations.
- [**HiGlass**](https://github.com/higlass/higlass) - [webtools] - A web-based multi-resolution viewer for Hi-C contact maps and genomic data with side-by-side comparison and track layering.
- [**CALDER**](https://github.com/CSOgroup/CALDER2) - [commandLine] - Computes chromatin compartment domains and TAD-like domains using graph-based approaches for hierarchical organization analysis.
- [**FAN-C**](https://github.com/vaquerizaslab/fanc) - [commandLine] - A command-line Hi-C analysis and visualization tool with matrix operations, normalization, and TAD calling.

---

## 19. General-Purpose & Workflow Management

Cross-domain bioinformatics utilities and workflow management systems.

- [**Galaxy**](https://galaxyproject.org/) - [webtools] - An open web-based platform for running bioinformatics tools, creating reproducible workflows, and sharing analyses.
- [**Bioconda**](https://bioconda.github.io/) - [commandLine] - A bioinformatics software distribution channel via Conda, providing dependency resolution across Linux, macOS, and Windows.
- [**Bioconductor**](https://www.bioconductor.org/) - [R] - An R-based repository of curated bioinformatics packages for high-throughput genomic data analysis.
- [**Biopython**](https://biopython.org/) - [commandLine] - Python tools for biological computation including sequence analysis, file parsing, and database access.
- [**Nextflow**](https://www.nextflow.io/) - [commandLine] - A workflow management system for scalable, reproducible bioinformatics pipelines with container and cloud/HPC support.
- [**Snakemake**](https://snakemake.readthedocs.io/) - [commandLine] - A Python-based workflow manager for reproducible data analysis with automatic job scheduling.
- [**nf-core**](https://nf-co.re/) - [commandLine] - A community-curated collection of best-practice Nextflow pipelines for standard bioinformatics analyses.
- [**Apptainer/Singularity**](https://apptainer.org/) - [commandLine] - A container platform for HPC environments providing portable, reproducible bioinformatics environments.
- [**IGV**](https://software.broadinstitute.org/software/igv/) - [commandLine] - A high-performance genomic data viewer supporting BAM, VCF, BigWig, and many other formats.
- [**UCSC Genome Browser**](https://genome.ucsc.edu/) - [webtools] - A web-based genome browser with curated assemblies, annotations, and custom track support.
- [**bedtools**](https://github.com/arq5x/bedtools2) - [commandLine] - A powerful suite of genomic interval arithmetic tools for intersect, merge, complement, and closest operations.
- [**htslib/tabix**](https://github.com/samtools/htslib) - [commandLine] - A C library for HTS data formats with tabix for indexing and rapid position-based retrieval.
- [**pysam**](https://github.com/pysam-developers/pysam) - [commandLine] - A Python wrapper for htslib providing access to SAM/BAM/CRAM and VCF/BCF files.
- [**SeqKit**](https://github.com/shenwei356/seqkit) - [commandLine] - A cross-platform, ultrafast toolkit for FASTA/Q manipulation with 30+ commands written in Go.

---

## License

MIT License

*Last updated: June 2026*
