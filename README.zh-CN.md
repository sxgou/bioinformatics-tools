# 生物信息学软件合集

![Bioinformatics](https://img.shields.io/badge/Bioinformatics-Tools-blue)
![PRs Welcome](https://img.shields.io/badge/PRs-welcome-brightgreen)

> 整理生物信息学各领域常用软件工具，按应用领域分类整理。

[**English version → README.md**](./README.md)

## 目录

- [1. 测序数据预处理](#1-测序数据预处理)
- [2. 序列比对](#2-序列比对)
- [3. 序列相似性搜索](#3-序列相似性搜索)
- [4. 基因组组装](#4-基因组组装)
- [5. 三代测序分析](#5-三代测序分析)
- [6. RNA-seq 分析](#6-rna-seq-分析)
- [7. 变异检测与注释](#7-变异检测与注释)
- [8. CRISPR 靶向深度测序分析](#8-crispr-靶向深度测序分析)
- [9. GWAS 全基因组关联分析](#9-gwas-全基因组关联分析)
- [10. 群体遗传学](#10-群体遗传学)
- [11. 结构变异检测](#11-结构变异检测)
- [12. 宏基因组学](#12-宏基因组学)
- [13. 表观基因组学](#13-表观基因组学)
- [14. 蛋白质组学](#14-蛋白质组学)
- [15. 单细胞组学分析](#15-单细胞组学分析)
- [16. 系统发育学](#16-系统发育学)
- [17. 蛋白结构与分子模拟](#17-蛋白结构与分子模拟)
- [18. Hi-C 与染色质构象分析](#18-hi-c-与染色质构象分析)
- [19. 通用工具与工作流](#19-通用工具与工作流)

---

## 1. 测序数据预处理

原始测序数据的质量控制、接头修剪、双端 reads 合并和过滤工具。

- [**FastQC**](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) - [commandLine] - 行业标准原始测序数据质控工具，生成模块化 HTML 报告，覆盖单碱基质量、GC 含量、N 含量、序列长度分布、过代表序列和接头含量等指标，帮助在下游分析前识别数据问题。
- [**MultiQC**](https://github.com/ewels/MultiQC) - [commandLine] - 将多个样本的生物信息学分析结果聚合为单个交互式 HTML 报告，自动搜索目录中的分析日志并编译为带交互式图表的汇总，支持 FastQC、STAR、HISAT2、Salmon 等多种工具。
- [**fastp**](https://github.com/OpenGene/fastp) - [commandLine] - 一站式 FASTQ 预处理工具，单次运行即可完成质控、接头修剪、质量过滤、逐 reads 质量修剪、poly-G 修剪等操作，并生成全面的 HTML/JSON 报告，替代多个传统工具的串联使用。
- [**Trimmomatic**](https://github.com/usadellab/Trimmomatic) - [commandLine] - 灵活的 Illumina NGS 数据修剪工具，支持接头裁剪、滑动窗口质量修剪、最小长度过滤和首尾低质量碱基去除，兼容单端和双端数据。
- [**Cutadapt**](https://cutadapt.readthedocs.io/) - [commandLine] - 从高通量测序 reads 中查找并去除接头序列、引物、poly-A 尾等不需要的序列，支持单端和双端数据，具有容错的接头匹配能力。
- [**FLASH**](https://ccb.jhu.edu/software/FLASH/) - [commandLine] - FLASH（Fast Length Adjustment of SHort reads）将重叠的双端 reads 快速准确地合并为完整序列，适用于原始 DNA 片段长度小于两倍读长的情况，合并后的更长 reads 可显著改善基因组和转录组组装。
- [**BBDuk**](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbduk-guide/) - [commandLine] - BBMap 套件组件，执行质量修剪、接头去除、基于 k-mer 的污染过滤和 reads 归一化，可过滤 Illumina 人工产物和已知污染物。
- [**PEAR**](https://github.com/tseemann/PEAR) - [commandLine] - 超快速、内存高效的双端 reads 合并工具，使用动态规划方法处理 Illumina 测序的重叠双端 reads。
- [**NGmerge**](https://github.com/noahp/NGmerge) - [commandLine] - PEAR 的现代替代方案，即使 reads 间存在大量错配也能有效合并，对高错误率 reads 具有改进的准确性。
- [**sickle**](https://github.com/najoshi/sickle) - [commandLine] - 基于滑动窗口的 FASTQ 自适应修剪工具，根据窗口内的平均质量从 3' 端修剪低质量碱基。
- [**Skewer**](https://github.com/relipmoc/skewer) - [commandLine] - 使用动态规划和 Needleman-Wunsch 算法快速准确地修剪双端接头。
- [**SOAPnuke**](https://github.com/BGI-flexlab/SOAPnuke) - [commandLine] - BGI 开发的 QC 和预处理工具，支持接头修剪、质量过滤、复杂度过滤和去重复，支持多线程大规模 NGS 数据处理。
- [**HTStream**](https://github.com/s4hts/HTStream) - [commandLine] - 基于流的 FASTQ 预处理套件，链式处理接头修剪、质量过滤和 reads 合并，无需中间文件同时保持配对信息。
- [**SeqFu**](https://github.com/telatin/seqfu2) - [commandLine] - Rust 编写的高性能 FASTQ/FASTA 操作和统计工具包，内存占用低、速度快，支持序列计数、过滤和质量指标计算。

---

## 2. 序列比对

将测序 reads 比对到参考基因组或进行多序列比对的工具。

- [**BWA**](https://github.com/lh3/bwa) - [commandLine] - 将低歧异度序列比对到大型参考基因组的软件包，包含三种算法：BWA-backtrack（≤100bp Illumina reads）、BWA-SW（70bp-1Mbp）和 BWA-MEM（70bp-10Mbp，最常用推荐）。
- [**BWA-MEM2**](https://github.com/bwa-mem2/bwa-mem2) - [commandLine] - BWA-MEM 的优化版本，使用 SSE/AVX2 向量化指令实现约 2 倍速度提升，输出完全一致，可直接替换 BWA-MEM。
- [**Bowtie 2**](https://bowtie-bio.sourceforge.net/bowtie2/) - [commandLine] - 基于 FM-index（Burrows-Wheeler 变换）的超快速、内存高效短读长比对器，支持有空位和无空位比对，适用于 50bp 到数千 bp 的 reads。
- [**STAR**](https://github.com/alexdobin/STAR) - [commandLine] - 超快速 RNA-seq 比对器，使用未压缩后缀数组中的顺序最大可映射种子搜索，无需注释文件即可直接从 reads 检测剪接位点、嵌合 RNA 和环状 RNA。
- [**HISAT2**](https://daehwankimlab.github.io/hisat2/) - [commandLine] - 使用层级 FM-index 和基于图的索引快速准确地进行 RNA-seq 剪接比对，内存占用低。
- [**Subread**](https://subread.sourceforge.net/) - [commandLine] - 使用创新 subread 索引策略的通用 reads 比对器，支持 DNA-seq 和 RNA-seq，同包还包括广泛使用的 featureCounts 定量工具。
- [**BBMap**](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbmap-guide/) - [commandLine] - 结合 k-mer 和带状比对的短读长比对器，对大基因组保持高灵敏度，支持 PacBio/ONT/Illumina，以处理大量错配和 indels 的高精度著称。
- [**LASTZ**](https://github.com/lastz/lastz) - [commandLine] - 大规模基因组序列成对比对工具，主要用于近缘物种间全基因组比对，支持全局和局部比对模式。
- [**BLAT**](https://genome.ucsc.edu/FAQ/FAQblat.html) - [webtools, commandLine] - 快速将 cDNA、mRNA 和蛋白质序列比对到基因组的工具，使用索引 k-mer 方法，比 BLAST 快数倍到数十倍，广泛用于 UCSC 基因组浏览器。
- [**MAFFT**](https://mafft.cbrc.jp/alignment/software/) - [commandLine] - 基于 FFT 的多序列比对程序，提供多种精炼方法，支持大量 DNA/RNA/蛋白质序列的比对。
- [**Clustal Omega**](https://www.ebi.ac.uk/Tools/msa/clustalo/) - [commandLine, webtools] - 使用种子引导树和 HMM 剖面技术的多序列比对程序，可扩展到数千条序列，产出高质量的蛋白质和 DNA 比对。
- [**MUSCLE**](https://www.drive5.com/muscle/) - [commandLine] - 通过迭代精炼实现快速准确的多序列比对，在中等规模数据集上提供最佳的速度/准确性平衡。

---

## 3. 序列相似性搜索

在大型数据库中搜索相似序列的工具，包括 BLAST 及其快速替代品。

- [**NCBI BLAST+**](https://blast.ncbi.nlm.nih.gov/) - [commandLine, webtools] - 序列相似性搜索的金标准工具，支持 blastn/blastp/blastx/tblastn/tblastx 多种模式，寻找序列间的局部相似区域并计算统计显著性。
- [**DIAMOND**](https://github.com/bbuchfink/diamond) - [commandLine] - 超快速蛋白比对工具，比 BLASTX 快 100-10,000 倍且保持可比灵敏度，使用双索引和种子扩展优化，广泛用于宏基因组功能注释。
- [**MMseqs2**](https://github.com/soedinglab/MMseqs2) - [commandLine] - 超快速灵敏的蛋白质序列搜索和聚类工具，蛋白序列聚类比 BLAST 快 10,000 倍，搜索比 UBLAST 快 400 倍，包含分类分类和剖面搜索模块。
- [**VSEARCH**](https://github.com/torognes/vsearch) - [commandLine] - 免费开源 64 位宏基因组序列处理工具，提供序列聚类、搜索、嵌合体检测、去冗余和洗牌等功能，是 USEARCH 的广泛使用的开源替代品。
- [**USEARCH**](https://www.drive5.com/usearch/) - [commandLine] - 高性能序列分析工具，提供序列搜索、UPARSE 聚类、OTU 挑选和嵌合体检测，广泛用于 16S rRNA 扩增子分析。
- [**LAST**](https://gitlab.com/mcfrith/last) - [commandLine] - 大规模序列比对工具，专为高效比对测序基因组设计，擅长处理含移码的核苷酸比对。
- [**FASTA**](https://fasta.bioch.virginia.edu/fasta_www2/fasta_list2.shtml) - [commandLine, webtools] - BLAST 出现前的经典序列相似性搜索工具套件，包含蛋白-蛋白（FASTA）、核苷酸-核苷酸（FASTN）和翻译搜索（FASTX/FASTY）。

---

## 4. 基因组组装

支持短读长、长读长或混合策略的 de novo 基因组组装工具。

- [**SPAdes**](https://github.com/ablab/spades) - [commandLine] - De Bruijn 图基因组组装器，适用于细菌、真菌和小型真核基因组，支持 Illumina/IonTorrent/PacBio/ONT，包含 metaSPAdes（宏基因组）、plasmidSPAdes（质粒）、rnaSPAdes（转录组）等模块。
- [**Flye**](https://github.com/fenderglass/Flye) - [commandLine] - 专为 PacBio/ONT 长读长设计的 de novo 组装器，使用 extend-repeat-merge 策略，无需错误校正即可产出高质量组装，支持宏基因组组装。
- [**Canu**](https://github.com/marbl/canu) - [commandLine] - 为 PacBio/ONT 读长设计的长读长组装器，遵循 Celera 组装器方法学三阶段：校正、修剪和组装（重叠-布局-共有），处理长读长典型的高错误率。
- [**Hifiasm**](https://github.com/chhylp123/hifiasm) - [commandLine] - 为 PacBio HiFi reads 设计的单倍型分辨 de novo 组装器，使用基于图的方法产出高度准确的二倍体分相组装。
- [**MEGAHIT**](https://github.com/voutcn/megahit) - [commandLine] - 超快速短读长宏基因组组装器，使用简洁 de Bruijn 图（SdBG）实现低内存消耗和高速运行，适合复杂宏基因组数据集。
- [**ABySS**](https://github.com/bcgsc/abyss) - [commandLine] - 适用于任意大小基因组的短读长 de novo 组装器，使用 de Bruijn 图，通过 MPI 支持并行计算，设计用于双端和 mate-pair 文库。
- [**Velvet**](https://www.ebi.ac.uk/~zerbino/velvet/) - [commandLine] - 使用 de Bruijn 图进行短读长组装的经典 de novo 基因组组装器，是将 de Bruijn 图组装引入真核基因组的开创性工具之一。
- [**SOAPdenovo2**](https://github.com/aquaskyline/SOAPdenovo2) - [commandLine] - BGI 开发的短读长 de novo 组装器，专为在有限内存条件下组装大型基因组设计，被用于大熊猫基因组首次 de novo 组装。
- [**Trinity**](https://github.com/trinityrnaseq/trinityrnaseq) - [commandLine] - RNA-seq 的 de novo 转录组组装器，通过三个模块（Inchworm、Chrysalis、Butterfly）重建全长转录本和异构体。
- [**Unicycler**](https://github.com/rrwick/Unicycler) - [commandLine] - 结合 SPAdes 短读长组装和长读长解决重复区域的细菌基因组混合组装器，输出完整环状染色体和质粒。
- [**MaSuRCA**](https://github.com/alekseyzimin/masurca) - [commandLine] - 从短读长生成 super-reads 后与长读长结合组装的混合组装器，适用于从小到大的各种基因组。
- [**Shasta**](https://github.com/paoloshasta/shasta) - [commandLine] - 专为 ONT 数据设计的快速 de novo 长读长组装器，使用运行长度编码和标记表示法，可在合理计算资源下组装人类规模基因组。
- [**Raven**](https://github.com/lbcb-sci/raven) - [commandLine] - 用于 ONT 和 PacBio CLR 读长的 de novo 长读长组装器，采用重叠-布局-共有方法（无错误校正），设计简洁快速。
- [**IDBA-UD**](https://github.com/loneknightpy/idba) - [commandLine] - 为深度不均的宏基因组短读长数据专门设计的迭代 de Bruijn 图组装器，通过多次迭代不同 k-mer 大小来恢复低丰度物种。

---

## 5. 三代测序分析

PacBio 和 Oxford Nanopore 测序数据的专业分析工具。

- [**minimap2**](https://github.com/lh3/minimap2) - [commandLine] - 多功能序列比对程序，支持 PacBio/ONT 读长比对、长读长重叠检测（高达 ~15% 错误率）、剪接感知比对、Illumina 读长比对、组装间比对和近缘物种全基因组比对。
- [**longshot**](https://github.com/pjedge/longshot) - [commandLine] - 针对长错误率高读长（PacBio/ONT）的二倍体基因组变异调用和分相工具，输入 BAM 输出分相 VCF，目前调用 SNV 并可对 indels 基因分型。
- [**Medaka**](https://github.com/nanoporetech/medaka) - [commandLine] - ONT 数据一致性序列生成工具，使用神经网络校正碱基识别错误，也可进行变异检测。
- [**Guppy**](https://community.nanoporetech.com/) - [commandLine] - ONT 官方碱基识别器，将原始电信号（FAST5）转换为核苷酸序列（FASTQ），支持 GPU 加速，提供多种速度/精度优化模型。
- [**Dorado**](https://github.com/nanoporetech/dorado) - [commandLine] - ONT 新一代碱基识别器，取代 Guppy，具有改进的准确度、POD5 格式支持和 simplex/duplex/碱基修饰检测模型。
- [**pbmm2**](https://github.com/PacificBiosciences/pbmm2) - [commandLine] - PacBio 官方 minimap2 封装，用于 SMRT 读长比对，提供一致的输出格式并针对 PacBio 数据特征优化。
- [**CCS**](https://github.com/PacificBiosciences/ccs) - [commandLine] - 循环一致性测序工具，从 PacBio subreads 计算多轮聚合酶循环的共有序列，生成精度 ≥Q20（通常 Q30）的 CCS/HiFi reads。
- [**IsoSeq**](https://github.com/PacificBiosciences/IsoSeq) - [commandLine] - PacBio 全长转录组分析流程，通过聚类和校正处理 Iso-Seq 数据生成全长转录本序列。
- [**Clair3**](https://github.com/HKU-BAL/Clair3) - [commandLine] - 基于深度学习的高精度 ONT 变异检测工具，使用 pileup 和全比对两阶段策略，准确度可媲美短读长流程。
- [**NanoPlot**](https://github.com/wdecoster/NanoPlot) - [commandLine] - 长读长测序数据可视化工具，生成读长分布、质量分数、产量和时间序列等图表。
- [**Pychopper**](https://github.com/nanoporetech/pychopper) - [commandLine] - 识别、定向和修剪全长 Nanopore cDNA reads，将 reads 分类为全长、截短或嵌合型，是转录组学预处理的重要工具。

---

## 6. RNA-seq 分析

转录组分析工具，包括比对、定量、差异表达和可变剪接检测。

- [**STAR**](https://github.com/alexdobin/STAR) - [commandLine] - 超快速 RNA-seq 比对器，使用后缀数组中的顺序最大可映射种子搜索，可直接检测剪接位点、嵌合 RNA 和环状 RNA。
- [**Salmon**](https://github.com/COMBINE-lab/salmon) - [commandLine] - 使用准映射技术绕过完整比对的转录本表达定量工具，自动校正 GC 偏好、片段长度偏好和序列特异性偏好。
- [**Kallisto**](https://github.com/pachterlab/kallisto) - [commandLine] - 基于 k-mer 匹配的伪比对技术快速定量 RNA-seq 转录本丰度，在标准笔记本上几分钟内处理 3000 万条 reads，支持 Bootstrap 不确定性量化。
- [**RSEM**](https://github.com/deweylab/RSEM) - [commandLine] - 使用期望最大化算法定量 RNA-seq 基因和转录本丰度，支持双端、链特异性和变长 reads，计算后验概率。
- [**featureCounts**](https://subread.sourceforge.net/) - [commandLine] - 高效的 reads 定量程序，统计比对到基因组特征的 reads 数量，以速度和低内存使用著称。
- [**StringTie**](https://ccb.jhu.edu/software/stringtie/) - [commandLine] - 从 RNA-seq 比对结果组装转录本并估算丰度，可发现新异构体，支持参考注释集成。
- [**DESeq2**](https://bioconductor.org/packages/DESeq2/) - [R] - 基于负二项分布 GLM 的差异表达分析 R/Bioconductor 包，使用收缩估计提高稳定性和可解释性，支持复杂实验设计。
- [**edgeR**](https://bioconductor.org/packages/edgeR/) - [R] - 使用经验贝叶斯估计和基于负二项模型精确检验的差异表达分析 R/Bioconductor 包，特别适合小样本量。
- [**limma**](https://bioconductor.org/packages/limma/) - [R] - 使用线性模型的差异表达分析 R 包，通过 voom 变换扩展到 RNA-seq，支持复杂设计。
- [**Cufflinks**](https://github.com/cole-trapnell-lab/cufflinks) - [commandLine] - 开创性 RNA-seq 分析套件，包含转录本组装（Cufflinks）、比较（Cuffcompare）和差异表达（Cuffdiff）。
- [**Arriba**](https://github.com/suhrig/arriba) - [commandLine] - 从 RNA-seq 检测基因融合的工具，使用 STAR 比对输出，包含全面过滤去除假象并提供验证可视化。
- [**STAR-Fusion**](https://github.com/STAR-Fusion/STAR-Fusion) - [commandLine] - 从 RNA-seq 检测基因融合的全面流程，广泛用于癌症基因组学。
- [**rMATS**](https://rnaseq-mats.sourceforge.net/) - [commandLine, webtools] - 从 RNA-seq 检测差异可变剪接事件（外显子跳跃、内含子保留、可变 5'/3' 剪接位点、互斥外显子）的统计模型工具。
- [**SUPPA2**](https://github.com/comprna/SUPPA2) - [commandLine] - 快速可变剪接分析工具，计算 PSI 值并进行差异剪接分析，可与多种定量工具集成。

---

## 7. 变异检测与注释

从测序数据中检测 SNV、indel 并注释功能效应的工具。

- [**samtools**](https://github.com/samtools/samtools) - [commandLine] - SAM/BAM/CRAM 文件操作套件，提供排序、索引、合并、过滤和 mpileup 功能生成基因型似然值。
- [**bcftools**](https://github.com/samtools/bcftools) - [commandLine] - VCF/BCF 文件变异调用和操作工具集，包含 SNP/indel 调用、过滤、比较和格式转换。
- [**GATK**](https://gatk.broadinstitute.org/) - [commandLine] - 行业标准变异检测工具包，包含 HaplotypeCaller（胚系）、Mutect2（体细胞）、联合基因分型、VQSR 过滤和 CNV 分析。
- [**DeepVariant**](https://github.com/google/deepvariant) - [commandLine] - 基于深度学习的变异检测工具，使用卷积神经网络将 pileup 图像分类识别变异，无需人工质量过滤。
- [**freebayes**](https://github.com/freebayes/freebayes) - [commandLine] - 基于贝叶斯单倍型的变异检测工具，检测 SNP、indel、MNP 和复杂事件，适合高倍体或混合群体样本。
- [**Strelka2**](https://github.com/Illumina/strelka) - [commandLine] - Illumina 开发的快速准确胚系/体细胞变异检测工具，使用单倍型方法和经验贝叶斯。
- [**Mutect2**](https://gatk.broadinstitute.org/) - [commandLine] - GATK 中的体细胞变异检测工具，专为肿瘤-正常配对设计，使用概率模型考虑纯度、错误和污染。
- [**VarScan2**](https://github.com/dkoboldt/varscan) - [commandLine] - 从肿瘤-正常配对测序数据检测体细胞突变、CNV 和 LOH。
- [**VEP**](https://www.ensembl.org/info/docs/tools/vep/) - [commandLine, webtools] - Ensembl 变异效应预测器，使用 Ensembl 注释确定变异对基因、转录本和蛋白质的功能影响。
- [**ANNOVAR**](https://annovar.openbioinformatics.org/) - [commandLine] - 遗传变异注释工具，支持基于基因、区域和过滤的注释，兼容 RefSeq、UCSC 和 ENCODE 数据库。
- [**SnpEff**](https://pcingola.github.io/SnpEff/) - [commandLine] - 快速变异注释工具，将效应分类为高/中/低/修饰影响，与 Galaxy 紧密集成。
- [**Pindel**](https://github.com/genome/pindel) - [commandLine] - 使用模式增长方法从双端短读长检测大缺失和中等插入的断点，对转座元件插入特别有效。

---

## 8. CRISPR 靶向深度测序分析

从深度测序和 Sanger 测序数据分析 CRISPR-Cas9 基因组编辑结果的工具。

- [**CRISPResso2**](https://github.com/pinellolab/CRISPResso2) - [commandLine, webtools] - 综合分析基因组编辑实验深度测序数据的套件，比对 reads 到参考、定量插入/缺失/突变、判断每条 reads 是否被编辑，以直观图表汇总编辑结果。
- [**CrispRVariants**](https://bioconductor.org/packages/CrispRVariants/) - [R] - R/Bioconductor 包，分析 CRISPR-Cas9 诱变实验结果，定位相对于切割位点的变异等位基因组合并计算突变率。
- [**CRISP-ID**](http://crispid.gbiomed.kuleuven.be) - [webtools] - 在线工具，基于 Sanger 测序检测 CRISPR 靶向区域精确 indel 大小和位置，识别单等位和双等位编辑事件。
- [**Microhomology-Predictor**](http://www.rgenome.net/mich-calculator/) - [webtools] - 在线工具，计算 ZFN/TALEN/Cas9 的微同源相关分数，预测 MMEJ 修复结果。
- [**TIDE**](https://tide.nki.nl/) - [webtools] - 通过分解 Sanger 色谱图定量 CRISPR 编辑效率，识别 indent 模式频率，无需深度测序。
- [**ICE**](https://ice.synthego.com/) - [webtools] - Synthego 在线工具，从 Sanger 数据分析 CRISPR 编辑结果，提供编辑效率、常见结果和敲除评分。
- [**CHOPCHOP**](https://chopchop.cbu.uib.no/) - [webtools] - 设计 CRISPR/Cas9/Cas12a gRNA，预测靶向效率和脱靶特异性，提供基因组可视化。
- [**CRISPick**](https://portals.broadinstitute.org/gppx/crispick) - [webtools] - Broad Institute gRNA 设计工具，结合效率和特异性评分，与全基因组 CRISPR 敲除文库集成。

---

## 9. GWAS 全基因组关联分析

全基因组关联分析、精细定位和多基因风险评分的统计工具。

- [**PLINK 1.9**](https://www.cog-genomics.org/plink/) - [commandLine] - 免费开源全基因组关联分析工具集，提供数据管理、质控、群体分层检测、关联检验和 LD 计算。
- [**PLINK 2.0**](https://www.cog-genomics.org/plink/2.0/) - [commandLine] - 下一代 PLINK，引入多线程、PGEN 格式支持，性能改进以支持生物样本库规模数据。
- [**BOLT-LMM**](https://alkesgroup.broadinstitute.org/BOLT-LMM/) - [commandLine] - 复杂性状 GWAS 的贝叶斯混合模型关联方法，使用效应大小先验提高功效，同时控制群体结构。
- [**SAIGE**](https://github.com/weizhouUMICH/SAIGE) - [commandLine] - 大规模队列的可扩展混合模型方法，使用鞍点逼近处理病例-对照不平衡，即使极端不平衡表型也能提供准确 p 值。
- [**GEMMA**](https://github.com/genetics-statistics/GEMMA) - [commandLine] - 使用标准线性混合模型进行全基因组关联分析，包括方差组分估计、关联检验和 BSLMM 分析。
- [**GCTA**](https://yanglab.westlake.edu.cn/software/gcta/) - [commandLine] - 估计 SNP 遗传力、遗传相关，进行条件/联合 GWAS 分析和多性状分析。
- [**MAGMA**](https://cncr.nl/research/magma/) - [commandLine, webtools] - 将 SNP 关联聚合到基因水平的基因-基关联分析，考虑 LD 并整合功能注释。
- [**FINEMAP**](http://www.christianbenner.com/) - [commandLine] - 从 GWAS 汇总统计识别因果变异的贝叶斯精细定位工具，计算后验概率。
- [**SuSiE**](https://github.com/stephenslab/susieR) - [R] - R 精细定位包，使用"单效应和"回归模型识别因果变异的置信集。
- [**FUMA**](https://fuma.ctglab.nl/) - [webtools] - GWAS 结果功能注释 Web 平台，识别独立显著 SNP、映射到基因并进行富集分析。
- [**Michigan Imputation Server**](https://imputationserver.sph.umich.edu/) - [webtools] - 免费云端基因型填充服务，使用 Minimac4 和 HRC/TOPMed/1000 Genomes 参考面板。
- [**Minimac4**](https://github.com/statgen/Minimac4) - [commandLine] - 快速基因型填充工具，使用隐马尔可夫模型，针对大规模填充高度优化。
- [**Beagle 5**](https://faculty.washington.edu/browning/beagle/beagle.html) - [commandLine] - Java 基因分相和填充程序，使用单倍型簇的图形模型，准确性高。
- [**PRSice-2**](https://github.com/choishingwan/PRSice) - [commandLine] - 多基因风险评分计算，使用 GWAS 效应大小加权，自动优化 p 值阈值并提供绘图。
- [**LDpred2**](https://github.com/privefl/bigsnpr) - [R] - R 贝叶斯 PRS 构建包，考虑 SNP 间 LD，提供多种模型提高预测准确性。

---

## 10. 群体遗传学

群体遗传学分析、模拟和推断工具。

- [**msprime**](https://github.com/tskit-dev/msprime) - [commandLine] - 可扩展溯祖模拟器，使用树序列数据结构高效模拟重组下的全基因组变异。
- [**SLiM**](https://github.com/MesserLab/SLiM) - [commandLine] - 正向进化模拟器，建模选择、突变、重组和群体结构变化，详细指定遗传架构。
- [**Fastsimcoal2**](http://cmpg.unibe.ch/software/fastsimcoal2/) - [commandLine] - 模拟复杂人口统计场景下的遗传多样性，通过最大似然从 SFS 推断参数。
- [**Dadi**](https://github.com/rajanil/dadi) - [commandLine] - 使用等位基因频谱扩散近似推断群体历史，含分化、迁移和规模变化。
- [**PSMC**](https://github.com/lh3/psmc) - [commandLine] - 从单个二倍体基因组推断历史有效群体规模，使用隐马尔可夫模型，广泛用于古群体动态研究。
- [**MSMC2**](https://github.com/stschiff/msmc2) - [commandLine] - 将 PSMC 扩展到多个分相基因组，提供更高分辨率的群体规模和分化时间推断。
- [**SMC++**](https://github.com/terhorst/smcpp) - [commandLine] - 使用复合似然方法联合建模 SFS 和 LD 信息估计群体规模历史。
- [**ADMIXTURE**](https://dalexander.github.io/admixture/) - [commandLine] - 最大似然法个体祖先成分估计，计算效率高，可处理大规模数据。
- [**STRUCTURE**](https://web.stanford.edu/group/pritchardlab/structure.html) - [commandLine] - 奠基性贝叶斯群体结构聚类方法，从基因型数据推断群体分配和混合。
- [**EIGENSOFT**](https://github.com/DReichLab/EIG) - [commandLine] - 包含 smartpca（PCA 群体结构检测）、Tracy-Widom 检验和 EIGENSTRAT 分层校正。
- [**TREEMIX**](https://github.com/joepickrell/treemix) - [commandLine] - 从等位基因频率推断群体分裂和混合事件，构建包含迁移边的树。
- [**AdmixTools**](https://github.com/DReichLab/AdmixTools) - [commandLine] - 使用 f-统计量分析群体混合，包含 D-统计量、qpWave/qpAdm 定量祖先比例。
- [**VCFtools**](https://vcftools.github.io/) - [commandLine] - VCF 过滤和群体遗传学计算工具，包括 Fst、Tajima D、核苷酸多样性和 LD 估计。
- [**ANGSD**](https://github.com/angsd/angsd) - [commandLine] - 直接使用基因型似然值分析 NGS 群体遗传学数据，估计等位基因频率、Fst 和 IBS 矩阵。
- [**SweepFinder2**](https://github.com/3hunter/SweepFinder2) - [commandLine] - 使用复合似然比检验检测选择性清除，无需知道选择等位基因。
- [**Selscan**](https://github.com/szpiech/selscan) - [commandLine] - 计算 iHS、XP-EHH 和 XP-CLR 统计量检测近期正选择。
- [**rehh**](https://cran.r-project.org/web/packages/rehh/) - [R] - R 包，通过 EHH 统计量检测选择，提供 iHS 和 Rsb 计算及可视化。
- [**pixy**](https://github.com/ksamuk/pixy) - [commandLine] - 从 VCF 文件计算核苷酸多样性（π）和分化（dxy），使用基因型似然法处理缺失数据。

---

## 11. 结构变异检测

专门检测和分析结构变异的工具。

- [**Snowman**](https://github.com/broadinstitute/SnowmanSV) - [commandLine] - 使用全基因组局部组装检测结构变异，高精度识别断点。
- [**seeksv**](https://github.com/qiukunlong/seeksv) - [commandLine] - 精确的结构变异和病毒整合检测，利用不一致 read 对和 split-read。
- [**Genomon SV**](https://github.com/Genomon-Project/GenomonSV) - [commandLine] - 从癌症基因组检测体细胞结构变异（缺失、倒位、易位）。
- [**BreakDancer**](https://github.com/genome/breakdancer) - [commandLine] - 全基因组 SV 检测，使用异常分离距离或方向的双端 reads 预测五种 SV 类型。
- [**Manta**](https://github.com/Illumina/manta) - [commandLine] - 从双端测序检测 SV 和中等 indels，针对大规模队列优化。
- [**GRIDSS**](https://github.com/PapenfussLab/gridss) - [commandLine] - 使用局部组装和机器学习检测 SV 断点，区分真实 SV 和假象。
- [**LUMPY**](https://github.com/arq5x/lumpy-sv) - [commandLine] - 集成不一致 read 对、split-read、读深度和先验知识的概率框架 SV 检测。
- [**Wham**](https://github.com/zeeev/wham) - [commandLine] - 对 read 映射特征聚类识别 SV，无需预定义 SV 模型。
- [**SVision**](https://github.com/xjtu-omics/SVision) - [commandLine] - 基于深度学习的 SV 检测和可视化，从长读长比对模式直接检测 SV。
- [**CuteSV**](https://github.com/tjiangHIT/cuteSV) - [commandLine] - 快速灵敏的长读长 SV 检测，使用聚类和集成学习，支持 PacBio 和 ONT。

---

## 12. 宏基因组学

微生物群落分类和功能分析工具。

- [**MetaPhlAn 4**](https://huttenhower.sph.harvard.edu/metaphlan/) - [commandLine] - 使用进化枝特异性标记基因进行菌株水平微生物谱分析。
- [**Kraken 2**](https://ccb.jhu.edu/software/kraken2/) - [commandLine] - 压缩数据库上的快速 k-mer 比对分类器。
- [**Bracken**](https://ccb.jhu.edu/software/bracken/) - [commandLine] - 使用贝叶斯概率从 Kraken 结果估计物种丰度。
- [**HUMAnN 3**](https://huttenhower.sph.harvard.edu/humann/) - [commandLine] - 估计基因家族和代谢通路丰度的宏基因组功能谱分析。
- [**MEGAN6**](https://ab.inf.uni-tuebingen.de/software/megan6/) - [commandLine] - 使用 LCA 算法的交互式分类和功能分析。
- [**Centrifuge**](https://ccb.jhu.edu/software/centrifuge/) - [commandLine] - 基于 FM-index 的快速分类器，低内存占用。
- [**MetaBAT 2**](https://bitbucket.org/berkeleylab/metabat) - [commandLine] - 使用四核苷酸频率和丰度的宏基因组 binning。
- [**CONCOCT**](https://github.com/BinPro/CONCOCT) - [commandLine] - 使用高斯混合模型的无监督宏基因组 binning。
- [**MaxBin 2**](https://sourceforge.net/projects/maxbin2/) - [commandLine] - 基于 EM 算法的自动化宏基因组 binning。
- [**CheckM**](https://github.com/Ecogenomics/CheckM) - [commandLine] - 使用谱系标记基因评估 MAG 质量。
- [**GTDB-Tk**](https://ecogenomics.github.io/GTDBTk/) - [commandLine] - 使用 GTDB 数据库进行细菌和古菌分类。
- [**Mash**](https://mash.readthedocs.io/) - [commandLine] - 使用 MinHash 进行快速基因组距离估计。
- [**Pavian**](https://github.com/fbreitwieser/pavian) - [webtools] - 交互式宏基因组分类结果 Web 可视化。
- [**QIIME 2**](https://qiime2.org/) - [commandLine] - 下一代微生物组生物信息学平台，支持扩增子和鸟枪法宏基因组学。

---

## 13. 表观基因组学

DNA 甲基化、组蛋白修饰和染色质可及性分析工具。

- [**MACS2/3**](https://github.com/macs3-project/MACS) - [commandLine] - 基于泊松模型的 ChIP-seq/ATAC-seq/CUT&Tag 峰检出工具。
- [**Bismark**](https://www.bioinformatics.babraham.ac.uk/projects/bismark/) - [commandLine] - 亚硫酸氢盐测序比对和单碱基分辨率 DNA 甲基化调用。
- [**MethylKit**](https://github.com/al2na/methylKit) - [R] - R/Bioconductor DNA 甲基化差异分析和注释包。
- [**MethylDackel**](https://github.com/dpryan79/MethylDackel) - [commandLine] - 快速甲基化提取器，从 Bismark 比对结果提取 CpG 甲基化指标。
- [**deepTools**](https://deeptools.readthedocs.io/) - [commandLine] - 表观基因组数据归一化、相关分析和可视化 Python 套件。
- [**HOMER**](https://homer.ucsd.edu/homer/) - [commandLine] - 峰检出、基序发现和注释套件。
- [**ChIPseeker**](https://github.com/YuLab-SMU/ChIPseeker) - [R] - R 峰注释和可视化包。
- [**bwa-meth**](https://github.com/brentp/bwa-meth) - [commandLine] - 专为全基因组亚硫酸氢盐测序设计的比对工具。
- [**GemBS**](https://github.com/heathsc/gemBS) - [commandLine] - 全面亚硫酸氢盐测序分析平台，集成质控报告。
- [**Picard Tools**](https://broadinstitute.github.io/picard/) - [commandLine] - Java 高通量测序数据处理工具集。

---

## 14. 蛋白质组学

基于质谱的蛋白质组学软件工具。

- [**MaxQuant**](https://www.maxquant.org/) - [commandLine] - 定量蛋白质组学软件，支持 LFQ/SILAC/TMT/iTRAQ 和 Andromeda 搜索引擎。
- [**FragPipe**](https://fragpipe.nesvilab.org/) - [commandLine] - 集成 MSFragger 鉴定、定量和验证的完整蛋白质组学流程。
- [**MSFragger**](https://msfragger.nesvilab.org/) - [commandLine] - 使用片段离子索引实现超快速肽段鉴定，支持开放质量搜索鉴定 PTM。
- [**OpenMS**](https://www.openms.de/) - [commandLine] - 开源 LC-MS 数据分析框架。
- [**Skyline**](https://skyline.ms/) - [commandLine] - Windows 靶向蛋白质组学方法创建和 SRM/PRM/DIA 数据分析。
- [**DIA-NN**](https://github.com/vdemichev/DiaNN) - [commandLine] - 基于深度神经网络的 DIA 蛋白质组学数据分析。
- [**ProteoWizard**](https://proteowizard.sourceforge.io/) - [commandLine] - 质谱数据格式转换工具。
- [**Comet**](https://uwpr.github.io/Comet/) - [commandLine] - 开源串联质谱搜索引擎。
- [**X!Tandem**](https://www.thegpm.org/tandem/) - [commandLine] - 流行开源串联质谱搜索引擎。
- [**Sage**](https://github.com/lazear/sage) - [commandLine] - 开源 GPU 加速肽搜索引擎，Rust 编写。
- [**pFind**](http://pfind.ict.ac.cn/software/pFind/) - [commandLine] - 高性能蛋白质组学数据库搜索工具。

---

## 15. 单细胞组学分析

单细胞 RNA-seq、ATAC-seq 和空间转录组学分析工具。

- [**Seurat**](https://satijalab.org/seurat/) - [R] - R 单细胞 RNA-seq 综合分析包，支持质控、聚类、差异表达、多样本整合和空间转录组学。
- [**Scanpy**](https://scanpy.readthedocs.io/) - [commandLine] - 可扩展 Python 单细胞基因表达分析工具包，专为大规模数据集设计。
- [**Cell Ranger**](https://support.10xgenomics.com/) - [commandLine] - 10x Genomics 官方单细胞数据处理流程。
- [**Monocle 3**](https://cole-trapnell-lab.github.io/monocle3/) - [R] - R 单细胞轨迹推断包，使用反向图嵌入学习发育轨迹。
- [**scvi-tools**](https://scvi-tools.org/) - [commandLine] - 单细胞组学深度生成建模框架，包含 scVI/scANVI/TotalVI。
- [**Harmony**](https://portals.broadinstitute.org/harmony/) - [commandLine] - 快速单细胞数据整合，通过 PCA 空间对齐校正批次效应。
- [**scVelo**](https://scvelo.readthedocs.io/) - [commandLine] - RNA 速度分析框架，使用剪接动力学动态模型。
- [**CellChat**](https://github.com/sqjin/CellChat) - [R] - R 细胞间通信网络推断包，使用配体-受体相互作用数据库。
- [**Squidpy**](https://squidpy.readthedocs.io/) - [commandLine] - Python 空间单细胞分析框架，整合转录组学和空间成像数据。

---

## 16. 系统发育学

从分子序列构建和分析系统发育树的工具。

- [**IQ-TREE 2**](http://www.iqtree.org/) - [commandLine] - 最大似然法系统发育推断，支持超快 Bootstrap、ModelFinder 自动模型选择和分区检验。
- [**RAxML-NG**](https://github.com/amkozlov/raxml-ng) - [commandLine] - 下一代 RAxML 最大似然法系统发育推断，针对大数据集优化。
- [**FastTree 2**](http://www.microbesonline.org/fasttree/) - [commandLine] - 使用启发式方法从数十万条序列的大规模比对快速构建近似最大似然树。
- [**MrBayes**](https://nbisweden.github.io/MrBayes/) - [commandLine] - 使用 MCMC 的贝叶斯系统发育推断，支持混合模型和后验概率估计。
- [**BEAST 2**](https://www.beast2.org/) - [commandLine] - 时间校准系统发育分析，用于分歧时间估计和群体动态。
- [**MEGA 11**](https://megasoftware.net/) - [commandLine] - 集成桌面软件，支持序列比对、系统发育树构建和分子进化分析。
- [**APE**](https://cran.r-project.org/package=ape) - [R] - R 基础系统发育学计算包。
- [**ggtree**](https://bioconductor.org/packages/ggtree/) - [R] - R 系统发育树可视化包，使用 ggplot2 语法。
- [**TrimAl**](https://github.com/scapella/trimal) - [commandLine] - 自动去除多序列比对中比对差的区域。

---

## 17. 蛋白结构与分子模拟

蛋白质结构预测、分子动力学模拟和计算药物发现工具。

- [**AlphaFold2**](https://github.com/google-deepmind/alphafold) - [commandLine] - DeepMind 革命性深度学习模型，从序列预测蛋白质 3D 结构达到接近实验精度。
- [**AlphaFold3**](https://github.com/google-deepmind/alphafold3) - [commandLine] - 扩展到蛋白-配体/DNA/RNA 复合物结构预测，使用基于扩散的架构。
- [**ColabFold**](https://github.com/sokrypton/ColabFold) - [webtools] - 通过 Google Colab 快速运行 AlphaFold2，无需大量计算资源。
- [**ESMFold**](https://github.com/facebookresearch/esm) - [commandLine] - Meta 的快速蛋白质结构预测模型，比 AlphaFold2 快约 60 倍。
- [**GROMACS**](https://www.gromacs.org/) - [commandLine] - 高性能分子动力学引擎，支持 GPU 加速和并行计算。
- [**AMBER**](https://ambermd.org/) - [commandLine] - 生物分子分子动力学模拟套件。
- [**OpenMM**](https://openmm.org/) - [commandLine] - GPU 加速分子模拟工具包，提供 Python API。
- [**AutoDock Vina**](https://vina.scripps.edu/) - [commandLine] - 开源分子对接和虚拟筛选工具。
- [**PyMOL**](https://pymol.org/) - [commandLine] - 分子可视化系统，支持高质量 3D 图像渲染。
- [**Rosetta**](https://www.rosettacommons.org/) - [commandLine] - 蛋白质结构预测、设计和蛋白-蛋白对接综合套件。

---

## 18. Hi-C 与染色质构象分析

从 Hi-C 数据分析 3D 染色质结构的工具。

- [**Juicer**](https://github.com/aidenlab/juicer) - [commandLine] - 全面 Hi-C 数据处理流程，从比对到接触矩阵生成。
- [**Juicebox**](https://github.com/aidenlab/Juicebox) - [commandLine] - Hi-C 接触图交互式二维热图可视化工具。
- [**HiC-Pro**](https://github.com/nservant/HiC-Pro) - [commandLine] - 优化的 Hi-C 数据处理流程，包含比对、过滤和质控。
- [**Cooler**](https://github.com/open2c/cooler) - [commandLine] - 稀疏压缩 Hi-C 接触矩阵存储格式和工具。
- [**HiGlass**](https://github.com/higlass/higlass) - [webtools] - 基于 Web 的多分辨率基因组数据查看器。
- [**CALDER**](https://github.com/CSOgroup/CALDER2) - [commandLine] - 使用基于图的方法计算染色质区室域和 TAD 样结构域。
- [**FAN-C**](https://github.com/vaquerizaslab/fanc) - [commandLine] - Hi-C 分析和可视化命令行工具。

---

## 19. 通用工具与工作流

跨领域生物信息学工具和工作流管理系统。

- [**Galaxy**](https://galaxyproject.org/) - [webtools] - 开放 Web 数据分析平台，支持工具运行、工作流创建和分析共享。
- [**Bioconda**](https://bioconda.github.io/) - [commandLine] - 基于 Conda 的生物信息学软件分发渠道，跨平台依赖解析。
- [**Bioconductor**](https://www.bioconductor.org/) - [R] - 基于 R 的策划生物信息学软件仓库。
- [**Biopython**](https://biopython.org/) - [commandLine] - Python 生物计算工具集。
- [**Nextflow**](https://www.nextflow.io/) - [commandLine] - 可扩展可复现工作流管理系统，支持容器和云/HPC。
- [**Snakemake**](https://snakemake.readthedocs.io/) - [commandLine] - 基于 Python 的工作流管理器。
- [**nf-core**](https://nf-co.re/) - [commandLine] - 社区策划的最佳实践 Nextflow 流程集合。
- [**Apptainer/Singularity**](https://apptainer.org/) - [commandLine] - 面向 HPC 的可复现容器平台。
- [**IGV**](https://software.broadinstitute.org/software/igv/) - [commandLine] - 高性能基因组数据可视化工具。
- [**UCSC Genome Browser**](https://genome.ucsc.edu/) - [webtools] - 基于 Web 的基因组浏览器。
- [**bedtools**](https://github.com/arq5x/bedtools2) - [commandLine] - 基因组区间运算工具套件。
- [**htslib/tabix**](https://github.com/samtools/htslib) - [commandLine] - HTS 数据 C 库和索引工具。
- [**pysam**](https://github.com/pysam-developers/pysam) - [commandLine] - htslib 的 Python 封装。
- [**SeqKit**](https://github.com/shenwei356/seqkit) - [commandLine] - 跨平台超快速 FASTA/Q 文件操作工具包。

---

## License

MIT License

*更新于: 2026年6月*
