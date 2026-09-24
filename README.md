[![testing](https://github.com/ctglab/excavator2/actions/workflows/test.yml/badge.svg)](https://github.com/ctglab/excavator2/actions/workflows/test.yml)
![GitHub docker image build](https://github.com/ctglab/excavator2/actions/workflows/docker.yml/badge.svg)
[![DOI](https://img.shields.io/badge/DOI-10.1093/nar/gkw718-blue.svg)](https://doi.org/10.1093/nar/gkw695)
# EXCAVATOR2

The Python-first port is being developed in `src/excavator2/`, with a small C++
extension under `cpp/`. See the [development setup](docs/development.md) and
[porting roadmap](docs/python-cpp-porting-roadmap.md). The development package supports
[FastCall and HSLM analysis from converted legacy inputs](docs/hslm-analysis.md).
Use the original workflow below for target generation and read preparation.

> The first **read count based** tool that exploits **all the reads** produced by **WES** experiments to detect **CNVs** with a **genome-wide resolution**.

**Copy Number Variants** (CNVs) are structural rearrangements contributing to phenotypic variation that have been proved to be **associated with many disease states**. Over the last years, the identification of CNVs from **whole-exome sequencing** (WES) has become a common practice for both research and clinical purpose as it represents a cost-effective alternative to whole-genome for the study of disease-associated variants affecting just the coding regions. Moreover, the sequencing of a smaller range of genomic regions—i.e., **Targeted Sequencing** (TS)—that are directly associated with disease and/or have a straight functional interpretation is widespread and extensively used for diagnostic purposes and treatment-response monitoring, especially in oncology.  

**EXCAVATOR2** is an Read Count-based tool that extends the functionalities of EXCAVATOR (Magi et al., 2013), its predecessor, to the identification of CNVs (overlapping or not overlapping exons) by integrating the analysis of **In-target** and **Off-target** reads produced by a WES experiment. Almost all commercial enrichment kits for WES (Samuels et al., 2013) indeed suffer from lack of specificity for the target regions and produce a significant amount of reads (up to 40%) mapping outside the target design. EXCAVATOR2 has two significant advantages: it makes it possible to investigate, from a WES experiment, the impact of CNVs also on the non-coding genome, and to recover with more precision most breakpoints, which usually fall outside of the targeted exons since they only encompass a sparse 1% of the genome.
Furthermore, EXCAVATOR2 can be exploited to investigate CNVs from TS experiments in which enriched regions correspond to few Mbs, usually a panel of few genes. Also in this case, a considerable portion of the sequenced reads are off-target reads originating from the rest of the genome. They provide a very low-coverage signal of the whole genome and can be an indicator of regions with copy-number alterations on a larger scale.
EXCAVATOR2 is tailored for the detection of both germline and somatic CNVs from different sequencing experiments (WES and TS) in various disease contexts and population genetic studies as well as in **small** and **large-scale re-sequencing population and cancer studies**.  

EXCAVATOR2 is a collection of bash, R, Perl, and Fortran scripts, it combines a three-step normalization procedure with a novel heterogeneous hidden Markov model algorithm and a calling method that classifies genomic regions into five copy number states.

## Article
Please **cite this publication** when using **EXCAVATOR2** software:

**[Enhanced copy number variants detection from whole-exome sequencing data using EXCAVATOR2](https://academic.oup.com/nar/article/44/20/e154/2607979/Enhanced-copy-number-variants-detection-from-whole)**

D'Aurizio R, Pippucci T, Tattini L, Giusti B, Pellegrini M, Magi A  
Nucleic Acids Res. (2016) 44 (20):e154 

## Installation

**EXCAVATOR2** could be either installed in a dedicated conda environment using the dedicated environment or using Docker. To start clone the repo with 

```
git clone https://github.com/ctglab/excavator2.git
```

If you wish to create a conda environment, use the file `excavator2.yml` as follows:

```
conda create -f excavator2.yml
```

Alternatively, if you want to use Docker, you could create the container from scratch using the `Dockerfile` within the repository or by pulling from the Dockerhub registry as follows

```
docker pull ctglab/excavator2:latest
```

Or Singularity/Apptainer (here showing Singularity, but the command line is the same)

```
singularity pull excavator2.sif docker://ctglab/excavator2:latest
```

## Testing the installation

We provided a minimal example of a scenario with Tumor/Normal WES samples downloaded from SRA, together with support files in the folder `.test`. 

Due to the size restriction in Github, two additional files are required (the human genome and the mappability track) that could be downloaded as follows:

```sh
curl -L https://hgdownload.soe.ucsc.edu/gbdb/hg38/hoffmanMappability/k100.Bismap.MultiTrackMappability.bw --output .test/ref/k100.Bismap.MultiTrackMappability.bw 
curl -L ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/references/GRCh38/GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta.gz --output .test/ref/GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta.gz
gzip -d .test/ref/GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta.gz
```
Then, in the conda environment scenario

```sh
conda activate excavator2
export PATH="/opt/excavator2:${PATH}"
TargetPerla.pl \
  -v -f -s .test/config.yaml \
  -o .test/output 
EXCAVATORDataPrepare.pl \
  -v -f \
  -s .test/sample_sheet.yaml \
  -t .test/output/hg38/SureSelectV7/w_30000 \
  -o .test/outputPrepare \
  -@ 3
EXCAVATORDataAnalysis.pl \
  -v -f \
  -s .test/sample_file_list.yaml \
  -i .test/outputPrepare \
  -t .test/output/hg38/SureSelectV7/w_30000 \
  -o .test/outputAnalysis \
  -@ 3 -e paired
```
The whole testing procedure takes around 4-5 minutes on a modern architecture CPU.

## Usage
**EXCAVATOR2** runs in three steps:
- [TargetPerla](#1-targetperla)
- [EXCAVATORDataPrepare](#2-excavatordataprepare)
- [EXCAVATORDataAnalysis](#3-excavatordataanalysis)

**IMPORTANT**: All file paths specified in CLI should be absolute paths. Using relative paths could lead to unexpected results.

### 1. TargetPerla
This module takes as **input** the **reference genome** and a **BED file** with the target regions of WES or genes panel to provide.

#### Command Line Interface

```sh
TargetPerla.pl \
  -v -f \
  -s <settings file> \
  -o <output folder> \
2> <log file>
```
Using singularity:
```
singularity run \
  --bind /home/workspace/src/exca2data/:srv/ \
  excavator2.sif TargetPerla.pl \
  -v -f \
  -s config.yaml \
  -o output \
2> logTarget.txt
```

|**Option**        |**Description**                                  |
|------------------|-------------------------------------------------|
|`-h, --help`      |Print help message.                              |
|`-m, --man`       |Print complete documentation.                    |
|`-v, --verbose`   |Use verbose output.                              |
|`-f, --force`     |Force output overwrite.                          |
|`-s, --settings`  |Path to YAML [settings file](#settings-file).    |
|`-o, --output`    |Path to output folder.                           |

**Log** messages can be saved by redirecting **STDERR** to a file.

#### Settings file
It's a [YAML](https://en.wikipedia.org/wiki/YAML) file containing the paths to **relevant files** and the **assembly name** (e.g. _hg19_ or _hg38_).  
```yaml
# Reference genome
Reference:
  Assembly: hg38 # assembly name
  FASTA: /srv/data/ucsc.hg38.fasta # path to reference genome FASTA file
  BigWig: /srv/data/ucsc.hg38.bw # path to reference genome BigWig file
  Chromosomes: /srv/data/ChromosomeCoordinate_hg38.txt # path to reference genome chromosome coordinates file
  Centromeres: /srv/data/CentromerePosition_hg38.txt # path to reference genome centromere positions file
  Gaps: /srv/data/GapHg38.UCSC.txt # path to reference genome gap positions file

# Target regions
Target:
  Name: 1000GP # target name
  BED: test/1000GPTarge.bed # path to target BED file
  Window: 30000 # sliding window size (bp)
```
The bigWig files are binary files reporting information about mappability for hg19 and hg38 assemblies. They were created by using the GEM mapper aligner (Derrien et al., 2012) belonging to the GEM suite (https://gemlibrary.sourceforge.net/), allowing up to two mismatches and considering sliding windows of 100-mer Targets

Target BED file can be downloaded from manufacturer's website and must follow this format ...

#### Output
**TargetPerla** module creates **three RData files per chromosome**, respectively containing its **mappability** (MAP), **local GCC content** (GCC) and (FRB) information from in-target and off-target reads.  
It **won't overwrite** output folder if already present, unless forced to do it with `-f` option.

> **R** version **>= 3.5.1** is needed to read RData files produced by each module. 


## 2. EXCAVATORDataPrepare
This module performs RC calculations, data normalization and data analysis on each .bam file listed in the .yaml file (e.g. samples.yaml).

#### Command Line Interface

> **ATTENTION**: Do **NOT** try to simultaneously execute two or more instances of **EXCAVATORDataPrepare** within the same output folder. Results would be unpredictable.

```sh
EXCAVATORDataPrepare.pl \
  -v -f \
  -s <samples file> \
  -t <target folder> \
  -o <output folder> \
  -@ <threads number> \
2> <log file>
```

Using singularity

```sh
singularity run \
  --bind home/workspace/src/exca2data/:srv/,\
  /home/workspace/db/:/db/,\
  /home/workspace/project1/file4exca2/:/file4exca2/,\
  /home/workspace/project1/bam/:/bam \
  excavator2.sif EXCAVATORDataPrepare.pl \
  -v -f \
  -s /file4exca2/ExpFilePrepare.yaml \
  -t output/hg19/1000GP/w_500000 \
  -o outputPrepare/ \
  -@ 4 \
2> logPrepare.txt
```

|**Option**        |**Description**                                  |
|------------------|-------------------------------------------------|
|`-h, --help`      |Print help message.                              |
|`-m, --man`       |Print complete documentation.                    |
|`-v, --verbose`   |Use verbose output.                              |
|`-f, --force`     |Force output overwrite.                          |
|`-s, --samples`   |Path to YAML [samples file](#samples-file).      |
|`-t, --target`    |Path to target folder.                           |
|`-o, --output`    |Path to output folder.                           |
|`-q, --mapq`      |Mapping quality threshold (default: 20).          |
|`-@, --threads`   |Number of threads to use (default: 1).           |

**Log** messages can be saved by redirecting **STDERR** to a file.

#### Samples file
It's a [YAML](https://en.wikipedia.org/wiki/YAML) file containing the paths to relevant **BAM/CRAM files** and their respective **label** (sample name). Sample labels must be unique.  
```yaml
# Test samples <Label>: <path to BAM/CRAM file>
Test1: /bam/test1.bam 
Test2: /bam/test2.bam
# ...
TestN: /bam/testN.bam
Control1: /bam/Control1.bam 
Control2: /bam/Control2.bam
# ...
ControlN: /bam/ControlN.bam
```

#### Output
**EXCAVATORDataPrepare** module creates a folder for each sample containing three subfolders (RC, RCNorm and Images) with, respectively, the calculated raw WMRC for In- and Off-target regions, the median normalized WMRC, and the plots showing the influence of GC content percentage and mappability on WMRC pre- and post-normalization.

## 3. EXCAVATORDataAnalysis
It performs the segmentation of the WMRC by means of the Heterogeneus Shifting Level Model algorithm and exploits FastCall algorithm to classify each segmented region as one of the five possible discrete states (2-copy deletion, 1-copy deletion, normal, 1-copy duplication and N-copy amplification). The FastCall calling procedure takes into account sample heterogeneity and use the Expectation Maximization algorithm to estimate the parameters of a five gaussian mixture model and to provide the probability that each segment belongs to a specific copy number state.

#### Command Line Interface

> **ATTENTION**: Do **NOT** try to simultaneously execute two or more instances of **EXCAVATORDataAnalysis** within the same output folder. Results would be unpredictable.

```sh
EXCAVATORDataAnalysis.pl \
  -v -f \
  -e <experimental design> \
  -s <samples file> \
  -i <input folder> \
  -t <target folder> \
  -o <output folder> \
  -@ <threads number> \
2> <log file>
```
Using singularity

```sh
singularity run \
  --bind home/workspace/src/exca2data/:srv/,\
  /home/workspace/db/:/db/,\
  /home/workspace/project1/file4exca2/:/file4exca2/,\
  /home/workspace/project1/bam/:/bam \
  excavator2.sif EXCAVATORDataAnalysis.pl \
  -vf -s /file4exca2/samples.yaml \
  -i outputPrepare/ \
  -t output/hg19/1000GP/w_500000 \
  -o outputAnalysis/ 
  -e paired 
  -@ 4 \
2> logAnalysis.txt
```

|**Option**        |**Description**                                      |
|------------------|-----------------------------------------------------|
|`-h, --help`      |Print help message.                                  |
|`-m, --man`       |Print complete documentation.                        |
|`-v, --verbose`   |Use verbose output.                                  |
|`-f, --force`     |Force output overwrite.                              |
|`-e, --experiment`|Experimental design ("pooled" or "paired").          |
|`-s, --samples`   |Path to YAML [samples file list](#samples-file-list).|
|`-t, --target`    |Path to target folder.                               |
|`-o, --output`    |Path to output folder.                               |
|`-p, --parameters`|Path to YAML [parameters file](#parameters-file) (default: EXACAVATORpath/parameters.yaml).    |
|`-@, --threads`   |Number of threads to use (default: 1).               |

**Log** messages can be saved by redirecting **STDERR** to a file.

#### Experimental settings 
EXCAVATOR2 can be used in experimental settings like **cancer genomic studies** in order to identify **somatic CNVs/CNAs** from pairs of tumor and matched normal samples. In this case, “paired” mode has to be set (-e option) and in *samples-file-list.yaml* test samples must be marked with a TX label (where X is an integer number) while control samples must be marked with CY (where Y is an integer number), whith X=Y. The analysis will be performed by comparing each test sample with its matched control sample.

In many applications, such as **germline CNV** detection or **population studies**, the control sample cannot be a proper matched sample. In this case, use “pooled” mode, so that each test sample will be compared with the same global control sample that results from pooling all samples labeled as controls by summing the WMRC region-by-region. 


#### Samples file list
It's a [YAML](https://en.wikipedia.org/wiki/YAML) file containing the **experimental group tags** and their respective **labels** (sample name) according to "-e option".   

```yaml
# Pooling control samples: 2 controls (C1, C2) - 1 test (T1)
T1: Test1 
C1: Control1
C2: Control2


# Paired control samples: 2 controls (C1, C2) - 2 test (T1, T2)
T1: Test1 
T2: Test2 
C1: Control1
C2: Control2

```

> Of note, to reduce technical biases, it is highly recommended to have control samples originating from sequencing libraries of the same type (e.g., paired-end, single-end) and processed by the same capture and sequencing protocols.


#### Parameters file
It's a YAML file storing running parameters **influencing the segmentation** as follows:  
- **Omega**: (default 0.1; range 0.0 to 1.0; suggested values 0.1–0.5) influences the evaluation of RC profile variation as real level shift, i.e., larger values favor the identification of slight variations over the stochastic noise signal;
- **Theta**: (default 1^−5; range 0.0 to 1.0) also influences the evaluation of profile level shift since it corresponds to the baseline probability to switch levels, thus controling specificity;
- **D_norm**: (default 10^5) modulates the ability of HSLM to detect both small and highly isolated coding regions of the genome and large and highly exon-covered genomic alterations.  

Parameters **influencing the classification** of each segmented region into five possible state (2-copy deletion, 1-copy deletion, normal, 1-copy duplication and N-copy amplification) and used by FastCall algorithm are as follows: 
- **Cellularity**: (default 1; range 0.0 to 1.0) is the fraction of tumor cells. This takes into account account the heterogeneity of tumor samples by accounting for the presence of normal cells at any proportion. It is 1 if no contamination of normal cells is present. The purity level of a tumor sample is usually measured experimentally, or can be computationally inferred from the data.  
- **d**: (default 0.5; suggested values 0.2 to 0.6) is the lower bound for the truncated gaussian calling the normal state (2 copies)
- **u**: (default 0.35; suggested values 0.1– 0.4) is the upper bound for the truncated gaussian calling the normal state (2 copies)
- **minExons**: (default 4) is the minimum number of exons to keep a segment in.

```yaml
# HSLM algorithm
HSLM:
  Omega:  0.1   # Omega parameter
  Theta:  1e-5  # Theta parameter (baseline probability m_i changes its value)
  D_norm: 10e5  # D_norm parameter

  
# FastCall Calling algorithm
FastCall:
  Cellularity: 1 # Cellularity parameter
  d: 0.5  # Threshold d for the truncated gaussian distribution
  u: 0.35 # Threshold u for the truncated gaussian distribution
  minExons: 4 # Segments with less exons than this threshold will be filtered out
```
