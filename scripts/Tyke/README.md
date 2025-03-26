# TykeVAR

## Installation

These instructions are valid for **Linux x86**. For other platforms (e.g. MacOS), instructions will need to be adapted.

Load conda, e.g.
```
module load miniconda/3
```
or (in some instances python has to be loaded to run the conda command)
```
module load anaconda3/2024.02
module load python
```

Installing two conda environments with python 3.10
```
conda create -n MosaicSim python=3.10
conda init
conda activate mosaicSim
```

Unload python if previously loaded so that python 3.10 from conda environment will be used
```
module unload python
```

Obtain the tool from github (replace $HOME with your preferred installation director)
```
cd  $HOME
git clone https://github.com/erikstricker/MosaicSim.git
```

To install the relevant python dependencies, run
```
REPO_ROOT="$HOME/SpikeVarTykeVar/"

pip install -r $REPO_ROOT/requirements.txt
```
Ensure to also load samtools >1.15.1 and bcftools>1.19

Once the requirements are installed, please install or load the following additional packages

_Installation_
```
conda install -c bioconda samtools bcftools
```
_Loading (e.g.)_
```
export PATH=/path/to/software/samtools/samtools-1.21/bin:$PATH
export PATH=/path/to/software/bcftools/bcftools-1.19/bin:$PATH
```
or
```
module load samtools-1.21
module load bcftools-1.19
```

## Dependencies
  
### TykeVar
- pysam (0.21.0) 
- numpy (1.25.2)
- biopython (1.81)
- samtools 1.15.1
- bcftools 1.19

## How to use it

### TykeVar

#### 1) Generate simulated VCF

The VCF simulator generates a random set of mosaic variants (SNVs and SVs). The variants
can be parameterized with VAF, number of variants to simulate and the size of the variations.
The generated file is in the Sniffles VCF format.

The variants generated here act as an input into the read editor step (described below) which
generates modified reads with the variants inserted into them. The same VCF file is also the
ground truth for evaluating mosaic variant callers.

```
python vcfgen.py <path_to_bam> <path_to_ref> <output_path_prefix>

e.g. python vcfgen.py chr22.bam hs37d5.fa chr22
The above generates a chr22SV.vcf and chr22SNV.vcf file
```

#### 2) Generate edited reads based on simulated VCF

```
python main.py -v <SIMULATED_VCF> -b <BAM> -r <REF> -o <OUTPUT_FASTQ>
```

This command above takes in the VCF which determines which variants to introduce into the reads.
The BAM file is used to find the reads which overlap with variant locations. Only a subset of the reads
corresponding to a particular variant location are edited. This is determined by the allele frequency.
The output FASTQ file has the edited reads. The query name of each read is kept the same.

## Example Implementation

### TykeVar

Here, we use the TykeVar workflow to modifiy reads of HG002 directly at their reference position by including artifical mutations to represent at variant allele frequency of 5%. In contrast to the above approach we do not introduce new haplotypes with this. However, more complex mutations (e.g. rearrangements, duplication or very long structural variants) will not be able to be introduced to the data itself, since the size of the reads is limited.


#### 1) Fetch data
In order to simulate and edit reads, the pipeline first needs an initial set of aligned reads and a reference. For our demonstration, we will use the GIAB datasets.

Reads - `ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/Ultralong_OxfordNanopore/guppy-V3.2.4_2020-01-22/HG002_hs37d5_ONT-UL_GIAB_20200122.phased.bam` and `ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/Ultralong_OxfordNanopore/guppy-V3.2.4_2020-01-22/HG002_hs37d5_ONT-UL_GIAB_20200122.phased.bam.bai`

Reference - `https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/references/GRCh37/hs37d5.fa.gz`

#### 2) Generate variants and modified reads

First we decompress the FASTA file.
```
gunzip hs37d5.fa.gz -c hs37d5.fa
```

Then we simulate variants
```
python vcfgen.py HG002_hs37d5_ONT-UL_GIAB_20200122.phased.bam hs37d5.fa hg002
```

Generate a set of modified reads with inserted variants.
```
python main.py -v hg002SV.vcf -b HG002_hs37d5_ONT-UL_GIAB_20200122.phased.bam -r hs37d5.fa -o hg002_modified_reads.fastq
```

#### 3) Re-align modified reads and merge them
Once the new reads are generated, they need to be re-aligned and re-inserted back into
the dataset by replacing the original reads.

Installation:
```console
# install samtools, bwa-mem2, minimap2
mkdir -p ~/tools && cd ~/tools
wget https://github.com/samtools/samtools/releases/download/1.18/samtools-1.18.tar.bz2
tar jxf samtools-1.18.tar.bz2 && cd samtools-1.18
./configure --prefix=$(pwd)
make && make install
export PATH=$HOME/tools/samtools-1.18/:$PATH

cd ~/tools
curl -L https://github.com/bwa-mem2/bwa-mem2/releases/download/v2.2.1/bwa-mem2-2.2.1_x64-linux.tar.bz2 | tar jxf -
export PATH=$HOME/tools/bwa-mem2-2.2.1_x64-linux/:$PATH

curl -L https://github.com/lh3/minimap2/releases/download/v2.26/minimap2-2.26_x64-linux.tar.bz2 | tar -jxvf -
export PATH=$HOME/tools/minimap2-2.26_x64-linux/:$PATH
```

Indexing and mapping:
```console
# indexing reference genome
bwa-mem2 index hs37d5.fa.gz
minimap2 -x map-ont -d mm2_hs37d5.mmi hs37d5.fa.gz

# mapping 
## short read
bwa-mem2 mem -t 14 ~/reference/hs37d5.fa.gz samp.fastq -o bwa_align/test.sam
samtools view -bS --no-PG bwa_align/test.sam | samtools sort -@ 12 --no-PG - > bwa_align/test.sorted.bam

## long read
minimap2 -a ~/reference/mm2_hs37d5.mmi -t 14 --secondary=no chr22.fastq | samtools view -bS --no-PG - > mod.chr22.bam
```

Postprocessing bam:
```
python filter_merge_bam.py -b chr22.HG002_hs37d5_ONT-UL_GIAB_20200122.phased.bam \
    -m mod.chr22.bam -o . --prefix mod_chr22 --primary
```

#### 4) Run your favorite mosaic variant caller

TODO:
