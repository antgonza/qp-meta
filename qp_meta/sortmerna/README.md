# rRNA removal

## SortMeRNA for Metatranscriptomics analysis

Updated  04/19/2020

SortMeRNA for removal of ribosomal reads from quality filtered metatranscriptomics data.

Latest SortMeRNA version: v2.1

SortMeRNA info available at https://bioinfo.lifl.fr/RNA/sortmerna/

Even if you are dealing with ribosomal RNA subtracted libraries, there will be some
residual ribosomal RNA in the libraries that you want to remove/separate from
the non ribosomal RNA sequences.

Input: Quality filtered metatranscriptomics reads (FASTA/FASTQ)

## Building Custom databases


Using the ARB package [2] to extract FASTA files for:
- 16S bacteria, 16S archaea and 18S eukarya using SSURef_NR99_119_SILVA_14_07_14_opt.arb
- 23S bacteria, 23S archaea and 28S eukarya using LSURef_119_SILVA_15_07_14_opt.arb.

Remove partial sequences from all databases using HMMER 3.1b1 [4] (meta_RNA [5] HMM profiles used)


Use SumaClust version 1.0.00 [6] to cluster sequences at various % id using the FASTA files generated in the previous step


Extract only the cluster centers



## Indexing custom databases

```indexdb_rna --ref db.fasta,db.idx ```

Multiple databases can be separated by ":"
indexdb_rna command available as a part of SortMeRNA package.


## Currently indexed databases

rRNA from bacteria, archaea, eukarya (SILVA and rfam database)

```
silva-bacterial-16s-id 90%
silva-bacterial-23s-id 98%
silva-archaeal-16s-id 95%
silva-archaeal-23s-id 98%
silva-eukarya-18s-id 95%
silva-eukarya-28s-id 98%
rfam-5s-database-id 98%
rfam-5.8s-database-id 98%
```


SILVA Databases can also be accessed at:
https://www.arb-silva.de/download/arb-files/

Rfam database:
https://rfam.xfam.org/
ftp://ftp.ebi.ac.uk/pub/databases/Rfam


## Usage

```sortmerna --ref $DB --reads $IN --aligned $RR --other $NR --paired_out --log --sam --fastx --blast 1 --num_alignments 1 -m 4096 -a 16```

## Description

```DIR="/opt/genome/qiita_dbs/sortmerna/rna-depletion/" ```

Path to where databases and their indexes are stored

```IN="/home/janedoe/data/DB/Test.fastq" ```

Input file, FASTA/FASTQ reads file, in this case interleaved fasta file

```RR="Test_ribosomal" ```

Prefix for output file conatining the identified ribosomal reads

```NR="Test_nonribosomal" ```

Prefix for output file conatining the identified non-ribosomal reads

```DB=./rRNA_databases/silva-bac-16s-id90.fasta,./index/silva-bac-16s-id90.idx:./rRNA_databases/silva-bac-23s-id98.fasta,./index/silva-bac-23s-id98.idx:./rRNA_databases/silva-arc-16s-id95.fasta,./index/silva-arc-16s-id95.idx:./rRNA_databases/silva-arc-23s-id98.fasta,./index/silva-arc-23s-id98.idx:./rRNA_databases/silva-euk-18s-id95.fasta,./index/silva-euk-18s-id95.idx:./rRNA_databases/silva-euk-28s-id98.fasta,./index/silva-euk-28s-id98.idx:./rRNA_databases/rfam-5s-database-id98.fasta,./index/rfam-5s-database-id98.idx:./rRNA_databases/rfam-5.8s-database-id98.fasta,./index/rfam-5.8s-database-id98.idx ```

Reference database(s) and their corresponding indexes separated by "," and multiple databases are separated by ":"

## Parameters ##


**--paired**  use only if you have an interleaved FASTA/Q file

**--sam**   option also outputs the SAM alignment files

**--aligned**   aligned reads filepath (in this case, ribosomal reads)

**--other**   rejected reads filepath + base file name (in this case, Non-ribosomal reads)

**--fastx**   output FASTA/FASTQ file

**-m**    Mbytes for loading reads into memory (default 1024)

**--num_alignments**  reports only the first n alignments per read reaching E-value

**--blast**   output alignments in various Blast-like formats                
   - 0   pairwise
   - 1   tabular (Blast -m 8 format)
   - 2   tabular + column for CIGAR
   - 3   tabular + columns for CIGAR and query coverage



## References

[1] Kopylova E., Noé L. and Touzet H., "SortMeRNA: Fast and accurate filtering of ribosomal RNAs in metatranscriptomic data", Bioinformatics (2012), doi:10.1093/bioinformatics/bts611.

[2] Wolfgang Ludwig et al., "ARB: a software environment for sequence data". Nucleic Acids Research (2004) 32:1363-1371

[3] Sarah W. Burge et al., "Rfam 11.0: 10 years of RNA families". Nucleic Acids Research (2012) doi: 10.1093/nar/gks1005

[4] R.D. Finn, J. Clements and S.R. Eddy, "HMMER web server: interactive sequence similarity searching". Nucleic Acids Research (2011) Web Server Issue 39:W29-W37

[5] Ying Huang, Paul Gilna and Weizhong Li, "Identification of ribosomal RNA genes in metagenomic fragments". Bioinformatics (2009) 25:1338-1340

[6] Céline Mercier et al., SUMATRA and SUMACLUST: fast and exact comparison and clustering of full-length barcode sequences.
