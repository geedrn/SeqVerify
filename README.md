# SeqVerify
SeqVerify is a Python-based command line tool for analysis of whole genome sequencing data for gene-editing verification. It performs insertion site detection, copy number variation (CNV) analysis through CNVPytor, and bacterial contamination detection through KRAKEN2 and BRACKEN.

## Install

### Dependencies
* Python >=3.9
  * Regex
  * Argparse
* SAMtools >=1.14
* BWA >=0.7
* CNVPytor >=1.3
* Kraken2 >= 2.0
* Matplotlib >= 2.2

### Install through bioconda

```
> hello, world!
```

## Usage

An example SeqVerify call can be found below:

```
seqverify.py --output sample --reads_1 sample_1.fastq --reads_2 sample_2.fastq --genome genome.fa --marker_sources marker.fa --database db8gb
```

### Anatomy of a SeqVerify call

SeqVerify has the following required arguments:
* ```--output``` (Type: String) Used as the identifying name for all the files and folders to do with the particular call. E.g. ```--output Sample1``` will cause the output folder to be named "seqverify_Sample1"
* ```--reads_1``` and ```--reads_2``` (Type: String/Path) The paired-read FASTA/FASTQ, or gzipped FASTA/FASTQ source files for the reads. Also accepts paths to the files if they're not in the working folder.
* ```--genome``` (Type: String/Path) Name or path of FASTA file to be used as reference genome for the reads (e.g. [CHM13](https://github.com/marbl/CHM13#downloads)).
* ```--markers``` (Type: String/Path) Names or paths of FASTA files containing the genomes of the markers to detect (transgened, unwanted plasmids, etc.). Accepts more than one if necessary, space-separated. Can also be left blank.

SeqVerify has the following optional arguments:
##### Performance
* ```--threads``` (Type: Integer) Determines how many CPU threads are used in the multithreaded portion of the pipeline. Set to 1 by default.
* ```--max_mem``` (Type: String) Determines what the maximum amount of memory to be used in the indexing of the reference genome should be. Expects an integer followed by "M" or "G" (case-sensitive) for "Megabytes" or "Gigabytes" respectively. Set by default to "16G", 16 Gigabytes.
* ```--start``` (Type: Integer) For core-usage optimization. If set to 1, runs the entire pipeline from the start. If set to 0, only runs the single-threaded parts of the pipeline (i.e. up to and including the indexing of the reference genome, but stopping after that). If set to 2, only runs the multi-threaded parts of the pipeline (i.e. assumes the reference genome is already indexed and continues from there). Set to 1 by default, the logic behind this argument being that, where core optimization may be desirable, the pipeline can be split up into two commands, one that just requests a single thread, and one that requests multiple. Additionally, if one is using the same reference chromosome and markers for multiple samples, it is useless to re-index the genome every time, so setting 0 can be run once for the entire cohort, and setting 2 can be run after that for each individual sample.
##### KRAKEN2
* ```--kraken``` (Type: String/Boolean) Determines if KRAKEN2 analysis is to be performed on the sample or not. Set to False by default, if set to True requires ```--database``` option in order to work correctly.
* ```--database``` (Type: String/Path) Path to valid KRAKEN2 database, or, if KRAKEN2 environmental variables are set, the name of the database. Only needed if ```--kraken``` set to True. No default setting.
##### Insertion Site Detection
* ```--granularity``` (Type: Integer) Determines how large (in bp) insertion site bins are, i.e. how far apart two insertions can be in order to count as the same insertion site. Set by default to 500, a value of 1 means all insertions on different coordinates will be counted as different insertion sites and show up separately on the readout.
* ```--min_matches``` (Type: Integer) Determines the minimum number of matches/insertions required for an insertion site to appear on the readout. Set by default to 1, i.e. shows any insertion, some of which may be false positives due to repetitive DNA or similar variables.
##### CNVPytor
* ```--bin_size``` (Type: Integer) Determines how many bp are binned together for the purposes of copy number variation detection. Default is 100000, the minimum value is 1, but anything below the original read length (1500 in the test data) will yield meaningless data.

## Output

SeqVerify will output two folders, ```seqverify_output``` and ```seqverify_temp_output``` where "output" is the variable set in the ```--output``` argument. 

```seqverify_temp_output``` contains files that can be deleted after the pipeline is run, that are generated and used by the pipeline (SAM headers, .BED files, etc.).

```seqverify_output``` contains the following files:

```seqverify_output_markers.bam``` and ```seqverify_output_markers.bam.bai``` A BAM file containing the reads given in ```--reads_1``` and ```--reads_2``` realigned to the reference genome augmented with the marker genomes and its corresponding index.

```seqverify_readout.txt``` The human-readable readout resulting from the insertion site detection portion of the pipeline. The outermost indentation represents the marker, the second indentation is the human chromosome it aligns to, and the third layer denotes the in-chromosome coordinate of the insertion site, and how many times the marker was aligned to it (with the minimum value being specified in ```--min_matches```).

```output.pytor``` The binaries output by CNVPytor. Can be used to generate further plots of specific regions if needed (refer to CNVPytor docs)

```output.global.0000.png``` The Manhattan plot of CNV across the sample as output from ```output.pytor```.

```fig_marker.png``` A CNV histogram of the marker in the genome. Default bin size is 30bp, and one histogram is output for every marker specified in ```--markers```. 

If KRAKEN2 analysis is being performed, the pipeline will also output the following files in ```seqverify_output```:

```classified_seqs_output.kreport``` A human-readable report of the contamination detected by KRAKEN2 in text format.

```classified_output.kraken``` The KRAKEN2 binaries used to generate the report.

```classified_seqs_output_1.fq``` and ```classified_seqs_output_2.fq```, FASTQ files containing the sequences that were classified in the KRAKEN2 database.




