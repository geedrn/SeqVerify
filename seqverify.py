#!/usr/bin/env python
 
import argparse
import os
from seqver_functions import *
from seqver_transgenes import *

parser = argparse.ArgumentParser()
#python magnify_commandline.py --output S04 --reads_1 S04_1.fastq --reads_2 S04_2.fastq --genome chm13v2.0.fa --markers T2A-mGreenLantern.fa T2A-tdTomato.fa transposons_in_S04.fa unwanted_plasmids.fa --granularity 500 --threads 15
parser.add_argument('--output', type=str, required=True)
parser.add_argument('--reads_1', type=str, required=True)
parser.add_argument('--reads_2', type=str, required=True)
parser.add_argument('--genome', type=str, required=True)
parser.add_argument('--markers', type=str, required=False, nargs='+')

parser.add_argument('--kraken', type=str, required=False,default="F")
parser.add_argument('--database',type=str,required=False)

parser.add_argument('--granularity',type=int, required=False, default=500)
parser.add_argument('--threads',type=int,required=False, default=1)
parser.add_argument('--max_mem', type=str, required=False, default='16G')
parser.add_argument('--min_matches', type=int, required=False,default=1)
parser.add_argument('--start', type=int, required=False, default=1)
parser.add_argument('--del_temp', type=int, required=False, default=1)

parser.add_argument('--bin_size',type=int, required=False, default=100000)

args = parser.parse_args()

output_name = args.output

#The pathFinder function detects whether arguments are paths or names and acts accordingly 
#such that it can handle both cases for ease of use.
reads_1_path = pathFinder(args.reads_1)
reads_2_path = pathFinder(args.reads_2)
genome_source_path = pathFinder(args.genome)
marker_sources_path = [pathFinder(i) for i in args.markers]

granularity_threshold = args.granularity
database = pathFinder(args.database)

min_matches = args.min_matches
bin_size = args.bin_size
threads = args.threads

if args.bin_size != 100000:
    if args.bin_size % 100 == 0:
        bin_size = args.bin_size
    else:
        raise ValueError("Custom bin size not divisible by 100")

#Deals with memory logic, converts max_mem argument into number of bytes of memory requested in order to build optimal BWA indexing block size
if args.max_mem[-1] == "M":
    max_mem = int(args.max_mem[:-1])*1000000
elif args.max_mem[-1] == "G":
    max_mem = int(args.max_mem[:-1])*1000000000
elif type(args.max_mem) == 'int':
    max_mem = int(args.max_mem)

base_block = max_mem / 8

#Stores Variable, Folder Names:

folder = f"seqverify_{output_name}"
temp_folder =f"seqverify_temp_{output_name}"
marker_list_file = f"seqverify_{output_name}_transgene_list.fa"
genome_with_markers = f"seqverify_{output_name}_collated.fa"
aligned_with_markers = f"seqverify_{output_name}_markers.sam"
aligned_diff_chr = f"seqverify_{output_name}_markers_diff_chr.sam"
aligned_diff_chr_bam = f"seqverify_{output_name}_markers_diff_chr.bam"
header = f"seqverify_{output_name}_markers_header.sam"
overall_coverage = f"seqverify_{output_name}_markers_coverage.cov"
pytor = f"{output_name}.pytor"
unaligned = f"seqverify_{output_name}_unaligned.sam"
unaligned_R1 = f"seqverify_{output_name}_unaligned_R1.fastq"
unaligned_R2 = f"seqverify_{output_name}_unaligned_R2.fastq"

if args.start == 1 or args.start == 0:
    #Preliminary Work - Generates the output folder given a specified --output string
    #                   as well as the temp folder for the pipeline to work within
    os.system(f'mkdir {folder}')
    os.system(f'mkdir {temp_folder}')

    #Prep Collated Genome


    #Collates all markers/transgenes into a single FASTA file for ease of handling
    for file in marker_sources_path: 
        os.system(f'echo "$(cat {file})" >> {temp_folder}/{marker_list_file}') 

    #Copies the reference genome into a new FASTA file that the transgenes will also be added to
    os.system(f'cp {genome_source_path} {temp_folder}/{genome_with_markers}')

    #Adds the transgenes to the new "collated" FASTA file
    os.system(f'echo "$(cat {temp_folder}/{marker_list_file})" >> {temp_folder}/{genome_with_markers}')

    #Builds the indices for the genome that has now been augmented with the transgenes, necessary for alignment
    os.system(f'bwa index -b {base_block} {temp_folder}/{genome_with_markers}')

if args.start == 1 or args.start == 2:
    #Beginning of multithreaded part -------------------------------------------------------------------------------

    #Aligns the new genome using BWA-MEM
    os.system(f'bwa mem -M -t {threads} {temp_folder}/{genome_with_markers} {reads_1_path} {reads_2_path} > {temp_folder}/{aligned_with_markers}')

    #Saves all names of transgenes/markers we want to find insertion sites of into a list for later usage and makes it into a string
    with open(f'{temp_folder}/{marker_list_file}') as f:
        marker_list = [i[1:].strip() for i in f if i.startswith('>')]
    markers = "|".join(marker_list)

    #Uses the "markers" string as a regular expression to find reads in the newly-aligned SAM file that align to a transgene and have a mate on a human chromosome
    #Exports these reads into a new SAM file, "diff_chr"
    os.system(f"samtools view -h -@ {threads} {temp_folder}/{aligned_with_markers} |gawk '{{if ((($3 ~ /^({markers})/)&&($7 ~ /chr([0-9]+|[XYM])/))||($1 ~ /^@/)) print $0}};' > {temp_folder}/{aligned_diff_chr}")

    #Functions imported from magnify_functions, refer to comments there for details on their functioning.
    #Takes "diff_chr", generates a .txt human-readable readout of the insertion sites, places it into the output folder.
    data = group(f'{temp_folder}/{aligned_diff_chr}')
    insertions = compress(data,granularity_threshold)
    readout(folder,insertions,marker_list,min_matches)

    os.system(f"samtools sort -@ {threads} {temp_folder}/{aligned_with_markers} > {folder}/{aligned_diff_chr_bam}")
    os.system(f"samtools index -@ {threads} {folder}/{aligned_diff_chr_bam}")

    os.system(f'cnvpytor -root {folder}/{pytor} -rd {folder}/{aligned_diff_chr_bam}')
    os.system(f'cnvpytor -root {folder}/{pytor} -his {bin_size}')
    os.system(f'cnvpytor -root {folder}/{pytor} -partition {bin_size}')
    os.system(f'cnvpytor -root {folder}/{pytor} -call {bin_size} > {folder}/calls.{bin_size}.tsv')
    os.system(f'cnvpytor -root {folder}/{pytor} -plot manhattan {bin_size} -o {folder}/{output_name}.png')

    os.system(f'samtools view -@ {threads} -H {folder}/{aligned_diff_chr_bam} > {temp_folder}/{header}')

    bed_file = region_bed(f"{temp_folder}/{header}",marker_list)

    os.system(f'samtools depth -b {bed_file} {folder}/{aligned_diff_chr_bam} > {temp_folder}/{overall_coverage}')

    chrHistograms(f'{temp_folder}/{overall_coverage}',marker_list)
    os.system(f'mv fig_* {folder}')
    os.system(f'mv *.cov {temp_folder}')

    #KRAKEN2/BRACKEN system
    if args.kraken[0].lower() == "t":
        os.system(f'samtools view -f 13 -@ {threads} {folder}/{aligned_diff_chr_bam}> {temp_folder}/{unaligned}')

        os.system(f'samtools fastq -@ {threads} {temp_folder}/{unaligned} -1 {temp_folder}/{unaligned_R1} -2 {temp_folder}/{unaligned_R2}')

        os.system(f'kraken2 --threads {threads} --db {database} --report {folder}/classified_seqs_{output_name}.kreport --paired --classified-out {folder}/classified_seqs_{output_name}#.fq {temp_folder}/{unaligned_R1} {temp_folder}/{unaligned_R2} > {folder}/classified_{output_name}.kraken')

        os.system(f'bracken -d {database} -i {folder}/classified_seqs_{output_name}.kreport -o {folder}/classified_seqs_{output_name}.bracken -r 150')

    if args.del_temp == 1:
        os.system(f'rm -r {temp_folder}')