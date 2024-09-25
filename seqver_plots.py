#Defines functions necessary for the plotting of transgene copy numbers, either through IGV or matplotlib.
import os
import matplotlib.pyplot as plt
from collections import OrderedDict
import logging
import traceback

import logging

def region_bed(temp_folder, sam_header, commands, chr_list, output, zoomout=200):
    logging.info(f"Starting region_bed function with parameters:")
    logging.info(f"temp_folder: {temp_folder}, sam_header: {sam_header}, output: {output}, zoomout: {zoomout}")
    logging.info(f"chr_list: {chr_list}")
    logging.info(f"commands: {commands}")

    try:
        with open(f"{temp_folder}/{output}", "w+") as bed:
            if chr_list:
                logging.info("Processing chromosome list from SAM header")
                with open(f"{temp_folder}/{sam_header}") as file:
                    for line in file:
                        tags = line.strip().split("\t")
                        if tags[0] == '@SQ':
                            chr_name = tags[1].split(":")[1]
                            if chr_name in chr_list:
                                end = int(tags[2].split(":")[1]) - 1
                                bed_line = f"{chr_name}\t1\t{end}\n"
                                bed.write(bed_line)
                                logging.info(f"Added chromosome line: {bed_line.strip()}")
            
            if commands is not None:
                logging.info("Processing commands")
                for command in commands:
                    fields = str(command).split("\t")
                    chr, start, seq = fields[0].split(":")[0], int(fields[0].split(":")[1].split("-")[0]), len(fields[1])
                    if start > zoomout:
                        start -= zoomout
                    seq += zoomout
                    bed_line = f"{chr}\t{str(start)}\t{int(start)+int(seq)}\n"
                    bed.write(bed_line)
                    logging.info(f"Added command line: {bed_line.strip()}")

        logging.info(f"BED file creation completed: {temp_folder}/{output}")
        
        # DEBUG: Read the content of the BED file
        with open(f"{temp_folder}/{output}", "r") as f:
            content = f.read()
            logging.info(f"BED file content:\n{content}")
        
        return True
    except Exception as e:
        logging.error(f"Error in region_bed function: {str(e)}")
        return False

def histogramData(coveragemap, chromosome, granularity=1): #Collects data to make a single histogram if IGV is not used
    os.system(f"gawk '{{if ($1 ~ /({chromosome})\>/) print $0}};' {coveragemap} > {chromosome}_coverage.cov") #gawks the output of samtools depth -b for the reads relevant to the chromosome
    with open(f"{chromosome}_coverage.cov") as file:
        bins, counts = [], [] #initializes an empty list of bins and an empty list of read numbers per bin 
        for line in file:
            features = line.split("\t")
            try:
                position, readnum = int(features[1].strip()), int(features[2].strip()) #takes the position and number of reads at that position for each line in the file
            except:
                continue
            try:
                if abs(bins[-1]-position) < granularity: #if the bin that was just analyzed falls within a previous bin, add the number of reads to the previous bin
                    counts[-1] += readnum
                else: #if not, make a new bin and a new corresponding number of reads
                    bins.append(position)
                    counts.append(readnum)
            except:
                    bins.append(position)
                    counts.append(readnum)
    counts = [int(i)/granularity for i in counts] #averages the readnums over the length of the bins as to indicate average read number
    return bins,counts

def histogram(bins,counts,chromosome,granularity=50,rounding=1, significant=1): #generates the plot if using matplotlib
    try:
        abnormal_bins, abnormal_counts = [], [] #initializes bins and read counts for any potentially abnormal bins

        mean = sum(counts)/len(counts) #calculates the average number of reads in the transgene overall

        counts = [int(round((i/mean)*rounding,0)) for i in counts] #calculates the mean copy number for all the bins 

        mean_cnv = sum(counts)/len(counts) #calculates the mean copy number overall

        bins = [str(int(round(i/granularity,0))) for i in bins] #calculates the new bins
        
        std_dev = (sum([(i-mean_cnv)**2 for i in counts])/len(counts))**0.5 #calculates the overall standard deviation in the copy numbers
        for index in range(len(counts)):
            count = counts[index]
            if abs(count-mean_cnv) >= significant*std_dev: #if the difference between a copy number and the mean copy number is larger than the threshold coefficient multiplied by the standard deviation, flags the bin as abnormal
                abnormal_bins.append(bins[index])
                abnormal_counts.append(counts[index])
                counts[index] = 0

        fig,ax = plt.subplots() #matplotlib logic
        ax.bar(bins,counts,color ='blue') #plots all normal bins in blue
        ax.bar(abnormal_bins,abnormal_counts,color ='green') #plots all abnormal/significant bins in green
        ax.axhline(y = mean_cnv, color = 'r', linestyle = '-') #plots a red average line

        #Defines labels, saves the image, closes the plot
        plt.xlabel(f"Base Coordinate (/{granularity})") 
        plt.ylabel(f"Avg. Copy Number per Bin (x{rounding})")
        plt.xticks([])
        plt.title("Chromosome Coverage")
        fig.savefig(f'fig_{chromosome}.png')
        plt.close('all')
    except:
        pass

def chrHistograms(coveragemap, chrList): #Generates histograms for every chromosome in chrList by calling the above two functions each time.
    for i in range(len(chrList)):
        try:
            bins, counts = histogramData(coveragemap, chrList[i])
            histogram(bins, counts, chrList[i])
        except ZeroDivisionError:
            continue

def igvScreenshot(temp_folder,folder,alignments,genome,bed_file,imageformat="png",gtf_file=None): #If using IGV, deals with IGV logic
    with open(f"{temp_folder}/seqverify_igv.bat","w+") as file: #Autogenerates an IGV-compatible bat file we will use later to screenshot the relevant parts
        file.write("new\n") #boilerplate code
        file.write('echo "very_beginning"\n')
        file.write(f"snapshotDirectory {folder}\n") #sets the folder for the screenshots
        file.write('echo "beginning"\n')
        file.write(f"genome {genome}\n") #loads the genome
        file.write(f'echo "genome {genome} loaded"\n')
        file.write(f"load {alignments}\n") #loads the bam file
        file.write(f'echo "alignments loaded"\n')
        if gtf_file is not None:
            file.write(f"load {gtf_file}\n")
            file.write(f'echo "gtf file loaded"\n')
        file.write(f"maxPanelHeight 500\n") #boilerplate code for adjustment of the screen size
        with open(f"{temp_folder}/{bed_file}","r") as bed: #writes instructions to take a screenshot of every transgene in the bed file
            for line in bed:
                name,begin,end = [i.strip() for i in line.split("\t")]
                file.write(f'echo "snapshot {name}"\n')
                file.write(f"goto {name}:{begin}-{end}\n") #makes IGV go to the entire transgene
                file.write(f'echo "goto {name}"\n')
                file.write(f"snapshot fig_{name}.{imageformat}\n") #takes screenshot and saves it
                file.write(f'echo "saved snapshot {name}"\n')
        file.write("exit") #boilerplate
    igv_cmd = f"xvfb-run --auto-servernum --server-args=\"-screen 0, 2048x1536x24\" igv -b {temp_folder}/seqverify_igv.bat"
    print(igv_cmd) # for debug
    os.system(igv_cmd) #runs XVFB, a headerless server emulator, to run IGV automatically without the need for a GUI.
    
def igvScreenshot_new(temp_folder, folder, alignments, genome, bed_file, imageformat="png", gtf_file=None):
    try:
        # DEBUG: bed file path
        if not os.path.exists(f"{temp_folder}/{bed_file}"):
            raise FileNotFoundError(f"BED file not found: {temp_folder}/{bed_file}")

        # DEBUG: genome file path
        if not os.path.exists(genome):
            raise FileNotFoundError(f"Genome file not found: {genome}")

        # DEBUG: alignment file path
        if not os.path.exists(alignments):
            raise FileNotFoundError(f"Alignment file not found: {alignments}")

        base_cmd = f"create_report {temp_folder}/{bed_file} --fasta {genome} --standalone --flanking 1000 --sequence 1 --begin 2 --end 3 --tracks {alignments}"

        if gtf_file is not None:
            if not os.path.exists(f"{folder}/{gtf_file}"):
                raise FileNotFoundError(f"GTF file not found: {folder}/{gtf_file}")
            cmd_str = f"{base_cmd} {folder}/{gtf_file} --output {folder}/igv_viewer.html"
        else:
            cmd_str = f"{base_cmd} --output {folder}/igv_viewer.html"

        logging.info(f"Executing command: {cmd_str}")
        print(f"Executing command: {cmd_str}")

        result = os.system(cmd_str)

        if result != 0:
            raise Exception(f"Command execution failed with exit code: {result}")

        if not os.path.exists(f"{folder}/igv_viewer.html"):
            raise FileNotFoundError(f"Output file not created: {folder}/igv_viewer.html")

        logging.info("IGV screenshot generation completed successfully")
        return True

    except Exception as e:
        logging.error(f"Error in igvScreenshot_new: {str(e)}")
        logging.error(f"Traceback: {traceback.format_exc()}")
        print(f"Error in igvScreenshot_new: {str(e)}")
        return False

def genome_configurator(temp_folder,pytor_conf,gc_name,genome,sam_header):
    special_chrs = {"chrX":"S","chrY":"S","chrM":"M"}
    genome_dict = {"seqverify_genome":{"name":"seqverify","species":"human","chromosomes":OrderedDict(),"gc_file":os.getcwd()+"/"+temp_folder+"/"+gc_name}}

    os.system(f"cnvpytor -root {temp_folder}/{gc_name} -gc {temp_folder}/{genome} -make_gc_file")
    
    with open(f"{temp_folder}/{sam_header}","r") as header:
        for line in header:
            tags = line.split("\t")
            if tags[0] == "@SQ":
                chr_name, chr_len = tags[1].split(":")[1], int(tags[2].split(":")[1].strip())
                if chr_name in special_chrs:
                    genome_dict["seqverify_genome"]["chromosomes"][chr_name] = tuple((chr_len,special_chrs[chr_name]))
                else:
                    genome_dict["seqverify_genome"]["chromosomes"][chr_name] = tuple((chr_len,"A"))
    
    with open(f"{temp_folder}/{pytor_conf}","w+") as conf:
        conf.write(f"import_reference_genomes = {repr(genome_dict)}")
