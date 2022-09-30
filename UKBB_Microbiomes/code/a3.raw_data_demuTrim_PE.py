# ========= import the module ======

import os
import argparse
import configparser
import subprocess

# ==================================
# read external argument
parser = argparse.ArgumentParser(description="Script for Demultiplexing")
parser.add_argument("-c", "--configfile", help="The config file.(default:./config.ini)",
                    default="./config.ini") #set the software location in this file
parser.add_argument("-o", "--outputDirectory", help="Output directory with results",
                    default="../results/1.Quality_Control/Demultiplexing")
parser.add_argument("-i", "--inputDirectory", help="input directory with fastq files",
                    default="../results/1.Quality_Control/joinedreads/")
parser.add_argument("-m", "--metadata", help="metadata",
                    default="./")    
parser.add_argument("-t", "--threads", help="Threads",
                    default="2")
parser.add_argument("-e", "--error", help="error rate",
                    default="1") # admit one error rate

# Might add something related to HPC

args = parser.parse_args()

# Process the command line arguments.
outputDirectory = os.path.abspath(args.outputDirectory)
configfile = os.path.abspath(args.configfile)
inputDirectory = os.path.abspath(args.inputDirectory)
threads = str(args.threads)
error =  str(args.error)
metadata =  os.path.abspath(args.metadata)


# get the software list by config file
config = configparser.ConfigParser()
config.read(configfile, encoding="utf-8")
python3 = config.get("software", "python3")
R = config.get("software", "R")
Trimmomatic = config.get("software", "Trimmomatic")
FastQC = config.get("software", "FastQC")


# check the directory
if os.path.exists(inputDirectory):
  print ("Input dir exists")
else :
  os.mkdir(inputDirectory)
if os.path.exists(outputDirectory):
  print ("Output dir exists")
else :
  os.mkdir(outputDirectory)

# write the command for QIIME2 preprocess

print ("============Start Importing============")
## import data
Command = "qiime tools import --type EMPPairedEndSequences \
    --input-path "+ inputDirectory + " \
    --output-path " + outputDirectory + "/emp-paired-end-sequences.qza"
# emp-paired-end-sequences.qza: barcodes and sequence fasta

#Command = Command + " && " + R + " ./Process.MetaData.R" #get the sample-metadata.csv
#subprocess.run(Command,shell=True,check=True)


print ("============Start Demultiplexing============")
## demultiplexing
Command1 = "qiime demux emp-paired \
    --i-seqs " + outputDirectory + "/emp-paired-end-sequences.qza \
    --m-barcodes-file " + metadata + "/sample-metadata.tsv \
    --m-barcodes-column BarcodeSequence   \
    --o-per-sample-sequences " + outputDirectory + "/demux.qza \
    --o-error-correction-details " +outputDirectory + "/demux-details.qza \
    --p-no-golay-error-correction"

## summarize and visualization
Command1 = Command1  + " && " + "qiime demux summarize \
    --i-data " + outputDirectory + "/demux.qza \
    --o-visualization " + outputDirectory + "/demux.qzv" #summarize
# demux.qza : a lot of demutiplexed fastq file in qza
# qiime tools view ../results/1.Quality_Control/Demultiplexing/demuz.qzv

#subprocess.run(Command1,shell=True,check=True)
# demux.qza : a lot of demutiplexed fastq file in qza 
# demux-details.qza: about the matching process detail of whole read record and barcode
# We donot care

## Trim the Primer // demux including detect the adapter/ However,not complete,therefore, we still need this
print ("============Start Trimming Adapter============")
Command2 = "qiime cutadapt trim-paired \
    --p-cores 4 \
    --i-demultiplexed-sequences " + outputDirectory + "/demux.qza \
    --p-anywhere-f AGAGTTTGATCCTGGCTCAG \
    --p-anywhere-r GCTGCCTCCCGTAGGAGT " + "--p-error-rate 0.1 \
    --o-trimmed-sequences " + outputDirectory + "/Primer_trimmed-seqs.qza"

Command2 = Command2  + " && " + "qiime demux summarize \
    --i-data " + outputDirectory + "/Primer_trimmed-seqs.qza \
    --o-visualization " + outputDirectory + "/Primer_trimmed-seqs.qzv" #summarize
subprocess.run(Command2,shell=True,check=True)




# finish - check - finish