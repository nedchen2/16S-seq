"""
This script provide dada2 denoise approach with Paired-end Data imported  
"""


# ========= import the module ======

import os
import argparse
import configparser
import subprocess

# ==================================
# read external argument
parser = argparse.ArgumentParser(description="Script for getting feature table")
parser.add_argument("-c", "--configfile", help="The config file.(default:./config.ini)",
                    default="./config.ini") #set the software location in this file
parser.add_argument("-o", "--outputDirectory", help="Output directory with results",
                    default="../results/1.Quality_Control/Demultiplexing")
parser.add_argument("-i", "--inputDirectory", help="input directory with fastq files",
                    default="../results/2.Feature_table/")
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
 
# get the software list by config file
config = configparser.ConfigParser()
config.read(configfile, encoding="utf-8")
python3 = config.get("software", "python3")
R = config.get("software", "R")
Trimmomatic = config.get("software", "Trimmomatic")
FastQC = config.get("software", "FastQC")

# write the command for QIIME2 preprocess

# check the directory
if os.path.exists(inputDirectory):
  print ("Input dir exists")
else :
  os.mkdir(inputDirectory)
if os.path.exists(outputDirectory):
  print ("Output dir exists")
else :
  os.mkdir(outputDirectory)

# ======================== start 

print ("============Start Demultiplexing Using dada2============")

Command = "qiime dada2 denoise-paired \
    --i-demultiplexed-seqs " + outputDirectory + "/Primer_trimmed-seqs.qza \
    --p-trunc-len-f 225 \
    --p-trunc-len-r 225 \
    --p-trim-left-f 0 \
    --p-trim-left-r 0 \
    --p-n-threads 4 \
    --o-representative-sequences  " + inputDirectory + "/rep-seqs.qza \
    --o-table " + inputDirectory + "/table.qza \
    --o-denoising-stats " + inputDirectory + "/denoising-stats.qza"

subprocess.run(Command,shell=True,check=True)