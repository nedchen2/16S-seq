# after preprocessing, we will do pairs join
# QIIME1 fastq-join
# Extract the barcode_in_label

# ========= import the module ======

import os
import argparse
import configparser
import subprocess

# ==================================
# read external argument
parser = argparse.ArgumentParser(description="Script for Join reads")
parser.add_argument("-c", "--configfile", help="The config file.(default:./config.ini)",
                    default="./config.ini") #set the software location in this file
parser.add_argument("-o", "--outputDirectory", help="Output directory with results",
                    default="../results/1.Quality_Control/joinedreads/")
parser.add_argument("-i", "--inputDirectory", help="input directory with fastq files",
                    default="../16sRaw/")
parser.add_argument("-t", "--threads", help="Threads",
                    default="2")

# Might add something related to HPC

args = parser.parse_args()

# Process the command line arguments.
outputDirectory = os.path.abspath(args.outputDirectory)
configfile = os.path.abspath(args.configfile)
inputDirectory = os.path.abspath(args.inputDirectory)
threads = str(args.threads)

# get the software list by config file
config = configparser.ConfigParser()
config.read(configfile, encoding="utf-8")
python3 = config.get("software", "python3")
R = config.get("software", "R")
Trimmomatic = config.get("software", "Trimmomatic")
Trimmomatic_adapter = config.get("software", "Trimmomatic_adapter")
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

# we need single file name;

R1 = "lane1_Undetermined_L001_R1_001.fastq.gz"
R2 = "lane1_Undetermined_L001_R2_001.fastq.gz"

# write the qiime1 to deal with the fastq
Command1 = "join_paired_ends.py" + " -f "+ inputDirectory + "/" + R1 + " -r " + inputDirectory + "/" + R2 +  " -m " + "fastq-join" + " -o " + outputDirectory
Command2 = "extract_barcodes.py" +  " --input_type barcode_in_label " + " -f " + outputDirectory + "/fastqjoin.join.fastq" + " --bc1_len 14 " + " -o " + outputDirectory + "/" # some of the qrgument here could be edit in other experiments
Command3 = "mv " + outputDirectory + "/fastqjoin.join.fastq " + outputDirectory + "/sequences.fastq" + " && " + " rm " + outputDirectory +"/fastqjoin.un*.fastq"  +" && " + " gzip -r " + outputDirectory #rename and zip


# run 
subprocess.run(Command1,shell=True,check=True)
subprocess.run(Command2,shell=True,check=True)
subprocess.run(Command3,shell=True,check=True)
