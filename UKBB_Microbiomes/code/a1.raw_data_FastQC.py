# Author : Congjia Chen

# ========= import the module ======

import os
import argparse
import configparser
import subprocess

# ==================================
# read external argument
parser = argparse.ArgumentParser(description="Script for FastQC")
parser.add_argument("-c", "--configfile", help="The config file.(default:./config.ini)",
                    default="./config.ini") #set the software location in this file
parser.add_argument("-o", "--outputDirectory", help="Output directory with results",
                    default="../results/1.Quality_Control/FastQC")
parser.add_argument("-i", "--inputDirectory", help="input directory with fastq files",
                    default="../16sRaw")
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

# write the command for Fastqc
Command = FastQC + " " + "-o" + " " + outputDirectory + " " + "-t" + " " + threads + " " + "-q" + " " + inputDirectory + "/*.fastq.gz"
# if do not want to zip the file add "--extract" to the command

# run 
subprocess.run(Command,shell=True,check=True)


## the quality in the front and behind is not good. Therefore, need to be trimmed