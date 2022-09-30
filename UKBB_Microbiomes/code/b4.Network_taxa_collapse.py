
# ========= import the module ======

from configparser import BasicInterpolation
from cProfile import run
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
                    default="../results/6.Network_analysis")
parser.add_argument("-i", "--inputDirectory", help="input directory with fastq files",
                    default="../results/4.Diversity_ana")
parser.add_argument("-m", "--metadata", help="metadata",
                    default="./")     
parser.add_argument("-t", "--threads", help="Threads",
                    default="3")
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
classifier =  os.path.abspath(args.metadata) # dir of the pre-trained classifier
 
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

#===========may be add some taxa collapse here 
print ("=============Start Taxa Collapsing=============")

Command1 = "qiime taxa collapse \
      --i-table " + inputDirectory + "/core-metrics-results/rarefied_table.qza\
      --i-taxonomy " + inputDirectory + "/taxonomy-corrected.qza\
      --p-level 6\
      --o-collapsed-table " + outputDirectory + "/table-l6.qza"

#subprocess.run(Command1,shell=True,check=True)


print ("=============Start Species Taxa Collapsing=============")

Command2 = "qiime taxa collapse \
      --i-table " + inputDirectory + "/core-metrics-results/rarefied_table.qza\
      --i-taxonomy " + inputDirectory + "/taxonomy-corrected.qza\
      --p-level 7\
      --o-collapsed-table " + outputDirectory + "/table-l7.qza"

#subprocess.run(Command2,shell=True,check=True)


print ("=============Start Species Novel Taxa Collapsing=============")

Command3 = "qiime taxa collapse \
      --i-table ../results/6.Network_analysis/core-metrics-results/rarefied_table.qza\
      --i-taxonomy ../results/6.Network_analysis/taxonomy-corrected.qza\
      --p-level 7\
      --o-collapsed-table " + outputDirectory + "/table-l7.qza"

subprocess.run(Command3,shell=True,check=True)



