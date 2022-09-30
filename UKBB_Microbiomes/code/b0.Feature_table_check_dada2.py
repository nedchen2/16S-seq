# ========= import the module ======

import os
import argparse
import configparser
import subprocess

# ==================================
# read external argument
parser = argparse.ArgumentParser(description="Script for feature table analysis")
parser.add_argument("-c", "--configfile", help="The config file.(default:./config.ini)",
                    default="./config.ini") #set the software location in this file
parser.add_argument("-o", "--outputDirectory", help="Output directory with results",
                    default="../results/1.Quality_Control/Demultiplexing")
parser.add_argument("-i", "--inputDirectory", help="input directory with fastq files",
                    default="../results/2.Feature_table/")
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

## Visualization
# The feature-table summarize command will give you information on 
# how many sequences are associated with each sample and with each feature, 
# histograms of those distributions, and some related summary statistics. 

# feature-table tabulate-seqs provide a mapping of feature IDs to sequences, 
# and provide links to easily BLAST each sequence against the NCBI nt database. 


# ========================== filtering ===================

Command0 = "qiime feature-table filter-features --i-table " + inputDirectory + "/table.qza  \
      --p-min-samples 1\
      --o-filtered-table" + inputDirectory + "/table_filter.qza"

Command0 = Command0 + "qiime feature-table filter-seqs --i-data " + inputDirectory + "/rep-seqs.qza  \
      --i-table " + inputDirectory + "/table_filter.qza"

#subprocess.run(Command0,shell=True,check=True)
# add some quality control here


# ==========================================================

# ========= test - chimera 
print ("=============Start Test Chimera=============")

Command2 = "qiime vsearch uchime-denovo \
  --i-table " + inputDirectory + "/table.qza \
  --i-sequences " + inputDirectory + "/rep-seqs.qza\
  --output-dir " + inputDirectory + "/uchime-dn-out2\
  --p-minh 0.20"

subprocess.run(Command2,shell=True,check=True)




## rename 

Command = "qiime  feature-table summarize \
  --i-table  " + inputDirectory + "/table.qza \
  --o-visualization  " + inputDirectory + "/table.qzv \
  --m-sample-metadata-file  " + metadata + "/sample-metadata.tsv && \
  qiime feature-table tabulate-seqs \
  --i-data  " + inputDirectory + "/rep-seqs.qza \
  --o-visualization  " + inputDirectory + "/rep-seqs.qzv" 

## extract the rep-seqs and table,combine them together

Command = Command + " && " + "qiime tools export \
      --input-path  " + inputDirectory + "/table.qza  \
      --output-path  " + inputDirectory + "/Feature-table-result"


Command = Command + " && " + "qiime tools export \
      --input-path  " + inputDirectory + "/table.qzv  \
      --output-path  " + inputDirectory + "/Feature-table-result"


Command = Command + " && " + "qiime tools export \
      --input-path  " + inputDirectory + "/rep-seqs.qza \
      --output-path  " + inputDirectory + "/Feature-table-result \
      " + " && " + "biom convert \
      -i " + inputDirectory + "/Feature-table-result/feature-table.biom \
      -o  " + inputDirectory + "/Feature-table-result/feature-table.tsv --to-tsv" #feature table

#os.system("cat ../results/2.Feature_table/Feature-table-result/feature-table.tsv | sed 's/#OTU ID/Feature.ID/g' > ../results/2.Feature_table/Feature-table-result/feature-table2.tsv ")

#subprocess.run(Command,shell=True,check=True)

#==================== combine the fasta with ASVs

Command1 = "python3 b0.z.fasta2tsv.py"

#subprocess.run(Command1,shell=True,check=True)