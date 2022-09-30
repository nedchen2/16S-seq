# In last script, we use some of the data from QIIME to do species construction analysis

# export the table
# qiime tools export   --input-path table.qza   --output-path exported-feature-table



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
                    default="../results/4.Diversity_ana")
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


# run 
#============
#===========may be add some taxa collapse here 
print ("=============Start Taxa Collapsing and Differential Abundance Analysis=============")

Command1 = "qiime taxa collapse \
      --i-table " + inputDirectory + "/feature-table-correct.qza\
      --i-taxonomy " + outputDirectory + "/taxonomy-corrected.qza\
      --p-level 6\
      --o-collapsed-table " + outputDirectory + "/table-l6.qza"

Command1 = Command1 + " && " + "qiime tools export \
      --input-path  " + outputDirectory + "/table-l6.qza \
      --output-path  " + outputDirectory + "/Feature-table-result \
      " + " && " + "biom convert \
      -i " + outputDirectory + "/Feature-table-result/feature-table.biom \
      -o  " + outputDirectory + "/Feature-table-result/feature-table-taxa.tsv --to-tsv \
       "

Command1 = Command1 + " && " + "qiime composition add-pseudocount\
    --i-table " + outputDirectory + "/table-l6.qza \
    --o-composition-table " + outputDirectory + "/comp-gut-table-l6.qza"

Metadatalist = ["Enterotype"]

for i in Metadatalist:
  Command1 = Command1 + " && " + "qiime composition ancom \
    --i-table " + outputDirectory + "/comp-gut-table-l6.qza \
    --m-metadata-file " + classifier + "/sample-metadata2.tsv \
    --m-metadata-column " + i + " \
    --o-visualization " + outputDirectory + "/Feature-table-result/l6-ancom-" + i +".qzv"

#subprocess.run(Command1,shell=True,check=True)

# ================================== esvS
Command4 = "qiime composition add-pseudocount\
    --i-table " + outputDirectory + "/feature-table-correct.qza \
    --o-composition-table " + outputDirectory + "/comp-gut-table-ESVs.qza"

Metadatalist = ["Species","CollectionSite","Infection"]

for i in Metadatalist:
  Command4 = Command4 + " && " + "qiime composition ancom \
    --i-table " + outputDirectory + "/comp-gut-table-ESVs.qza \
    --m-metadata-file " + classifier + "/sample-metadata2.tsv \
    --m-metadata-column " + i + " \
    --o-visualization " + outputDirectory + "/Feature-table-result/ESVs-ancom-" + i +".qzv"

#subprocess.run(Command4,shell=True,check=True)




#according to the information so far, there are no difference taxa according to the current data.
#only difference exists between species



