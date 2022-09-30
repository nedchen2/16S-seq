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
                    default="../results/3.Taxonomy_ana")
parser.add_argument("-i", "--inputDirectory", help="input directory with fastq files",
                    default="../results/2.Feature_table/")
parser.add_argument("-m", "--metadata", help="metadata",
                    default="./")     
parser.add_argument("-t", "--threads", help="Threads",
                    default="3")
parser.add_argument("-e", "--error", help="error rate",
                    default="1") # admit one error rate
parser.add_argument("-u", "--uncultured", help="Exclude Uncultured or not",
                    default="Exclude") # admit one error rate                    

# Might add something related to HPC

args = parser.parse_args()

# Process the command line arguments.
outputDirectory = os.path.abspath(args.outputDirectory)
configfile = os.path.abspath(args.configfile)
inputDirectory = os.path.abspath(args.inputDirectory)
threads = str(args.threads)
error =  str(args.error)
classifier =  os.path.abspath(args.metadata) # dir of the pre-trained classifier
uncultured = str(args.uncultured)
 
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

# Region of 16S sequence
PRIMER1 = "AGAGTTTGATCCTGGCTCAG" #8F
PRIMER2 = "GCTGCCTCCCGTAGGAGT"   #338R


# Later could add some arguments (Outputdirectory etc)
# add some function to check and download the sequence and taxanomy

# silva data base sequence
# source and license
# https://docs.qiime2.org/2022.2/data-resources/
# Output the sequence into code

print ("=============Start Training the Classifier=============")

Command0 = "qiime feature-classifier extract-reads \
--i-sequences " + classifier + "/silva-138-99-seqs.qza \
--p-f-primer " +  PRIMER1 + " \
--p-r-primer " + PRIMER2 + " \
--o-reads "+ classifier + "/silva-138-99-8F338R.qza \
--verbose"

# unidentified
# gut_metagenome
# human_gut
# soil_bacterium
if uncultured != "Exclude":
  Command0 = Command0 + " && qiime feature-classifier fit-classifier-naive-bayes \
    --i-reference-reads "+ classifier + "/silva-138-99-8F338R.qza \
    --i-reference-taxonomy "+ classifier + "/silva-138-99-tax.qza \
    --o-classifier "+ classifier + "/silva-138-99-8F338R-classifier.qza"  # Train the classifer # Output the classifer to code 
  classifier_name = "silva-138-99-8F338R-classifier.qza"
else:
  Command0 = Command0 + " && qiime rescript filter-taxa \
    --i-taxonomy ./silva-138-99-tax.qza \
    --p-exclude \"g__uncultured\" \
    --o-filtered-taxonomy ./silva-138-99-tax-ExUncultured.qza"# or we exclude the uncultured species  # exclude the taxa
  Command0 = Command0 + " && qiime taxa filter-seqs --i-sequences silva-138-99-8F338R.qza \
    --i-taxonomy silva-138-99-tax.qza \
       --p-exclude  \"g__uncultured\" \
         --o-filtered-sequences ./silva-138-99-8F338R-ExUncultured.qza"   # exclude the sequence
  Command0 =  Command0 + " && qiime feature-classifier fit-classifier-naive-bayes \
    --i-reference-reads "+ classifier + "/silva-138-99-8F338R-ExUncultured.qza \
    --i-reference-taxonomy "+ classifier + "/silva-138-99-tax-ExUncultured.qza \
    --o-classifier "+ classifier + "/silva-138-99-8F338R-ExUncultured-classifier.qza"  # Train the uncultured classifer
  classifier_name = "silva-138-99-8F338R-ExUncultured-classifier.qza" #uncultured classifier

subprocess.run(Command0,shell=True,check=True)
# qiime feature-classifier fit-classifier-naive-bayes --i-reference-reads ./silva-139-99-8F338R-ExUncultured.qza --i-reference-taxonomy ./silva-138-99-tax-ExUncultured.qza --o-classifier silva-138-99-8F338R-ExUncultured-classifier.qza


## TAXONOMY analysis
# the region we amplified is 8F-338R--- V1,V2 region
# Given the sequence, I want to know the taxonomy  information
# Use silva pre-trained  database q2-feature-classifier 
# silva may better than greengene
# however, the trained database is much bigger than GG, which will cause memory inadequate,
# therefore,here we use gg full-length data

# We tried the pre-trained GG Database,silva database classifer

# We will use HPC to do silva classification

# We will try classify the taxanomy with self-trained classifier which is recommanded by QIIME
print ("=============Start Classification=============")


Command ="qiime feature-classifier classify-sklearn \
  --i-classifier " + classifier + "/" + classifier_name + "\
  --i-reads " + inputDirectory + "/rep-seqs.qza \
  --o-classification " + outputDirectory + "/taxonomy.qza\
  --p-n-jobs " + threads

Command = Command + " && " + "qiime tools export \
      --input-path " + outputDirectory + "/taxonomy.qza \
      --output-path " + outputDirectory + "/Taxonomy_export\
      "
# visualize the barplot of taxonomy distribution
# ================ Check the blank sample here

Command = Command + " && " + "qiime taxa barplot \
   --i-table " + inputDirectory + "/table.qza \
   --i-taxonomy " + outputDirectory + "/taxonomy.qza \
   --m-metadata-file " + classifier + "/sample-metadata.tsv \
   --o-visualization " + outputDirectory + "/taxa-bar-plots.qzv"
# taxa-bar-plots.qzv: the counts of taxa in different level
# could use the data in it for example the level-2.csv to visualize

Command = Command + " && " + "qiime tools export \
      --input-path " + outputDirectory + "/taxa-bar-plots.qzv \
      --output-path " + outputDirectory + "/Taxonomy_export \
      "

subprocess.run(Command,shell=True,check=True)



# Do blast consensus analysis

Command2 = "python3 b1.w.Feature_taxa_Doblast.py"

#subprocess.run(Command2,shell=True,check=True)
# ================= remove the possible contaminant


