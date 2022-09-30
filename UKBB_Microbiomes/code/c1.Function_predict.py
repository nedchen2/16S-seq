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
                    default="../results/8.DEG")
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


print ("=============Start Functional Profiling=============")

Command2 = "qiime picrust2 full-pipeline \
   --i-table " + inputDirectory + "/core-metrics-results/rarefied_table.qza \
   --i-seq " + inputDirectory + "/rep-seqs-correct.qza \
   --output-dir " + outputDirectory + "/q2-picrust2_output \
   --p-placement-tool sepp \
   --p-threads 2 \
   --p-hsp-method pic \
   --p-max-nsti 1\
   --p-min-align 0.97\
   --verbose"


Command2 = Command2 + " && " +  "qiime tools export \
      --input-path  " + outputDirectory + "/q2-picrust2_output/pathway_abundance.qza \
      --output-path  " + outputDirectory + "/q2-picrust2_output/Functional_profiling_pathway \
      " + " && " + "biom convert \
      -i " + outputDirectory + "/q2-picrust2_output/Functional_profiling_pathway/feature-table.biom \
      -o  " + outputDirectory + "/q2-picrust2_output/Functional_profiling_pathway/Pathway_abundance.tsv --to-tsv \
       "

Command2 = Command2 + " && " + "qiime tools export \
      --input-path  " + outputDirectory + "/q2-picrust2_output/ko_metagenome.qza \
      --output-path  " + outputDirectory + "/q2-picrust2_output/Functional_profiling_ko \
      " + " && " + "biom convert \
      -i " + outputDirectory + "/q2-picrust2_output/Functional_profiling_ko/feature-table.biom \
      -o  " + outputDirectory + "/q2-picrust2_output/Functional_profiling_ko/ko_metagenome.tsv --to-tsv \
       "

Command2 = Command2 + " && " + "qiime tools export \
      --input-path  " + outputDirectory + "/q2-picrust2_output/ec_metagenome.qza \
      --output-path  " + outputDirectory + "/q2-picrust2_output/Functional_profiling_ec \
      " + " && " + "biom convert \
      -i " + outputDirectory + "/q2-picrust2_output/Functional_profiling_ec/feature-table.biom \
      -o  " + outputDirectory + "/q2-picrust2_output/Functional_profiling_ec/ec_metagenome.tsv --to-tsv \
       "

subprocess.run(Command2,shell=True,check=True)



# use PICRUST2 TO ADD THE DISCRIPTION
# =============== Enrichment analysis

# conda activate picrust2
# Command = "add_descriptions.py -i ../results/3.Taxonomy_ana/q2-picrust2_output/Functional_profiling_pathway/Pathway_abundance.tsv -o ../results/3.Taxonomy_ana/q2-picrust2_output/Functional_profiling_pathway/pred_metagenome_unstrat_descrip.tsv -m METACYC"
# add_descriptions.py -i ../results/3.Taxonomy_ana/q2-picrust2_output/Functional_profiling_ec/ec_metagenome.tsv -o ../results/3.Taxonomy_ana/q2-picrust2_output/Functional_profiling_ec/pred_metagenome_unstrat_descrip.tsv -m EC
# add_descriptions.py -i ../results/3.Taxonomy_ana/q2-picrust2_output/Functional_profiling_ko/ko_metagenome.tsv -o ../results/3.Taxonomy_ana/q2-picrust2_output/Functional_profiling_ko/pred_metagenome_unstrat_descrip.tsv -m KO

# add_descriptions.py -i ../results/8.DEG/q2-picrust2_output/Functional_profiling_ko/ko_metagenome.tsv -o ../results/8.DEG/q2-picrust2_output/Functional_profiling_ko/pred_metagenome_unstrat_descrip.tsv -m KO


# ============== differential Analysis of KO and Metacyc 
Metadatalist = ["Enterotype","Crithidia_binomial","Nosema_binomial"]

Command3 = "echo '=============Start Differential Functional Analysis=============' "

Command3 = Command3 + " && " + "qiime composition add-pseudocount\
    --i-table " + outputDirectory + "/q2-picrust2_output/ko_metagenome.qza \
    --o-composition-table " + outputDirectory + "/q2-picrust2_output/comp-gut-table-ko.qza"

for i in Metadatalist:
  Command3 = Command3  + " && qiime composition ancom \
    --i-table " + outputDirectory + "/q2-picrust2_output/comp-gut-table-ko.qza\
    --m-metadata-file " + classifier + "/sample-metadata2.tsv \
    --m-metadata-column " + i + " \
    --o-visualization " + outputDirectory + "/q2-picrust2_output/ko-ancom-" + i +".qzv"

#subprocess.run(Command3,shell=True,check=True)
