
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
                    default="../results/4.Diversity_ana/")
parser.add_argument("-i", "--inputDirectory", help="input directory with fastq files",
                    default="../results/2.Feature_table/")
parser.add_argument("-t", "--threads", help="Threads",
                    default="2")
parser.add_argument("-m", "--metadata", help="metadata",
                    default="./")   
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



# re import the data into qiime
# add "# Constructed from biom file" to the first column first row

#Command0 = "Rscript b2.g.reimport_and_filter.R"

Command0 = "pwd"

Command0 = Command0 + " && biom convert -i ../results/4.Diversity_ana/corrected_abundance_table.tsv \
  -o ../results/4.Diversity_ana/table-correct.biom\
  --table-type=\"OTU table\" --to-hdf5"

Command0 = Command0 + " && qiime tools import \
    --input-path ../results/4.Diversity_ana/table-correct.biom  --type 'FeatureTable[Frequency]' \
    --input-format BIOMV210Format  \
    --output-path ../results/4.Diversity_ana/feature-table.qza"

Command0 = Command0 + " && qiime tools import \
  --input-path ../results/4.Diversity_ana/taxonomy-corrected.tsv \
  --type 'FeatureData[Taxonomy]' \
  --output-path ../results/4.Diversity_ana/taxonomy-corrected.qza"

# filter feature table , remove the sequence error

Command0 = Command0 + " && qiime feature-table filter-features \
  --i-table ../results/4.Diversity_ana/feature-table.qza \
  --p-min-frequency 10 \
  --o-filtered-table ../results/4.Diversity_ana/feature-frequency-filtered-table.qza"

# 50% prevalence

Command0 = Command0 + " && qiime feature-table filter-features \
  --i-table ../results/4.Diversity_ana/feature-frequency-filtered-table.qza \
  --p-min-samples 45 \
  --o-filtered-table ../results/4.Diversity_ana/feature-table-correct.qza"

"""
qiime feature-table filter-samples \
  --i-table ../results/4.Diversity_ana/feature-table.qza \
  --p-min-frequency 1000 \
  --o-filtered-table ../results/4.Diversity_ana/sample-frequency-filtered-table.qza
"""


# filter sequence
Command0 =  Command0 + " && qiime feature-table filter-seqs \
  --i-data  ../results/2.Feature_table/rep-seqs.qza \
  --i-table ../results/4.Diversity_ana/feature-table-correct.qza\
  --o-filtered-data ../results/4.Diversity_ana/rep-seqs-correct.qza"


Command0 = Command0 + " && qiime  feature-table summarize \
  --i-table  ../results/4.Diversity_ana/feature-table-correct.qza \
  --o-visualization  ../results/4.Diversity_ana/feature-table-correct.qzv \
  --m-sample-metadata-file  ./sample-metadata.tsv"

Command0 = Command0 + " && qiime diversity alpha-rarefaction \
  --i-table ../results/4.Diversity_ana/feature-table-correct.qza \
  --p-max-depth 10000 \
  --p-min-depth 10 \
  --m-metadata-file ./sample-metadata.tsv \
  --o-visualization ../results/4.Diversity_ana/alpha-rarefaction.qzv"


Command = "qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences ../results/4.Diversity_ana/rep-seqs-correct.qza \
  --o-alignment  " + outputDirectory + "/aligned-rep-seqs.qza \
  --o-masked-alignment  " + outputDirectory + "/masked-aligned-rep-seqs.qza \
  --o-tree  " + outputDirectory + "/unrooted-tree.qza \
  --o-rooted-tree  " + outputDirectory + "/rooted-tree.qza"

# 5000
# decided from feature-table-correct.qzv

Command = Command + " && " + "\
qiime diversity core-metrics-phylogenetic \
  --i-phylogeny " + outputDirectory + "/rooted-tree.qza \
  --i-table " +  outputDirectory + "/feature-table-correct.qza \
  --p-sampling-depth 1010\
  --m-metadata-file " + metadata + "/sample-metadata.tsv \
  --output-dir " + outputDirectory + "/core-metrics-results"
# this core matrix give out certain outpt
# No difference in alpha diversity

#3355 

#


#core-metrics-results/unweighted_unifrac_emperor.qzv
#core-metrics-results/jaccard_emperor.qzv
#core-metrics-results/bray_curtis_emperor.qzv
#core-metrics-results/weighted_unifrac_emperor.qzv




#=====================Alternative==========================

#Alpha_index_list = ["observed_features","ace","simpson","shannon"]

#Command = " echo '=============Start Alpha Diversity Analysis=============' "

#for i in Alpha_index_list:
#  Command = Command + " && qiime diversity alpha \
#  --i-table  "  + inputDirectory + "/table.qza \
#  --p-metric " + i + " \
#  --o-alpha-diversity " + outputDirectory + "/" + i + "_vector.qza"
#  Command = Command + "&& qiime tools export \
#      --input-path  " + outputDirectory + "/" + i + "_vector.qza \
#      --output-path  " + outputDirectory + "/alpha_diversity_vector/" + i  

#subprocess.run(Command,shell=True,check=True)

#=================================================


# ================== Beta Diversity ===========

# stat and visualize
Command1 =  "qiime diversity alpha-group-significance \
  --i-alpha-diversity " + outputDirectory + "/core-metrics-results/faith_pd_vector.qza \
  --m-metadata-file sample-metadata.tsv \
  --o-visualization " + outputDirectory + "/core-metrics-results/faith-pd-group-significance.qzv"

Command1 = Command1 + " && " + "qiime diversity alpha-group-significance \
  --i-alpha-diversity " + outputDirectory + "/core-metrics-results/shannon_vector.qza \
  --m-metadata-file sample-metadata.tsv \
  --o-visualization " + outputDirectory + "/core-metrics-results/shannon-group-significance.qzv"

Command1 = Command1 + " && " + "qiime diversity alpha-group-significance \
  --i-alpha-diversity " + outputDirectory + "/core-metrics-results/evenness_vector.qza \
  --m-metadata-file " + metadata + "/sample-metadata.tsv \
  --o-visualization " + outputDirectory + "/core-metrics-results/evenness-group-significance.qzv"

Command1 = Command1 + " && " + "qiime diversity alpha-group-significance \
  --i-alpha-diversity " + outputDirectory + "/core-metrics-results/observed_features_vector.qza \
  --m-metadata-file sample-metadata.tsv \
  --o-visualization " + outputDirectory + "/core-metrics-results/observed_features_vector-group-significance.qzv"







Command2 = "cd ../results/4.Diversity_ana" 

list_of_metadata = ["CollectionSite","Species","Interaction","gut_parasite_richness"]
list_of_metadata = ["Crithidia_binomial","Nosema_binomial"]
matrix_name = ["weighted_unifrac_distance_matrix","unweighted_unifrac_distance_matrix","bray_curtis_distance_matrix","jaccard_distance_matrix"]

for j in matrix_name:
  for i in list_of_metadata:
    Command2 = Command2 + " && qiime diversity beta-group-significance \
    --i-distance-matrix core-metrics-results/" + j + ".qza \
    --m-metadata-file ../../code/sample-metadata.tsv \
    --m-metadata-column " + i + " \
    --o-visualization core-metrics-results/" + j + "-" + i + "-significance.qzv \
    --p-pairwise "



"""
  qiime diversity beta-group-significance \
  metadata-file ../../code/sample-metadata.tsv \
  --m-metadata-column CollectionSite \
  --o-visualization core-metrics-results/weighted-unifrac-CollectionSite-significance.qzv \
  --p-pairwise--i-distance-matrix core-metrics-results/weighted_unifrac_distance_matrix.qza \
  --m-

qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results/weighted_unifrac_distance_matrix.qza \
  --m-metadata-file ../../code/sample-metadata.tsv \
  --m-metadata-column Species \
  --o-visualization core-metrics-results/weighted-unifrac-Species-significance.qzv \
  --p-pairwise

"""



#subprocess.run(Command0,shell=True,check=True)
#subprocess.run(Command,shell=True,check=True)
#subprocess.run(Command1,shell=True,check=True)
subprocess.run(Command2,shell=True,check=True)
