# Author : Congjia Chen

# ========= import the module ======

import os
import glob
import sys            #mainly leave the running
import pandas as pd
import argparse
import configparser
import pickle

# ==================================
# read external argument
parser = argparse.ArgumentParser(description="Script for confirm the number of raw_data and output sample id")
parser.add_argument("-c", "--configfile", help="The config file.(default:./config.ini)",
                    default="./config.ini") #set the software location in this file
parser.add_argument("-o", "--outputDirectory", help="Output directory with results",
                    default="../results")
parser.add_argument("-i", "--inputDirectory", help="input directory with fastq files",
                    default="../16sRaw")

# Might add something related to HPC

args = parser.parse_args()

# Process the command line arguments.
outputDirectory = os.path.abspath(args.outputDirectory)
configfile = os.path.abspath(args.configfile)
inputDirectory = os.path.abspath(args.inputDirectory)

# get the software list by config file
config = configparser.ConfigParser()
config.read(configfile, encoding="utf-8")
python3 = config.get("software", "python3")
R = config.get("software", "R")
Trimmomatic = config.get("software", "Trimmomatic")

# check the directory
if os.path.exists(inputDirectory):
  print ("Input dir exists")
else :
  os.mkdir(inputDirectory)
if os.path.exists(outputDirectory):
  print ("Output dir exists")
else :
  os.mkdir(outputDirectory)

# get the file list by glob

raw_data_list = glob.glob(inputDirectory + "/*.fastq.gz")

Reads_name_list = [os.path.basename(i).split(sep=".")[0] for i in raw_data_list]

Sample_name = [i.split(sep ="_")[-1] for i in Reads_name_list]

# check the Sample number

if len(set(Sample_name)) == 0 :
    print ("[ERROR] No Sample in given input Directory")
    sys.exit(0)
elif len(Reads_name_list) % 2 != 0 or len(Reads_name_list)/2 != len(set(Sample_name)) : # check if the sample is paired
    print ("[ERROR] Please check the sample with the paired reads")
    sys.exit(0)
else:
    print ("We are PROCESSING " + str(len(set(Sample_name))) + " Samples")

# Output the dataframe of Sample list

sample_frame = pd.DataFrame({"Reads":Reads_name_list, "Sample": Sample_name})

with open('../results/sample_frame.pickle', 'wb') as f:
    print ("=====================storing the result=====================")
    pickle.dump(sample_frame, f)
    print ("The result has been stored in the ../results")