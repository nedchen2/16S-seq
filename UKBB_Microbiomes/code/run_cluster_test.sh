#!/bin/bash
#PBS -l walltime=8:00:00
#PBS -l select=1:ncpus=8:mem=16gb

module load anaconda3/personal

#activate the env
source activate qiime2-2021.11

echo " QIIME2 is about to run"

# this is where my submit directory is
cd $PBS_O_WORKDIR

#pipeline for deblur-PE

#python3 $HOME/UKBB_Microbiomes/code/a3.raw_data_demuTrim.py
#python3 $HOME/UKBB_Microbiomes/code/a4.raw_data_2_Feature_deblur.py
#python3 $HOME/UKBB_Microbiomes/code/b0.Feature_table_check.py
python3 $HOME/UKBB_Microbiomes/code/b1.Feature_taxa_classification_back.py 

echo "Mission complete" 