import os
import argparse

parser = argparse.ArgumentParser(description="format_transform.")
parser.add_argument("-fa", "--fasta", help="Input .fa file [default: ./fasta]",
                    default="../results/2.Feature_table/Feature-table-result/dna-sequences.fasta")
parser.add_argument("-r", "--resultfile", help="Output file with selected gene [default:./Seq_selected.txt].",
                    default="../results/2.Feature_table/Feature-table-result/dna-sequence.tsv")
args = parser.parse_args()

fasta = os.path.abspath(args.fasta)
outputDirectory = os.path.abspath(args.resultfile)

#check file
if os.path.exists(outputDirectory):
    os.remove(outputDirectory)

#create file
fa = open(fasta, 'r')
fr = open(outputDirectory, 'a')

dict={}

fr.write("Feature.ID" + "\t" + "Sequence" + "\n") # add the title

for line in fa:
    if line.startswith('>'):
        line = line.strip('\n')
        line = line.split(">", 1)
        name = str(line[1])
        dict[name] = ''
    else:
        dict[name] += line.replace('\n','')

for ID in dict.keys():
        fr.write(ID + "\t" + dict[ID] + "\n")
fr.close()
fa.close()

# ========= merge the sequence with the feature table
import pandas as pd

df = pd.read_table("../results/2.Feature_table/Feature-table-result/feature-table.tsv",skiprows = 0,header = 1)
df = df.rename(columns={"#OTU ID":"Feature.ID"}).set_index("Feature.ID")

df2 = df2 = pd.read_table("../results/2.Feature_table/Feature-table-result/dna-sequence.tsv",index_col = "Feature.ID")

df3 = df.merge(df2, on="Feature.ID")

df3.to_csv("../results/2.Feature_table/Feature-table-result/feature-table-with-seq.tsv",sep = "\t")


