# Auto convert id to taxonomy

from Bio import Entrez
import glob 
import os
import pandas as pd
Entrez.email = "nedchen2@126.com" 
file_list = glob.glob("../results/blast/*.csv")

df = pd.DataFrame()
for i in range(0,len(file_list)):
    Featureid = os.path.basename(file_list[i]).split(".")[0]
    df2 = pd.read_csv(file_list[i])
    string = df2.loc[0,"Accession  "]
    AccessionNum = string.split(",")[1].replace("\"","").replace(")","") 
    print ("*******ID -> Taxonomy*********")
    print (i)
    print (string)
    print (AccessionNum)
    data = Entrez.efetch(db = "nucleotide", id = AccessionNum,rettype = "gb",retmode="xml" )
    record = Entrez.read(data)  
    Taxonomy = record[0]["GBSeq_taxonomy"]
    df.loc[i,"Row.names"] = Featureid 
    df.loc[i,"AccessionNum"] = AccessionNum
    df.loc[i,"Taxonomy"] = Taxonomy
    df.loc[i,"Organism"] = record[0]["GBSeq_organism"] 

df.to_csv("../results/blast/Manual_blast_result.tsv",sep="\t")
