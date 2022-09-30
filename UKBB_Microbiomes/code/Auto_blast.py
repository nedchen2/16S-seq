from Bio.Blast import NCBIWWW
import pandas as pd
from Bio.Blast import NCBIXML

df = pd.read_csv("../results/4.Diversity_ana/table_contaminant.csv")
n=1
Topn = df.iloc[[42,56,58,71,97,106,114,173,169,145,146]].reset_index()

fasta_string = Topn["Sequence"]

E_VALUE_THRESH = 0.04

for i in range(0,len(fasta_string)):
    result_handle = NCBIWWW.qblast("blastn", "nt", fasta_string[i],megablast= "TRUE",entrez_query = "all [filter] NOT(environmental samples[organism] OR metagenomes[orgn])", hitlist_size=1)
    blast_record = NCBIXML.read(result_handle)
    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            if hsp.expect < E_VALUE_THRESH:
                Topn.loc[i,"title"] = alignment.title
                Topn.loc[i,"length"] = alignment.length
                Topn.loc[i,"evalue"] = hsp.expect
                Topn.loc[i,"MimicPercIdentity"] = hsp.identities/hsp.align_length
                print ('********Alignment*********')
                print ('sequence:', alignment.title)
                print ('length:', alignment.length)
                print ('e value:', hsp.expect)
                print ('MimicPercIdentity', hsp.identities/hsp.align_length)
                print (hsp.query[0:75] + '...')
                print (hsp.match[0:75] + '...')
                print (hsp.sbjct[0:75] + '...')

filename = "../results/4.Diversity_ana/top" + str(n) + "_unassigned_blast_result.tsv"

Topn.to_csv(filename,sep = "\t")


