# do blast
import os

os.system("mkdir ../results/blast")
os.system("mkdir ../results/blast/0.99-1hit/")
os.system("mkdir ../results/blast/0.97-5hit/")


Command = "cd ../results/blast"
#download NCBI
Command = Command + " &&  qiime rescript get-ncbi-data \
    --p-query '33175[BioProject] OR 33317[BioProject]' \
    --o-sequences ncbi-refseqs-unfiltered.qza \
    --o-taxonomy ncbi-refseqs-taxonomy-unfiltered.qza"
#filter length
Command = Command + " && qiime rescript filter-seqs-length-by-taxon \
    --i-sequences ncbi-refseqs-unfiltered.qza \
    --i-taxonomy ncbi-refseqs-taxonomy-unfiltered.qza \
    --p-labels Archaea Bacteria \
    --p-min-lens 900 1200 \
    --o-filtered-seqs ncbi-refseqs.qza \
    --o-discarded-seqs ncbi-refseqs-tooshort.qza"

Command = Command + " && qiime rescript filter-taxa \
    --i-taxonomy ncbi-refseqs-taxonomy-unfiltered.qza \
    --m-ids-to-keep-file ncbi-refseqs.qza \
    --o-filtered-taxonomy ncbi-refseqs-taxonomy.qza"

#filter uncultured
Command = Command + " && qiime taxa filter-seqs \
    --i-sequences ncbi-refseqs.qza \
    --i-taxonomy ncbi-refseqs-taxonomy.qza \
    --p-exclude  \"uncultured\",\"Uncultured\"  \
    --o-filtered-sequences ncbi-refseqs-filtered.qza"

Command = Command + " && qiime rescript filter-taxa \
    --i-taxonomy ncbi-refseqs-taxonomy.qza \
    --m-ids-to-keep-file  ncbi-refseqs-filtered.qza\
    --o-filtered-taxonomy ncbi-refseqs-taxonomy-filtered.qza"

#os.system(Command)



#do blast

Command2 = "cd ../results/blast && qiime feature-classifier classify-consensus-blast \
    --i-query ../2.Feature_table/rep-seqs.qza\
    --i-reference-reads ncbi-refseqs-filtered.qza\
    --i-reference-taxonomy ncbi-refseqs-taxonomy-filtered.qza\
    --p-perc-identity 0.99\
    --p-query-cov 0.90\
    --p-maxaccepts 1\
    --o-classification ./0.99-1hit/taxonomy.qza"


Command2 = Command2 + " && qiime feature-classifier classify-consensus-blast \
    --i-query ../2.Feature_table/rep-seqs.qza\
    --i-reference-reads ncbi-refseqs-filtered.qza\
    --i-reference-taxonomy ncbi-refseqs-taxonomy-filtered.qza\
    --p-perc-identity 0.97\
    --p-query-cov 0.90\
    --p-maxaccepts 5\
    --o-classification ./0.97-5hit/taxonomy.qza"  

Command2 = Command2 + " && qiime tools export \
     --input-path  ./0.99-1hit/taxonomy.qza \
     --output-path ./0.99-1hit/" 

Command2 =  Command2 + " && qiime tools export \
     --input-path  ./0.97-5hit/taxonomy.qza \
     --output-path ./0.97-5hit/"

os.system(Command2)

#vsearch 
"""
use vsearch to look at the assignment
os.system("qiime feature-classifier classify-consensus-vsearch \
    --i-query ../2.Feature_table/rep-seqs.qza  \
    --i-reference-reads ncbi-refseqs-filtered.qza \
    --i-reference-taxonomy ncbi-refseqs-taxonomy-filtered.qza \
    --p-perc-identity 0.90 \
    --p-query-cov 0.3 \
    --p-top-hits-only \
    --p-maxaccepts 10 \
    --p-strand 'both' \
    --p-unassignable-label 'Unassigned' \
    --p-threads 2 \
    --o-classification taxonomy.qza")
"""    