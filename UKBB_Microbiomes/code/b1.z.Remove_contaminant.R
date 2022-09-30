require(tidyverse)

# ======== Summary of the OTU table including confidence, taxa annotation and sequence,and total frequency

taxonomy <- read.delim ("../results/3.Taxonomy_ana/Taxonomy_export/taxonomy.tsv")


#feature_table <- read.table("../results/2.Feature_table/Feature-table-result/feature-table2.tsv",header=T,check.names = F)
#feature_table_seq <-  read.table("../results/2.Feature_table/Feature-table-result/feature-table-with-seq.tsv",header=T,check.names = F) 

NUMBER <- ncol(read.table("../results/2.Feature_table/Feature-table-result/feature-table-with-seq.tsv",header=T,check.names = F))
test = read.table("../results/2.Feature_table/Feature-table-result/feature-table-with-seq.tsv",header=T,check.names = F) %>%
  column_to_rownames(var = "Feature.ID")

feature_table_with_seq <- read.table("../results/2.Feature_table/Feature-table-result/feature-table-with-seq.tsv",header=T,check.names = F)  %>%
  column_to_rownames(var = "Feature.ID") %>%
  mutate(Frequency = rowSums(.[1:(NUMBER-2)])) %>%  # 2 is the Non-sample column number
  arrange(desc(Frequency)) %>% rownames_to_column(var = "Feature.ID")


Join_Three_table <-left_join(feature_table_with_seq,taxonomy,by = "Feature.ID") %>% arrange(desc(Frequency))

Join_Three_table[,"Suspectable"] = ""

# first load the possible contaminant genus to the R

list_of_contaminate <- read.csv("../code/List_of_contaminent.tsv",sep = "\t")

vector_of_contaminate <- str_remove(unlist(str_split(list_of_contaminate$List.of.constituent.contaminant.genera,","))," " )


for (i in seq(length(vector_of_contaminate))){
  index <- grepl(paste0("g__",vector_of_contaminate[i]),
                 Join_Three_table$Taxon,
                 ignore.case = F)
  Join_Three_table[index,"Suspectable"] = paste0(Join_Three_table[index,"Suspectable"],"; ",vector_of_contaminate[i]) 
}

# See what kind of contaminant we have 
unique(Join_Three_table[Join_Three_table$Suspectable != "","Suspectable"])
Join_Three_table[seq(1:10),"Taxon"]

write.table(Join_Three_table,file="../results/3.Taxonomy_ana/Feature_abundance_table_with_seq_and_taxa_annotation.tsv",row.names = F)


# ===== deal with unassigned and add URL to the fasta file

idx1 = !grepl("g__",Join_Three_table$Taxon,ignore.case = T)
idx2 = grepl("uncultured",Join_Three_table$Taxon,ignore.case = T)

idx = idx1 | idx2

add_blastn <- function(data){
  blasteURL = paste0("http://www.ncbi.nlm.nih.gov/BLAST/Blast.cgi?ALIGNMENT_VIEW=Pairwise&PROGRAM=blastn&DATABASE=nt&CMD=Put&QUERY=",data$Sequence)
  return(blasteURL)
}

print (paste0("Ratio of unassigned:",sum(idx)/nrow(Join_Three_table)*100,"%"))
unassigned = Join_Three_table[idx,] %>% mutate(blastURL = add_blastn(.))

# output for autoblast
write.table(unassigned,file="../results/3.Taxonomy_ana/Unassigned_blastn-result.tsv",row.names = F)



# =============== test for remove the contaminate
df_test <- feature_table_with_seq %>% dplyr::select(!c("Sequence","Frequency")) %>% column_to_rownames(var = "Feature.ID")

# Here, we assume sample 2 as the blank sample
# subtract the feature from blank sample

df_subtraction = df_test - df_test$Blank 

df_subtraction[df_subtraction < 0 ] = 0 

# read the taxonomy information
df_taxonomy <- Join_Three_table %>% dplyr::select("Feature.ID","Sequence","Taxon","Confidence") %>% column_to_rownames(var = "Feature.ID")

df_join <- merge(df_subtraction,df_taxonomy,by="row.names")

df_join[,"Suspectable"] = ""

# first load the possible contaminant genus to the R

list_of_contaminate <- read.csv("../code/List_of_contaminent.tsv",sep = "\t")

vector_of_contaminate <- str_remove(unlist(str_split(list_of_contaminate$List.of.constituent.contaminant.genera,","))," " )

for (i in seq(length(vector_of_contaminate))){
  index <- grepl(paste0("g__",vector_of_contaminate[i]),
                 df_join$Taxon,
                 ignore.case = F)
  df_join[index,"Suspectable"] = paste0(df_join[index,"Suspectable"],"; ",vector_of_contaminate[i]) 
}

# See what kind of contaminant we have 
unique(df_join[df_join$Suspectable != "","Suspectable"])


df_final = df_join  %>% mutate(Frequency = rowSums(.[2:92])) %>%  # 2 is the Non-sample column number
  arrange(desc(Frequency))

write.table(df_final,file="../results/3.Taxonomy_ana/Feature_abundance_table_without_blank_sample.tsv",row.names = F)

# ===== deal with unassigned and add URL to the fasta file

idx3 = !grepl("g__",df_final$Taxon,ignore.case = T)
idx4 = grepl("uncultured",df_final$Taxon,ignore.case = T)

idx_noblank = idx3 | idx4

add_blastn <- function(data){
  blasteURL = paste0("http://www.ncbi.nlm.nih.gov/BLAST/Blast.cgi?ALIGNMENT_VIEW=Pairwise&PROGRAM=blastn&DATABASE=nt&CMD=Put&QUERY=",data$Sequence)
  return(blasteURL)
}

unassigned_noblank = df_final[idx_noblank,] %>% mutate(blastURL = add_blastn(.))

# output for autoblast
write.table(unassigned_noblank,file="../results/3.Taxonomy_ana/Unassigned_blastn-result-noblank.tsv",row.names = F)


# ============ Combine the consensus blast result with the unassigned table ==============

test_taxa <-  "k__Bacteria; p__Firmicutes; c__Bacilli; o__Lactobacillales; f__Lactobacillaceae; g__Lactobacillus; s__bombicola"

uninomial_to_bionomial <- function(test_taxa){
  species <- str_extract(test_taxa,"(?<=s__)(.*)")
  
  print("Species assigned")
  
  genus <- str_extract(test_taxa,"(?<=g__)(.*)(?=; s__?)")
  
  bionomial_nomen <- paste0(genus,"_",species)
  
  taxonomy_before <- str_extract(test_taxa,"(.*)(?=; s__?)")
  
  new_taxonomy <- paste0(taxonomy_before,"; s__",bionomial_nomen)
  
  return(new_taxonomy)
}
  

Deal_with_GreenGene <- function(test_taxa){
  species <- str_extract(test_taxa,"(?<=s__)(.*)")
  genus   <- str_extract(test_taxa,"(?<=g__)(.*)(?=; s__?)")
  family  <- str_extract(test_taxa,"(?<=f__)(.*)(?=; g__?)")
  order   <- str_extract(test_taxa,"(?<=o__)(.*)(?=; f__?)")
  class   <- str_extract(test_taxa,"(?<=c__)(.*)(?=; o__?)") 
  phylum  <- str_extract(test_taxa,"(?<=p__)(.*)(?=; c__?)") 
  kindom  <- str_extract(test_taxa,"(?<=k__)(.*)(?=; p__?)")
  
  
 if (test_taxa == "Unassigned"){
   return(test_taxa)
  }else if (kindom == "" | is.na(kindom)){
    
    print("Kindom unassigned")
    
    taxonomy_before <- str_extract(test_taxa,"k__.*")
    return(taxonomy_before)
     
  }else if (phylum == "" | is.na(phylum)){
    
    print("phylum unassigned")
    
    taxonomy_before <- str_extract(test_taxa,"(.*)(?=; p__?)")
    return(taxonomy_before)
    
  }else if (class  == ""  | is.na(class)){
    
    print("class unassigned")
    
    taxonomy_before <- str_extract(test_taxa,"(.*)(?=; c__?)")
    return(taxonomy_before)
    
  }else if (order  == ""  | is.na(order)){
    
    print("order unassigned")
    
    taxonomy_before <- str_extract(test_taxa,"(.*)(?=; o__?)")
    return(taxonomy_before)
    
  }else if (family == ""  | is.na(family)){
    
    print("family unassigned")
  
    taxonomy_before <- str_extract(test_taxa,"(.*)(?=; f__?)")
    return(taxonomy_before)
    
  }else if (genus == "" | is.na(genus)) {
    
    print("genus unassigned")
    
    taxonomy_before <- str_extract(test_taxa,"(.*)(?=; g__?)")
    return(taxonomy_before)
    
  }else if (species == ""  | is.na(species)){
    
    print("Species unassigned")
    
    taxonomy_before <- str_extract(test_taxa,"(.*)(?=; s__?)")
    return(taxonomy_before)
  
  }else{
    print("Species assigned")
    
    genus <- str_extract(test_taxa,"(?<=g__)(.*)(?=; s__?)")
    
    bionomial_nomen <- paste0(genus,"_",species)
    
    taxonomy_before <- str_extract(test_taxa,"(.*)(?=; s__?)")
    
    new_taxonomy <- paste0(taxonomy_before,"; s__",bionomial_nomen)
    
    return(new_taxonomy)
  }
}

GreenGene_2_SILVA <- function(df_taxa){
  new_taxonomy <- c()
  for (i in seq(nrow(df_taxa))){
  new_taxonomy  <- c(new_taxonomy,Deal_with_GreenGene(df_taxa[i,"Taxon"]))
  }
  return(new_taxonomy)
}


Combine_unassigned_0.99 <- function(){
  taxonomy_blast <- read.delim ("../results/blast/0.99-1hit/taxonomy.tsv")
  colnames(taxonomy_blast) <- c("Row.names","Taxon_blast","Consensus")
  idx_need_binomial <- str_detect(taxonomy_blast$Taxon_blast,"s__")
  taxonomy_blast[idx_need_binomial,"Taxon_blast"] = uninomial_to_bionomial(taxonomy_blast$Taxon_blast[idx_need_binomial])
  
  
  df_test2 <- left_join(unassigned_noblank,taxonomy_blast,by="Row.names") 
  
  df_final2 <- df_test2 %>% mutate(Taxon_blast_revised =str_replace(str_replace(str_replace(str_replace(Taxon_blast,"k__","d__"),"c__Betaproteobacteria; o__Burkholderiales","c__Gammaproteobacteria; o__Burkholderiales"),"g__Escherichia","g__Escherichia-Shigella"),"g__Shigella","g__Escherichia-Shigella")) %>% 
    mutate(Taxon_blast_revised=str_replace(Taxon_blast_revised,"s__Snodgrassella_alvi wkB2","s__Snodgrassella_alvi"))%>%
    mutate(Taxon_blast_revised=str_replace(Taxon_blast_revised,"s__Phytobacter_massiliensis JC163","s__Phytobacter_massiliensis")) %>%
    mutate(Taxon_blast_revised=str_replace(Taxon_blast_revised,"s__Rhizorhapis_suberifaciens NBRC 15211","s__Rhizorhapis_suberifaciens"))%>%
    dplyr::select("Row.names","Taxon","Confidence","Taxon_blast_revised","Consensus","Frequency")
  
  # index of possibly improved
  idx5 = !grepl("g__",df_test2$Taxon_blast,ignore.case = T)
  idx6 = grepl("uncultured",df_test2$Taxon_blast,ignore.case = T)
  
  idx_blast = !(idx5 | idx6)
  
  # index of whether the result from silva is identical to result the blast
  index_blast2 =c()
  for (i in seq(nrow(df_final2))){
    result <- grepl(df_final2$Taxon[i],df_final2$Taxon_blast_revised[i],ignore.case = T)
    index_blast2 = c(index_blast2,result)
  }
  
  idx_final <- idx_blast & index_blast2
  idx_dissimilar <- idx_blast & !index_blast2
  idx_residue <- !(idx_final | idx_dissimilar)
  
  df_final2[idx_final,"Improve"] = "Improve(0.99-1hit)"
  df_final2[idx_dissimilar,"Improve"] = "Improve but not identical with silva"
  df_final2[idx_residue,"Improve"] = "No improve"
  
  write.table(df_final2,file="../results/3.Taxonomy_ana/Unassigned_blastn-result-noblank-corrected_blast_99.tsv",row.names = F)
  return(df_final2)
}




Combine_unassigned_0.97 <- function(){
  taxonomy_blast <- read.delim ("../results/blast/0.97-5hit/taxonomy.tsv")
  colnames(taxonomy_blast) <- c("Row.names","Taxon_blast","Consensus")
  idx_need_binomial <- str_detect(taxonomy_blast$Taxon_blast,"s__")
  taxonomy_blast[idx_need_binomial,"Taxon_blast"] = uninomial_to_bionomial(taxonomy_blast$Taxon_blast[idx_need_binomial])
  
  df_test2 <- left_join(unassigned_noblank,taxonomy_blast,by="Row.names") 
  
  df_final2 <- df_test2 %>% mutate(Taxon_blast_revised =str_replace(str_replace(str_replace(str_replace(Taxon_blast,"k__","d__"),"c__Betaproteobacteria; o__Burkholderiales","c__Gammaproteobacteria; o__Burkholderiales"),"g__Escherichia","g__Escherichia-Shigella"),"g__Shigella","g__Escherichia-Shigella")) %>% 
    mutate(Taxon_blast_revised=str_replace(Taxon_blast_revised,"s__Snodgrassella_alvi wkB2","s__Snodgrassella_alvi"))%>%
    mutate(Taxon_blast_revised=str_replace(Taxon_blast_revised,"s__Phytobacter_massiliensis JC163","s__Phytobacter_massiliensis")) %>%
    mutate(Taxon_blast_revised=str_replace(Taxon_blast_revised,"s__Rhizorhapis_suberifaciens NBRC 15211","s__Rhizorhapis_suberifaciens"))%>%
    dplyr::select("Row.names","Taxon","Confidence","Taxon_blast_revised","Consensus","Frequency")
  
  # index of possibly improved
  idx5 = !grepl("g__",df_test2$Taxon_blast,ignore.case = T)
  idx6 = grepl("uncultured",df_test2$Taxon_blast,ignore.case = T)
  
  idx_blast = !(idx5 | idx6)
  
  # index of whether the result from silva is identical to result the blast
  index_blast2 =c()
  for (i in seq(nrow(df_final2))){
    result <- grepl(df_final2$Taxon[i],df_final2$Taxon_blast_revised[i],ignore.case = T)
    index_blast2 = c(index_blast2,result)
  }
  
  idx_final <- idx_blast & index_blast2
  idx_dissimilar <- idx_blast & !index_blast2
  idx_residue <- !(idx_final | idx_dissimilar)
  
  df_final2[idx_final,"Improve"] = "Improve(0.97-5hit)"
  df_final2[idx_dissimilar,"Improve"] = "Improve but not identical with silva"
  df_final2[idx_residue,"Improve"] = "No improve"
  
  write.table(df_final2,file="../results/3.Taxonomy_ana/Unassigned_blastn-result-noblank-corrected_blast_97.tsv",row.names = F)
  return(df_final2)
}


unassigned_0.97 <- Combine_unassigned_0.97()
unassigned_0.99 <- Combine_unassigned_0.99()



idx_need_correct <- !str_detect(unassigned_0.97$Taxon_blast_revised,"g__") # 0.97 unassigned
idx_need_correct2 <- str_detect(unassigned_0.99$Taxon_blast_revised,"g__") # 0.99 assigned

idx_need_correct_final <- idx_need_correct & idx_need_correct2 


unassigned_0.97[idx_need_correct_final,"Taxon_blast_revised"] = unassigned_0.99[idx_need_correct_final,"Taxon_blast_revised"]

unassigned_0.97[idx_need_correct_final,"Improve"] = unassigned_0.99[idx_need_correct_final,"Improve"]

unassigned_final_correct <- unassigned_0.97

# =============================================


Ratio_of_unassigned<- function(data){
  idx_new = !grepl("g__",data$Taxon,ignore.case = T)
  idx_new_species = !grepl("s__",data$Taxon,ignore.case = T)
  
  print (paste0("Ratio of unassigned genus:",sum(idx_new)/nrow(data)*100,"%"))
  print (paste0("Ratio of unassigned species:",sum(idx_new_species)/nrow(data)*100,"%"))
  return("==============")
}

# fix some of the content by manual blast

final_fix <- function(taxa_blast_final){
  TOP20 = taxa_blast_final %>% 
    dplyr::select(c("Row.names","Taxon_blast_revised","Consensus","Improve"))
  
  # We decide to only correct the top20 unassigned feature
  
  df_correct = TOP20[TOP20$Improve == "Improve(0.97-5hit)"|TOP20$Improve == "Improve(0.99-1hit)",]
  
  print(paste0("the final corrected number of top20 unassigned ",nrow(df_correct)))
  print(df_correct[,"Taxon_blast_revised"])
  
  idx_replace = df_final$Row.names%in%df_correct$Row.names
  
  df_final <- left_join(df_final,df_correct,by="Row.names")
  
  df_final[idx_replace,"Taxon"] <- df_final[idx_replace,"Taxon_blast_revised"]
  
  
  df_final <- df_final %>% dplyr::select(!"Taxon_blast_revised") %>% arrange(desc(Frequency)) %>% 
    replace_na(list("Improve"="No change")) %>%
    dplyr::select(-Blank)
  
  a = Ratio_of_unassigned(df_final[1:500,])
  
  write.table(df_final,file="../results/3.Taxonomy_ana/Corrected_Abundance_taxa_table.csv",row.names = F,sep=",")
}

final_fix(unassigned_final_correct)


# ==============Deal with different database NCBI,Greengene,SILVA

# taxonomy_greengene   # the dataframe of GreeenGene
taxonomy_greengene =  read.delim ("../results/greengene/Taxonomy_export/taxonomy.tsv")%>%
  mutate(Taxon_GreenGene = str_replace(GreenGene_2_SILVA(.),"k__","d__")) %>% 
  dplyr::select("Feature.ID", "Taxon_GreenGene", "Confidence") %>% rename(Taxon=`Taxon_GreenGene`)


taxonomy_NCBI = read.delim ("../results/blast/0.97-5hit/taxonomy.tsv") %>%
  mutate(Taxon_NCBI = str_replace(GreenGene_2_SILVA(.),"k__","d__")) %>% 
  dplyr::select("Feature.ID", "Taxon_NCBI", "Consensus") %>% rename(Taxon=`Taxon_NCBI`,
                                                             Confidence=`Consensus`)


# ============= In last par, we have combine the NCBI result with SILVA and GreenGene

# if the lowest taxa of SILVA is identical with those in the NCBI and GreenGene

# slice()
# pull()

xtab_set <- function(A,B){
  both    <-  union(A,B)
  inA     <-  both %in% A
  inB     <-  both %in% B
  return(table(inA,inB))
}

score_the_SILVA <- function(df_taxa,df_taxa2){
  #df_taxa  : silva
  #df_taxa2 : NCBI or GreenGene
  df_taxa_final <- df_taxa
  stopifnot(sum(row.names(df_taxa) == row.names(df_taxa2)) == nrow(df_taxa))#check if the index is the same
  
  df_taxa <- df_taxa %>% rownames_to_column("Order") %>%
    separate(col = Taxon,into = c("kindom","phylum","class","order","family","genus","species"),sep = "; ",fill = "right" ) %>% 
    dplyr::select(-Confidence,-Order) %>% arrange("Feature.ID") %>% column_to_rownames("Feature.ID")
  
  print (df_taxa %>% slice(1:10) %>% rownames_to_column("Feature.ID") %>% pull("Feature.ID"))
  
  df_taxa2 <- df_taxa2 %>% rownames_to_column("Order") %>%
    separate(col = Taxon,into = c("kindom","phylum","class","order","family","genus","species"),sep = "; ",fill = "right" ) %>% 
    dplyr::select(-Confidence,-Order)  %>% arrange("Feature.ID") %>% column_to_rownames("Feature.ID")
  
  print (df_taxa2 %>% slice(1:10) %>% rownames_to_column("Feature.ID") %>% pull("Feature.ID"))
  
  score = c()
  rank_depth = c()
  Evaluation <- c()
  for (i in seq(nrow(df_taxa2))){
    A = unlist(df_taxa[i,])
    B = unlist(df_taxa2[i,])
    score = c(score,sum(diag(xtab_set(A,B))))
    rank_depth = c(rank_depth,length(A))
    
    deepest_of_A = length(na.omit(A))
    deepest_of_B = length(na.omit(B))
    
    # compare the depth first
    if (length(na.omit(A)) >= length(na.omit(B))){
      #print ("depth of the SILVA is better than NCBI or GreenGene" )
      check_depth =  deepest_of_B
      evaluation =  unlist(df_taxa[i,])[check_depth] == unlist(df_taxa2[i,])[check_depth]
      
    } else {
      #print ("depth of NCBI or GreenGene is better than the SILVA" )
      check_depth =  deepest_of_A
      evaluation = unlist(df_taxa[i,])[check_depth] == unlist(df_taxa2[i,])[check_depth]
    }
    
    if(evaluation){Evaluation <- c(Evaluation,"Consistent")}else{Evaluation <- c(Evaluation,"Inconsistent")}
  }
  
  df_taxa_final[,"score"]=score
  df_taxa_final[,"rank_depth"]=rank_depth
  df_taxa_final[,"Evaluation"]=Evaluation
  
  df_taxa_final <- df_taxa_final %>% dplyr::select(!c("Taxon","Confidence"))
  
  return (df_taxa_final)
}


#score the dissimilarity of greengene and SILVA

df_evaluation <- score_the_SILVA(taxonomy,taxonomy_greengene)
table(df_evaluation$Evaluation)
#Consistent Inconsistent 
#2189         1395 
mean(df_evaluation$score)
#4.054408

df_evaluation2 <- score_the_SILVA(taxonomy,taxonomy_NCBI)
table(df_evaluation2$Evaluation)
#Consistent Inconsistent 
#107         3477 
mean(df_evaluation2$score) 
#1.09375





#Rename
colnames(taxonomy_greengene) <- c("Feature.ID", "Taxon_GreenGene", "Confidence_GreenGene") 

# taxonomy        # the dataframe of SILVA

colnames(taxonomy) <- c("Feature.ID", "Taxon_SILVA", "Confidence_SILVA") 


# taxonomy_blast  # the dataframe of blast

colnames(taxonomy_NCBI) <- c("Feature.ID", "Taxon_NCBI", "Consensus_NCBI") 



combine_taxo = left_join(taxonomy_greengene,left_join(taxonomy,taxonomy_NCBI,by="Feature.ID"),by="Feature.ID")




#  ============ Deal with the blast result downloaded from the Website ========


#label = unassigned %>% select("Feature.ID","Frequency","Taxon","Confidence")

#blast_dir = "../results/blast/"

#blast_result = c()
#for (i in 1:length(list.files(blast_dir))){
#  df = read.csv(paste0(blast_dir, list.files(blast_dir)[i])) %>% mutate(Feature.ID = unlist(strsplit(list.files(blast_dir)[i],split = ".",fixed = T))[1])
#  blast_result =  rbind(blast_result,df[1,])
#}


#Top10_unassigned_genus_blast = left_join(blast_result,label,by = "Feature.ID") %>% 
#  select("Taxon", "Frequency" ,  "Confidence"  ,"Description","Scientific.Name", "Max.Score" ,
#         "Total.Score","Query.Cover","E.value","Per..ident",     
#        "Acc..Len", "Accession","Feature.ID"   )

#write.table(Top10_unassigned_genus_blast ,file="../results/blast/Top10_unassigned_genus_blast.xls",row.names = F)




