# clean

if(!require("dplyr")) {install.packages("dplyr")}
require(dplyr)
require(tidyverse)

#previous metadata

df_old= read.table("map_UKBB.txt",sep = "\t",header = T)

df_old = df_old[(df_old$Sample.ID != ""),]

# We found that the barcode of Index2i5Sequence TAGATCGC
index = df_old$Index2i5Sequence == "TAGATCGC"
replace = "GCGTAAGN"

df_old$Index2.modified = str_replace(df_old$Index2i5Sequence,pattern = "TAGATCGC",replacement = "GCGTAAGN")

df.old.revised  <- df_old %>% mutate(BarcodeSequence=paste0(substr(df_old$Index1i7Sequence,1,7),substr(df_old$Index2.modified,1,7)),
                     LinkerPrimerSequence=Linker1PrimerSequence,
                     sample_name=Sample.ID) %>% dplyr::select("sample_name","BarcodeSequence","LinkerPrimerSequence","CollectionSite","Species","Infection","Infection.intensity","Description")  %>% 
                      replace_na(list(CollectionSite = "Blank",Species = "Blank",Infection = "Blank",Infection.intensity = "Blank"))
  


write.table(df.old.revised,file = "../results/1.Quality_Control/Demultiplexing/sample-metadata.tsv",row.names = F,sep="\t")

write.table(df.old.revised,file = "./sample-metadata.tsv",row.names = F,sep="\t")

#to discover the pattern of barcode sequence
#zcat lane1_Undetermined_L001_R2_1.fastq.gz | grep '^@' | cut -d: -f10 | sort | uniq -c | sort - -k1rn | head -n120 > result.txt

df_top120Barcode = read.delim("../16sRaw/result.txt",header = F) %>% 
  mutate(repeats = extract_numeric(V1),
  BarcodeSequence = str_sub(V1,-14)) %>%
  select("repeats","BarcodeSequence" ) %>% 
  full_join(df.old.revised,by = c("BarcodeSequence"))

write.table(df_top120Barcode,file = "./Top120Barcode.tsv",row.names = F,sep="\t")





#updated metadata
df_new = read.table("micro_map.txt",sep="\t",header = T)

df_new$Index2.modified = str_replace(df_new$Index2i5Sequence,pattern = "TAGATCGC",replacement = "GCGTAAGN")

df.new.revised <-  df_new %>% mutate(BarcodeSequence=paste0(substr(df_new$Index1i7Sequence,1,7),substr(df_new$Index2.modified,1,7)),
                                LinkerPrimerSequence=Linker1PrimerSequence,
                                sample_name=as.character(Bee_ID)) %>% dplyr::select( "sample_name" , "BarcodeSequence",  "LinkerPrimerSequence"  ,"radseq_ID_list", "sibling_sets", "Site", "conopid_larvae","Apicystis_binomial",    
                                                                                     "mean_Apicystis_sp_ul" , "Crithidia_binomial" ,"mean_Crithidia_cell_ul", "Nosema_binomial" , "mean_Nosema_sp_ul"  ,  
                                                                                     "para_richness")

write.table(df.new.revised,file = "./sample-metadata-NewMappingFile.tsv",row.names = F,sep="\t")



df.merged = df.new.revised %>%  left_join(df.old.revised,by = "sample_name")

df.merged$BarcodeSequence == df.merged$BarcodeSequence_old
str(df.merged)

#write.table(df.merged,file = "../results/1.Quality_Control/Demultiplexing/sample-metadata.tsv",row.names = F,sep="\t")

#write.table(df.merged,file = "./sample-metadata.tsv",row.names = F,sep="\t")


df_final = read.table("combined_mapping_file.csv",sep=",",header = T) %>% 
  mutate(Subgenus  = ifelse(test = finalised_species == "B.terrestris", yes = "Bombus.sensu.stricto", no= "Bombus.Megabombus"))%>% 
  rename(sample_name =`Bee_ID` ,
         Species = `finalised_species`,
         CollectionSite = `Site`)   %>% 
  mutate(Interaction = paste0(Species,".",CollectionSite))

write.table(df_final,file = "./sample-metadata.tsv",row.names = F,sep="\t")


df2 <- read.csv("../code/metadata_Microbiome.csv") %>% mutate(sample_name = as.character(SampleID)) %>% dplyr::select(sample_name, Enterotype)

df<-df_final %>% mutate(Infection = ifelse(gut_parasite_richness != 0, "infected","Notinfected")) %>% left_join(df2)

write.table(df,file = "./sample-metadata2.tsv",row.names = F,sep="\t")




# =========== test if the final one is consistent with the previous one

df_test <- df.old.revised %>% dplyr::select(sample_name,BarcodeSequence) %>% rename(BarcodeSequence_old=BarcodeSequence)

df_final %>% left_join(df_test) %>% mutate(test= BarcodeSequence == BarcodeSequence_old) %>% pull(test) %>% sum()
#91


# ========== Draw a parasite infection

metadata = read.csv("../code/sample-metadata.tsv",sep ="\t") %>% rename(SampleID = sample_name)
metadata <- metadata[-91,]
# ========== parasite prevalence and species/location

list_of_species = subset(metadata,select =  c("SampleID" ,"Species")) # the name is enterotype but actuallt it is Species
colnames(list_of_species) = c("SampleID" ,"group")

list_of_sites = subset(metadata,select =  c("SampleID" ,"CollectionSite"))
colnames(list_of_sites) = c("SampleID" ,"group")

list_of_Crithidia = subset(metadata,select = c("SampleID" , "Crithidia_binomial"))
list_of_Apicystis = subset(metadata,select = c("SampleID" ,"Apicystis_binomial")) 
list_of_Nosema = subset(metadata,select = c("SampleID" ,"Nosema_binomial" ))


Enterotype_Percentage <- function(list_of_group = list_of_species  ,Name = "Species"){
  
  Crithidia_stack_plot <-left_join(list_of_group, list_of_Crithidia) %>% 
    group_by(group) %>% 
    summarise(Parasite.present = sum(Crithidia_binomial),frequency = n()) %>% 
    mutate(Percent =Parasite.present/frequency,
           no.Parasite.present = frequency - Parasite.present, 
           parasite = "Crithidia") %>% 
    pivot_longer(cols = c("Parasite.present", "no.Parasite.present"), names_to = "Parasite.Status", values_to = "Parasite.count")
  
  Apicystis_stack_plot <- left_join(list_of_group, list_of_Apicystis) %>% 
    group_by(group) %>% 
    summarise(Parasite.present = sum(Apicystis_binomial),frequency = n()) %>% 
    mutate(Percent =Parasite.present/frequency,
           no.Parasite.present = frequency - Parasite.present, 
           parasite = "Apicystis") %>% 
    pivot_longer(cols = c("Parasite.present", "no.Parasite.present"), names_to = "Parasite.Status", values_to = "Parasite.count")
  
  Nosema_stack_plot  <- left_join(list_of_group, list_of_Nosema) %>% 
    group_by(group) %>% 
    summarise(Parasite.present = sum(Nosema_binomial),frequency = n()) %>% 
    mutate(Percent =Parasite.present/frequency,
           no.Parasite.present = frequency - Parasite.present, 
           parasite = "Nosema") %>% 
    pivot_longer(cols = c("Parasite.present", "no.Parasite.present"), names_to = "Parasite.Status", values_to = "Parasite.count")
  
  data = rbind(Nosema_stack_plot,Apicystis_stack_plot,Crithidia_stack_plot)
  
  p = ggplot(data, aes(x = group, y = Parasite.count , fill = Parasite.Status)) +
    geom_bar(stat = "identity",
             position = "fill",
             width = 0.7)+
    scale_y_continuous(labels = scales::percent) +
    facet_grid(~ parasite, scales = "free_x", switch = "x") +
    theme(strip.background = element_blank()) +
    theme(axis.ticks.x = element_blank(), axis.text.x = element_blank()) +
    xlab("") + ylab("Percentage (%)") + labs(fill = "Status") +
    theme_classic() + 
    theme(text = element_text(family = "sans", size = 10)) +
    scale_fill_ordinal()
  
  ggsave(plot = p,filename = paste0("../results/7.Final_graph/parasite_prev_",Name,"_stack_plot.pdf"), height=4, width=11, device="pdf") # save a PDF 3 inches by 4 inches
}

Enterotype_Percentage()

Enterotype_Percentage(list_of_group = list_of_sites, Name = "Sites")




# box plot 

Parasite_loads_Enterotype_boxplot <- function( lm_genus = metadata ,group =  "Crithidia_total_per_bee" , group2 = "Species" ){
  
  data <- lm_genus %>% 
    dplyr::select("Nosema_total_per_bee" ,"Crithidia_total_per_bee","Apicystis_total_per_bee", group2) %>% 
    rownames_to_column("SampleID") %>%
    pivot_longer(cols = c("Nosema_total_per_bee" ,"Crithidia_total_per_bee","Apicystis_total_per_bee"), names_to = "Parasite.Status",values_to = "Parasite.Loads") %>% 
    mutate(log.Parasite.Load = log10(Parasite.Loads)) %>% 
    mutate_if(is.numeric, list(~na_if(.,-Inf))) %>% 
    mutate_if(is.numeric, list(~na_if(.,Inf))) %>% 
    replace_na(., replace = list( log.Parasite.Load = 0)) %>%
    subset(Parasite.Status == group) %>% subset( log.Parasite.Load != 0 )
  
  colnames(data)[2] <- "Enterotype"
  
  colnames(data)[3] <- "group"
  
  Infected_size = data %>% count(Enterotype)
  
  data <- data %>% left_join(Infected_size) %>% mutate(Enterotype_n = paste0(Enterotype,"(n=",n,")")) 
  
  model = aov(log.Parasite.Load ~ Enterotype, data=data)
  # Tukey-significance test
  Tukey_HSD = TukeyHSD(model, ordered = TRUE, conf.level = 0.95)
  
  Tukey_HSD_table = as.data.frame(Tukey_HSD$Enterotype)
  
  print (Tukey_HSD)
  
  All_group_test = kruskal.test(log.Parasite.Load ~ Enterotype, data=data)
  
  #pairwise.wilcox.test(data$log.Parasite.Load,data$Enterotype,p.adjust.method = "BH",paired = F,)
  
  p = ggplot(data, aes( x = Enterotype_n, y = log.Parasite.Load, color = Enterotype)) +
    geom_boxplot(alpha=1, outlier.shape = NA, outlier.size=0, size=0.7, width=0.5, fill="transparent") +
    labs(x= group, y= bquote(paste("Parasite Loads",(log[10]))), color =  group2,subtitle = paste("kruskal.test: ",All_group_test$p.value)) + theme_classic() + scale_fill_gdocs()+
    geom_jitter(position=position_jitter(0.17), size=1, alpha=0.7)+
    theme(text=element_text(family="sans", size=14)) + scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07","yellow"))
  
  ggsave(plot = p,filename = paste("../results/7.Final_graph/Parasite_loads_",group2,"_", group,"_boxplot.pdf"), height=4, width=11, device="pdf") # save a PDF 3 inches by 4 inches
}


Parasite_loads_Enterotype_boxplot()

Parasite_loads_Enterotype_boxplot(group = "Nosema_total_per_bee")

Parasite_loads_Enterotype_boxplot(group = "Apicystis_total_per_bee")

Parasite_loads_Enterotype_boxplot(group2 = "CollectionSite")


Parasite_loads_Enterotype_boxplot(group = "Nosema_total_per_bee",group2 = "CollectionSite")


Parasite_loads_Enterotype_boxplot(group = "Apicystis_total_per_bee",group2 = "CollectionSite")











