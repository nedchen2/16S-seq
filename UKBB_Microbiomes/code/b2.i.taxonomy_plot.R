require(tidyverse)
require(ggthemes)
#require(pheatmap)
require(reshape2)
require(RColorBrewer)


# ==========

library(qiime2R)

metadata<- read.table("./sample-metadata.tsv",header = T) %>% dplyr::rename(SampleID = `sample_name`)
feature_table<-read_qza("../results/4.Diversity_ana/core-metrics-results/rarefied_table.qza")$data
taxonomy<-read_qza("../results/4.Diversity_ana/taxonomy-corrected.qza")$data
relative.feature.table <- apply(feature_table, 2, function(x) x/sum(x)*100) #convert to percent

#problem here: use the table after taxa collapse might be better 
plot_percent_heatmap <- function(SVs = relative.feature.table,topn = 30){

  SVsToPlot <- data.frame(MeanAbundance=rowMeans(SVs)) %>% #find the average abundance of a SV
    rownames_to_column("Feature.ID") %>%
    arrange(desc(MeanAbundance)) %>%
    top_n(topn, MeanAbundance) %>%
    pull(Feature.ID) #extract only the names from the table
  
  data <- SVs %>% as.data.frame() %>%
    rownames_to_column("Feature.ID") %>%
    gather(-Feature.ID, key="SampleID", value="Abundance") %>%
    mutate(Feature.ID=if_else(Feature.ID %in% SVsToPlot,  Feature.ID, "Remainder")) %>% #flag features to be collapsed
    group_by(SampleID, Feature.ID) %>%
    summarize(Abundance=sum(Abundance))%>%
    left_join(metadata) %>%
    mutate(NormAbundance=log10(Abundance+0.01)) %>% # do a log10 transformation after adding a 0.01% pseudocount. Could also add 1 read before transformation to percent
    left_join(taxonomy) %>%
    mutate(Feature=gsub("[dpcofgs]__", "", Taxon)) %>% 
    replace_na(.,list(Feature = "others")) # trim out leading text from taxonomy string
    
  p <- ggplot(data = data, aes(x=SampleID, y=Feature, fill=NormAbundance)) +
    geom_tile() +
    facet_grid(~`Species`, scales="free_x") +
    theme_q2r() +
    theme(axis.text.x=element_text(angle=90, hjust=1,size = 4,face = "bold")) + 
    labs(title = paste0("Top ", topn, " Abundance Heatmap"), y = "Taxonomy") +
    scale_fill_viridis_c(name="log10(% Abundance)") 
  p
  ggsave(plot = p,filename = "../results/4.Diversity_ana/top_features_heatmap.pdf", height=4, width=11, device="pdf") # save a PDF 3 inches by 4 inches
}

#plot_percent_heatmap(topn = 10)

# for taxon
plot_percent_heatmap_taxon <- function(SVs = relative.feature.table,topn = 30){
  
  SVsToPlot <- data.frame(MeanAbundance=rowMeans(SVs)) %>% #find the average abundance of a SV
    rownames_to_column("Feature.ID") %>%
    arrange(desc(MeanAbundance)) %>%
    top_n(topn, MeanAbundance) %>%
    pull(Feature.ID) #extract only the names from the table
  
  data <- SVs %>% as.data.frame() %>%
      rownames_to_column("Feature.ID") %>%
      gather(-Feature.ID, key="SampleID", value="Abundance") %>%
      mutate(Feature.ID=if_else(Feature.ID %in% SVsToPlot,  Feature.ID, "Remainder")) %>%  #flag features to be collapsed
      left_join(taxonomy) %>%
      group_by(SampleID,Taxon) %>%
      summarize(Abundance=sum(Abundance))%>% 
    left_join(metadata) %>%
    mutate(NormAbundance=log10(Abundance+0.01)) %>% # do a log10 transformation after adding a 0.01% pseudocount. Could also add 1 read before transformation to percent
    mutate(Feature=gsub("[dpcofgs]__", "", Taxon)) %>% 
    replace_na(.,list(Feature = "others")) # trim out leading text from taxonomy string
  
  Taxon_num <- length(unique(data$Taxon)) - 1
  
  p <- ggplot(data = data, aes(x=SampleID, y=Feature, fill=NormAbundance)) +
    geom_tile() +
    facet_grid(~`Species`, scales="free_x") +
    theme_q2r() +
    theme(axis.text.x=element_text(angle=90, hjust=1,size = 4,face = "bold")) + 
    labs(title = paste0("Top ", topn," (",Taxon_num," unique taxon)", " Abundance Heatmap"), y = "Taxonomy") +
    scale_fill_viridis_c(name="log10(% Abundance)") 
  p
  ggsave(plot = p,filename = paste0("../results/4.Diversity_ana/top",topn,"_features_heatmap.pdf"), height=4, width=11, device="pdf") # save a PDF 3 inches by 4 inches
}
plot_percent_heatmap_taxon(topn = 10)
plot_percent_heatmap_taxon(topn = 50)

# ============ stackbarplot 
colors = brewer.pal(11,"Paired")
colors = 

plot_stack_barplot <- function(SVs = relative.feature.table, subgroup = FALSE,topn = 10){ 
  
  SVsToPlot <- data.frame(MeanAbundance=rowMeans(SVs)) %>% #find the average abundance of a SV
    rownames_to_column("Feature.ID") %>%
    arrange(desc(MeanAbundance)) %>%
    top_n(topn, MeanAbundance) %>%
    pull(Feature.ID) #extract only the names from the table
  
  data <- SVs %>% as.data.frame() %>%
    rownames_to_column("Feature.ID") %>%
    gather(-Feature.ID, key="SampleID", value="Abundance") %>%
    mutate(Feature.ID=if_else(Feature.ID %in% SVsToPlot,  Feature.ID, "Remainder")) %>% #flag features to be collapsed
    group_by(SampleID, Feature.ID) %>%
    summarize(Abundance=sum(Abundance))%>%
    left_join(metadata) %>%
    mutate(NormAbundance=log10(Abundance+0.01)) %>% # do a log10 transformation after adding a 0.01% pseudocount. Could also add 1 read before transformation to percent
    left_join(taxonomy) %>%
    mutate(Feature=gsub("[dpcofgs]__", "", Taxon)) %>% 
    tidyr::replace_na(.,list(Feature = "others")) # trim out leading text from taxonomy string
  
  if (!subgroup){
    p=ggplot(data, aes(x = SampleID, y = Abundance , fill = Feature)) +
      geom_bar(stat = "identity",
               position = "fill",
               width = 0.7) +
      scale_y_continuous(labels = scales::percent) +
      #facet_grid(~ group, scales = "free_x", switch = "x") +
      theme(strip.background = element_blank()) +
      theme(axis.ticks.x = element_blank(), axis.text.x = element_blank()) +
      xlab("Samples") + ylab("Percentage (%)") +
      theme_classic() + theme(axis.text.x = element_text(
        angle = 90,
        vjust = 1,
        hjust = 1,
        size = 5
      )) +
      theme(text = element_text(family = "sans", size = 10)) +
      scale_fill_manual(values = colors)
    #scale_fill_brewer(palette="Set3")
  }else{
    p=ggplot(data, aes(x = SampleID, y = Abundance , fill = Feature)) +
      geom_bar(stat = "identity",
               position = "fill",
               width = 0.7) +
      scale_y_continuous(labels = scales::percent) +
      facet_grid(~ Species, scales = "free_x", switch = "x") +
      theme(strip.background = element_blank()) +
      theme(axis.ticks.x = element_blank(), axis.text.x = element_blank()) +
      xlab("Samples") + ylab("Percentage (%)") + labs(fill = "Taxonomy") +
      theme_classic() + theme(axis.text.x = element_text(
        angle = 90,
        vjust = 1,
        hjust = 1,
        size = 4
      )) +
      theme(text = element_text(family = "sans", size = 10)) +
      scale_fill_manual(values = colors)
  }
 
  return(p)
}

plot_stack_barplot(SVs= relative.feature.table)

colors = rainbow(11)

plot_stack_barplot2 <- function(SVs =  relative_genus_abundances * 100, subgroup = FALSE,topn = 10){ 
  
  SVsToPlot <- data.frame(MeanAbundance=rowMeans(SVs)) %>% #find the average abundance of a SV
    rownames_to_column("Feature.ID") %>%
    arrange(desc(MeanAbundance)) %>%
    top_n(topn, MeanAbundance) %>%
    pull(Feature.ID) #extract only the names from the table
  
  data <- SVs %>% as.data.frame() %>%
    rownames_to_column("Feature.ID") %>%
    gather(-Feature.ID, key="SampleID", value="Abundance") %>%
    mutate(Feature.ID=if_else(Feature.ID %in% SVsToPlot,  Feature.ID, "others")) %>% #flag features to be collapsed
    group_by(SampleID, Feature.ID) %>%
    summarize(Abundance=sum(Abundance))%>%
    left_join(metadata) %>%
    mutate(NormAbundance=log10(Abundance+0.01)) # do a log10 transformation after adding a 0.01% pseudocount. Could also add 1 read before transformation to percent
    # trim out leading text from taxonomy string
  
  if (!subgroup){
    p=ggplot(data, aes(x = SampleID, y = Abundance , fill = Feature.ID)) +
      geom_bar(stat = "identity",
               position = "fill",
               width = 0.7) +
      scale_y_continuous(labels = scales::percent) +
      #facet_grid(~ group, scales = "free_x", switch = "x") +
      theme(strip.background = element_blank()) +
      theme(axis.ticks.x = element_blank(), axis.text.x = element_blank()) +
      xlab("Samples") + ylab("Percentage (%)") +
      theme_classic() + theme(axis.text.x = element_text(
        angle = 90,
        vjust = 1,
        hjust = 1,
        size = 5
      )) +
      theme(text = element_text(family = "sans", size = 10)) +
      scale_fill_manual(values = colors)
    #scale_fill_brewer(palette="Set3")
  }else{
    p=ggplot(data, aes(x = SampleID, y = Abundance , fill = Feature.ID)) +
      geom_bar(stat = "identity",
               position = "fill",
               width = 0.7) +
      scale_y_continuous(labels = scales::percent) +
      facet_grid(~ Species, scales = "free_x", switch = "x") +
      theme(strip.background = element_blank()) +
      theme(axis.ticks.x = element_blank(), axis.text.x = element_blank()) +
      xlab("Samples") + ylab("Percentage (%)") + labs(fill = "Taxonomy") +
      theme_classic() + theme(axis.text.x = element_text(
        angle = 90,
        vjust = 1,
        hjust = 1,
        size = 4
      )) +
      theme(text = element_text(family = "sans", size = 20)) +
      scale_fill_manual(values = colors) + coord_cartesian(expand = FALSE)
  }
  
  return(p)
}


p = plot_stack_barplot2(SVs= relative_genus_abundances * 100, subgroup = T)

ggsave(plot = p,filename = "../results/7.Final_graph//top10_features_Genus.pdf", height=4, width=11, device="pdf") # save a PDF 3 inches by 4 inches


# =========== phyloseq and enterotypes =============
library("phyloseq")

metadata<- read.table("./sample-metadata.tsv",header = T) %>% dplyr::rename(SampleID = `sample_name`)

#filtered  - rarefied    --- literature
physeq<-qza_to_phyloseq(
  features="../results/4.Diversity_ana/core-metrics-results/rarefied_table.qza",
  tree="../results/4.Diversity_ana/rooted-tree.qza",
  "../results/4.Diversity_ana/taxonomy-corrected.qza"
)

sample_data <- sample_data(metadata %>% column_to_rownames("SampleID") )

physeq1 = merge_phyloseq(physeq, sample_data)
physeq1


otu_table_export0 <- as.data.frame(otu_table(physeq1))
taxa_table_export0 <- as.data.frame(tax_table(physeq1))
merged_table0 <- merge(otu_table_export0,taxa_table_export0,by.x="row.names",by.y = "row.names") 

write.csv(merged_table0,"../../Data_Available//Original_ASV_abundance_data.csv",row.names = T)


test <- as.data.frame(otu_table(physeq1))
sum(rowSums(test) == 1)

# plot_richness
plot_richness(physeq1, x="Species")

# plot network
ig <- make_network(physeq1, "samples",max.dist = 0.8)
plot_network(ig, physeq1, color="Species", shape="CollectionSite", line_weight=0.3)





# ================================Enterotypes ==================
test_phyloseq2 = tax_glom(physeq1, "Phylum")
otu_table_export2 <- as.data.frame(otu_table(test_phyloseq2))
taxa_table_export2 <- as.data.frame(tax_table(test_phyloseq2))
merged_table2 <- merge(otu_table_export2,taxa_table_export2,by.x="row.names",by.y = "row.names") %>% 
  dplyr::select(!c("Row.names", "Kingdom", "Genus"  ,"Class" ,"Order", "Family","Species")) %>% column_to_rownames("Phylum")

write.csv(merged_table2,"../results/5.Enterotype/phylum_abundance_data.csv",row.names = T,col.names = T)



# Phyloseq to get the genus level table
test_phyloseq = tax_glom(physeq1, "Genus")
test_phyloseq

plot_richness(test_phyloseq, x="Species",measures = "shannon")

otu_table_export <- as.data.frame(otu_table(test_phyloseq))
taxa_table_export <- as.data.frame(tax_table(test_phyloseq))

merged_table <- merge(otu_table_export,taxa_table_export,by.x="row.names",by.y = "row.names") %>% 
  dplyr::select(!c("Row.names", "Kingdom", "Phylum"  ,"Class" ,"Order", "Family","Species")) %>% column_to_rownames("Genus")

write.csv(merged_table,"../results/5.Enterotype/genus_abundance_data.csv",row.names = T,col.names = T)

relative_genus_abundances <- apply(merged_table, 2, function(x) x/sum(x)) 
# Bombus terrestis 
# ====================== just a test

sampleInHorRod <-  metadata %>% filter(Species != "B.terrestris") %>% pull("SampleID")

sampleInTerrestris <- metadata %>% filter(Species == "B.terrestris") %>% pull("SampleID")




idx = colnames(relative_genus_abundances) %in% sampleInTerrestris

idx = colnames(relative_genus_abundances) %in% sampleInHorRod 

relative_genus_abundances <- relative_genus_abundances[,idx]
#===========================

colSums(relative_genus_abundances)
ncol(relative_genus_abundances)

# ============ TEST FOR USING COLLAPSE RESULT FROM QIIME
re

df_taxa_test <-read_qza("../results/6.Network_analysis/table-l6.qza")$data


#relative_genus_abundances <- apply(df_taxa_test , 2, function(x) x/sum(x)) 



# ==== denoise
noise.removal <- function(dataframe, percent=0.01, top=NULL){
  dataframe->Matrix
  bigones <- rowSums(Matrix)*100/(sum(rowSums(Matrix))) > percent 
  Matrix_1 <- Matrix[bigones,]
  print(percent)
  return(Matrix_1)
}

data.denoized=noise.removal(relative_genus_abundances, percent=0.01)



#========= dist matrix

dist.JSD <- function(inMatrix, pseudocount=0.000001, ...) {
  KLD <- function(x,y) sum(x *log(x/y))
  JSD<- function(x,y) sqrt(0.5 * KLD(x, (x+y)/2) + 0.5 * KLD(y, (x+y)/2))
  matrixColSize <- length(colnames(inMatrix))
  matrixRowSize <- length(rownames(inMatrix))
  colnames <- colnames(inMatrix)
  resultsMatrix <- matrix(0, matrixColSize, matrixColSize)
  
  inMatrix = apply(inMatrix,1:2,function(x) ifelse (x==0,pseudocount,x))
  
  for(i in 1:matrixColSize) {
    for(j in 1:matrixColSize) { 
      resultsMatrix[i,j]=JSD(as.vector(inMatrix[,i]),
                             as.vector(inMatrix[,j]))
    }
  }
  colnames -> colnames(resultsMatrix) -> rownames(resultsMatrix)
  as.dist(resultsMatrix)->resultsMatrix
  attr(resultsMatrix, "method") <- "dist"
  return(resultsMatrix) 
}

data.dist=dist.JSD(data.denoized)

pam.clustering=function(x,k) { # x is a distance matrix and k the number of clusters
  require(cluster)
  cluster = as.vector(pam(as.dist(x), k, diss=TRUE)$clustering)
  return(cluster)
}

hierarchical.clustering <- function(dist_object = data.dist, k = 3){
  upgma <- hclust(dist_object)
  cluster =   cutree(tree = upgma,k = k)
  return(cluster)
}




# === evaluate the enterotypes number
require(clusterSim)
#nclusters = index.G1(t(data.denoized), data.cluster, d = data.dist, centrotypes = "centroids")
nclusters=NULL

for (k in 1:20) { 
  if (k==1) {
    nclusters[k]=NA 
  } else {
    data.cluster_temp=pam.clustering(data.dist, k)
    nclusters[k]=index.G1(t(data.denoized),data.cluster_temp,  d = data.dist,
                          centrotypes = "medoids")
  }
}
plot(nclusters, type="h", xlab="k clusters", ylab="CH index")



# ==== use silhouette method to see the clusters performance
nclusters=NULL
for (k in 1:20) { 
  if (k==1) {
    nclusters[k]=NA 
  } else {
    data.cluster_temp=pam.clustering(data.dist, k)
    nclusters[k]= mean(silhouette(data.cluster_temp, data.dist)[,3])
  }
}
plot(nclusters, type="h", xlab="k clusters", ylab="silhouette")



# ==== cluster


  data.cluster=pam.clustering(data.dist, k=3)
  
  #data.cluster = hierarchical.clustering(k = 4)
  
  table(data.cluster)
  
  metadata<- read.table("./sample-metadata.tsv",header = T) %>% dplyr::rename(SampleID = `sample_name`)
  
  Enterotypes <- as.data.frame(t(rbind(colnames(data.denoized),data.cluster)))
  
  colnames(Enterotypes) <- c("SampleID","Enterotype")
  
  Enterotypes <- Enterotypes %>% mutate(Enterotype = paste0("Enterotype_",Enterotype))
  
  metadata <- merge(metadata,Enterotypes,by.x = "SampleID",by.y = "SampleID")
  
  metadata[metadata$Enterotype == "Enterotype_1","SampleID"]
  
  
  write.csv(metadata, file = "../results/7.Final_graph/Enterotype_HorRud.csv",row.names = F)
  
  write.csv(metadata, file = "../results/7.Final_graph/Enterotype_terrestris.csv",row.names = F)
  
  Enterotype_n <- metadata %>% group_by(Enterotype) %>% count() %>% dplyr::rename(total = n)
  
  Enterotype_species <- metadata %>% group_by(Enterotype) %>% count(Species)%>% left_join(Enterotype_n) %>% mutate(Percent = n/total)
  write.csv(Enterotype_species, file = "../results/7.Final_graph/Enterotype_species_table.csv")
  
  table(metadata$Species,metadata$Enterotype)
  
  Enterotype_sites <- metadata %>% group_by(Enterotype) %>% count(CollectionSite)%>% left_join(Enterotype_n) %>% mutate(Percent = n/total)
  write.csv(Enterotype_sites, file = "../results/7.Final_graph/Enterotype_sites_table.csv")
  
  table(metadata$CollectionSite,metadata$Enterotype)

  chisq.test(metadata$Species,metadata$Enterotype,correct = T,simulate.p.value = T)
  
  table(metadata$Species,metadata$Enterotype)
  
  chisq.test(metadata$CollectionSite,metadata$Enterotype)
  
  
Arsenophonus_abundance = as.data.frame(t(relative_genus_abundances)*100) %>% rownames_to_column("sample_name") %>% 
  dplyr::select("sample_name","Arsenophonus","Gilliamella","Snodgrassella","Apibacter") %>% rename( SampleID = sample_name)


metadata <- metadata %>% left_join(Arsenophonus_abundance)

write.csv(metadata, file = "../results/5.Enterotype/metadata_Microbiome.csv")
write.csv(metadata, file = "./metadata_Microbiome.csv",row.names = F)


obs.silhouette=mean(silhouette(data.cluster, data.dist)[,3])

obs.silhouette

# terrestris 0.4235



# ================ have a look based on the enterotypes

# cluster the sample



# ===seeThename 


require(ade4)

obs.pca=dudi.pca(data.frame(t(data.denoized)), scannf=F, nf=10)
obs.bet=bca(obs.pca, fac=as.factor(data.cluster), scannf=F, nf=10) 
s.class(obs.bet$ls, fac=as.factor(data.cluster), grid=F)
s.class(obs.bet$ls, fac=as.factor(data.cluster), grid=F, cell=0, cstar=0, col=c(4,2,3))

obs.pcoa=dudi.pco(data.dist, scannf=F, nf=5)


# === draw PCOA by my script

ordination_result = obs.pcoa$li %>% rownames_to_column("SampleID") %>% dplyr::rename(PC1 = A1,
                                                                              PC2 = A2)

plot_pcoA_with_beta <- function(ordination_result,groupID = "Enterotype" ){
  #faith_pd_vector<-read_qza("../results/4.Diversity_ana/core-metrics-results/faith_pd_vector.qza")$data %>% rownames_to_column("SampleID") 
  #ordination_result[,"SampleID"] <- as.character(ordination_result[,"SampleID"] )
  
  combin_matrix <- ordination_result%>%
    dplyr::select(SampleID, PC1, PC2) %>%
    left_join(metadata) 
  
  col_vector = str_replace(colnames(combin_matrix),groupID,replacement = "group")
  
  colnames(combin_matrix) = col_vector
  
  p = ggplot(data = combin_matrix, aes(x=PC1, y=PC2)) +
    geom_point(mapping = aes( color= Species,shape = CollectionSite),size = 5,alpha = 0.5 ) + #alpha controls transparency and helps when points are overlapping
    theme_q2r() + 
    scale_shape_manual(values=c(15,16,17,18), name="CollectionSite") + #see http://www.sthda.com/sthda/RDoc/figure/graphs/r-plot-pch-symbols-points-in-r.png for numeric shape codes
    #scale_size_continuous(name="Phylogenetic Diversity") +
    scale_color_gdocs(name = "Species") +
    stat_ellipse(aes(linetype=group),level = 0.95) + 
    labs(x="PCo1",
         y="PCo2") +
    geom_hline(yintercept = 0,lwd = 0.5)+
    geom_vline(xintercept = 0,lwd = 0.5)+ 
    theme(text = element_text(size = 20))
  p
  
  return(p)
}

p = plot_pcoA_with_beta(ordination_result = ordination_result)

ggsave(p,filename = "../results/7.Final_graph/Enterotype_PCoA_All.png",dpi = "print",width = 8)

s.class(obs.pcoa$li, fac=as.factor(data.cluster), grid=T,label = c(1,2,3,4,5))


#sample_data <- sample_data(metadata %>% column_to_rownames("SampleID"))
#physeq2 = merge_phyloseq(test_phyloseq, sample_data)
#GPr<-transform_sample_counts(physeq2 , function(x) x / sum(x) )
#gpsfb = subset_taxa(GPr, Genus=="Snodgrassella")
#plot_bar(gpsfb, "Enterotype", "Abundance", title="Snodgrassella")
#gpsfb = subset_taxa(physeq2, Genus=="Lactobacillus")
#plot_bar(gpsfb, "Enterotype", "Abundance", title="Lactobacillus")
#gpsfb = subset_taxa(physeq2, Genus=="Serratia")
#plot_bar(gpsfb, "Enterotype", "Abundance", title="Serratia")
#gpsfb = subset_taxa(physeq2, Genus=="Hafnia-Obesumbacterium")
#plot_bar(gpsfb, "Enterotype", "Abundance", title="Hafnia-Obesumbacterium")


# =============== box plot for single genus

build_data_frame <- function( genus.name){
  
  return(as.data.frame(t(relative_genus_abundances)) %>% 
           rownames_to_column("SampleID") %>% 
           pivot_longer(cols = rownames(relative_genus_abundances), names_to = "Genus",values_to = "Relative.abundance") %>%
           left_join(metadata) %>% subset(Genus == genus.name))
  
}

list_of_genrea = c("Snodgrassella","Lactobacillus", "Hafnia-Obesumbacterium","Serratia","Gilliamella")

Enterotype_box <- function(genus_name,group = "Enterotype") {
  
  df=build_data_frame(genus.name = genus_name )
  
  df<-df %>% dplyr::select("Relative.abundance",all_of(group)) 
  
  colnames(df) <- c("Relative.abundance","Group")
  
  model = aov(Relative.abundance ~ Group, data=df)
  # Tukey-significance test
  Tukey_HSD = TukeyHSD(model, ordered = TRUE, conf.level = 0.95)
  
  Tukey_HSD_table = as.data.frame(Tukey_HSD$Group)
  
  print (Tukey_HSD)
  
  All_group_test = kruskal.test(Relative.abundance ~ Group, data=df)
  
  PAIRWISE_RESULT = pairwise.wilcox.test(x = df$Relative.abundance,g=df$Group,p.adjust.method = "BH")
  print(PAIRWISE_RESULT)
  
  p = ggplot(df, aes(x=Group, y=Relative.abundance, color=Group)) +
    geom_boxplot(alpha=1, outlier.shape = NA, outlier.size=0, size=0.7, width=0.5, fill="transparent") +
    labs(x="Groups", y="relative abundance",title = paste0(genus_name,"| p-value:",All_group_test$p.value) , color = group) + theme_classic() + scale_fill_gdocs()+
    geom_jitter(position=position_jitter(0.17), size=1, alpha=0.7)+
    theme(text=element_text(family="sans", size=14))
  p
  return(p)
}

Enterotype_box(genus_name = "Gilliamella")
Enterotype_box(genus_name = "Snodgrassella")
Enterotype_box(genus_name = "Lactobacillus")
Enterotype_box(genus_name = "Apibacter")
Enterotype_box(genus_name = "Hafnia-Obesumbacterium")
Enterotype_box(genus_name = "Arsenophonus") # siginificant while clustering through B.R and B.H
Enterotype_box(genus_name = "Fructobacillus")
Enterotype_box(genus_name = "Escherichia-Shigella")
Enterotype_box(genus_name = "Pseudomonas")
Enterotype_box(genus_name = "Bifidobacterium")
Enterotype_box(genus_name = "Pantoea")
Enterotype_box(genus_name = "Ralstonia")
Enterotype_box(genus_name = "Vibrio")

as.data.frame(rowSums(merged_table)) %>% top_n(30) %>% arrange(desc(`rowSums(merged_table)`))

Enterotype_box(genus_name = "Gilliamella",group = "Species" )
Enterotype_box(genus_name = "Snodgrassella",group = "Species" )
Enterotype_box(genus_name = "Lactobacillus",group = "Species" )
Enterotype_box(genus_name = "Apibacter",group = "Species" )

# "Hafnia-Obesumbacterium","Serratia" Not significant

Enterotype_box(genus_name = "Fructobacillus",group = "Species" ) # No difference
Enterotype_box(genus_name = "Arsenophonus",group = "Species" )  # No difference
Enterotype_box(genus_name = "Pseudomonas",group = "Species" ) #No
Enterotype_box(genus_name = "Pantoea",group = "Species" )#No
Enterotype_box(genus_name = "Bifidobacterium",group = "Species" )#No
Enterotype_box(genus_name = "Dubosiella",group = "Species" )#No
Enterotype_box(genus_name = "Neokomagataea",group = "Species" )#No
Enterotype_box(genus_name = "Ralstonia",group = "Species" ) # No difference

Enterotype_box(genus_name = "Pseudomonas",group = "CollectionSite" )

metadata2 <- metadata 
metadata <- metadata %>% mutate(infection = paste0(gut_parasite_richness,"_infection")) %>% dplyr::select(-Enterotype) %>% rename(Enterotype=infection)
metadata <- metadata2


metadata <- read.csv(file =  "./metadata_Microbiome.csv") %>% mutate(SampleID = as.character(SampleID))

Enterotype_heatmap_taxon <- function(SVs = relative_genus_abundances,topn = 50){
  
  SVs = SVs * 100
  
  SVsToPlot <- data.frame(MeanAbundance=rowMeans(SVs)) %>% #find the average abundance of a SV
    rownames_to_column("Feature.ID") %>%
    arrange(desc(MeanAbundance)) %>%
    top_n(topn, MeanAbundance) %>%
    pull(Feature.ID) #extract only the names from the table
  
  data <- SVs %>% as.data.frame() %>%
    rownames_to_column("Feature.ID") %>%
    gather(-Feature.ID, key="SampleID", value="Abundance") %>%
    mutate(Genus=if_else(Feature.ID %in% SVsToPlot,  Feature.ID, "Others")) %>%  #flag features to be collapsed
    group_by(SampleID,Genus) %>%
    summarize(Abundance=sum(Abundance))%>% 
    left_join(metadata) %>%
    mutate(NormAbundance=log10(Abundance+0.01)) %>% # do a log10 transformation after adding a 0.01% pseudocount. Could also add 1 read before transformation to percent
    tidyr::replace_na(.,list(Genus = "others")) # trim out leading text from taxonomy string
  

  p <- ggplot(data = data, aes(x=SampleID, y=Genus, fill=NormAbundance)) +
    geom_tile() +
    facet_grid(~`Species`, scales="free_x") +
    theme_q2r() +
    theme(axis.text.x=element_text(angle=90, hjust=1,size = 4,face = "bold")) + 
    labs( y = "Genus") +
    scale_fill_viridis_c(name="log10(% Abundance)") 
  p
  return(p)
}

p = Enterotype_heatmap_taxon(topn = 10)

ggsave(p, filename = "../results/7.Final_graph/Heatmap_top10_genus.png", dpi = "print",width =  8, height = 3, scale= 0.8)


# ===== calculate the mean of the enterotype

substr_relative_genus_abundances_heatmap <- as.data.frame(substr_table) %>% rownames_to_column("Genus") %>%
  pivot_longer(cols = colnames(substr_table),names_to = "SampleID",values_to = "relative_abundance") %>% 
  left_join(Enterotypes) %>% group_by(Genus,Enterotype) %>% summarise(mean_abundance = mean(relative_abundance)) %>%
  pivot_wider(names_from = Enterotype,values_from = mean_abundance) %>% column_to_rownames("Genus")

write.csv(substr_relative_genus_abundances_heatmap,"../results/7.Final_graph/Core_bacteria_Enterotype_mean_relative_abundance.csv")

data <- substr_relative_genus_abundances_heatmap %>% rownames_to_column("Genus") %>% pivot_longer(cols = colnames(substr_relative_genus_abundances_heatmap),names_to = "Enterotype",values_to = "Abundance") %>% mutate(NormAbundance=log10(Abundance+0.01)) %>% 
  mutate(Percent_Abundance = round(Abundance * 100,2),
         Enterotype_display = str_remove(Enterotype,pattern = "Enterotype_")) 

p <- ggplot(data = data, aes(x= Enterotype_display, y=Genus, fill=NormAbundance)) +
  geom_tile(color = "white",lwd = 1.5,linetype = 1) +
  coord_fixed() + 
  geom_text(aes(label =  Percent_Abundance),color = "white",size = 4) + 
  theme(axis.text.x=element_text(size = 8,face = "bold")) + 
  scale_fill_gradient(low = "navy",high = "firebrick") +
  labs( y = "Genus", x = "") +
  guides(fill= guide_colorbar(title="log10(% Abundance)")) + theme(legend.position = "bottom")
  

ggsave(p,filename = "../results/7.Final_graph/Enterotype_MeanAbundance_Heatmap.png",width = 8, height = 6)


# ===== Top Abundant

Enterotype_heatmap_taxon(topn = 10)

Enterotype_heatmap_taxon(topn = 25)




# ====== if the infection will influence 

# prevalence
chisq.test(metadata$Apicystis_binomial,metadata$Enterotype)
chisq.test(metadata$Crithidia_binomial,metadata$Enterotype)
chisq.test(metadata$Nosema_binomial,metadata$Enterotype)
chisq.test(metadata$gut_parasite_richness,metadata$Enterotype)



# parasite load



# =========== Enterotype Percentage






# =========== 

# ===========SIMPER Analysis

library("vegan")

merged_table = df_taxa_test

test_simper = as.data.frame(t(merged_table)) 

metadata2 <- metadata[colnames(merged_table)%in%metadata$SampleID,] %>% column_to_rownames("SampleID")

metadata2 <- metadata2[row.names(test_simper),] %>% rownames_to_column("SampleID")

test_simper <- test_simper %>% rownames_to_column("SampleID") %>% dplyr::select(!"SampleID")

simper(test_simper,metadata2$Enterotype)

# ======== list of biomarker from simper

list_of_biomarker_simper = c("Lactobacillus", "Snodgrassella", "Arsenophonus", "Apibacter", "Pseudomonas", "Fructobacillus", "Gilliamella")
#list_of_biomarker = c("Lactobacillus", "Snodgrassella", "Arsenophonus", "Apibacter", "Fructobacillus", "Gilliamella")

list_of_biomarker = c("Lactobacillus", "Snodgrassella","Bifidobacterium", "Arsenophonus", "Apibacter", "Fructobacillus", "Gilliamella","Pantoea","Ralstonia","Vibrio")

#Pantoea 

# =========== Random Forest
require(randomForest)
library(caret)
require(tidyverse)

data.denoized=noise.removal(relative_genus_abundances, percent=0.01)

rownames(data.denoized) <- str_remove(rownames(data.denoized),pattern = "\\[")
rownames(data.denoized) <- str_remove(rownames(data.denoized),pattern = "\\]")
rownames(data.denoized) <- str_replace(rownames(data.denoized),pattern = "-",replacement = "_")
rownames(data.denoized) <- str_replace(rownames(data.denoized),pattern = "-",replacement = "_")
rownames(data.denoized) <- str_replace(rownames(data.denoized),pattern = "-",replacement = "_")

df_Crithidia <- df_meta %>% rownames_to_column("SampleID")%>% dplyr::select(SampleID,CrithidiaInfection)
df_Nosema <- df_meta %>% rownames_to_column("SampleID")%>% dplyr::select(SampleID,NosemaInfection)
Enterotypes <- df_meta%>% rownames_to_column("SampleID")%>% dplyr::select(SampleID,Enterotype)

#======= Crithidia infection
data <- as.data.frame(t(data.denoized)) %>% 
  rownames_to_column("SampleID") %>% 
  left_join(df_Crithidia) %>%  column_to_rownames("SampleID") %>% dplyr::rename( Group = CrithidiaInfection) %>% mutate(Group = as.factor(Group))

#=======  Nosema infection
data <- as.data.frame(t(data.denoized)) %>% 
  rownames_to_column("SampleID") %>% 
  left_join(df_Nosema) %>%  column_to_rownames("SampleID") %>% dplyr::rename( Group = NosemaInfection) %>% mutate(Group = as.factor(Group))

#======= Enterotype
data <- as.data.frame(t(data.denoized)) %>% 
  rownames_to_column("SampleID") %>% 
  left_join(Enterotypes) %>%  column_to_rownames("SampleID") %>% dplyr::rename( Group = Enterotype) %>% mutate(Group = as.factor(Group))

set.seed(222)
ind <- sample(2, nrow(data), replace = TRUE, prob = c(0.8, 0.2))
train <- data[ind==1,]
test <- data[ind==2,]

rf <- randomForest(Group~., data=train, proximity=TRUE) 
print(rf)


# PREDICT AND CONFUSION MATRIX 
p2 <- predict(rf, test)
confusionMatrix(p2, test$Group)

#importance(rf)

varImpPlot(rf,
           sort = T,
           n.var = 10,
           main = "Top 10 - Variable Importance")

list_of_biomarker_rf <- varImp(rf) %>% arrange(desc(Overall)) %>% top_n(10) %>% rownames_to_column("Genus") %>% pull(Genus)

list_of_biomarker <- intersect(list_of_biomarker_simper,list_of_biomarker_rf)


# =====  SELECT THE DRIVING MICROBIOME


substr_relative_genus_abundances <- apply(merged_table, 2, function(x) x/sum(x)) 

substr_table <- substr_relative_genus_abundances[list_of_biomarker,]
hist(colSums(substr_table))
Unexpected_samples <- colnames(substr_table)[colSums(substr_table) < 0.9]

infected_sample <- metadata[metadata$gut_parasite_richness != 0,"SampleID"]
length(infected_sample)
length(Unexpected_samples)

sum(Unexpected_samples %in% infected_sample)/length(Unexpected_samples)
# 0.7142857

Enterotype_heatmap_taxon(SVs =substr_table, topn = length(list_of_biomarker))

# ===== 







# =========== Venn Graph

data <- metadata %>% dplyr::select(Nosema_binomial,Crithidia_binomial,Apicystis_binomial,Conopid_larvae)
data <- metadata %>% dplyr::select(Nosema_binomial,Crithidia_binomial,Apicystis_binomial)

data <- metadata %>% filter(SampleID != "5",SampleID !="Blank") %>% dplyr::select(parasite_counts_available,RADseq_available,amplicon_16S_available) %>%
  mutate(Parasite.counts = ifelse(parasite_counts_available == "YES", 1, 0),
         RADseq = ifelse(RADseq_available == "YES", 1, 0),
         AmpliconSeq = ifelse(amplicon_16S_available == "YES", 1, 0)) %>%  dplyr::select(Parasite.counts,RADseq,AmpliconSeq)


library("venn")
library("ggplot2")
pdf("../results/7.Final_graph/ExperVennGraph.pdf", width=4, height=4, onefile=F)
venn(data,zcolor="red,yellow,blue,green,orange,purple",box =F,ggplot = T,ilabels = T, plotsize = 10,ellipse = T)
dev.off()

png("./VennGraph.png", width=3000, height=3000, res=300)
venn(data,zcolor="red,yellow,blue,green,orange,purple",box = T,plotsize = 10,ggplot = T,sncs = 6)
dev.off()





# ============== Deeper in parasite

faith_pd_vector<- read_qza("../results/4.Diversity_ana/core-metrics-results/faith_pd_vector.qza")$data %>% rownames_to_column("SampleID") 
shannon_vector <- read_qza("../results/4.Diversity_ana/core-metrics-results/shannon_vector.qza")$data %>% rownames_to_column("SampleID") 

lm_genus <- metadata %>% dplyr::select("SampleID","Nosema_total_per_bee" ,"Crithidia_total_per_bee","Apicystis_total_per_bee","Enterotype")

Genus <- log10(as.data.frame(t(substr_table))) %>% rownames_to_column("SampleID") 

lm_genus <- left_join(lm_genus,Genus) %>% left_join(.,shannon_vector) %>% left_join(.,faith_pd_vector) %>% column_to_rownames("SampleID")


# ===== low - parasite - load sample and particular taxa
suppressPackageStartupMessages(library(ggrepel))
sample <- data[data$log.Parasite.Load < 4.5,] %>% dplyr::select("SampleID","Enterotype")
colnames(sample) <- c("SampleID","Enterotype_label")

sample <- data %>% dplyr::select("SampleID","Enterotype") %>% subset(Enterotype == "Enterotype_2")

Genus_2 <- as.data.frame(t(substr_table)) %>% rownames_to_column("SampleID") %>% left_join(sample) %>% dplyr::select("Gilliamella","Apibacter","Snodgrassella","Enterotype_label")

# ====== 2d plot
p1 = ggplot(Genus_2,aes(x=log10(Gilliamella),y=log10(Apibacter))) + 
  geom_point(size = 3) + 
  geom_label_repel(aes(label = Enterotype_label), size = 2,, force = 1,color = "black") + 
  theme_clean() 

p2 = ggplot(Genus_2,aes(x=log10(Gilliamella),y=log10(Snodgrassella))) + 
  geom_point(size = 3) + 
  geom_label_repel(aes(label = Enterotype_label), size = 2,, force = 1,color = "black") + 
  theme_clean() 

p3 = ggplot(Genus_2,aes(x=log10(Apibacter),y=log10(Snodgrassella))) + 
  geom_point(size = 3) + 
  geom_label_repel(aes(label = Enterotype_label), size = 2,, force = 1,color = "black") + 
  theme_clean() 


# ====== 3d plot
library("plotly")

Genus_2 <- Genus_2 %>% replace_na(replace = list(Enterotype_label="Others(Low Parasite Loads)"))

plot_ly(Genus_2,x=~Apibacter,y=~Snodgrassella,z=~Gilliamella,color = ~Enterotype_label,colors = c('#BF382A', '#0C4B8E',"green")) %>% add_markers(size=1)

# =========== Parasite_loads_Enterotype_boxplot

Parasite_loads_Enterotype_boxplot <- function( group =  "Crithidia_total_per_bee" ){
  
  data <- lm_genus %>% 
    dplyr::select("Nosema_total_per_bee" ,"Crithidia_total_per_bee","Apicystis_total_per_bee","Enterotype") %>% 
    rownames_to_column("SampleID") %>%
    pivot_longer(cols = c("Nosema_total_per_bee" ,"Crithidia_total_per_bee","Apicystis_total_per_bee"), names_to = "Parasite.Status",values_to = "Parasite.Loads") %>% 
    mutate(log.Parasite.Load = log10(Parasite.Loads)) %>% 
    mutate_if(is.numeric, list(~na_if(.,-Inf))) %>% 
    mutate_if(is.numeric, list(~na_if(.,Inf))) %>% 
    replace_na(., replace = list( log.Parasite.Load = 0)) %>%
    subset(Parasite.Status == group) %>% subset( log.Parasite.Load != 0 )
  
  Infected_size = data %>%  count(Enterotype)
  
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
    labs(x= group, y= bquote(paste("Parasite Loads",(log[10]))), color =  "Enterotype",subtitle = paste("kruskal.test: ",All_group_test$p.value)) + theme_classic() + scale_fill_gdocs()+
    geom_jitter(position=position_jitter(0.17), size=1, alpha=0.7)+
    theme(text=element_text(family="sans", size=14)) + scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07","yellow"))
  
  ggsave(plot = p,filename = paste("../results/5.Enterotype/Parasite_loads_Enterotype_", group,"_boxplot.pdf"), height=4, width=11, device="pdf") # save a PDF 3 inches by 4 inches
}


Parasite_loads_Enterotype_boxplot()

Parasite_loads_Enterotype_boxplot(group = "Nosema_total_per_bee")

Parasite_loads_Enterotype_boxplot(group = "Apicystis_total_per_bee")

# =========== Parasite load heatmap

Enterotype_heatmap_parasite <- function(){
  
  data <- lm_genus %>% 
    dplyr::select("Nosema_total_per_bee" ,"Crithidia_total_per_bee" ,"Apicystis_total_per_bee" ,"Enterotype" ) %>% 
    rownames_to_column("SampleID") %>%
    pivot_longer(cols = c("Nosema_total_per_bee" ,"Crithidia_total_per_bee","Apicystis_total_per_bee"), names_to = "Parasite.Status",values_to = "Parasite.Loads") %>% 
    mutate(log.Parasite.Load = log10(Parasite.Loads)) %>% 
    mutate_if(is.numeric, list(~na_if(.,-Inf))) %>% 
    mutate_if(is.numeric, list(~na_if(.,Inf))) %>% 
    replace_na(., replace = list( log.Parasite.Load = 0)) 
  
  
  p1 <- ggplot(data = data, aes(x = SampleID , y = Parasite.Status, fill =  log.Parasite.Load)) +
    geom_tile() +
    facet_grid(~`Enterotype`, scales="free_x") +
    theme_q2r() +
    theme(axis.text.x=element_text(angle=90, hjust=1,size = 4,face = "bold")) + 
    scale_fill_viridis_c() 
  
  p1
  
  data <- data %>% subset( log.Parasite.Load != 0 ) %>% group_by(Enterotype,Parasite.Status) %>% 
    summarise(Sum.Parasite.Load = sum(Parasite.Loads),
              SampleSize = n()) %>% 
    mutate(log.Mean.load = log10(Sum.Parasite.Load/SampleSize)) 
  
  
  p2 <- ggplot(data = data, aes(x = Enterotype , y = Parasite.Status, fill = log.Mean.load)) +
    geom_tile() +
    facet_grid(~`Enterotype`, scales="free_x") +
    theme_q2r() +
    theme(axis.text.x=element_text(angle=90, hjust=1,size = 4,face = "bold")) + 
    scale_fill_viridis_c() 
  p2
  
}

# ======== Upset plot of Enterptype and Species and Collection Sites
#,"gut_parasite_richness"
list_of_Species = subset(metadata,select = c("SampleID" ,"Species"))
colnames(list_of_Species) <- c("SampleID","Set")

list_of_CollectionSite = subset(metadata,select = c("SampleID" ,"CollectionSite"))
colnames(list_of_CollectionSite) <- c("SampleID","Set")

list_of_Enterotype =  subset(metadata,select = c("SampleID" ,"Enterotype"))      
colnames(list_of_Enterotype ) <- c("SampleID","Set")

list_of_Nosema = subset(metadata,select = c("SampleID" ,"Nosema_binomial" )) %>% mutate(Nosema_binomial = paste0("Nosema_richness(",Nosema_binomial,")"))  
colnames(list_of_Nosema ) <- c("SampleID","Set")

list_of_Crithidia = subset(metadata,select = c("SampleID" , "Crithidia_binomial")) %>% mutate(Crithidia_binomial = paste0("Crithidia_richness(",Crithidia_binomial,")"))  
colnames(list_of_Crithidia) <- c("SampleID","Set")

list_of_Apicystis = subset(metadata,select = c("SampleID" ,"Apicystis_binomial")) %>% mutate(Apicystis_binomial = paste0("Apicystis_richness(",Apicystis_binomial,")"))  
colnames(list_of_Apicystis ) <- c("SampleID","Set")

list_of_parasite = subset(metadata,select = c("SampleID" ,"gut_parasite_richness")) %>% mutate(gut_parasite_richness = paste0("Parasite_richness(",gut_parasite_richness,")"))  
colnames(list_of_parasite ) <- c("SampleID","Set")



list_of_Crithidia = subset(metadata,select = c("SampleID" , "Crithidia_binomial")) 
list_of_Apicystis = subset(metadata,select = c("SampleID" ,"Apicystis_binomial")) 
list_of_Nosema = subset(metadata,select = c("SampleID" ,"Nosema_binomial" ))


UpSET_compare <- function(lista = list_of_Enterotype , listb= list_of_Species ,listc = NULL,listd=NULL,list_binomial1 = NULL,list_binomial2=NULL,list_binomial3=NULL){
  
  test_for_upset <- rbind.data.frame(lista,listb,listc,listd)
  
  data_for_upset <- as.data.frame.matrix(table(test_for_upset)) 
  
  data_for_upset[data_for_upset>0]=1
  
  data_for_upset <- data_for_upset %>% rownames_to_column("SampleID")
  
  if (is.null(list_binomial1) ){
    data <- data_for_upset%>% column_to_rownames("SampleID")
    
  }else{
    data <- left_join(left_join(left_join(data_for_upset,list_binomial1),list_binomial2),list_binomial3) %>% column_to_rownames("SampleID")
  }
  
  library("UpSetR")
  datnum <- length(colnames(data))
  pdf("../results/5.Enterotype/VennGraph.pdf", width=8, height=6, onefile=F)
  p = upset(data, nsets=datnum, nintersects=NA, number.angles=0,point.size=1, line.size=0.5, mainbar.y.label="Intersection sample number", sets.x.label="Total sample number", text.scale=c(1.3,1.3,1,1,1.3,1), mb.ratio=c(0.55,0.45), order.by="freq", show.numbers="yes", sets.bar.color=c("black"),main.bar.color = "navy")
  print (p)
  dev.off()
}

UpSET_compare(listc =list_of_CollectionSite )
UpSET_compare(lista = list_of_Enterotype , listb= list_of_Species )
UpSET_compare(listb = list_of_CollectionSite)
UpSET_compare(listb = list_of_parasite)

UpSET_compare(lista =  list_of_Enterotype, listb = NULL,list_binomial1 = list_of_Apicystis,list_binomial2=list_of_Crithidia,list_binomial3=list_of_Nosema)


# =========== linear regression

plot_lm_scatter <- function(data = lm_genus, response = "Crithidia_total_per_bee",explantory = "shannon_entropy",Group = "Enterotype", confidence_interval = T,group = F) {
  library(ggpubr)
  library(ggplot2)
  library(ggthemes)
  
  df_subset <- subset(data, select=c(response,explantory,Group))
  
  
  #log10 transform to response variable
  df_subset[,response] <- log10(df_subset[,response]) 
  colnames(df_subset) <- c("response","explantory","group")
  
  #replace the "-/Inf"
  if( "Inf" %in% df_subset[,"response"] | "-Inf" %in% df_subset[,"response"]){
    tmp <- df_subset[which(df_subset[,"response"] != "Inf" & df_subset[,"response"] != "-Inf"),]
    max = max(tmp[,"response"])
    min = min(tmp[,"response"])
    df_subset[,"response"] <-  as.numeric(sub("-Inf", 0, df_subset[,"response"]))
    df_subset[,"response"] <- as.numeric(sub("Inf", max, df_subset[,"response"]))
  }
  

  if( "Inf" %in% df_subset[,"explantory"] | "-Inf" %in% df_subset[,"explantory"]){
    tmp <- df_subset[which(df_subset[,"explantory"] != "Inf" & df_subset[,"explantory"] != "-Inf"),]
    max = max(tmp[,"explantory"])
    min = min(tmp[,"explantory"])
    df_subset[,"explantory"] <- as.numeric(sub("-Inf", min,  df_subset[,"explantory"]))
    df_subset[,"explantory"] <- as.numeric(sub("Inf", max, df_subset[,"explantory"]))
  }
  # remove uninfected samples
  df_subset <- subset(df_subset, response !=0)
  
  fit_model <- try(lm(data = df_subset,response ~  explantory),silent = T)
  if (as.vector(summary(fit_model))[2] != "try-error"){
  
    report_tmp <- as.data.frame(coef(summary(fit_model)))
    report_tmp["Parameter"] <- rownames(report_tmp) #add parameter to one col
    # convert to longer sheet
    report_tmp2 <- report_tmp 
    report_tmp2["Explantory"] <- explantory
    report_tmp2["Response"]   <- response
    report_tmp2["SampleSize"] <- nrow(df_subset)
    report_tmp2["Adjusted_R"] <- summary(fit_model)$adj.r.squared
    report_tmp2["p_value"]    <- pf(summary(fit_model)$fstatistic[1],summary(fit_model)$fstatistic[2],summary(fit_model)$fstatistic[3],lower.tail = F)
  }else{
    cat("========================")
    stop("please check the table")
  }
  
  if (group){
   p = ggplot(df_subset,aes(x=explantory,y=response,color=group)) + 
      geom_point(aes(color=group ),size = 3) + 
      geom_smooth(method = "lm",formula = y~x, se = confidence_interval) +
      stat_cor(data = df_subset, method = "pearson",label.sep = " ",label.x.npc = 0.8) + 
      theme_clean() + 
      labs(y = bquote(paste(.(response),(Log[10]))),x = explantory, title = "linear regression")
    

  }else{
    p = ggplot(df_subset,aes(x=explantory,y=response)) + 
      geom_point(size = 3) + 
      geom_smooth(method = "lm",formula = y~x, se = confidence_interval) +
      stat_cor(data = df_subset, method = "pearson",label.sep = " ",label.x.npc = 0.8) + 
      theme_clean() + 
      labs(y =  bquote(paste(.(response),(Log[10]))),x = explantory, title = "linear regression")
  }
  
  return(list(picture = p, report = report_tmp2))
}


plot_lm_scatter(explantory = "faith_pd",group = T)
plot_lm_scatter(explantory = "faith_pd",group = F)

plot_lm_scatter(group = T)
plot_lm_scatter(group = F)

plot_lm_scatter(response = "Nosema_total_per_bee" ,group = T)

# Apibacter and Correaltion

#metadata <- 




# ============ Correlation


list_of_biomarker = c( "Snodgrassella","Bifidobacterium", "Apibacter", "Gilliamella")

substr_table <- substr_relative_genus_abundances[list_of_biomarker,]

require(PerformanceAnalytics)
Genus <- log10(as.data.frame(t(substr_table))+0.0001)

Genus <- Genus[list_of_biomarker]

chart.Correlation.linear <-function (R, histogram = TRUE, method=c("pearson", "kendall", "spearman"), ...)
  { # @author R Development Core Team
    # @author modified by Peter Carl & Marek Lahoda
    # Visualization of a Correlation Matrix. On top the (absolute) value of the correlation plus the result 
    # of the cor.test as stars. On botttom, the bivariate scatterplots, with a linear regression fit. 
    # On diagonal, the histograms with probability, density and normal density (gaussian) distribution.
    
    x = checkData(R, method="matrix")
    
    if(missing(method)) method=method[1] #only use one
    cormeth <- method
    
    # Published at http://addictedtor.free.fr/graphiques/sources/source_137.R
    panel.cor <- function(x, y, digits=2, prefix="", use="pairwise.complete.obs", method=cormeth, cex.cor, ...)
    {
      usr <- par("usr"); on.exit(par(usr))
      par(usr = c(0, 1, 0, 1))
      r <- cor(x, y, use=use, method=method) # MG: remove abs here
      txt <- format(c(r, 0.123456789), digits=digits)[1]
      txt <- paste(prefix, txt, sep="")
      if(missing(cex.cor)) cex <- 2.5
      
      test <- cor.test(as.numeric(x),as.numeric(y), method=method)
      # borrowed from printCoefmat
      Signif <- symnum(test$p.value, corr = FALSE, na = FALSE,
                       cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                       symbols = c("***", "**", "*", ".", " "))
      # MG: add abs here and also include a 30% buffer for small numbers
      text(0.5, 0.5, txt, cex = cex)
      text(.8, .8, Signif, cex=cex, col=2)
    }
    
    #remove method from dotargs
    dotargs <- list(...)
    dotargs$method <- NULL
    rm(method)
    
    hist.panel = function (x, ...=NULL ) {
      par(new = TRUE)
      hist(x,
           col = "light gray",
           probability = TRUE,
           axes = FALSE,
           main = "",
           breaks = "FD")
      lines(density(x, na.rm=TRUE),
            col = "red",
            lwd = 1)
      # adding line representing density of normal distribution with parameters correponding to estimates of mean and standard deviation from the data 
      ax.x = seq(min(x), max(x), 0.1)                                                  # ax.x containts points corresponding to data range on x axis
      density.est = dnorm(ax.x, mean = mean(x), sd = sd(x))   # density corresponding to points stored in vector ax.x 
      lines(ax.x, density.est, col = "blue", lwd = 1, lty = 1)                                # adding line representing density into histogram
      rug(x)
    }
    
    # Linear regression line fit over points
    reg <- function(x, y, ...) {
      points(x,y, cex=0.5, pch = 19,...)
      abline(lm(y~x), col = "red") 
    }
    
    # Draw the chart
    if(histogram)
      pairs(x, gap=0, lower.panel=reg, upper.panel=panel.cor, diag.panel=hist.panel)
    else
      pairs(x, gap=1, lower.panel=reg, upper.panel=panel.cor) 
  }


pdf(file = "../results/7.Final_graph/CorrelationPlot.pdf",width = 9)
chart.Correlation.linear(Genus, histogram = F)
dev.off()


# ================
# install.packages("devtools")
devtools::install_github("Hy4m/linkET", force = TRUE)

library(linkET)
library(dplyr)
## 
## Attaching package: 'dplyr'
## The following objects are masked from 'package:stats':
## 
##     filter, lag
## The following objects are masked from 'package:base':
## 
##     intersect, setdiff, setequal, union
library(ggplot2)

environmental <- metadata %>% dplyr:: select("SampleID","Crithidia_total_per_bee","Nosema_total_per_bee","Apicystis_total_per_bee") %>% column_to_rownames("SampleID") 
idx = row.names(environmental)

#list_of_biomarker = c( "Snodgrassella","Bifidobacterium","Gilliamella", "Apibacter","Pantoea","Ralstonia","Vibrio","Lactobacillus","Fructobacillus" )

#idx_of_bact = colnames(spec)[order(colnames(spec) %in% list_of_biomarker)]

spec <- as.data.frame(t(data.denoized))[idx,] * 100


mantel <- mantel_test(spec , environmental,spec_select = list(Bac01 = c("Lactobacillus","Fructobacillus","Arsenophonus"),
                                                             Bac02 = c("Snodgrassella","Bifidobacterium","Gilliamella","Apibacter"),
                                                             Bac03 = c("Pantoea","Ralstonia","Vibrio","Lactobacillus","Fructobacillus"))) %>% 
  mutate(rd = cut(r, breaks = c(-Inf, 0, 0.2, Inf), # 对相关系数进行分割，便于映射大小
                  labels = c("< 0", "0 - 0.2", ">= 0.2")),
         pd = cut(p, breaks = c(-Inf, 0.01, 0.05, Inf), # 对P值进行分割，便于映射颜色
                  labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))
## `mantel_test()` using 'bray' dist method for 'spec'.
## `mantel_test()` using 'euclidean' dist method for 'env'.

# 画图
p = qcorrplot(correlate(environmental), type = "lower", diag = FALSE) +
  geom_square(size = 1) +
  geom_couple(aes(colour = pd, size = rd), # 这行代码是关键
              data = mantel, 
              curvature = nice_curvature(),nudge_x = 0.1) +
  # 下面就是各种颜色和名称设置
  scale_fill_gradientn(colours = RColorBrewer::brewer.pal(11, "RdBu")) +
  scale_size_manual(values = c(0.5, 1, 2)) +
  scale_colour_manual(values = color_pal(3)) +
  guides(size = guide_legend(title = "Mantel's r",
                             override.aes = list(colour = "grey35"), 
                             order = 2),
         colour = guide_legend(title = "Mantel's p", 
                               override.aes = list(size = 3), 
                               order = 1),
         fill = guide_colorbar(title = "Pearson's r", order = 3))  


ggsave(p, filename = "../../UKBB_Microbiomes/results/7.Final_graph/Mantel_p.png",dpi = "print", width = 8)


# =====================================================
# Enterotype-specific test

Genus <- log10(as.data.frame(t(substr_table))+0.0001) %>% rownames_to_column("SampleID") %>% left_join(Enterotypes) 

list_of_biomarker = c( "Snodgrassella","Bifidobacterium","Gilliamella", "Apibacter", "Lactobacillus","Fructobacillus", "Arsenophonus","Pantoea","Ralstonia","Vibrio")

Enterotype1_Genus <- Genus %>% filter(Enterotype == "Enterotype_1")
                                      
list_of_biomarker = c( "Snodgrassella","Bifidobacterium","Gilliamella", "Apibacter", "Lactobacillus","Fructobacillus", "Arsenophonus","Pantoea","Ralstonia","Vibrio")
Enterotype1_Genus<- Enterotype1_Genus[list_of_biomarker]

pdf(file = "../results/7.Final_graph/CorrelationPlotEnterotype1.pdf",width = 9)
chart.Correlation.linear(Enterotype1_Genus, histogram = F)
dev.off()


Enterotype2_Genus <- Genus %>% filter(Enterotype == "Enterotype_2")
list_of_biomarker = c("Snodgrassella","Bifidobacterium","Gilliamella", "Apibacter")
Enterotype2_Genus<- Enterotype2_Genus[list_of_biomarker]

pdf(file = "../results/7.Final_graph/CorrelationPlotEnterotype2.pdf",width = 9)


chart.Correlation.linear(Enterotype2_Genus, histogram = F)

dev.off()

Enterotype3_Genus <- Genus %>% filter(Enterotype == "Enterotype_3")
list_of_biomarker = c("Lactobacillus","Pantoea","Ralstonia","Vibrio")

Enterotype3_Genus<- Enterotype3_Genus[list_of_biomarker]


pdf(file = "../results/7.Final_graph/CorrelationPlotEnterotype3.pdf",width = 9)
chart.Correlation.linear(Enterotype3_Genus, histogram = F)
dev.off()


# ================ differential abundance analysis of genus

library("DESeq2")

# revise the metadata

sample_data <- sample_data(metadata %>% column_to_rownames("SampleID") )

# add pseudo number
otu_table <- otu_table(transform_sample_counts(test_phyloseq, function(x) x + 0.01 ))

test_phyloseq = merge_phyloseq(test_phyloseq, sample_data)

test_phyloseq = merge_phyloseq(test_phyloseq, otu_table)

# DEseq
#diagdds = phyloseq_to_deseq2(test_phyloseq, ~ Species) # this might be error-prone . We may have a much more precise base-control group

diagdds = phyloseq_to_deseq2(test_phyloseq, ~ Enterotype) 

diagdds = DESeq(diagdds, test="Wald", fitType="parametric")

res = results(diagdds, cooksCutoff = FALSE)

sigtab = cbind(as(res, "data.frame"), as(tax_table(test_phyloseq)[rownames(res), ], "matrix"))
head(sigtab)

# ===== function of 

library("ggplot2")

# point plot
phylum_and_genus <-function() {
  theme_set(theme_bw())
  scale_fill_discrete <- function(palname = "Set1", ...) {
    scale_fill_brewer(palette = palname, ...)
  }
  sigtab <- sigtab[sigtab$pvalue < 0.01,]
  # Phylum order
  x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
  x = sort(x, TRUE)
  sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
  # Genus order
  x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
  x = sort(x, TRUE)
  sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))
  ggplot(sigtab, aes(x=Genus, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
    theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))
}



# volcano


plot_volcano <- function(){
  logFC_cutoff = log(2, 2)
  pval_cutoff <- 0.05
  
  sigtab$change = as.factor(ifelse(sigtab$pvalue < 0.05 & abs(sigtab$log2FoldChange) > logFC_cutoff, ifelse(sigtab$log2FoldChange > logFC_cutoff ,'Up','Down'),'Filtered'))
  
  
  l <- list(a = "p-value")
  print (l)
  eq <- substitute(-log[10]~a, l)
  
  #找到范围
  x=ifelse(max(c(floor(max(sigtab$log2FoldChange)),abs(floor(min(sigtab[,"log2FoldChange"])))))>15,15,max(c(floor(max(sigtab$log2FoldChange)),abs(floor(min(sigtab[,"log2FoldChange"]))))))
  print (xlim)
  if ( is.na(x) ){stop("请检查表格")}
  
  
  labels = c("Snodgrassella","Lactobacillus","Gilliamella")
  
  sigtab$label <-''
  
  sigtab[sigtab$Genus %in% labels, "label"] = as.character(sigtab[sigtab$Genus %in% labels,]$Genus)
  
  suppressPackageStartupMessages(library(ggrepel))
  
  g = ggplot(data=sigtab, aes(x=sigtab$log2FoldChange, y=-log10(sigtab$pvalue), color=change)) + 
    geom_point( size=2,alpha=0.5) + 
    theme_set(theme_set(theme_bw(base_size=15))) + 
    theme(legend.title=element_blank()) +   
    xlab(expression(paste(log[2], " Fold change"))) + 
    ylab(as.expression(eq)) + 
    theme(plot.title = element_text(hjust = 0.5)) + 
    guides(shape=guide_legend(override.aes=list(size=4))) +
    scale_colour_manual(values = c("Up"=c("#8B0000"), "Down"=c("darkblue"),"Filtered"=c("grey")),na.translate=FALSE) + 
    geom_text_repel(aes(label = label), size = 4, force = 1,color = "black") + 
    theme(text=element_text(size=14,family="ArialMT"))+
    scale_x_continuous(limits = c(-x,x))+
    theme(panel.grid = element_blank()) + 
    geom_hline(yintercept = -log(pval_cutoff, 10), linetype = "solid", color = c("grey"), size = 0.5) + 
    geom_vline(xintercept = c(-logFC_cutoff, logFC_cutoff), linetype = "solid", color = c("grey"), size = 0.5)
  
  g
  return(g)
}


# ======== co occurrance 

library(vegan)
library(bipartite)

C.score(relative_genus_abundances,normalise = F)

oecosimu(relative_genus_abundances,bipartite::C.score,"swap",burnin = 100,thin = 10,statistic = "evals",nsimul = 1000,parallel = 3)


