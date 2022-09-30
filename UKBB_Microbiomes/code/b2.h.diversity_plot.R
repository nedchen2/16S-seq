

# =============================Or use vegan package to see the diversity
require(ggthemes)
require(vegan)
require(tidyverse)
require(qiime2R)

if(!require(qiime2R)){
  install.packages("remotes")
  remotes::install_github("jbisanz/qiime2R")
}


path <- "../results/4.Diversity_ana/plot"
if (!dir.exists(path )){
  dir.create(path)
}

# ================= Use Vegan to calculate 
feature_abundance_with_seq <- read.csv("../results/3.Taxonomy_ana/Corrected_Abundance_taxa_table.csv",check.names = F)
feature_abundance_matrix = feature_abundance_with_seq %>% select(!c("Sequence","Taxon","Confidence","Suspectable","Frequency","Consensus","Improve")) %>% 
  column_to_rownames(var = "Row.names")
alpha_matrix = matrix()
richness <- specnumber(feature_abundance_matrix,MARGIN = 2)
invsimpson  <- diversity(feature_abundance_matrix,
                               MARGIN = 2,
                               index = "invsimpson")

simpson  <-  diversity(feature_abundance_matrix,
                                    MARGIN = 2,
                                    index = "simpson")

shannon  <- diversity(feature_abundance_matrix,
                                  MARGIN = 2,
                                  index = "shannon")

PieLou.evenness <- shannon/log(richness)

Fisher_alpha <- fisher.alpha(feature_abundance_matrix,MARGIN = 2)

alpha_matrix = as.data.frame (cbind(richness,shannon,simpson,invsimpson,PieLou.evenness,Fisher_alpha ))

# ================= import the data from qiime

faith_pd_vector<-read_qza("../results/4.Diversity_ana/core-metrics-results/faith_pd_vector.qza")$data %>% rownames_to_column("SampleID") 
shannon_vector <- read_qza("../results/4.Diversity_ana/core-metrics-results/shannon_vector.qza")$data %>% rownames_to_column("SampleID") 

alpha_matrix = as.data.frame(left_join(faith_pd_vector,shannon_vector)) %>% column_to_rownames("SampleID")

alpha_matrix %>% rownames_to_column("sample_name") %>% write_csv("../results/4.Diversity_ana/Alpha_Diversity_Matrix.csv",)

plot_alpha_barplot <- function (alpha_div = alpha_matrix ,index="shannon",groupID = "Infection"){
  
  # "CollectionSite"       "Species"              "Infection"
  
  # "richness"        "shannon"         "simpson"         "invsimpson"      "PieLou.evenness" "Fisher_alpha"
  
  # "faith_pd"        "shannon_entropy"
  
  metadata <- read.table("./sample-metadata.tsv",header = T,row.names = 1)
  
  # cross-filtering
  idx = rownames(metadata) %in% rownames(alpha_div)
  metadata = metadata[idx,,drop=F]
  
  #alpha_div = as.data.frame(alpha_div[rownames(metadata),],row.names = rownames(metadata))
  
  sampFile = as.data.frame(metadata[, groupID],row.names = row.names(metadata))
  colnames(sampFile) <- c(groupID)
  
  
  # combine alpha diversity and meta table
  df = cbind(alpha_div[rownames(sampFile),index], sampFile)
  colnames(df) = c(index,"group")
  
  # aov
  model = aov(df[[index]] ~ group, data=df)
  # Tukey-significance test
  Tukey_HSD = TukeyHSD(model, ordered = TRUE, conf.level = 0.95)
  
  Tukey_HSD_table = as.data.frame(Tukey_HSD$group)
  
  #save the result from the tukey test
  write.table(paste(date(), "\nGroup\t", groupID, "\n\t", sep=""), file=paste("alpha_boxplot_TukeyHSD.txt",sep=""),append = T, quote = F, eol = "", row.names = F, col.names = F)
  suppressWarnings(write.table(Tukey_HSD_table, file=paste("alpha_boxplot_TukeyHSD.txt",sep=""), append = T, quote = F, sep="\t", eol = "\n", na = "NA", dec = ".", row.names = T, col.names = T))
  
  
  # 函数：将Tukey检验结果P值转换为显著字母分组
  # 输入文件为图基检验结果和分组
  generate_label_df = function(TUKEY, variable){
    library(multcompView)
    # 转换P值为字母分组
    ## 提取图基检验中分组子表的第4列P adjust值
    Tukey.levels = TUKEY[[variable]][,4]
    ## multcompLetters函数将两两p值转换为字母，data.frame并生成列名为Letters的数据框
    Tukey.labels = data.frame(multcompLetters(Tukey.levels)['Letters'])
    
    # 按分组名字母顺序
    ## 提取字母分组行名为group组名
    Tukey.labels$group = rownames(Tukey.labels)
    # 按组名的字母顺序排列，默认的Levels
    Tukey.labels=Tukey.labels[order(Tukey.labels$group), ]
    return(Tukey.labels)
  }
  
  # 当只有两组时，用LSD标注字母
  if (length(unique(df$group)) == 2){
    # LSD检验，添加差异组字母
    library(agricolae)
    out = LSD.test(model, "group", p.adj="none")
    stat = out$groups
    # 分组结果添入Index
    df$stat=stat[as.character(df$group),]$groups
    # 当大于两组时，用multcompView标注字母
  }else{
    # library(multcompView)
    LABELS = generate_label_df(Tukey_HSD , "group")
    df$stat=LABELS[as.character(df$group),]$Letters
  }
  
  # 设置分组位置为各组y最大值+高的5%
  max=max(df[,c(index)])
  min=min(df[,index])
  x = df[,c("group",index)]
  y = x %>% group_by(group) %>% summarise_(Max=paste('max(',index,')',sep=""))
  y=as.data.frame(y)
  rownames(y)=y$group
  df$y=y[as.character(df$group),]$Max + (max-min)*0.05
  # head(df)
  
  # ================== calculate the sd and mean
  
  colnames(df) = c("alpha","group","stat","y")
  went = df %>% group_by(group) %>% summarize(mean = mean(alpha),SD = sd(alpha),label = unique(stat)) 
  #or
  #wen1 = as.data.frame(tapply(as.vector(as.matrix(df[1])), df$group, mean, na.rm=TRUE))
  #wen2 = as.data.frame(tapply(as.vector(as.matrix(df[1])), df$group, sd, na.rm=TRUE))
  #went = cbind(wen1, wen2)
  #colnames(went) = c("mean","SD")
  #aa = distinct(df, group, .keep_all = TRUE)
  #went$label = aa$stat[match(row.names(went),aa$group)]
  #went$Groups = row.names(went)
  
  a = max(went$mean + went$SD)*1.1
  
  # =========== barplot
  p = ggplot(went, aes(x = group, y = mean, colour= group)) +
    geom_bar(aes(colour= group, fill = group), stat = "identity", width = 0.4, position = "dodge") +
    scale_y_continuous(expand = c(0,0),limits = c(0,a)) +
    geom_errorbar(aes(ymin = mean-SD, ymax=mean+SD), colour="black", width=0.1, size = 1) +
    geom_text(aes(label = label,y = mean+SD, x = group, vjust = -0.3), color = "black") +
    theme_classic() + scale_color_gdocs()+ scale_fill_gdocs()+
    theme(text=element_text(family="sans", size=14)) + labs(y = paste0("mean of ", index ))
  p
  return(p)
  
}


# test

plot_alpha_barplot(groupID = "CollectionSite")

plot_alpha_barplot(groupID = "Species")

# ========== boxplot
#df[[index]]  can be used to aes variable to the ggplot

df_meta =  read.table("./sample-metadata.tsv",header = T,row.names = 1)

df_meta = read.csv("./metadata_Microbiome.csv",header = T,row.names = 1) %>% mutate(CrithidiaInfection = ifelse(Crithidia_binomial == 1, "infected","not infected"),
                                                                                    NosemaInfection = ifelse(Nosema_binomial == 1, "infected","not infected"))


plot_alpha_boxplot <- function (alpha_div = alpha_matrix ,index="shannon_entropy",groupID = "Species", input_metadata = df_meta , ylab = "alpha diversity"){
  
  metadata <- input_metadata
  # cross-filtering
  idx = rownames(metadata) %in% rownames(alpha_div)
  metadata = metadata[idx,,drop=F]
  
  #alpha_div = as.data.frame(alpha_div[rownames(metadata),],row.names = rownames(metadata))
  
  sampFile = as.data.frame(metadata[, groupID],row.names = row.names(metadata))
  colnames(sampFile) <- c(groupID)
  
  
  # combine alpha diversity and meta table
  df = cbind(alpha_div[rownames(sampFile),index], sampFile)
  colnames(df) = c(index,"group")
  
  # aov
  model = aov(df[[index]] ~ group, data=df)
  # Tukey-significance test
  Tukey_HSD = TukeyHSD(model, ordered = TRUE, conf.level = 0.95)
  
  Tukey_HSD_table = as.data.frame(Tukey_HSD$group)
  
  print(Tukey_HSD_table)
  #save the result from the tukey test
  #write.table(paste(date(), "\nGroup\t", groupID, "\n\t", sep=""), file=paste("../results/4.Diversity_ana/alpha_boxplot_TukeyHSD.txt",sep=""),append = T, quote = F, eol = "", row.names = F, col.names = F)
  #suppressWarnings(write.table(Tukey_HSD_table, file=paste("../results/4.Diversity_ana/alpha_boxplot_TukeyHSD.txt",sep=""), append = T, quote = F, sep="\t", eol = "\n", na = "NA", dec = ".", row.names = T, col.names = T))
  
  
  # 函数：将Tukey检验结果P值转换为显著字母分组
  # 输入文件为图基检验结果和分组
  generate_label_df = function(TUKEY, variable){
    library(multcompView)
    # 转换P值为字母分组
    ## 提取图基检验中分组子表的第4列P adjust值
    Tukey.levels = TUKEY[[variable]][,4]
    ## multcompLetters函数将两两p值转换为字母，data.frame并生成列名为Letters的数据框
    Tukey.labels = data.frame(multcompLetters(Tukey.levels,reversed = T)['Letters'])
    
    # 按分组名字母顺序
    ## 提取字母分组行名为group组名
    Tukey.labels$group = rownames(Tukey.labels)
    # 按组名的字母顺序排列，默认的Levels
    Tukey.labels=Tukey.labels[order(Tukey.labels$group), ]
    return(Tukey.labels)
  }
  
  # 当只有两组时，用LSD标注字母
  if (length(unique(df$group)) == 2){
    # LSD检验，添加差异组字母
    library(agricolae)
    out = LSD.test(model, "group", p.adj="none")
    stat = out$groups
    # 分组结果添入Index
    df$stat=stat[as.character(df$group),]$groups
    # 当大于两组时，用multcompView标注字母
  }else{
    library(multcompView)
    LABELS = generate_label_df(Tukey_HSD , "group")
    df$stat=LABELS[as.character(df$group),]$Letters
  }
  
  # 设置分组位置为各组y最大值+高的5%
  max=max(df[,c(index)])
  min=min(df[,index])
  x = df[,c("group",index)]
  y = x %>% group_by(group) %>% summarise_(Max=paste('max(',index,')',sep=""))
  y=as.data.frame(y)
  rownames(y)=y$group
  df$y=y[as.character(df$group),]$Max + (max-min)*0.05
  # head(df)
  
  # ================== calculate the sd and mean
  
  colnames(df) = c("alpha","group","stat","y")
  went = df %>% group_by(group) %>% summarize(mean = mean(alpha),SD = sd(alpha),label = unique(stat)) 
  #or
  #wen1 = as.data.frame(tapply(as.vector(as.matrix(df[1])), df$group, mean, na.rm=TRUE))
  #wen2 = as.data.frame(tapply(as.vector(as.matrix(df[1])), df$group, sd, na.rm=TRUE))
  #went = cbind(wen1, wen2)
  #colnames(went) = c("mean","SD")
  #aa = distinct(df, group, .keep_all = TRUE)
  #went$label = aa$stat[match(row.names(went),aa$group)]
  #went$Groups = row.names(went)
  
  a = max(went$mean + went$SD)*1.1
  
  # =========== boxplot
  p = ggplot(df, aes(x=group, y=alpha, color=group)) +
    geom_boxplot(alpha=1, outlier.shape = NA, outlier.size=0, size=0.7, width=0.5, fill="transparent")  + theme_classic() + scale_fill_gdocs() +
    geom_text(data=df, aes(x=group, y=y, color=group, label=stat)) + 
    labs(x="Groups", y=paste(ylab), color=groupID) +
    geom_jitter(position=position_jitter(0.17), size=1, alpha=0.7)+
    theme(text=element_text(family="sans", size=20)) + labs(y = ylab, x = "") + 
    scale_color_manual(values = c("navy","firebrick","gold3","grey")) +
    theme(legend.position = c(0.8,0.8))
  p
  return(p)
  
}


p = plot_alpha_boxplot(index = "shannon_entropy", groupID = "Species")
ggsave(p,filename = "../../UKBB_Microbiomes/results/7.Final_graph/shannon_species_boxplot.png")

p = plot_alpha_boxplot(index = "faith_pd", groupID = "Species")
ggsave(p,filename = "../../UKBB_Microbiomes/results/7.Final_graph/faith_pd_species_boxplot.png",width = 7, height = 4,dpi = "print")

p = plot_alpha_boxplot(index = "faith_pd", groupID = "CrithidiaInfection")
ggsave(p,filename = "../../UKBB_Microbiomes/results/7.Final_graph/faith_pd_Crithidia_boxplot.png",width = 7, height = 5,dpi = "print")

p = plot_alpha_boxplot(index = "faith_pd", groupID = "NosemaInfection")
ggsave(p,filename = "../../UKBB_Microbiomes/results/7.Final_graph/faith_pd_Nosema_boxplot.png",width = 7, height = 5,dpi = "print")

p = plot_alpha_boxplot(index = "faith_pd", groupID = "CollectionSite")
ggsave(p,filename = "../../UKBB_Microbiomes/results/7.Final_graph/faith_pd_collectionSite_boxplot.png",width = 7, height = 5,dpi = "print")




# =============== beta diversity

library("qiime2R")

metadata<-read.table("./sample-metadata.tsv",header = T) %>% rename(SampleID = `sample_name`) 
metadata <- read.csv("./metadata_Microbiome.csv",header = T) %>% mutate(SampleID = as.character(SampleID))

metadata <- df_meta %>% rownames_to_column(var = "SampleID")

uwunifrac<-read_qza("../results/4.Diversity_ana/core-metrics-results/unweighted_unifrac_pcoa_results.qza")$data$Vectors
wunifrac <-read_qza("../results/4.Diversity_ana/core-metrics-results/weighted_unifrac_pcoa_results.qza")$data$Vectors
bc_dissimilar <-read_qza("../results/4.Diversity_ana/core-metrics-results/bray_curtis_pcoa_results.qza")$data$Vectors
jaccard_dis <- read_qza("../results/4.Diversity_ana/core-metrics-results/jaccard_pcoa_results.qza")$data$Vectors

proportion_uwunifra<-read_qza("../results/4.Diversity_ana/core-metrics-results/unweighted_unifrac_pcoa_results.qza")$data$ProportionExplained
proportion_wunifra<-read_qza("../results/4.Diversity_ana/core-metrics-results/weighted_unifrac_pcoa_results.qza")$data$ProportionExplained
proportion_bc <- read_qza("../results/4.Diversity_ana/core-metrics-results/bray_curtis_pcoa_results.qza")$data$ProportionExplained
proportion_jac <-  read_qza("../results/4.Diversity_ana/core-metrics-results/jaccard_pcoa_results.qza")$data$ProportionExplained

plot_pcoA_with_alpha_beta <- function(ordination_result,proportion_result,groupID = "Enterotype" ){
  faith_pd_vector<-read_qza("../results/4.Diversity_ana/core-metrics-results/faith_pd_vector.qza")$data %>% rownames_to_column("SampleID") 
  ordination_result[,"SampleID"] <- as.character(ordination_result[,"SampleID"] )
  
  combin_matrix <- ordination_result%>%
    dplyr::select(SampleID, PC1, PC2) %>%
    left_join(metadata) %>%
    left_join(faith_pd_vector) 
  
  col_vector = str_replace(colnames(combin_matrix),groupID,replacement = "group")
  
  colnames(combin_matrix) = col_vector
  
  p = ggplot(data = combin_matrix, aes(x=PC1, y=PC2)) +
    geom_point(mapping = aes( color=group),alpha=0.5) + #alpha controls transparency and helps when points are overlapping
    theme_classic() + 
    #scale_shape_manual(values=c(16,1,6,5), name="CollectionSite") + #see http://www.sthda.com/sthda/RDoc/figure/graphs/r-plot-pch-symbols-points-in-r.png for numeric shape codes
    #scale_size_continuous(name="Phylogenetic Diversity") +
    scale_color_discrete(name=groupID) + scale_color_manual(values = c("navy","firebrick","gold3","grey")) + stat_ellipse(aes(color=group),linetype = 2,level = 0.68) + 
    labs(x=paste("PCo1 (", format(100 * proportion_result[1] , digits=4), "%)", sep=""),
         y=paste("PCo2 (", format(100 * proportion_result[2] , digits=4), "%)", sep=""),color = groupID) +
    geom_hline(yintercept = 0,size = 0.1) +
    geom_vline(xintercept = 0,size = 0.1) + 
    theme(legend.position = c(0.85,0.85),text=element_text(family="sans", size=20))
  p
  
  return(p)
}

p = plot_pcoA_with_alpha_beta(ordination_result=wunifrac,proportion_result=proportion_wunifra, groupID = "Species")
ggsave(p, filename = "../results/7.Final_graph/PCoA_beta_wei_unifra_dis_Species.png",width = 10, height = 9, dpi = "print",scale=0.7)

p = plot_pcoA_with_alpha_beta(ordination_result=wunifrac,proportion_result=proportion_wunifra, groupID = "CollectionSite")
ggsave(p, filename = "../results/7.Final_graph/PCoA_beta_wei_unifra_dis_CollectionSite.png",width = 10, height = 9, dpi = "print",scale=0.7)

p = plot_pcoA_with_alpha_beta(ordination_result=wunifrac,proportion_result=proportion_wunifra, groupID = "CrithidiaInfection")
ggsave(p, filename = "../results/7.Final_graph/PCoA_beta_wei_unifra_dis_Crithidia.png",width = 10, height = 9, dpi = "print",scale = 0.7)

p = plot_pcoA_with_alpha_beta(ordination_result=wunifrac,proportion_result=proportion_wunifra, groupID = "NosemaInfection")
ggsave(p, filename = "../results/7.Final_graph/PCoA_beta_wei_unifra_dis_Nosema.png",width = 10, height = 9, dpi = "print",scale = 0.7)




p = plot_pcoA_with_alpha_beta(ordination_result=uwunifrac,proportion_result=proportion_uwunifra, groupID = "Species" )
ggsave(p, filename = "../results/7.Final_graph/PCoA_beta_unwei_unifra_dis.pdf")


p = plot_pcoA_with_alpha_beta(bc_dissimilar,proportion_result = proportion_bc, groupID = "Species")
ggsave(p, filename = "../results/7.Final_graph/PCoA_beta_BC_dis.pdf")

p = plot_pcoA_with_alpha_beta(jaccard_dis,proportion_result = proportion_jac, groupID = "Species")
ggsave(p, filename = "../results/7.Final_graph/PCoA_beta_JAC_dis.pdf")



# ============= other method fro beta diversity


wunifra <-read_qza("../results/4.Diversity_ana/core-metrics-results/weighted_unifrac_distance_matrix.qza")
wunifra2 <-as.matrix(wunifra$data)

unwunifra <- read_qza("../results/4.Diversity_ana/core-metrics-results/unweighted_unifrac_distance_matrix.qza")
unwunifra2 <- as.matrix(unwunifra$data) # distance matrix

bc_distance <- read_qza("../results/4.Diversity_ana/core-metrics-results/bray_curtis_distance_matrix.qza")

jaccard_dis <- read_qza("../results/4.Diversity_ana/core-metrics-results/jaccard_distance_matrix.qza")


#b_c_dissimilar <- read_qza("../results/4.Diversity_ana/core-metrics-results/bray_curtis_distance_matrix.qza")
#b_c_dissimilar2  <- as.matrix(b_c_dissimilar$data) # distance matrix

beta_pcoa <- function(dis_mat, groupID="CollectionSite", ellipse=T, label=F, PCo=12) {
  # 依赖关系检测与安装
  p_list=c("ggplot2", "vegan", "ggrepel")
  for(p in p_list){
    if (!requireNamespace(p)){
      install.packages(p)}
    suppressWarnings(suppressMessages(library(p,character.only=T)))}
  
  # 测试默认参数
  # dis_mat=beta_unifrac
  # metadata=metadata
  # groupID="Group"
  # ellipse=T
  # label=F
  # PCo=12
  
  metadata <- read.table("./sample-metadata.tsv",header = T,row.names = 1)
  # 交叉筛选
  idx=rownames(metadata) %in% rownames(dis_mat)
  metadata=metadata[idx,,drop=F]
  
  # 提取样品组信息,默认为group可指定
  sampFile=as.data.frame(metadata[, groupID],row.names=row.names(metadata))
  # colnames(sampFile)[1]="group"
  
  # PCoA
  pcoa=cmdscale(dis_mat, k=3, eig=T) # k is dimension, 3 is recommended; eig is eigenvalues
  points=as.data.frame(pcoa$points) # get coordinate string, format to dataframme
  eig=pcoa$eig
  points=cbind(points, sampFile[rownames(points),])
  colnames(points)=c("x", "y", "z","group")
  
  # 按1、2轴绘图
  if (PCo == 12){
    p=ggplot(points, aes(x=x, y=y, color=group))  +
      labs(x=paste("PCo 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
           y=paste("PCo 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep=""), color=groupID)
  }
  # 按1、3轴绘图
  if (PCo == 13){
    p=ggplot(points, aes(x=x, y=z, color=group))  +
      labs(x=paste("PCo 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
           y=paste("PCo 3 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep=""), color=groupID)
  }
  # 按2、3轴绘图
  if (PCo == 23){
    p=ggplot(points, aes(x=y, y=z, color=group))  +
      labs(x=paste("PCo 2 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
           y=paste("PCo 3 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep=""), color=groupID)
  }
  p=p + geom_point(alpha=.7, size=2) + theme_classic() + theme(text=element_text(family="sans", size=7))
  # 是否添加置信椭圆
  if (ellipse == T){
    p=p + stat_ellipse(level=0.68)
  }
  # 是否显示样本标签
  if (label == T){
    p=p + geom_text_repel(label=paste(rownames(points)), colour="black", size=3)
  }
  p
}

beta_pcoa(dis_mat = unwunifra2)
#beta_pcoa(dis_mat = b_c_dissimilar2, groupID = "Species")

# ===================================== Proceed the feature table into alpha diversity ===========
# read the file from the alpha diversity result
#alpha_div <- read.table("../results/4.Alpha_result/alpha_diversity_vector/observed_features/alpha-diversity.tsv")


# =============== other significant visualization format
# Not necessary

#p_value_one <- data.frame(
#  group = c("Infected", "Infected", "Uninfected", "Uninfected"),
#  alpha = c(50, 60, 60, 50))
#p + geom_line(data = p_value_one, aes(x = group, y =alpha, group = 1)) +
#  annotate("text", x = 1.5, y = 60, label = "*",size = 8, color = "#22292F")



# ================= plotting dendrogram

plot_dendrogram <- function(dist_object = unwunifra$data){
  library("ape")
  upgma <- hclust(dist_object)
  
  # label the dendogram with 
  cluster <- as.data.frame(as.matrix(dist_object)) %>% merge(metadata,by.x = "row.names",by.y="SampleID") %>% 
    rename( SampleID = `Row.names`) %>%
    dplyr:: select("SampleID","Species")
  
  x = factor(cluster$Species)
  
  colors = c("navy","firebrick","gold3")
  test_setname <- setNames(as.numeric(x),cluster$SampleID)
  plot(as.phylo(upgma),cex = 0.8,type = "fan",tip.color = colors[test_setname],label.offset = 0.01,edge.color ="steelblue")
  legend("left",legend = levels(x),fill = colors,border = F)
}


pdf("../results/7.Final_graph/dendrogram_unwei_unifrac.pdf", width = 10)
plot_dendrogram(dist_object = unwunifra$data)
dev.off()

pdf("../results/7.Final_graph/dendrogram_wei_unifrac.pdf", width = 10)
plot_dendrogram(dist_object = wunifra$data)
dev.off()

pdf("../results/7.Final_graph/dendrogram_jac.pdf", width = 10)
plot_dendrogram(dist_object = jaccard_dis$data)
dev.off()

pdf("../results/7.Final_graph/dendrogram_bc.pdf", width = 10)
plot_dendrogram(dist_object = bc_distance$data)
dev.off()




#plot_dendrogram(dist_object = b_c_dissimilar$data)


# ================== Statistical analysis

# ======== alpha diversity (kruskal-wallis-pairwise-CollectionSite.csv)

# ==== Species and CollectionSites could be found in qiime 

# Hoever, the interaction between them could not be found

metadata <- read.table("./sample-metadata.tsv",header = T) %>% 
  mutate(Interaction = paste0(Species,CollectionSite)) %>% rename(SampleID = sample_name)
table(metadata$Interaction)

joined_table <- left_join(shannon_vector,metadata) %>% left_join(faith_pd_vector)
table(joined_table$Interaction)

kruskal.test(shannon_entropy~Interaction,data = joined_table)

# ======== beta diversity (permanova-pairwise.csv)
library("vegan")

a = as.dist(wunifra$data)
dispersion<-betadisper(a,group = joined_table$Interaction)
permutest(dispersion)







