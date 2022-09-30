# ======================== Use vegan here
# ======================== alpha rarefraction plot

set.seed(100)
col <- sample(colors(),24)
# we need 24 color for our project
lty <- c("solid", "dashed", "longdash", "dotdash")
pars <- expand.grid(col = col, lty = lty, stringsAsFactors = FALSE)
head(pars)

out <- with(pars[1:ncol(feature_abundance_matrix), ],
            rarecurve(feature_abundance_matrix, step = 20, col = col,
                      lty = lty, label = FALSE,xlab = "Sampling depth per sample",ylab = "Observed feature numbers",MARGIN =2 ))


# ======= plot rarefied result
# Two plot function
source("./b2.h.alpha_analysis.R")

filter_matrix_5000 <- feature_abundance_matrix[,colSums(feature_abundance_matrix)>5000]

set.seed(100)
Srare <- as.matrix(rarefy(filter_matrix_5000, sample = 5000,MARGIN = 2))
colnames(Srare) <- "RarefiedNum"


# ======= plot unrarefied result


Grouplist = c( "CollectionSite","Species")
alpha_matric = c("richness","shannon","simpson","invsimpson","PieLou.evenness","Fisher_alpha")


for (matric in alpha_matric){
  for (groupid in Grouplist){
    pic_name = paste0("../results/4.Diversity_ana/plot/",groupid,"_",matric,".png")
    p = plot_alpha_boxplot(index = matric, groupID = groupid)
    ggsave(filename = pic_name,device = "png",width = 8,height = 6)
  }
}
