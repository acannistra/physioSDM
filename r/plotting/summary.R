suppressPackageStartupMessages({
  library(dplyr)
  library(reshape2)
  library(ggplot2)
  require(gridExtra)
})

#args = commandArgs(trailingOnly = T)
results_dir = "../results/remote_results/706cb7aa-7891-11e7-b49f-0242ac110002/"
results_fname = "results.csv"
results = read.csv(paste0(results_dir, results_fname))

mean_graf = mean(results$graf_cur_auc)
mean_phys = mean(results$phys_cur_auc)
mean_me   = mean(results$maxent_cur_auc)
print(nrow(results))
#results = sample_frac(results[!duplicated(results$species),], size=0.44) ## why are there duplicate species??


melted = melt(results[,c("species", 'graf_cur_auc', 'phys_cur_auc', 'maxent_cur_auc')])




# ggplot(data = results_nodup[1:10,], mapping=aes(x=species, y='variable'), position=position_dodge()) + 
#   geom_col(mapping=aes(y=phys_cur_auc),  color='red', fill=NA) +
#   geom_col(mapping=aes(y=graf_cur_auc), color='blue', fill=NA) +
#   geom_col(mapping=aes(y=maxent_cur_auc), fill=NA, color='yellow')+
#   theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
#   
p1 = ggplot(melted, aes(x=species, y=value)) + 
  geom_col(aes(fill=variable), position='dodge') + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  labs(title="AUC performance comparision - Current Distribution Recovery")
  

means = melted %>% group_by(variable) %>% summarise(mean = mean(value), std = sd(value))


p2 = ggplot(melted, aes(x=variable)) + 
  geom_bar(data=means, stat='identity', aes(y=mean, fill=variable))  + 
  geom_errorbar(data=means, aes(ymin = mean-std, ymax = mean+std), width=0.2, position=position_dodge(0.9)) + 
  labs(title=sprintf("Algorithm Means (n = %d)", length(unique(melted$species)))) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) 

png(paste0(results_dir, "summaryPlot.png"), units = 'in', width=11, height=8.5, res=300)
grid.arrange(p1, p2)
dev.off()
                     