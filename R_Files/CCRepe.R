# Script by Leire Garate.
# This script is used to generate networks with ccrepe, can be modified by using diferent inputs to generate different networks.


library(tidyverse)
library(vegan)
library(ccrepe)


es.data <- PrevFilteEstuaryData
es.data <- column_to_rownames(es.data,var="...1")
es.data.t <- as.data.frame(t(es.data))
es.data.t.rowsum <- apply(es.data.t,1,sum)
es.data.norm <- es.data.t/es.data.t.rowsum
apply(es.data.norm,1,sum)

es.data.t.ccrp <- ccrepe(x=es.data.norm, iterations=250)

es.data.ccrp.filt <- es.data.t.ccrp$q.values %>% as.data.frame %>% tibble::rownames_to_column() %>% 
  tidyr::pivot_longer(-rowname)
es.data.ccrp.filt <- na.omit(es.data.ccrp.filt)

es.data.ccrp.sim <- es.data.t.ccrp$sim.score %>% as.data.frame %>% tibble::rownames_to_column() %>% 
  tidyr::pivot_longer(-rowname)
es.data.ccrp.sim <- na.omit(es.data.ccrp.sim)
colnames(es.data.ccrp.sim)[3] <- "sim.score"

es.data.mrg <- merge(es.data.ccrp.filt, es.data.ccrp.sim, by.x = c("rowname", "name"), by.y = c("rowname", "name"))
colnames(es.data.mrg)[1] <- "node1"
colnames(es.data.mrg)[2] <- "node2"



# Filter the results to keep only those with q-value < 0.001 or 0.0001
es.data.mrg.0.001 <- es.data.mrg[es.data.mrg$value < 0.001,]
#es.data.mrg.0.0001 <- es.data.mrg[es.data.mrg$value < 0.0001,]
es.data.mrg.0.001 <- es.data.mrg.0.001 %>% distinct(value, sim.score, .keep_all = TRUE)

#Add a column of abs value of sim.score to ebiMCr
es.data.mrg.0.001$score_abs <- abs(es.data.mrg.0.001$sim.score)



#If need be, filter to get the percentile
es.data.mrg.0.001Sorted <- es.data.mrg.0.001 %>% arrange(desc(score_abs))
es.data.mrg.0.001Edges <- nrow(es.data.mrg.0.001)
es.data.mrgTopXEdges <- ceiling(0.01 * es.data.mrg.0.001Edges)

# Select the 1% 
es.data.mrgTop1 <- es.data.mrg.0.001Sorted[1:es.data.mrgTopXEdges, ]
write_tsv(es.data.mrgTop1, "es.data.mrgTop1perc.txt")



# Write the results to a file change estuary and method to Bidasoa or Oiartzun and MAGs, 16S or KOs
write.table(es.data.mrg.0.001, file = "estuary_method_ccrp.txt", append = FALSE, quote = F, sep = "\t", dec = ".", row.names = F, col.names = T)


