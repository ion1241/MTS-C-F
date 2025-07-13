#Bin abundance on TS
#Load Libraries####
library(tidyverse)
library(vegan)
library(wTO)
#library(Hmisc)
library(igraph)
#library(clusterProfiler) #Not avail in this version of R

#Load functions####
#FilterCoocurring Removes associations between taxa that do not occur > X% together (default 50).
## (i.e. 50% of the times that at least one of them occur)
## Accepts object from e.g. wTO w Results (The 2 first columns being node names)
## and an OTU table w taxa as rows


getwd()
setwd("../CoverM/")

#Load Files
dates <- read_csv2("Sample_date_readcount.csv")

# Define the suffix to be removed
suffix <- "output.tsv"



temp<- list.files(pattern = "*output.tsv")
#Leer todas las muestras de Oiartzun
for (i in 1:length(temp)) assign(temp[i], readr::read_tsv(temp[i],  name_repair = "check_unique", skip = 1))

# List all objects in the environment
all_objects <- ls()

# Identify dataframes whose names start with "S"
s_dataframes <- all_objects[grepl("*output.tsv", all_objects) & sapply(all_objects, function(x) is.data.frame(get(x)))]

# Store these dataframes into a list
s_dataframe_list <- lapply(s_dataframes, get)

# Name the list elements with the dataframe names
names(s_dataframe_list) <- s_dataframes


# Function to rename columns
rename_columns <- function(df, filename) {
  colnames(df) <- c("MAG", filename)
  return(df)
}

# Iterate over the list and rename columns
for (name in names(s_dataframe_list)) {
  s_dataframe_list[[name]] <- rename_columns(s_dataframe_list[[name]], name)
}

# Print the renamed dataframes
print(s_dataframe_list)
#retrieve the dataframes of s_dataframe_list to the environment
list2env(s_dataframe_list, envir = .GlobalEnv)

#Juntar todas por MAG
#from s_dataframe_list merge all dataframes by the first column
MAG_TS_All <- Reduce(function(x, y) merge(x, y, by = "MAG", all = TRUE), s_dataframe_list)

#change colnames from MAG_TS_All that appear in dates$SerieNumber for the ones in dates$'Date Info'
MAG_TS_All<-as.data.frame(MAG_TS_All)

rownames(MAG_TS_All) <- MAG_TS_All$MAG
MAG_TS_All <- MAG_TS_All[,-1]
#Order colnames of MAG_TS_All
MAG_TS_All <- MAG_TS_All[,order(colnames())]
names(MAG_TS_All) %in% dates$SerieNumber

#Remove from the colnames the output.tsv part
colnames(MAG_TS_All) <- gsub(suffix, "", colnames(MAG_TS_All))

MAG_TS_Allnames = MAG_TS_All[,dates$SerieNumber]
table(names(MAG_TS_Allnames) == dates$SerieNumber)
colnames(MAG_TS_Allnames) = dates$'Date Info'
summary(MAG_TS_Allnames)

#Make rownames the first column, called Genome
MAG_TS_Allnames$Genome <- rownames(MAG_TS_Allnames)
#Make a DF with the median, mean max and min for each column of MAG_TS_Allnames
MAG_TS_Allnames_stats <- data.frame(
  mean = colMeans(MAG_TS_Allnames[, -ncol(MAG_TS_Allnames)], na.rm = TRUE),
  median = apply(MAG_TS_Allnames[, -ncol(MAG_TS_Allnames)], 2, median, na.rm = TRUE),
  max = apply(MAG_TS_Allnames[, -ncol(MAG_TS_Allnames)], 2, max, na.rm = TRUE),
  min = apply(MAG_TS_Allnames[, -ncol(MAG_TS_Allnames)], 2, min, na.rm = TRUE)
)

#Join Dates with MAG_TS_Allnames
MAG_TS_Allnames_stats <- merge(MAG_TS_Allnames_stats, dates, by.x = "row.names", by.y = "Date Info")
#From MAG_TS_Allnames_stats Drop columns 6 and 7, rename the first column to Sample
MAG_TS_Allnames_stats <- MAG_TS_Allnames_stats[, -c(6, 7)]
colnames(MAG_TS_Allnames_stats)[1] <- "Sample"
write.table(MAG_TS_Allnames_stats, "MAG_TS_Allnames_stats.tsv", col.names = T, row.names = F, sep = "\t", quote=F, dec = ",")

#Filter so only the rows with the same  Genome as user_genome in Taxonomy_clust_compTop are retained
Taxonomy_clust_compTop<-read_tsv("../../FunNet/ORF_mapping/Taxonomy_clust_compTop.tsv")
MAG_TS_Allnames<-MAG_TS_Allnames[MAG_TS_Allnames$Genome %in% Taxonomy_clust_compTop$user_genome,]

#Filter MAG_TS_All a solo columnas que contengan EBI20*
MAG_TS_EBI20<-MAG_TS_Allnames[,grepl("EBI20",colnames(MAG_TS_Allnames))]
write.table(MAG_TS_EBI20, "MAG_TS_EBI20.tsv", col.names = NA, row.names = TRUE, sep = "\t", quote=F)
MAG_TS_EOI15<-MAG_TS_Allnames[,grepl("EOI15",colnames(MAG_TS_Allnames))]
write.table(MAG_TS_EOI15, "MAG_TS_EOI15.tsv", col.names = NA, row.names = TRUE, sep = "\t", quote=F)
colSums(MAG_TS_EBI20)


####BIDASOA####
#normalizar MAG_TS_EBI20 usando vegan decostand normalize
MAG_TS_EBI20_norm<-decostand(MAG_TS_EBI20, method = "total", MARGIN = 2)
colSums(MAG_TS_EBI20_norm)

BidasoaMAG.pa <- decostand(MAG_TS_EBI20_norm, method = "pa")
# 50% prevalence filter
MAG_TS_EBI20_prev <- MAG_TS_EBI20_norm[rowSums(BidasoaMAG.pa)>=10,]
head(MAG_TS_EBI20_prev)
write.csv(MAG_TS_EBI20_prev, file = "BidasoaTS_MAGs_prev.csv", dec = ".")


####OIARTZUN####
#normalizar MAG_TS_EOI15 usando vegan decostand normalize
MAG_TS_EOI15_norm<-decostand(MAG_TS_EOI15, method = "total", MARGIN = 2)
colSums(MAG_TS_EOI15_norm)

OiartzunMAG.pa <- decostand(MAG_TS_EOI15_norm, method = "pa")
# 50% prevalence filter
MAG_TS_EOI15_prev <- MAG_TS_EOI15_norm[rowSums(OiartzunMAG.pa)>=10,]
head(MAG_TS_EOI15_prev)
write.csv(MAG_TS_EOI15_prev, file = "OiartzunTS_MAG_prev.csv", dec = ".")

