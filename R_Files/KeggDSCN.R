
library(tidyverse)
library(vegan)
library(igraph)


#Sampling Data and readcount number
SampleData<-read.csv(file = "Sample_date_readcount.csv", sep = ";")
BidasoaData<-SampleData[c(1:20),]
OiartzunData<-SampleData[-c(1:20),]
ORFInfo<-read_delim("ProkORF_Tax.txt")


keggIDs<-read_csv2("ko_tree.csv")

keggIDs<-droplevels(subset(keggIDs, catA_name != "Human Diseases"))
keggIDs<-droplevels(subset(keggIDs, catA_name != "Genetic Information Processing"))
keggIDs<-droplevels(subset(keggIDs, catA_name != "Human Diseases"))
keggIDs<-droplevels(subset(keggIDs, catA_name != "Not Included in Pathway or Brite"))
keggIDs<-droplevels(subset(keggIDs, catA_name != "Not Included in Pathway or Brite"))



#Global Ocean Metatranscriptome single copy genes

Singlecopygenes<- list("K06942", "K01889", "K01887", "K01875", "K01883", "K01869", "K01873", "K01409", "K03106", "K03110")


#Oiartzun TS Data is in reads per kilobsae Change to TPM#### Change dierectory to Oiartzun mapping files####
setwd("Oiartzun/")
temp<- list.files(pattern = "*_clean.txt")
for (i in 1:length(temp)) assign(temp[i], readr::read_tsv(temp[i],  name_repair = "check_unique"))

OiartzunTS_All <- Reduce(function(x, y) merge(x, y, by = "#ID", all = TRUE), mget(ls(pattern = "*_clean.txt"))) #Numero de  115.982.915 ORFs

rm(list = grep("*clean.txt", ls(), value = TRUE))

gc()
write_tsv(OiartzunTS_All, "Oiartzun/Oiartzun_AllTS.tsv")

#OiartzunTS_All <- read_tsv("Oiartzun/Oiartzun_AllTS.tsv")

#Order by colname and get TPM:
#We have reads/lenght of ORF in bp
OiartzunTS_All<-as.data.frame(OiartzunTS_All)

rownames(OiartzunTS_All)<-OiartzunTS_All$`#ID`
OiartzunTS_All<-OiartzunTS_All[,-1]

#Set Sampling dates as Colnames instead of SXX codes
match(OiartzunData[,"SerieNumber"], names((OiartzunTS_All)))
names(OiartzunTS_All)[match(OiartzunData[,"SerieNumber"], names((OiartzunTS_All)))] = OiartzunData[,"Date.Info"]


#Normalizar a TPM
OiartzunTS_All[is.na(OiartzunTS_All)] <- 0
OiartzunTS_All<-OiartzunTS_All[,order(names(OiartzunTS_All))] #Ordenar columnas
OiartzunTS_All<-OiartzunTS_All/1000#Pasar de reads/lenght of gene  a reads/kb lenght gene = RPK
OiartzunScaling<-as.data.frame(colSums(OiartzunTS_All)/1000000) #Scaling factor
OiartzunScalingW<-as.data.frame(t(OiartzunScaling))
OiartzunScalingW<-OiartzunScalingW[,order(names(OiartzunScalingW))] #Ordenar columnas


OiartzunTPM <- OiartzunTS_All / OiartzunScalingW[rep(1, nrow(OiartzunTS_All)), ]

colSums(OiartzunTPM) #Tiene que dar 10⁶
OiartzunTPM$ID<-rownames(OiartzunTPM) #Uncoment this only for writing table
head(OiartzunTPM)
write_tsv(OiartzunTPM, "Oiartzun/Oiartzun_AllTS_TPM.tsv")


#Añadir los datos de KEGG a la tabla
head(OiartzunTPM)
colnames(OiartzunTPM)
colnames(ORFInfo)
OiartzunTSData<-inner_join(OiartzunTPM, ORFInfo, by = c("ID"="ORF ID")) # 36M de ORFs con info taxonomica de prok
#Checks
colnames(OiartzunTSData)


##Normalizar por single copy genes####
#Eliminar Info de taxonomia y ORFdata
OiartzunTSProk<-OiartzunTSData[,-c(22,23,24,25)] #Se borra to salvo valores aqui 

#Group by KEGG IDs
OiartzunTSKegg <- OiartzunTSProk %>%
  group_by(KEGG) %>%
  summarise_if(is.numeric, sum, na.rm = TRUE)
colnames(OiartzunTSKegg)
colSums(OiartzunTSKegg[,-1])
write_tsv(OiartzunTSKegg, "Oiartzun/Oiartzun_AllTS_TPM_KD.tsv")
OiartzunTSKegg <- read_tsv("Oiartzun/Oiartzun_AllTS_TPM_KD.tsv")



#Normalize by singlecopygene KEGGs
#Filter to get only the single copy genes
singlecopy_df <- OiartzunTSKegg %>% filter(KEGG %in% Singlecopygenes)

#Calculate the scaling factor for each column
scaling_factors <- colMeans(singlecopy_df[,-1], na.rm = TRUE)

# Normalize a new table
OiartzunTSKeggSCN <- OiartzunTSKegg

for (col in names(OiartzunTSKeggSCN)[-1]) {
  OiartzunTSKeggSCN[[col]] <- OiartzunTSKeggSCN[[col]] / scaling_factors[col]
}

#Set the single copy genes to 1 TPM
for (gene in Singlecopygenes) {
  OiartzunTSKeggSCN[OiartzunTSKeggSCN$KEGG == gene, -1] <- 1
}



write_tsv(OiartzunTSKeggSCN, "Oiartzun/OiartzunTSKEGGD_SCN.tsv")
OiartzunTSKeggSCN <- as.data.frame(read_tsv("Oiartzun/OiartzunTSKEGGD_SCN.tsv"))
row.names(OiartzunTSKeggSCN) <- OiartzunTSKeggSCN$KEGG
OiartzunTSKeggSCN <- OiartzunTSKeggSCN[,-1]


#Filter less than 1% fold of scg
OiartzunTSKeggSCN[OiartzunTSKeggSCN < 0.01 & OiartzunTSKeggSCN != 0.01] <- 0


#Filter by prevalence

Oiartzun.pa <- decostand(OiartzunTSKeggSCN, method = "pa")
#50% prevalence filter
OiartzunTS_KD_SCNprev <- OiartzunTSKeggSCN[rowSums(Oiartzun.pa)>=10,]
head(OiartzunTS_KD_SCNprev)
colSums(OiartzunTS_KD_SCNprev[,-1])



#Order the columns by date to avoid problems with network software
Ocol_names <- colnames(OiartzunTS_KD_SCNprev)

# Extract the date part from the column names
Odates <- sub("EOI15_S_", "", Ocol_names)

# Convert the extracted dates to Date objects
Odate_objects <- as.Date(Odates, format = "%d.%m.%Y")

# Order the columns based on the Date objects
Oordered_indices <- order(Odate_objects)

# Reorder the dataframe columns
OiartzunTS_KD_SCNprev <- OiartzunTS_KD_SCNprev[, Oordered_indices]

# Verify the column names are ordered
colnames(OiartzunTS_KD_SCNprev)

write.table(OiartzunTS_KD_SCNprev, "OiartzunTS_KeggD_SCN.txt",row.names = T)




#Bidasoa TS Data is in reads per kilobsae Change to TPM#### Change dierectory to Bidasoa mapping files####
setwd("Bidasoa/")
temp<- list.files(pattern = "*_clean.txt")
for (i in 1:length(temp)) assign(temp[i], readr::read_tsv(temp[i],  name_repair = "check_unique"))

BidasoaTS_All <- Reduce(function(x, y) merge(x, y, by = "#ID", all = TRUE), mget(ls(pattern = "*_clean.txt"))) #Numero de  115.982.915 ORFs

rm(list = grep("*clean.txt", ls(), value = TRUE))

gc()
write_tsv(BidasoaTS_All, "Bidasoa/Bidasoa_AllTS.tsv")

#BidasoaTS_All <- read_tsv("Bidasoa/Bidasoa_AllTS.tsv")

#Order by colname and get TPM:
#We have reads/lenght of ORF in bp
BidasoaTS_All<-as.data.frame(BidasoaTS_All)

rownames(BidasoaTS_All)<-BidasoaTS_All$`#ID`
BidasoaTS_All<-BidasoaTS_All[,-1]

#Set Sampling dates as Colnames instead of SXX codes
match(BidasoaData[,"SerieNumber"], names((BidasoaTS_All)))
names(BidasoaTS_All)[match(BidasoaData[,"SerieNumber"], names((BidasoaTS_All)))] = BidasoaData[,"Date.Info"]


#Normalizar a TPM
BidasoaTS_All[is.na(BidasoaTS_All)] <- 0
BidasoaTS_All<-BidasoaTS_All[,order(names(BidasoaTS_All))] #Ordenar columnas
BidasoaTS_All<-BidasoaTS_All/1000#Pasar de reads/lenght of gene  a reads/kb lenght gene = RPK
BidasoaScaling<-as.data.frame(colSums(BidasoaTS_All)/1000000) #Scaling factor
BidasoaScalingW<-as.data.frame(t(BidasoaScaling))
BidasoaScalingW<-BidasoaScalingW[,order(names(BidasoaScalingW))] #Ordenar columnas


BidasoaTPM <- BidasoaTS_All / BidasoaScalingW[rep(1, nrow(BidasoaTS_All)), ]

colSums(BidasoaTPM) #Tiene que dar 10⁶
BidasoaTPM$ID<-rownames(BidasoaTPM) #Uncoment this only for writing table
head(BidasoaTPM)
write_tsv(BidasoaTPM, "Bidasoa/Bidasoa_AllTS_TPM.tsv")


#Añadir los datos de KEGG a la tabla
head(BidasoaTPM)
colnames(BidasoaTPM)
colnames(ORFInfo)
BidasoaTSData<-inner_join(BidasoaTPM, ORFInfo, by = c("ID"="ORF ID")) # 36M de ORFs con info taxonomica de prok
#Checks
colnames(BidasoaTSData)


##Normalizar por single copy genes####
#Eliminar Info de taxonomia y ORFdata
BidasoaTSProk<-BidasoaTSData[,-c(22,23,24,25)] #Se borra to salvo valores aqui 

#Group by KEGG IDs
BidasoaTSKegg <- BidasoaTSProk %>%
  group_by(KEGG) %>%
  summarise_if(is.numeric, sum, na.rm = TRUE)
colnames(BidasoaTSKegg)
colSums(BidasoaTSKegg[,-1])
write_tsv(BidasoaTSKegg, "Bidasoa/Bidasoa_AllTS_TPM_KD.tsv")
BidasoaTSKegg <- read_tsv("Bidasoa/Bidasoa_AllTS_TPM_KD.tsv")



#Normalize by singlecopygene KEGGs
#Filter to get only the single copy genes
singlecopy_df <- BidasoaTSKegg %>% filter(KEGG %in% Singlecopygenes)

#Calculate the scaling factor for each column
scaling_factors <- colMeans(singlecopy_df[,-1], na.rm = TRUE)

# Normalize a new table
BidasoaTSKeggSCN <- BidasoaTSKegg

for (col in names(BidasoaTSKeggSCN)[-1]) {
  BidasoaTSKeggSCN[[col]] <- BidasoaTSKeggSCN[[col]] / scaling_factors[col]
}

#Set the single copy genes to 1 TPM
for (gene in Singlecopygenes) {
  BidasoaTSKeggSCN[BidasoaTSKeggSCN$KEGG == gene, -1] <- 1
}



write_tsv(BidasoaTSKeggSCN, "Bidasoa/BidasoaTSKEGGD_SCN.tsv")
BidasoaTSKeggSCN <- as.data.frame(read_tsv("Bidasoa/BidasoaTSKEGGD_SCN.tsv"))
row.names(BidasoaTSKeggSCN) <- BidasoaTSKeggSCN$KEGG
BidasoaTSKeggSCN <- BidasoaTSKeggSCN[,-1]


#Filter less than 1% fold of scg
BidasoaTSKeggSCN[BidasoaTSKeggSCN < 0.01 & BidasoaTSKeggSCN != 0.01] <- 0


#Filter by prevalence

Bidasoa.pa <- decostand(BidasoaTSKeggSCN, method = "pa")
#50% prevalence filter
BidasoaTS_KD_SCNprev <- BidasoaTSKeggSCN[rowSums(Bidasoa.pa)>=10,]
head(BidasoaTS_KD_SCNprev)
colSums(BidasoaTS_KD_SCNprev[,-1])



#Order the columns by date to avoid problems with network software
Ocol_names <- colnames(BidasoaTS_KD_SCNprev)

# Extract the date part from the column names
Odates <- sub("EOI15_S_", "", Ocol_names)

# Convert the extracted dates to Date objects
Odate_objects <- as.Date(Odates, format = "%d.%m.%Y")

# Order the columns based on the Date objects
Oordered_indices <- order(Odate_objects)

# Reorder the dataframe columns
BidasoaTS_KD_SCNprev <- BidasoaTS_KD_SCNprev[, Oordered_indices]

# Verify the column names are ordered
colnames(BidasoaTS_KD_SCNprev)

write.table(BidasoaTS_KD_SCNprev, "BidasoaTS_KeggD_SCN.txt",row.names = T)
