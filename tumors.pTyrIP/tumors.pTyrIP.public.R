#_______________________________________________________________________________________________________________________
# 20.01.2022
# 
# Project: proteomics.truncated.FGFR2.is.oncogene.cancer.Zingg.et.al 
# 
# Script name: tumors.pTyrIP.public
#
# Purpose of script: analysis protein phosphorylation data tumors truncated FGFR2 is oncogene cancer project Zingg et al 
#
# Author: Frank Rolfs 
#
# Notes:
#
# 
#_______________________________________________________________________________________________________________________


### set working directory
try(setwd(dirname(rstudioapi::getActiveDocumentContext()$path)))
getwd()

### load packages
library(tidyverse)
library(data.table)
library(pbapply)
library(ComplexHeatmap)
library(circlize)
library(viridis)
library(reshape2)
library(splitstackshape)
library(readxl)
library(scico)
library(forcats)
library(limma)
library(cmapR)
library(gtools)
library(paletteer)

###################################################################################################################################
##############################################################################################################################################################
#####################################################################################################################################################################
###### define Perseus expand function ###############################################################################################################################


### define Perseus-like expand function
PerseusPSexpand.17.02.21 <- function(DF1){
  print("Wait a moment... ...")
  
  #get all column names data
  allDFcolNAMES <- colnames(DF1) 
  
  # get column names with Intensities for 1x, 2x, 3x phosphorylated candidates
  allINT123names <- unlist(str_extract_all(string=allDFcolNAMES, pattern= "Intensity.+___\\d" )) 
  
  # prepare names for columns after expansion
  eachsamplesIntname <- unique(str_replace_all(allINT123names, pattern="___\\d", replacement="")) 
  
  # make DF with Intensity data only
  DFonlyINT123data <- DF1 %>% select(all_of(allINT123names)) 
  
  #count columns in DFonlyINT123data
  NCOLS_DFonlyINT123data <- ncol(DFonlyINT123data) 
  
  #make temporary DF with after expansion column names and row length of original DF
  tempDF1 <-matrix(NA, nrow=nrow(DF1), ncol=length(eachsamplesIntname))
  tempDF1 <- as.data.frame(tempDF1)
  colnames(tempDF1) <- eachsamplesIntname
  
  
  # concatenate data for __1, __2, __3 per sample (3 cols each) and write results in temporary DF prepared before
  z=0
  for(i in seq(1,ncol(DFonlyINT123data), by=3 ) ){
    #print(i)
    z=z+1
    
    everythirdCol <- seq(i,i+2, by=1) 
    tempDATAconcatenated <- DFonlyINT123data %>% unite(collapsed.result, all_of(everythirdCol), sep = ';') %>% select(collapsed.result)  %>% pull()
    
    tempDF1[, z] <- tempDATAconcatenated 
  }
  
  #combine original DF with temporary DF; leave out __1, __2, __3 columns
  DF1expanded <- cbind(DF1, tempDF1)
  DF1expanded <- DF1expanded %>% select(-contains("___"))
  
  #split concatenated DF & add multiplicity information
  DF1expanded <- cSplit(DF1expanded,                                  
                        splitCols    = eachsamplesIntname, 
                        sep          = ";", 
                        direction    = "long", 
                        type.convert = TRUE)
  
  DF1expanded <- as.data.frame(DF1expanded)
  
  DF1expanded$PS_Multiplicity <- c("_1", "_2", "_3")
  
  print("... ... done.")
  return(DF1expanded)
}


#####################################################################################################################################################################
### load data tumors IMAC DDA from MaxQuant Phospho(STY)Sites.txt

# load tumors pTyrIP sample labels
load("DZ.IP.tum.labels.Rdata")
glimpse(DZ.IP.tum.labels)


## phosphodata MQ _ Version FR _ no GFP and no FGFR2 variant FASTA  _ newer 2021/101 mouse FASTA file, MQ v2.0.3.0
#IP.tum.DZ <- fread("/Users/frankrolfs/NLPostDR/VUmc OPL/Phosphoproteomics with Daniel Zingg NKI/mouse tumors/pTyrIP.mouse.tumors.DZ/QE3_210824_OPL1025_FR_pTyrIP_FGFR2_mouse_tumor_byFR_noGFP_noVariants 01.12.2021/Phospho (STY)Sites.txt",integer64 = "numeric") #

IP.tum.DZ <- fread("tumors.pTyrIP.Phospho(STY)Sites.txt",integer64 = "numeric") #
glimpse(IP.tum.DZ)

nrow(IP.tum.DZ) #5382
ncol(IP.tum.DZ) #76
colnames(IP.tum.DZ)
length(unique(IP.tum.DZ$id)) #5382
glimpse(IP.tum.DZ)


### select columns of interest
IP.tum.DZ.2 <- IP.tum.DZ %>% select(
  
  "id",                          "Proteins",                    "Positions within proteins",   "Leading proteins",            "Protein",                    
  "Protein names",               "Gene names",                  "Fasta headers",               "Localization prob",           "Number of Phospho (STY)",    
  "Amino acid",                  "Sequence window",             "Phospho (STY) Probabilities", "Position in peptide",         "Reverse",                    
  "Potential contaminant",       "Positions",                   "Position",                    "Peptide IDs",                 "Mod. peptide IDs",      
  
  matches("___[0-9]"), -Intensity___1, -Intensity___2, -Intensity___3
  
)

glimpse(IP.tum.DZ.2)


### remove spaces from column names
column.names <- colnames(IP.tum.DZ.2)
column.names.changed <- unname(sapply(column.names , function(x) str_replace_all(string=x, pattern=" ", replacement=".")))
column.names.changed
colnames(IP.tum.DZ.2) <- column.names.changed

glimpse(IP.tum.DZ.2)
colnames(IP.tum.DZ.2)


#change column classes using mutate_at
IP.tum.DZ.2 <- mutate_at(IP.tum.DZ.2,  21:ncol(IP.tum.DZ.2), list(as.numeric) )
glimpse(IP.tum.DZ.2)


### remove reverse hits
table(IP.tum.DZ.2$Reverse) #52+
IP.tum.DZ.3 <- IP.tum.DZ.2[IP.tum.DZ.2$Reverse =="",]
nrow(IP.tum.DZ.3) #5330


### remove contaminants
table(IP.tum.DZ.3$Potential.contaminant) #67+
IP.tum.DZ.3 <- IP.tum.DZ.3[IP.tum.DZ.3$Potential.contaminant == "",]
nrow(IP.tum.DZ.3) #5263


### remove rows completely zero
ncol(IP.tum.DZ.3) #95
IP.tum.DZ.3$zero.check <- pbapply(IP.tum.DZ.3[, c(21:95)], 1, function(x) sum(x) ) #
table(IP.tum.DZ.3$zero.check >0) #F799 T4464

IP.tum.DZ.3.complete.zero.rows <- filter(IP.tum.DZ.3, zero.check <=0)
glimpse(IP.tum.DZ.3.complete.zero.rows) #

IP.tum.DZ.3 <- filter(IP.tum.DZ.3, zero.check > 0)
nrow(IP.tum.DZ.3) #4464


length(unique(IP.tum.DZ.3$id)) ##4464 = total number PS
table(IP.tum.DZ.3$Amino.acid) #S513 T380 Y3571 
length(unique(IP.tum.DZ.3$Sequence.window)) #4463

### remove column zero.check
IP.tum.DZ.3 <- IP.tum.DZ.3 %>% select(-zero.check)
glimpse(IP.tum.DZ.3)


### create variable with short sample names with order as in Phospho(STY)Sites.txt
samples.in.experiment <- colnames(IP.tum.DZ.3)
samples.in.experiment <- samples.in.experiment[21:ncol(IP.tum.DZ.3)]
samples.in.experiment <- str_remove(string=samples.in.experiment, pattern="Intensity.")
samples.in.experiment <- str_remove(string=samples.in.experiment, pattern="___[0-9]")
samples.in.experiment <- unique(samples.in.experiment)
samples.in.experiment

### loop over data to count number of PS per sample ###############################################################################
colnames(IP.tum.DZ.3)
ncol(IP.tum.DZ.3) #95

result.PS.count.IP.tum.DZ.3 <- c()
sample.name.count <- 1

for(i.1 in seq(from=21, to=95, by=3)){ 
  print(i.1)
  
  i.2 <- i.1+1
  i.3 <- i.2+1
  
  sample.name <- samples.in.experiment[sample.name.count]
  print(sample.name)
  col.names <- colnames(IP.tum.DZ.3)
  col.names.select <- col.names[c(i.1, i.2, i.3)]
  
  temp <- select(IP.tum.DZ.3,Amino.acid, col.names.select)
  temp$sum_1_2_3 <- pbapply(temp[, 2:4], 1, function(x) sum(x))
  
  numberPS.A <- filter(temp, sum_1_2_3 > 0)
  numberPS <- nrow(numberPS.A)
  
  nMVs <- filter(temp, sum_1_2_3 <= 0)
  nMVs <- nrow(nMVs)
  
  percentageSTY <- table(numberPS.A$Amino.acid)
  
  temp.result <- tibble(sample=sample.name, numberPS=numberPS, nMVs=nMVs, nS=percentageSTY[1], nT=percentageSTY[2], nY=percentageSTY[3])
  
  result.PS.count.IP.tum.DZ.3  <- bind_rows(result.PS.count.IP.tum.DZ.3, temp.result)
  
  sample.name.count <- sample.name.count +1
  
}
result.PS.count.IP.tum.DZ.3

#change result data frame
result.PS.count.IP.tum.DZ.3$sample2 <- as.numeric(result.PS.count.IP.tum.DZ.3$sample)
result.PS.count.IP.tum.DZ.3 <- result.PS.count.IP.tum.DZ.3 %>% arrange(sample2)
glimpse(result.PS.count.IP.tum.DZ.3)

mean(result.PS.count.IP.tum.DZ.3$numberPS) #1343.68

#merge results counting with additional sample label information
glimpse(DZ.IP.tum.labels)

result.PS.count.IP.tum.DZ.3 <- merge(result.PS.count.IP.tum.DZ.3, DZ.IP.tum.labels, by.x = "sample2", by.y = "sample", all.x = T, all.y = T, sort = F)
result.PS.count.IP.tum.DZ.3$variant.group.3 <- paste0(result.PS.count.IP.tum.DZ.3$sample2, "_", result.PS.count.IP.tum.DZ.3$variant.group.2)
glimpse(result.PS.count.IP.tum.DZ.3)


#define sample order
IP.tum.DZ.3.samples.order <- c(
  "2_Fgfr2",
  
  "14_Ate1",
  
  "7_Bicc1",
  "9_Bicc1",
  
  "17_Tacc2",
  "19_Tacc2",
  
  "4_dE18",
  "5_dE18",
  "6_dE18",
  
  "16_dE18-Ate1",
  
  "10_dE18-Bicc1",
  "11_dE18-Bicc1",
  "12_dE18-Bicc1",
  
  "20_dE18-Tacc2",
  "22_dE18-Tacc2",
  
  "23_dE18-IGR1",
  "24_dE18-IGR1",
  
  "25_dE18-IGR2",
  "26_dE18-IGR2",
  
  "3_E18-C2",
  "27_E18-C2",
  "28_E18-C2",
  
  "29_E18-C3",
  "30_E18-C3",
  
  "32_E18-C4" )
IP.tum.DZ.3.samples.order

# calculate percentage MV PS all classes level
nrow(IP.tum.DZ.3)  #464 = total number PS

result.PS.count.IP.tum.DZ.3$percentageMV <- pbsapply(result.PS.count.IP.tum.DZ.3$nMVs, function(x) (100*x)/4464)
result.PS.count.IP.tum.DZ.3

# calculate percentage pS/pT/pY PS all classes level
for(i in 1:nrow(result.PS.count.IP.tum.DZ.3)){
  print(i)
  
  total.PS.sample <- result.PS.count.IP.tum.DZ.3$numberPS[i]
  n.pS <- result.PS.count.IP.tum.DZ.3$nS[i]
  n.pT <- result.PS.count.IP.tum.DZ.3$nT[i]
  n.pY <- result.PS.count.IP.tum.DZ.3$nY[i]
  
  result.PS.count.IP.tum.DZ.3$percentage.pS[i] <- (100*n.pS)/total.PS.sample
  result.PS.count.IP.tum.DZ.3$percentage.pT[i] <- (100*n.pT)/total.PS.sample
  result.PS.count.IP.tum.DZ.3$percentage.pY[i] <- (100*n.pY)/total.PS.sample
}

glimpse(result.PS.count.IP.tum.DZ.3)


#reorder samples
result.PS.count.IP.tum.DZ.3 <- result.PS.count.IP.tum.DZ.3[match(IP.tum.DZ.3.samples.order, result.PS.count.IP.tum.DZ.3$variant.group.3 ),]
glimpse(result.PS.count.IP.tum.DZ.3)

# calcuate means
mean.PS.allclassPS.per.sample <- mean(result.PS.count.IP.tum.DZ.3$numberPS); mean.PS.allclassPS.per.sample #1343.68

mean.PS.allclassPS.per.sample.pS <- mean(result.PS.count.IP.tum.DZ.3$percentage.pS);round(mean.PS.allclassPS.per.sample.pS, 1) #%pS 5.1
mean.PS.allclassPS.per.sample.pT <- mean(result.PS.count.IP.tum.DZ.3$percentage.pT);round(mean.PS.allclassPS.per.sample.pT, 1) #%pT 4.7
mean.PS.allclassPS.per.sample.pY <- mean(result.PS.count.IP.tum.DZ.3$percentage.pY);round(mean.PS.allclassPS.per.sample.pY, 1) #%pY 90.2


## plot number of PS(all classes)
ggplot(data = result.PS.count.IP.tum.DZ.3)+
  geom_col(aes(x=variant.group.3, y=numberPS, fill=variant.group.2), color="black")+
  theme(legend.position="none") +
  scale_x_discrete(limits= IP.tum.DZ.3.samples.order)+
  scale_fill_viridis(discrete=T, option="inferno")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1.0, vjust = 0.5))+
  ylab("number PS")+
  xlab(NULL)+
  scale_y_continuous(breaks=seq(from=0, to=2400, by=200), limits = c(0,2400))+
  ggtitle("# PS (all classes) per sample" )+
  theme(plot.title = element_text(hjust = 0.5, face= "bold"))+
  theme(axis.text.x = element_text(face= "plain", colour="black", size=14))+
  theme(axis.text.y = element_text(face= "plain", colour="black", size=14))+
  geom_vline(xintercept = c(1.5,2.5,  4.5, 6.5, 9.5, 10.5, 13.5, 15.5, 17.5, 19.5, 22.5, 24.5) , size=0.25, linetype = 2)



### loop over data to count number of PS per sample class 1 PS  ###################################################################################################################################

glimpse(IP.tum.DZ.3)

# Class1 PS filter for PS with localization probability > 0.75
table(IP.tum.DZ.3$Localization.prob >= 0.75) 
IP.tum.DZ.3.class1.COUNT <- IP.tum.DZ.3[which(IP.tum.DZ.3$Localization.prob >= 0.75),]
nrow(IP.tum.DZ.3.class1.COUNT) #3532
table(IP.tum.DZ.3.class1.COUNT$Localization.prob >= 0.75) 

### calculate number of class 1 phosphosite per sample
glimpse(IP.tum.DZ.3.class1.COUNT)
colnames(IP.tum.DZ.3.class1.COUNT)
ncol(IP.tum.DZ.3.class1.COUNT) #95

#loop over data to count number of class 1 PS per sample
result.class1.PS.count.IP.tum.DZ.3 <- c()
sample.name.count <- 1

for(i.1 in seq(from=21, to=95, by=3)){ #53
  print(i.1)
  
  i.2 <- i.1+1
  i.3 <- i.2+1
  
  sample.name <- samples.in.experiment[sample.name.count]
  col.names <- colnames(IP.tum.DZ.3.class1.COUNT)
  col.names.select <- col.names[c(i.1, i.2, i.3)]
  
  temp <- select(IP.tum.DZ.3.class1.COUNT, Amino.acid, col.names.select)
  temp$sum_1_2_3 <- pbapply(temp[, 2:4], 1, function(x) sum(x))
  
  numberPS.A <- filter(temp, sum_1_2_3 > 0)
  numberPS <- nrow(numberPS.A)
  
  nMVs <- filter(temp, sum_1_2_3 <= 0)
  nMVs <- nrow(nMVs)
  
  percentageSTY <- table(numberPS.A$Amino.acid)
  
  temp.result <- tibble(class1.sample=sample.name, class1.numberPS=numberPS, class1.nMVs=nMVs, class1.nS=percentageSTY[1], class1.nT=percentageSTY[2], class1.nY=percentageSTY[3])
  
  result.class1.PS.count.IP.tum.DZ.3  <- bind_rows(result.class1.PS.count.IP.tum.DZ.3 , temp.result)
  
  sample.name.count <- sample.name.count +1
}

result.class1.PS.count.IP.tum.DZ.3 


# calculate percentage MV PS class1
result.class1.PS.count.IP.tum.DZ.3$percentageMV <- pbsapply(result.class1.PS.count.IP.tum.DZ.3 $class1.nMVs, function(x) (100*x)/3532) #version FR 01.12.2021
result.class1.PS.count.IP.tum.DZ.3 


# calculate percentage pS/pT/pY PS class1
for(i in 1:nrow(result.class1.PS.count.IP.tum.DZ.3 )){
  print(i)
  
  total.PS.sample <- result.class1.PS.count.IP.tum.DZ.3$class1.numberPS[i]
  n.pS <- result.class1.PS.count.IP.tum.DZ.3 $class1.nS[i]
  n.pT <- result.class1.PS.count.IP.tum.DZ.3 $class1.nT[i]
  n.pY <- result.class1.PS.count.IP.tum.DZ.3 $class1.nY[i]
  
  result.class1.PS.count.IP.tum.DZ.3 $class1.percentage.pS[i] <- (100*n.pS)/total.PS.sample
  result.class1.PS.count.IP.tum.DZ.3 $class1.percentage.pT[i] <- (100*n.pT)/total.PS.sample
  result.class1.PS.count.IP.tum.DZ.3 $class1.percentage.pY[i] <- (100*n.pY)/total.PS.sample
}

glimpse(result.class1.PS.count.IP.tum.DZ.3 )

#change result data frame and merge results counting with additional sample label information
result.class1.PS.count.IP.tum.DZ.3$class1.sample <- as.numeric(result.class1.PS.count.IP.tum.DZ.3$class1.sample)
glimpse(result.class1.PS.count.IP.tum.DZ.3)

result.class1.PS.count.IP.tum.DZ.3 <- merge(result.class1.PS.count.IP.tum.DZ.3, DZ.IP.tum.labels, by.x = "class1.sample", by.y = "sample", all.x = T, all.y = T, sort = F)
glimpse(result.class1.PS.count.IP.tum.DZ.3)

result.class1.PS.count.IP.tum.DZ.3$variant.group.3 <- paste0(result.class1.PS.count.IP.tum.DZ.3$class1.sample, "_", result.class1.PS.count.IP.tum.DZ.3$variant.group.2)
glimpse(result.class1.PS.count.IP.tum.DZ.3)


#reorder samples
result.class1.PS.count.IP.tum.DZ.3 <- result.class1.PS.count.IP.tum.DZ.3[match(IP.tum.DZ.3.samples.order, result.class1.PS.count.IP.tum.DZ.3$variant.group.3),]
glimpse(result.class1.PS.count.IP.tum.DZ.3)

#calculate means
mean.result.class1.PS.count.IP.tum.DZ.3  <- mean(result.class1.PS.count.IP.tum.DZ.3$class1.numberPS); mean.result.class1.PS.count.IP.tum.DZ.3   #1229.48




### plot number of class 1 PS
ggplot(data = result.class1.PS.count.IP.tum.DZ.3 )+
  geom_col(aes(x=  variant.group.3 , y= class1.numberPS, fill=variant.group.2), color="black")+
  theme(legend.position="none") +
  scale_x_discrete(limits= IP.tum.DZ.3.samples.order)+
  scale_fill_viridis(discrete=T, option="inferno")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1.0, vjust = 0.5))+
  ylab("number PS")+
  xlab(NULL)+
  scale_y_continuous(breaks=seq(from=0, to=2400, by=200), limits = c(0,2400))+
  ggtitle("# PS class1 per sample" )+
  theme(plot.title = element_text(hjust = 0.5, face= "bold"))+
  theme(axis.text.x = element_text(face= "plain", colour="black", size=14))+
  theme(axis.text.y = element_text(face= "plain", colour="black", size=14))+
  geom_vline(xintercept = c(1.5,2.5,  4.5, 6.5, 9.5, 10.5, 13.5, 15.5, 17.5, 19.5, 22.5, 24.5) , size=0.25, linetype = 2)


### plot number class 1: Ser, Thr, Tyr sites
ggplot()+
  geom_col(data= result.class1.PS.count.IP.tum.DZ.3, aes(x=variant.group.3 , y=class1.nY), color="black", fill="yellow")+
  geom_col(data= result.class1.PS.count.IP.tum.DZ.3, aes(x=variant.group.3 , y=class1.nS), color="black", fill="blue")+
  geom_col(data= result.class1.PS.count.IP.tum.DZ.3, aes(x=variant.group.3 , y=class1.nT), color="black", fill="black")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1.0, vjust = 0.5))+
  ylab("class1.nY")+
  xlab(NULL)+
  ggtitle(" # PS class 1 Tyr - yellow \n # PS class 1 Ser - blue  \n # PS class 1 Thr - black") +
  theme(plot.title = element_text(hjust = 0.5, face= "bold"))+  #center title
  theme(axis.text.x = element_text(face= "plain", colour="black", size=14))+
  theme(axis.text.y = element_text(face= "plain", colour="black", size=14))+
  geom_vline(xintercept = c(1.5,2.5,  4.5, 6.5, 9.5, 10.5, 13.5, 15.5, 17.5, 19.5, 22.5, 24.5) , size=0.25, linetype = 2)




#### tumors pTyrIP expand with PerseusPSexpand function (see above) ################################################################################################################
#### tumors pTyrIP expand with PerseusPSexpand function (see above) ################################################################################################################

glimpse(IP.tum.DZ.3)


DZ.tum.IP_expanded <- PerseusPSexpand.17.02.21(IP.tum.DZ.3)

glimpse(DZ.tum.IP_expanded)



#remove rows that have all zero intensities
colnames(DZ.tum.IP_expanded)

DZ.tum.IP_expanded$ZEROcheckINT <- pbapply(DZ.tum.IP_expanded[, c(21:45)], MARGIN=1, function(x) sum(x)) #
table(DZ.tum.IP_expanded$ZEROcheckINT > 0)

DZ.tum.IP_expanded <- filter(DZ.tum.IP_expanded, ZEROcheckINT > 0)
nrow(DZ.tum.IP_expanded) #4759


#remove column ZEROcheckINT
DZ.tum.IP_expanded <- select(DZ.tum.IP_expanded,  -ZEROcheckINT)
glimpse(DZ.tum.IP_expanded)




# tumors pTyrIP normalization & histogram standard median shift version #########################################################################################################################################
# tumors pTyrIP normalization & histogram standard median shift version #########################################################################################################################################

colnames(DZ.tum.IP_expanded)

#log2 transform & replace -Inf with NA
Log2Trafo <- DZ.tum.IP_expanded[, c(21:45)] #
Log2Trafo <- log2(Log2Trafo)
Log2Trafo[Log2Trafo  == -Inf] <- NA

glimpse(Log2Trafo)

DZ.tum.IP_expanded_log2 <- DZ.tum.IP_expanded
DZ.tum.IP_expanded_log2[, c(21:45)] <- Log2Trafo #

glimpse(DZ.tum.IP_expanded_log2, list.len=900)


### calculate median per column RawINT

#reorder
DZ.tum.IP_expanded_log2 <- DZ.tum.IP_expanded_log2 %>% select("id",                          "Proteins",                    "Positions.within.proteins",   "Leading.proteins",            "Protein",                    
                                                              "Protein.names",               "Gene.names",                  "Fasta.headers",               "Localization.prob",           "Number.of.Phospho.(STY)",    
                                                              "Amino.acid",                  "Sequence.window",             "Phospho.(STY).Probabilities", "Position.in.peptide",         "Reverse",                    
                                                              "Potential.contaminant",       "Positions",                   "Position",                    "Peptide.IDs",                 "Mod..peptide.IDs",   
                                                              
                                                              "Intensity.2",
                                                              
                                                              "Intensity.14",
                                                              
                                                              "Intensity.7",
                                                              "Intensity.9",
                                                              
                                                              "Intensity.17",
                                                              "Intensity.19",
                                                              
                                                              "Intensity.4",
                                                              "Intensity.5",
                                                              "Intensity.6",
                                                              
                                                              "Intensity.16",
                                                              
                                                              "Intensity.10",
                                                              "Intensity.11",
                                                              "Intensity.12",
                                                              
                                                              "Intensity.20",
                                                              "Intensity.22",
                                                              
                                                              "Intensity.23",
                                                              "Intensity.24",
                                                              
                                                              "Intensity.25",
                                                              "Intensity.26",
                                                              
                                                              "Intensity.3",
                                                              "Intensity.27",
                                                              "Intensity.28",
                                                              
                                                              "Intensity.29",
                                                              "Intensity.30",
                                                              
                                                              "Intensity.32",
                                                              
                                                              "PS_Multiplicity"
)
glimpse(DZ.tum.IP_expanded_log2)
colnames(DZ.tum.IP_expanded_log2)


ColMedianfindDF <- glimpse(DZ.tum.IP_expanded_log2[, c(21:45)]) #
vectormediancolumnsINT <- c()
for(i in 1:ncol(ColMedianfindDF)){
  print(i)
  #i=5
  temp <- median(ColMedianfindDF[, i], na.rm = TRUE)
  vectormediancolumnsINT <- append(vectormediancolumnsINT, temp)
}
vectormediancolumnsINT

median(vectormediancolumnsINT)


# data frame to plot medians as lines
Mdf <- data.frame(variable=colnames(DZ.tum.IP_expanded_log2[, c(21:45)]), median=vectormediancolumnsINT)

Mdf$variable2 <- c( "2_Fgfr2",
                    
                    "14_Ate1",
                    
                    "7_Bicc1",
                    "9_Bicc1",
                    
                    "17_Tacc2",
                    "19_Tacc2",
                    
                    "4_dE18",
                    "5_dE18",
                    "6_dE18",
                    
                    "16_dE18-Ate1",
                    
                    "10_dE18-Bicc1",
                    "11_dE18-Bicc1",
                    "12_dE18-Bicc1",
                    
                    "20_dE18-Tacc2",
                    "22_dE18-Tacc2",
                    
                    "23_dE18-IGR1",
                    "24_dE18-IGR1",
                    
                    "25_dE18-IGR2",
                    "26_dE18-IGR2",
                    
                    "3_E18-C2",
                    "27_E18-C2",
                    "28_E18-C2",
                    
                    "29_E18-C3",
                    "30_E18-C3",
                    
                    "32_E18-C4" )

Mdf$group <- c("Fgfr2",      "Ate1",     "Bicc1",      "Bicc1",      "Tacc2",      "Tacc2",      "dE18",       "dE18",       "dE18",       "dE18-Ate1",  "dE18-Bicc1",
               "dE18-Bicc1", "dE18-Bicc1", "dE18-Tacc2", "dE18-Tacc2", "dE18-IGR1",  "dE18-IGR1",  "dE18-IGR2",  "dE18-IGR2",  "E18-C2",     "E18-C2",     "E18-C2",    
               "E18-C3",     "E18-C3",     "E18-C4"    )
Mdf

## plain density plot
DensityPlotData <- DZ.tum.IP_expanded_log2[, c(21:45)]
glimpse(DensityPlotData)

DensityPlotDataMelted <- reshape2::melt(DensityPlotData)
glimpse(DensityPlotDataMelted)
head(DensityPlotDataMelted)

ggplot(DensityPlotDataMelted, aes(x=value, colour = variable) ) + 
  geom_density() +
  theme(legend.position="bottom") +
  geom_vline(data=Mdf, aes(xintercept=median, color=variable))

# plot with additional sample names
glimpse(DensityPlotDataMelted)
DensityPlotDataMelted.2 <- merge(DensityPlotDataMelted, Mdf[, c("variable", "variable2", "group")], by.x = "variable", by.y = "variable", all.x = T, sort = F)

ggplot(DensityPlotDataMelted.2, aes(x=value, colour = variable2 ) ) + 
  geom_density() +
  scale_color_viridis(discrete=T, option="magma")+
  theme(legend.position="bottom")+
  ggtitle("mouse tumors DZ pTyr IP\n data distribution") +
  theme(plot.title = element_text(hjust = 0.5, face= "bold"))+
  theme(axis.text.x = element_text(face= "plain", colour="black", size=14))+
  theme(axis.text.y = element_text(face= "plain", colour="black", size=14))


## Start define function to calculate median shift factors 
MedianShiftFactorFinder <- function(arg){
  tempfactor <- 23.2535/arg # log2(1e07) = 23.2535
  return(tempfactor)
}
## Start define function to calculate median shift factors 

## find median shift factors
shiftfactors <- c()
for(i in vectormediancolumnsINT ){
  print(i)
  temp <- MedianShiftFactorFinder(arg = i)
  shiftfactors <- append(shiftfactors, temp)
}
shiftfactors 


## apply median shift factors
temp_for_median_shift <-  DZ.tum.IP_expanded_log2[, c(21:45)]
glimpse(temp_for_median_shift)
glimpse(shiftfactors)

for(q in 1:ncol(temp_for_median_shift)){
  print(q)
  temp_for_median_shift[, q] <- temp_for_median_shift[, q]*shiftfactors[q]
}

glimpse(temp_for_median_shift)


###CHECK
# calculate median per column RawINT
colnames(temp_for_median_shift)

ColMedianfindDF <- temp_for_median_shift
vectormediancolumnsINTCHECK <- c()
for(i in 1:ncol(ColMedianfindDF)){
  print(i)
  temp <- median(ColMedianfindDF[, i], na.rm = TRUE)
  vectormediancolumnsINTCHECK <- append(vectormediancolumnsINTCHECK, temp)
}
vectormediancolumnsINTCHECK

median(vectormediancolumnsINT)

# data frame to plot medians as lines
MdfCHECK <- data.frame(variable=colnames(temp_for_median_shift), median=vectormediancolumnsINTCHECK); MdfCHECK
MdfCHECK

DensityPlotDataCHECK <- temp_for_median_shift
glimpse(DensityPlotDataCHECK)


DensityPlotDataCHECKMelted <- reshape2::melt(DensityPlotDataCHECK)

#plain plot check
ggplot(DensityPlotDataCHECKMelted, aes(x=value, colour = variable) ) + 
  geom_density() +
  theme(legend.position="bottom") +
  geom_vline(data=MdfCHECK, aes(xintercept=median, color=variable))

#plot check with additional sample names
glimpse(DensityPlotDataCHECKMelted)

DensityPlotDataCHECKMelted.2 <- merge(DensityPlotDataCHECKMelted, Mdf[, c("variable", "variable2", "group")], by.x = "variable", by.y = "variable", all.x = T, sort = F)
glimpse(DensityPlotDataCHECKMelted.2)

ggplot(DensityPlotDataCHECKMelted.2, aes(x=value, colour = variable2) ) + 
  geom_density() +
  scale_color_viridis(discrete=T, option="magma")+
  theme(legend.position="bottom")+
  ggtitle("mouse tumors DZ pTyr IP \n data distribution - normalized") +
  theme(plot.title = element_text(hjust = 0.5, face= "bold"))+
  theme(axis.text.x = element_text(face= "plain", colour="black", size=14))+
  theme(axis.text.y = element_text(face= "plain", colour="black", size=14))




### rename columns after normalization
newCOLNAMESnorm <- colnames(temp_for_median_shift); newCOLNAMESnorm

newCOLNAMESnorm <- paste0("Norm.",newCOLNAMESnorm); newCOLNAMESnorm

colnames(temp_for_median_shift) <- newCOLNAMESnorm
glimpse(temp_for_median_shift)


### add norm INT to original Dataframe 
DZ.tum.IP_expanded_log2_norm <- DZ.tum.IP_expanded_log2

DZ.tum.IP_expanded_log2_norm  <- cbind(DZ.tum.IP_expanded_log2_norm , temp_for_median_shift)
glimpse(DZ.tum.IP_expanded_log2_norm, list.len=999)


### box plot after normalization
data.for.boxplot.after.norm.DZ.pTyr.IP.MBR.ON <- DZ.tum.IP_expanded_log2_norm %>% select(                 
  "Norm.Intensity.2",
  
  "Norm.Intensity.14",
  
  "Norm.Intensity.7",
  "Norm.Intensity.9",
  
  "Norm.Intensity.17",
  "Norm.Intensity.19",
  
  "Norm.Intensity.4",
  "Norm.Intensity.5",
  "Norm.Intensity.6",
  
  "Norm.Intensity.16",
  
  "Norm.Intensity.10",
  "Norm.Intensity.11",
  "Norm.Intensity.12",
  
  "Norm.Intensity.20",
  "Norm.Intensity.22",
  
  "Norm.Intensity.23",
  "Norm.Intensity.24",
  
  "Norm.Intensity.25",
  "Norm.Intensity.26",
  
  "Norm.Intensity.3",
  "Norm.Intensity.27",
  "Norm.Intensity.28",
  
  "Norm.Intensity.29",
  "Norm.Intensity.30",
  
  "Norm.Intensity.32")

glimpse(data.for.boxplot.after.norm.DZ.pTyr.IP.MBR.ON)
colnames(data.for.boxplot.after.norm.DZ.pTyr.IP.MBR.ON) <- c("2_Fgfr2",
                                                             
                                                             "14_Ate1",
                                                             
                                                             "7_Bicc1",
                                                             "9_Bicc1",
                                                             
                                                             "17_Tacc2",
                                                             "19_Tacc2",
                                                             
                                                             "4_dE18",
                                                             "5_dE18",
                                                             "6_dE18",
                                                             
                                                             "16_dE18-Ate1",
                                                             
                                                             "10_dE18-Bicc1",
                                                             "11_dE18-Bicc1",
                                                             "12_dE18-Bicc1",
                                                             
                                                             "20_dE18-Tacc2",
                                                             "22_dE18-Tacc2",
                                                             
                                                             "23_dE18-IGR1",
                                                             "24_dE18-IGR1",
                                                             
                                                             "25_dE18-IGR2",
                                                             "26_dE18-IGR2",
                                                             
                                                             "3_E18-C2",
                                                             "27_E18-C2",
                                                             "28_E18-C2",
                                                             
                                                             "29_E18-C3",
                                                             "30_E18-C3",
                                                             
                                                             "32_E18-C4" )
glimpse(data.for.boxplot.after.norm.DZ.pTyr.IP.MBR.ON)

melted.data.for.boxplot.after.norm.DZ.pTyr.IP.MBR.ON <- reshape2::melt(data.for.boxplot.after.norm.DZ.pTyr.IP.MBR.ON)
glimpse(melted.data.for.boxplot.after.norm.DZ.pTyr.IP.MBR.ON )

#plot
ggplot(data = melted.data.for.boxplot.after.norm.DZ.pTyr.IP.MBR.ON)+
  geom_boxplot(mapping = aes(x= variable, y = value), outlier.size = 0.5, notch = TRUE, notchwidth = 0.1)+ 
  scale_x_discrete(limits= c("2_Fgfr2",
                             
                             "14_Ate1",
                             
                             "7_Bicc1",
                             "9_Bicc1",
                             
                             "17_Tacc2",
                             "19_Tacc2",
                             
                             "4_dE18",
                             "5_dE18",
                             "6_dE18",
                             
                             "16_dE18-Ate1",
                             
                             "10_dE18-Bicc1",
                             "11_dE18-Bicc1",
                             "12_dE18-Bicc1",
                             
                             "20_dE18-Tacc2",
                             "22_dE18-Tacc2",
                             
                             "23_dE18-IGR1",
                             "24_dE18-IGR1",
                             
                             "25_dE18-IGR2",
                             "26_dE18-IGR2",
                             
                             "3_E18-C2",
                             "27_E18-C2",
                             "28_E18-C2",
                             
                             "29_E18-C3",
                             "30_E18-C3",
                             
                             "32_E18-C4" ))+
  scale_fill_viridis(discrete=T, option="inferno")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1.0, vjust = 0.5))+

  theme(legend.position="none", legend.justification = "center")+
  ggtitle("tumors DZ pTyr IP - normalized data - PS level") +
  theme(plot.title = element_text(hjust = 0.5, face= "bold"))+
  theme(axis.text.x = element_text(face= "plain", colour="black", size=14))+
  theme(axis.text.y = element_text(face= "plain", colour="black", size=14))+
  xlab(NULL)+ 
  ylab("log2(Int.)") 


### box plot WITHOUT normalization
glimpse(DZ.tum.IP_expanded_log2_norm)

data.for.boxplot.after.norm.DZ.pTyr.IP.MBR.ON <- DZ.tum.IP_expanded_log2_norm %>% select( "Intensity.2",
                                                                                          
                                                                                          "Intensity.14",
                                                                                          
                                                                                          "Intensity.7",
                                                                                          "Intensity.9",
                                                                                          
                                                                                          "Intensity.17",
                                                                                          "Intensity.19",
                                                                                          
                                                                                          "Intensity.4",
                                                                                          "Intensity.5",
                                                                                          "Intensity.6",
                                                                                          
                                                                                          "Intensity.16",
                                                                                          
                                                                                          "Intensity.10",
                                                                                          "Intensity.11",
                                                                                          "Intensity.12",
                                                                                          
                                                                                          "Intensity.20",
                                                                                          "Intensity.22",
                                                                                          
                                                                                          "Intensity.23",
                                                                                          "Intensity.24",
                                                                                          
                                                                                          "Intensity.25",
                                                                                          "Intensity.26",
                                                                                          
                                                                                          "Intensity.3",
                                                                                          "Intensity.27",
                                                                                          "Intensity.28",
                                                                                          
                                                                                          "Intensity.29",
                                                                                          "Intensity.30",
                                                                                          
                                                                                          "Intensity.32")
glimpse(data.for.boxplot.after.norm.DZ.pTyr.IP.MBR.ON)
colnames(data.for.boxplot.after.norm.DZ.pTyr.IP.MBR.ON) <- c("2_Fgfr2",
                                                             
                                                             "14_Ate1",
                                                             
                                                             "7_Bicc1",
                                                             "9_Bicc1",
                                                             
                                                             "17_Tacc2",
                                                             "19_Tacc2",
                                                             
                                                             "4_dE18",
                                                             "5_dE18",
                                                             "6_dE18",
                                                             
                                                             "16_dE18-Ate1",
                                                             
                                                             "10_dE18-Bicc1",
                                                             "11_dE18-Bicc1",
                                                             "12_dE18-Bicc1",
                                                             
                                                             "20_dE18-Tacc2",
                                                             "22_dE18-Tacc2",
                                                             
                                                             "23_dE18-IGR1",
                                                             "24_dE18-IGR1",
                                                             
                                                             "25_dE18-IGR2",
                                                             "26_dE18-IGR2",
                                                             
                                                             "3_E18-C2",
                                                             "27_E18-C2",
                                                             "28_E18-C2",
                                                             
                                                             "29_E18-C3",
                                                             "30_E18-C3",
                                                             
                                                             "32_E18-C4" )
glimpse(data.for.boxplot.after.norm.DZ.pTyr.IP.MBR.ON)

melted.data.for.boxplot.after.norm.DZ.pTyr.IP.MBR.ON <- reshape2::melt(data.for.boxplot.after.norm.DZ.pTyr.IP.MBR.ON)
glimpse(melted.data.for.boxplot.after.norm.DZ.pTyr.IP.MBR.ON )

#plot
ggplot(data = melted.data.for.boxplot.after.norm.DZ.pTyr.IP.MBR.ON)+
  geom_boxplot(mapping = aes(x= variable, y = value), outlier.size = 0.5, notch = TRUE, notchwidth = 0.1)+ 
  scale_x_discrete(limits= c("2_Fgfr2",
                             
                             "14_Ate1",
                             
                             "7_Bicc1",
                             "9_Bicc1",
                             
                             "17_Tacc2",
                             "19_Tacc2",
                             
                             "4_dE18",
                             "5_dE18",
                             "6_dE18",
                             
                             "16_dE18-Ate1",
                             
                             "10_dE18-Bicc1",
                             "11_dE18-Bicc1",
                             "12_dE18-Bicc1",
                             
                             "20_dE18-Tacc2",
                             "22_dE18-Tacc2",
                             
                             "23_dE18-IGR1",
                             "24_dE18-IGR1",
                             
                             "25_dE18-IGR2",
                             "26_dE18-IGR2",
                             
                             "3_E18-C2",
                             "27_E18-C2",
                             "28_E18-C2",
                             
                             "29_E18-C3",
                             "30_E18-C3",
                             
                             "32_E18-C4" ))+
  
  theme(axis.text.x = element_text(angle = 90, hjust = 1.0, vjust = 0.5))+
  theme(legend.position="none", legend.justification = "center")+
  ggtitle("tumors DZ pTyr IP - NOT normalized data - PS level") +
  theme(plot.title = element_text(hjust = 0.5, face= "bold"))+ 
  theme(axis.text.x = element_text(face= "plain", colour="black", size=14))+
  theme(axis.text.y = element_text(face= "plain", colour="black", size=14))+
  xlab(NULL)+ 
  ylab("log2(Int.)") 



#####################################################################################################################################################################
#### analysis via class 1 PS ##################################################################################################################################
glimpse(DZ.tum.IP_expanded_log2_norm)



#filter class 1
DZ.tum.IP_expanded_log2_normclass1.PS <- DZ.tum.IP_expanded_log2_norm %>% filter(Localization.prob>= 0.75)
glimpse(DZ.tum.IP_expanded_log2_normclass1.PS) #3796
length(unique(DZ.tum.IP_expanded_log2_normclass1.PS$Sequence.window)) #3532

# add additional column to better calculate number PS 
# take Proteins and add Psite without multiplicity info ###
for(i in 1:nrow(DZ.tum.IP_expanded_log2_normclass1.PS)){
  print(i)
  temp_prot_name <- DZ.tum.IP_expanded_log2_normclass1.PS$Proteins[i]
  tempPSite <- DZ.tum.IP_expanded_log2_normclass1.PS$Positions[i]
  tempAA <- DZ.tum.IP_expanded_log2_normclass1.PS$Amino.acid[i]
  
  DZ.tum.IP_expanded_log2_normclass1.PS$Prot.Name_AApos[i] <- paste0(temp_prot_name, "_", tempAA, tempPSite)
}
length(unique(DZ.tum.IP_expanded_log2_normclass1.PS$Prot.Name_AApos)) #3532 == the number of class 1 PS in exp

glimpse(DZ.tum.IP_expanded_log2_normclass1.PS)

#  prepare new column
# take genename and add Psite ###
for(i in 1:nrow(DZ.tum.IP_expanded_log2_normclass1.PS)){
  print(i)
  temp_gene_name <- DZ.tum.IP_expanded_log2_normclass1.PS$Gene.names[i]
  tempPSite <- DZ.tum.IP_expanded_log2_normclass1.PS$Positions[i]
  tempAA <- DZ.tum.IP_expanded_log2_normclass1.PS$Amino.acid[i]
  tempPSmultiplicity <- DZ.tum.IP_expanded_log2_normclass1.PS$PS_Multiplicity[i]
  
  DZ.tum.IP_expanded_log2_normclass1.PS$Gene.Name_AApos[i] <- paste0(temp_gene_name, "_", tempAA, tempPSite, "x", tempPSmultiplicity)
}
glimpse(DZ.tum.IP_expanded_log2_normclass1.PS)


#add identifier for plotting reason e.g. PTMSEA later
DZ.tum.IP_expanded_log2_normclass1.PS$ID.FR.all.C1.PS <- 1:nrow(DZ.tum.IP_expanded_log2_normclass1.PS)
glimpse(DZ.tum.IP_expanded_log2_normclass1.PS)
m

#save data
#save(DZ.tum.IP_expanded_log2_normclass1.PS, file="DZ.tum.IP_expanded_log2_normclass1.PS.Rdata")
#load("DZ.tum.IP_expanded_log2_normclass1.PS.Rdata")
#glimpse(DZ.tum.IP_expanded_log2_normclass1.PS)


#############################################################################################################################################################
#############################################################################################################################################################
# correlation tumors pTyrIP class 1 PS

glimpse( DZ.tum.IP_expanded_log2_normclass1.PS)


#select columns and order
DF_for_correlation_plot.pTyr.IP.class1PS <-  DZ.tum.IP_expanded_log2_normclass1.PS %>% select(contains("Norm.Intensity.")               )
glimpse(DF_for_correlation_plot.pTyr.IP.class1PS)



##check if there are all NA rows 
DF_for_correlation_plot.pTyr.IP.class1PS$na.rows <- pbapply(DF_for_correlation_plot.pTyr.IP.class1PS, 1, function(x) sum(is.na(x)))
glimpse(DF_for_correlation_plot.pTyr.IP.class1PS)
table(DF_for_correlation_plot.pTyr.IP.class1PS$na.rows)

# change sample names and remove column na.rows
DF_for_correlation_plot.pTyr.IP.class1PS <- DF_for_correlation_plot.pTyr.IP.class1PS %>% select(-na.rows )
glimpse(DF_for_correlation_plot.pTyr.IP.class1PS)
colnames(DF_for_correlation_plot.pTyr.IP.class1PS)
colnames(DF_for_correlation_plot.pTyr.IP.class1PS) <- c("2_Fgfr2",
                                                        
                                                        "14_Fgfr2-Ate1",
                                                        
                                                        "7_Fgfr2-Bicc1",
                                                        "9_Fgfr2-Bicc1",
                                                        
                                                        "17_Fgfr2-Tacc2",
                                                        "19_Fgfr2-Tacc2",
                                                        
                                                        "4_Fgfr2-dE18",
                                                        "5_Fgfr2-dE18",
                                                        "6_Fgfr2-dE18",
                                                        
                                                        "16_Fgfr2-dE18-Ate1",
                                                        
                                                        "10_Fgfr2-dE18-Bicc1",
                                                        "11_Fgfr2-dE18-Bicc1",
                                                        "12_Fgfr2-dE18-Bicc1",
                                                        
                                                        "20_Fgfr2-dE18-Tacc2",
                                                        "22_Fgfr2-dE18-Tacc2",
                                                        
                                                        "23_Fgfr2-dE18-IGR1",
                                                        "24_Fgfr2-dE18-IGR1",
                                                        
                                                        "25_Fgfr2-dE18-IGR2",
                                                        "26_Fgfr2-dE18-IGR2",
                                                        
                                                        "3_Fgfr2-E18-C2",
                                                        "27_Fgfr2-E18-C2",
                                                        "28_Fgfr2-E18-C2",
                                                        
                                                        "29_Fgfr2-E18-C3",
                                                        "30_Fgfr2-E18-C3",
                                                        
                                                        "32_Fgfr2-E18-C4" )
glimpse(DF_for_correlation_plot.pTyr.IP.class1PS)

#prepare additional sample names
tum.pTyrIP.name.change <- tibble(orig.col.name = colnames(DF_for_correlation_plot.pTyr.IP.class1PS))
tum.pTyrIP.name.change$new.col.name.1 <- pbsapply(tum.pTyrIP.name.change$orig.col.name, function(x) unlist(str_split(string=x, pattern="_"))[2]   )
tum.pTyrIP.name.change$new.col.name.2 <- paste0(tum.pTyrIP.name.change$new.col.name.1, ".", c(2,2,  1,3,  1,3,  1,2,3,   2,  1,2,3,   1,3,  1,2,   1,2,  1,2,3,   1,2,   2))
print(tum.pTyrIP.name.change , n=45)
colnames(DF_for_correlation_plot.pTyr.IP.class1PS) <- tum.pTyrIP.name.change$new.col.name.2
glimpse(DF_for_correlation_plot.pTyr.IP.class1PS)


###correlate
corr2 <- cor(DF_for_correlation_plot.pTyr.IP.class1PS, method = "pearson", use = "na.or.complete")
glimpse(corr2)
head(corr2)
min(corr2) #0.3650069
max(corr2)
round(min(corr2), 2)



# prepare annotation data
annot_df_for_heatmap <- data.frame(samples = colnames(DF_for_correlation_plot.pTyr.IP.class1PS),
                                   group = c("Fgfr2",      "Fgfr2-Ate1",     "Fgfr2-Bicc1",      "Fgfr2-Bicc1",      "Fgfr2-Tacc2",      "Fgfr2-Tacc2",      "Fgfr2-dE18",       "Fgfr2-dE18",       "Fgfr2-dE18",       "Fgfr2-dE18-Ate1",  "Fgfr2-dE18-Bicc1",
                                             "Fgfr2-dE18-Bicc1", "Fgfr2-dE18-Bicc1", "Fgfr2-dE18-Tacc2", "Fgfr2-dE18-Tacc2", "Fgfr2-dE18-IGR1",  "Fgfr2-dE18-IGR1",  "Fgfr2-dE18-IGR2",  "Fgfr2-dE18-IGR2",  "Fgfr2-E18-C2",     "Fgfr2-E18-C2",     "Fgfr2-E18-C2",    
                                             "Fgfr2-E18-C3",     "Fgfr2-E18-C3",     "Fgfr2-E18-C4"    )
                                   )
annot_df_for_heatmap 


# make shorter annotation data frame
annot_df_for_heatmap.short <- data.frame( 
  group  = annot_df_for_heatmap$group
)

annot_df_for_heatmap.short

# definecolors for annotation bar
annot.colors = list(group = c("Fgfr2"="blue",
                              "Fgfr2-Ate1"="white",
                              "Fgfr2-Bicc1"="black",
                              "Fgfr2-Tacc2"="deeppink",
                              "Fgfr2-dE18"="red",
                              "Fgfr2-dE18-Ate1"="lavender",
                              "Fgfr2-dE18-Bicc1"="grey40",
                              "Fgfr2-dE18-Tacc2"="lightskyblue",
                              
                              "Fgfr2-dE18-IGR1"="gold1",
                              "Fgfr2-dE18-IGR2"="gold2",
                              
                              "Fgfr2-E18-C2"="chocolate2",
                              "Fgfr2-E18-C3"="chocolate3",
                              "Fgfr2-E18-C4"="chocolate4"))
annot.colors

# Create the heatmap annotation
ha.fr <- HeatmapAnnotation(group=annot_df_for_heatmap.short$group,
                           treatment=annot_df_for_heatmap.short$treatment,
                           col = annot.colors,
                           annotation_legend_param = list(grid_height = unit(8, "mm"))
)
ha.fr


# define colors
FR.heatmap.colors.2 <-colorRamp2(c(0.4, 1.0), c("grey90", "grey10"))



#pdf("DZ.tumors.pTyrIP.correlation.pdf", width=25/2.54, height=20/2.54, useDingbats=FALSE) 
Heatmap(corr2, 
        name = "corr.coeff",
        
        top_annotation = ha.fr,
        col = FR.heatmap.colors.2,
        
        clustering_distance_rows = function(m) as.dist((1-m)/2), 
        clustering_method_rows = "ward.D2",
        
        clustering_distance_columns = function(m) as.dist((1-m)/2), 
        clustering_method_columns = "ward.D2",
        
        column_dend_height = unit(40, "mm"),
        row_dend_width = unit(40, "mm"),
        
        heatmap_legend_param = list(ncol = 1, nrow = 1, legend_height = unit(60, "mm"))
        
)
#dev.off()


####################################################################################################################################
#############################################################################################################################################################

## heatmap for all Fgfr2 sites in dataset for all groups with collapsing _1/_2/_3 via sum 
### PS numbers also for NP_963895.2

#load("DZ.tum.IP_expanded_log2_normclass1.PS.Rdata")
glimpse( DZ.tum.IP_expanded_log2_normclass1.PS)

# filter for protein candidate of interest, pick columns of interest and rename
HM.plotdata.pTyrIP.tum <- DZ.tum.IP_expanded_log2_normclass1.PS %>% filter(stringr::str_detect(Proteins , "P21803")) 
glimpse(HM.plotdata.pTyrIP.tum)

HM.plotdata.pTyrIP.tum <- HM.plotdata.pTyrIP.tum %>% select("Gene.Name_AApos", 
                                                            "Norm.Intensity.2",
                                                            
                                                            "Norm.Intensity.14",
                                                            
                                                            "Norm.Intensity.7",
                                                            "Norm.Intensity.9",
                                                            
                                                            "Norm.Intensity.17",
                                                            "Norm.Intensity.19",
                                                            
                                                            "Norm.Intensity.4",
                                                            "Norm.Intensity.5",
                                                            "Norm.Intensity.6",
                                                            
                                                            "Norm.Intensity.16",
                                                            
                                                            "Norm.Intensity.10",
                                                            "Norm.Intensity.11",
                                                            "Norm.Intensity.12",
                                                            
                                                            "Norm.Intensity.20",
                                                            "Norm.Intensity.22",
                                                            
                                                            "Norm.Intensity.23",
                                                            "Norm.Intensity.24",
                                                            
                                                            "Norm.Intensity.25",
                                                            "Norm.Intensity.26",
                                                            
                                                            "Norm.Intensity.3",
                                                            "Norm.Intensity.27",
                                                            "Norm.Intensity.28",
                                                            
                                                            "Norm.Intensity.29",
                                                            "Norm.Intensity.30",
                                                            
                                                            "Norm.Intensity.32",
                                                            "Proteins",
                                                            "Positions.within.proteins",
                                                            "Amino.acid",
                                                            "Gene.names" ,
                                                            "PS_Multiplicity")
glimpse(HM.plotdata.pTyrIP.tum)

colnames(HM.plotdata.pTyrIP.tum) <- c("Gene.Name_AApos", 
                                      "2_Fgfr2",
                                      
                                      "14_Fgfr2-Ate1",
                                      
                                      "7_Fgfr2-Bicc1",
                                      "9_Fgfr2-Bicc1",
                                      
                                      "17_Fgfr2-Tacc2",
                                      "19_Fgfr2-Tacc2",
                                      
                                      "4_Fgfr2-dE18",
                                      "5_Fgfr2-dE18",
                                      "6_Fgfr2-dE18",
                                      
                                      "16_Fgfr2-dE18-Ate1",
                                      
                                      "10_Fgfr2-dE18-Bicc1",
                                      "11_Fgfr2-dE18-Bicc1",
                                      "12_Fgfr2-dE18-Bicc1",
                                      
                                      "20_Fgfr2-dE18-Tacc2",
                                      "22_Fgfr2-dE18-Tacc2",
                                      
                                      "23_Fgfr2-dE18-IGR1",
                                      "24_Fgfr2-dE18-IGR1",
                                      
                                      "25_Fgfr2-dE18-IGR2",
                                      "26_Fgfr2-dE18-IGR2",
                                      
                                      "3_Fgfr2-E18-C2",
                                      "27_Fgfr2-E18-C2",
                                      "28_Fgfr2-E18-C2",
                                      
                                      "29_Fgfr2-E18-C3",
                                      "30_Fgfr2-E18-C3",
                                      
                                      "32_Fgfr2-E18-C4",
                                      "Proteins",
                                      "Positions.within.proteins",
                                      "Amino.acid",
                                      "Gene.names" ,
                                      "PS_Multiplicity")
glimpse(HM.plotdata.pTyrIP.tum)

# define function to get isoform specific site position from MaxQuant output
position.in.mq.string <- function(mq.string, mq.position.string, Amino.acid, Gene.names, PS_Multiplicity, prot.of.interest){
  temp.df <- tibble(proteins = unlist(str_split(mq.string, pattern = ";")),
                    positions = unlist(str_split(mq.position.string, pattern = ";")),
                    Amino.acid = Amino.acid, 
                    Gene.names = Gene.names,
                    PS_Multiplicity = PS_Multiplicity
  )
  temp.df
  temp.df <- temp.df %>% filter( proteins %in% c(prot.of.interest))
  temp.df
  isoform.corrected.PS.position <-  temp.df$positions
  return(isoform.corrected.PS.position )
}


### add isoform site numbers

# isoform 1  ="P21803"
for(i in 1:nrow(HM.plotdata.pTyrIP.tum)){
  print(i)
  HM.plotdata.pTyrIP.tum$isoform.1.corrected.PS.position[i] <- position.in.mq.string(prot.of.interest ="P21803", #isoform 1
                                                                                     Gene.names=HM.plotdata.pTyrIP.tum$Gene.names[i],
                                                                                     Amino.acid=HM.plotdata.pTyrIP.tum$Amino.acid[i], 
                                                                                     mq.position.string=HM.plotdata.pTyrIP.tum$Positions.within.proteins[i],
                                                                                     mq.string = HM.plotdata.pTyrIP.tum$Proteins[i],
                                                                                     PS_Multiplicity=HM.plotdata.pTyrIP.tum$PS_Multiplicity[i])
}
HM.plotdata.pTyrIP.tum$isoform.1.corrected.Gene.Name_AApos = paste0(HM.plotdata.pTyrIP.tum$Gene.names, "_", HM.plotdata.pTyrIP.tum$Amino.acid, HM.plotdata.pTyrIP.tum$isoform.1.corrected.PS.position, "x", HM.plotdata.pTyrIP.tum$PS_Multiplicity)
glimpse(HM.plotdata.pTyrIP.tum)


# isoform 2 ="P21803-2"
for(i in 1:nrow(HM.plotdata.pTyrIP.tum)){
  print(i)
  HM.plotdata.pTyrIP.tum$isoform.2.corrected.PS.position[i] <- position.in.mq.string(prot.of.interest ="P21803-2", #isoform 2
                                                                                     Gene.names=HM.plotdata.pTyrIP.tum$Gene.names[i],
                                                                                     Amino.acid=HM.plotdata.pTyrIP.tum$Amino.acid[i], 
                                                                                     mq.position.string=HM.plotdata.pTyrIP.tum$Positions.within.proteins[i],
                                                                                     mq.string = HM.plotdata.pTyrIP.tum$Proteins[i],
                                                                                     PS_Multiplicity=HM.plotdata.pTyrIP.tum$PS_Multiplicity[i])
}
HM.plotdata.pTyrIP.tum$isoform.2.corrected.Gene.Name_AApos = paste0(HM.plotdata.pTyrIP.tum$Gene.names, "_", HM.plotdata.pTyrIP.tum$Amino.acid, HM.plotdata.pTyrIP.tum$isoform.2.corrected.PS.position, "x", HM.plotdata.pTyrIP.tum$PS_Multiplicity)
glimpse(HM.plotdata.pTyrIP.tum)


#prepare helper data frame for later order 
HM.plotdata.pTyrIP.tum.melted <- reshape2::melt(HM.plotdata.pTyrIP.tum)
glimpse(HM.plotdata.pTyrIP.tum.melted)
colnames(HM.plotdata.pTyrIP.tum.melted) <- c("Gene.Name_AApos", "Proteins","Positions.within.proteins", "Amino.acid","Gene.names", "PS_Multiplicity", "isoform.1.corrected.PS.position", "isoform.1.corrected.Gene.Name_AApos","isoform.2.corrected.PS.position", "isoform.2.corrected.Gene.Name_AApos",  "variable","log2.int")
glimpse(HM.plotdata.pTyrIP.tum.melted)


order.df.HM.plotdata.pTyrIP.tum.melted <- unique(tibble(Gene.Name_AApos = HM.plotdata.pTyrIP.tum.melted$Gene.Name_AApos,
                                                        Proteins = HM.plotdata.pTyrIP.tum.melted$Proteins, 
                                                        Positions.within.proteins = HM.plotdata.pTyrIP.tum.melted$Positions.within.proteins,
                                                        Amino.acid = HM.plotdata.pTyrIP.tum.melted$Amino.acid,
                                                        Gene.names = HM.plotdata.pTyrIP.tum.melted$Gene.names, 
                                                        PS_Multiplicity =HM.plotdata.pTyrIP.tum.melted$PS_Multiplicity,
                                                        isoform.1.corrected.PS.position =HM.plotdata.pTyrIP.tum.melted$isoform.1.corrected.PS.position,
                                                        isoform.1.corrected.Gene.Name_AApos = paste0(HM.plotdata.pTyrIP.tum.melted$Gene.names, "_", HM.plotdata.pTyrIP.tum.melted$Amino.acid, HM.plotdata.pTyrIP.tum.melted$isoform.1.corrected.PS.position, "x", HM.plotdata.pTyrIP.tum.melted$PS_Multiplicity),
                                                        isoform.2.corrected.PS.position =HM.plotdata.pTyrIP.tum.melted$isoform.2.corrected.PS.position,
                                                        isoform.2.corrected.Gene.Name_AApos = paste0(HM.plotdata.pTyrIP.tum.melted$Gene.names, "_", HM.plotdata.pTyrIP.tum.melted$Amino.acid, HM.plotdata.pTyrIP.tum.melted$isoform.2.corrected.PS.position, "x", HM.plotdata.pTyrIP.tum.melted$PS_Multiplicity)
))
glimpse(order.df.HM.plotdata.pTyrIP.tum.melted)

### for better visualization/plotting => collapse _1/_2/_3 via  sum
### note: sum is based on _1/_2/_3 intensities in dataframe at hand => unfiltered data might contain more _1/_2/_3 PS that are not included in the sum here

HM.plotdata.pTyrIP.tum.collapsed <- HM.plotdata.pTyrIP.tum
glimpse(HM.plotdata.pTyrIP.tum.collapsed )

HM.plotdata.pTyrIP.tum.collapsed$isoform.1.corrected.Gene.Name_AApos.collapsed <- pbsapply(HM.plotdata.pTyrIP.tum.collapsed$isoform.1.corrected.Gene.Name_AApos, function(x) str_remove_all(string=x, pattern="x_\\d") )
HM.plotdata.pTyrIP.tum.collapsed$isoform.2.corrected.Gene.Name_AApos.collapsed <- pbsapply(HM.plotdata.pTyrIP.tum.collapsed$isoform.2.corrected.Gene.Name_AApos, function(x) str_remove_all(string=x, pattern="x_\\d") )
glimpse(HM.plotdata.pTyrIP.tum.collapsed)
## unlog, do collapsing, relog
HM.plotdata.pTyrIP.tum.collapsed [, c(   "2_Fgfr2",
                                         
                                         "14_Fgfr2-Ate1",
                                         
                                         "7_Fgfr2-Bicc1",
                                         "9_Fgfr2-Bicc1",
                                         
                                         "17_Fgfr2-Tacc2",
                                         "19_Fgfr2-Tacc2",
                                         
                                         "4_Fgfr2-dE18",
                                         "5_Fgfr2-dE18",
                                         "6_Fgfr2-dE18",
                                         
                                         "16_Fgfr2-dE18-Ate1",
                                         
                                         "10_Fgfr2-dE18-Bicc1",
                                         "11_Fgfr2-dE18-Bicc1",
                                         "12_Fgfr2-dE18-Bicc1",
                                         
                                         "20_Fgfr2-dE18-Tacc2",
                                         "22_Fgfr2-dE18-Tacc2",
                                         
                                         "23_Fgfr2-dE18-IGR1",
                                         "24_Fgfr2-dE18-IGR1",
                                         
                                         "25_Fgfr2-dE18-IGR2",
                                         "26_Fgfr2-dE18-IGR2",
                                         
                                         "3_Fgfr2-E18-C2",
                                         "27_Fgfr2-E18-C2",
                                         "28_Fgfr2-E18-C2",
                                         
                                         "29_Fgfr2-E18-C3",
                                         "30_Fgfr2-E18-C3",
                                         
                                         "32_Fgfr2-E18-C4")] <- 2^HM.plotdata.pTyrIP.tum.collapsed [, c( "2_Fgfr2",
                                                                                                         
                                                                                                         "14_Fgfr2-Ate1",
                                                                                                         
                                                                                                         "7_Fgfr2-Bicc1",
                                                                                                         "9_Fgfr2-Bicc1",
                                                                                                         
                                                                                                         "17_Fgfr2-Tacc2",
                                                                                                         "19_Fgfr2-Tacc2",
                                                                                                         
                                                                                                         "4_Fgfr2-dE18",
                                                                                                         "5_Fgfr2-dE18",
                                                                                                         "6_Fgfr2-dE18",
                                                                                                         
                                                                                                         "16_Fgfr2-dE18-Ate1",
                                                                                                         
                                                                                                         "10_Fgfr2-dE18-Bicc1",
                                                                                                         "11_Fgfr2-dE18-Bicc1",
                                                                                                         "12_Fgfr2-dE18-Bicc1",
                                                                                                         
                                                                                                         "20_Fgfr2-dE18-Tacc2",
                                                                                                         "22_Fgfr2-dE18-Tacc2",
                                                                                                         
                                                                                                         "23_Fgfr2-dE18-IGR1",
                                                                                                         "24_Fgfr2-dE18-IGR1",
                                                                                                         
                                                                                                         "25_Fgfr2-dE18-IGR2",
                                                                                                         "26_Fgfr2-dE18-IGR2",
                                                                                                         
                                                                                                         "3_Fgfr2-E18-C2",
                                                                                                         "27_Fgfr2-E18-C2",
                                                                                                         "28_Fgfr2-E18-C2",
                                                                                                         
                                                                                                         "29_Fgfr2-E18-C3",
                                                                                                         "30_Fgfr2-E18-C3",
                                                                                                         
                                                                                                         "32_Fgfr2-E18-C4")] 
glimpse(HM.plotdata.pTyrIP.tum.collapsed )
##collapse
HM.plotdata.pTyrIP.tum.collapsed.2  <- HM.plotdata.pTyrIP.tum.collapsed %>% 
  group_by(isoform.1.corrected.Gene.Name_AApos.collapsed, isoform.2.corrected.Gene.Name_AApos.collapsed) %>% 
  summarise_at(c( "2_Fgfr2",
                  
                  "14_Fgfr2-Ate1",
                  
                  "7_Fgfr2-Bicc1",
                  "9_Fgfr2-Bicc1",
                  
                  "17_Fgfr2-Tacc2",
                  "19_Fgfr2-Tacc2",
                  
                  "4_Fgfr2-dE18",
                  "5_Fgfr2-dE18",
                  "6_Fgfr2-dE18",
                  
                  "16_Fgfr2-dE18-Ate1",
                  
                  "10_Fgfr2-dE18-Bicc1",
                  "11_Fgfr2-dE18-Bicc1",
                  "12_Fgfr2-dE18-Bicc1",
                  
                  "20_Fgfr2-dE18-Tacc2",
                  "22_Fgfr2-dE18-Tacc2",
                  
                  "23_Fgfr2-dE18-IGR1",
                  "24_Fgfr2-dE18-IGR1",
                  
                  "25_Fgfr2-dE18-IGR2",
                  "26_Fgfr2-dE18-IGR2",
                  
                  "3_Fgfr2-E18-C2",
                  "27_Fgfr2-E18-C2",
                  "28_Fgfr2-E18-C2",
                  
                  "29_Fgfr2-E18-C3",
                  "30_Fgfr2-E18-C3",
                  
                  "32_Fgfr2-E18-C4"), sum, na.rm = TRUE)
glimpse(HM.plotdata.pTyrIP.tum.collapsed.2)
HM.plotdata.pTyrIP.tum.collapsed.2[HM.plotdata.pTyrIP.tum.collapsed.2  == 0] <- NA
##relog
HM.plotdata.pTyrIP.tum.collapsed.2[c( "2_Fgfr2",
                                      
                                      "14_Fgfr2-Ate1",
                                      
                                      "7_Fgfr2-Bicc1",
                                      "9_Fgfr2-Bicc1",
                                      
                                      "17_Fgfr2-Tacc2",
                                      "19_Fgfr2-Tacc2",
                                      
                                      "4_Fgfr2-dE18",
                                      "5_Fgfr2-dE18",
                                      "6_Fgfr2-dE18",
                                      
                                      "16_Fgfr2-dE18-Ate1",
                                      
                                      "10_Fgfr2-dE18-Bicc1",
                                      "11_Fgfr2-dE18-Bicc1",
                                      "12_Fgfr2-dE18-Bicc1",
                                      
                                      "20_Fgfr2-dE18-Tacc2",
                                      "22_Fgfr2-dE18-Tacc2",
                                      
                                      "23_Fgfr2-dE18-IGR1",
                                      "24_Fgfr2-dE18-IGR1",
                                      
                                      "25_Fgfr2-dE18-IGR2",
                                      "26_Fgfr2-dE18-IGR2",
                                      
                                      "3_Fgfr2-E18-C2",
                                      "27_Fgfr2-E18-C2",
                                      "28_Fgfr2-E18-C2",
                                      
                                      "29_Fgfr2-E18-C3",
                                      "30_Fgfr2-E18-C3",
                                      
                                      "32_Fgfr2-E18-C4")] <- log2(HM.plotdata.pTyrIP.tum.collapsed.2[c( "2_Fgfr2",
                                                                                                        
                                                                                                        "14_Fgfr2-Ate1",
                                                                                                        
                                                                                                        "7_Fgfr2-Bicc1",
                                                                                                        "9_Fgfr2-Bicc1",
                                                                                                        
                                                                                                        "17_Fgfr2-Tacc2",
                                                                                                        "19_Fgfr2-Tacc2",
                                                                                                        
                                                                                                        "4_Fgfr2-dE18",
                                                                                                        "5_Fgfr2-dE18",
                                                                                                        "6_Fgfr2-dE18",
                                                                                                        
                                                                                                        "16_Fgfr2-dE18-Ate1",
                                                                                                        
                                                                                                        "10_Fgfr2-dE18-Bicc1",
                                                                                                        "11_Fgfr2-dE18-Bicc1",
                                                                                                        "12_Fgfr2-dE18-Bicc1",
                                                                                                        
                                                                                                        "20_Fgfr2-dE18-Tacc2",
                                                                                                        "22_Fgfr2-dE18-Tacc2",
                                                                                                        
                                                                                                        "23_Fgfr2-dE18-IGR1",
                                                                                                        "24_Fgfr2-dE18-IGR1",
                                                                                                        
                                                                                                        "25_Fgfr2-dE18-IGR2",
                                                                                                        "26_Fgfr2-dE18-IGR2",
                                                                                                        
                                                                                                        "3_Fgfr2-E18-C2",
                                                                                                        "27_Fgfr2-E18-C2",
                                                                                                        "28_Fgfr2-E18-C2",
                                                                                                        
                                                                                                        "29_Fgfr2-E18-C3",
                                                                                                        "30_Fgfr2-E18-C3",
                                                                                                        
                                                                                                        "32_Fgfr2-E18-C4")] )
glimpse(HM.plotdata.pTyrIP.tum.collapsed.2)


### reshape and reorder, define order for plot
HM.plotdata.pTyrIP.tum.collapsed.melted <- reshape2::melt(HM.plotdata.pTyrIP.tum.collapsed.2)
glimpse(HM.plotdata.pTyrIP.tum.collapsed.melted)
colnames(HM.plotdata.pTyrIP.tum.collapsed.melted ) <- c("isoform.1.corrected.Gene.Name_AApos.collapsed","isoform.2.corrected.Gene.Name_AApos.collapsed", "variable","log2.int")
glimpse(HM.plotdata.pTyrIP.tum.collapsed.melted )


#define order for plot, get gene names, corrected positions, aminco acid etc. and order according to number
order.df.HM.plotdata.pTyrIP.tum.collapsed.melted <- unique(tibble(isoform.1.corrected.Gene.Name_AApos.collapsed = HM.plotdata.pTyrIP.tum.collapsed.melted$isoform.1.corrected.Gene.Name_AApos.collapsed,
                                                                  isoform.2.corrected.Gene.Name_AApos.collapsed = HM.plotdata.pTyrIP.tum.collapsed.melted$isoform.2.corrected.Gene.Name_AApos.collapsed))
order.df.HM.plotdata.pTyrIP.tum.collapsed.melted


order.df.HM.plotdata.pTyrIP.tum.collapsed.melted$Gene <- pbsapply(order.df.HM.plotdata.pTyrIP.tum.collapsed.melted$isoform.1.corrected.Gene.Name_AApos.collapsed , function(x) unlist(str_split(string=x,  pattern="_") )[1])

order.df.HM.plotdata.pTyrIP.tum.collapsed.melted$isoform.1.corrected.AAPos <- pbsapply(order.df.HM.plotdata.pTyrIP.tum.collapsed.melted$isoform.1.corrected.Gene.Name_AApos.collapsed, function(x) unlist(str_split(string=x,  pattern="_") )[2])
order.df.HM.plotdata.pTyrIP.tum.collapsed.melted$isoform.2.corrected.AAPos <- pbsapply(order.df.HM.plotdata.pTyrIP.tum.collapsed.melted$isoform.2.corrected.Gene.Name_AApos.collapsed, function(x) unlist(str_split(string=x,  pattern="_") )[2])

order.df.HM.plotdata.pTyrIP.tum.collapsed.melted$isoform.corrected.Amino.acid <- pbsapply(order.df.HM.plotdata.pTyrIP.tum.collapsed.melted$isoform.1.corrected.AAPos, function(x) unlist((str_extract_all(string=x,  pattern="[A-Z]") )))

order.df.HM.plotdata.pTyrIP.tum.collapsed.melted$isoform.1.corrected.Pos <- pbsapply(order.df.HM.plotdata.pTyrIP.tum.collapsed.melted$isoform.1.corrected.AAPos, function(x) unlist((str_extract_all(string=x,  pattern="[0-9]+") )))
order.df.HM.plotdata.pTyrIP.tum.collapsed.melted$isoform.2.corrected.Pos <- pbsapply(order.df.HM.plotdata.pTyrIP.tum.collapsed.melted$isoform.2.corrected.AAPos, function(x) unlist((str_extract_all(string=x,  pattern="[0-9]+") )))

order.df.HM.plotdata.pTyrIP.tum.collapsed.melted <- order.df.HM.plotdata.pTyrIP.tum.collapsed.melted %>% arrange(desc(isoform.2.corrected.Pos, isoform.2.corrected.Amino.acid))
order.df.HM.plotdata.pTyrIP.tum.collapsed.melted


isoformcorrected.collapes.order.29.11.21 <- order.df.HM.plotdata.pTyrIP.tum.collapsed.melted$isoform.2.corrected.Gene.Name_AApos.collapsed
isoformcorrected.collapes.order.29.11.21


glimpse(HM.plotdata.pTyrIP.tum.collapsed.melted)


#plot the heatmap with site numbers for P21803-2
ggplot(
  HM.plotdata.pTyrIP.tum.collapsed.melted , aes(x=variable, y=isoform.2.corrected.Gene.Name_AApos.collapsed, fill=log2.int)) +
  geom_tile() + 
  scale_fill_viridis( option = "inferno", na.value = "grey30",name="log2(intensity)") + 
  coord_equal()+
  scale_y_discrete(limits=  isoformcorrected.collapes.order.29.11.21)+
  theme(legend.position="bottom", legend.justification = "center")+
  theme(axis.line = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  theme(legend.title.align=0.5)+
  theme(plot.title = element_text(hjust = 0.5, face="bold", size=16)) +
  theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1))+
  theme(axis.text.y= element_text(size=8))+ 
  xlab(NULL) + 
  xlab(NULL) + 
  ylab(NULL) +
  geom_vline(xintercept = c(1.5,2.5,  4.5, 6.5, 9.5, 10.5, 13.5, 15.5, 17.5, 19.5, 22.5, 24.5) , size=1.0, linetype = "solid", color="white")


### add corresponding human FGFR2 phosphosite numbers to the heatmap plot ######################################################

glimpse(HM.plotdata.pTyrIP.tum.collapsed.2)
order.df.HM.plotdata.pTyrIP.tum.collapsed.melted

# online Clustal Omega alignment showed that P21803-2 (707 aa) and NP_963895.2 (726 aa) are very similar. The difference is N-terminal: 
# 19aa  (MGLPSTWRYGRGPGIGTVT) are added for  NP_963895.2. Everything else is identical.
# Thus P21803-2 position +19 should be the positions for NP_963895.2

#prepare additional data frame and calculate the new site numbers
HM.plotdata.pTyrIP.tum.collapsed.3 <- HM.plotdata.pTyrIP.tum.collapsed.2

HM.plotdata.pTyrIP.tum.collapsed.3$isoform.1.corrected.AApos <- pbsapply(HM.plotdata.pTyrIP.tum.collapsed.3$isoform.1.corrected.Gene.Name_AApos.collapsed, function(x) unlist(str_split(string=x, pattern="_"))[2])
HM.plotdata.pTyrIP.tum.collapsed.3$isoform.1.corrected.pos <- pbsapply(HM.plotdata.pTyrIP.tum.collapsed.3$isoform.1.corrected.AApos, function(x) unlist(str_extract_all(string=x, pattern="[0-9]+")))

HM.plotdata.pTyrIP.tum.collapsed.3$isoform.2.corrected.AApos <- pbsapply(HM.plotdata.pTyrIP.tum.collapsed.3$isoform.2.corrected.Gene.Name_AApos.collapsed, function(x) unlist(str_split(string=x, pattern="_"))[2])
HM.plotdata.pTyrIP.tum.collapsed.3$isoform.2.corrected.pos <- pbsapply(HM.plotdata.pTyrIP.tum.collapsed.3$isoform.2.corrected.AApos, function(x) unlist(str_extract_all(string=x, pattern="[0-9]+")))
glimpse(HM.plotdata.pTyrIP.tum.collapsed.3)

HM.plotdata.pTyrIP.tum.collapsed.3$pos.NP_963895.2 <- as.numeric(HM.plotdata.pTyrIP.tum.collapsed.3$isoform.2.corrected.pos) +19
glimpse(HM.plotdata.pTyrIP.tum.collapsed.3)


# go to www.phosphosite.org/ => search for Fgfr2 => look a t site table => enable other species/isoforms mouse => copy list to Excle for lookup of human phosphosite numbers starting from mouse

human.PSP.site.number.lookup <- read_excel("mouse.Fgfr2.modifications.entry.PSP.with.iso.2.mouse.and.human.02.12.2021.xlsx")
human.PSP.site.number.lookup$site.Fgfr2.mouse.iso2_AA <- pbsapply(human.PSP.site.number.lookup$site.Fgfr2.mouse.iso2, function(x) str_extract(string=x, pattern="[A-Z]"))
human.PSP.site.number.lookup$site.Fgfr2.mouse.iso1 <- pbsapply(human.PSP.site.number.lookup$site.Fgfr2.mouse.iso1, function(x) str_remove_all(string=x, pattern="-p"))
human.PSP.site.number.lookup$site.Fgfr2.mouse.iso2 <- pbsapply(human.PSP.site.number.lookup$site.Fgfr2.mouse.iso2, function(x) str_remove_all(string=x, pattern="-p"))
human.PSP.site.number.lookup$site.FGFR2.human.iso1<- pbsapply(human.PSP.site.number.lookup$site.FGFR2.human.iso1, function(x) str_remove_all(string=x, pattern="-p"))

glimpse(human.PSP.site.number.lookup)

# filter for phosphosites only
human.PSP.site.number.lookup.2 <- human.PSP.site.number.lookup %>% filter(site.Fgfr2.mouse.iso2_AA %in% c("S", "T", "Y"))
glimpse(human.PSP.site.number.lookup.2)

#merge PSP mouse/human site information with heatmap data
HM.plotdata.pTyrIP.tum.collapsed.4 <- merge(HM.plotdata.pTyrIP.tum.collapsed.3, human.PSP.site.number.lookup.2[, c("site.Fgfr2.mouse.iso1", "site.FGFR2.human.iso1" )], by.x = "isoform.1.corrected.AApos", by.y = "site.Fgfr2.mouse.iso1", all.x = T, all.y = F, sort = F)
HM.plotdata.pTyrIP.tum.collapsed.4$Gene <- pbsapply(HM.plotdata.pTyrIP.tum.collapsed.4$isoform.1.corrected.Gene.Name_AApos.collapsed, function(x) unlist(str_split(string=x,  pattern="_") )[1])
HM.plotdata.pTyrIP.tum.collapsed.4 <- HM.plotdata.pTyrIP.tum.collapsed.4  %>% arrange(HM.plotdata.pTyrIP.tum.collapsed.4$isoform.1.corrected.pos)

glimpse(HM.plotdata.pTyrIP.tum.collapsed.4)

#some human sites have to be entered manually since not all sites are in PSP yet #they correspond mostly to numbers from P21803-1
unname(HM.plotdata.pTyrIP.tum.collapsed.4$isoform.1.corrected.AApos)
HM.plotdata.pTyrIP.tum.collapsed.4$site.FGFR2.human.iso1

HM.plotdata.pTyrIP.tum.collapsed.4$AA.site.FGFR2.human.iso1_manual.change <- c("Y466","Y561","Y575","Y586","S587","Y588","Y608","Y616","Y656","Y657","T660")
HM.plotdata.pTyrIP.tum.collapsed.4$site.FGFR2.human.iso1_manual.change <- pbsapply(HM.plotdata.pTyrIP.tum.collapsed.4$AA.site.FGFR2.human.iso1_manual.change, function(x) str_remove_all(string=x, pattern="[A-Z]"))
glimpse(HM.plotdata.pTyrIP.tum.collapsed.4)

#make.new.plot.label
#e.g. Fgfr2_Y466|352|371|466 => 1st number: P21803-1; 2nd number: P21803-2; 3rd number: NP_963895.2; 4th number: P21802-1
HM.plotdata.pTyrIP.tum.collapsed.4$combined.plot.label <- paste0(HM.plotdata.pTyrIP.tum.collapsed.4$isoform.1.corrected.Gene.Name_AApos.collapsed, "|", HM.plotdata.pTyrIP.tum.collapsed.4$isoform.2.corrected.pos, "|",  HM.plotdata.pTyrIP.tum.collapsed.4$pos.NP_963895.2, "|", HM.plotdata.pTyrIP.tum.collapsed.4$site.FGFR2.human.iso1_manual.change)
glimpse(HM.plotdata.pTyrIP.tum.collapsed.4)

isoform.corrected.collapsed.order.02.12.2021 <- rev(HM.plotdata.pTyrIP.tum.collapsed.4$combined.plot.label)


#change samples names and reorder for plot
tum.pTyrIP.name.change <- tibble(orig.col.name = colnames(HM.plotdata.pTyrIP.tum.collapsed.4 %>% select(`2_Fgfr2` : `32_Fgfr2-E18-C4`)) ) 
tum.pTyrIP.name.change$new.col.name.1 <- pbsapply(tum.pTyrIP.name.change$orig.col.name, function(x) unlist(str_split(string=x, pattern="_"))[2]   )
tum.pTyrIP.name.change$new.col.name.2 <- paste0(tum.pTyrIP.name.change$new.col.name.1, ".", c(2,2,  1,3,  1,3,  1,2,3,   2,  1,2,3,   1,3,  1,2,   1,2,  1,2,3,   1,2,   2))
print(tum.pTyrIP.name.change , n=45)



colnames(HM.plotdata.pTyrIP.tum.collapsed.4) <- c("isoform.1.corrected.AApos","isoform.1.corrected.Gene.Name_AApos.collapsed","isoform.2.corrected.Gene.Name_AApos.collapsed",
                                                  tum.pTyrIP.name.change$new.col.name.2,                      
                                                  "isoform.1.corrected.pos","isoform.2.corrected.AApos","isoform.2.corrected.pos","pos.NP_963895.2","site.FGFR2.human.iso1","Gene","AA.site.FGFR2.human.iso1_manual.change",       
                                                  "site.FGFR2.human.iso1_manual.change","combined.plot.label")
glimpse(HM.plotdata.pTyrIP.tum.collapsed.4)



# reshape for plot
melted.HM.plotdata.pTyrIP.tum.collapsed.4 <- reshape2::melt(HM.plotdata.pTyrIP.tum.collapsed.4 %>% select(combined.plot.label, "Fgfr2.2" : "Fgfr2-E18-C4.2"))
glimpse(melted.HM.plotdata.pTyrIP.tum.collapsed.4)

custom.sample.order.DZ.tum.pTyrIP <- c("Fgfr2.2", 
                                       "Fgfr2-dE18.1", "Fgfr2-dE18.2", "Fgfr2-dE18.3", 
                                       "Fgfr2-Bicc1.1", "Fgfr2-Bicc1.3", 
                                       "Fgfr2-dE18-Bicc1.1", "Fgfr2-dE18-Bicc1.2", "Fgfr2-dE18-Bicc1.3", 
                                       "Fgfr2-Ate1.2",
                                       "Fgfr2-dE18-Ate1.2",
                                       "Fgfr2-Tacc2.1", "Fgfr2-Tacc2.3",
                                       "Fgfr2-dE18-Tacc2.1", "Fgfr2-dE18-Tacc2.3",
                                       "Fgfr2-dE18-IGR1.1", "Fgfr2-dE18-IGR1.2", 
                                       "Fgfr2-dE18-IGR2.1", "Fgfr2-dE18-IGR2.2",
                                       "Fgfr2-E18-C2.1", "Fgfr2-E18-C2.2", 
                                       "Fgfr2-E18-C3.1", "Fgfr2-E18-C3.2",
                                       "Fgfr2-E18-C4.2"   
)
custom.sample.order.DZ.tum.pTyrIP
length(custom.sample.order.DZ.tum.pTyrIP)

### plot
ggplot(
  melted.HM.plotdata.pTyrIP.tum.collapsed.4 , aes(x=variable, y=combined.plot.label, fill=value)) +
  geom_tile() + 
  scale_fill_viridis( option = "inferno", na.value = "grey30",name="log2(intensity)") + 
  coord_equal()+
  scale_y_discrete(limits=  isoform.corrected.collapsed.order.02.12.2021 )+
  scale_x_discrete(limits=  custom.sample.order.DZ.tum.pTyrIP)+
  theme(legend.position="bottom", legend.justification = "center")+
  theme(axis.line = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  theme(legend.title.align=0.5)+
  theme(plot.title = element_text(hjust = 0.5, face="bold", size=16)) + 
  theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1))+
  theme(axis.text.y= element_text(size=8))+
  xlab(NULL) + 
  xlab(NULL) + 
  ylab(NULL) +
  theme(axis.text.x = element_blank()) +
  theme(axis.ticks.x = element_blank())+
  
  geom_vline(xintercept = c(1.5, 4.5, 6.5, 9.5,10.5, 11.5, 13.5, 15.5, 17.5, 19.5, 21.5, 23.5) , size=1.0, linetype = "solid", color="white")+
  annotate("text", x = 1, y = nrow(HM.plotdata.pTyrIP.tum.collapsed.2)+1, label = "Fgfr2", fontface = "bold", angle='90', hjust = 0, size=4)+
  annotate("text", x = 3, y = nrow(HM.plotdata.pTyrIP.tum.collapsed.2)+1, label = "Fgfr2-dE18", fontface = "bold", angle='90', hjust = 0, size=4)+
  annotate("text", x = 5.5, y = nrow(HM.plotdata.pTyrIP.tum.collapsed.2)+1, label = "Fgfr2-Bicc1", fontface = "bold", angle='90', hjust = 0, size=4)+
  annotate("text", x = 8, y = nrow(HM.plotdata.pTyrIP.tum.collapsed.2)+1, label = "Fgfr2-dE18-Bicc1", fontface = "bold", angle='90', hjust = 0, size=4)+
  annotate("text", x = 10, y = nrow(HM.plotdata.pTyrIP.tum.collapsed.2)+1, label = "Fgfr2-Ate1", fontface = "bold", angle='90', hjust = 0, size=4)+
  annotate("text", x = 11, y = nrow(HM.plotdata.pTyrIP.tum.collapsed.2)+1, label = "Fgfr2-dE18-Ate1", fontface = "bold", angle='90', hjust = 0, size=4)+
  annotate("text", x = 12.5, y = nrow(HM.plotdata.pTyrIP.tum.collapsed.2)+1, label = "Fgfr2-Tacc2", fontface = "bold", angle='90', hjust = 0, size=4)+
  annotate("text", x = 14.5, y = nrow(HM.plotdata.pTyrIP.tum.collapsed.2)+1, label = "Fgfr2-dE18-Tacc2", fontface = "bold", angle='90', hjust = 0, size=4)+
  annotate("text", x = 16.5, y = nrow(HM.plotdata.pTyrIP.tum.collapsed.2)+1, label = "Fgfr2-dE18-IGR1", fontface = "bold", angle='90', hjust = 0, size=4)+
  annotate("text", x = 18.5, y = nrow(HM.plotdata.pTyrIP.tum.collapsed.2)+1, label = "Fgfr2-dE18-IGR2", fontface = "bold", angle='90', hjust = 0, size=4)+
  annotate("text", x = 20.5, y = nrow(HM.plotdata.pTyrIP.tum.collapsed.2)+1, label = "Fgfr2-E18-C2", fontface = "bold", angle='90', hjust = 0, size=4)+
  annotate("text", x = 22.5, y = nrow(HM.plotdata.pTyrIP.tum.collapsed.2)+1, label = "Fgfr2-E18-C3", fontface = "bold", angle='90', hjust = 0, size=4)+
  annotate("text", x = 24, y = nrow(HM.plotdata.pTyrIP.tum.collapsed.2)+1, label = "Fgfr2-E18-C4", fontface = "bold", angle='90', hjust = 0, size=4)+
  annotate("text", x = 24, y = nrow(HM.plotdata.pTyrIP.tum.collapsed.2)+10, label = "",  color = "transparent")

#ggsave("DZ.tumors.pTyrIP.Fgfr2.sites.pdf", useDingbats=FALSE,  width = 14, height =14, units = "cm")




















