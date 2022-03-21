#_______________________________________________________________________________________________________________________
# 20.01.2022
# 
# Project: proteomics.truncated.FGFR2.is.oncogene.cancer.Zingg.et.al 
# 
# Script name: tumors.global.phospho.IMAC.public
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
  allDFcolNAMES <- colnames(DF1) #; allDFcolNAMES
  
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
### load data tumors IMAC DDA from MaxQuant Phospho(STY)Sites.txt (Rdata)

### sample overview
#A	Fgfr2	1
#B	Fgfr2-E18-C2	3
#C	Fgfr2-dE18	4
#D	Fgfr2-Bicc1	7
#E	Fgfr2-dE18-Bicc1	10
#F	Fgfr2-Ate1	13
#G	Fgfr2-dE18-Ate1	15
#H	Fgfr2-Tacc2	17
#I	Fgfr2-dE18-Tacc2	20
#J	Fgfr2-dE18-IGR1	23
#K	Fgfr2-dE18-IGR2	25
#L	Fgfr2-E18-C3	29
#M	Fgfr2-E18-C4	31
#N	Fgfr2	2
#O	Fgfr2-dE18	5
#P	Fgfr2-Bicc1	8
#Q	Fgfr2-dE18-Bicc1	11
#R	Fgfr2-Ate1	14
#S	Fgfr2-dE18-Ate1	16
#T	Fgfr2-Tacc2	18
#U	Fgfr2-dE18-Tacc2	21
#V	Fgfr2-dE18-IGR1	24
#W	Fgfr2-dE18-IGR2	26
#X	Fgfr2-E18-C2	27
#Y	Fgfr2-E18-C3	30
#Z	Fgfr2-E18-C4	32
#AA	Fgfr2-dE18	6
#AB	Fgfr2-Bicc1	9
#AC	Fgfr2-dE18-Bicc1	12
#AD	Fgfr2-Tacc2	19
#AE	Fgfr2-dE18-Tacc2	22
#AF	Fgfr2-E18-C2	28

# load Zingg tumor global phospho IMAC sample information
load("DZ.IMAC.tum.labels.Rdata")
glimpse(DZ.IMAC.tum.labels)

# load data
load("tumors.IMAC.Phospho(STY)Sites.Rdata")
glimpse(IMAC.tum.DZ)

nrow(IMAC.tum.DZ) #29329
ncol(IMAC.tum.DZ) #467
colnames(IMAC.tum.DZ)
length(unique(IMAC.tum.DZ$id)) #29329

### select columns of interest
IMAC.tum.DZ.2 <- IMAC.tum.DZ %>% select(
  
  "id",                          "Proteins",                    "Positions within proteins",   "Leading proteins",            "Protein",                    
  "Protein names",               "Gene names",                  "Fasta headers",               "Localization prob",           "Number of Phospho (STY)",    
  "Amino acid",                  "Sequence window",             "Phospho (STY) Probabilities", "Position in peptide",         "Reverse",                    
  "Potential contaminant",       "Positions",                   "Position",                    "Peptide IDs",                 "Mod. peptide IDs",      
  
  matches("___[0-9]"), -Intensity___1, -Intensity___2, -Intensity___3
  
)

glimpse(IMAC.tum.DZ.2)


### remove spaces from column names
column.names <- colnames(IMAC.tum.DZ.2)
column.names.changed <- unname(sapply(column.names , function(x) str_replace_all(string=x, pattern=" ", replacement=".")))
colnames(IMAC.tum.DZ.2) <- column.names.changed

glimpse(IMAC.tum.DZ.2)
colnames(IMAC.tum.DZ.2)


#change column classes
IMAC.tum.DZ.2 <- mutate_at(IMAC.tum.DZ.2,  21:ncol(IMAC.tum.DZ.2), list(as.numeric) )

glimpse(IMAC.tum.DZ.2)



### remove reverse hits
table(IMAC.tum.DZ.2$Reverse) #289+
IMAC.tum.DZ.3 <- IMAC.tum.DZ.2[IMAC.tum.DZ.2$Reverse =="",]
nrow(IMAC.tum.DZ.3) #29040


### remove contaminants
table(IMAC.tum.DZ.3$Potential.contaminant) #312+
IMAC.tum.DZ.3 <- IMAC.tum.DZ.3[IMAC.tum.DZ.3$Potential.contaminant == "",]
nrow(IMAC.tum.DZ.3) #28728


### remove rows completely zero
glimpse(IMAC.tum.DZ.3)
colnames(IMAC.tum.DZ.3)
ncol(IMAC.tum.DZ.3) #116
IMAC.tum.DZ.3$zero.check <- pbapply(IMAC.tum.DZ.3[, c(21:116)], 1, function(x) sum(x) ) #
table(IMAC.tum.DZ.3$zero.check >0) #F6801 T21927 

IMAC.tum.DZ.3.complete.zero.rows <- filter(IMAC.tum.DZ.3, zero.check <=0)
glimpse(IMAC.tum.DZ.3.complete.zero.rows) #

IMAC.tum.DZ.3 <- filter(IMAC.tum.DZ.3, zero.check > 0)
nrow(IMAC.tum.DZ.3) #21927
table(IMAC.tum.DZ.3$zero.check >0)


length(unique(IMAC.tum.DZ.3$id))#21927 = total number PS
table(IMAC.tum.DZ.3$Amino.acid) #S18291  T3160   Y476 
length(unique(IMAC.tum.DZ.3$Sequence.window)) #21909

### remove column zero.check
IMAC.tum.DZ.3 <- IMAC.tum.DZ.3 %>% select(-zero.check)
glimpse(IMAC.tum.DZ.3)


### create variable with short sample names with order as in Phospho(STY)Sites.txt
samples.in.experiment <- colnames(IMAC.tum.DZ.3)
samples.in.experiment <- samples.in.experiment[21:ncol(IMAC.tum.DZ.3)]
samples.in.experiment <- str_remove(string=samples.in.experiment, pattern="Intensity.")
samples.in.experiment <- str_remove(string=samples.in.experiment, pattern="___[0-9]")
samples.in.experiment <- unique(samples.in.experiment)
samples.in.experiment

### loop over data to count number of PS per sample ##############################################################################
colnames(IMAC.tum.DZ.3)
ncol(IMAC.tum.DZ.3) #16

result.PS.count.IMAC.tum.DZ.3 <- c()
sample.name.count <- 1

for(i.1 in seq(from=21, to=116, by=3)){ 
  print(i.1)
  
  i.2 <- i.1+1
  i.3 <- i.2+1
  
  sample.name <- samples.in.experiment[sample.name.count]
  print(sample.name)
  col.names <- colnames(IMAC.tum.DZ.3); col.names
  col.names.select <- col.names[c(i.1, i.2, i.3)]
  
  temp <- select(IMAC.tum.DZ.3,Amino.acid, col.names.select)
  temp$sum_1_2_3 <- pbapply(temp[, 2:4], 1, function(x) sum(x))
  
  numberPS.A <- filter(temp, sum_1_2_3 > 0)
  numberPS <- nrow(numberPS.A)
  
  nMVs <- filter(temp, sum_1_2_3 <= 0)
  nMVs <- nrow(nMVs)
  
  percentageSTY <- table(numberPS.A$Amino.acid); percentageSTY
  
  temp.result <- tibble(sample=sample.name, numberPS=numberPS, nMVs=nMVs, nS=percentageSTY[1], nT=percentageSTY[2], nY=percentageSTY[3])
  
  result.PS.count.IMAC.tum.DZ.3  <- bind_rows(result.PS.count.IMAC.tum.DZ.3, temp.result)
  
  sample.name.count <- sample.name.count +1
  
}
result.PS.count.IMAC.tum.DZ.3

glimpse(result.PS.count.IMAC.tum.DZ.3)

mean(result.PS.count.IMAC.tum.DZ.3$numberPS) #7157.5


#merge results counting with additional sample label information
glimpse(DZ.IMAC.tum.labels)

result.PS.count.IMAC.tum.DZ.3 <- merge(result.PS.count.IMAC.tum.DZ.3, DZ.IMAC.tum.labels, by.x = "sample", by.y = "Proteomics.Sample.Label", all.x = T, all.y = T, sort = F)
glimpse(result.PS.count.IMAC.tum.DZ.3)

#define sample order
IMAC.tum.DZ.3.samples.order <- c(
  "1_Fgfr2",
  "2_Fgfr2",
  
  "13_Fgfr2-Ate1",
  "14_Fgfr2-Ate1",
  
  "7_Fgfr2-Bicc1",
  "8_Fgfr2-Bicc1",
  "9_Fgfr2-Bicc1",
  
  "17_Fgfr2-Tacc2",
  "18_Fgfr2-Tacc2",
  "19_Fgfr2-Tacc2",
  
  "4_Fgfr2-dE18",
  "5_Fgfr2-dE18",
  "6_Fgfr2-dE18",
  
  "15_Fgfr2-dE18-Ate1",
  "16_Fgfr2-dE18-Ate1",
  
  "10_Fgfr2-dE18-Bicc1",
  "11_Fgfr2-dE18-Bicc1",
  "12_Fgfr2-dE18-Bicc1",
  
  "20_Fgfr2-dE18-Tacc2",
  "21_Fgfr2-dE18-Tacc2",
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
  
  "31_Fgfr2-E18-C4" ,
  "32_Fgfr2-E18-C4"    
)
IMAC.tum.DZ.3.samples.order

## calculate percentage MV PS all classes level
nrow(IMAC.tum.DZ.3)  #21927 = total number PS

result.PS.count.IMAC.tum.DZ.3$percentageMV <- pbsapply(result.PS.count.IMAC.tum.DZ.3$nMVs, function(x) (100*x)/21927) #
result.PS.count.IMAC.tum.DZ.3


# calculate percentage pS/pT/pY PS all classes level
for(i in 1:nrow(result.PS.count.IMAC.tum.DZ.3)){
  print(i)
  
  total.PS.sample <- result.PS.count.IMAC.tum.DZ.3$numberPS[i]
  n.pS <- result.PS.count.IMAC.tum.DZ.3$nS[i]
  n.pT <- result.PS.count.IMAC.tum.DZ.3$nT[i]
  n.pY <- result.PS.count.IMAC.tum.DZ.3$nY[i]
  
  result.PS.count.IMAC.tum.DZ.3$percentage.pS[i] <- (100*n.pS)/total.PS.sample
  result.PS.count.IMAC.tum.DZ.3$percentage.pT[i] <- (100*n.pT)/total.PS.sample
  result.PS.count.IMAC.tum.DZ.3$percentage.pY[i] <- (100*n.pY)/total.PS.sample
}
result.PS.count.IMAC.tum.DZ.3
glimpse(result.PS.count.IMAC.tum.DZ.3)

#reorder samples
result.PS.count.IMAC.tum.DZ.3 <- result.PS.count.IMAC.tum.DZ.3[match(IMAC.tum.DZ.3.samples.order, result.PS.count.IMAC.tum.DZ.3$Proteomics.Sample.Label_2),]
result.PS.count.IMAC.tum.DZ.3

# calcuate means
mean.PS.allclassPS.per.sample <- mean(result.PS.count.IMAC.tum.DZ.3$numberPS); mean.PS.allclassPS.per.sample #7157.5

mean.PS.allclassPS.per.sample.pS <- mean(result.PS.count.IMAC.tum.DZ.3$percentage.pS);round(mean.PS.allclassPS.per.sample.pS, 1) #%pS 87.6
mean.PS.allclassPS.per.sample.pT <- mean(result.PS.count.IMAC.tum.DZ.3$percentage.pT);round(mean.PS.allclassPS.per.sample.pT, 1) #%pT 11.1
mean.PS.allclassPS.per.sample.pY <- mean(result.PS.count.IMAC.tum.DZ.3$percentage.pY);round(mean.PS.allclassPS.per.sample.pY, 1) #%pY 1.3


### plot number of PS (all classes) per sample
ggplot(data = result.PS.count.IMAC.tum.DZ.3)+
  geom_col(aes(x=Proteomics.Sample.Label_2, y=numberPS, fill=variant.group), color="black")+
  theme(legend.position="none") +
  scale_x_discrete(limits= IMAC.tum.DZ.3.samples.order)+
  scale_fill_viridis(discrete=T, option="inferno")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1.0, vjust = 0.5))+
  ylab("number PS")+
  xlab(NULL)+
  scale_y_continuous(breaks=seq(from=0, to=9200, by=400), limits = c(0,9200))+
  ggtitle("# PS (all classes) per sample" )+
  theme(plot.title = element_text(hjust = 0.5, face= "bold"))+
  theme(axis.text.x = element_text(face= "plain", colour="black", size=14))+
  theme(axis.text.y = element_text(face= "plain", colour="black", size=14)) +
  geom_vline(xintercept = c(2.5, 4.5, 7.5, 10.5, 13.5, 15.5, 18.5, 21.5, 23.5, 25.5, 28.5, 30.5) , size=0.25, linetype = 2)




### loop over data to count number of PS per sample class 1 PS ################################################################################################################################### ###################################################################################################################################

glimpse(IMAC.tum.DZ.3)

# Class1 PS: filter for PS with localization probability > 0.75
table(IMAC.tum.DZ.3$Localization.prob >= 0.75) 
IMAC.tum.DZ.3.class1.COUNT <- IMAC.tum.DZ.3[which(IMAC.tum.DZ.3$Localization.prob >= 0.75),]
nrow(IMAC.tum.DZ.3.class1.COUNT) #16945
table(IMAC.tum.DZ.3.class1.COUNT$Localization.prob >= 0.75) 

### calculate number of class 1 phosphosite per sample
glimpse(IMAC.tum.DZ.3.class1.COUNT)
colnames(IMAC.tum.DZ.3.class1.COUNT)
ncol(IMAC.tum.DZ.3.class1.COUNT) #116

#loop over data to count number of class 1 PS per sample
result.class1.PS.count.IMAC.tum.DZ.3 <- c()
sample.name.count <- 1

for(i.1 in seq(from=21, to=116, by=3)){
  print(i.1)
 
  i.2 <- i.1+1
  i.3 <- i.2+1
  
  sample.name <- samples.in.experiment[sample.name.count]
  col.names <- colnames(IMAC.tum.DZ.3.class1.COUNT)
  col.names.select <- col.names[c(i.1, i.2, i.3)]
  
  temp <- select(IMAC.tum.DZ.3.class1.COUNT, Amino.acid, col.names.select)
  temp$sum_1_2_3 <- pbapply(temp[, 2:4], 1, function(x) sum(x))
  
  numberPS.A <- filter(temp, sum_1_2_3 > 0)
  numberPS <- nrow(numberPS.A)
  
  nMVs <- filter(temp, sum_1_2_3 <= 0)
  nMVs <- nrow(nMVs); nMVs
  
  percentageSTY <- table(numberPS.A$Amino.acid)
  
  temp.result <- tibble(class1.sample=sample.name, class1.numberPS=numberPS, class1.nMVs=nMVs, class1.nS=percentageSTY[1], class1.nT=percentageSTY[2], class1.nY=percentageSTY[3])
  
  result.class1.PS.count.IMAC.tum.DZ.3  <- bind_rows(result.class1.PS.count.IMAC.tum.DZ.3 , temp.result)
  
  sample.name.count <- sample.name.count +1
}
result.class1.PS.count.IMAC.tum.DZ.3 


# calculate percentage MV PS class1
result.class1.PS.count.IMAC.tum.DZ.3$percentageMV <- pbsapply(result.class1.PS.count.IMAC.tum.DZ.3 $class1.nMVs, function(x) (100*x)/16945) 
result.class1.PS.count.IMAC.tum.DZ.3 


# calculate percentage pS/pT/pY PS class1
for(i in 1:nrow(result.class1.PS.count.IMAC.tum.DZ.3 )){
  print(i)
  
  total.PS.sample <- result.class1.PS.count.IMAC.tum.DZ.3 $class1.numberPS[i]
  n.pS <- result.class1.PS.count.IMAC.tum.DZ.3 $class1.nS[i]
  n.pT <- result.class1.PS.count.IMAC.tum.DZ.3 $class1.nT[i]
  n.pY <- result.class1.PS.count.IMAC.tum.DZ.3 $class1.nY[i]
  
  result.class1.PS.count.IMAC.tum.DZ.3 $class1.percentage.pS[i] <- (100*n.pS)/total.PS.sample
  result.class1.PS.count.IMAC.tum.DZ.3 $class1.percentage.pT[i] <- (100*n.pT)/total.PS.sample
  result.class1.PS.count.IMAC.tum.DZ.3 $class1.percentage.pY[i] <- (100*n.pY)/total.PS.sample
}
glimpse(result.class1.PS.count.IMAC.tum.DZ.3 )


#merge results counting with additional sample label information
glimpse(DZ.IMAC.tum.labels)

result.class1.PS.count.IMAC.tum.DZ.3 <- merge(result.class1.PS.count.IMAC.tum.DZ.3, DZ.IMAC.tum.labels, by.x = "class1.sample", by.y = "Proteomics.Sample.Label", all.x = T, all.y = T, sort = F)
glimpse(result.class1.PS.count.IMAC.tum.DZ.3)

#reorder samples
result.class1.PS.count.IMAC.tum.DZ.3 <- result.class1.PS.count.IMAC.tum.DZ.3[match(IMAC.tum.DZ.3.samples.order, result.class1.PS.count.IMAC.tum.DZ.3$Proteomics.Sample.Label_2),]
result.class1.PS.count.IMAC.tum.DZ.3

#calculate means
mean.result.class1.PS.count.IMAC.tum.DZ.3  <- mean(result.class1.PS.count.IMAC.tum.DZ.3$class1.numberPS); mean.result.class1.PS.count.IMAC.tum.DZ.3   #6604



## plot number of class 1 PS
ggplot(data = result.class1.PS.count.IMAC.tum.DZ.3 )+
  geom_col(aes(x=  Proteomics.Sample.Label_2 , y= class1.numberPS, fill=variant.group), color="black")+
  theme(legend.position="none") +
  scale_x_discrete(limits= IMAC.tum.DZ.3.samples.order)+
  scale_fill_viridis(discrete=T, option="inferno")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1.0, vjust = 0.5))+
  ylab("number PS")+
  xlab(NULL)+
  scale_y_continuous(breaks=seq(from=0, to=9200, by=400), limits = c(0,9200))+
  ggtitle("# PS class 1 per sample" )+
  theme(plot.title = element_text(hjust = 0.5, face= "bold"))+  
  theme(axis.text.x = element_text(face= "plain", colour="black", size=14))+
  theme(axis.text.y = element_text(face= "plain", colour="black", size=14))+
  geom_vline(xintercept = c(2.5, 4.5, 7.5, 10.5, 13.5, 15.5, 18.5, 21.5, 23.5, 25.5, 28.5, 30.5) , size=0.25, linetype = 2)


### plot number class 1: Ser, Thr, Tyr sites
glimpse(result.class1.PS.count.IMAC.tum.DZ.3)
result.class1.PS.count.IMAC.tum.DZ.3$variant.group.2

#normal plot
ggplot()+
  geom_col(data= result.class1.PS.count.IMAC.tum.DZ.3, aes(x=Proteomics.Sample.Label_2 , y=class1.nS), color="black", fill="blue")+
  geom_col(data= result.class1.PS.count.IMAC.tum.DZ.3, aes(x=Proteomics.Sample.Label_2 , y=class1.nT), color="black", fill="black")+
  geom_col(data= result.class1.PS.count.IMAC.tum.DZ.3, aes(x=Proteomics.Sample.Label_2 , y=class1.nY), color="black", fill="yellow")+
  scale_x_discrete(limits= c(IMAC.tum.DZ.3.samples.order))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1.0, vjust = 0.5))+
  ylab("number of sites")+
  xlab(NULL)+
  ggtitle(" # PS class 1 Tyr - yellow \n # PS class 1 Ser - blue  \n # PS class 1 Thr - black") +
  theme(plot.title = element_text(hjust = 0.5, face= "bold"))+  #center title
  theme(axis.text.x = element_text(face= "plain", colour="black", size=14))+
  theme(axis.text.y = element_text(face= "plain", colour="black", size=14))+
  geom_vline(xintercept = c(2.5, 4.5, 7.5, 10.5, 13.5, 15.5, 18.5, 21.5, 23.5, 25.5, 28.5, 30.5) , size=0.25, linetype = 2)

#plot with zoomed y-axis
ggplot()+
  geom_col(data= result.class1.PS.count.IMAC.tum.DZ.3, aes(x=Proteomics.Sample.Label_2 , y=class1.percentage.pS), color="black", fill="blue")+
  geom_col(data= result.class1.PS.count.IMAC.tum.DZ.3, aes(x=Proteomics.Sample.Label_2 , y=class1.percentage.pT), color="black", fill="black")+
  geom_col(data= result.class1.PS.count.IMAC.tum.DZ.3, aes(x=Proteomics.Sample.Label_2 , y=class1.percentage.pY), color="black", fill="yellow")+
  scale_x_discrete(limits= c(IMAC.tum.DZ.3.samples.order))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1.0, vjust = 0.5))+
  ylab("percentage")+
  xlab(NULL)+
  ggtitle(" # PS class 1 Tyr - yellow \n # PS class 1 Ser - blue  \n # PS class 1 Thr - black") +
  theme(plot.title = element_text(hjust = 0.5, face= "bold"))+  #center title
  theme(axis.text.x = element_text(face= "plain", colour="black", size=14))+
  theme(axis.text.y = element_text(face= "plain", colour="black", size=14))+
  coord_cartesian(ylim=c(0, 10.5))+
  geom_vline(xintercept = c(2.5, 4.5, 7.5, 10.5, 13.5, 15.5, 18.5, 21.5, 23.5, 25.5, 28.5, 30.5) , size=0.25, linetype = 2)



#### tumors IMAC MBR ON expand with PerseusPSexpand function (see above) ################################################################################################################
#### tumors IMAC MBR ON expand with PerseusPSexpand function (see above) ################################################################################################################

glimpse(IMAC.tum.DZ.3)

DZ.tum.IMAC_expanded <- PerseusPSexpand.17.02.21(IMAC.tum.DZ.3)


glimpse(DZ.tum.IMAC_expanded)



#remove rows that have all zero intensities
colnames(DZ.tum.IMAC_expanded)

DZ.tum.IMAC_expanded$ZEROcheckINT <- pbapply(DZ.tum.IMAC_expanded[, c(21:52)], MARGIN=1, function(x) sum(x)) #
table(DZ.tum.IMAC_expanded$ZEROcheckINT > 0)
nrow(DZ.tum.IMAC_expanded) #65781

DZ.tum.IMAC_expanded <- filter(DZ.tum.IMAC_expanded, ZEROcheckINT > 0)
nrow(DZ.tum.IMAC_expanded) #25464

#remove column ZEROcheckINT
DZ.tum.IMAC_expanded <- select(DZ.tum.IMAC_expanded,  -ZEROcheckINT)
glimpse(DZ.tum.IMAC_expanded)




# tumors IMAC  MBR ON normalization & histogram standard median shift version #########################################################################################################################################
# tumors IMAC  MBR ON normalization & histogram standard median shift version #########################################################################################################################################

glimpse(DZ.tum.IMAC_expanded)
colnames(DZ.tum.IMAC_expanded)

#log2 transform & replace -Inf with NA
colnames(DZ.tum.IMAC_expanded)
Log2Trafo <- DZ.tum.IMAC_expanded[, c(21:52)] #
Log2Trafo <- log2(Log2Trafo)
Log2Trafo[Log2Trafo  == -Inf] <- NA

glimpse(Log2Trafo)

DZ.tum.IMAC_expanded_log2 <- DZ.tum.IMAC_expanded
DZ.tum.IMAC_expanded_log2[, c(21:52)] <- Log2Trafo #

glimpse(DZ.tum.IMAC_expanded_log2, list.len=900)






# calculate median per column RawINT
colnames(DZ.tum.IMAC_expanded_log2)

glimpse(DZ.IMAC.tum.labels)
DZ.IMAC.tum.labels$Proteomics.Sample.Label_3 <- paste0("Intensity.", DZ.IMAC.tum.labels$Proteomics.Sample.Label)
DZ.IMAC.tum.labels <-  DZ.IMAC.tum.labels[match(IMAC.tum.DZ.3.samples.order,  DZ.IMAC.tum.labels$Proteomics.Sample.Label_2),]
DZ.IMAC.tum.labels$Proteomics.Sample.Label_3
DZ.IMAC.tum.labels$variant.group

#reorder
DZ.tum.IMAC_expanded_log2 <- DZ.tum.IMAC_expanded_log2 %>% select("id",                          "Proteins",                    "Positions.within.proteins",   "Leading.proteins",            "Protein",                    
                                                                  "Protein.names",               "Gene.names",                  "Fasta.headers",               "Localization.prob",           "Number.of.Phospho.(STY)",    
                                                                  "Amino.acid",                  "Sequence.window",             "Phospho.(STY).Probabilities", "Position.in.peptide",         "Reverse",                    
                                                                  "Potential.contaminant",       "Positions",                   "Position",                    "Peptide.IDs",                 "Mod..peptide.IDs",   
                                                                  
                                                                  "Intensity.A",  "Intensity.N",  "Intensity.F",  "Intensity.R",  "Intensity.D",  "Intensity.P",  "Intensity.AB", "Intensity.H",  "Intensity.T",  "Intensity.AD",
                                                                  "Intensity.C",  "Intensity.O",  "Intensity.AA", "Intensity.G",  "Intensity.S",  "Intensity.E",  "Intensity.Q",  "Intensity.AC", "Intensity.I",  "Intensity.U", 
                                                                  "Intensity.AE", "Intensity.J",  "Intensity.V",  "Intensity.K",  "Intensity.W",  "Intensity.B",  "Intensity.X",  "Intensity.AF", "Intensity.L",  "Intensity.Y", 
                                                                  "Intensity.M",  "Intensity.Z", 
                                                                  
                                                                  "PS_Multiplicity"
)

glimpse(DZ.tum.IMAC_expanded_log2)
colnames(DZ.tum.IMAC_expanded_log2)


ColMedianfindDF <- glimpse(DZ.tum.IMAC_expanded_log2[, c(21:52)]) #
vectormediancolumnsINT <- c()
for(i in 1:ncol(ColMedianfindDF)){
  print(i)
  #i=5
  temp <- median(ColMedianfindDF[, i], na.rm = TRUE); temp
  vectormediancolumnsINT <- append(vectormediancolumnsINT, temp)
}
vectormediancolumnsINT

median(vectormediancolumnsINT)


# data frame to plot medians as lines
Mdf <- data.frame(variable=colnames(DZ.tum.IMAC_expanded_log2[, c(21:52)]), median=vectormediancolumnsINT); Mdf

Mdf$variable2 <- c("1_Fgfr2",
                   "2_Fgfr2",
                   
                   "13_Fgfr2-Ate1",
                   "14_Fgfr2-Ate1",
                   
                   "7_Fgfr2-Bicc1",
                   "8_Fgfr2-Bicc1",
                   "9_Fgfr2-Bicc1",
                   
                   "17_Fgfr2-Tacc2",
                   "18_Fgfr2-Tacc2",
                   "19_Fgfr2-Tacc2",
                   
                   "4_Fgfr2-dE18",
                   "5_Fgfr2-dE18",
                   "6_Fgfr2-dE18",
                   
                   "15_Fgfr2-dE18-Ate1",
                   "16_Fgfr2-dE18-Ate1",
                   
                   "10_Fgfr2-dE18-Bicc1",
                   "11_Fgfr2-dE18-Bicc1",
                   "12_Fgfr2-dE18-Bicc1",
                   
                   "20_Fgfr2-dE18-Tacc2",
                   "21_Fgfr2-dE18-Tacc2",
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
                   
                   "31_Fgfr2-E18-C4" ,
                   "32_Fgfr2-E18-C4")

Mdf$group <- c("Fgfr2",            "Fgfr2",            "Fgfr2-Ate1",       "Fgfr2-Ate1",       "Fgfr2-Bicc1",      "Fgfr2-Bicc1",      "Fgfr2-Bicc1",      "Fgfr2-Tacc2",     
               "Fgfr2-Tacc2",      "Fgfr2-Tacc2",      "Fgfr2-dE18",       "Fgfr2-dE18",       "Fgfr2-dE18",       "Fgfr2-dE18-Ate1",  "Fgfr2-dE18-Ate1",  "Fgfr2-dE18-Bicc1",
               "Fgfr2-dE18-Bicc1", "Fgfr2-dE18-Bicc1", "Fgfr2-dE18-Tacc2", "Fgfr2-dE18-Tacc2", "Fgfr2-dE18-Tacc2", "Fgfr2-dE18-IGR1",  "Fgfr2-dE18-IGR1",  "Fgfr2-dE18-IGR2", 
               "Fgfr2-dE18-IGR2",  "Fgfr2-E18-C2",     "Fgfr2-E18-C2",     "Fgfr2-E18-C2",     "Fgfr2-E18-C3",     "Fgfr2-E18-C3",     "Fgfr2-E18-C4",     "Fgfr2-E18-C4")
Mdf

## plain density plot
DensityPlotData <- DZ.tum.IMAC_expanded_log2[, c(21:52)]

glimpse(DensityPlotData)

DensityPlotDataMelted <- reshape2::melt(DensityPlotData)
glimpse(DensityPlotDataMelted)


DensityPlot  <- ggplot(DensityPlotDataMelted, aes(x=value, colour = variable) ) + 
  geom_density() +
  theme(legend.position="bottom") +
  geom_vline(data=Mdf, aes(xintercept=median, color=variable))

DensityPlot 

# plot with additional sample names
glimpse(DensityPlotDataMelted)
DensityPlotDataMelted.2 <- merge(DensityPlotDataMelted, Mdf[, c("variable", "variable2", "group")], by.x = "variable", by.y = "variable", all.x = T, sort = F)
glimpse(DensityPlotDataMelted.2)

ggplot(DensityPlotDataMelted.2, aes(x=value, colour = variable2 ) ) + 
  geom_density() +
  scale_color_viridis(discrete=T, option="magma")+
  theme(legend.position="bottom")+
  ggtitle("mouse tumors DZ IMAC\n data distribution") +
  theme(plot.title = element_text(hjust = 0.5, face= "bold"))+  #center title
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
temp_for_median_shift <-  DZ.tum.IMAC_expanded_log2[, c(21:52)]
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
  #i=5
  temp <- median(ColMedianfindDF[, i], na.rm = TRUE); temp
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
DensityPlotCHECK  <- ggplot(DensityPlotDataCHECKMelted, aes(x=value, colour = variable) ) + 
  geom_density() +
  theme(legend.position="bottom") +
  geom_vline(data=MdfCHECK, aes(xintercept=median, color=variable))

DensityPlotCHECK


#plot check with additional sample names
DensityPlotDataCHECKMelted.2 <- merge(DensityPlotDataCHECKMelted, Mdf[, c("variable", "variable2", "group")], by.x = "variable", by.y = "variable", all.x = T, sort = F)
glimpse(DensityPlotDataCHECKMelted.2)

ggplot(DensityPlotDataCHECKMelted.2, aes(x=value, colour = variable2) ) + 
  geom_density() +
  scale_color_viridis(discrete=T, option="magma")+
  theme(legend.position="bottom")+
  ggtitle("mouse tumors DZ IMAC\n data distribution - normalized") +
  theme(plot.title = element_text(hjust = 0.5, face= "bold"))+  #center title
  theme(axis.text.x = element_text(face= "plain", colour="black", size=14))+
  theme(axis.text.y = element_text(face= "plain", colour="black", size=14))


### rename columns after normalization
newCOLNAMESnorm <- colnames(temp_for_median_shift); newCOLNAMESnorm

newCOLNAMESnorm <- paste0("Norm.",newCOLNAMESnorm); newCOLNAMESnorm

colnames(temp_for_median_shift) <- newCOLNAMESnorm
glimpse(temp_for_median_shift)


### add norm INT to original Dataframe 

DZ.tum.IMAC_expanded_log2_norm <- DZ.tum.IMAC_expanded_log2

DZ.tum.IMAC_expanded_log2_norm  <- cbind(DZ.tum.IMAC_expanded_log2_norm , temp_for_median_shift)
glimpse(DZ.tum.IMAC_expanded_log2_norm, list.len=999)

#prepare additional sample label
glimpse(DZ.IMAC.tum.labels)
DZ.IMAC.tum.labels$Proteomics.Sample.Label_4 <- paste0("Norm.Intensity.", DZ.IMAC.tum.labels$Proteomics.Sample.Label)
DZ.IMAC.tum.labels$Proteomics.Sample.Label_4

## box plot after normalization
data.for.boxplot.after.norm.DZ.pTyr.IP.MBR.ON <- DZ.tum.IMAC_expanded_log2_norm %>% select(                 
  "Norm.Intensity.A",  "Norm.Intensity.N",  "Norm.Intensity.F",  "Norm.Intensity.R",  "Norm.Intensity.D",  "Norm.Intensity.P",  "Norm.Intensity.AB",
  "Norm.Intensity.H",  "Norm.Intensity.T",  "Norm.Intensity.AD", "Norm.Intensity.C", "Norm.Intensity.O",  "Norm.Intensity.AA", "Norm.Intensity.G", 
  "Norm.Intensity.S",  "Norm.Intensity.E", "Norm.Intensity.Q",  "Norm.Intensity.AC", "Norm.Intensity.I",  "Norm.Intensity.U",  "Norm.Intensity.AE",
  "Norm.Intensity.J",  "Norm.Intensity.V",  "Norm.Intensity.K",  "Norm.Intensity.W",  "Norm.Intensity.B", "Norm.Intensity.X",  "Norm.Intensity.AF",
  "Norm.Intensity.L",  "Norm.Intensity.Y",  "Norm.Intensity.M",  "Norm.Intensity.Z" )

glimpse(data.for.boxplot.after.norm.DZ.pTyr.IP.MBR.ON)
colnames(data.for.boxplot.after.norm.DZ.pTyr.IP.MBR.ON) <- c("1_Fgfr2",
                                                             "2_Fgfr2",
                                                             
                                                             "13_Fgfr2-Ate1",
                                                             "14_Fgfr2-Ate1",
                                                             
                                                             "7_Fgfr2-Bicc1",
                                                             "8_Fgfr2-Bicc1",
                                                             "9_Fgfr2-Bicc1",
                                                             
                                                             "17_Fgfr2-Tacc2",
                                                             "18_Fgfr2-Tacc2",
                                                             "19_Fgfr2-Tacc2",
                                                             
                                                             "4_Fgfr2-dE18",
                                                             "5_Fgfr2-dE18",
                                                             "6_Fgfr2-dE18",
                                                             
                                                             "15_Fgfr2-dE18-Ate1",
                                                             "16_Fgfr2-dE18-Ate1",
                                                             
                                                             "10_Fgfr2-dE18-Bicc1",
                                                             "11_Fgfr2-dE18-Bicc1",
                                                             "12_Fgfr2-dE18-Bicc1",
                                                             
                                                             "20_Fgfr2-dE18-Tacc2",
                                                             "21_Fgfr2-dE18-Tacc2",
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
                                                             
                                                             "31_Fgfr2-E18-C4" ,
                                                             "32_Fgfr2-E18-C4")
glimpse(data.for.boxplot.after.norm.DZ.pTyr.IP.MBR.ON)

melted.data.for.boxplot.after.norm.DZ.pTyr.IP.MBR.ON <- reshape2::melt(data.for.boxplot.after.norm.DZ.pTyr.IP.MBR.ON)
glimpse(melted.data.for.boxplot.after.norm.DZ.pTyr.IP.MBR.ON )
unique(melted.data.for.boxplot.after.norm.DZ.pTyr.IP.MBR.ON$variable)

#plot
ggplot(data = melted.data.for.boxplot.after.norm.DZ.pTyr.IP.MBR.ON)+
  geom_boxplot(mapping = aes(x= variable, y = value), outlier.size = 0.5, notch = TRUE, notchwidth = 0.1)+ 
  scale_x_discrete(limits= c("1_Fgfr2",
                             "2_Fgfr2",
                             
                             "13_Fgfr2-Ate1",
                             "14_Fgfr2-Ate1",
                             
                             "7_Fgfr2-Bicc1",
                             "8_Fgfr2-Bicc1",
                             "9_Fgfr2-Bicc1",
                             
                             "17_Fgfr2-Tacc2",
                             "18_Fgfr2-Tacc2",
                             "19_Fgfr2-Tacc2",
                             
                             "4_Fgfr2-dE18",
                             "5_Fgfr2-dE18",
                             "6_Fgfr2-dE18",
                             
                             "15_Fgfr2-dE18-Ate1",
                             "16_Fgfr2-dE18-Ate1",
                             
                             "10_Fgfr2-dE18-Bicc1",
                             "11_Fgfr2-dE18-Bicc1",
                             "12_Fgfr2-dE18-Bicc1",
                             
                             "20_Fgfr2-dE18-Tacc2",
                             "21_Fgfr2-dE18-Tacc2",
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
                             
                             "31_Fgfr2-E18-C4" ,
                             "32_Fgfr2-E18-C4"))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1.0, vjust = 0.5))+
  theme(legend.position="none", legend.justification = "center")+
  ggtitle("tumors DZ IMAC - normalized data - PS level") +
  theme(plot.title = element_text(hjust = 0.5, face= "bold"))+
  theme(axis.text.x = element_text(face= "plain", colour="black", size=14))+
  theme(axis.text.y = element_text(face= "plain", colour="black", size=14))+
  xlab(NULL)+ 
  ylab("log2(Int.)") 


## box plot WITHOUT normalization
glimpse(DZ.tum.IMAC_expanded_log2_norm)

data.for.boxplot.after.norm.DZ.pTyr.IP.MBR.ON <- DZ.tum.IMAC_expanded_log2_norm %>% select(  "Intensity.A",  "Intensity.N",  "Intensity.F",  "Intensity.R",  "Intensity.D",  "Intensity.P",  "Intensity.AB", "Intensity.H",  "Intensity.T",  "Intensity.AD",
                                                                                             "Intensity.C",  "Intensity.O",  "Intensity.AA", "Intensity.G",  "Intensity.S",  "Intensity.E",  "Intensity.Q",  "Intensity.AC", "Intensity.I",  "Intensity.U", 
                                                                                             "Intensity.AE", "Intensity.J",  "Intensity.V",  "Intensity.K",  "Intensity.W",  "Intensity.B",  "Intensity.X",  "Intensity.AF", "Intensity.L",  "Intensity.Y", 
                                                                                             "Intensity.M",  "Intensity.Z")
glimpse(data.for.boxplot.after.norm.DZ.pTyr.IP.MBR.ON)
colnames(data.for.boxplot.after.norm.DZ.pTyr.IP.MBR.ON) <- c("1_Fgfr2",
                                                             "2_Fgfr2",
                                                             
                                                             "13_Fgfr2-Ate1",
                                                             "14_Fgfr2-Ate1",
                                                             
                                                             "7_Fgfr2-Bicc1",
                                                             "8_Fgfr2-Bicc1",
                                                             "9_Fgfr2-Bicc1",
                                                             
                                                             "17_Fgfr2-Tacc2",
                                                             "18_Fgfr2-Tacc2",
                                                             "19_Fgfr2-Tacc2",
                                                             
                                                             "4_Fgfr2-dE18",
                                                             "5_Fgfr2-dE18",
                                                             "6_Fgfr2-dE18",
                                                             
                                                             "15_Fgfr2-dE18-Ate1",
                                                             "16_Fgfr2-dE18-Ate1",
                                                             
                                                             "10_Fgfr2-dE18-Bicc1",
                                                             "11_Fgfr2-dE18-Bicc1",
                                                             "12_Fgfr2-dE18-Bicc1",
                                                             
                                                             "20_Fgfr2-dE18-Tacc2",
                                                             "21_Fgfr2-dE18-Tacc2",
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
                                                             
                                                             "31_Fgfr2-E18-C4" ,
                                                             "32_Fgfr2-E18-C4")
glimpse(data.for.boxplot.after.norm.DZ.pTyr.IP.MBR.ON)

melted.data.for.boxplot.after.norm.DZ.pTyr.IP.MBR.ON <- reshape2::melt(data.for.boxplot.after.norm.DZ.pTyr.IP.MBR.ON)
glimpse(melted.data.for.boxplot.after.norm.DZ.pTyr.IP.MBR.ON )
unique(melted.data.for.boxplot.after.norm.DZ.pTyr.IP.MBR.ON$variable)

ggplot(data = melted.data.for.boxplot.after.norm.DZ.pTyr.IP.MBR.ON)+
  geom_boxplot(mapping = aes(x= variable, y = value), outlier.size = 0.5, notch = TRUE, notchwidth = 0.1)+ 
  scale_x_discrete(limits= c("1_Fgfr2",
                             "2_Fgfr2",
                             
                             "13_Fgfr2-Ate1",
                             "14_Fgfr2-Ate1",
                             
                             "7_Fgfr2-Bicc1",
                             "8_Fgfr2-Bicc1",
                             "9_Fgfr2-Bicc1",
                             
                             "17_Fgfr2-Tacc2",
                             "18_Fgfr2-Tacc2",
                             "19_Fgfr2-Tacc2",
                             
                             "4_Fgfr2-dE18",
                             "5_Fgfr2-dE18",
                             "6_Fgfr2-dE18",
                             
                             "15_Fgfr2-dE18-Ate1",
                             "16_Fgfr2-dE18-Ate1",
                             
                             "10_Fgfr2-dE18-Bicc1",
                             "11_Fgfr2-dE18-Bicc1",
                             "12_Fgfr2-dE18-Bicc1",
                             
                             "20_Fgfr2-dE18-Tacc2",
                             "21_Fgfr2-dE18-Tacc2",
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
                             
                             "31_Fgfr2-E18-C4" ,
                             "32_Fgfr2-E18-C4"))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1.0, vjust = 0.5))+
  theme(legend.position="none", legend.justification = "center")+
  ggtitle("tumors DZ IMAC - NOT normalized data - PS level") +
  theme(plot.title = element_text(hjust = 0.5, face= "bold"))+  
  theme(axis.text.x = element_text(face= "plain", colour="black", size=14))+
  theme(axis.text.y = element_text(face= "plain", colour="black", size=14))+
  xlab(NULL)+ 
  ylab("log2(Int.)") 



#####################################################################################################################################################################
#### analysis via class 1 PS ##################################################################################################################################

glimpse(DZ.tum.IMAC_expanded_log2_norm)

#filter class 1
DZ.tum.IMAC_expanded_log2_normclass1.PS <- DZ.tum.IMAC_expanded_log2_norm %>% filter(Localization.prob>= 0.75)
glimpse(DZ.tum.IMAC_expanded_log2_normclass1.PS) #20196
length(unique(DZ.tum.IMAC_expanded_log2_normclass1.PS$Sequence.window)) #16930

# add additional column to better calculate number PS 
# take Proteins and add Psite without multiplicity info ###
for(i in 1:nrow(DZ.tum.IMAC_expanded_log2_normclass1.PS)){
  print(i)
  temp_prot_name <- DZ.tum.IMAC_expanded_log2_normclass1.PS$Proteins[i];temp_prot_name
  tempPSite <- DZ.tum.IMAC_expanded_log2_normclass1.PS$Positions[i]; tempPSite
  tempAA <- DZ.tum.IMAC_expanded_log2_normclass1.PS$Amino.acid[i];  tempAA
  DZ.tum.IMAC_expanded_log2_normclass1.PS$Prot.Name_AApos[i] <- paste0(temp_prot_name, "_", tempAA, tempPSite)
}
length(unique(DZ.tum.IMAC_expanded_log2_normclass1.PS$Prot.Name_AApos)) #16945 == the number of class 1 PS in exp

glimpse(DZ.tum.IMAC_expanded_log2_normclass1.PS)

#  prepare new column
# take genename and add Psite ###
for(i in 1:nrow(DZ.tum.IMAC_expanded_log2_normclass1.PS)){
  print(i)
  temp_gene_name <- DZ.tum.IMAC_expanded_log2_normclass1.PS$Gene.names[i];temp_gene_name
  tempPSite <- DZ.tum.IMAC_expanded_log2_normclass1.PS$Positions[i]; tempPSite
  tempAA <- DZ.tum.IMAC_expanded_log2_normclass1.PS$Amino.acid[i];  tempAA
  tempPSmultiplicity <- DZ.tum.IMAC_expanded_log2_normclass1.PS$PS_Multiplicity[i]; tempPSmultiplicity
  DZ.tum.IMAC_expanded_log2_normclass1.PS$Gene.Name_AApos[i] <- paste0(temp_gene_name, "_", tempAA, tempPSite, "x", tempPSmultiplicity)
}
glimpse(DZ.tum.IMAC_expanded_log2_normclass1.PS)


#add identifier for plotting reason e.g. PTMSEA later
glimpse(DZ.tum.IMAC_expanded_log2_normclass1.PS)

DZ.tum.IMAC_expanded_log2_normclass1.PS$ID.FR.all.C1.PS <- 1:nrow(DZ.tum.IMAC_expanded_log2_normclass1.PS)
glimpse(DZ.tum.IMAC_expanded_log2_normclass1.PS)

#save data
#save(DZ.tum.IMAC_expanded_log2_normclass1.PS, file="DZ.tum.IMAC_expanded_log2_normclass1.PS.Rdata")
#load("DZ.tum.IMAC_expanded_log2_normclass1.PS")
#glimpse(DZ.tum.IMAC_expanded_log2_normclass1.PS)



#############################################################################################################################################################
#############################################################################################################################################################
# correlation tumors IMAC class 1 PS

glimpse( DZ.tum.IMAC_expanded_log2_normclass1.PS)


#select columns and order
DF_for_correlation_plot.tum.IMAC.class1PS <-  DZ.tum.IMAC_expanded_log2_normclass1.PS %>% select(contains("Norm.Intensity.")               )
glimpse(DF_for_correlation_plot.tum.IMAC.class1PS)



##check if there are all NA rows 
DF_for_correlation_plot.tum.IMAC.class1PS$na.rows <- pbapply(DF_for_correlation_plot.tum.IMAC.class1PS, 1, function(x) sum(is.na(x)))
glimpse(DF_for_correlation_plot.tum.IMAC.class1PS)
table(DF_for_correlation_plot.tum.IMAC.class1PS$na.rows)

# change sample names and remove column na.rows
DF_for_correlation_plot.tum.IMAC.class1PS <- DF_for_correlation_plot.tum.IMAC.class1PS %>% select(-na.rows)
colnames(DF_for_correlation_plot.tum.IMAC.class1PS)
colnames(DF_for_correlation_plot.tum.IMAC.class1PS) <- c("1_Fgfr2",
                                                         "2_Fgfr2",
                                                         
                                                         "13_Fgfr2-Ate1",
                                                         "14_Fgfr2-Ate1",
                                                         
                                                         "7_Fgfr2-Bicc1",
                                                         "8_Fgfr2-Bicc1",
                                                         "9_Fgfr2-Bicc1",
                                                         
                                                         "17_Fgfr2-Tacc2",
                                                         "18_Fgfr2-Tacc2",
                                                         "19_Fgfr2-Tacc2",
                                                         
                                                         "4_Fgfr2-dE18",
                                                         "5_Fgfr2-dE18",
                                                         "6_Fgfr2-dE18",
                                                         
                                                         "15_Fgfr2-dE18-Ate1",
                                                         "16_Fgfr2-dE18-Ate1",
                                                         
                                                         "10_Fgfr2-dE18-Bicc1",
                                                         "11_Fgfr2-dE18-Bicc1",
                                                         "12_Fgfr2-dE18-Bicc1",
                                                         
                                                         "20_Fgfr2-dE18-Tacc2",
                                                         "21_Fgfr2-dE18-Tacc2",
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
                                                         
                                                         "31_Fgfr2-E18-C4" ,
                                                         "32_Fgfr2-E18-C4")


glimpse(DF_for_correlation_plot.tum.IMAC.class1PS)


#change sample names again
tum.IMAC.name.change <- tibble(orig.col.name = colnames(DF_for_correlation_plot.tum.IMAC.class1PS))
tum.IMAC.name.change$new.col.name.1 <- pbsapply(tum.IMAC.name.change$orig.col.name, function(x) unlist(str_split(string=x, pattern="_"))[2]   )
tum.IMAC.name.change$new.col.name.2 <- paste0(tum.IMAC.name.change$new.col.name.1, ".", c(1, 2,   1, 2,   1, 2, 3,   1, 2, 3,    1, 2, 3,   1, 2,   1, 2, 3,    1, 2, 3,    1, 2,  1, 2,   1, 2, 3 ,   1, 2,   1, 2))
print(tum.IMAC.name.change , n=45)
colnames(DF_for_correlation_plot.tum.IMAC.class1PS) <- tum.IMAC.name.change$new.col.name.2
glimpse(DF_for_correlation_plot.tum.IMAC.class1PS)


#correlate
corr2 <- cor(DF_for_correlation_plot.tum.IMAC.class1PS, method = "pearson", use = "na.or.complete")
glimpse(corr2)
head(corr2)
min(corr2) #0.6705535
max(corr2)
round(min(corr2), 2)

# prepare annotation data
annot_df_for_heatmap <- data.frame(samples = colnames(DF_for_correlation_plot.tum.IMAC.class1PS),
                                   group = c("Fgfr2",            "Fgfr2",            "Fgfr2-Ate1",       "Fgfr2-Ate1",       "Fgfr2-Bicc1",      "Fgfr2-Bicc1",      "Fgfr2-Bicc1",      "Fgfr2-Tacc2",     
                                             "Fgfr2-Tacc2",      "Fgfr2-Tacc2",      "Fgfr2-dE18",       "Fgfr2-dE18",       "Fgfr2-dE18",       "Fgfr2-dE18-Ate1",  "Fgfr2-dE18-Ate1",  "Fgfr2-dE18-Bicc1",
                                             "Fgfr2-dE18-Bicc1", "Fgfr2-dE18-Bicc1", "Fgfr2-dE18-Tacc2", "Fgfr2-dE18-Tacc2", "Fgfr2-dE18-Tacc2", "Fgfr2-dE18-IGR1",  "Fgfr2-dE18-IGR1",  "Fgfr2-dE18-IGR2", 
                                             "Fgfr2-dE18-IGR2",  "Fgfr2-E18-C2",     "Fgfr2-E18-C2",     "Fgfr2-E18-C2",     "Fgfr2-E18-C3",     "Fgfr2-E18-C3",     "Fgfr2-E18-C4",     "Fgfr2-E18-C4") #
                                  )

annot_df_for_heatmap 


# make shorter annotation data frame
annot_df_for_heatmap.short <- data.frame( 
  group  = annot_df_for_heatmap$group
)

annot_df_for_heatmap.short

# define colors for annotation bar
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
FR.heatmap.colors.2 <-colorRamp2(c(min(corr2) , 1.0), c("grey90", "grey10"))



#pdf("DZ.tumors.IMAC.correlation.pdf", width=25/2.54, height=20/2.54, useDingbats=FALSE)
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



############################################################################################################################################################
############################################################################################################################################################

### heatmap for all Fgfr2 sites in dataset for all groups with collapsing _1/_2/_3 via sum 
### PS numbers also for NP_963895.2


glimpse(DZ.tum.IMAC_expanded_log2_normclass1.PS)

# filter for protein candidate of interest
HM.plotdata.IMAC.tum <- DZ.tum.IMAC_expanded_log2_normclass1.PS %>% filter(stringr::str_detect(Proteins , "P21803")) #

glimpse(HM.plotdata.IMAC.tum)


#pick columns of interest and rename
HM.plotdata.IMAC.tum <- HM.plotdata.IMAC.tum %>% select("Gene.Name_AApos", 
                                                        "Norm.Intensity.A",  "Norm.Intensity.N",  "Norm.Intensity.F",  "Norm.Intensity.R",  "Norm.Intensity.D",  "Norm.Intensity.P",  "Norm.Intensity.AB",
                                                        "Norm.Intensity.H",  "Norm.Intensity.T",  "Norm.Intensity.AD", "Norm.Intensity.C", "Norm.Intensity.O",  "Norm.Intensity.AA", "Norm.Intensity.G", 
                                                        "Norm.Intensity.S",  "Norm.Intensity.E", "Norm.Intensity.Q",  "Norm.Intensity.AC", "Norm.Intensity.I",  "Norm.Intensity.U",  "Norm.Intensity.AE",
                                                        "Norm.Intensity.J",  "Norm.Intensity.V",  "Norm.Intensity.K",  "Norm.Intensity.W",  "Norm.Intensity.B", "Norm.Intensity.X",  "Norm.Intensity.AF",
                                                        "Norm.Intensity.L",  "Norm.Intensity.Y",  "Norm.Intensity.M",  "Norm.Intensity.Z" ,
                                                        "Proteins",
                                                        "Positions.within.proteins",
                                                        "Amino.acid",
                                                        "Gene.names" ,
                                                        "PS_Multiplicity")
glimpse(HM.plotdata.IMAC.tum)

colnames(HM.plotdata.IMAC.tum) <- c("Gene.Name_AApos", 
                                    "1_Fgfr2",
                                    "2_Fgfr2",
                                    
                                    "13_Fgfr2-Ate1",
                                    "14_Fgfr2-Ate1",
                                    
                                    "7_Fgfr2-Bicc1",
                                    "8_Fgfr2-Bicc1",
                                    "9_Fgfr2-Bicc1",
                                    
                                    "17_Fgfr2-Tacc2",
                                    "18_Fgfr2-Tacc2",
                                    "19_Fgfr2-Tacc2",
                                    
                                    "4_Fgfr2-dE18",
                                    "5_Fgfr2-dE18",
                                    "6_Fgfr2-dE18",
                                    
                                    "15_Fgfr2-dE18-Ate1",
                                    "16_Fgfr2-dE18-Ate1",
                                    
                                    "10_Fgfr2-dE18-Bicc1",
                                    "11_Fgfr2-dE18-Bicc1",
                                    "12_Fgfr2-dE18-Bicc1",
                                    
                                    "20_Fgfr2-dE18-Tacc2",
                                    "21_Fgfr2-dE18-Tacc2",
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
                                    
                                    "31_Fgfr2-E18-C4" ,
                                    "32_Fgfr2-E18-C4",
                                    "Proteins",
                                    "Positions.within.proteins",
                                    "Amino.acid",
                                    "Gene.names" ,
                                    "PS_Multiplicity")
glimpse(HM.plotdata.IMAC.tum)

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
for(i in 1:nrow(HM.plotdata.IMAC.tum)){
  print(i)
  HM.plotdata.IMAC.tum$isoform.1.corrected.PS.position[i] <- position.in.mq.string(prot.of.interest ="P21803", #isoform 1
                                                                                   Gene.names=HM.plotdata.IMAC.tum$Gene.names[i],
                                                                                   Amino.acid=HM.plotdata.IMAC.tum$Amino.acid[i], 
                                                                                   mq.position.string=HM.plotdata.IMAC.tum$Positions.within.proteins[i],
                                                                                   mq.string = HM.plotdata.IMAC.tum$Proteins[i],
                                                                                   PS_Multiplicity=HM.plotdata.IMAC.tum$PS_Multiplicity[i])
}
HM.plotdata.IMAC.tum$isoform.1.corrected.Gene.Name_AApos = paste0(HM.plotdata.IMAC.tum$Gene.names, "_", HM.plotdata.IMAC.tum$Amino.acid, HM.plotdata.IMAC.tum$isoform.1.corrected.PS.position, "x", HM.plotdata.IMAC.tum$PS_Multiplicity)
glimpse(HM.plotdata.IMAC.tum)

# isoform 2 ="P21803-2"
for(i in 1:nrow(HM.plotdata.IMAC.tum)){
  print(i)
  HM.plotdata.IMAC.tum$isoform.2.corrected.PS.position[i] <- position.in.mq.string(prot.of.interest ="P21803-2", #isoform 2
                                                                                   Gene.names=HM.plotdata.IMAC.tum$Gene.names[i],
                                                                                   Amino.acid=HM.plotdata.IMAC.tum$Amino.acid[i], 
                                                                                   mq.position.string=HM.plotdata.IMAC.tum$Positions.within.proteins[i],
                                                                                   mq.string = HM.plotdata.IMAC.tum$Proteins[i],
                                                                                   PS_Multiplicity=HM.plotdata.IMAC.tum$PS_Multiplicity[i])
}
HM.plotdata.IMAC.tum$isoform.2.corrected.Gene.Name_AApos = paste0(HM.plotdata.IMAC.tum$Gene.names, "_", HM.plotdata.IMAC.tum$Amino.acid, HM.plotdata.IMAC.tum$isoform.2.corrected.PS.position, "x", HM.plotdata.IMAC.tum$PS_Multiplicity)
glimpse(HM.plotdata.IMAC.tum)


#prepare helper data frame for later order 
HM.plotdata.IMAC.tum.melted <- reshape2::melt(HM.plotdata.IMAC.tum)
glimpse(HM.plotdata.IMAC.tum.melted)
colnames(HM.plotdata.IMAC.tum.melted) <- c("Gene.Name_AApos", "Proteins","Positions.within.proteins", "Amino.acid","Gene.names", "PS_Multiplicity", "isoform.1.corrected.PS.position", "isoform.1.corrected.Gene.Name_AApos","isoform.2.corrected.PS.position", "isoform.2.corrected.Gene.Name_AApos",  "variable","log2.int")
glimpse(HM.plotdata.IMAC.tum.melted)

order.df.HM.plotdata.IMAC.tum.melted <- unique(tibble(Gene.Name_AApos = HM.plotdata.IMAC.tum.melted$Gene.Name_AApos,
                                                      Proteins = HM.plotdata.IMAC.tum.melted$Proteins, 
                                                      Positions.within.proteins = HM.plotdata.IMAC.tum.melted$Positions.within.proteins,
                                                      Amino.acid = HM.plotdata.IMAC.tum.melted$Amino.acid,
                                                      Gene.names = HM.plotdata.IMAC.tum.melted$Gene.names, 
                                                      PS_Multiplicity =HM.plotdata.IMAC.tum.melted$PS_Multiplicity,
                                                      isoform.1.corrected.PS.position =HM.plotdata.IMAC.tum.melted$isoform.1.corrected.PS.position,
                                                      isoform.1.corrected.Gene.Name_AApos = paste0(HM.plotdata.IMAC.tum.melted$Gene.names, "_", HM.plotdata.IMAC.tum.melted$Amino.acid, HM.plotdata.IMAC.tum.melted$isoform.1.corrected.PS.position, "x", HM.plotdata.IMAC.tum.melted$PS_Multiplicity),
                                                      isoform.2.corrected.PS.position =HM.plotdata.IMAC.tum.melted$isoform.2.corrected.PS.position,
                                                      isoform.2.corrected.Gene.Name_AApos = paste0(HM.plotdata.IMAC.tum.melted$Gene.names, "_", HM.plotdata.IMAC.tum.melted$Amino.acid, HM.plotdata.IMAC.tum.melted$isoform.2.corrected.PS.position, "x", HM.plotdata.IMAC.tum.melted$PS_Multiplicity)
))
glimpse(order.df.HM.plotdata.IMAC.tum.melted)


### for better visualization/plotting => collapse _1/_2/_3 via  sum
### note: sum is based on _1/_2/_3 intensities in dataframe at hand => unfiltered data might contain more _1/_2/_3 PS that are not included in the sum here

HM.plotdata.IMAC.tum.collapsed <- HM.plotdata.IMAC.tum
glimpse(HM.plotdata.IMAC.tum.collapsed )

HM.plotdata.IMAC.tum.collapsed$isoform.1.corrected.Gene.Name_AApos.collapsed <- pbsapply(HM.plotdata.IMAC.tum.collapsed$isoform.1.corrected.Gene.Name_AApos, function(x) str_remove_all(string=x, pattern="x_\\d") )
HM.plotdata.IMAC.tum.collapsed$isoform.2.corrected.Gene.Name_AApos.collapsed <- pbsapply(HM.plotdata.IMAC.tum.collapsed$isoform.2.corrected.Gene.Name_AApos, function(x) str_remove_all(string=x, pattern="x_\\d") )
glimpse(HM.plotdata.IMAC.tum.collapsed)
## unlog, do collapsing, relog
HM.plotdata.IMAC.tum.collapsed [, c("1_Fgfr2",
                                    "2_Fgfr2",
                                    
                                    "13_Fgfr2-Ate1",
                                    "14_Fgfr2-Ate1",
                                    
                                    "7_Fgfr2-Bicc1",
                                    "8_Fgfr2-Bicc1",
                                    "9_Fgfr2-Bicc1",
                                    
                                    "17_Fgfr2-Tacc2",
                                    "18_Fgfr2-Tacc2",
                                    "19_Fgfr2-Tacc2",
                                    
                                    "4_Fgfr2-dE18",
                                    "5_Fgfr2-dE18",
                                    "6_Fgfr2-dE18",
                                    
                                    "15_Fgfr2-dE18-Ate1",
                                    "16_Fgfr2-dE18-Ate1",
                                    
                                    "10_Fgfr2-dE18-Bicc1",
                                    "11_Fgfr2-dE18-Bicc1",
                                    "12_Fgfr2-dE18-Bicc1",
                                    
                                    "20_Fgfr2-dE18-Tacc2",
                                    "21_Fgfr2-dE18-Tacc2",
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
                                    
                                    "31_Fgfr2-E18-C4" ,
                                    "32_Fgfr2-E18-C4")] <- 2^HM.plotdata.IMAC.tum.collapsed [, c("1_Fgfr2",
                                                                                                 "2_Fgfr2",
                                                                                                 
                                                                                                 "13_Fgfr2-Ate1",
                                                                                                 "14_Fgfr2-Ate1",
                                                                                                 
                                                                                                 "7_Fgfr2-Bicc1",
                                                                                                 "8_Fgfr2-Bicc1",
                                                                                                 "9_Fgfr2-Bicc1",
                                                                                                 
                                                                                                 "17_Fgfr2-Tacc2",
                                                                                                 "18_Fgfr2-Tacc2",
                                                                                                 "19_Fgfr2-Tacc2",
                                                                                                 
                                                                                                 "4_Fgfr2-dE18",
                                                                                                 "5_Fgfr2-dE18",
                                                                                                 "6_Fgfr2-dE18",
                                                                                                 
                                                                                                 "15_Fgfr2-dE18-Ate1",
                                                                                                 "16_Fgfr2-dE18-Ate1",
                                                                                                 
                                                                                                 "10_Fgfr2-dE18-Bicc1",
                                                                                                 "11_Fgfr2-dE18-Bicc1",
                                                                                                 "12_Fgfr2-dE18-Bicc1",
                                                                                                 
                                                                                                 "20_Fgfr2-dE18-Tacc2",
                                                                                                 "21_Fgfr2-dE18-Tacc2",
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
                                                                                                 
                                                                                                 "31_Fgfr2-E18-C4" ,
                                                                                                 "32_Fgfr2-E18-C4")] 
glimpse(HM.plotdata.IMAC.tum.collapsed )
##collapse
HM.plotdata.IMAC.tum.collapsed.2  <- HM.plotdata.IMAC.tum.collapsed %>% 
  group_by(isoform.1.corrected.Gene.Name_AApos.collapsed, isoform.2.corrected.Gene.Name_AApos.collapsed) %>% 
  summarise_at(c("1_Fgfr2",
                 "2_Fgfr2",
                 
                 "13_Fgfr2-Ate1",
                 "14_Fgfr2-Ate1",
                 
                 "7_Fgfr2-Bicc1",
                 "8_Fgfr2-Bicc1",
                 "9_Fgfr2-Bicc1",
                 
                 "17_Fgfr2-Tacc2",
                 "18_Fgfr2-Tacc2",
                 "19_Fgfr2-Tacc2",
                 
                 "4_Fgfr2-dE18",
                 "5_Fgfr2-dE18",
                 "6_Fgfr2-dE18",
                 
                 "15_Fgfr2-dE18-Ate1",
                 "16_Fgfr2-dE18-Ate1",
                 
                 "10_Fgfr2-dE18-Bicc1",
                 "11_Fgfr2-dE18-Bicc1",
                 "12_Fgfr2-dE18-Bicc1",
                 
                 "20_Fgfr2-dE18-Tacc2",
                 "21_Fgfr2-dE18-Tacc2",
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
                 
                 "31_Fgfr2-E18-C4" ,
                 "32_Fgfr2-E18-C4"), sum, na.rm = TRUE)
glimpse(HM.plotdata.IMAC.tum.collapsed.2)
HM.plotdata.IMAC.tum.collapsed.2[HM.plotdata.IMAC.tum.collapsed.2  == 0] <- NA
##relog
HM.plotdata.IMAC.tum.collapsed.2[c("1_Fgfr2",
                                   "2_Fgfr2",
                                   
                                   "13_Fgfr2-Ate1",
                                   "14_Fgfr2-Ate1",
                                   
                                   "7_Fgfr2-Bicc1",
                                   "8_Fgfr2-Bicc1",
                                   "9_Fgfr2-Bicc1",
                                   
                                   "17_Fgfr2-Tacc2",
                                   "18_Fgfr2-Tacc2",
                                   "19_Fgfr2-Tacc2",
                                   
                                   "4_Fgfr2-dE18",
                                   "5_Fgfr2-dE18",
                                   "6_Fgfr2-dE18",
                                   
                                   "15_Fgfr2-dE18-Ate1",
                                   "16_Fgfr2-dE18-Ate1",
                                   
                                   "10_Fgfr2-dE18-Bicc1",
                                   "11_Fgfr2-dE18-Bicc1",
                                   "12_Fgfr2-dE18-Bicc1",
                                   
                                   "20_Fgfr2-dE18-Tacc2",
                                   "21_Fgfr2-dE18-Tacc2",
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
                                   
                                   "31_Fgfr2-E18-C4" ,
                                   "32_Fgfr2-E18-C4")] <- log2(HM.plotdata.IMAC.tum.collapsed.2[c("1_Fgfr2",
                                                                                                  "2_Fgfr2",
                                                                                                  
                                                                                                  "13_Fgfr2-Ate1",
                                                                                                  "14_Fgfr2-Ate1",
                                                                                                  
                                                                                                  "7_Fgfr2-Bicc1",
                                                                                                  "8_Fgfr2-Bicc1",
                                                                                                  "9_Fgfr2-Bicc1",
                                                                                                  
                                                                                                  "17_Fgfr2-Tacc2",
                                                                                                  "18_Fgfr2-Tacc2",
                                                                                                  "19_Fgfr2-Tacc2",
                                                                                                  
                                                                                                  "4_Fgfr2-dE18",
                                                                                                  "5_Fgfr2-dE18",
                                                                                                  "6_Fgfr2-dE18",
                                                                                                  
                                                                                                  "15_Fgfr2-dE18-Ate1",
                                                                                                  "16_Fgfr2-dE18-Ate1",
                                                                                                  
                                                                                                  "10_Fgfr2-dE18-Bicc1",
                                                                                                  "11_Fgfr2-dE18-Bicc1",
                                                                                                  "12_Fgfr2-dE18-Bicc1",
                                                                                                  
                                                                                                  "20_Fgfr2-dE18-Tacc2",
                                                                                                  "21_Fgfr2-dE18-Tacc2",
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
                                                                                                  
                                                                                                  "31_Fgfr2-E18-C4" ,
                                                                                                  "32_Fgfr2-E18-C4")] )
glimpse(HM.plotdata.IMAC.tum.collapsed.2)


### reshape and reorder, define order for plot
HM.plotdata.IMAC.tum.collapsed.melted <- reshape2::melt(HM.plotdata.IMAC.tum.collapsed.2)
glimpse(HM.plotdata.IMAC.tum.collapsed.melted)
colnames(HM.plotdata.IMAC.tum.collapsed.melted ) <- c("isoform.1.corrected.Gene.Name_AApos.collapsed","isoform.2.corrected.Gene.Name_AApos.collapsed", "variable","log2.int")
glimpse(HM.plotdata.IMAC.tum.collapsed.melted )


#define order for plot, get gene names, corrected positions, aminco acid etc. and order accor4ding to number
order.df.HM.plotdata.IMAC.tum.collapsed.melted <- unique(tibble(isoform.1.corrected.Gene.Name_AApos.collapsed = HM.plotdata.IMAC.tum.collapsed.melted$isoform.1.corrected.Gene.Name_AApos.collapsed,
                                                                isoform.2.corrected.Gene.Name_AApos.collapsed = HM.plotdata.IMAC.tum.collapsed.melted$isoform.2.corrected.Gene.Name_AApos.collapsed))


order.df.HM.plotdata.IMAC.tum.collapsed.melted$Gene <- pbsapply(order.df.HM.plotdata.IMAC.tum.collapsed.melted$isoform.1.corrected.Gene.Name_AApos.collapsed , function(x) unlist(str_split(string=x,  pattern="_") )[1])

order.df.HM.plotdata.IMAC.tum.collapsed.melted$isoform.1.corrected.AAPos <- pbsapply(order.df.HM.plotdata.IMAC.tum.collapsed.melted$isoform.1.corrected.Gene.Name_AApos.collapsed, function(x) unlist(str_split(string=x,  pattern="_") )[2])
order.df.HM.plotdata.IMAC.tum.collapsed.melted$isoform.2.corrected.AAPos <- pbsapply(order.df.HM.plotdata.IMAC.tum.collapsed.melted$isoform.2.corrected.Gene.Name_AApos.collapsed, function(x) unlist(str_split(string=x,  pattern="_") )[2])

order.df.HM.plotdata.IMAC.tum.collapsed.melted$isoform.corrected.Amino.acid <- pbsapply(order.df.HM.plotdata.IMAC.tum.collapsed.melted$isoform.1.corrected.AAPos, function(x) unlist((str_extract_all(string=x,  pattern="[A-Z]") )))

order.df.HM.plotdata.IMAC.tum.collapsed.melted$isoform.1.corrected.Pos <- pbsapply(order.df.HM.plotdata.IMAC.tum.collapsed.melted$isoform.1.corrected.AAPos, function(x) unlist((str_extract_all(string=x,  pattern="[0-9]+") )))
order.df.HM.plotdata.IMAC.tum.collapsed.melted$isoform.2.corrected.Pos <- pbsapply(order.df.HM.plotdata.IMAC.tum.collapsed.melted$isoform.2.corrected.AAPos, function(x) unlist((str_extract_all(string=x,  pattern="[0-9]+") )))

order.df.HM.plotdata.IMAC.tum.collapsed.melted <- order.df.HM.plotdata.IMAC.tum.collapsed.melted %>% arrange(desc(isoform.2.corrected.Pos, isoform.2.corrected.Amino.acid))
order.df.HM.plotdata.IMAC.tum.collapsed.melted

isoformcorrected.collapes.order.29.11.21 <- order.df.HM.plotdata.IMAC.tum.collapsed.melted$isoform.2.corrected.Gene.Name_AApos.collapsed
isoformcorrected.collapes.order.29.11.21

glimpse(HM.plotdata.IMAC.tum.collapsed.melted)

#plot the heatmap with site numbers for P21803-2
ggplot(
  HM.plotdata.IMAC.tum.collapsed.melted , aes(x=variable, y=isoform.2.corrected.Gene.Name_AApos.collapsed, fill=log2.int)) +
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
  geom_vline(xintercept = c(2.5, 4.5, 7.5, 10.5, 13.5, 15.5, 18.5, 21.5, 23.5, 25.5, 28.5, 30.5) , size=1.0, linetype = "solid", color="white")

### add corresponding human FGFR2 phosphosite numbers to the heatmap plot ######################################################

glimpse(HM.plotdata.IMAC.tum.collapsed.2)
order.df.HM.plotdata.IMAC.tum.collapsed.melted

# online Clustal Omega alignment showed that P21803-2 (707 aa) and NP_963895.2 (726 aa) are very similar. The difference is N-terminal: 
# 19aa  (MGLPSTWRYGRGPGIGTVT) are added for  NP_963895.2. Everything else is identical.
# Thus P21803-2 position +19 should be the positions for NP_963895.2

#prepare additional data frame and calculate the new site numbers
HM.plotdata.IMAC.tum.collapsed.3 <- HM.plotdata.IMAC.tum.collapsed.2

HM.plotdata.IMAC.tum.collapsed.3$isoform.1.corrected.AApos <- pbsapply(HM.plotdata.IMAC.tum.collapsed.3$isoform.1.corrected.Gene.Name_AApos.collapsed, function(x) unlist(str_split(string=x, pattern="_"))[2])
HM.plotdata.IMAC.tum.collapsed.3$isoform.1.corrected.pos <- pbsapply(HM.plotdata.IMAC.tum.collapsed.3$isoform.1.corrected.AApos, function(x) unlist(str_extract_all(string=x, pattern="[0-9]+")))

HM.plotdata.IMAC.tum.collapsed.3$isoform.2.corrected.AApos <- pbsapply(HM.plotdata.IMAC.tum.collapsed.3$isoform.2.corrected.Gene.Name_AApos.collapsed, function(x) unlist(str_split(string=x, pattern="_"))[2])
HM.plotdata.IMAC.tum.collapsed.3$isoform.2.corrected.pos <- pbsapply(HM.plotdata.IMAC.tum.collapsed.3$isoform.2.corrected.AApos, function(x) unlist(str_extract_all(string=x, pattern="[0-9]+")))
glimpse(HM.plotdata.IMAC.tum.collapsed.3)


HM.plotdata.IMAC.tum.collapsed.3$pos.NP_963895.2 <- as.numeric(HM.plotdata.IMAC.tum.collapsed.3$isoform.2.corrected.pos) +19
glimpse(HM.plotdata.IMAC.tum.collapsed.3)

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
HM.plotdata.IMAC.tum.collapsed.4 <- merge(HM.plotdata.IMAC.tum.collapsed.3, human.PSP.site.number.lookup.2[, c("site.Fgfr2.mouse.iso1", "site.FGFR2.human.iso1" )], by.x = "isoform.1.corrected.AApos", by.y = "site.Fgfr2.mouse.iso1", all.x = T, all.y = F, sort = F)
HM.plotdata.IMAC.tum.collapsed.4$Gene <- pbsapply(HM.plotdata.IMAC.tum.collapsed.4$isoform.1.corrected.Gene.Name_AApos.collapsed, function(x) unlist(str_split(string=x,  pattern="_") )[1])
HM.plotdata.IMAC.tum.collapsed.4 <- HM.plotdata.IMAC.tum.collapsed.4  %>% arrange(HM.plotdata.IMAC.tum.collapsed.4$isoform.1.corrected.pos)

glimpse(HM.plotdata.IMAC.tum.collapsed.4)

#some human sites have to be entered manually since not all sites are in PSP yet (looked up in online Clustal Omega alignment)
#they correspond mostly to numbers from P21803-1
unname(HM.plotdata.IMAC.tum.collapsed.4$isoform.1.corrected.AApos)
HM.plotdata.IMAC.tum.collapsed.4$site.FGFR2.human.iso1

HM.plotdata.IMAC.tum.collapsed.4$AA.site.FGFR2.human.iso1_manual.change <- c("T419", "T429", "S452", "S453", "Y561", "Y586", "S587", "Y608", "Y656", "Y657")
HM.plotdata.IMAC.tum.collapsed.4$site.FGFR2.human.iso1_manual.change <- pbsapply(HM.plotdata.IMAC.tum.collapsed.4$AA.site.FGFR2.human.iso1_manual.change, function(x) str_remove_all(string=x, pattern="[A-Z]"))
glimpse(HM.plotdata.IMAC.tum.collapsed.4)

#make.new.plot.label
#e.g. Fgfr2_T419|305|324|419 => 1st number: P21803-1; 2nd number: P21803-2; 3rd number: NP_963895.2; 4th number: P21802-1
HM.plotdata.IMAC.tum.collapsed.4$combined.plot.label <- paste0(HM.plotdata.IMAC.tum.collapsed.4$isoform.1.corrected.Gene.Name_AApos.collapsed, "|", 
                                                               HM.plotdata.IMAC.tum.collapsed.4$isoform.2.corrected.pos, "|",  
                                                               HM.plotdata.IMAC.tum.collapsed.4$pos.NP_963895.2, "|", 
                                                               HM.plotdata.IMAC.tum.collapsed.4$site.FGFR2.human.iso1_manual.change)
glimpse(HM.plotdata.IMAC.tum.collapsed.4)

#order for plot
isoform.corrected.collapsed.order.02.12.2021 <- rev(HM.plotdata.IMAC.tum.collapsed.4$combined.plot.label)


#change samples names and reorder for plot
name.change.DZ.tum.IMAC.dia <- tibble(orig.col.name = colnames(HM.plotdata.IMAC.tum.collapsed.4 %>% select(`1_Fgfr2` : `32_Fgfr2-E18-C4`)) )
name.change.DZ.tum.IMAC.dia$new.col.name.1 <- pbsapply(name.change.DZ.tum.IMAC.dia$orig.col.name, function(x) unlist(str_split(string=x, pattern="_"))[2]   )
name.change.DZ.tum.IMAC.dia$new.col.name.2 <- paste0(name.change.DZ.tum.IMAC.dia$new.col.name.1, ".", c(1, 2,   1, 2,   1, 2, 3,   1, 2, 3,    1, 2, 3,   1, 2,   1, 2, 3,    1, 2, 3,    1, 2,  1, 2,   1, 2, 3 ,   1, 2,   1, 2))
print(name.change.DZ.tum.IMAC.dia, n=32)

colnames(HM.plotdata.IMAC.tum.collapsed.4)
colnames(HM.plotdata.IMAC.tum.collapsed.4) <- c("isoform.1.corrected.AApos","isoform.1.corrected.Gene.Name_AApos.collapsed",
                                                "isoform.2.corrected.Gene.Name_AApos.collapsed",
                                                name.change.DZ.tum.IMAC.dia$new.col.name.2,                      
                                                "isoform.1.corrected.pos","isoform.2.corrected.AApos",
                                                "isoform.2.corrected.pos","pos.NP_963895.2","site.FGFR2.human.iso1","Gene","AA.site.FGFR2.human.iso1_manual.change",       
                                                "site.FGFR2.human.iso1_manual.change","combined.plot.label")
glimpse(HM.plotdata.IMAC.tum.collapsed.4)

# reshape for plot
melted.HM.plotdata.IMAC.tum.collapsed.4 <- reshape2::melt(HM.plotdata.IMAC.tum.collapsed.4 %>% select(combined.plot.label, "Fgfr2.1" : "Fgfr2-E18-C4.2"))
glimpse(melted.HM.plotdata.IMAC.tum.collapsed.4)

custom.sample.order.DZ.tum.IMAC <- c("Fgfr2.1", "Fgfr2.2", 
                                     "Fgfr2-dE18.1", "Fgfr2-dE18.2", "Fgfr2-dE18.3", 
                                     "Fgfr2-Bicc1.1", "Fgfr2-Bicc1.2" , "Fgfr2-Bicc1.3", 
                                     "Fgfr2-dE18-Bicc1.1", "Fgfr2-dE18-Bicc1.2", "Fgfr2-dE18-Bicc1.3", 
                                     "Fgfr2-Ate1.1", "Fgfr2-Ate1.2",
                                     "Fgfr2-dE18-Ate1.1", "Fgfr2-dE18-Ate1.2",
                                     "Fgfr2-Tacc2.1", "Fgfr2-Tacc2.2", "Fgfr2-Tacc2.3",
                                     "Fgfr2-dE18-Tacc2.1", "Fgfr2-dE18-Tacc2.2", "Fgfr2-dE18-Tacc2.3",
                                     "Fgfr2-dE18-IGR1.1", "Fgfr2-dE18-IGR1.2", 
                                     "Fgfr2-dE18-IGR2.1", "Fgfr2-dE18-IGR2.2",
                                     "Fgfr2-E18-C2.1", "Fgfr2-E18-C2.2", "Fgfr2-E18-C2.3", 
                                     "Fgfr2-E18-C3.1", "Fgfr2-E18-C3.2",
                                     "Fgfr2-E18-C4.1", "Fgfr2-E18-C4.2"   
)
custom.sample.order.DZ.tum.IMAC
length(custom.sample.order.DZ.tum.IMAC)

#plot
ggplot(
  melted.HM.plotdata.IMAC.tum.collapsed.4 , aes(x=variable, y=combined.plot.label, fill=value)) + 
  geom_tile() + 
  scale_fill_viridis( option = "inferno", na.value = "grey30",name="log2(intensity)") + 
  coord_equal()+
  scale_y_discrete(limits=  isoform.corrected.collapsed.order.02.12.2021 )+
  scale_x_discrete(limits=  custom.sample.order.DZ.tum.IMAC)+
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
  geom_vline(xintercept = c(2.5, 5.5, 8.5, 11.5, 13.5, 15.5, 18.5, 21.5, 23.5, 25.5, 28.5, 30.5) , size=1.0, linetype = "solid", color="white")+
  annotate("text", x = 1.5, y = nrow(HM.plotdata.IMAC.tum.collapsed.2)+1, label = "Fgfr2", fontface = "bold", angle='90', hjust = 0, size=5)+
  annotate("text", x = 4, y = nrow(HM.plotdata.IMAC.tum.collapsed.2)+1, label = "Fgfr2-dE18", fontface = "bold", angle='90', hjust = 0, size=5)+
  annotate("text", x = 7, y = nrow(HM.plotdata.IMAC.tum.collapsed.2)+1, label = "Fgfr2-Bicc1", fontface = "bold", angle='90', hjust = 0, size=5)+
  annotate("text", x = 10, y = nrow(HM.plotdata.IMAC.tum.collapsed.2)+1, label = "Fgfr2-dE18-Bicc1", fontface = "bold", angle='90', hjust = 0, size=5)+
  annotate("text", x = 12.5, y = nrow(HM.plotdata.IMAC.tum.collapsed.2)+1, label = "Fgfr2-Ate1", fontface = "bold", angle='90', hjust = 0, size=5)+
  annotate("text", x = 14.5, y = nrow(HM.plotdata.IMAC.tum.collapsed.2)+1, label = "Fgfr2-dE18-Ate1", fontface = "bold", angle='90', hjust = 0, size=5)+
  annotate("text", x = 17, y = nrow(HM.plotdata.IMAC.tum.collapsed.2)+1, label = "Fgfr2-Tacc2", fontface = "bold", angle='90', hjust = 0, size=5)+
  annotate("text", x = 20, y = nrow(HM.plotdata.IMAC.tum.collapsed.2)+1, label = "Fgfr2-dE18-Tacc2", fontface = "bold", angle='90', hjust = 0, size=5)+
  annotate("text", x = 22.5, y = nrow(HM.plotdata.IMAC.tum.collapsed.2)+1, label = "Fgfr2-dE18-IGR1", fontface = "bold", angle='90', hjust = 0, size=5)+
  annotate("text", x = 24.5, y = nrow(HM.plotdata.IMAC.tum.collapsed.2)+1, label = "Fgfr2-dE18-IGR2", fontface = "bold", angle='90', hjust = 0, size=5)+
  annotate("text", x = 27, y = nrow(HM.plotdata.IMAC.tum.collapsed.2)+1, label = "Fgfr2-E18-C2", fontface = "bold", angle='90', hjust = 0, size=5)+
  annotate("text", x = 29.5, y = nrow(HM.plotdata.IMAC.tum.collapsed.2)+1, label = "Fgfr2-E18-C3", fontface = "bold", angle='90', hjust = 0, size=5)+
  annotate("text", x = 31.5, y = nrow(HM.plotdata.IMAC.tum.collapsed.2)+1, label = "Fgfr2-E18-C4", fontface = "bold", angle='90', hjust = 0, size=5)+
  annotate("text", x = 32, y = nrow(HM.plotdata.IMAC.tum.collapsed.2)+12, label = "",  color = "transparent")

#ggsave("DZ.tumors.IMAC.Fgfr2.sites.pdf", useDingbats=FALSE,  width = 18, height =14, units = "cm") #


#############################################################################################################################################################
#############################################################################################################################################################

##############################################################################################################################################################
#tumors IMAC:  single sample PTMSEA (ssPTMSEA) using mouse PTMSigDB

#load("DZ.tum.IMAC_expanded_log2_normclass1.PS")
glimpse(DZ.tum.IMAC_expanded_log2_normclass1.PS) #20196

###reshape data for ssPTMSEA
DZ.tumors.IMAC.MBR.ON_expanded_log2_norm.class1.PS_ssPTMSEA <- DZ.tum.IMAC_expanded_log2_normclass1.PS
glimpse(DZ.tumors.IMAC.MBR.ON_expanded_log2_norm.class1.PS_ssPTMSEA )

## start define function to make aa15.window
make.aa.15.window.from.31.sequence.window <- function(input){
  nchar(input)
  input <- unlist(str_split(string=input, pattern="")); input
  return(paste0(input[9:23], collapse="") )
}
make.aa.15.window.from.31.sequence.window(input="RLETSTSCFYQPQRRsVILDGRSGRQIE___")
make.aa.15.window.from.31.sequence.window(input="RLETSTSCFYQPQRRsVILDGRSGRQIE___")
## stop define function to make aa15.window

### prepare  aa15.window from 31 aa Sequence.window MaxQuant
# take first entry sequence window if there are two entries
DZ.tumors.IMAC.MBR.ON_expanded_log2_norm.class1.PS_ssPTMSEA$aa15.window <- pbsapply(DZ.tumors.IMAC.MBR.ON_expanded_log2_norm.class1.PS_ssPTMSEA$Sequence.window, function(x) unlist(str_split(string =x, pattern = ";"))[1]  )
DZ.tumors.IMAC.MBR.ON_expanded_log2_norm.class1.PS_ssPTMSEA$aa15.window <- pbsapply(DZ.tumors.IMAC.MBR.ON_expanded_log2_norm.class1.PS_ssPTMSEA$aa15.window, function(x) make.aa.15.window.from.31.sequence.window(input = x)  )

glimpse(DZ.tumors.IMAC.MBR.ON_expanded_log2_norm.class1.PS_ssPTMSEA)

table(duplicated(DZ.tumors.IMAC.MBR.ON_expanded_log2_norm.class1.PS_ssPTMSEA$aa15.window))#3318
which(duplicated(DZ.tumors.IMAC.MBR.ON_expanded_log2_norm.class1.PS_ssPTMSEA$aa15.window))


### prepare new data frame and remove aa15-window duplicates

DZ.tumors.IMAC.MBR.ON_expanded_log2_norm.class1.PS_ssPTMSEA_2 <- DZ.tumors.IMAC.MBR.ON_expanded_log2_norm.class1.PS_ssPTMSEA
glimpse(DZ.tumors.IMAC.MBR.ON_expanded_log2_norm.class1.PS_ssPTMSEA_2 )

# find unique aa15.windows
unique.15aa.windows <- unique(DZ.tumors.IMAC.MBR.ON_expanded_log2_norm.class1.PS_ssPTMSEA_2$aa15.window)
glimpse(unique.15aa.windows) #16878


# PTMSEA does not accept duplicate aa15.windows
# if there are duplicate entries (e.g. because of _1_2_3) here, use the entry with the highest sum intensity over all samples, 
# this will favour entries with low NAs, good signal an can possibly distinguish duplicate situations with the same number of data points present
DZ.tumors.IMAC.MBR.ON_expanded_log2_norm.class1.PS_ssPTMSEA_3 <- c()
for (i in unique.15aa.windows){
  print(i)
  
  temp <-DZ.tumors.IMAC.MBR.ON_expanded_log2_norm.class1.PS_ssPTMSEA_2 %>% filter(aa15.window %in% c(i))
  temp$temp.sum <- apply(temp[, c( "Norm.Intensity.A",  "Norm.Intensity.N",  "Norm.Intensity.F",  "Norm.Intensity.R",  "Norm.Intensity.D",  "Norm.Intensity.P",  "Norm.Intensity.AB",
                                   "Norm.Intensity.H",  "Norm.Intensity.T",  "Norm.Intensity.AD", "Norm.Intensity.C", "Norm.Intensity.O",  "Norm.Intensity.AA", "Norm.Intensity.G", 
                                   "Norm.Intensity.S",  "Norm.Intensity.E", "Norm.Intensity.Q",  "Norm.Intensity.AC", "Norm.Intensity.I",  "Norm.Intensity.U",  "Norm.Intensity.AE",
                                   "Norm.Intensity.J",  "Norm.Intensity.V",  "Norm.Intensity.K",  "Norm.Intensity.W",  "Norm.Intensity.B", "Norm.Intensity.X",  "Norm.Intensity.AF",
                                   "Norm.Intensity.L",  "Norm.Intensity.Y",  "Norm.Intensity.M",  "Norm.Intensity.Z" 
  )], 1, function(x) sum(x, na.rm=T))

  temp <- temp %>% arrange(desc(temp.sum)) #
  temp <- temp[1,]
  
  DZ.tumors.IMAC.MBR.ON_expanded_log2_norm.class1.PS_ssPTMSEA_3 <- bind_rows(DZ.tumors.IMAC.MBR.ON_expanded_log2_norm.class1.PS_ssPTMSEA_3 , temp)
}

glimpse(DZ.tumors.IMAC.MBR.ON_expanded_log2_norm.class1.PS_ssPTMSEA_3)
table(duplicated(DZ.tumors.IMAC.MBR.ON_expanded_log2_norm.class1.PS_ssPTMSEA_3$aa15.window))


# add_p to end of 15aa window
DZ.tumors.IMAC.MBR.ON_expanded_log2_norm.class1.PS_ssPTMSEA_3$aa15.window.2 <- paste0(DZ.tumors.IMAC.MBR.ON_expanded_log2_norm.class1.PS_ssPTMSEA_3$aa15.window, "-p")
glimpse(DZ.tumors.IMAC.MBR.ON_expanded_log2_norm.class1.PS_ssPTMSEA_3) 


### prepare input for PTMSEA via https://cloud.genepattern.org
gene.pattern.PTM.SEA.input <- DZ.tumors.IMAC.MBR.ON_expanded_log2_norm.class1.PS_ssPTMSEA_3 %>% dplyr::select(aa15.window.2, 
                                                                                                              "Norm.Intensity.A",  "Norm.Intensity.N",  "Norm.Intensity.F",  "Norm.Intensity.R",  "Norm.Intensity.D",  "Norm.Intensity.P",  "Norm.Intensity.AB",
                                                                                                              "Norm.Intensity.H",  "Norm.Intensity.T",  "Norm.Intensity.AD", "Norm.Intensity.C", "Norm.Intensity.O",  "Norm.Intensity.AA", "Norm.Intensity.G", 
                                                                                                              "Norm.Intensity.S",  "Norm.Intensity.E", "Norm.Intensity.Q",  "Norm.Intensity.AC", "Norm.Intensity.I",  "Norm.Intensity.U",  "Norm.Intensity.AE",
                                                                                                              "Norm.Intensity.J",  "Norm.Intensity.V",  "Norm.Intensity.K",  "Norm.Intensity.W",  "Norm.Intensity.B", "Norm.Intensity.X",  "Norm.Intensity.AF",
                                                                                                              "Norm.Intensity.L",  "Norm.Intensity.Y",  "Norm.Intensity.M",  "Norm.Intensity.Z" 
)
glimpse(gene.pattern.PTM.SEA.input)
nrow(gene.pattern.PTM.SEA.input) #16878

#impute a zero (e.g. for black and white regulation situations) 
gene.pattern.PTM.SEA.input[is.na(gene.pattern.PTM.SEA.input)] <- 0
glimpse(gene.pattern.PTM.SEA.input)


#prepare gct input file with cmapR package
gene.pattern.PTM.SEA.input_GCT <- as.matrix(gene.pattern.PTM.SEA.input[,c( "Norm.Intensity.A",  "Norm.Intensity.N",  "Norm.Intensity.F",  "Norm.Intensity.R",  "Norm.Intensity.D",  "Norm.Intensity.P",  "Norm.Intensity.AB",
                                                                           "Norm.Intensity.H",  "Norm.Intensity.T",  "Norm.Intensity.AD", "Norm.Intensity.C", "Norm.Intensity.O",  "Norm.Intensity.AA", "Norm.Intensity.G", 
                                                                           "Norm.Intensity.S",  "Norm.Intensity.E", "Norm.Intensity.Q",  "Norm.Intensity.AC", "Norm.Intensity.I",  "Norm.Intensity.U",  "Norm.Intensity.AE",
                                                                           "Norm.Intensity.J",  "Norm.Intensity.V",  "Norm.Intensity.K",  "Norm.Intensity.W",  "Norm.Intensity.B", "Norm.Intensity.X",  "Norm.Intensity.AF",
                                                                           "Norm.Intensity.L",  "Norm.Intensity.Y",  "Norm.Intensity.M",  "Norm.Intensity.Z" )])
glimpse(gene.pattern.PTM.SEA.input_GCT)
colnames(gene.pattern.PTM.SEA.input_GCT) <- c( "Norm.Intensity.A",  "Norm.Intensity.N",  "Norm.Intensity.F",  "Norm.Intensity.R",  "Norm.Intensity.D",  "Norm.Intensity.P",  "Norm.Intensity.AB",
                                               "Norm.Intensity.H",  "Norm.Intensity.T",  "Norm.Intensity.AD", "Norm.Intensity.C", "Norm.Intensity.O",  "Norm.Intensity.AA", "Norm.Intensity.G", 
                                               "Norm.Intensity.S",  "Norm.Intensity.E", "Norm.Intensity.Q",  "Norm.Intensity.AC", "Norm.Intensity.I",  "Norm.Intensity.U",  "Norm.Intensity.AE",
                                               "Norm.Intensity.J",  "Norm.Intensity.V",  "Norm.Intensity.K",  "Norm.Intensity.W",  "Norm.Intensity.B", "Norm.Intensity.X",  "Norm.Intensity.AF",
                                               "Norm.Intensity.L",  "Norm.Intensity.Y",  "Norm.Intensity.M",  "Norm.Intensity.Z" )
rownames(gene.pattern.PTM.SEA.input_GCT) <- gene.pattern.PTM.SEA.input$aa15.window.2
glimpse(gene.pattern.PTM.SEA.input_GCT)

gene.pattern.PTM.SEA.input_GCT.2 <- new("GCT", mat=gene.pattern.PTM.SEA.input_GCT)
glimpse(gene.pattern.PTM.SEA.input_GCT.2)

#save GCT file for PTMSEA input
write_gct(gene.pattern.PTM.SEA.input_GCT.2, "tumors.IMAC.mouse.windows.zero.imputed.ssPTMSEA.input.gct")

### run PTMSEA via Gene Pattern at https://cloud.genepattern.org ################################################################
# pathwaydb: ptm.sig.db.all.flanking.mouse.v1.9.0.gmt see https://github.com/broadinstitute/ssGSEA2.0
# all settings default
# save resulsta and ...combined.gct file
# change ...combined.gct file type to .txt and remove first two rows


###############################################################################################################################################################
###############################################################################################################################################################
### combine ssPTMSEA results mouse tumors IMAC and KB1P mouse tumors


#load ssPTMSEA results FGFR2 tumors
ss.PTMSEA.result.DZ.IMAC.tum <- fread("FGFR2.tumors.IMAC.ssPTMSEA.result.combined.txt")
glimpse(ss.PTMSEA.result.DZ.IMAC.tum )


#load ssPTMSEA results KB1P tumors
ss.PTMSEA.result.EG.Tiox.tum.NT <- fread("KB1P.tumors.ssPTMSEA.result.combined.txt")
glimpse(ss.PTMSEA.result.EG.Tiox.tum.NT)

#combine sPTMSEA results
ss.PTMSE.global.phospho.DZ.EG.comb <- merge(ss.PTMSEA.result.DZ.IMAC.tum %>% select(id, Norm.Intensity.A : Norm.Intensity.Z), ss.PTMSEA.result.EG.Tiox.tum.NT %>% select(id, Int.1_2 : Int.51_52), by.x = "id", by.y = "id", all.x = T, all.y = T, sort = F)
glimpse(ss.PTMSE.global.phospho.DZ.EG.comb)

### count data presence and calculate some measures
#samples DZ (FGFR2) 32; samples EG (KB1P): 24

ss.PTMSE.global.phospho.DZ.EG.comb$data.presence.DZ <- apply(ss.PTMSE.global.phospho.DZ.EG.comb[, c("Norm.Intensity.A",  "Norm.Intensity.N",  "Norm.Intensity.F",  "Norm.Intensity.R",  "Norm.Intensity.D",  "Norm.Intensity.P",  "Norm.Intensity.AB",
                                                                                                    "Norm.Intensity.H",  "Norm.Intensity.T",  "Norm.Intensity.AD", "Norm.Intensity.C", "Norm.Intensity.O",  "Norm.Intensity.AA", "Norm.Intensity.G",
                                                                                                    "Norm.Intensity.S",  "Norm.Intensity.E", "Norm.Intensity.Q",  "Norm.Intensity.AC", "Norm.Intensity.I",  "Norm.Intensity.U",  "Norm.Intensity.AE",
                                                                                                    "Norm.Intensity.J",  "Norm.Intensity.V",  "Norm.Intensity.K",  "Norm.Intensity.W",  "Norm.Intensity.B", "Norm.Intensity.X",  "Norm.Intensity.AF",
                                                                                                    "Norm.Intensity.L",  "Norm.Intensity.Y",  "Norm.Intensity.M",  "Norm.Intensity.Z")],1, function(x) sum(!is.na(x)))

ss.PTMSE.global.phospho.DZ.EG.comb$data.presence.EG <- apply(ss.PTMSE.global.phospho.DZ.EG.comb[, c( "Int.1_2", "Int.3_4", "Int.9_10","Int.11_12","Int.17_18","Int.19_20","Int.25_26","Int.27_28", "Int.33_34","Int.35_36", "Int.41_42","Int.43_44",
                                                                                                     "Int.65_66","Int.67_68","Int.73_74","Int.75_76", "Int.81_82","Int.83_84","Int.89_90","Int.91_92", "Int.57_58","Int.59_60",
                                                                                                     "Int.49_50","Int.51_52")],1, function(x) sum(!is.na(x)))

ss.PTMSE.global.phospho.DZ.EG.comb$median.NES.DZ<- apply(ss.PTMSE.global.phospho.DZ.EG.comb[, c("Norm.Intensity.A",  "Norm.Intensity.N",  "Norm.Intensity.F",  "Norm.Intensity.R",  "Norm.Intensity.D",  "Norm.Intensity.P",  "Norm.Intensity.AB",
                                                                                                "Norm.Intensity.H",  "Norm.Intensity.T",  "Norm.Intensity.AD", "Norm.Intensity.C", "Norm.Intensity.O",  "Norm.Intensity.AA", "Norm.Intensity.G",
                                                                                                "Norm.Intensity.S",  "Norm.Intensity.E", "Norm.Intensity.Q",  "Norm.Intensity.AC", "Norm.Intensity.I",  "Norm.Intensity.U",  "Norm.Intensity.AE",
                                                                                                "Norm.Intensity.J",  "Norm.Intensity.V",  "Norm.Intensity.K",  "Norm.Intensity.W",  "Norm.Intensity.B", "Norm.Intensity.X",  "Norm.Intensity.AF",
                                                                                                "Norm.Intensity.L",  "Norm.Intensity.Y",  "Norm.Intensity.M",  "Norm.Intensity.Z")],1, function(x) mean(x, na.rm=T))


ss.PTMSE.global.phospho.DZ.EG.comb$median.NES.EG<- apply(ss.PTMSE.global.phospho.DZ.EG.comb[, c("Int.1_2", "Int.3_4", "Int.9_10","Int.11_12","Int.17_18","Int.19_20","Int.25_26","Int.27_28", "Int.33_34","Int.35_36", "Int.41_42","Int.43_44",
                                                                                                "Int.65_66","Int.67_68","Int.73_74","Int.75_76", "Int.81_82","Int.83_84","Int.89_90","Int.91_92", "Int.57_58","Int.59_60",
                                                                                                "Int.49_50","Int.51_52")],1, function(x) mean(x, na.rm=T))


ss.PTMSE.global.phospho.DZ.EG.comb <- ss.PTMSE.global.phospho.DZ.EG.comb %>% arrange(median.NES.DZ)
glimpse(ss.PTMSE.global.phospho.DZ.EG.comb)

### determine shared and unique pathway terms for FGFR2 tumors and KB1P tumors
ss.PTMSE.global.phospho.DZ.EG.comb_DZ.unique <- ss.PTMSE.global.phospho.DZ.EG.comb %>% filter(data.presence.EG == 0) %>% arrange(median.NES.DZ)
ss.PTMSE.global.phospho.DZ.EG.comb_DZ.unique

ss.PTMSE.global.phospho.DZ.EG.comb_EG.unique <- ss.PTMSE.global.phospho.DZ.EG.comb %>% filter(data.presence.DZ == 0) %>% arrange(median.NES.DZ)
ss.PTMSE.global.phospho.DZ.EG.comb_EG.unique

ss.PTMSE.global.phospho.DZ.EG.comb_EG.shared <- ss.PTMSE.global.phospho.DZ.EG.comb %>% filter(data.presence.DZ == 32, data.presence.EG == 24) %>% arrange(median.NES.DZ)
ss.PTMSE.global.phospho.DZ.EG.comb_EG.shared

#combine unique and shared data frames
ss.PTMSE.global.phospho.DZ.EG.comb_2 <- bind_rows(ss.PTMSE.global.phospho.DZ.EG.comb_DZ.unique, ss.PTMSE.global.phospho.DZ.EG.comb_EG.shared, ss.PTMSE.global.phospho.DZ.EG.comb_EG.unique)
glimpse(ss.PTMSE.global.phospho.DZ.EG.comb_2)


#reshape data for plot
melted.ss.PTMSE.global.phospho.DZ.EG.comb <- reshape2::melt(ss.PTMSE.global.phospho.DZ.EG.comb_2 %>% select(-data.presence.DZ, -data.presence.EG, -median.NES.DZ, -median.NES.EG))
melted.ss.PTMSE.global.phospho.DZ.EG.comb$group <- pbsapply(melted.ss.PTMSE.global.phospho.DZ.EG.comb$variable, function(x) if(x %in% c(
  "Norm.Intensity.A",  "Norm.Intensity.N",  "Norm.Intensity.F",  "Norm.Intensity.R",  "Norm.Intensity.D",  "Norm.Intensity.P",  "Norm.Intensity.AB",
                                                                                                                                        "Norm.Intensity.H",  "Norm.Intensity.T",  "Norm.Intensity.AD", "Norm.Intensity.C", "Norm.Intensity.O",  "Norm.Intensity.AA", "Norm.Intensity.G", 
                                                                                                                                        "Norm.Intensity.S",  "Norm.Intensity.E", "Norm.Intensity.Q",  "Norm.Intensity.AC", "Norm.Intensity.I",  "Norm.Intensity.U",  "Norm.Intensity.AE",
                                                                                                                                        "Norm.Intensity.J",  "Norm.Intensity.V",  "Norm.Intensity.K",  "Norm.Intensity.W",  "Norm.Intensity.B", "Norm.Intensity.X",  "Norm.Intensity.AF",
                                                                                                                                        "Norm.Intensity.L",  "Norm.Intensity.Y",  "Norm.Intensity.M",  "Norm.Intensity.Z")) {"Fgfr2 variant tumors"} else if (x %in% c("Int.1_2", "Int.3_4", "Int.9_10","Int.11_12","Int.17_18","Int.19_20","Int.25_26","Int.27_28", "Int.33_34","Int.35_36", "Int.41_42","Int.43_44",
                                                                                                                                                                                                                                                                       "Int.65_66","Int.67_68","Int.73_74","Int.75_76", "Int.81_82","Int.83_84","Int.89_90","Int.91_92", "Int.57_58","Int.59_60",
                                                                                                                                                                                                                                                                       "Int.49_50","Int.51_52")) {"KB1P tumors"} else {"xxx"})

glimpse(melted.ss.PTMSE.global.phospho.DZ.EG.comb)
table(melted.ss.PTMSE.global.phospho.DZ.EG.comb$group)


# define order in plot
order.plot.ssptmses.tums.180122 <- ss.PTMSE.global.phospho.DZ.EG.comb_2$id
order.plot.ssptmses.tums.180122

#plot results
ggplot()+
  theme_classic()+
  theme(panel.grid.major = element_line(),
        panel.grid.minor = element_line()) +
  geom_boxplot(data = melted.ss.PTMSE.global.phospho.DZ.EG.comb, aes(x = value, y = fct_inorder(id), fill= group) , outlier.size = 0.25,  outlier.stroke = 0.25)+
  scale_fill_manual(values=c("red", "grey75"))+
  theme(legend.position="top", legend.justification = "center")+
  theme(axis.text.x=element_text(angle=0,vjust=0.,hjust=0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+ 
  xlab("NES")+
  ylab(NULL)+
  scale_y_discrete(limits=  rev(order.plot.ssptmses.tums.180122))+
  theme(axis.text.y= element_text(size=10))+  
  theme(axis.text.x= element_text(size=12)) 

#ggsave("ss.PTMSEA.tumors.FGFR2.KB1P.combined.pdf", useDingbats=FALSE,  width = 16, height =18, units = "cm") #


###################################################################################################################################
###################################################################################################################################
##################################################################################################################################
##################################################################################################################################

### tumors IMAC: select and heatmap visualize candidate phosphosites and their upstream kinases 
### perform several different two group analyses and subject results to RoKAI tool (https://rokai.io/)
### combine results and visualize selection

### sample and group overview
#Group 1 (G1)
#"1_Fgfr2","2_Fgfr2", "13_Fgfr2-Ate1","14_Fgfr2-Ate1",

#Group 2 (G2)
#"4_Fgfr2-dE18","5_Fgfr2-dE18","6_Fgfr2-dE18","23_Fgfr2-dE18-IGR1","24_Fgfr2-dE18-IGR1","25_Fgfr2-dE18-IGR2","26_Fgfr2-dE18-IGR2","29_Fgfr2-E18-C3","30_Fgfr2-E18-C3","31_Fgfr2-E18-C4","32_Fgfr2-E18-C4",

#Group 3 (G3)
#"7_Fgfr2-Bicc1","8_Fgfr2-Bicc1","9_Fgfr2-Bicc1","17_Fgfr2-Tacc2","18_Fgfr2-Tacc2","19_Fgfr2-Tacc2",

#Group 4 (G4)
#"15_Fgfr2-dE18-Ate1","16_Fgfr2-dE18-Ate1","10_Fgfr2-dE18-Bicc1","11_Fgfr2-dE18-Bicc1","12_Fgfr2-dE18-Bicc1","20_Fgfr2-dE18-Tacc2","21_Fgfr2-dE18-Tacc2","22_Fgfr2-dE18-Tacc2",

#Group 5 (G5)
#"3_Fgfr2-E18-C2","27_Fgfr2-E18-C2","28_Fgfr2-E18-C2"


### two group comparisons to perform (G5 is excluded here but also heatmap plotted later)
#G1+G3 vs. G2
#G1+G3 vs. G4
#G2 vs G4

####start define FC function FR_FC_LOGData_2.3.17 ######
FR_FC_LOGData_2.3.17 <- function(control, treatment, LOGtype) {
  
  control <- control[which(!is.infinite(control))]; control
  control <- control[which(!is.na(control))]; control
  control <- control[which(control != 0)]; control
  control
  
  treatment <- treatment[which(!is.infinite(treatment))]; treatment
  treatment <- treatment[which(!is.na(treatment))]; treatment
  treatment <- treatment[which(treatment != 0)]; treatment
  treatment
  
  length.control <- length(control);length.control
  length.treatment <- length(treatment); length.treatment
  
  mean.control <- mean(control); mean.control
  mean.treatment <- mean(treatment); mean.treatment
  
  if (length.control == 0 & length.treatment == 0) {
    FC <- 1
  } else if ( length.control >= 1 & length.treatment == 0)  {
    FC <- -10000
  } else if (length.control == 0 & length.treatment >= 1) {
    FC <- 10000
  } else if  (mean.control == mean.treatment) {
    FC <- 1
  } else {
    treatment.up.is.plus <- mean.control < mean.treatment; treatment.up.is.plus
    treatment.down.is.minus <- mean.control > mean.treatment; treatment.down.is.minus
    
    tempdf <- data.frame(group = c("control", "treatment"), mean= c(mean.control, mean.treatment)); tempdf
    tempdf <-  tempdf[order(tempdf$mean, decreasing = TRUE),];  tempdf
    
    FC <- tempdf$mean[1] - tempdf$mean[2]; FC
    FC <- LOGtype^FC
    FC
    if(treatment.up.is.plus) {FC <- FC} else {FC <- FC*(-1)}}
  
  return(FC)
}
#### stop define FC function FR_FC_LOGData_2.3.17 ######


## start define function to make aa15.window
make.aa.15.window.from.31.sequence.window <- function(input){
  nchar(input)
  input <- unlist(str_split(string=input, pattern="")); input
  return(paste0(input[9:23], collapse="") )
}
make.aa.15.window.from.31.sequence.window(input="RLETSTSCFYQPQRRsVILDGRSGRQIE___")
make.aa.15.window.from.31.sequence.window(input="RLETSTSCFYQPQRRsVILDGRSGRQIE___")
## stop define function to make aa15.window



#### get data for two group comparisons #########################################################################################
glimpse(DZ.tum.IMAC_expanded_log2_normclass1.PS) #20,196

DZ.tum.IMAC_expanded_log2_normclass1.PS_for.L.sep.2.group.comp <- DZ.tum.IMAC_expanded_log2_normclass1.PS %>% select(id : Mod..peptide.IDs, PS_Multiplicity, "Norm.Intensity.A",  "Norm.Intensity.N",  "Norm.Intensity.F",  "Norm.Intensity.R",  "Norm.Intensity.D",  "Norm.Intensity.P",  "Norm.Intensity.AB",
                                                                                                                     "Norm.Intensity.H",  "Norm.Intensity.T",  "Norm.Intensity.AD", "Norm.Intensity.C", "Norm.Intensity.O",  "Norm.Intensity.AA", "Norm.Intensity.G", 
                                                                                                                     "Norm.Intensity.S",  "Norm.Intensity.E", "Norm.Intensity.Q",  "Norm.Intensity.AC", "Norm.Intensity.I",  "Norm.Intensity.U",  "Norm.Intensity.AE",
                                                                                                                     "Norm.Intensity.J",  "Norm.Intensity.V",  "Norm.Intensity.K",  "Norm.Intensity.W",  "Norm.Intensity.B", "Norm.Intensity.X",  "Norm.Intensity.AF",
                                                                                                                     "Norm.Intensity.L",  "Norm.Intensity.Y",  "Norm.Intensity.M",  "Norm.Intensity.Z",  Prot.Name_AApos :  ID.FR.all.C1.PS )
glimpse(DZ.tum.IMAC_expanded_log2_normclass1.PS_for.L.sep.2.group.comp)

#change column names
colnames(DZ.tum.IMAC_expanded_log2_normclass1.PS_for.L.sep.2.group.comp) <- c("id",                          "Proteins",                    "Positions.within.proteins",   "Leading.proteins",            "Protein",                    
                                                                              "Protein.names",               "Gene.names",                  "Fasta.headers",               "Localization.prob",           "Number.of.Phospho.(STY)",    
                                                                              "Amino.acid",                  "Sequence.window",             "Phospho.(STY).Probabilities", "Position.in.peptide",         "Reverse",                    
                                                                              "Potential.contaminant",       "Positions",                   "Position",                    "Peptide.IDs",                 "Mod..peptide.IDs",           
                                                                              "PS_Multiplicity",
                                                                              "1_Fgfr2",
                                                                              "2_Fgfr2",
                                                                              
                                                                              "13_Fgfr2-Ate1",
                                                                              "14_Fgfr2-Ate1",
                                                                              
                                                                              "7_Fgfr2-Bicc1",
                                                                              "8_Fgfr2-Bicc1",
                                                                              "9_Fgfr2-Bicc1",
                                                                              
                                                                              "17_Fgfr2-Tacc2",
                                                                              "18_Fgfr2-Tacc2",
                                                                              "19_Fgfr2-Tacc2",
                                                                              
                                                                              "4_Fgfr2-dE18",
                                                                              "5_Fgfr2-dE18",
                                                                              "6_Fgfr2-dE18",
                                                                              
                                                                              "15_Fgfr2-dE18-Ate1",
                                                                              "16_Fgfr2-dE18-Ate1",
                                                                              
                                                                              "10_Fgfr2-dE18-Bicc1",
                                                                              "11_Fgfr2-dE18-Bicc1",
                                                                              "12_Fgfr2-dE18-Bicc1",
                                                                              
                                                                              "20_Fgfr2-dE18-Tacc2",
                                                                              "21_Fgfr2-dE18-Tacc2",
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
                                                                              
                                                                              "31_Fgfr2-E18-C4" ,
                                                                              "32_Fgfr2-E18-C4",
                                                                              "Prot.Name_AApos",
                                                                              "Gene.Name_AApos",            
                                                                              "ID.FR.all.C1.PS")
glimpse(DZ.tum.IMAC_expanded_log2_normclass1.PS_for.L.sep.2.group.comp)

#reorder columns
DZ.tum.IMAC_expanded_log2_normclass1.PS_for.L.sep.2.group.comp <- DZ.tum.IMAC_expanded_log2_normclass1.PS_for.L.sep.2.group.comp %>% select(
  "id",                          "Proteins",                    "Positions.within.proteins",   "Leading.proteins",            "Protein",                    
                                                                                                                                            "Protein.names",               "Gene.names",                  "Fasta.headers",               "Localization.prob",           "Number.of.Phospho.(STY)",    
                                                                                                                                            "Amino.acid",                  "Sequence.window",             "Phospho.(STY).Probabilities", "Position.in.peptide",         "Reverse",                    
                                                                                                                                            "Potential.contaminant",       "Positions",                   "Position",                    "Peptide.IDs",                 "Mod..peptide.IDs",           
                                                                                                                                            "PS_Multiplicity",
                                                                                                                                            #Group 1
                                                                                                                                            "1_Fgfr2",
                                                                                                                                            "2_Fgfr2",
                                                                                                                                            
                                                                                                                                            "13_Fgfr2-Ate1",
                                                                                                                                            "14_Fgfr2-Ate1",
                                                                                                                                            
                                                                                                                                            
                                                                                                                                            #Group2
                                                                                                                                            "4_Fgfr2-dE18",
                                                                                                                                            "5_Fgfr2-dE18",
                                                                                                                                            "6_Fgfr2-dE18",
                                                                                                                                            
                                                                                                                                            "23_Fgfr2-dE18-IGR1",
                                                                                                                                            "24_Fgfr2-dE18-IGR1",
                                                                                                                                            
                                                                                                                                            "25_Fgfr2-dE18-IGR2",
                                                                                                                                            "26_Fgfr2-dE18-IGR2",
                                                                                                                                            
                                                                                                                                            "29_Fgfr2-E18-C3",
                                                                                                                                            "30_Fgfr2-E18-C3",
                                                                                                                                            
                                                                                                                                            "31_Fgfr2-E18-C4" ,
                                                                                                                                            "32_Fgfr2-E18-C4",
                                                                                                                                            
                                                                                                                                            #Group3
                                                                                                                                            "7_Fgfr2-Bicc1",
                                                                                                                                            "8_Fgfr2-Bicc1",
                                                                                                                                            "9_Fgfr2-Bicc1",
                                                                                                                                            
                                                                                                                                            "17_Fgfr2-Tacc2",
                                                                                                                                            "18_Fgfr2-Tacc2",
                                                                                                                                            "19_Fgfr2-Tacc2",
                                                                                                                                            
                                                                                                                                            #Group4
                                                                                                                                            "15_Fgfr2-dE18-Ate1",
                                                                                                                                            "16_Fgfr2-dE18-Ate1",
                                                                                                                                            
                                                                                                                                            "10_Fgfr2-dE18-Bicc1",
                                                                                                                                            "11_Fgfr2-dE18-Bicc1",
                                                                                                                                            "12_Fgfr2-dE18-Bicc1",
                                                                                                                                            
                                                                                                                                            "20_Fgfr2-dE18-Tacc2",
                                                                                                                                            "21_Fgfr2-dE18-Tacc2",
                                                                                                                                            "22_Fgfr2-dE18-Tacc2",
                                                                                                                                            
                                                                                                                                            ##Group 5  (exclude)
                                                                                                                                            #"3_Fgfr2-E18-C2",
                                                                                                                                            #"27_Fgfr2-E18-C2",
                                                                                                                                            #"28_Fgfr2-E18-C2"
                                                                                                                                            
                                                                                                                                            "Prot.Name_AApos",
                                                                                                                                            "Gene.Name_AApos",            
                                                                                                                                            "ID.FR.all.C1.PS")


glimpse(DZ.tum.IMAC_expanded_log2_normclass1.PS_for.L.sep.2.group.comp)

#count data presence per group
DZ.tum.IMAC_expanded_log2_normclass1.PS_for.L.sep.2.group.comp$data.presence.Group.1 <- apply(DZ.tum.IMAC_expanded_log2_normclass1.PS_for.L.sep.2.group.comp[, c("1_Fgfr2","2_Fgfr2","13_Fgfr2-Ate1","14_Fgfr2-Ate1")],1, function(x) sum(!is.na(x)))
DZ.tum.IMAC_expanded_log2_normclass1.PS_for.L.sep.2.group.comp$data.presence.Group.2 <- apply(DZ.tum.IMAC_expanded_log2_normclass1.PS_for.L.sep.2.group.comp[, c("4_Fgfr2-dE18","5_Fgfr2-dE18","6_Fgfr2-dE18","23_Fgfr2-dE18-IGR1","24_Fgfr2-dE18-IGR1","25_Fgfr2-dE18-IGR2","26_Fgfr2-dE18-IGR2","29_Fgfr2-E18-C3","30_Fgfr2-E18-C3","31_Fgfr2-E18-C4","32_Fgfr2-E18-C4")],1, function(x) sum(!is.na(x)))
DZ.tum.IMAC_expanded_log2_normclass1.PS_for.L.sep.2.group.comp$data.presence.Group.3 <- apply(DZ.tum.IMAC_expanded_log2_normclass1.PS_for.L.sep.2.group.comp[, c("7_Fgfr2-Bicc1","8_Fgfr2-Bicc1","9_Fgfr2-Bicc1","17_Fgfr2-Tacc2","18_Fgfr2-Tacc2","19_Fgfr2-Tacc2")],1, function(x) sum(!is.na(x)))
DZ.tum.IMAC_expanded_log2_normclass1.PS_for.L.sep.2.group.comp$data.presence.Group.4 <- apply(DZ.tum.IMAC_expanded_log2_normclass1.PS_for.L.sep.2.group.comp[, c("15_Fgfr2-dE18-Ate1","16_Fgfr2-dE18-Ate1","10_Fgfr2-dE18-Bicc1","11_Fgfr2-dE18-Bicc1","12_Fgfr2-dE18-Bicc1","20_Fgfr2-dE18-Tacc2","21_Fgfr2-dE18-Tacc2","22_Fgfr2-dE18-Tacc2")],1, function(x) sum(!is.na(x)))
#DZ.tum.IMAC_expanded_log2_normclass1.PS_for.L.sep.2.group.comp$data.presence.Group.5 <- apply(DZ.tum.IMAC_expanded_log2_normclass1.PS_for.L.sep.2.group.comp[, c("3_Fgfr2-E18-C2","27_Fgfr2-E18-C2","28_Fgfr2-E18-C2")],1, function(x) sum(!is.na(x)))

glimpse(DZ.tum.IMAC_expanded_log2_normclass1.PS_for.L.sep.2.group.comp)


##################################################################################################################################
##################################################################################################################################
### G1G3_G2 DZ tumors IMAC Limma seperate 2 group comparison for RokaiTwo ########################################################

#Group 1
#"1_Fgfr2","2_Fgfr2", "13_Fgfr2-Ate1","14_Fgfr2-Ate1",

#Group 2
#"4_Fgfr2-dE18","5_Fgfr2-dE18","6_Fgfr2-dE18","23_Fgfr2-dE18-IGR1","24_Fgfr2-dE18-IGR1","25_Fgfr2-dE18-IGR2","26_Fgfr2-dE18-IGR2","29_Fgfr2-E18-C3","30_Fgfr2-E18-C3","31_Fgfr2-E18-C4","32_Fgfr2-E18-C4",

#Group 3
#"7_Fgfr2-Bicc1","8_Fgfr2-Bicc1","9_Fgfr2-Bicc1","17_Fgfr2-Tacc2","18_Fgfr2-Tacc2","19_Fgfr2-Tacc2",

#Group 4
#"15_Fgfr2-dE18-Ate1","16_Fgfr2-dE18-Ate1","10_Fgfr2-dE18-Bicc1","11_Fgfr2-dE18-Bicc1","12_Fgfr2-dE18-Bicc1","20_Fgfr2-dE18-Tacc2","21_Fgfr2-dE18-Tacc2","22_Fgfr2-dE18-Tacc2",

#Group 5
#"3_Fgfr2-E18-C2","27_Fgfr2-E18-C2","28_Fgfr2-E18-C2"



## select columns of interest
glimpse(DZ.tum.IMAC_expanded_log2_normclass1.PS_for.L.sep.2.group.comp)

G1and3_G2_Limma.2gc <- DZ.tum.IMAC_expanded_log2_normclass1.PS_for.L.sep.2.group.comp %>% select("ID.FR.all.C1.PS",
                                                                                                 "id",                          "Proteins",                    "Positions.within.proteins",   "Leading.proteins",            "Protein",                    
                                                                                                 "Protein.names",               "Gene.names",                  "Fasta.headers",               "Localization.prob",           "Number.of.Phospho.(STY)",    
                                                                                                 "Amino.acid",                  "Sequence.window",             "Phospho.(STY).Probabilities", "Position.in.peptide",         "Reverse",                    
                                                                                                 "Potential.contaminant",       "Positions",                   "Position",                    "Peptide.IDs",                 "Mod..peptide.IDs",   
                                                                                                 "PS_Multiplicity", "Prot.Name_AApos", "Gene.Name_AApos",
                                                                                                 
                                                                                                 #Group 1
                                                                                                 "1_Fgfr2","2_Fgfr2", "13_Fgfr2-Ate1","14_Fgfr2-Ate1",
                                                                                                 #Group 3
                                                                                                 "7_Fgfr2-Bicc1","8_Fgfr2-Bicc1","9_Fgfr2-Bicc1","17_Fgfr2-Tacc2","18_Fgfr2-Tacc2","19_Fgfr2-Tacc2",
                                                                                                 
                                                                                                 #Group 2
                                                                                                 "4_Fgfr2-dE18","5_Fgfr2-dE18","6_Fgfr2-dE18","23_Fgfr2-dE18-IGR1","24_Fgfr2-dE18-IGR1","25_Fgfr2-dE18-IGR2","26_Fgfr2-dE18-IGR2","29_Fgfr2-E18-C3","30_Fgfr2-E18-C3","31_Fgfr2-E18-C4","32_Fgfr2-E18-C4",
                                                                                                 
                                                                                                 data.presence.Group.1,
                                                                                                 data.presence.Group.3,
                                                                                                 data.presence.Group.2
                                                                                                 
)

#extend data apresence calculation
G1and3_G2_Limma.2gc$data.presence.all <- G1and3_G2_Limma.2gc$data.presence.Group.1 + G1and3_G2_Limma.2gc$data.presence.Group.3 + G1and3_G2_Limma.2gc$data.presence.Group.2
G1and3_G2_Limma.2gc$data.presence.Group.1and3 <- G1and3_G2_Limma.2gc$data.presence.Group.1 + G1and3_G2_Limma.2gc$data.presence.Group.3
glimpse(G1and3_G2_Limma.2gc)

#remove complete zero rows
table(G1and3_G2_Limma.2gc$data.presence.all == 0) #F17831 T2365 

G1and3_G2_Limma.2gc <- G1and3_G2_Limma.2gc %>% filter(data.presence.all > 0)
nrow(G1and3_G2_Limma.2gc) #17831

## X out of Y filter 
#Group 1+3: 4+6=10 samples => 50% => 5 samples => 6 samples
#Group 2: 11 samples => 50% => 5.5 samples => 6 samples
glimpse(G1and3_G2_Limma.2gc)
temp.left <- filter(G1and3_G2_Limma.2gc, data.presence.Group.1and3 >= 6) ; nrow(temp.left)    #5594
temp.right <- filter(G1and3_G2_Limma.2gc, data.presence.Group.2 >= 6) ; nrow(temp.right) #6828

temp.left.and.right <- bind_rows(temp.left , temp.right) ; nrow(temp.left.and.right)
temp.left.and.right <- distinct(temp.left.and.right) ; nrow(temp.left.and.right) #7592

G1and3_G2_Limma.2gc.XYfilter <- temp.left.and.right
nrow(G1and3_G2_Limma.2gc.XYfilter) #6756
length(unique(G1and3_G2_Limma.2gc.XYfilter$Prot.Name_AApos)) #6868

#add FRID
G1and3_G2_Limma.2gc.XYfilter$FRID <- 1:nrow(G1and3_G2_Limma.2gc.XYfilter) #here FRID is also row number
glimpse(G1and3_G2_Limma.2gc.XYfilter)

#find likely Limma NA pval situations 
#Group 1and3 90% => 9 samples
#Group 2 90% => 10 samples
possible.NApvals.A <- G1and3_G2_Limma.2gc.XYfilter %>% filter((data.presence.Group.1and3  >= 9 & data.presence.Group.2  == 0) | (data.presence.Group.1and3  == 0 & data.presence.Group.2  >= 10)  )

glimpse(possible.NApvals.A)
length(unique(possible.NApvals.A$FRID))
FRIDS.possible.NApvals <- possible.NApvals.A$FRID
FRIDS.possible.NApvals 

# for these filtered  black and white situation add zero so that p.val can be calculated
W <- G1and3_G2_Limma.2gc.XYfilter[FRIDS.possible.NApvals, ]
glimpse(W)
W <- W %>% mutate_at(vars(
  #Group 1
  "1_Fgfr2","2_Fgfr2", "13_Fgfr2-Ate1","14_Fgfr2-Ate1",
  #Group 3
  "7_Fgfr2-Bicc1","8_Fgfr2-Bicc1","9_Fgfr2-Bicc1","17_Fgfr2-Tacc2","18_Fgfr2-Tacc2","19_Fgfr2-Tacc2",
  
  #Group 2
  "4_Fgfr2-dE18","5_Fgfr2-dE18","6_Fgfr2-dE18","23_Fgfr2-dE18-IGR1","24_Fgfr2-dE18-IGR1","25_Fgfr2-dE18-IGR2","26_Fgfr2-dE18-IGR2","29_Fgfr2-E18-C3","30_Fgfr2-E18-C3","31_Fgfr2-E18-C4","32_Fgfr2-E18-C4"),
  ~replace(., is.na(.), 0))
glimpse(W)

#combine again
G1and3_G2_Limma.2gc.XYfilter.2 <- G1and3_G2_Limma.2gc.XYfilter %>% filter(! FRID %in% FRIDS.possible.NApvals ) 
G1and3_G2_Limma.2gc.XYfilter.2 <- bind_rows(G1and3_G2_Limma.2gc.XYfilter.2, W)
G1and3_G2_Limma.2gc.XYfilter.2 <- G1and3_G2_Limma.2gc.XYfilter.2 %>% arrange(FRID)
glimpse(G1and3_G2_Limma.2gc.XYfilter.2)
table(duplicated(G1and3_G2_Limma.2gc.XYfilter.2$FRID))

##Limma statistics
#get data 
data.for.limma <- G1and3_G2_Limma.2gc.XYfilter.2 %>% select(
  #Group 1
  "1_Fgfr2","2_Fgfr2", "13_Fgfr2-Ate1","14_Fgfr2-Ate1",
  #Group 3
  "7_Fgfr2-Bicc1","8_Fgfr2-Bicc1","9_Fgfr2-Bicc1","17_Fgfr2-Tacc2","18_Fgfr2-Tacc2","19_Fgfr2-Tacc2",
  
  #Group 2
  "4_Fgfr2-dE18","5_Fgfr2-dE18","6_Fgfr2-dE18","23_Fgfr2-dE18-IGR1","24_Fgfr2-dE18-IGR1","25_Fgfr2-dE18-IGR2","26_Fgfr2-dE18-IGR2","29_Fgfr2-E18-C3","30_Fgfr2-E18-C3","31_Fgfr2-E18-C4","32_Fgfr2-E18-C4"
)
rownames(data.for.limma) <- paste0(G1and3_G2_Limma.2gc.XYfilter.2$ID.FR.all.C1.PS, "__", G1and3_G2_Limma.2gc.XYfilter.2$Gene.Name_AApos) 
nrow(data.for.limma)
glimpse(data.for.limma)


design <- cbind(left=c(rep(1, 10), rep(0, 11)), right=c(rep(0,10), rep(1, 11))); design
fit <- lmFit(data.for.limma, design) # Warning: Partial NA coefficients for 210 probe(s) 
cont.matrix <- makeContrasts( left.vs.right =  right - left, levels=design) 
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
topTable(fit2, adjust="BH")

Limma.p.value <- as.numeric(fit2$p.value)
glimpse(Limma.p.value)
Limma.adj.p.valueBH <- p.adjust(Limma.p.value, "BH") 
glimpse(Limma.adj.p.valueBH)
Limma.t <- as.numeric(fit2$t)
glimpse(Limma.t)

G1and3_G2_Limma.2gc.XYfilter.2.B <- G1and3_G2_Limma.2gc.XYfilter.2
G1and3_G2_Limma.2gc.XYfilter.2.B$Limma.p.value <- as.numeric(Limma.p.value)
G1and3_G2_Limma.2gc.XYfilter.2.B$Limma.adj.p.valueBH <- as.numeric(Limma.adj.p.valueBH)
G1and3_G2_Limma.2gc.XYfilter.2.B$Limma.t <- as.numeric(Limma.t)

glimpse(G1and3_G2_Limma.2gc.XYfilter.2.B)

###calculate FC
for(i in 1:nrow(G1and3_G2_Limma.2gc.XYfilter.2.B)){
  print(i)
  leftvals <- as.numeric(G1and3_G2_Limma.2gc.XYfilter.2.B[i,c("1_Fgfr2","2_Fgfr2", "13_Fgfr2-Ate1","14_Fgfr2-Ate1",     "7_Fgfr2-Bicc1","8_Fgfr2-Bicc1","9_Fgfr2-Bicc1","17_Fgfr2-Tacc2","18_Fgfr2-Tacc2","19_Fgfr2-Tacc2")]); leftvals #
  rightvals <- as.numeric(G1and3_G2_Limma.2gc.XYfilter.2.B[i,c("4_Fgfr2-dE18","5_Fgfr2-dE18","6_Fgfr2-dE18","23_Fgfr2-dE18-IGR1","24_Fgfr2-dE18-IGR1","25_Fgfr2-dE18-IGR2","26_Fgfr2-dE18-IGR2","29_Fgfr2-E18-C3","30_Fgfr2-E18-C3","31_Fgfr2-E18-C4","32_Fgfr2-E18-C4")]); rightvals #
  
  FC <- FR_FC_LOGData_2.3.17(control = leftvals, treatment = rightvals, LOGtype = 2); FC
  
  G1and3_G2_Limma.2gc.XYfilter.2.B$FC_FR[i] <- FC
}
glimpse(G1and3_G2_Limma.2gc.XYfilter.2.B)

#check numbers unique Sequence.window
glimpse(G1and3_G2_Limma.2gc.XYfilter.2.B)#7592
G1and3_G2_Limma.2gc.XYfilter.2.B %>% filter(Limma.p.value < 0.05) %>% select(Sequence.window) %>% distinct() %>% nrow() #1870 total sign
G1and3_G2_Limma.2gc.XYfilter.2.B %>% filter(Limma.p.value < 0.05) %>% filter(FC_FR > 1) %>% select(Sequence.window) %>% distinct() %>% nrow() #1089 up
G1and3_G2_Limma.2gc.XYfilter.2.B %>% filter(Limma.p.value < 0.05) %>% filter(FC_FR < 1) %>% select(Sequence.window) %>% distinct() %>% nrow() #781 down

# filter na pvalues
table(is.na(G1and3_G2_Limma.2gc.XYfilter.2.B$Limma.p.value))
G1and3_G2_Limma.2gc.XYfilter.2.B <- G1and3_G2_Limma.2gc.XYfilter.2.B %>% filter(!is.na(Limma.p.value))
glimpse(G1and3_G2_Limma.2gc.XYfilter.2.B) #7382
table(is.na(G1and3_G2_Limma.2gc.XYfilter.2.B$Limma.p.value))

#check +/- 10000 FCs
table(G1and3_G2_Limma.2gc.XYfilter.2.B$FC_FR >= 10000) #4
table(G1and3_G2_Limma.2gc.XYfilter.2.B$FC_FR <= -10000) #5
summary(G1and3_G2_Limma.2gc.XYfilter.2.B$FC_FR)

#find second min or max and add +/-5 for +/- 10000 FCs
temp <- G1and3_G2_Limma.2gc.XYfilter.2.B %>% filter(!FC_FR %in% c(10000, -10000)) %>% arrange(desc(FC_FR))
head(temp$FC_FR)
tail(temp$FC_FR)
max(temp$FC_FR)
min(temp$FC_FR)

G1and3_G2_Limma.2gc.XYfilter.2.B$FC_FR <- replace(G1and3_G2_Limma.2gc.XYfilter.2.B$FC_FR, G1and3_G2_Limma.2gc.XYfilter.2.B$FC_FR==10000, max(temp$FC_FR)+5)
G1and3_G2_Limma.2gc.XYfilter.2.B$FC_FR <- replace(G1and3_G2_Limma.2gc.XYfilter.2.B$FC_FR, G1and3_G2_Limma.2gc.XYfilter.2.B$FC_FR==-10000, min(temp$FC_FR)-5)
summary(G1and3_G2_Limma.2gc.XYfilter.2.B$FC_FR)
glimpse(G1and3_G2_Limma.2gc.XYfilter.2.B) #7382



##G1G3_G2_Limma.2gc.XYfilter.2.B reshape for RoKAI

G1G3_G2_RokaiTwo <- G1and3_G2_Limma.2gc.XYfilter.2.B

# make aa15.window
# take first entry sequence window if there are two entries
G1G3_G2_RokaiTwo$aa15.window <- pbsapply(G1G3_G2_RokaiTwo$Sequence.window, function(x) unlist(str_split(string =x, pattern = ";"))[1]  )
G1G3_G2_RokaiTwo$aa15.window <- pbsapply(G1G3_G2_RokaiTwo$aa15.window, function(x) make.aa.15.window.from.31.sequence.window(input = x)  )
glimpse(G1G3_G2_RokaiTwo)

table(duplicated(G1G3_G2_RokaiTwo$aa15.window))
which(duplicated(G1G3_G2_RokaiTwo$aa15.window))

# if there are duplicate entries (e.g. because of _1_2_3) take the entry with the best p-value
unique.15aa.windows <- unique(G1G3_G2_RokaiTwo$aa15.window)
length(unique.15aa.windows) #66657
glimpse(G1G3_G2_RokaiTwo) #7382


G1G3_G2_RokaiTwo.2 <- c()
for (i in unique.15aa.windows){
  print(i)
  temp <- G1G3_G2_RokaiTwo %>% filter(aa15.window %in% c(i))
  temp <- temp %>% arrange(Limma.p.value)
  temp <- temp[1,]
  G1G3_G2_RokaiTwo.2 <- bind_rows(G1G3_G2_RokaiTwo.2, temp)
}
table(duplicated(G1G3_G2_RokaiTwo.2$aa15.window))

#calculate log2(FC)
G1G3_G2_RokaiTwo.2$log2FC <- log2(abs(G1G3_G2_RokaiTwo.2$FC_FR))
G1G3_G2_RokaiTwo.2$sign.FC_FR <- sign(G1G3_G2_RokaiTwo.2$FC_FR)
G1G3_G2_RokaiTwo.2$log2FC <- sign(G1G3_G2_RokaiTwo.2$FC_FR) * log2(abs(G1G3_G2_RokaiTwo.2$FC_FR))
glimpse(G1G3_G2_RokaiTwo.2) 

G1G3_G2_RokaiTwo.2 <- G1G3_G2_RokaiTwo.2 %>% arrange(desc(log2FC))
head(G1G3_G2_RokaiTwo.2$FC_FR)
tail(G1G3_G2_RokaiTwo.2$FC_FR)

### prepare RoKAI input and save data
glimpse(G1G3_G2_RokaiTwo.2) 
G1G3_G2_RokaiTwo.2.06.01.2022 <- G1G3_G2_RokaiTwo.2 %>% select(Protein, Position, log2FC) 
table(is.na(G1G3_G2_RokaiTwo.2.06.01.2022$Protein))
G1G3_G2_RokaiTwo.2.06.01.2022 <- G1G3_G2_RokaiTwo.2.06.01.2022 %>% filter(!is.na(Protein))
glimpse(G1G3_G2_RokaiTwo.2.06.01.2022)
table(is.na(G1G3_G2_RokaiTwo.2.06.01.2022$Protein))
table(is.na(G1G3_G2_RokaiTwo.2.06.01.2022$Position))
table(is.na(G1G3_G2_RokaiTwo.2.06.01.2022$log2FC))
colnames(G1G3_G2_RokaiTwo.2.06.01.2022) <- c("Protein", "Position", "Quantification")


#write_csv(G1G3_G2_RokaiTwo.2.06.01.2022,"G1G3_G2.RokaiTWO.input.limma.sep.2.group.comp.csv", col_names=TRUE)
#save(G1G3_G2_RokaiTwo  , file="G1G3_G2_RokaiTwo.limma.sep.2.group.comp.Rdata")


# perform analysis via https://rokai.io/
# tool settings was default: Reference Proteome: Uniprot Mouse / Fold changes: normalized / Kinase substrate dataset: PSP + Signor / RoKAI network: KS+PPI+SD+CoEv
# save RoKAI output kinase_table.csv and kinase_targets.csv and load these later




#################################################################################################################################
## G1G3_G4 DZ tumors IMAC Limma seperate 2 group comparison for RokaiTwo ##############################################################################

#Group 1
#"1_Fgfr2","2_Fgfr2", "13_Fgfr2-Ate1","14_Fgfr2-Ate1",

#Group 2
#"4_Fgfr2-dE18","5_Fgfr2-dE18","6_Fgfr2-dE18","23_Fgfr2-dE18-IGR1","24_Fgfr2-dE18-IGR1","25_Fgfr2-dE18-IGR2","26_Fgfr2-dE18-IGR2","29_Fgfr2-E18-C3","30_Fgfr2-E18-C3","31_Fgfr2-E18-C4","32_Fgfr2-E18-C4",

#Group 3
#"7_Fgfr2-Bicc1","8_Fgfr2-Bicc1","9_Fgfr2-Bicc1","17_Fgfr2-Tacc2","18_Fgfr2-Tacc2","19_Fgfr2-Tacc2",

#Group 4
#"15_Fgfr2-dE18-Ate1","16_Fgfr2-dE18-Ate1","10_Fgfr2-dE18-Bicc1","11_Fgfr2-dE18-Bicc1","12_Fgfr2-dE18-Bicc1","20_Fgfr2-dE18-Tacc2","21_Fgfr2-dE18-Tacc2","22_Fgfr2-dE18-Tacc2",


#Group 5
#"3_Fgfr2-E18-C2","27_Fgfr2-E18-C2","28_Fgfr2-E18-C2"

## select columns of interest
glimpse(DZ.tum.IMAC_expanded_log2_normclass1.PS_for.L.sep.2.group.comp)

G1andG3_G4_Limma.2gc <- DZ.tum.IMAC_expanded_log2_normclass1.PS_for.L.sep.2.group.comp %>% select("ID.FR.all.C1.PS",
                                                                                                  "id",                          "Proteins",                    "Positions.within.proteins",   "Leading.proteins",            "Protein",                    
                                                                                                  "Protein.names",               "Gene.names",                  "Fasta.headers",               "Localization.prob",           "Number.of.Phospho.(STY)",    
                                                                                                  "Amino.acid",                  "Sequence.window",             "Phospho.(STY).Probabilities", "Position.in.peptide",         "Reverse",                    
                                                                                                  "Potential.contaminant",       "Positions",                   "Position",                    "Peptide.IDs",                 "Mod..peptide.IDs",   
                                                                                                  "PS_Multiplicity", "Prot.Name_AApos", "Gene.Name_AApos",
                                                                                                  
                                                                                                  #Group 1
                                                                                                  "1_Fgfr2","2_Fgfr2", "13_Fgfr2-Ate1","14_Fgfr2-Ate1",
                                                                                                  #Group 3
                                                                                                  "7_Fgfr2-Bicc1","8_Fgfr2-Bicc1","9_Fgfr2-Bicc1","17_Fgfr2-Tacc2","18_Fgfr2-Tacc2","19_Fgfr2-Tacc2",
                                                                                                  
                                                                                                  
                                                                                                  #Group 4
                                                                                                  "15_Fgfr2-dE18-Ate1","16_Fgfr2-dE18-Ate1","10_Fgfr2-dE18-Bicc1","11_Fgfr2-dE18-Bicc1","12_Fgfr2-dE18-Bicc1","20_Fgfr2-dE18-Tacc2","21_Fgfr2-dE18-Tacc2","22_Fgfr2-dE18-Tacc2",
                                                                                                  
                                                                                                  data.presence.Group.1,
                                                                                                  data.presence.Group.3,
                                                                                                  
                                                                                                  data.presence.Group.4
                                                                                                  
)

#extend data presence calculation
G1andG3_G4_Limma.2gc$data.presence.all <- G1andG3_G4_Limma.2gc$data.presence.Group.1 + G1andG3_G4_Limma.2gc$data.presence.Group.3 + G1andG3_G4_Limma.2gc$data.presence.Group.4
G1andG3_G4_Limma.2gc$data.presence.Group.1and3 <- G1andG3_G4_Limma.2gc$data.presence.Group.1 + G1andG3_G4_Limma.2gc$data.presence.Group.3
glimpse(G1andG3_G4_Limma.2gc)

#remove complete zero rows
table(G1andG3_G4_Limma.2gc$data.presence.all == 0) #F18133 T2063

G1andG3_G4_Limma.2gc <- G1andG3_G4_Limma.2gc %>% filter(data.presence.all > 0)
nrow(G1andG3_G4_Limma.2gc) #18133

## X out of Y filter 
#Group 1+3: 4+6=10 samples => 50% => 5 samples => 6 samples
#Group 4: 8 samples => 50% => 4 samples => 5 samples
glimpse(G1andG3_G4_Limma.2gc)
temp.left <- filter(G1andG3_G4_Limma.2gc, data.presence.Group.1and3 >= 6) ; nrow(temp.left)    #5594
temp.right <- filter(G1andG3_G4_Limma.2gc, data.presence.Group.4 >= 5) ; nrow(temp.right) #6058

temp.left.and.right <- bind_rows(temp.left , temp.right) ; nrow(temp.left.and.right)
temp.left.and.right <- distinct(temp.left.and.right) ; nrow(temp.left.and.right) #6771

G1andG3_G4_Limma.2gc.XYfilter <- temp.left.and.right
nrow(G1andG3_G4_Limma.2gc.XYfilter) #6771
length(unique(G1andG3_G4_Limma.2gc.XYfilter$Prot.Name_AApos)) #6135

#add FRID
G1andG3_G4_Limma.2gc.XYfilter$FRID <- 1:nrow(G1andG3_G4_Limma.2gc.XYfilter) #here FRID is also row number
glimpse(G1andG3_G4_Limma.2gc.XYfilter)

#find likely Limma NA pval situations 
#Group 1and3 90% => 9 samples
#Group 4 90% => 7 samples
possible.NApvals.A <- G1andG3_G4_Limma.2gc.XYfilter %>% filter((data.presence.Group.1and3  >= 9 & data.presence.Group.4  == 0) | (data.presence.Group.1and3  == 0 & data.presence.Group.4  >= 7)  )

glimpse(possible.NApvals.A)
length(unique(possible.NApvals.A$FRID))
FRIDS.possible.NApvals <- possible.NApvals.A$FRID
FRIDS.possible.NApvals 

# for these filtered  black and white situation add zero so that p.val can be calculated
W <- G1andG3_G4_Limma.2gc.XYfilter[FRIDS.possible.NApvals, ]
glimpse(W)
W <- W %>% mutate_at(vars(
  #Group 1
  "1_Fgfr2","2_Fgfr2", "13_Fgfr2-Ate1","14_Fgfr2-Ate1",
  #Group 3
  "7_Fgfr2-Bicc1","8_Fgfr2-Bicc1","9_Fgfr2-Bicc1","17_Fgfr2-Tacc2","18_Fgfr2-Tacc2","19_Fgfr2-Tacc2",
  
  #Group 4
  "15_Fgfr2-dE18-Ate1","16_Fgfr2-dE18-Ate1","10_Fgfr2-dE18-Bicc1","11_Fgfr2-dE18-Bicc1","12_Fgfr2-dE18-Bicc1","20_Fgfr2-dE18-Tacc2","21_Fgfr2-dE18-Tacc2","22_Fgfr2-dE18-Tacc2"),
  ~replace(., is.na(.), 0))#
glimpse(W)

#combine again
G1andG3_G4_Limma.2gc.XYfilter.2 <- G1andG3_G4_Limma.2gc.XYfilter %>% filter(! FRID %in% FRIDS.possible.NApvals ) 
G1andG3_G4_Limma.2gc.XYfilter.2 <- bind_rows(G1andG3_G4_Limma.2gc.XYfilter.2, W)
G1andG3_G4_Limma.2gc.XYfilter.2 <- G1andG3_G4_Limma.2gc.XYfilter.2 %>% arrange(FRID)
glimpse(G1andG3_G4_Limma.2gc.XYfilter.2)
table(duplicated(G1andG3_G4_Limma.2gc.XYfilter.2$FRID))

##Limma statistics
#get data 
data.for.limma <- G1andG3_G4_Limma.2gc.XYfilter.2 %>% select(
  #Group 1
  "1_Fgfr2","2_Fgfr2", "13_Fgfr2-Ate1","14_Fgfr2-Ate1",
  #Group 3
  "7_Fgfr2-Bicc1","8_Fgfr2-Bicc1","9_Fgfr2-Bicc1","17_Fgfr2-Tacc2","18_Fgfr2-Tacc2","19_Fgfr2-Tacc2",
  
  #Group 4
  "15_Fgfr2-dE18-Ate1","16_Fgfr2-dE18-Ate1","10_Fgfr2-dE18-Bicc1","11_Fgfr2-dE18-Bicc1","12_Fgfr2-dE18-Bicc1","20_Fgfr2-dE18-Tacc2","21_Fgfr2-dE18-Tacc2","22_Fgfr2-dE18-Tacc2"
  
)
rownames(data.for.limma) <- paste0(G1andG3_G4_Limma.2gc.XYfilter.2$ID.FR.all.C1.PS, "__", G1andG3_G4_Limma.2gc.XYfilter.2$Gene.Name_AApos) 
nrow(data.for.limma)
glimpse(data.for.limma)


design <- cbind(left=c(rep(1, 10), rep(0, 8)), right=c(rep(0,10), rep(1, 8))); design
fit <- lmFit(data.for.limma, design) # Warning: Partial NA coefficients for 33 probe(s) 
cont.matrix <- makeContrasts( left.vs.right =  right - left, levels=design) 
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
topTable(fit2, adjust="BH")

Limma.p.value <- as.numeric(fit2$p.value)
glimpse(Limma.p.value)
Limma.adj.p.valueBH <- p.adjust(Limma.p.value, "BH") #
glimpse(Limma.adj.p.valueBH)
Limma.t <- as.numeric(fit2$t)
glimpse(Limma.t)

G1andG3_G4_Limma.2gc.XYfilter.2.B <- G1andG3_G4_Limma.2gc.XYfilter.2
G1andG3_G4_Limma.2gc.XYfilter.2.B$Limma.p.value <- as.numeric(Limma.p.value)
G1andG3_G4_Limma.2gc.XYfilter.2.B$Limma.adj.p.valueBH <- as.numeric(Limma.adj.p.valueBH)
G1andG3_G4_Limma.2gc.XYfilter.2.B$Limma.t <- as.numeric(Limma.t)

glimpse(G1andG3_G4_Limma.2gc.XYfilter.2.B)

###calculate FC
for(i in 1:nrow(G1andG3_G4_Limma.2gc.XYfilter.2.B)){
  print(i)
  #i=250
  leftvals <- as.numeric(G1andG3_G4_Limma.2gc.XYfilter.2.B[i,c("1_Fgfr2","2_Fgfr2", "13_Fgfr2-Ate1","14_Fgfr2-Ate1",     "7_Fgfr2-Bicc1","8_Fgfr2-Bicc1","9_Fgfr2-Bicc1","17_Fgfr2-Tacc2","18_Fgfr2-Tacc2","19_Fgfr2-Tacc2")]); leftvals #
  rightvals <- as.numeric(G1andG3_G4_Limma.2gc.XYfilter.2.B[i,c("15_Fgfr2-dE18-Ate1","16_Fgfr2-dE18-Ate1","10_Fgfr2-dE18-Bicc1","11_Fgfr2-dE18-Bicc1","12_Fgfr2-dE18-Bicc1","20_Fgfr2-dE18-Tacc2","21_Fgfr2-dE18-Tacc2","22_Fgfr2-dE18-Tacc2")]); rightvals #
  
  FC <- FR_FC_LOGData_2.3.17(control = leftvals, treatment = rightvals, LOGtype = 2); FC
  
  G1andG3_G4_Limma.2gc.XYfilter.2.B$FC_FR[i] <- FC
}
glimpse(G1andG3_G4_Limma.2gc.XYfilter.2.B)

#check numbers unique Sequence.window
glimpse(G1andG3_G4_Limma.2gc.XYfilter.2.B)#6771
G1andG3_G4_Limma.2gc.XYfilter.2.B %>% filter(Limma.p.value < 0.05) %>% select(Sequence.window) %>% distinct() %>% nrow() #446 total sign
G1andG3_G4_Limma.2gc.XYfilter.2.B %>% filter(Limma.p.value < 0.05) %>% filter(FC_FR > 1) %>% select(Sequence.window) %>% distinct() %>% nrow() #252 up
G1andG3_G4_Limma.2gc.XYfilter.2.B %>% filter(Limma.p.value < 0.05) %>% filter(FC_FR < 1) %>% select(Sequence.window) %>% distinct() %>% nrow() #194 down

# filter na pvalues
table(is.na(G1andG3_G4_Limma.2gc.XYfilter.2.B$Limma.p.value))
G1andG3_G4_Limma.2gc.XYfilter.2.B <- G1andG3_G4_Limma.2gc.XYfilter.2.B %>% filter(!is.na(Limma.p.value))
glimpse(G1andG3_G4_Limma.2gc.XYfilter.2.B) #6738
table(is.na(G1andG3_G4_Limma.2gc.XYfilter.2.B$Limma.p.value))

#check +/- 10000 FCs
table(G1andG3_G4_Limma.2gc.XYfilter.2.B$FC_FR >= 10000) #0
table(G1andG3_G4_Limma.2gc.XYfilter.2.B$FC_FR <= -10000) #0
summary(G1andG3_G4_Limma.2gc.XYfilter.2.B$FC_FR)

#find second min or max and add +/-5 for +/- 10000 FCs
temp <- G1andG3_G4_Limma.2gc.XYfilter.2.B %>% filter(!FC_FR %in% c(10000, -10000)) %>% arrange(desc(FC_FR))
head(temp$FC_FR)
tail(temp$FC_FR)
max(temp$FC_FR)
min(temp$FC_FR)

G1andG3_G4_Limma.2gc.XYfilter.2.B$FC_FR <- replace(G1andG3_G4_Limma.2gc.XYfilter.2.B$FC_FR, G1andG3_G4_Limma.2gc.XYfilter.2.B$FC_FR==10000, max(temp$FC_FR)+5)
G1andG3_G4_Limma.2gc.XYfilter.2.B$FC_FR <- replace(G1andG3_G4_Limma.2gc.XYfilter.2.B$FC_FR, G1andG3_G4_Limma.2gc.XYfilter.2.B$FC_FR==-10000, min(temp$FC_FR)-5)
summary(G1andG3_G4_Limma.2gc.XYfilter.2.B$FC_FR)
glimpse(G1andG3_G4_Limma.2gc.XYfilter.2.B) #6738


##G1G3_G4_Limma.2gc.XYfilter.2.B reshape for RoKAI

G1G3_G4_RokaiTwo <- G1andG3_G4_Limma.2gc.XYfilter.2.B

# make aa15.window
# take first entry sequence window if there are two entries
G1G3_G4_RokaiTwo$aa15.window <- pbsapply(G1G3_G4_RokaiTwo$Sequence.window, function(x) unlist(str_split(string =x, pattern = ";"))[1]  )
G1G3_G4_RokaiTwo$aa15.window <- pbsapply(G1G3_G4_RokaiTwo$aa15.window, function(x) make.aa.15.window.from.31.sequence.window(input = x)  )
glimpse(G1G3_G4_RokaiTwo)

table(duplicated(G1G3_G4_RokaiTwo$aa15.window))#656
which(duplicated(G1G3_G4_RokaiTwo$aa15.window))

# if there are duplicate entries (e.g. because of _1_2_3) take the entry with the best p-value
unique.15aa.windows <- unique(G1G3_G4_RokaiTwo$aa15.window)
length(unique.15aa.windows) #6082
glimpse(G1G3_G4_RokaiTwo) #6738


G1G3_G4_RokaiTwo.2 <- c()
for (i in unique.15aa.windows){
  print(i)
  temp <- G1G3_G4_RokaiTwo %>% filter(aa15.window %in% c(i))
  temp <- temp %>% arrange(Limma.p.value)
  temp <- temp[1,]
  G1G3_G4_RokaiTwo.2 <- bind_rows(G1G3_G4_RokaiTwo.2, temp)
}
table(duplicated(G1G3_G4_RokaiTwo.2$aa15.window))

#calculate log2(FC)
G1G3_G4_RokaiTwo.2$log2FC <- log2(abs(G1G3_G4_RokaiTwo.2$FC_FR))
G1G3_G4_RokaiTwo.2$sign.FC_FR <- sign(G1G3_G4_RokaiTwo.2$FC_FR)
G1G3_G4_RokaiTwo.2$log2FC <- sign(G1G3_G4_RokaiTwo.2$FC_FR) * log2(abs(G1G3_G4_RokaiTwo.2$FC_FR))
glimpse(G1G3_G4_RokaiTwo.2) 

G1G3_G4_RokaiTwo.2 <- G1G3_G4_RokaiTwo.2 %>% arrange(desc(log2FC))
head(G1G3_G4_RokaiTwo.2$FC_FR)
tail(G1G3_G4_RokaiTwo.2$FC_FR)

### prepare RoKAI input and save data
glimpse(G1G3_G4_RokaiTwo.2) 
G1G3_G4_RokaiTwo.2.06.01.2022 <- G1G3_G4_RokaiTwo.2 %>% select(Protein, Position, log2FC) 
table(is.na(G1G3_G4_RokaiTwo.2.06.01.2022$Protein))
G1G3_G4_RokaiTwo.2.06.01.2022 <- G1G3_G4_RokaiTwo.2.06.01.2022 %>% filter(!is.na(Protein))
glimpse(G1G3_G4_RokaiTwo.2.06.01.2022)
table(is.na(G1G3_G4_RokaiTwo.2.06.01.2022$Protein))
table(is.na(G1G3_G4_RokaiTwo.2.06.01.2022$Position))
table(is.na(G1G3_G4_RokaiTwo.2.06.01.2022$log2FC))
colnames(G1G3_G4_RokaiTwo.2.06.01.2022) <- c("Protein", "Position", "Quantification")

#write_csv(G1G3_G4_RokaiTwo.2.06.01.2022,"G1G3_G4.RokaiTWO.input.limma.sep.2.group.comp.csv", col_names=TRUE)
#save(G1G3_G4_RokaiTwo  , file="G1G3_G4_RokaiTwo.limma.sep.2.group.comp.Rdata")


# perform analysis via https://rokai.io/
# tool settings was default: Reference Proteome: Uniprot Mouse / Fold changes: normalized / Kinase substrate dataset: PSP + Signor / RoKAI network: KS+PPI+SD+CoEv
# save RoKAI output kinase_table.csv and kinase_targets.csv and load these later



###################################################################################################################################
## G2_G4 DZ tumors IMAC Limma seperate 2 group comparison for RokaiTwo ##############################################################################

#Group 1
#"1_Fgfr2","2_Fgfr2", "13_Fgfr2-Ate1","14_Fgfr2-Ate1",

#Group 2
#"4_Fgfr2-dE18","5_Fgfr2-dE18","6_Fgfr2-dE18","23_Fgfr2-dE18-IGR1","24_Fgfr2-dE18-IGR1","25_Fgfr2-dE18-IGR2","26_Fgfr2-dE18-IGR2","29_Fgfr2-E18-C3","30_Fgfr2-E18-C3","31_Fgfr2-E18-C4","32_Fgfr2-E18-C4",

#Group 3
#"7_Fgfr2-Bicc1","8_Fgfr2-Bicc1","9_Fgfr2-Bicc1","17_Fgfr2-Tacc2","18_Fgfr2-Tacc2","19_Fgfr2-Tacc2",

#Group 4
#"15_Fgfr2-dE18-Ate1","16_Fgfr2-dE18-Ate1","10_Fgfr2-dE18-Bicc1","11_Fgfr2-dE18-Bicc1","12_Fgfr2-dE18-Bicc1","20_Fgfr2-dE18-Tacc2","21_Fgfr2-dE18-Tacc2","22_Fgfr2-dE18-Tacc2",

#Group 5
#"3_Fgfr2-E18-C2","27_Fgfr2-E18-C2","28_Fgfr2-E18-C2"

## select columns of interest
glimpse(DZ.tum.IMAC_expanded_log2_normclass1.PS_for.L.sep.2.group.comp)

G4_G2_Limma.2gc <- DZ.tum.IMAC_expanded_log2_normclass1.PS_for.L.sep.2.group.comp %>% select("ID.FR.all.C1.PS",
                                                                                             "id",                          "Proteins",                    "Positions.within.proteins",   "Leading.proteins",            "Protein",                    
                                                                                             "Protein.names",               "Gene.names",                  "Fasta.headers",               "Localization.prob",           "Number.of.Phospho.(STY)",    
                                                                                             "Amino.acid",                  "Sequence.window",             "Phospho.(STY).Probabilities", "Position.in.peptide",         "Reverse",                    
                                                                                             "Potential.contaminant",       "Positions",                   "Position",                    "Peptide.IDs",                 "Mod..peptide.IDs",   
                                                                                             "PS_Multiplicity", "Prot.Name_AApos", "Gene.Name_AApos",
                                                                                             
                                                                                             #Group 2
                                                                                             "4_Fgfr2-dE18","5_Fgfr2-dE18","6_Fgfr2-dE18","23_Fgfr2-dE18-IGR1","24_Fgfr2-dE18-IGR1","25_Fgfr2-dE18-IGR2","26_Fgfr2-dE18-IGR2","29_Fgfr2-E18-C3","30_Fgfr2-E18-C3","31_Fgfr2-E18-C4","32_Fgfr2-E18-C4",
                                                                                             
                                                                                             #Group 4
                                                                                             "15_Fgfr2-dE18-Ate1","16_Fgfr2-dE18-Ate1","10_Fgfr2-dE18-Bicc1","11_Fgfr2-dE18-Bicc1","12_Fgfr2-dE18-Bicc1","20_Fgfr2-dE18-Tacc2","21_Fgfr2-dE18-Tacc2","22_Fgfr2-dE18-Tacc2",
                                                                                             
                                                                                             data.presence.Group.2,
                                                                                             data.presence.Group.4
                                                                                             
)

#extend data presence calculation
G4_G2_Limma.2gc$data.presence.all <- G4_G2_Limma.2gc$data.presence.Group.2 + G4_G2_Limma.2gc$data.presence.Group.4
glimpse(G4_G2_Limma.2gc)

#remove complete zero rows
table(G4_G2_Limma.2gc$data.presence.all == 0) #F17555 T2641

G4_G2_Limma.2gc <- G4_G2_Limma.2gc %>% filter(data.presence.all > 0)
nrow(G4_G2_Limma.2gc) #17555

## X out of Y filter 
#Group 2: 11 samples => 50% => 5.5 samples => 6 samples
#Group 4: 8 samples => 50% => 4 samples => 5 samples
glimpse(G4_G2_Limma.2gc)
temp.left <- filter(G4_G2_Limma.2gc, data.presence.Group.2 >= 6) ; nrow(temp.left)    #6828
temp.right <- filter(G4_G2_Limma.2gc, data.presence.Group.4 >= 5) ; nrow(temp.right) #6058

temp.left.and.right <- bind_rows(temp.left , temp.right) ; nrow(temp.left.and.right)
temp.left.and.right <- distinct(temp.left.and.right) ; nrow(temp.left.and.right) #7905

G4_G2_Limma.2gc.XYfilter <- temp.left.and.right
nrow(G4_G2_Limma.2gc.XYfilter) #7905
length(unique(G4_G2_Limma.2gc.XYfilter$Prot.Name_AApos)) #7158

#add FRID
G4_G2_Limma.2gc.XYfilter$FRID <- 1:nrow(G4_G2_Limma.2gc.XYfilter) #here FRID is also row number
glimpse(G4_G2_Limma.2gc.XYfilter)

#find likely Limma NA pval situations 
#Group 2 90% => 10 samples
#Group 4 90% => 7 samples
possible.NApvals.A <- G4_G2_Limma.2gc.XYfilter %>% filter((data.presence.Group.2  >= 10 & data.presence.Group.4  == 0) | (data.presence.Group.2  == 0 & data.presence.Group.4  >= 7)  )

glimpse(possible.NApvals.A)
length(unique(possible.NApvals.A$FRID))
FRIDS.possible.NApvals <- possible.NApvals.A$FRID
FRIDS.possible.NApvals 

# for these filtered  black and white situation add zero so that p.val can be calculated
W <- G4_G2_Limma.2gc.XYfilter[FRIDS.possible.NApvals, ]
glimpse(W)
W <- W %>% mutate_at(vars(
  #Group 2
  "4_Fgfr2-dE18","5_Fgfr2-dE18","6_Fgfr2-dE18","23_Fgfr2-dE18-IGR1","24_Fgfr2-dE18-IGR1","25_Fgfr2-dE18-IGR2","26_Fgfr2-dE18-IGR2","29_Fgfr2-E18-C3","30_Fgfr2-E18-C3","31_Fgfr2-E18-C4","32_Fgfr2-E18-C4",
  
  #Group 4
  "15_Fgfr2-dE18-Ate1","16_Fgfr2-dE18-Ate1","10_Fgfr2-dE18-Bicc1","11_Fgfr2-dE18-Bicc1","12_Fgfr2-dE18-Bicc1","20_Fgfr2-dE18-Tacc2","21_Fgfr2-dE18-Tacc2","22_Fgfr2-dE18-Tacc2"),
  ~replace(., is.na(.), 0))
glimpse(W)

#combine again
G4_G2_Limma.2gc.XYfilter.2 <- G4_G2_Limma.2gc.XYfilter %>% filter(! FRID %in% FRIDS.possible.NApvals ) 
G4_G2_Limma.2gc.XYfilter.2 <- bind_rows(G4_G2_Limma.2gc.XYfilter.2, W)
G4_G2_Limma.2gc.XYfilter.2 <- G4_G2_Limma.2gc.XYfilter.2 %>% arrange(FRID)
glimpse(G4_G2_Limma.2gc.XYfilter.2)
table(duplicated(G4_G2_Limma.2gc.XYfilter.2$FRID))

##Limma statistics
#get data 
data.for.limma <- G4_G2_Limma.2gc.XYfilter.2 %>% select(
  #Group 2
  "4_Fgfr2-dE18","5_Fgfr2-dE18","6_Fgfr2-dE18","23_Fgfr2-dE18-IGR1","24_Fgfr2-dE18-IGR1","25_Fgfr2-dE18-IGR2","26_Fgfr2-dE18-IGR2","29_Fgfr2-E18-C3","30_Fgfr2-E18-C3","31_Fgfr2-E18-C4","32_Fgfr2-E18-C4",
  
  #Group 4
  "15_Fgfr2-dE18-Ate1","16_Fgfr2-dE18-Ate1","10_Fgfr2-dE18-Bicc1","11_Fgfr2-dE18-Bicc1","12_Fgfr2-dE18-Bicc1","20_Fgfr2-dE18-Tacc2","21_Fgfr2-dE18-Tacc2","22_Fgfr2-dE18-Tacc2"
)
rownames(data.for.limma) <- paste0(G4_G2_Limma.2gc.XYfilter.2$ID.FR.all.C1.PS, "__", G4_G2_Limma.2gc.XYfilter.2$Gene.Name_AApos) 
nrow(data.for.limma)
glimpse(data.for.limma)


design <- cbind(left=c(rep(1, 11), rep(0, 8)), right=c(rep(0,11), rep(1, 8))); design
fit <- lmFit(data.for.limma, design) # Warning: Partial NA coefficients for 283 probe(s) 
cont.matrix <- makeContrasts( left.vs.right =  right - left, levels=design) 
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
topTable(fit2, adjust="BH")

Limma.p.value <- as.numeric(fit2$p.value)
glimpse(Limma.p.value)
Limma.adj.p.valueBH <- p.adjust(Limma.p.value, "BH") #
glimpse(Limma.adj.p.valueBH)
Limma.t <- as.numeric(fit2$t)
glimpse(Limma.t)

G4_G2_Limma.2gc.XYfilter.2.B <- G4_G2_Limma.2gc.XYfilter.2
G4_G2_Limma.2gc.XYfilter.2.B$Limma.p.value <- as.numeric(Limma.p.value)
G4_G2_Limma.2gc.XYfilter.2.B$Limma.adj.p.valueBH <- as.numeric(Limma.adj.p.valueBH)
G4_G2_Limma.2gc.XYfilter.2.B$Limma.t <- as.numeric(Limma.t)

glimpse(G4_G2_Limma.2gc.XYfilter.2.B)

###calculate FC
for(i in 1:nrow(G4_G2_Limma.2gc.XYfilter.2.B)){
  print(i)
  leftvals <- as.numeric(G4_G2_Limma.2gc.XYfilter.2.B[i,c("4_Fgfr2-dE18","5_Fgfr2-dE18","6_Fgfr2-dE18","23_Fgfr2-dE18-IGR1","24_Fgfr2-dE18-IGR1","25_Fgfr2-dE18-IGR2","26_Fgfr2-dE18-IGR2","29_Fgfr2-E18-C3","30_Fgfr2-E18-C3","31_Fgfr2-E18-C4","32_Fgfr2-E18-C4")]); leftvals #
  rightvals <- as.numeric(G4_G2_Limma.2gc.XYfilter.2.B[i,c("15_Fgfr2-dE18-Ate1","16_Fgfr2-dE18-Ate1","10_Fgfr2-dE18-Bicc1","11_Fgfr2-dE18-Bicc1","12_Fgfr2-dE18-Bicc1","20_Fgfr2-dE18-Tacc2","21_Fgfr2-dE18-Tacc2","22_Fgfr2-dE18-Tacc2")]); rightvals #
  
  FC <- FR_FC_LOGData_2.3.17(control = leftvals, treatment = rightvals, LOGtype = 2); FC
  
  G4_G2_Limma.2gc.XYfilter.2.B$FC_FR[i] <- FC
}
glimpse(G4_G2_Limma.2gc.XYfilter.2.B)

#check numbers unique Sequence.window
glimpse(G4_G2_Limma.2gc.XYfilter.2.B)#7905
G4_G2_Limma.2gc.XYfilter.2.B %>% filter(Limma.p.value < 0.05) %>% select(Sequence.window) %>% distinct() %>% nrow() #2438 total sign
G4_G2_Limma.2gc.XYfilter.2.B %>% filter(Limma.p.value < 0.05) %>% filter(FC_FR > 1) %>% select(Sequence.window) %>% distinct() %>% nrow() #1204 up
G4_G2_Limma.2gc.XYfilter.2.B %>% filter(Limma.p.value < 0.05) %>% filter(FC_FR < 1) %>% select(Sequence.window) %>% distinct() %>% nrow() #1236 down

# filter na pvalues
table(is.na(G4_G2_Limma.2gc.XYfilter.2.B$Limma.p.value))
G4_G2_Limma.2gc.XYfilter.2.B <- G4_G2_Limma.2gc.XYfilter.2.B %>% filter(!is.na(Limma.p.value))
glimpse(G4_G2_Limma.2gc.XYfilter.2.B) #7622
table(is.na(G4_G2_Limma.2gc.XYfilter.2.B$Limma.p.value))

#check +/- 10000 FCs
table(G4_G2_Limma.2gc.XYfilter.2.B$FC_FR >= 10000) #27
table(G4_G2_Limma.2gc.XYfilter.2.B$FC_FR <= -10000) #18
summary(G4_G2_Limma.2gc.XYfilter.2.B$FC_FR)

#find second min or max and add +/-5 for +/- 10000 FCs
temp <- G4_G2_Limma.2gc.XYfilter.2.B %>% filter(!FC_FR %in% c(10000, -10000)) %>% arrange(desc(FC_FR))
head(temp$FC_FR)
tail(temp$FC_FR)
max(temp$FC_FR)
min(temp$FC_FR)

G4_G2_Limma.2gc.XYfilter.2.B$FC_FR <- replace(G4_G2_Limma.2gc.XYfilter.2.B$FC_FR, G4_G2_Limma.2gc.XYfilter.2.B$FC_FR==10000, max(temp$FC_FR)+5)
G4_G2_Limma.2gc.XYfilter.2.B$FC_FR <- replace(G4_G2_Limma.2gc.XYfilter.2.B$FC_FR, G4_G2_Limma.2gc.XYfilter.2.B$FC_FR==-10000, min(temp$FC_FR)-5)
summary(G4_G2_Limma.2gc.XYfilter.2.B$FC_FR)
glimpse(G4_G2_Limma.2gc.XYfilter.2.B) #7622

##G4_G2_Rokai reshape for RoKAI
G2_G4_RokaiTwo <- G4_G2_Limma.2gc.XYfilter.2.B

# make aa15.window
# take first entry sequence window if there are two entries
G2_G4_RokaiTwo$aa15.window <- pbsapply(G2_G4_RokaiTwo$Sequence.window, function(x) unlist(str_split(string =x, pattern = ";"))[1]  )
G2_G4_RokaiTwo$aa15.window <- pbsapply(G2_G4_RokaiTwo$aa15.window, function(x) make.aa.15.window.from.31.sequence.window(input = x)  )
glimpse(G2_G4_RokaiTwo)

table(duplicated(G2_G4_RokaiTwo$aa15.window))
which(duplicated(G2_G4_RokaiTwo$aa15.window))

# if there are duplicate entries (e.g. because of _1_2_3) take the entry with the best p-value
unique.15aa.windows <- unique(G2_G4_RokaiTwo$aa15.window)
length(unique.15aa.windows) #6864
glimpse(G2_G4_RokaiTwo) #7622

G2_G4_RokaiTwo.2 <- c()
for (i in unique.15aa.windows){
  print(i)
  temp <- G2_G4_RokaiTwo %>% filter(aa15.window %in% c(i))
  temp <- temp %>% arrange(Limma.p.value)
  temp <- temp[1,]
  G2_G4_RokaiTwo.2 <- bind_rows(G2_G4_RokaiTwo.2, temp)
}
table(duplicated(G2_G4_RokaiTwo.2$aa15.window))

#calculate log2(FC)
G2_G4_RokaiTwo.2$log2FC <- log2(abs(G2_G4_RokaiTwo.2$FC_FR))
G2_G4_RokaiTwo.2$sign.FC_FR <- sign(G2_G4_RokaiTwo.2$FC_FR)
G2_G4_RokaiTwo.2$log2FC <- sign(G2_G4_RokaiTwo.2$FC_FR) * log2(abs(G2_G4_RokaiTwo.2$FC_FR))
glimpse(G2_G4_RokaiTwo.2) 

G2_G4_RokaiTwo.2 <- G2_G4_RokaiTwo.2 %>% arrange(desc(log2FC))
head(G2_G4_RokaiTwo.2$FC_FR)
tail(G2_G4_RokaiTwo.2$FC_FR)

### prepare RoKAI input and save data
glimpse(G2_G4_RokaiTwo.2) 
G2_G4_RokaiTwo.2.06.01.2022 <- G2_G4_RokaiTwo.2 %>% select(Protein, Position, log2FC) 
table(is.na(G2_G4_RokaiTwo.2.06.01.2022$Protein))
G2_G4_RokaiTwo.2.06.01.2022 <- G2_G4_RokaiTwo.2.06.01.2022 %>% filter(!is.na(Protein))
glimpse(G2_G4_RokaiTwo.2.06.01.2022)
table(is.na(G2_G4_RokaiTwo.2.06.01.2022$Protein))
table(is.na(G2_G4_RokaiTwo.2.06.01.2022$Position))
table(is.na(G2_G4_RokaiTwo.2.06.01.2022$log2FC))
colnames(G2_G4_RokaiTwo.2.06.01.2022) <- c("Protein", "Position", "Quantification")


#write_csv(G2_G4_RokaiTwo.2.06.01.2022,"G2_G4.RokaiTWO.input.limma.sep.2.group.comp.csv", col_names=TRUE)
#save(G2_G4_RokaiTwo  , file="G2_G4_RokaiTwo.limma.sep.2.group.comp.Rdata")


#http://rokai.io
#tool settings was default: Reference Proteome: Uniprot Mouse / Fold changes: normalized / Kinase substrate dataset: PSP + Signor / RoKAI network: KS+PPI+SD+CoEv.

# perform analysis via https://rokai.io/
# tool settings was default: Reference Proteome: Uniprot Mouse / Fold changes: normalized / Kinase substrate dataset: PSP + Signor / RoKAI network: KS+PPI+SD+CoEv
# save RoKAI output kinase_table.csv and kinase_targets.csv and load these later


############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################
### tumors IMAC: get RoKAI outputs and combine RoKAI results kinase table
### filter for number of substrates >= 2  for each comparison

#
G1G3_G2_rokai_kinase.table <- fread("G1G3_G2_RokaiTwo_kinase_table.csv")
G1G3_G2_rokai_kinase.table <- G1G3_G2_rokai_kinase.table  %>% filter(NumSubs >= 2)
colnames(G1G3_G2_rokai_kinase.table) <- c("KinaseID","Name","Gene",  paste0("G1G3_G2", colnames(G1G3_G2_rokai_kinase.table[, c(4:9)])))
glimpse(G1G3_G2_rokai_kinase.table)

#
G1G3_G4_rokai_kinase.table <- fread("G1G3_G4_RokaiTwo_kinase_table.csv")
G1G3_G4_rokai_kinase.table <- G1G3_G4_rokai_kinase.table  %>% filter(NumSubs >= 2)
colnames(G1G3_G4_rokai_kinase.table) <- c("KinaseID","Name","Gene",  paste0("G1G3_G4", colnames(G1G3_G4_rokai_kinase.table[, c(4:9)])))
glimpse(G1G3_G4_rokai_kinase.table)

#
G2_G4_rokai_kinase.table <- fread("G2_G4_RokaiTwo_kinase_table.csv")
G2_G4_rokai_kinase.table <- G2_G4_rokai_kinase.table  %>% filter(NumSubs >= 2)
colnames(G2_G4_rokai_kinase.table) <- c("KinaseID","Name","Gene",  paste0("G2_G4", colnames(G2_G4_rokai_kinase.table[, c(4:9)])))
glimpse(G2_G4_rokai_kinase.table)


#get all KinaseIDs corresponding Names and Genes and remove duplicates
comb.ROKAI.KinaseID <- bind_rows(G1G3_G2_rokai_kinase.table %>% select(KinaseID, Name , Gene),
                                 G1G3_G4_rokai_kinase.table  %>% select(KinaseID, Name , Gene),
                                 G2_G4_rokai_kinase.table %>% select(KinaseID, Name , Gene)
) %>% distinct()
glimpse(comb.ROKAI.KinaseID) #89
comb.ROKAI.KinaseID$Gene

#merge data 
comb.ROKAI.result <- merge(comb.ROKAI.KinaseID, G1G3_G2_rokai_kinase.table[, c(1, 4:9)], by.x = "KinaseID", by.y = "KinaseID", all.x = T, all.y = T, sort = F)
comb.ROKAI.result <- merge(comb.ROKAI.result, G1G3_G4_rokai_kinase.table[, c(1, 4:9)], by.x = "KinaseID", by.y = "KinaseID", all.x = T, all.y = T, sort = F)
comb.ROKAI.result <- merge(comb.ROKAI.result, G2_G4_rokai_kinase.table[, c(1, 4:9)], by.x = "KinaseID", by.y = "KinaseID", all.x = T, all.y = T, sort = F)
glimpse(comb.ROKAI.result)

#plot kinase activity scores
plot.comb.ROKAI.result <- comb.ROKAI.result %>% arrange(desc(G1G3_G2ZScore))  %>% select(Gene,  G1G3_G2ZScore, G1G3_G4ZScore , G2_G4ZScore)
plot.comb.ROKAI.result$sum.z.scores <- apply(plot.comb.ROKAI.result[, c("G1G3_G2ZScore", "G1G3_G4ZScore" , "G2_G4ZScore")], 1, function(x) sum(x, na.rm=T))
plot.comb.ROKAI.result <- plot.comb.ROKAI.result %>%  arrange(desc(sum.z.scores)) %>% select(-sum.z.scores)
glimpse(plot.comb.ROKAI.result)

#order for plot
order.plot.comb.ROKAI.result <- plot.comb.ROKAI.result$Gene
order.plot.comb.ROKAI.result #

#reshape for plot 
melted.plot.comb.ROKAI.result <- reshape2::melt(plot.comb.ROKAI.result) 
glimpse(melted.plot.comb.ROKAI.result)

#plot heatmap kinase activity scores
ggplot(melted.plot.comb.ROKAI.result  , aes(x=variable, y=Gene, fill=value)) + 
  geom_tile(color = "white") + 
  scale_fill_viridis( option = "turbo", na.value = "grey60", name="z.score") + 
  coord_equal()+
  scale_y_discrete(limits= rev(order.plot.comb.ROKAI.result), expand=c(0,0))+
  theme(legend.position="right", legend.justification = "center")+
  ggtitle("IMAC tumors - RoKAITwo Kinase activity z-score ") +
  theme(plot.title = element_text(hjust = 0.5, face="bold", size=16)) + 
  theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1))+
  xlab(NULL) + 
  ylab(NULL) +
  theme(axis.text.x = element_text(face= "plain", colour="black", size=10)) 


################################################################################################################################################################
################################################################################################################################################################
### tumors IMAC, RoKAI outputs: combine and visualize kinase targets

glimpse(DZ.tum.IMAC_expanded_log2_normclass1.PS_for.L.sep.2.group.comp)

kinases.of.interest <- order.plot.comb.ROKAI.result
length(kinases.of.interest) #89
kinases.of.interest

#load Rokai Kinase target results
G1G3_G2_rokai_kinase.targets <- fread("G1G3_G2_RokaiTwo_kinase_targets.csv")
colnames(G1G3_G2_rokai_kinase.targets) <- c("KinID","KinName","KinGene","SubsProtein","SubsGene","Position","DataSource","Flanking" ,  paste0("G1G3_G2_", colnames(G1G3_G2_rokai_kinase.targets[, c(9:12)])))
glimpse(G1G3_G2_rokai_kinase.targets)


G1G3_G4_rokai_kinase.targets <- fread("G1G3_G4_RokaiTwo_kinase_targets.csv")
colnames(G1G3_G4_rokai_kinase.targets ) <- c("KinID","KinName","KinGene","SubsProtein","SubsGene","Position","DataSource","Flanking" ,  paste0("G1G3_G4_", colnames(G1G3_G4_rokai_kinase.targets [, c(9:12)])))
glimpse(G1G3_G4_rokai_kinase.targets )


G2_G4_rokai_kinase.targets <- fread("G2_G4_RokaiTwo_kinase_targets.csv")
colnames(G2_G4_rokai_kinase.targets) <- c("KinID","KinName","KinGene","SubsProtein","SubsGene","Position","DataSource","Flanking" ,  paste0("G2_G4_", colnames(G2_G4_rokai_kinase.targets[, c(9:12)])))
glimpse(G2_G4_rokai_kinase.targets)



#get all KinaseIDs, targtes corresponding Names and Genes and remove duplicates
comb.ROKAI.targets <- bind_rows(G1G3_G2_rokai_kinase.targets %>% select(KinID, KinName, KinGene, SubsProtein, SubsGene, Position, DataSource, Flanking),
                                G1G3_G4_rokai_kinase.targets %>% select(KinID, KinName, KinGene, SubsProtein, SubsGene, Position, DataSource, Flanking),
                                G2_G4_rokai_kinase.targets   %>% select(KinID, KinName, KinGene, SubsProtein, SubsGene, Position, DataSource, Flanking)) %>% distinct()

glimpse(comb.ROKAI.targets)

#make flanking sequence all capital letters
comb.ROKAI.targets$Flanking <- str_to_upper(comb.ROKAI.targets$Flanking)
glimpse(comb.ROKAI.targets)


# filter for Kinases of interest and get phosphosite sequence
comb.ROKAI.targets_2 <- comb.ROKAI.targets %>% filter(KinGene %in% kinases.of.interest)
glimpse(comb.ROKAI.targets_2)

comb.ROKAI.targets_2_unique.flanking <- unique(comb.ROKAI.targets_2 %>% select(Flanking) %>% pull())
glimpse(comb.ROKAI.targets_2_unique.flanking)


###load two group comparison data (see above) and select sign sites with FC >=1.5 or FC <= -1.5 ###

#load("G1G3_G2_RokaiTwo.limma.sep.2.group.comp.Rdata")
glimpse(G1G3_G2_RokaiTwo) 
G1G3_G2_RokaiTwo_2group.filter.sign.FC <- G1G3_G2_RokaiTwo %>% filter(Limma.p.value < 0.05) %>% filter(FC_FR >=1.5 | FC_FR <= -1.5) %>% select(ID.FR.all.C1.PS) %>% pull()
glimpse(G1G3_G2_RokaiTwo_2group.filter.sign.FC )

#load("G1G3_G4_RokaiTwo.limma.sep.2.group.comp.Rdata")
glimpse(G1G3_G4_RokaiTwo)
G1G3_G4_RokaiTwo_2group.filter.sign.FC <- G1G3_G4_RokaiTwo  %>% filter(Limma.p.value < 0.05) %>% filter(FC_FR >=1.5 | FC_FR <= -1.5) %>% select(ID.FR.all.C1.PS) %>% pull()
glimpse(G1G3_G4_RokaiTwo_2group.filter.sign.FC)

#load("G2_G4_RokaiTwo.limma.sep.2.group.comp.Rdata")
glimpse(G2_G4_RokaiTwo)
G2_G4_RokaiTwo_2group.filter.sign.FC <- G2_G4_RokaiTwo %>% filter(Limma.p.value < 0.05) %>% filter(FC_FR >=1.5 | FC_FR <= -1.5) %>% select(ID.FR.all.C1.PS) %>% pull()
glimpse(G2_G4_RokaiTwo_2group.filter.sign.FC)


### get unique candidate IDs 
all.2group.filter.sign.FC.C1.IDs <- unique( c(G1G3_G2_RokaiTwo_2group.filter.sign.FC,
                                              G1G3_G4_RokaiTwo_2group.filter.sign.FC,
                                              G2_G4_RokaiTwo_2group.filter.sign.FC
) )
glimpse(all.2group.filter.sign.FC.C1.IDs)



###get all the intensity data for all samples and the selected sites (all.2group.filter.sign.FC.C1.IDs)
glimpse(DZ.tum.IMAC_expanded_log2_normclass1.PS_for.L.sep.2.group.comp)

DZ.TUMORS.IMAC.MBR.ON_expanded_log2_norm.class1.PS_filter.sign.FC <- DZ.tum.IMAC_expanded_log2_normclass1.PS_for.L.sep.2.group.comp %>% filter(ID.FR.all.C1.PS %in% all.2group.filter.sign.FC.C1.IDs)

glimpse(DZ.TUMORS.IMAC.MBR.ON_expanded_log2_norm.class1.PS_filter.sign.FC ) #2834

# make aa15.window
# take first entry sequence window if there are two entries
DZ.TUMORS.IMAC.MBR.ON_expanded_log2_norm.class1.PS_filter.sign.FC$aa15.window <- pbsapply(DZ.TUMORS.IMAC.MBR.ON_expanded_log2_norm.class1.PS_filter.sign.FC$Sequence.window, function(x) unlist(str_split(string =x, pattern = ";"))[1]  )
DZ.TUMORS.IMAC.MBR.ON_expanded_log2_norm.class1.PS_filter.sign.FC$aa15.window <- pbsapply(DZ.TUMORS.IMAC.MBR.ON_expanded_log2_norm.class1.PS_filter.sign.FC$aa15.window, function(x) make.aa.15.window.from.31.sequence.window(input = x)  )

glimpse(DZ.TUMORS.IMAC.MBR.ON_expanded_log2_norm.class1.PS_filter.sign.FC)


# filter for the combined RoKAI target sites
glimpse(comb.ROKAI.targets_2_unique.flanking)

DZ.TUMORS.IMAC.MBR.ON_expanded_log2_norm.class1.PS_filter.sign.FC_rokai.filter <- DZ.TUMORS.IMAC.MBR.ON_expanded_log2_norm.class1.PS_filter.sign.FC %>% filter(aa15.window %in% comb.ROKAI.targets_2_unique.flanking)
glimpse(DZ.TUMORS.IMAC.MBR.ON_expanded_log2_norm.class1.PS_filter.sign.FC_rokai.filter)


#merge information on upstream kinase
glimpse(comb.ROKAI.targets_2)

DZ.TUMORS.IMAC.MBR.ON_expanded_log2_norm.class1.PS_filter.sign.FC_rokai.filter_merge.rokai.upstream.Kins <- merge(DZ.TUMORS.IMAC.MBR.ON_expanded_log2_norm.class1.PS_filter.sign.FC_rokai.filter, comb.ROKAI.targets_2, by.x = "aa15.window", by.y = "Flanking", all.x = T, all.y = F, sort = F)
glimpse(DZ.TUMORS.IMAC.MBR.ON_expanded_log2_norm.class1.PS_filter.sign.FC_rokai.filter_merge.rokai.upstream.Kins)

#merge Limma p values from seperate two group comparisons #
glimpse(DZ.tum.IMAC_expanded_log2_normclass1.PS_for.L.sep.2.group.comp)
glimpse(G1G3_G2_RokaiTwo)#
glimpse(G1G3_G4_RokaiTwo)#
glimpse(G2_G4_RokaiTwo)#



DZ.TUMORS.IMAC.MBR.ON_expanded_log2_norm.class1.PS_filter.sign.FC_rokai.filter_merge.rokai.upstream.Kins_2 <- merge(DZ.TUMORS.IMAC.MBR.ON_expanded_log2_norm.class1.PS_filter.sign.FC_rokai.filter_merge.rokai.upstream.Kins, 
                                                                                                                    G1G3_G2_RokaiTwo %>% select("ID.FR.all.C1.PS", "Limma.p.value", "Limma.adj.p.valueBH", "FC_FR") %>% rename(Limma.p.value_G1G3_G2 = Limma.p.value, Limma.adj.p.valueBH_G1G3_G2 = Limma.adj.p.valueBH, FC_FR_G1G3_G2 = FC_FR), 
                                                                                                                    by.x = "ID.FR.all.C1.PS", by.y = "ID.FR.all.C1.PS", all.x = T, all.y = F, sort = F)

DZ.TUMORS.IMAC.MBR.ON_expanded_log2_norm.class1.PS_filter.sign.FC_rokai.filter_merge.rokai.upstream.Kins_2 <- merge(DZ.TUMORS.IMAC.MBR.ON_expanded_log2_norm.class1.PS_filter.sign.FC_rokai.filter_merge.rokai.upstream.Kins_2, 
                                                                                                                    G1G3_G4_RokaiTwo %>% select("ID.FR.all.C1.PS", "Limma.p.value", "Limma.adj.p.valueBH", "FC_FR") %>% rename(Limma.p.value_G1G3_G4 = Limma.p.value, Limma.adj.p.valueBH_G1G3_G4 = Limma.adj.p.valueBH, FC_FR_G1G3_G4 = FC_FR), 
                                                                                                                    by.x = "ID.FR.all.C1.PS", by.y = "ID.FR.all.C1.PS", all.x = T, all.y = F, sort = F)

DZ.TUMORS.IMAC.MBR.ON_expanded_log2_norm.class1.PS_filter.sign.FC_rokai.filter_merge.rokai.upstream.Kins_2 <- merge(DZ.TUMORS.IMAC.MBR.ON_expanded_log2_norm.class1.PS_filter.sign.FC_rokai.filter_merge.rokai.upstream.Kins_2, 
                                                                                                                    G2_G4_RokaiTwo %>% select("ID.FR.all.C1.PS", "Limma.p.value", "Limma.adj.p.valueBH", "FC_FR") %>% rename(Limma.p.value_G2_G4 = Limma.p.value, Limma.adj.p.valueBH_G2_G4 = Limma.adj.p.valueBH, FC_FR_G2_G4 = FC_FR), 
                                                                                                                    by.x = "ID.FR.all.C1.PS", by.y = "ID.FR.all.C1.PS", all.x = T, all.y = F, sort = F)

#add data group 5
glimpse(DZ.tum.IMAC_expanded_log2_normclass1.PS) #20,196

df.to.add.group5 <- DZ.tum.IMAC_expanded_log2_normclass1.PS %>% select(id : Mod..peptide.IDs, PS_Multiplicity, "Norm.Intensity.A",  "Norm.Intensity.N",  "Norm.Intensity.F",  "Norm.Intensity.R",  "Norm.Intensity.D",  "Norm.Intensity.P",  "Norm.Intensity.AB",
                                                                       "Norm.Intensity.H",  "Norm.Intensity.T",  "Norm.Intensity.AD", "Norm.Intensity.C", "Norm.Intensity.O",  "Norm.Intensity.AA", "Norm.Intensity.G", 
                                                                       "Norm.Intensity.S",  "Norm.Intensity.E", "Norm.Intensity.Q",  "Norm.Intensity.AC", "Norm.Intensity.I",  "Norm.Intensity.U",  "Norm.Intensity.AE",
                                                                       "Norm.Intensity.J",  "Norm.Intensity.V",  "Norm.Intensity.K",  "Norm.Intensity.W",  "Norm.Intensity.B", "Norm.Intensity.X",  "Norm.Intensity.AF",
                                                                       "Norm.Intensity.L",  "Norm.Intensity.Y",  "Norm.Intensity.M",  "Norm.Intensity.Z",  Prot.Name_AApos :  ID.FR.all.C1.PS )
glimpse(df.to.add.group5 )
colnames(df.to.add.group5 ) <- c("id",                          "Proteins",                    "Positions.within.proteins",   "Leading.proteins",            "Protein",                    
                                 "Protein.names",               "Gene.names",                  "Fasta.headers",               "Localization.prob",           "Number.of.Phospho.(STY)",    
                                 "Amino.acid",                  "Sequence.window",             "Phospho.(STY).Probabilities", "Position.in.peptide",         "Reverse",                    
                                 "Potential.contaminant",       "Positions",                   "Position",                    "Peptide.IDs",                 "Mod..peptide.IDs",           
                                 "PS_Multiplicity",
                                 "1_Fgfr2",
                                 "2_Fgfr2",
                                 
                                 "13_Fgfr2-Ate1",
                                 "14_Fgfr2-Ate1",
                                 
                                 "7_Fgfr2-Bicc1",
                                 "8_Fgfr2-Bicc1",
                                 "9_Fgfr2-Bicc1",
                                 
                                 "17_Fgfr2-Tacc2",
                                 "18_Fgfr2-Tacc2",
                                 "19_Fgfr2-Tacc2",
                                 
                                 "4_Fgfr2-dE18",
                                 "5_Fgfr2-dE18",
                                 "6_Fgfr2-dE18",
                                 
                                 "15_Fgfr2-dE18-Ate1",
                                 "16_Fgfr2-dE18-Ate1",
                                 
                                 "10_Fgfr2-dE18-Bicc1",
                                 "11_Fgfr2-dE18-Bicc1",
                                 "12_Fgfr2-dE18-Bicc1",
                                 
                                 "20_Fgfr2-dE18-Tacc2",
                                 "21_Fgfr2-dE18-Tacc2",
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
                                 
                                 "31_Fgfr2-E18-C4" ,
                                 "32_Fgfr2-E18-C4",
                                 "Prot.Name_AApos",
                                 "Gene.Name_AApos",            
                                 "ID.FR.all.C1.PS")
glimpse(df.to.add.group5)

df.to.add.group5 <- df.to.add.group5 %>% select("ID.FR.all.C1.PS", "3_Fgfr2-E18-C2","27_Fgfr2-E18-C2","28_Fgfr2-E18-C2")
glimpse(df.to.add.group5)




DZ.TUMORS.IMAC.MBR.ON_expanded_log2_norm.class1.PS_filter.sign.FC_rokai.filter_merge.rokai.upstream.Kins_2 <- merge(DZ.TUMORS.IMAC.MBR.ON_expanded_log2_norm.class1.PS_filter.sign.FC_rokai.filter_merge.rokai.upstream.Kins_2, 
                                                                                                                    df.to.add.group5 %>% select("ID.FR.all.C1.PS", "3_Fgfr2-E18-C2","27_Fgfr2-E18-C2","28_Fgfr2-E18-C2"), 
                                                                                                                    by.x = "ID.FR.all.C1.PS", by.y = "ID.FR.all.C1.PS", all.x = T, all.y = F, sort = F)



glimpse(DZ.TUMORS.IMAC.MBR.ON_expanded_log2_norm.class1.PS_filter.sign.FC_rokai.filter_merge.rokai.upstream.Kins_2)


#prepare a unique rowname to be able to plot everything
DZ.TUMORS.IMAC.MBR.ON_expanded_log2_norm.class1.PS_filter.sign.FC_rokai.filter_merge.rokai.upstream.Kins_2$Gene.Name_AApos.collapsed <- pbsapply(DZ.TUMORS.IMAC.MBR.ON_expanded_log2_norm.class1.PS_filter.sign.FC_rokai.filter_merge.rokai.upstream.Kins_2$Gene.Name_AApos, function(x) str_remove_all(string=x, pattern="x_\\d") )
glimpse(DZ.TUMORS.IMAC.MBR.ON_expanded_log2_norm.class1.PS_filter.sign.FC_rokai.filter_merge.rokai.upstream.Kins_2)

DZ.TUMORS.IMAC.MBR.ON_expanded_log2_norm.class1.PS_filter.sign.FC_rokai.filter_merge.rokai.upstream.Kins_2$row_Gene.Name_AApos <- paste0(paste0(rownames(DZ.TUMORS.IMAC.MBR.ON_expanded_log2_norm.class1.PS_filter.sign.FC_rokai.filter_merge.rokai.upstream.Kins_2), "_", DZ.TUMORS.IMAC.MBR.ON_expanded_log2_norm.class1.PS_filter.sign.FC_rokai.filter_merge.rokai.upstream.Kins_2$Gene.Name_AApos))
glimpse(DZ.TUMORS.IMAC.MBR.ON_expanded_log2_norm.class1.PS_filter.sign.FC_rokai.filter_merge.rokai.upstream.Kins_2)





## prepare heatmap facet plot for results combined RoKAI #################################################################################

facet.hm.plot <- DZ.TUMORS.IMAC.MBR.ON_expanded_log2_norm.class1.PS_filter.sign.FC_rokai.filter_merge.rokai.upstream.Kins_2 %>% select(
  #"row_Gene.Name_AApos", 
  "Gene.Name_AApos", 
  "KinGene",
  #Group 1
  "1_Fgfr2","2_Fgfr2", "13_Fgfr2-Ate1","14_Fgfr2-Ate1",
  
  #Group 2
  "4_Fgfr2-dE18","5_Fgfr2-dE18","6_Fgfr2-dE18","23_Fgfr2-dE18-IGR1","24_Fgfr2-dE18-IGR1","25_Fgfr2-dE18-IGR2","26_Fgfr2-dE18-IGR2","29_Fgfr2-E18-C3","30_Fgfr2-E18-C3","31_Fgfr2-E18-C4","32_Fgfr2-E18-C4",
  
  #Group 3
  "7_Fgfr2-Bicc1","8_Fgfr2-Bicc1","9_Fgfr2-Bicc1","17_Fgfr2-Tacc2","18_Fgfr2-Tacc2","19_Fgfr2-Tacc2",
  
  #Group 4
  "15_Fgfr2-dE18-Ate1","16_Fgfr2-dE18-Ate1","10_Fgfr2-dE18-Bicc1","11_Fgfr2-dE18-Bicc1","12_Fgfr2-dE18-Bicc1","20_Fgfr2-dE18-Tacc2","21_Fgfr2-dE18-Tacc2","22_Fgfr2-dE18-Tacc2",
  
  #Group 5
  "3_Fgfr2-E18-C2","27_Fgfr2-E18-C2","28_Fgfr2-E18-C2"
)

glimpse(facet.hm.plot)

#get overview how many observations per upstream Kinase  (some sites duplicated per kinase term)
print(facet.hm.plot %>% group_by(KinGene) %>% summarize(count = n()) %>% arrange(desc(count)), n=68 )#

#define order via number of substrates per upstream kinase
order.upstreamKin.via.num.subs <- facet.hm.plot %>% group_by(KinGene) %>% summarize(count = n()) %>% arrange(desc(count)) %>% select(KinGene) %>% pull()
order.upstreamKin.via.num.subs
order.upstreamKin.via.num.subs[1] #the kinase with the highes substrates here: Mapk1

# remove duplicate PS per upstream Kinase start with Kinase that has highest substrate number
glimpse(facet.hm.plot)
facet.hm.plot_2 <- facet.hm.plot %>% filter(KinGene %in% c(order.upstreamKin.via.num.subs[1]))
glimpse(facet.hm.plot_2)
temp.substrates.start <- facet.hm.plot %>% filter(KinGene %in% c(order.upstreamKin.via.num.subs[1])) %>% select(Gene.Name_AApos) %>% pull()
temp.substrates.start

for(i in order.upstreamKin.via.num.subs[2:length(order.upstreamKin.via.num.subs)]){
  print(i)
  temp.substrates.next.kin <- facet.hm.plot %>% filter(KinGene %in% c(i)) %>% select(Gene.Name_AApos) %>% pull()
  temp.substrates.next.kin.include <- setdiff(temp.substrates.next.kin, temp.substrates.start)
  
  temp.df.include <- facet.hm.plot %>% filter(KinGene %in% c(i)) %>% filter(Gene.Name_AApos %in% temp.substrates.next.kin.include)
  
  facet.hm.plot_2 <- bind_rows(facet.hm.plot_2, temp.df.include)
  
  temp.substrates.next.kin.exclude <- setdiff( temp.substrates.start, temp.substrates.next.kin)
  
  temp.substrates.start <- c(temp.substrates.start, temp.substrates.next.kin.include)
}
glimpse(facet.hm.plot_2)

glimpse(facet.hm.plot_2)
table(duplicated(facet.hm.plot_2$Gene.Name_AApos))

###for better visualization collapse _1/_2/_3 sites
###collapsing _1/_2/_3 via using sum
# note: sum is based on _1/_2/_3 intensities in the dataframe at hand => unfiltered data might contain more _1/_2/_3 PS that are not included in the sum here
# unlog, collapse, relog


facet.hm.plot_2$Gene.Name_AApos.collapsed <- pbsapply(facet.hm.plot_2$Gene.Name_AApos, function(x) str_remove_all(string=x, pattern="x_\\d") )
glimpse(facet.hm.plot_2)
## unlog, do collapsing, relog
facet.hm.plot_2[, c( 
  #Group 1
  "1_Fgfr2","2_Fgfr2", "13_Fgfr2-Ate1","14_Fgfr2-Ate1",
  #Group 2
  "4_Fgfr2-dE18","5_Fgfr2-dE18","6_Fgfr2-dE18","23_Fgfr2-dE18-IGR1","24_Fgfr2-dE18-IGR1","25_Fgfr2-dE18-IGR2","26_Fgfr2-dE18-IGR2","29_Fgfr2-E18-C3","30_Fgfr2-E18-C3","31_Fgfr2-E18-C4","32_Fgfr2-E18-C4",
  #Group 3
  "7_Fgfr2-Bicc1","8_Fgfr2-Bicc1","9_Fgfr2-Bicc1","17_Fgfr2-Tacc2","18_Fgfr2-Tacc2","19_Fgfr2-Tacc2",
  #Group 4
  "15_Fgfr2-dE18-Ate1","16_Fgfr2-dE18-Ate1","10_Fgfr2-dE18-Bicc1","11_Fgfr2-dE18-Bicc1","12_Fgfr2-dE18-Bicc1","20_Fgfr2-dE18-Tacc2","21_Fgfr2-dE18-Tacc2","22_Fgfr2-dE18-Tacc2",
  #Group 5
  "3_Fgfr2-E18-C2","27_Fgfr2-E18-C2","28_Fgfr2-E18-C2")] <- 2^facet.hm.plot_2[, c(
    #Group 1
    "1_Fgfr2","2_Fgfr2", "13_Fgfr2-Ate1","14_Fgfr2-Ate1",
    #Group 2
    "4_Fgfr2-dE18","5_Fgfr2-dE18","6_Fgfr2-dE18","23_Fgfr2-dE18-IGR1","24_Fgfr2-dE18-IGR1","25_Fgfr2-dE18-IGR2","26_Fgfr2-dE18-IGR2","29_Fgfr2-E18-C3","30_Fgfr2-E18-C3","31_Fgfr2-E18-C4","32_Fgfr2-E18-C4",
    #Group 3
    "7_Fgfr2-Bicc1","8_Fgfr2-Bicc1","9_Fgfr2-Bicc1","17_Fgfr2-Tacc2","18_Fgfr2-Tacc2","19_Fgfr2-Tacc2",
    #Group 4
    "15_Fgfr2-dE18-Ate1","16_Fgfr2-dE18-Ate1","10_Fgfr2-dE18-Bicc1","11_Fgfr2-dE18-Bicc1","12_Fgfr2-dE18-Bicc1","20_Fgfr2-dE18-Tacc2","21_Fgfr2-dE18-Tacc2","22_Fgfr2-dE18-Tacc2",
    #Group 5
    "3_Fgfr2-E18-C2","27_Fgfr2-E18-C2","28_Fgfr2-E18-C2")] 
glimpse(facet.hm.plot_2)
##collapse
facet.hm.plot_2_B <-facet.hm.plot_2 %>% group_by(Gene.Name_AApos.collapsed) %>% summarise_at(c(
  #Group 1
  "1_Fgfr2","2_Fgfr2", "13_Fgfr2-Ate1","14_Fgfr2-Ate1",
  #Group 2
  "4_Fgfr2-dE18","5_Fgfr2-dE18","6_Fgfr2-dE18","23_Fgfr2-dE18-IGR1","24_Fgfr2-dE18-IGR1","25_Fgfr2-dE18-IGR2","26_Fgfr2-dE18-IGR2","29_Fgfr2-E18-C3","30_Fgfr2-E18-C3","31_Fgfr2-E18-C4","32_Fgfr2-E18-C4",
  #Group 3
  "7_Fgfr2-Bicc1","8_Fgfr2-Bicc1","9_Fgfr2-Bicc1","17_Fgfr2-Tacc2","18_Fgfr2-Tacc2","19_Fgfr2-Tacc2",
  #Group 4
  "15_Fgfr2-dE18-Ate1","16_Fgfr2-dE18-Ate1","10_Fgfr2-dE18-Bicc1","11_Fgfr2-dE18-Bicc1","12_Fgfr2-dE18-Bicc1","20_Fgfr2-dE18-Tacc2","21_Fgfr2-dE18-Tacc2","22_Fgfr2-dE18-Tacc2",
  #Group 5
  "3_Fgfr2-E18-C2","27_Fgfr2-E18-C2","28_Fgfr2-E18-C2"), sum, na.rm = TRUE)
glimpse(facet.hm.plot_2_B)
facet.hm.plot_2_B[facet.hm.plot_2_B == 0] <- NA
##relog
facet.hm.plot_2_B[c(
  #Group 1
  "1_Fgfr2","2_Fgfr2", "13_Fgfr2-Ate1","14_Fgfr2-Ate1",
  #Group 2
  "4_Fgfr2-dE18","5_Fgfr2-dE18","6_Fgfr2-dE18","23_Fgfr2-dE18-IGR1","24_Fgfr2-dE18-IGR1","25_Fgfr2-dE18-IGR2","26_Fgfr2-dE18-IGR2","29_Fgfr2-E18-C3","30_Fgfr2-E18-C3","31_Fgfr2-E18-C4","32_Fgfr2-E18-C4",
  #Group 3
  "7_Fgfr2-Bicc1","8_Fgfr2-Bicc1","9_Fgfr2-Bicc1","17_Fgfr2-Tacc2","18_Fgfr2-Tacc2","19_Fgfr2-Tacc2",
  #Group 4
  "15_Fgfr2-dE18-Ate1","16_Fgfr2-dE18-Ate1","10_Fgfr2-dE18-Bicc1","11_Fgfr2-dE18-Bicc1","12_Fgfr2-dE18-Bicc1","20_Fgfr2-dE18-Tacc2","21_Fgfr2-dE18-Tacc2","22_Fgfr2-dE18-Tacc2",
  #Group 5
  "3_Fgfr2-E18-C2","27_Fgfr2-E18-C2","28_Fgfr2-E18-C2")] <- log2(facet.hm.plot_2_B[c(
    #Group 1
    "1_Fgfr2","2_Fgfr2", "13_Fgfr2-Ate1","14_Fgfr2-Ate1",
    #Group 2
    "4_Fgfr2-dE18","5_Fgfr2-dE18","6_Fgfr2-dE18","23_Fgfr2-dE18-IGR1","24_Fgfr2-dE18-IGR1","25_Fgfr2-dE18-IGR2","26_Fgfr2-dE18-IGR2","29_Fgfr2-E18-C3","30_Fgfr2-E18-C3","31_Fgfr2-E18-C4","32_Fgfr2-E18-C4",
    #Group 3
    "7_Fgfr2-Bicc1","8_Fgfr2-Bicc1","9_Fgfr2-Bicc1","17_Fgfr2-Tacc2","18_Fgfr2-Tacc2","19_Fgfr2-Tacc2",
    #Group 4
    "15_Fgfr2-dE18-Ate1","16_Fgfr2-dE18-Ate1","10_Fgfr2-dE18-Bicc1","11_Fgfr2-dE18-Bicc1","12_Fgfr2-dE18-Bicc1","20_Fgfr2-dE18-Tacc2","21_Fgfr2-dE18-Tacc2","22_Fgfr2-dE18-Tacc2",
    #Group 5
    "3_Fgfr2-E18-C2","27_Fgfr2-E18-C2","28_Fgfr2-E18-C2")] )
glimpse(facet.hm.plot_2_B)
#get KinGene column back
facet.hm.plot_2_B <- merge(facet.hm.plot_2_B, facet.hm.plot_2[, c("Gene.Name_AApos.collapsed", "KinGene")], by.x="Gene.Name_AApos.collapsed", by.y="Gene.Name_AApos.collapsed", all.x=T, all.y=F, sort=F)
glimpse(facet.hm.plot_2_B)


#make new data frame
facet.hm.plot_3 <- facet.hm.plot_2_B

### reorder data frame (collapsed sites) and calculate FC for G1+G3 vs G2
second.order.upstreamKin.via.num.subs <- facet.hm.plot_3 %>% group_by(KinGene) %>% summarize(count = n()) %>% arrange(desc(count)) %>% select(KinGene) %>% pull()
second.order.upstreamKin.via.num.subs

#calculate FC G1+G3 vs G2
glimpse(facet.hm.plot_3)
for(i in 1:nrow(facet.hm.plot_3)){
  print(i)
  leftvals <- as.numeric(facet.hm.plot_3[i,c("1_Fgfr2","2_Fgfr2", "13_Fgfr2-Ate1","14_Fgfr2-Ate1",  "7_Fgfr2-Bicc1","8_Fgfr2-Bicc1","9_Fgfr2-Bicc1","17_Fgfr2-Tacc2","18_Fgfr2-Tacc2","19_Fgfr2-Tacc2")]); leftvals #
  
  rightvals <- as.numeric(facet.hm.plot_3[i,c("4_Fgfr2-dE18","5_Fgfr2-dE18","6_Fgfr2-dE18","23_Fgfr2-dE18-IGR1","24_Fgfr2-dE18-IGR1","25_Fgfr2-dE18-IGR2","26_Fgfr2-dE18-IGR2","29_Fgfr2-E18-C3","30_Fgfr2-E18-C3","31_Fgfr2-E18-C4","32_Fgfr2-E18-C4")]); rightvals #
  
  FC <- FR_FC_LOGData_2.3.17(control = leftvals, treatment = rightvals, LOGtype = 2); FC
  
  facet.hm.plot_3$FC_FR[i] <- FC
}
glimpse(facet.hm.plot_3)


#reorder
output <- c()
for(i in second.order.upstreamKin.via.num.subs){
  print(i)
  temp <- facet.hm.plot_3 %>% filter(KinGene %in% i) %>% arrange(FC_FR)
  
  output <-bind_rows(output, temp )
}

facet.hm.plot_3 <- output
glimpse(facet.hm.plot_3)

# select columns of interest
facet.hm.plot_3 <- facet.hm.plot_3 %>% select(Gene.Name_AApos.collapsed, KinGene, 
                                              #Group 1
                                              "1_Fgfr2","2_Fgfr2", "13_Fgfr2-Ate1","14_Fgfr2-Ate1",
                                              #Group 2
                                              "4_Fgfr2-dE18","5_Fgfr2-dE18","6_Fgfr2-dE18","23_Fgfr2-dE18-IGR1","24_Fgfr2-dE18-IGR1","25_Fgfr2-dE18-IGR2","26_Fgfr2-dE18-IGR2","29_Fgfr2-E18-C3","30_Fgfr2-E18-C3","31_Fgfr2-E18-C4","32_Fgfr2-E18-C4",
                                              #Group 3
                                              "7_Fgfr2-Bicc1","8_Fgfr2-Bicc1","9_Fgfr2-Bicc1","17_Fgfr2-Tacc2","18_Fgfr2-Tacc2","19_Fgfr2-Tacc2",
                                              #Group 4
                                              "15_Fgfr2-dE18-Ate1","16_Fgfr2-dE18-Ate1","10_Fgfr2-dE18-Bicc1","11_Fgfr2-dE18-Bicc1","12_Fgfr2-dE18-Bicc1","20_Fgfr2-dE18-Tacc2","21_Fgfr2-dE18-Tacc2","22_Fgfr2-dE18-Tacc2",
                                              #Group 5
                                              "3_Fgfr2-E18-C2","27_Fgfr2-E18-C2","28_Fgfr2-E18-C2")
glimpse(facet.hm.plot_3)


###rescale: zscore 
#start define function
zscorescalebaseR<- function(x){
  #arg <- c(1, 2, 3, 4, 5)
  temp <- (x-mean(x, na.rm = T))/sd(x, na.rm = T)
  return(temp)
}
x <- c(1, 2, 3, 4, 5); x
zscorescalebaseR(x = x)
#stop define functions

facet.hm.plot_3_Z <- facet.hm.plot_3 %>% select(
  #Group 1
  "1_Fgfr2","2_Fgfr2", "13_Fgfr2-Ate1","14_Fgfr2-Ate1",
  #Group 2
  "4_Fgfr2-dE18","5_Fgfr2-dE18","6_Fgfr2-dE18","23_Fgfr2-dE18-IGR1","24_Fgfr2-dE18-IGR1","25_Fgfr2-dE18-IGR2","26_Fgfr2-dE18-IGR2","29_Fgfr2-E18-C3","30_Fgfr2-E18-C3","31_Fgfr2-E18-C4","32_Fgfr2-E18-C4",
  #Group 3
  "7_Fgfr2-Bicc1","8_Fgfr2-Bicc1","9_Fgfr2-Bicc1","17_Fgfr2-Tacc2","18_Fgfr2-Tacc2","19_Fgfr2-Tacc2",
  #Group 4
  "15_Fgfr2-dE18-Ate1","16_Fgfr2-dE18-Ate1","10_Fgfr2-dE18-Bicc1","11_Fgfr2-dE18-Bicc1","12_Fgfr2-dE18-Bicc1","20_Fgfr2-dE18-Tacc2","21_Fgfr2-dE18-Tacc2","22_Fgfr2-dE18-Tacc2",
  #Group 5
  "3_Fgfr2-E18-C2","27_Fgfr2-E18-C2","28_Fgfr2-E18-C2")
glimpse(facet.hm.plot_3_Z)


facet.hm.plot_3_Z <- t(apply(facet.hm.plot_3_Z, 1, function(x) zscorescalebaseR(x = x)))
glimpse(facet.hm.plot_3_Z)
sd(facet.hm.plot_3_Z[1,], na.rm = T) #should be 1

facet.hm.plot_3_Z <- as.data.frame(facet.hm.plot_3_Z)
glimpse(facet.hm.plot_3_Z)
facet.hm.plot_3_Z$Gene.Name_AApos.collapsed <- facet.hm.plot_3$Gene.Name_AApos.collapsed
facet.hm.plot_3_Z$KinGene <- facet.hm.plot_3$KinGene
facet.hm.plot_3_Z <- facet.hm.plot_3_Z %>% select(Gene.Name_AApos.collapsed, KinGene, 
                                                  #Group 1
                                                  "1_Fgfr2","2_Fgfr2", "13_Fgfr2-Ate1","14_Fgfr2-Ate1",
                                                  #Group 2
                                                  "4_Fgfr2-dE18","5_Fgfr2-dE18","6_Fgfr2-dE18","23_Fgfr2-dE18-IGR1","24_Fgfr2-dE18-IGR1","25_Fgfr2-dE18-IGR2","26_Fgfr2-dE18-IGR2","29_Fgfr2-E18-C3","30_Fgfr2-E18-C3","31_Fgfr2-E18-C4","32_Fgfr2-E18-C4",
                                                  #Group 3
                                                  "7_Fgfr2-Bicc1","8_Fgfr2-Bicc1","9_Fgfr2-Bicc1","17_Fgfr2-Tacc2","18_Fgfr2-Tacc2","19_Fgfr2-Tacc2",
                                                  #Group 4
                                                  "15_Fgfr2-dE18-Ate1","16_Fgfr2-dE18-Ate1","10_Fgfr2-dE18-Bicc1","11_Fgfr2-dE18-Bicc1","12_Fgfr2-dE18-Bicc1","20_Fgfr2-dE18-Tacc2","21_Fgfr2-dE18-Tacc2","22_Fgfr2-dE18-Tacc2",
                                                  #Group 5
                                                  "3_Fgfr2-E18-C2","27_Fgfr2-E18-C2","28_Fgfr2-E18-C2")
glimpse(facet.hm.plot_3_Z)



# melt data frame for plot and change column names
melted.facet.hm.plot <- reshape2::melt(facet.hm.plot_3_Z)

colnames(melted.facet.hm.plot ) <- c("PG.Genes.3", "pathway",    "variable",   "log2.int")    
glimpse(melted.facet.hm.plot)




#plot all data found
ggplot()+
  geom_tile(data=melted.facet.hm.plot, aes(y = fct_inorder(PG.Genes.3), x = variable, fill = log2.int))+ 
  scale_x_discrete(limits=c(
    #Group 1
    "1_Fgfr2","2_Fgfr2", "13_Fgfr2-Ate1","14_Fgfr2-Ate1",
    #Group 3
    "7_Fgfr2-Bicc1","8_Fgfr2-Bicc1","9_Fgfr2-Bicc1","17_Fgfr2-Tacc2","18_Fgfr2-Tacc2","19_Fgfr2-Tacc2",
    #Group 5
    "3_Fgfr2-E18-C2","27_Fgfr2-E18-C2","28_Fgfr2-E18-C2",
    #Group 2
    "4_Fgfr2-dE18","5_Fgfr2-dE18","6_Fgfr2-dE18","23_Fgfr2-dE18-IGR1","24_Fgfr2-dE18-IGR1","25_Fgfr2-dE18-IGR2","26_Fgfr2-dE18-IGR2","29_Fgfr2-E18-C3","30_Fgfr2-E18-C3","31_Fgfr2-E18-C4","32_Fgfr2-E18-C4",
    #Group 4
    "15_Fgfr2-dE18-Ate1","16_Fgfr2-dE18-Ate1","10_Fgfr2-dE18-Bicc1","11_Fgfr2-dE18-Bicc1","12_Fgfr2-dE18-Bicc1","20_Fgfr2-dE18-Tacc2","21_Fgfr2-dE18-Tacc2","22_Fgfr2-dE18-Tacc2"))+
  scale_fill_scico(palette = 'vikO',name="row z-score \n log2(int.)", na.value = "grey50")+ #
  xlab(NULL) + 
  ylab(NULL)+
  theme(strip.background=element_rect(colour="transparent", fill="transparent"))+
  theme(plot.title = element_text(hjust = 0.5, face="bold"))+  
  theme(legend.position = "right")+
  theme(axis.text.x = element_text(size = 10)) +
  theme(axis.text.y = element_text(size = 12))+
  facet_grid(fct_inorder(pathway) ~ ., scales = "free", space = "free")+
  theme(strip.text.y = element_text(angle=0, face="bold", colour="black", size=10, lineheight = 5))+ 
  theme(axis.text.x  = element_text(angle=90, hjust = 1.0, vjust = 0.5, face="bold"))+
  theme(axis.text.y = element_text(size = 7))+
  geom_vline(xintercept = c(      10.5,    13.5, 24.5) , size=1.5, linetype = "solid", color="white") #G1+3, 5, 2,4

#ggsave("tumors.IMAC.global.phospho.active.kinase.sites.all.pdf", useDingbats=FALSE,  width = 20, height = 65, units = "cm")


### visualize only selected site candidates, group/adjust names upstream kinases ################################################

## define selected sites
DZ.tum.IMAC.candidates.hm_B <- tibble(Gene.Name_AApos.collapsed=c(
  #AKT
  "Gab2_T388",              "Akt1s1_S184",            "Akt1s1_T247",            "Acly_S455",              "Bad_S112",              
  "Bad_S136","Pfkfb2_S486", "Stx7_S126",             "Tbc1d1_S231",                                       "Foxo1_S295",             
  
  #CDK / cell cycle
  "Ccdc6_S237",            
  "Cdk7_S164",              "Cdk7_T170",              "Coro1a_T418",                      
  "Rb1_S243",               "Rb1_T367",               "Runx1_S276",            "Smad3_T8",              
  "Tp53bp1_S382",           
  
  #CK2
  "Ctnnb1_S191",            "Ctnnb1_S552",            "Ctnnb1_S675",            "Ctnnd1_S252;252",       
  "Ctnnd1_S268;268",                    "Hdac2_S394",             "Hdac2_S422",            
  "Hdac2_S424",             
  
  
  #MAPK
  "Pak2_S141",              "Ptk2_S910","Ripk1_S321",             
  "Map3k2_S153",            "Map3k3_S337",            "Map2k1;Map2k2_S218;222", "Map2k1;Map2k2_S222;226", "Map2k1_T386",           
  "Map2k4_S78",             "Mapk1_T183",             "Mapk1_Y185",             "Mapk3_T203",             "Mapk3_Y205","Ksr1_T274;274","Ksr1_S392;392",            
  "Jun_S246",               "Junb_S256",              "Fosl2_S200","Eps8_S624",             "Sp1_S61",                "Sp3_S72",               
  "Gsk3b_S9",               "Map1b_T1784",            "Map2_S1783",             
  
  #mTOR
  "Rps6ka1_S352",           "Rps6ka3_S386",          
  "Rps6_S236",              "Rps6_S240",              "Rps6_S244",              "Eif4b_S422",             "Eif4ebp1_T36",          
  "Eif4ebp1_S64","Crtc2_S70",                      "Igf2bp1_S181",           "Pdcd4_S457", "Trim28_S473", "Usp20_S132"),
  
  Group.DZ =c( rep("AKT", 10), rep("CDK / cell cycle", 9), rep("CK2", 8), rep("MAPK", 24), rep("mTOR", 13)   ),
  order.FR = c(1:64),
  order.FR.rev = c(64:1)
)
print(DZ.tum.IMAC.candidates.hm_B, n=64)

#filter for selected candidates
glimpse(facet.hm.plot_3_Z)

facet.hm.plot_filter.custom.DZ <- facet.hm.plot_3_Z %>% filter(Gene.Name_AApos.collapsed %in% c(DZ.tum.IMAC.candidates.hm_B$Gene.Name_AApos.collapsed, preserve = TRUE)) %>% unique()
glimpse(facet.hm.plot_filter.custom.DZ)
table(duplicated(facet.hm.plot_filter.custom.DZ$Gene.Name_AApos.collapsed))

#merge order and adjust names upstream kinases arrange to custom order
facet.hm.plot_filter.custom.DZ <- merge(facet.hm.plot_filter.custom.DZ, DZ.tum.IMAC.candidates.hm_B, by.x = "Gene.Name_AApos.collapsed", by.y = "Gene.Name_AApos.collapsed", all.x=T, all.y = T, sort=FALSE)
facet.hm.plot_filter.custom.DZ <- facet.hm.plot_filter.custom.DZ  %>% arrange(match(Group.DZ, c("AKT", "CDK / cell cycle", "CK2", "MAPK", "mTOR")), order.FR.rev) %>% unique()
glimpse(facet.hm.plot_filter.custom.DZ)

#add alternative name
facet.hm.plot_filter.custom.DZ$Gene.Name <- pbsapply(facet.hm.plot_filter.custom.DZ$Gene.Name_AApos.collapsed, function(x) unlist(str_split(string=x, pattern="_"))[1])


## reshape for plot
melted.facet.hm.plot_filter.custom.DZ <- reshape2::melt(facet.hm.plot_filter.custom.DZ  %>% select(Gene.Name_AApos.collapsed, Group.DZ, 
                                                                                                   #Group 1
                                                                                                   "1_Fgfr2","2_Fgfr2", "13_Fgfr2-Ate1","14_Fgfr2-Ate1",
                                                                                                   #Group 2
                                                                                                   "4_Fgfr2-dE18","5_Fgfr2-dE18","6_Fgfr2-dE18","23_Fgfr2-dE18-IGR1","24_Fgfr2-dE18-IGR1","25_Fgfr2-dE18-IGR2","26_Fgfr2-dE18-IGR2","29_Fgfr2-E18-C3","30_Fgfr2-E18-C3","31_Fgfr2-E18-C4","32_Fgfr2-E18-C4",
                                                                                                   #Group 3
                                                                                                   "7_Fgfr2-Bicc1","8_Fgfr2-Bicc1","9_Fgfr2-Bicc1","17_Fgfr2-Tacc2","18_Fgfr2-Tacc2","19_Fgfr2-Tacc2",
                                                                                                   #Group 4
                                                                                                   "15_Fgfr2-dE18-Ate1","16_Fgfr2-dE18-Ate1","10_Fgfr2-dE18-Bicc1","11_Fgfr2-dE18-Bicc1","12_Fgfr2-dE18-Bicc1","20_Fgfr2-dE18-Tacc2","21_Fgfr2-dE18-Tacc2","22_Fgfr2-dE18-Tacc2",
                                                                                                   #Group 5
                                                                                                   "3_Fgfr2-E18-C2","27_Fgfr2-E18-C2","28_Fgfr2-E18-C2"))
glimpse(melted.facet.hm.plot_filter.custom.DZ)

#add Gene.Name 
melted.facet.hm.plot_filter.custom.DZ$Gene.Name <- pbsapply(melted.facet.hm.plot_filter.custom.DZ$Gene.Name_AApos.collapsed, function(x) unlist(str_split(string=x, pattern="_"))[1])
glimpse(melted.facet.hm.plot_filter.custom.DZ)


### facet plot selected candidates DZ
ggplot()+
  geom_tile(data=melted.facet.hm.plot_filter.custom.DZ, aes(y = fct_inorder(Gene.Name_AApos.collapsed), x = variable, fill = value))+
  scale_x_discrete(limits=c(
    #Group 1
    "1_Fgfr2","2_Fgfr2", "13_Fgfr2-Ate1","14_Fgfr2-Ate1",
    #Group 3
    "7_Fgfr2-Bicc1","8_Fgfr2-Bicc1","9_Fgfr2-Bicc1","17_Fgfr2-Tacc2","18_Fgfr2-Tacc2","19_Fgfr2-Tacc2",
    #Group 5
    "3_Fgfr2-E18-C2","27_Fgfr2-E18-C2","28_Fgfr2-E18-C2",
    #Group 2
    "4_Fgfr2-dE18","5_Fgfr2-dE18","6_Fgfr2-dE18","23_Fgfr2-dE18-IGR1","24_Fgfr2-dE18-IGR1","25_Fgfr2-dE18-IGR2","26_Fgfr2-dE18-IGR2","29_Fgfr2-E18-C3","30_Fgfr2-E18-C3","31_Fgfr2-E18-C4","32_Fgfr2-E18-C4",
    #Group 4
    "15_Fgfr2-dE18-Ate1","16_Fgfr2-dE18-Ate1","10_Fgfr2-dE18-Bicc1","11_Fgfr2-dE18-Bicc1","12_Fgfr2-dE18-Bicc1","20_Fgfr2-dE18-Tacc2","21_Fgfr2-dE18-Tacc2","22_Fgfr2-dE18-Tacc2"))+
  scale_fill_scico(palette = 'vikO',name="row z-score \n log2(int.)", na.value = "grey50")+ #
  xlab(NULL) + 
  ylab(NULL)+
  theme(strip.background=element_rect(colour="transparent", fill="transparent"))+ 
  theme(plot.title = element_text(hjust = 0.5, face="bold"))+ 
  theme(legend.position = "right")+
  theme(axis.text.x = element_text(size = 10)) +
  theme(axis.text.y = element_text(size = 12))+
  facet_grid(fct_inorder(Group.DZ) ~ ., scales = "free", space = "free")+
  theme(strip.text.y = element_text(angle=0, face="bold", colour="black", size=10, lineheight = 5))+
  theme(axis.text.x  = element_text(angle=90, hjust = 1.0, vjust = 0.5, face="bold"))+
  theme(axis.text.y = element_text(size = 7))+
  geom_vline(xintercept = c(10.5, 13.5, 24.5) , size=1.5, linetype = "solid", color="white") #G1+3,5, 2,4

#ggsave("tumors.IMAC.global.phospho.custom.candidates.pdf", useDingbats=FALSE,  width =17, height =19, units = "cm")


#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
### tumors IMAC: seperate two group comparisons for PTMSEA

##see two group comparisons above
#G1+G3 vs. G2
#G1+G3 vs. G4
#G2 vs G4

#################################################################################################################################
#G1+G3 vs G2 tumors IMAC Limma seperate 2 group comparison for classic PTMSEA ##############################################################################

glimpse(G1and3_G2_Limma.2gc.XYfilter.2.B)

##G1and3_G2_Limma.2gc.XYfilter.2.B reshape for classic PTMSEA
#mouse windows
G1and3_G2_PTMSEA <- G1and3_G2_Limma.2gc.XYfilter.2.B
glimpse(G1and3_G2_PTMSEA)

# make aa15.window
# take first entry sequence window if there are two entries
G1and3_G2_PTMSEA$aa15.window <- pbsapply(G1and3_G2_PTMSEA$Sequence.window, function(x) unlist(str_split(string =x, pattern = ";"))[1]  )
G1and3_G2_PTMSEA$aa15.window <- pbsapply(G1and3_G2_PTMSEA$aa15.window, function(x) make.aa.15.window.from.31.sequence.window(input = x)  )
glimpse(G1and3_G2_PTMSEA)

table(duplicated(G1and3_G2_PTMSEA$aa15.window))
which(duplicated(G1and3_G2_PTMSEA$aa15.window))


#remove duplicated amino acid windows 
##if there are duplicate entries (e.g. because of _1_2_3) take the entry with the best p-value
unique.15aa.windows <- unique(G1and3_G2_PTMSEA$aa15.window)
glimpse(unique.15aa.windows) #6657

G1and3_G2_PTMSEA.2 <- c()
for (i in unique.15aa.windows){
  print(i)
  temp <- G1and3_G2_PTMSEA %>% filter(aa15.window %in% c(i))
  temp <- temp %>% arrange(Limma.p.value)
  temp <- temp[1,]
  G1and3_G2_PTMSEA.2 <- bind_rows(G1and3_G2_PTMSEA.2, temp)
}

glimpse(G1and3_G2_PTMSEA.2)
table(duplicated(G1and3_G2_PTMSEA.2$aa15.window))

#order according to signed , (log-transformed) pvalue
G1and3_G2_PTMSEA.2$signFC_FR <- sign(G1and3_G2_PTMSEA.2$FC_FR)
G1and3_G2_PTMSEA.2$neg.log10.p <- -log10(G1and3_G2_PTMSEA.2$Limma.p.value)
G1and3_G2_PTMSEA.2$rank.metric <- G1and3_G2_PTMSEA.2$neg.log10.p * G1and3_G2_PTMSEA.2$signFC_FR
glimpse(G1and3_G2_PTMSEA.2) 
table(is.na(G1and3_G2_PTMSEA.2$rank.metric))

G1and3_G2_PTMSEA.2 <- G1and3_G2_PTMSEA.2 %>% arrange(desc(rank.metric))
head(G1and3_G2_PTMSEA.2$FC_FR)
tail(G1and3_G2_PTMSEA.2$FC_FR)

# add_p to end of 15aa window
G1and3_G2_PTMSEA.2$aa15.window.2 <- paste0(G1and3_G2_PTMSEA.2$aa15.window, "-p")
glimpse(G1and3_G2_PTMSEA.2) 


#save data
save(G1and3_G2_PTMSEA.2, file = "G1and3_G2.classic.PTMSEA.Rdata")

#prepare input for gene pattern/PTMSEA 
gene.pattern.PTM.SEA.input <- G1and3_G2_PTMSEA.2 %>% select(aa15.window.2, 
                                                            #Group 1
                                                            "1_Fgfr2","2_Fgfr2", "13_Fgfr2-Ate1","14_Fgfr2-Ate1",
                                                            #Group 3
                                                            "7_Fgfr2-Bicc1","8_Fgfr2-Bicc1","9_Fgfr2-Bicc1","17_Fgfr2-Tacc2","18_Fgfr2-Tacc2","19_Fgfr2-Tacc2",
                                                            
                                                            #Group 2
                                                            "4_Fgfr2-dE18","5_Fgfr2-dE18","6_Fgfr2-dE18","23_Fgfr2-dE18-IGR1","24_Fgfr2-dE18-IGR1","25_Fgfr2-dE18-IGR2","26_Fgfr2-dE18-IGR2","29_Fgfr2-E18-C3","30_Fgfr2-E18-C3","31_Fgfr2-E18-C4","32_Fgfr2-E18-C4",
                                                            
                                                            "rank.metric")
glimpse(gene.pattern.PTM.SEA.input)
nrow(gene.pattern.PTM.SEA.input) #6657

nrow(filter(gene.pattern.PTM.SEA.input, !is.na(rank.metric)))#6657


#make GCT object with cmapR
gene.pattern.PTM.SEA.input_GCT <- as.matrix(gene.pattern.PTM.SEA.input[,"rank.metric"])
glimpse(gene.pattern.PTM.SEA.input_GCT)
colnames(gene.pattern.PTM.SEA.input_GCT) <- "rank.metric"
rownames(gene.pattern.PTM.SEA.input_GCT) <- gene.pattern.PTM.SEA.input$aa15.window.2
glimpse(gene.pattern.PTM.SEA.input_GCT)

gene.pattern.PTM.SEA.input_GCT.2 <- new("GCT", mat=gene.pattern.PTM.SEA.input_GCT)
glimpse(gene.pattern.PTM.SEA.input_GCT.2)

#save GCT file for PTMSEA input
write_gct(gene.pattern.PTM.SEA.input_GCT.2, "DZ.tumors.IMAC.mouse.windows.G1and3_G2.PTMSEA.input.gct")



#################################################################################################################################
#G1+G3 vs G4 tumors IMAC Limma seperate 2 group comparison for classic PTMSEA ##############################################################################

glimpse(G1andG3_G4_Limma.2gc.XYfilter.2.B)

##G1andG3_G4_Limma.2gc.XYfilter.2.B reshape for classic PTMSEA
#mouse windows
G1andG3_G4_PTMSEA <- G1andG3_G4_Limma.2gc.XYfilter.2.B
glimpse(G1andG3_G4_PTMSEA)

# make aa15.window
# take first entry sequence window if there are two entries
G1andG3_G4_PTMSEA$aa15.window <- pbsapply(G1andG3_G4_PTMSEA$Sequence.window, function(x) unlist(str_split(string =x, pattern = ";"))[1]  )
G1andG3_G4_PTMSEA$aa15.window <- pbsapply(G1andG3_G4_PTMSEA$aa15.window, function(x) make.aa.15.window.from.31.sequence.window(input = x)  )
glimpse(G1andG3_G4_PTMSEA)

table(duplicated(G1andG3_G4_PTMSEA$aa15.window))
which(duplicated(G1andG3_G4_PTMSEA$aa15.window))


#remove duplicated amino acid windows 
##if there are duplicate entries (e.g. because of _1_2_3) take the entry with the best p-value
unique.15aa.windows <- unique(G1andG3_G4_PTMSEA$aa15.window)
glimpse(unique.15aa.windows) #6082


G1andG3_G4_PTMSEA.2 <- c()
for (i in unique.15aa.windows){
  print(i)
  temp <- G1andG3_G4_PTMSEA %>% filter(aa15.window %in% c(i))
  temp <- temp %>% arrange(Limma.p.value)
  temp <- temp[1,]
  G1andG3_G4_PTMSEA.2 <- bind_rows(G1andG3_G4_PTMSEA.2, temp)
}

glimpse(G1andG3_G4_PTMSEA.2)
table(duplicated(G1andG3_G4_PTMSEA.2$aa15.window))

#order according to signed , (log-transformed) pvalue
G1andG3_G4_PTMSEA.2$signFC_FR <- sign(G1andG3_G4_PTMSEA.2$FC_FR)
G1andG3_G4_PTMSEA.2$neg.log10.p <- -log10(G1andG3_G4_PTMSEA.2$Limma.p.value)
G1andG3_G4_PTMSEA.2$rank.metric <- G1andG3_G4_PTMSEA.2$neg.log10.p * G1andG3_G4_PTMSEA.2$signFC_FR
glimpse(G1andG3_G4_PTMSEA.2) 
table(is.na(G1andG3_G4_PTMSEA.2$rank.metric))

G1andG3_G4_PTMSEA.2 <- G1andG3_G4_PTMSEA.2 %>% arrange(desc(rank.metric))
head(G1andG3_G4_PTMSEA.2$FC_FR)
tail(G1andG3_G4_PTMSEA.2$FC_FR)

# add_p to end of 15aa window
G1andG3_G4_PTMSEA.2$aa15.window.2 <- paste0(G1andG3_G4_PTMSEA.2$aa15.window, "-p")
glimpse(G1andG3_G4_PTMSEA.2) 


# sava data
save(G1andG3_G4_PTMSEA.2, file = "G1andG3_G4.classic.PTMSEA.Rdata")



#prepare input for gene pattern/PTMSEA
gene.pattern.PTM.SEA.input <- G1andG3_G4_PTMSEA.2 %>% select(aa15.window.2, 
                                                             #Group 1
                                                             "1_Fgfr2","2_Fgfr2", "13_Fgfr2-Ate1","14_Fgfr2-Ate1",
                                                             #Group 3
                                                             "7_Fgfr2-Bicc1","8_Fgfr2-Bicc1","9_Fgfr2-Bicc1","17_Fgfr2-Tacc2","18_Fgfr2-Tacc2","19_Fgfr2-Tacc2",
                                                             
                                                             
                                                             #Group 4
                                                             "15_Fgfr2-dE18-Ate1","16_Fgfr2-dE18-Ate1","10_Fgfr2-dE18-Bicc1","11_Fgfr2-dE18-Bicc1","12_Fgfr2-dE18-Bicc1","20_Fgfr2-dE18-Tacc2","21_Fgfr2-dE18-Tacc2","22_Fgfr2-dE18-Tacc2",
                                                             
                                                             "rank.metric")
glimpse(gene.pattern.PTM.SEA.input)
nrow(gene.pattern.PTM.SEA.input) #6082

#make GCT object with cmapR
gene.pattern.PTM.SEA.input_GCT <- as.matrix(gene.pattern.PTM.SEA.input[,"rank.metric"])
glimpse(gene.pattern.PTM.SEA.input_GCT)
colnames(gene.pattern.PTM.SEA.input_GCT) <- "rank.metric"
rownames(gene.pattern.PTM.SEA.input_GCT) <- gene.pattern.PTM.SEA.input$aa15.window.2
glimpse(gene.pattern.PTM.SEA.input_GCT)

gene.pattern.PTM.SEA.input_GCT.2 <- new("GCT", mat=gene.pattern.PTM.SEA.input_GCT)
glimpse(gene.pattern.PTM.SEA.input_GCT.2)

#save GCT file for PTMSEA input
write_gct(gene.pattern.PTM.SEA.input_GCT.2, "DZ.tumors.IMAC.mouse.windows.G1andG3_G4.PTMSEA.input.gct")




#####################################################################################################################################
#####################################################################################################################################
###tumors IMAC PTMSEA classic Group 2 vs Group 4

glimpse(G4_G2_Limma.2gc.XYfilter.2.B) #7,622

##reshape for PTMSEA

#mouse windows
G2_G4_PTMSEA <- G4_G2_Limma.2gc.XYfilter.2.B
glimpse(G2_G4_PTMSEA)


# make aa15.window
# take first entry sequence window if there are two entries
G2_G4_PTMSEA$aa15.window <- pbsapply(G2_G4_PTMSEA$Sequence.window, function(x) unlist(str_split(string =x, pattern = ";"))[1]  )
G2_G4_PTMSEA$aa15.window <- pbsapply(G2_G4_PTMSEA$aa15.window, function(x) make.aa.15.window.from.31.sequence.window(input = x)  )
glimpse(G2_G4_PTMSEA)


table(duplicated(G2_G4_PTMSEA$aa15.window))
which(duplicated(G2_G4_PTMSEA$aa15.window))


#remove duplicated amino acid windows 
##if there are duplicate entries (e.g. because of _1_2_3) take the entry with the best p-value
unique.15aa.windows <- unique(G2_G4_PTMSEA$aa15.window)
glimpse(unique.15aa.windows) #6864


G2_G4_PTMSEA.2 <- c()
for (i in unique.15aa.windows){
  print(i)
  temp <- G2_G4_PTMSEA %>% filter(aa15.window %in% c(i))
  temp <- temp %>% arrange(Limma.p.value)
  temp <- temp[1,]
  G2_G4_PTMSEA.2 <- bind_rows(G2_G4_PTMSEA.2, temp)
}

glimpse(G2_G4_PTMSEA.2)
table(duplicated(G2_G4_PTMSEA.2$aa15.window))

#order according to signed , (log-transformed) pvalue
G2_G4_PTMSEA.2$signFC_FR <- sign(G2_G4_PTMSEA.2$FC_FR)
G2_G4_PTMSEA.2$neg.log10.p <- -log10(G2_G4_PTMSEA.2$Limma.p.value)
G2_G4_PTMSEA.2$rank.metric <- G2_G4_PTMSEA.2$neg.log10.p * G2_G4_PTMSEA.2$signFC_FR
glimpse(G2_G4_PTMSEA.2) 
table(is.na(G2_G4_PTMSEA.2$rank.metric))

G2_G4_PTMSEA.2 <- G2_G4_PTMSEA.2 %>% arrange(desc(rank.metric))
head(G2_G4_PTMSEA.2$FC_FR)
tail(G2_G4_PTMSEA.2$FC_FR)

# add_p to end of 15aa window
G2_G4_PTMSEA.2$aa15.window.2 <- paste0(G2_G4_PTMSEA.2$aa15.window, "-p")
glimpse(G2_G4_PTMSEA.2) 


### save data
save(G2_G4_PTMSEA.2, file = "G2_G4.classic.PTMSEA.Rdata")

#prepare gene pattern/PTMSEA input
gene.pattern.PTM.SEA.input <- G2_G4_PTMSEA.2 %>% select(aa15.window.2, 
                                                        #Group 2
                                                        "4_Fgfr2-dE18","5_Fgfr2-dE18","6_Fgfr2-dE18","23_Fgfr2-dE18-IGR1","24_Fgfr2-dE18-IGR1","25_Fgfr2-dE18-IGR2","26_Fgfr2-dE18-IGR2","29_Fgfr2-E18-C3","30_Fgfr2-E18-C3","31_Fgfr2-E18-C4","32_Fgfr2-E18-C4",
                                                        
                                                        #Group 4
                                                        "15_Fgfr2-dE18-Ate1","16_Fgfr2-dE18-Ate1","10_Fgfr2-dE18-Bicc1","11_Fgfr2-dE18-Bicc1","12_Fgfr2-dE18-Bicc1","20_Fgfr2-dE18-Tacc2","21_Fgfr2-dE18-Tacc2","22_Fgfr2-dE18-Tacc2",
                                                        
                                                        "rank.metric")
glimpse(gene.pattern.PTM.SEA.input)
nrow(gene.pattern.PTM.SEA.input) #6864


#make GCT object with cmapR
gene.pattern.PTM.SEA.input_GCT <- as.matrix(gene.pattern.PTM.SEA.input[,"rank.metric"])
glimpse(gene.pattern.PTM.SEA.input_GCT)
colnames(gene.pattern.PTM.SEA.input_GCT) <- "rank.metric"
rownames(gene.pattern.PTM.SEA.input_GCT) <- gene.pattern.PTM.SEA.input$aa15.window.2
glimpse(gene.pattern.PTM.SEA.input_GCT)

gene.pattern.PTM.SEA.input_GCT.2 <- new("GCT", mat=gene.pattern.PTM.SEA.input_GCT)
glimpse(gene.pattern.PTM.SEA.input_GCT.2)

#save GCT file for PTMSEA input
write_gct(gene.pattern.PTM.SEA.input_GCT.2, "DZ.tumors.IMAC.mouse.windows.G2_G4.PTMSEA.input.gct")


#####################################################################################################################################################################
#####################################################################################################################################################################
### tumor IMAC: analyze and result classic PTM-SEA obtained via gene pattern

## load corresponding data an plot, repeat with remaining results, change plot titles etc. accordingly

##mouse windows 
##DZ tumors IMAC G1+G3 vs G2 classic PTMSEA  
#dz.tum..PTMSEA <- fread("DZ.tumors.mouse.G1and3_G2.PTMSEA-combined.txt")
#nrow(dz.tum..PTMSEA) #30
#glimpse(dz.tum..PTMSEA)

##mouse windows 
##DZ tumors IMAC G1+G3 vs G4 classic PTMSEA  
#dz.tum..PTMSEA <- fread("/Users/frankrolfs/NLPostDR/VUmc OPL/Phosphoproteomics with Daniel Zingg NKI/mouse tumors/IMAC.mouse.tumors.DZ/DZ.tumors.IMAC.classicPTMSEA/G1andG3_G4/403670/11.01.22.DZ.tumors.IMAC.mouse.windows.G1andG3_G4.PTMSEA.input-combined.txt")
#nrow(dz.tum..PTMSEA)
#glimpse(dz.tum..PTMSEA)


##mouse windows 
##DZ tumors IMAC G2 vs G4 classic PTMSEA  
dz.tum..PTMSEA <- fread("/Users/frankrolfs/NLPostDR/VUmc OPL/Phosphoproteomics with Daniel Zingg NKI/mouse tumors/IMAC.mouse.tumors.DZ/DZ.tumors.IMAC.classicPTMSEA/G2_G4/403614/11.01.22.DZ.tumors.IMAC.mouse.G2_G4.PTMSEA-combined.txt")
nrow(dz.tum..PTMSEA)
glimpse(dz.tum..PTMSEA)


# reshape data
data.result.ptmsea <- dz.tum..PTMSEA 
glimpse(data.result.ptmsea)


#calculate  -log10 p-value
data.result.ptmsea$neg.log.pval <- -log10(data.result.ptmsea$pvalue.rank.metric)
glimpse(data.result.ptmsea)


###LOLLIPOP Plot PTMSEA result  
glimpse(data.result.ptmsea)
data.result.ptmsea$fdr.pvalue <- pbsapply(data.result.ptmsea$fdr.pvalue.rank.metric, function(x)x<0.05) #
data.result.ptmsea$pvalue.below.0.05 <- pbsapply(data.result.ptmsea$pvalue.rank.metric, function(x)x<0.05) #

ggplot() +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  geom_hline(yintercept=0, linetype="solid", color="black", size=0.5)+
  geom_segment(data=data.result.ptmsea,  aes(x=reorder(id, rank.metric), xend=reorder(id, rank.metric), y=0, yend=rank.metric)) +
  geom_point(data=data.result.ptmsea, aes(x=reorder(id, rank.metric),  y=rank.metric, shape=pvalue.below.0.05),size=3.5)+ # 
  scale_shape_manual(values=c(17,19),name="p-value < 0.05")+
  theme(legend.position="right")+
  labs(x="", y="Normalized Enrichment Score")+
  coord_flip() +
  theme(axis.text.x = element_text(size = 16))+
  theme(axis.text.y = element_text(size = 16))+ 
 
  #ggtitle("tumors IMAC \n G1andG3 vs. G2 PTM-SEA \n mouse windows")+
  #ggtitle("tumors IMAC \n G1+G3 vs. G4 PTM-SEA \n mouse windows")+
  ggtitle("tumors IMAC \n G2 vs. G4 PTM-SEA \n mouse windows")+
 
  theme(legend.position = "right")+
  theme(axis.text.x = element_text(size = 14)) +
  theme(axis.text.y = element_text(size = 14)) +
  theme(axis.ticks.y=element_blank() ) +
  theme(axis.line.y =element_blank() ) +
  theme(plot.title = element_text(hjust = 0.5,face = "bold", size = 16)) 

#ggsave("Fgfr2.tumors.G1andG3.vs.G2.PTMSEA.mouse.version.no.filters.pdf", useDingbats=FALSE,  width = 20, height = 15, units = "cm") 

#ggsave("Fgfr2.tumors.G1andG3.vs.G4.PTMSEA.mouse.version.no.filters.pdf", useDingbats=FALSE,  width = 20, height = 15, units = "cm") 

#ggsave("Fgfr2.tumors.G2.vs.G4.PTMSEA.mouse.version.no.filters.pdf", useDingbats=FALSE,  width = 20, height = 15, units = "cm") 






















