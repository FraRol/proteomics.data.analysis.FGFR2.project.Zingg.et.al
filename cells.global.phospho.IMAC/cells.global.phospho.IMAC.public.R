#_______________________________________________________________________________________________________________________
# 17.01.2022
# 
# Project: proteomics.truncated.FGFR2.is.oncogene.cancer.Zingg.et.al 
# 
# Script name: cells.global.phospho.IMAC.public
#
# Purpose of script: analysis protein phosphorylation data cells truncated FGFR2 is oncogene cancer project Zingg et al 
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
### load data cells IMAC DDA from MaxQuant Phospho(STY)Sites.txt

### sample overview
# A	1	  GFP
# F	2	  GFP
# K	3	  GFP
# P	4	  GFP

# B	5	  Fgfr2-FL
# G	6	  Fgfr2-FL
# L	7	  Fgfr2-FL
# Q	8	  Fgfr2-FL

# C	9	  Fgfr2-dE18
# H	10	Fgfr2-dE18
# M	11	Fgfr2-dE18
# R	12	Fgfr2-dE18

# D	13	Fgfr2-FL-Bicc1
# I	14	Fgfr2-FL-Bicc1
# N	15	Fgfr2-FL-Bicc1
# S	16	Fgfr2-FL-Bicc1

# E	17	Fgfr2-dE18-Bicc1
# J	18	Fgfr2-dE18-Bicc1
# O	19	Fgfr2-dE18-Bicc1
# T	20	Fgfr2-dE18-Bicc1



load("cells.IMAC_Phospho(STY)Sites.txt.Rdata")

glimpse(DZ.IMAC.MBR.ON)

nrow(DZ.IMAC.MBR.ON) #20871
ncol(DZ.IMAC.MBR.ON) #311
colnames(DZ.IMAC.MBR.ON)
length(unique(DZ.IMAC.MBR.ON$id)) #20871


### select columns of interest
DZ.IMAC.MBR.ON.2 <- DZ.IMAC.MBR.ON %>% select(
  
  "id",                          "Proteins",                    "Positions within proteins",   "Leading proteins",            "Protein",                    
  "Protein names",               "Gene names",                  "Fasta headers",               "Localization prob",           "Number of Phospho (STY)",    
  "Amino acid",                  "Sequence window",             "Phospho (STY) Probabilities", "Position in peptide",         "Reverse",                    
  "Potential contaminant",       "Positions",                   "Position",                    "Peptide IDs",                 "Mod. peptide IDs",      
  
  matches("___[0-9]"), -Intensity___1, -Intensity___2, -Intensity___3
  
)

glimpse(DZ.IMAC.MBR.ON.2)


### remove spaces from column names
column.names <- colnames(DZ.IMAC.MBR.ON.2)
column.names.changed <- unname(sapply(column.names , function(x) str_replace_all(string=x, pattern=" ", replacement=".")))
colnames(DZ.IMAC.MBR.ON.2) <- column.names.changed

glimpse(DZ.IMAC.MBR.ON.2)
colnames(DZ.IMAC.MBR.ON.2)

#change column classes
DZ.IMAC.MBR.ON.2 <- mutate_at(DZ.IMAC.MBR.ON.2,  21:ncol(DZ.IMAC.MBR.ON.2), list(as.numeric) )

glimpse(DZ.IMAC.MBR.ON.2)

### remove reverse hits
table(DZ.IMAC.MBR.ON.2$Reverse) #205+
DZ.IMAC.MBR.ON.3 <- DZ.IMAC.MBR.ON.2[DZ.IMAC.MBR.ON.2$Reverse =="",]
nrow(DZ.IMAC.MBR.ON.3) #20666

###remove contaminants
table(DZ.IMAC.MBR.ON.3$Potential.contaminant) #59+
DZ.IMAC.MBR.ON.3 <- DZ.IMAC.MBR.ON.3[DZ.IMAC.MBR.ON.3$Potential.contaminant == "",]
nrow(DZ.IMAC.MBR.ON.3) #20607 

### remove rows completely zero
glimpse(DZ.IMAC.MBR.ON.3)
colnames(DZ.IMAC.MBR.ON.3)
ncol(DZ.IMAC.MBR.ON.3) #80
DZ.IMAC.MBR.ON.3$zero.check <- pbapply(DZ.IMAC.MBR.ON.3[, c(21:80)], 1, function(x) sum(x) ) #
table(DZ.IMAC.MBR.ON.3$zero.check >0) #F3459 T17148

DZ.IMAC.MBR.ON.3.complete.zero.rows <- filter(DZ.IMAC.MBR.ON.3, zero.check <= 0)
glimpse(DZ.IMAC.MBR.ON.3.complete.zero.rows) #

DZ.IMAC.MBR.ON.3 <- filter(DZ.IMAC.MBR.ON.3, zero.check > 0)
nrow(DZ.IMAC.MBR.ON.3) #17148
table(DZ.IMAC.MBR.ON.3$zero.check >0)


length(unique(DZ.IMAC.MBR.ON.3$id))#17148 = total number PS
table(DZ.IMAC.MBR.ON.3$Amino.acid) #S14055 T2481 Y612
length(unique(DZ.IMAC.MBR.ON.3$Sequence.window)) #17135

### remove column zero.check
DZ.IMAC.MBR.ON.3 <- DZ.IMAC.MBR.ON.3 %>% select(-zero.check)
glimpse(DZ.IMAC.MBR.ON.3)


### create variable with short sample names with order as in Phospho(STY)Sites.txt
samples.in.experiment <- colnames(DZ.IMAC.MBR.ON.3)
samples.in.experiment <- samples.in.experiment[21:ncol(DZ.IMAC.MBR.ON.3)]
samples.in.experiment <- str_remove(string=samples.in.experiment, pattern="Intensity.")
samples.in.experiment <- str_remove(string=samples.in.experiment, pattern="___[0-9]")
samples.in.experiment <- unique(samples.in.experiment)
samples.in.experiment



### loop over data to count number of PS per sample ############################################################################
colnames(DZ.IMAC.MBR.ON.3)
ncol(DZ.IMAC.MBR.ON.3) #80

result.PS.count.DZ.IMAC.MBR.ON <- c()
sample.name.count <- 1

for(i.1 in seq(from=21, to=80, by=3)){ 
  print(i.1)
  
  i.2 <- i.1+1
  i.3 <- i.2+1
  
  sample.name <- samples.in.experiment[sample.name.count]
  print(sample.name)
  col.names <- colnames(DZ.IMAC.MBR.ON.3)
  col.names.select <- col.names[c(i.1, i.2, i.3)]
  
  temp <- select(DZ.IMAC.MBR.ON.3,Amino.acid, col.names.select)
  temp$sum_1_2_3 <- pbapply(temp[, 2:4], 1, function(x) sum(x))
  
  numberPS.A <- filter(temp, sum_1_2_3 > 0)
  numberPS <- nrow(numberPS.A); numberPS
  
  nMVs <- filter(temp, sum_1_2_3 <= 0)
  nMVs <- nrow(nMVs); nMVs
  
  percentageSTY <- table(numberPS.A$Amino.acid); percentageSTY
  
  temp.result <- tibble(sample=sample.name, numberPS=numberPS, nMVs=nMVs, nS=percentageSTY[1], nT=percentageSTY[2], nY=percentageSTY[3])
  
  result.PS.count.DZ.IMAC.MBR.ON  <- bind_rows(result.PS.count.DZ.IMAC.MBR.ON, temp.result)
  
  sample.name.count <- sample.name.count +1
  
}
result.PS.count.DZ.IMAC.MBR.ON

mean(result.PS.count.DZ.IMAC.MBR.ON$numberPS) #8849.8

### sample overview
# A	1	GFP
# F	2	GFP
# K	3	GFP
# P	4	GFP

# B	5	Fgfr2-FL
# G	6	Fgfr2-FL
# L	7	Fgfr2-FL
# Q	8	Fgfr2-FL

# C	9	Fgfr2-dE18
# H	10	Fgfr2-dE18
# M	11	Fgfr2-dE18
# R	12	Fgfr2-dE18

# D	13	Fgfr2-FL-Bicc1
# I	14	Fgfr2-FL-Bicc1
# N	15	Fgfr2-FL-Bicc1
# S	16	Fgfr2-FL-Bicc1

# E	17	Fgfr2-dE18-Bicc1
# J	18	Fgfr2-dE18-Bicc1
# O	19	Fgfr2-dE18-Bicc1
# T	20	Fgfr2-dE18-Bicc1


DZ.IMAC.MBR.ON.samples.order <- c(
  "A", "F", "K", "P", #GFP
  
  
  "B", "G", "L", "Q",  #FL
  
  "C", "H", "M", "R",  #dE18
  
  "D", "I", "N", "S",  #FL-Bicc1
  
  "E", "J", "O", "T"   #dE18-Bicc1
  
  
)
DZ.IMAC.MBR.ON.samples.order

# calculate percentage MV PS all classes level
nrow(DZ.IMAC.MBR.ON.3)  #17148 = total number PS

result.PS.count.DZ.IMAC.MBR.ON$percentageMV <- pbsapply(result.PS.count.DZ.IMAC.MBR.ON$nMVs, function(x) (100*x)/17148) #
result.PS.count.DZ.IMAC.MBR.ON

# calculate percentage pS/pT/pY PS all classes level
for(i in 1:nrow(result.PS.count.DZ.IMAC.MBR.ON)){
  print(i)
  
  total.PS.sample <- result.PS.count.DZ.IMAC.MBR.ON$numberPS[i]
  n.pS <- result.PS.count.DZ.IMAC.MBR.ON$nS[i]
  n.pT <- result.PS.count.DZ.IMAC.MBR.ON$nT[i]
  n.pY <- result.PS.count.DZ.IMAC.MBR.ON$nY[i]
  
  result.PS.count.DZ.IMAC.MBR.ON$percentage.pS[i] <- (100*n.pS)/total.PS.sample
  result.PS.count.DZ.IMAC.MBR.ON$percentage.pT[i] <- (100*n.pT)/total.PS.sample
  result.PS.count.DZ.IMAC.MBR.ON$percentage.pY[i] <- (100*n.pY)/total.PS.sample
}
result.PS.count.DZ.IMAC.MBR.ON
glimpse(result.PS.count.DZ.IMAC.MBR.ON)

# add alternative sample name
result.PS.count.DZ.IMAC.MBR.ON$sample2 <- c("GFP.1",
                                            "FL.1",
                                            "dE18.1",
                                            "FL-Bicc1.1",
                                            "dE18-Bicc1.1",
                                            
                                            "GFP.2",
                                            "FL.2",
                                            "dE18.2",
                                            "FL-Bicc1.2",
                                            "dE18-Bicc1.2",
                                            
                                            "GFP.3",
                                            "FL.3",
                                            "dE18.3",
                                            "FL-Bicc1.3",
                                            "dE18-Bicc1.3",
                                            
                                            "GFP.4",
                                            "FL.4",
                                            "dE18.4",
                                            "FL-Bicc1.4",
                                            "dE18-Bicc1.4"
)
result.PS.count.DZ.IMAC.MBR.ON

# calcuate means
mean.PS.allclassPS.per.sample <- mean(result.PS.count.DZ.IMAC.MBR.ON$numberPS); mean.PS.allclassPS.per.sample #8849.8

mean.PS.allclassPS.per.sample.pS <- mean(result.PS.count.DZ.IMAC.MBR.ON$percentage.pS);round(mean.PS.allclassPS.per.sample.pS, 1) #%pS 85.9
mean.PS.allclassPS.per.sample.pT <- mean(result.PS.count.DZ.IMAC.MBR.ON$percentage.pT);round(mean.PS.allclassPS.per.sample.pT, 1) #%pT 11.8
mean.PS.allclassPS.per.sample.pY <- mean(result.PS.count.DZ.IMAC.MBR.ON$percentage.pY);round(mean.PS.allclassPS.per.sample.pY, 1) #%pY 2.3

# mean count all class PS  GFP
result.PS.count.DZ.IMAC.MBR.ON %>% filter(sample %in% c("A", "F", "K", "P")) %>% select(numberPS) %>% pull() %>% mean() #9181
#mean count all class PS  FL
result.PS.count.DZ.IMAC.MBR.ON %>% filter(sample %in% c("B", "G", "L", "Q")) %>% select(numberPS) %>% pull() %>% mean() #8999
#mean count all class PS  dE18
result.PS.count.DZ.IMAC.MBR.ON %>% filter(sample %in% c("C", "H", "M", "R")) %>% select(numberPS) %>% pull() %>% mean() #9026.25
#mean count all class PS  FL-Bicc1.1
result.PS.count.DZ.IMAC.MBR.ON %>% filter(sample %in% c("D", "I", "N", "S")) %>% select(numberPS) %>% pull() %>% mean() #9434.5
#mean count all class PS  dE18-Bicc1.1
result.PS.count.DZ.IMAC.MBR.ON %>% filter(sample %in% c("E", "J", "O", "T")) %>% select(numberPS) %>% pull() %>% mean() #7608.25


### loop over data to count number of PS per sample class 1 PS ###################################################################################################################################

glimpse(DZ.IMAC.MBR.ON.3)

### Class1 PS: filter for PS with localization probability > 0.75
table(DZ.IMAC.MBR.ON.3$Localization.prob >= 0.75) 
DZ.IMAC.MBR.ON.3.class1.COUNT <- DZ.IMAC.MBR.ON.3[which(DZ.IMAC.MBR.ON.3$Localization.prob >= 0.75),]
nrow(DZ.IMAC.MBR.ON.3.class1.COUNT) #13683
table(DZ.IMAC.MBR.ON.3.class1.COUNT$Localization.prob >= 0.75) 

### calculate number of class 1 phosphosite per sample
glimpse(DZ.IMAC.MBR.ON.3.class1.COUNT)
colnames(DZ.IMAC.MBR.ON.3.class1.COUNT)
ncol(DZ.IMAC.MBR.ON.3.class1.COUNT) #80

# loop over data to count number of class 1 PS per sample
result.class1.PS.count.DZ.IMAC.MBR.ON <- c()
sample.name.count <- 1

for(i.1 in seq(from=21, to=80, by=3)){
  print(i.1)
 
  i.2 <- i.1+1
  i.3 <- i.2+1
  
  sample.name <- samples.in.experiment[sample.name.count]
  print(sample.name)
  col.names <- colnames(DZ.IMAC.MBR.ON.3.class1.COUNT)
  col.names.select <- col.names[c(i.1, i.2, i.3)]
  
  temp <- select(DZ.IMAC.MBR.ON.3.class1.COUNT, Amino.acid, col.names.select)
  temp$sum_1_2_3 <- pbapply(temp[, 2:4], 1, function(x) sum(x))
  
  numberPS.A <- filter(temp, sum_1_2_3 > 0)
  numberPS <- nrow(numberPS.A)
  
  nMVs <- filter(temp, sum_1_2_3 <= 0)
  nMVs <- nrow(nMVs); nMVs
  
  percentageSTY <- table(numberPS.A$Amino.acid); percentageSTY
  
  temp.result <- tibble(class1.sample=sample.name, class1.numberPS=numberPS, class1.nMVs=nMVs, class1.nS=percentageSTY[1], class1.nT=percentageSTY[2], class1.nY=percentageSTY[3])
  
  result.class1.PS.count.DZ.IMAC.MBR.ON  <- bind_rows(result.class1.PS.count.DZ.IMAC.MBR.ON , temp.result)
  
  sample.name.count <- sample.name.count +1
}
result.class1.PS.count.DZ.IMAC.MBR.ON 


# calculate percentage MV PS class1
result.class1.PS.count.DZ.IMAC.MBR.ON$percentageMV <- pbsapply(result.class1.PS.count.DZ.IMAC.MBR.ON $class1.nMVs, function(x) (100*x)/13683) 
result.class1.PS.count.DZ.IMAC.MBR.ON 


# calculate percentage pS/pT/pY PS class1
for(i in 1:nrow(result.class1.PS.count.DZ.IMAC.MBR.ON )){
  print(i)
  
  total.PS.sample <- result.class1.PS.count.DZ.IMAC.MBR.ON $class1.numberPS[i]
  n.pS <- result.class1.PS.count.DZ.IMAC.MBR.ON $class1.nS[i]
  n.pT <- result.class1.PS.count.DZ.IMAC.MBR.ON $class1.nT[i]
  n.pY <- result.class1.PS.count.DZ.IMAC.MBR.ON $class1.nY[i]
  
  result.class1.PS.count.DZ.IMAC.MBR.ON $class1.percentage.pS[i] <- (100*n.pS)/total.PS.sample
  result.class1.PS.count.DZ.IMAC.MBR.ON $class1.percentage.pT[i] <- (100*n.pT)/total.PS.sample
  result.class1.PS.count.DZ.IMAC.MBR.ON $class1.percentage.pY[i] <- (100*n.pY)/total.PS.sample
}
glimpse(result.class1.PS.count.DZ.IMAC.MBR.ON)

# add alternative sample name
result.class1.PS.count.DZ.IMAC.MBR.ON$sample2 <- c("GFP.1",
                                                   "FL.1",
                                                   "dE18.1",
                                                   "FL-Bicc1.1",
                                                   "dE18-Bicc1.1",
                                                   
                                                   "GFP.2",
                                                   "FL.2",
                                                   "dE18.2",
                                                   "FL-Bicc1.2",
                                                   "dE18-Bicc1.2",
                                                   
                                                   "GFP.3",
                                                   "FL.3",
                                                   "dE18.3",
                                                   "FL-Bicc1.3",
                                                   "dE18-Bicc1.3",
                                                   
                                                   "GFP.4",
                                                   "FL.4",
                                                   "dE18.4",
                                                   "FL-Bicc1.4",
                                                   "dE18-Bicc1.4")
glimpse(result.class1.PS.count.DZ.IMAC.MBR.ON)

# calcuate means
mean.result.class1.PS.count.DZ.IMAC.MBR.ON  <- mean(result.class1.PS.count.DZ.IMAC.MBR.ON $class1.numberPS); mean.result.class1.PS.count.DZ.IMAC.MBR.ON   #8205.65

#mean count all class PS  GFP
result.class1.PS.count.DZ.IMAC.MBR.ON %>% filter(class1.sample %in% c("A", "F", "K", "P")) %>% select(class1.numberPS) %>% pull() %>% mean() #8490
#mean count all class PS  FL
result.class1.PS.count.DZ.IMAC.MBR.ON %>% filter(class1.sample %in% c("B", "G", "L", "Q")) %>% select(class1.numberPS) %>% pull() %>% mean() #8356
#mean count all class PS  dE18
result.class1.PS.count.DZ.IMAC.MBR.ON %>% filter(class1.sample %in% c("C", "H", "M", "R")) %>% select(class1.numberPS) %>% pull() %>% mean() #8389
#mean count all class PS  FL-Bicc1
result.class1.PS.count.DZ.IMAC.MBR.ON %>% filter(class1.sample %in% c("D", "I", "N", "S")) %>% select(class1.numberPS) %>% pull() %>% mean() #8758.75
#mean count all class PS  dE18-Bicc1
result.class1.PS.count.DZ.IMAC.MBR.ON %>% filter(class1.sample %in% c("E", "J", "O", "T")) %>% select(class1.numberPS) %>% pull() %>% mean() #7034.5 with T #
result.class1.PS.count.DZ.IMAC.MBR.ON %>% filter(class1.sample %in% c("E", "J", "O"     )) %>% select(class1.numberPS) %>% pull() %>% mean() #8716 w/o T #


### plot number of PS CLASS1 AND number of PS (all classes)
ggplot()+
  geom_col(data = result.PS.count.DZ.IMAC.MBR.ON, aes(x=sample2, y=numberPS))+
  geom_col(data= result.class1.PS.count.DZ.IMAC.MBR.ON, aes(x=sample2, y=class1.numberPS), color="black", fill="yellow")+
  scale_x_discrete(limits= c("GFP.1","GFP.2","GFP.3","GFP.4",  "FL.1","FL.2","FL.3","FL.4",  "dE18.1","dE18.2","dE18.3","dE18.4",  "FL-Bicc1.1","FL-Bicc1.2","FL-Bicc1.3","FL-Bicc1.4",   "dE18-Bicc1.1","dE18-Bicc1.2","dE18-Bicc1.3","dE18-Bicc1.4"))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1.0, vjust = 0.5))+
  ylab("number PS")+
  xlab(NULL)+
  scale_y_continuous(breaks=seq(from=0, to=11000, by=500), limits = c(0,11000))+
  ggtitle("# PS (all classes) per sample - grey \n # class 1PS per sample - yellow") +
  theme(plot.title = element_text(hjust = 0.5, face= "bold"))+  
  theme(axis.text.x = element_text(face= "plain", colour="black", size=14))+
  theme(axis.text.y = element_text(face= "plain", colour="black", size=14))+
  geom_vline(xintercept = c(4.5, 8.5, 12.5, 16.5) , size=0.25, linetype = 2)


### plot percentages of pTyr PS
ggplot()+
  geom_col(data = result.PS.count.DZ.IMAC.MBR.ON, aes(x=sample2, y=percentage.pY))+
  geom_col(data= result.class1.PS.count.DZ.IMAC.MBR.ON, aes(x=sample2, y=class1.percentage.pY), color="black", fill="yellow")+
  scale_x_discrete(limits= c("GFP.1","GFP.2","GFP.3","GFP.4",  "FL.1","FL.2","FL.3","FL.4",  "dE18.1","dE18.2","dE18.3","dE18.4",  "FL-Bicc1.1","FL-Bicc1.2","FL-Bicc1.3","FL-Bicc1.4",   "dE18-Bicc1.1","dE18-Bicc1.2","dE18-Bicc1.3","dE18-Bicc1.4"))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1.0, vjust = 0.5))+
  ylab("percentage")+
  xlab(NULL)+
  scale_y_continuous(breaks=seq(from=0, to=5, by=0.5), limits = c(0,5))+
  ggtitle("% pTyr sites (all classes) per sample - grey \n % class 1 pTyr sites per sample - yellow") +
  theme(plot.title = element_text(hjust = 0.5, face= "bold"))+  
  theme(axis.text.x = element_text(face= "plain", colour="black", size=14))+
  theme(axis.text.y = element_text(face= "plain", colour="black", size=14))+
  geom_vline(xintercept = c(4.5, 8.5, 12.5, 16.5) , size=0.25, linetype = 2)


#mean count all class PS  GFP
result.class1.PS.count.DZ.IMAC.MBR.ON %>% filter(class1.sample %in% c("A", "F", "K", "P")) %>% select(class1.percentage.pY) %>% pull() %>% mean() #0.7899692
#mean count all class PS  FL
result.class1.PS.count.DZ.IMAC.MBR.ON %>% filter(class1.sample %in% c("B", "G", "L", "Q")) %>% select(class1.percentage.pY) %>% pull() %>% mean() #2.060383
#mean count all class PS  dE18
result.class1.PS.count.DZ.IMAC.MBR.ON %>% filter(class1.sample %in% c("C", "H", "M", "R")) %>% select(class1.percentage.pY) %>% pull() %>% mean() #2.748839
#mean count all class PS  FL-Bicc1
result.class1.PS.count.DZ.IMAC.MBR.ON %>% filter(class1.sample %in% c("D", "I", "N", "S")) %>% select(class1.percentage.pY) %>% pull() %>% mean() #2.259006
#mean count all class PS  dE18-Bicc1
result.class1.PS.count.DZ.IMAC.MBR.ON %>% filter(class1.sample %in% c("E", "J", "O", "T")) %>% select(class1.percentage.pY) %>% pull() %>% mean() #3.470072 with T #
result.class1.PS.count.DZ.IMAC.MBR.ON %>% filter(class1.sample %in% c("E", "J", "O"     )) %>% select(class1.percentage.pY) %>% pull() %>% mean() #3.102474 w/o T #

result.class1.PS.count.DZ.IMAC.MBR.ON  %>% arrange(match(class1.sample , c("A", "F", "K", "P",    "B", "G", "L", "Q",    "C", "H", "M", "R",    "D", "I", "N", "S",  "E", "J", "O", "T")))


#### cells IMAC MBR ON expand with PerseusPSexpand function (see above) ################################################################################################################
#### cells IMAC MBR ON expand with PerseusPSexpand function (see above) ################################################################################################################


glimpse(DZ.IMAC.MBR.ON.3)

DZ.IMAC.MBR.ON_expanded <- PerseusPSexpand.17.02.21(DZ.IMAC.MBR.ON.3)


glimpse(DZ.IMAC.MBR.ON_expanded)


#remove  sample T (low spectra) & remove rows that have all zero intensities
DZ.IMAC.MBR.ON_expanded <- DZ.IMAC.MBR.ON_expanded %>% select(-c("Intensity.T"))
glimpse(DZ.IMAC.MBR.ON_expanded)
colnames(DZ.IMAC.MBR.ON_expanded)

DZ.IMAC.MBR.ON_expanded$ZEROcheckINT <- pbapply(DZ.IMAC.MBR.ON_expanded[, c(21:39)], MARGIN=1, function(x) sum(x)) #
table(DZ.IMAC.MBR.ON_expanded$ZEROcheckINT > 0)
nrow(DZ.IMAC.MBR.ON_expanded) #51444

DZ.IMAC.MBR.ON_expanded <- filter(DZ.IMAC.MBR.ON_expanded, ZEROcheckINT > 0)
nrow(DZ.IMAC.MBR.ON_expanded) #20097

#remove column ZEROcheckINT
DZ.IMAC.MBR.ON_expanded <- select(DZ.IMAC.MBR.ON_expanded,  -ZEROcheckINT)
glimpse(DZ.IMAC.MBR.ON_expanded)


# cells IMAC  MBR ON normalization & histogram standard median shift version #########################################################################################################################################
# cells IMAC  MBR ON normalization & histogram standard median shift version #########################################################################################################################################

glimpse(DZ.IMAC.MBR.ON_expanded)
colnames(DZ.IMAC.MBR.ON_expanded)

#log2 transform & replace -Inf with NA
colnames(DZ.IMAC.MBR.ON_expanded)
Log2Trafo <- DZ.IMAC.MBR.ON_expanded[, c(21:39)] #
Log2Trafo <- log2(Log2Trafo)
Log2Trafo[Log2Trafo  == -Inf] <- NA

glimpse(Log2Trafo)

DZ.IMAC.MBR.ON_expanded_log2 <- DZ.IMAC.MBR.ON_expanded
DZ.IMAC.MBR.ON_expanded_log2[, c(21:39)] <- Log2Trafo #

glimpse(DZ.IMAC.MBR.ON_expanded_log2, list.len=900)


# calculate median per column RawINT
colnames(DZ.IMAC.MBR.ON_expanded_log2)

ColMedianfindDF <- glimpse(DZ.IMAC.MBR.ON_expanded_log2[, c(21:39)]) #
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
Mdf <- data.frame(variable=colnames(DZ.IMAC.MBR.ON_expanded_log2[, c(21:39)]), median=vectormediancolumnsINT); Mdf
Mdf$variable2 <- c("GFP.1",
                   "FL.1",
                   "dE18.1",
                   "FL-Bicc1.1",
                   "dE18-Bicc1.1",
                   
                   "GFP.2",
                   "FL.2",
                   "dE18.2",
                   "FL-Bicc1.2",
                   "dE18-Bicc1.2",
                   
                   "GFP.3",
                   "FL.3",
                   "dE18.3",
                   "FL-Bicc1.3",
                   "dE18-Bicc1.3",
                   
                   "GFP.4",
                   "FL.4",
                   "dE18.4",
                   "FL-Bicc1.4")

Mdf$group <- c("GFP",
               "FL",
               "dE18",
               "FL-Bicc1",
               "dE18-Bicc1",
               
               "GFP",
               "FL",
               "dE18",
               "FL-Bicc1",
               "dE18-Bicc1",
               
               "GFP",
               "FL",
               "dE18",
               "FL-Bicc1",
               "dE18-Bicc1",
               
               "GFP",
               "FL",
               "dE18",
               "FL-Bicc1")
Mdf

## Plot
DensityPlotData <- DZ.IMAC.MBR.ON_expanded_log2[, c(21:39)]

glimpse(DensityPlotData)

DensityPlotDataMelted <- reshape2::melt(DensityPlotData)
glimpse(DensityPlotDataMelted)

DensityPlot  <- ggplot(DensityPlotDataMelted, aes(x=value, colour = variable) ) + 
  geom_density() +
  theme(legend.position="bottom") +
  geom_vline(data=Mdf, aes(xintercept=median, color=variable))
DensityPlot 



## Start define function to calculate median shift factors 
MedianShiftFactorFinder <- function(arg){
  tempfactor <- 23.2535/arg # log2(1e07) = 23.2535
  return(tempfactor)
}
## Stop define function to calculate median shift factors 

## find median shift factors
shiftfactors <- c()
for(i in vectormediancolumnsINT ){
  print(i)
  temp <- MedianShiftFactorFinder(arg = i)
  shiftfactors <- append(shiftfactors, temp)
}
shiftfactors 


## apply median shift factors
temp_for_median_shift <-  DZ.IMAC.MBR.ON_expanded_log2[, c(21:39)]
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

DensityPlotCHECK  <- ggplot(DensityPlotDataCHECKMelted, aes(x=value, colour = variable) ) + 
  geom_density() +
  theme(legend.position="bottom") +
  geom_vline(data=MdfCHECK, aes(xintercept=median, color=variable))
DensityPlotCHECK


### rename columns after normalization
newCOLNAMESnorm <- colnames(temp_for_median_shift); newCOLNAMESnorm

newCOLNAMESnorm <- paste0("Norm.",newCOLNAMESnorm); newCOLNAMESnorm

colnames(temp_for_median_shift) <- newCOLNAMESnorm
glimpse(temp_for_median_shift)


### add norm INT to original Dataframe 

DZ.IMAC.MBR.ON_expanded_log2_norm <- DZ.IMAC.MBR.ON_expanded_log2

DZ.IMAC.MBR.ON_expanded_log2_norm  <- cbind(DZ.IMAC.MBR.ON_expanded_log2_norm , temp_for_median_shift)
glimpse(DZ.IMAC.MBR.ON_expanded_log2_norm, list.len=999)


## box plot after normalization
data.for.boxplot.after.norm.DZ.IMAC.MBR.ON <- DZ.IMAC.MBR.ON_expanded_log2_norm %>% select(contains("Norm.Intensity."))
glimpse(data.for.boxplot.after.norm.DZ.IMAC.MBR.ON)
colnames(data.for.boxplot.after.norm.DZ.IMAC.MBR.ON) <- c("GFP.1",
                                                          "FL.1",
                                                          "dE18.1",
                                                          "FL-Bicc1.1",
                                                          "dE18-Bicc1.1",
                                                          
                                                          "GFP.2",
                                                          "FL.2",
                                                          "dE18.2",
                                                          "FL-Bicc1.2",
                                                          "dE18-Bicc1.2",
                                                          
                                                          "GFP.3",
                                                          "FL.3",
                                                          "dE18.3",
                                                          "FL-Bicc1.3",
                                                          "dE18-Bicc1.3",
                                                          
                                                          "GFP.4",
                                                          "FL.4",
                                                          "dE18.4",
                                                          "FL-Bicc1.4")
glimpse(data.for.boxplot.after.norm.DZ.IMAC.MBR.ON)

melted.data.for.boxplot.after.norm.DZ.IMAC.MBR.ON <- reshape2::melt(data.for.boxplot.after.norm.DZ.IMAC.MBR.ON)
glimpse(melted.data.for.boxplot.after.norm.DZ.IMAC.MBR.ON )
unique(melted.data.for.boxplot.after.norm.DZ.IMAC.MBR.ON$variable)

ggplot(data = melted.data.for.boxplot.after.norm.DZ.IMAC.MBR.ON)+
  geom_boxplot(mapping = aes(x= variable, y = value, fill=variable), outlier.size = 0.5, notch = TRUE, notchwidth = 0.1)+
  scale_x_discrete(limits= c("GFP.1","GFP.2","GFP.3","GFP.4",  "FL.1","FL.2","FL.3","FL.4",  "dE18.1","dE18.2","dE18.3","dE18.4",  "FL-Bicc1.1","FL-Bicc1.2","FL-Bicc1.3","FL-Bicc1.4",   "dE18-Bicc1.1","dE18-Bicc1.2","dE18-Bicc1.3"))+
  scale_fill_manual(values = rep(c("darkgreen", "blue", "red", "black", "grey40"),4) )+ 
  theme(axis.text.x = element_text(angle = 90, hjust = 1.0, vjust = 0.5))+
  theme(legend.position="none", legend.justification = "center")+
  ggtitle("cells 1 DZ IMAC - normalized data - PS level") +
  theme(plot.title = element_text(hjust = 0.5, face= "bold"))+  
  theme(axis.text.x = element_text(face= "plain", colour="black", size=14))+
  theme(axis.text.y = element_text(face= "plain", colour="black", size=14))+
  xlab(NULL)+ 
  ylab("log2(Int.)") 



#####################################################################################################################################################################
#### analysis via class 1 PS ##################################################################################################################################
glimpse(DZ.IMAC.MBR.ON_expanded_log2_norm) #

DZ.IMAC.MBR.ON_expanded_log2_norm.class1.PS <- DZ.IMAC.MBR.ON_expanded_log2_norm %>% filter(Localization.prob>= 0.75)
glimpse(DZ.IMAC.MBR.ON_expanded_log2_norm.class1.PS) # 16430


# add additional column to better calculate number PS 
# take Proteins and add Psite without multiplicity info ###
for(i in 1:nrow(DZ.IMAC.MBR.ON_expanded_log2_norm.class1.PS)){
  print(i)
  temp_prot_name <- DZ.IMAC.MBR.ON_expanded_log2_norm.class1.PS$Proteins[i];temp_prot_name
  tempPSite <- DZ.IMAC.MBR.ON_expanded_log2_norm.class1.PS$Positions[i]; tempPSite
  tempAA <- DZ.IMAC.MBR.ON_expanded_log2_norm.class1.PS$Amino.acid[i];  tempAA
  DZ.IMAC.MBR.ON_expanded_log2_norm.class1.PS$Prot.Name_AApos[i] <- paste0(temp_prot_name, "_", tempAA, tempPSite)
}
length(unique(DZ.IMAC.MBR.ON_expanded_log2_norm.class1.PS$Prot.Name_AApos))

glimpse(DZ.IMAC.MBR.ON_expanded_log2_norm.class1.PS)

#  prepare new column
#  take genename and add Psite ###
for(i in 1:nrow(DZ.IMAC.MBR.ON_expanded_log2_norm.class1.PS)){
  print(i)
  #i=60
  temp_gene_name <- DZ.IMAC.MBR.ON_expanded_log2_norm.class1.PS$Gene.names[i];temp_gene_name
  tempPSite <- DZ.IMAC.MBR.ON_expanded_log2_norm.class1.PS$Positions[i]; tempPSite
  tempAA <- DZ.IMAC.MBR.ON_expanded_log2_norm.class1.PS$Amino.acid[i];  tempAA
  tempPSmultiplicity <- DZ.IMAC.MBR.ON_expanded_log2_norm.class1.PS$PS_Multiplicity[i]; tempPSmultiplicity
  DZ.IMAC.MBR.ON_expanded_log2_norm.class1.PS$Gene.Name_AApos[i] <- paste0(temp_gene_name, "_", tempAA, tempPSite, "x", tempPSmultiplicity)
}
glimpse(DZ.IMAC.MBR.ON_expanded_log2_norm.class1.PS)

#add identifier for plotting reasons e.g. PTMSEA later
glimpse(DZ.IMAC.MBR.ON_expanded_log2_norm.class1.PS)

DZ.IMAC.MBR.ON_expanded_log2_norm.class1.PS$ID.FR.all.C1.PS <- 1:nrow(DZ.IMAC.MBR.ON_expanded_log2_norm.class1.PS)

#save data
#save(DZ.IMAC.MBR.ON_expanded_log2_norm.class1.PS, file="DZ.IMAC.MBR.ON_expanded_log2_norm.class1.PS.Rdata")
#load("DZ.IMAC.MBR.ON_expanded_log2_norm.class1.PS.Rdata")
#glimpse(DZ.IMAC.MBR.ON_expanded_log2_norm.class1.PS) #16430

############################################################################################################################################################
############################################################################################################################################################

### heatmap for all Fgfr2 sites in dataset for all groups with collapsing _1/_2/_3 via sum 
### PS numbers also for NP_963895.2

glimpse( DZ.IMAC.MBR.ON_expanded_log2_norm.class1.PS)

HM.plot.data <- DZ.IMAC.MBR.ON_expanded_log2_norm.class1.PS
HM.plot.data <- HM.plot.data %>% filter(stringr::str_detect(Gene.Name_AApos, "Fgfr2"))
glimpse(HM.plot.data)

#pick columns of interest and rename
HM.plotdata <- HM.plot.data %>% select("Gene.Name_AApos", 
                                       "Norm.Intensity.A", "Norm.Intensity.F", "Norm.Intensity.K", "Norm.Intensity.P", 
                                       "Norm.Intensity.B", "Norm.Intensity.G", "Norm.Intensity.L", "Norm.Intensity.Q",    
                                       "Norm.Intensity.C", "Norm.Intensity.H", "Norm.Intensity.M", "Norm.Intensity.R", 
                                       "Norm.Intensity.D", "Norm.Intensity.I", "Norm.Intensity.N", "Norm.Intensity.S",    
                                       "Norm.Intensity.E", "Norm.Intensity.J", "Norm.Intensity.O",
                                       "Proteins",
                                       "Positions.within.proteins",
                                       "Amino.acid" , 
                                       "Gene.names",
                                       "PS_Multiplicity"
) 
glimpse(HM.plotdata)

colnames(HM.plotdata) <- c("Gene.Name_AApos", 
                           "GFP.1","GFP.2","GFP.3","GFP.4",  
                           "FL.1","FL.2","FL.3","FL.4",  
                           "dE18.1","dE18.2","dE18.3","dE18.4",  
                           "FL-Bicc1.1","FL-Bicc1.2","FL-Bicc1.3","FL-Bicc1.4",   
                           "dE18-Bicc1.1","dE18-Bicc1.2","dE18-Bicc1.3",
                           "Proteins",
                           "Positions.within.proteins",
                           "Amino.acid",
                           "Gene.names" ,
                           "PS_Multiplicity")
glimpse(HM.plotdata)



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
for(i in 1:nrow(HM.plotdata)){
  print(i)
  HM.plotdata$isoform.1.corrected.PS.position[i] <- position.in.mq.string(prot.of.interest ="P21803", #isoform 1
                                                                          Gene.names=HM.plotdata$Gene.names[i],
                                                                          Amino.acid=HM.plotdata$Amino.acid[i], 
                                                                          mq.position.string=HM.plotdata$Positions.within.proteins[i],
                                                                          mq.string = HM.plotdata$Proteins[i],
                                                                          PS_Multiplicity=HM.plotdata$PS_Multiplicity[i])
}
HM.plotdata$isoform.1.corrected.Gene.Name_AApos = paste0(HM.plotdata$Gene.names, "_", HM.plotdata$Amino.acid, HM.plotdata$isoform.1.corrected.PS.position, "x", HM.plotdata$PS_Multiplicity)
glimpse(HM.plotdata)

# isoform 2 ="P21803-2"
for(i in 1:nrow(HM.plotdata)){
  print(i)
  HM.plotdata$isoform.2.corrected.PS.position[i] <- position.in.mq.string(prot.of.interest ="P21803-2", #isoform 2
                                                                          Gene.names=HM.plotdata$Gene.names[i],
                                                                          Amino.acid=HM.plotdata$Amino.acid[i], 
                                                                          mq.position.string=HM.plotdata$Positions.within.proteins[i],
                                                                          mq.string = HM.plotdata$Proteins[i],
                                                                          PS_Multiplicity=HM.plotdata$PS_Multiplicity[i])
}
HM.plotdata$isoform.2.corrected.Gene.Name_AApos = paste0(HM.plotdata$Gene.names, "_", HM.plotdata$Amino.acid, HM.plotdata$isoform.2.corrected.PS.position, "x", HM.plotdata$PS_Multiplicity)

glimpse(HM.plotdata)

#prepare helper data frame for later order
HM.plotdata.melted <- reshape2::melt(HM.plotdata)
glimpse(HM.plotdata.melted)
colnames(HM.plotdata.melted) <- c("Gene.Name_AApos", "Proteins","Positions.within.proteins", "Amino.acid","Gene.names", "PS_Multiplicity", "isoform.1.corrected.PS.position", "isoform.1.corrected.Gene.Name_AApos","isoform.2.corrected.PS.position", "isoform.2.corrected.Gene.Name_AApos",  "variable","log2.int")
glimpse(HM.plotdata.melted)


order.df.HM.plotdata.melted <- unique(tibble(Gene.Name_AApos = HM.plotdata.melted$Gene.Name_AApos,
                                             Proteins = HM.plotdata.melted$Proteins,
                                             Positions.within.proteins = HM.plotdata.melted$Positions.within.proteins,
                                             Amino.acid = HM.plotdata.melted$Amino.acid,
                                             Gene.names = HM.plotdata.melted$Gene.names,
                                             PS_Multiplicity =HM.plotdata.melted$PS_Multiplicity,
                                             isoform.1.corrected.PS.position =HM.plotdata.melted$isoform.1.corrected.PS.position,
                                             isoform.1.corrected.Gene.Name_AApos = paste0(HM.plotdata.melted$Gene.names, "_", HM.plotdata.melted$Amino.acid, HM.plotdata.melted$isoform.1.corrected.PS.position, "x", HM.plotdata.melted$PS_Multiplicity),
                                             isoform.2.corrected.PS.position =HM.plotdata.melted$isoform.2.corrected.PS.position,
                                             isoform.2.corrected.Gene.Name_AApos = paste0(HM.plotdata.melted$Gene.names, "_", HM.plotdata.melted$Amino.acid, HM.plotdata.melted$isoform.2.corrected.PS.position, "x", HM.plotdata.melted$PS_Multiplicity)
))
glimpse(order.df.HM.plotdata.melted)


### for better visualization/plotting => collapse _1/_2/_3 via  sum
### note: sum is based on _1/_2/_3 intensities in dataframe at hand => unfiltered data might contain more _1/_2/_3 PS that are not included in the sum here

HM.plotdata.collapsed <- HM.plotdata
glimpse(HM.plotdata.collapsed)

HM.plotdata.collapsed$isoform.1.corrected.Gene.Name_AApos.collapsed <- pbsapply(HM.plotdata.collapsed$isoform.1.corrected.Gene.Name_AApos, function(x) str_remove_all(string=x, pattern="x_\\d") )
HM.plotdata.collapsed$isoform.2.corrected.Gene.Name_AApos.collapsed <- pbsapply(HM.plotdata.collapsed$isoform.2.corrected.Gene.Name_AApos, function(x) str_remove_all(string=x, pattern="x_\\d") )
glimpse(HM.plotdata.collapsed)
## unlog, do collapsing, relog
HM.plotdata.collapsed[, c("GFP.1","GFP.2","GFP.3","GFP.4","FL.1","FL.2","FL.3","FL.4","dE18.1","dE18.2","dE18.3","dE18.4","FL-Bicc1.1","FL-Bicc1.2","FL-Bicc1.3","FL-Bicc1.4","dE18-Bicc1.1","dE18-Bicc1.2","dE18-Bicc1.3")] <- 2^HM.plotdata.collapsed [, c("GFP.1","GFP.2","GFP.3","GFP.4","FL.1","FL.2","FL.3","FL.4","dE18.1","dE18.2","dE18.3","dE18.4","FL-Bicc1.1","FL-Bicc1.2","FL-Bicc1.3","FL-Bicc1.4","dE18-Bicc1.1","dE18-Bicc1.2","dE18-Bicc1.3")] 
glimpse(HM.plotdata.collapsed )
# collapse
HM.plotdata.collapsed.2  <- HM.plotdata.collapsed %>% group_by(isoform.1.corrected.Gene.Name_AApos.collapsed, isoform.2.corrected.Gene.Name_AApos.collapsed) %>% summarise_at(c("GFP.1","GFP.2","GFP.3","GFP.4","FL.1","FL.2","FL.3","FL.4","dE18.1","dE18.2","dE18.3","dE18.4","FL-Bicc1.1","FL-Bicc1.2","FL-Bicc1.3","FL-Bicc1.4","dE18-Bicc1.1","dE18-Bicc1.2","dE18-Bicc1.3"), sum, na.rm = TRUE)
glimpse(HM.plotdata.collapsed.2)
HM.plotdata.collapsed.2[HM.plotdata.collapsed.2  == 0] <- NA
# relog
HM.plotdata.collapsed.2[c("GFP.1","GFP.2","GFP.3","GFP.4","FL.1","FL.2","FL.3","FL.4","dE18.1","dE18.2","dE18.3","dE18.4","FL-Bicc1.1","FL-Bicc1.2","FL-Bicc1.3","FL-Bicc1.4","dE18-Bicc1.1","dE18-Bicc1.2","dE18-Bicc1.3")] <- log2(HM.plotdata.collapsed.2[c("GFP.1","GFP.2","GFP.3","GFP.4","FL.1","FL.2","FL.3","FL.4","dE18.1","dE18.2","dE18.3","dE18.4","FL-Bicc1.1","FL-Bicc1.2","FL-Bicc1.3","FL-Bicc1.4","dE18-Bicc1.1","dE18-Bicc1.2","dE18-Bicc1.3")] )
glimpse(HM.plotdata.collapsed.2)


### reshape and reorder, define order for plot
HM.plotdata.collapsed.melted <- reshape2::melt(HM.plotdata.collapsed.2)
glimpse(HM.plotdata.collapsed.melted)
colnames(HM.plotdata.collapsed.melted ) <- c("isoform.1.corrected.Gene.Name_AApos.collapsed","isoform.2.corrected.Gene.Name_AApos.collapsed", "variable","log2.int")
glimpse(HM.plotdata.collapsed.melted )



#define order for plot, get gene names, corrected positions, aminco acid etc. and order accor4ding to number
order.df.HM.plotdata.collapsed.melted <- unique(tibble(isoform.1.corrected.Gene.Name_AApos.collapsed = HM.plotdata.collapsed.melted$isoform.1.corrected.Gene.Name_AApos.collapsed,
                                                       isoform.2.corrected.Gene.Name_AApos.collapsed = HM.plotdata.collapsed.melted$isoform.2.corrected.Gene.Name_AApos.collapsed))
order.df.HM.plotdata.collapsed.melted

order.df.HM.plotdata.collapsed.melted$Gene <- pbsapply(order.df.HM.plotdata.collapsed.melted$isoform.1.corrected.Gene.Name_AApos.collapsed , function(x) unlist(str_split(string=x,  pattern="_") )[1])

order.df.HM.plotdata.collapsed.melted$isoform.1.corrected.AAPos <- pbsapply(order.df.HM.plotdata.collapsed.melted$isoform.1.corrected.Gene.Name_AApos.collapsed, function(x) unlist(str_split(string=x,  pattern="_") )[2])
order.df.HM.plotdata.collapsed.melted$isoform.2.corrected.AAPos <- pbsapply(order.df.HM.plotdata.collapsed.melted$isoform.2.corrected.Gene.Name_AApos.collapsed, function(x) unlist(str_split(string=x,  pattern="_") )[2])

order.df.HM.plotdata.collapsed.melted$isoform.corrected.Amino.acid <- pbsapply(order.df.HM.plotdata.collapsed.melted$isoform.1.corrected.AAPos, function(x) unlist((str_extract_all(string=x,  pattern="[A-Z]") )))

order.df.HM.plotdata.collapsed.melted$isoform.1.corrected.Pos <- pbsapply(order.df.HM.plotdata.collapsed.melted$isoform.1.corrected.AAPos, function(x) unlist((str_extract_all(string=x,  pattern="[0-9]+") )))
order.df.HM.plotdata.collapsed.melted$isoform.2.corrected.Pos <- pbsapply(order.df.HM.plotdata.collapsed.melted$isoform.2.corrected.AAPos, function(x) unlist((str_extract_all(string=x,  pattern="[0-9]+") )))


order.df.HM.plotdata.collapsed.melted <- order.df.HM.plotdata.collapsed.melted %>% arrange(desc(isoform.2.corrected.Pos, isoform.2.corrected.Amino.acid))
order.df.HM.plotdata.collapsed.melted


isoformcorrected.collapes.order.29.11.21 <- order.df.HM.plotdata.collapsed.melted$isoform.2.corrected.Gene.Name_AApos.collapsed
isoformcorrected.collapes.order.29.11.21


glimpse(HM.plotdata.collapsed.melted)

#plot the heatmap with site numbers for P21803-2
ggplot(
  HM.plotdata.collapsed.melted , aes(x=variable, y=isoform.2.corrected.Gene.Name_AApos.collapsed, fill=log2.int)) +  #fct_inorder(PG.ProteinNames)
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
  theme(axis.text.x = element_blank()) +
  theme(axis.ticks.x = element_blank())+
  geom_vline(xintercept = c(4.5, 8.5, 12.5, 16.5) , size=1.0, linetype = "solid", color="white")+
  annotate("text", x = 2.5, y = nrow(HM.plotdata.collapsed.2)+2, label = "GFP", fontface = "bold", size=2.5)+
  annotate("text", x = 6.5, y = nrow(HM.plotdata.collapsed.2)+2, label = "FL", fontface = "bold", size=2.5)+
  annotate("text", x = 10.5, y = nrow(HM.plotdata.collapsed.2)+2, label = "dE18", fontface = "bold", size=2.5)+
  annotate("text", x = 14.5, y = nrow(HM.plotdata.collapsed.2)+2, label = "FL\nBicc1", fontface = "bold", size=2.5)+
  annotate("text", x = 18.0, y = nrow(HM.plotdata.collapsed.2)+2, label = "dE18\nBicc1", fontface = "bold", size=2.5)+
  annotate("text", x = 22.5, y = nrow(HM.plotdata.collapsed.2)+3.5, label = "",  color = "transparent")

### add corresponding human FGFR2 phosphosite numbers to the heatmap plot ######################################################

glimpse(HM.plotdata.collapsed.2)
order.df.HM.plotdata.collapsed.melted

# online Clustal Omega alignment showed that P21803-2 (707 aa) and NP_963895.2 (726 aa) are very similar. The difference is N-terminal: 
# 19aa  (MGLPSTWRYGRGPGIGTVT) are added for  NP_963895.2. Everything else is identical.
# Thus P21803-2 position +19 should be the positions for NP_963895.2

#prepare additional data frame and calculate the new site numbers
HM.plotdata.collapsed.3 <- HM.plotdata.collapsed.2

HM.plotdata.collapsed.3$isoform.1.corrected.AApos <- pbsapply(HM.plotdata.collapsed.3$isoform.1.corrected.Gene.Name_AApos.collapsed, function(x) unlist(str_split(string=x, pattern="_"))[2])
HM.plotdata.collapsed.3$isoform.1.corrected.pos <- pbsapply(HM.plotdata.collapsed.3$isoform.1.corrected.AApos, function(x) unlist(str_extract_all(string=x, pattern="[0-9]+")))

HM.plotdata.collapsed.3$isoform.2.corrected.AApos <- pbsapply(HM.plotdata.collapsed.3$isoform.2.corrected.Gene.Name_AApos.collapsed, function(x) unlist(str_split(string=x, pattern="_"))[2])
HM.plotdata.collapsed.3$isoform.2.corrected.pos <- pbsapply(HM.plotdata.collapsed.3$isoform.2.corrected.AApos, function(x) unlist(str_extract_all(string=x, pattern="[0-9]+")))
glimpse(HM.plotdata.collapsed.3)


HM.plotdata.collapsed.3$pos.NP_963895.2 <- as.numeric(HM.plotdata.collapsed.3$isoform.2.corrected.pos) +19
glimpse(HM.plotdata.collapsed.3)

# go to www.phosphosite.org/ => search for Fgfr2 => look a t site table => enable other species/isoforms mouse => copy list to Excel for lookup of human phosphosite numbers starting from mouse
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
HM.plotdata.collapsed.4 <- merge(HM.plotdata.collapsed.3, human.PSP.site.number.lookup.2[, c("site.Fgfr2.mouse.iso1", "site.FGFR2.human.iso1" )], by.x = "isoform.1.corrected.AApos", by.y = "site.Fgfr2.mouse.iso1", all.x = T, all.y = F, sort = F)
HM.plotdata.collapsed.4$Gene <- pbsapply(HM.plotdata.collapsed.4$isoform.1.corrected.Gene.Name_AApos.collapsed, function(x) unlist(str_split(string=x,  pattern="_") )[1])
HM.plotdata.collapsed.4 <- HM.plotdata.collapsed.4  %>% arrange(HM.plotdata.collapsed.4$isoform.1.corrected.pos)

glimpse(HM.plotdata.collapsed.4)

#some human sites have to be entered manually since not all sites are in PSP yet (looked up in online Clustal Omega alignment)
#they correspond mostly to numbers from P21803-1
unname(HM.plotdata.collapsed.4$isoform.1.corrected.AApos)
HM.plotdata.collapsed.4$site.FGFR2.human.iso1

HM.plotdata.collapsed.4$AA.site.FGFR2.human.iso1_manual.change <- c("T419","T429","S452","S453","T457","S464","Y466","Y561","Y575","Y586","S587","Y588","Y608","Y616","Y656","Y657","Y733" )
HM.plotdata.collapsed.4$site.FGFR2.human.iso1_manual.change <- pbsapply(HM.plotdata.collapsed.4$AA.site.FGFR2.human.iso1_manual.change, function(x) str_remove_all(string=x, pattern="[A-Z]"))
glimpse(HM.plotdata.collapsed.4)

#make.new.plot.label
#e.g. Fgfr2_T419|305|324|419 => 1st number: P21803-1; 2nd number: P21803-2; 3rd number: NP_963895.2; 4th number: P21802-1
HM.plotdata.collapsed.4$combined.plot.label <- paste0(HM.plotdata.collapsed.4$isoform.1.corrected.Gene.Name_AApos.collapsed, "|", 
                                                      HM.plotdata.collapsed.4$isoform.2.corrected.pos, "|",  
                                                      HM.plotdata.collapsed.4$pos.NP_963895.2, "|", 
                                                      HM.plotdata.collapsed.4$site.FGFR2.human.iso1_manual.change)
glimpse(HM.plotdata.collapsed.4)

#order for plot
isoform.corrected.collapsed.order.02.12.2021 <- rev(HM.plotdata.collapsed.4$combined.plot.label)

# reshape for plot
melted.HM.plotdata.collapsed.4 <- reshape2::melt(HM.plotdata.collapsed.4 %>% select(combined.plot.label, GFP.1    : `dE18-Bicc1.3`))
glimpse(melted.HM.plotdata.collapsed.4)

# plot
ggplot(
  melted.HM.plotdata.collapsed.4 , aes(x=variable, y=combined.plot.label, fill=value)) +  
  geom_tile() + 
  scale_fill_viridis( option = "inferno", na.value = "grey30",name="log2(intensity)") + 
  coord_equal()+
  scale_y_discrete(limits=  isoform.corrected.collapsed.order.02.12.2021 )+
  theme(legend.position="bottom", legend.justification = "center")+
  theme(axis.line = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  theme(legend.title.align=0.5)+
  theme(plot.title = element_text(hjust = 0.5, face="bold", size=16)) + #center title
  theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1))+  ## Vertical text on x axis
  theme(axis.text.y= element_text(size=8))+  ## Vertical text on x axis
  xlab(NULL) + 
  xlab(NULL) + 
  ylab(NULL) +
  theme(axis.text.x = element_blank()) +
  theme(axis.ticks.x = element_blank())+
  
  geom_vline(xintercept = c(4.5, 8.5, 12.5, 16.5) , size=1.0, linetype = "solid", color="white")+
  annotate("text", x = 2.5, y = nrow(HM.plotdata.collapsed.4)+2, label = "GFP", fontface = "bold", size=2.5)+
  annotate("text", x = 6.5, y = nrow(HM.plotdata.collapsed.4)+2, label = "FL", fontface = "bold", size=2.5)+
  annotate("text", x = 10.5, y = nrow(HM.plotdata.collapsed.4)+2, label = "dE18", fontface = "bold", size=2.5)+
  annotate("text", x = 14.5, y = nrow(HM.plotdata.collapsed.4)+2, label = "FL\nBicc1", fontface = "bold", size=2.5)+
  annotate("text", x = 18.0, y = nrow(HM.plotdata.collapsed.4)+2, label = "dE18\nBicc1", fontface = "bold", size=2.5)+
  annotate("text", x = 22.5, y = nrow(HM.plotdata.collapsed.4)+3.5, label = "",  color = "transparent")

#ggsave("cells.IMAC.global.phospho.Fgfr2.sites.collapsed.with.isoform.numbers.pdf", useDingbats=FALSE,  width = 10, height =12, units = "cm")







### cells IMAC correlation #################################################################################################################################
### cells IMAC correlation #################################################################################################################################

### correlation without Bicc1 sites 

glimpse(DZ.IMAC.MBR.ON_expanded_log2_norm.class1.PS)

#define sites to exclude
sites_excluded_for_correlation_plot <- DZ.IMAC.MBR.ON_expanded_log2_norm.class1.PS %>% filter(stringr::str_detect(Gene.Name_AApos, "Bicc1"))
glimpse(sites_excluded_for_correlation_plot) #60 rows

#exclude Bicc1 sites
DF_for_correlation_plot <- DZ.IMAC.MBR.ON_expanded_log2_norm.class1.PS %>% filter(!stringr::str_detect(Gene.Name_AApos, "Bicc1"))
glimpse(DF_for_correlation_plot) #16,370 rows

#select columns and order
DF_for_correlation_plot <- DF_for_correlation_plot %>% select("Norm.Intensity.A", "Norm.Intensity.F", "Norm.Intensity.K", "Norm.Intensity.P", 
                                                              "Norm.Intensity.B", "Norm.Intensity.G", "Norm.Intensity.L", "Norm.Intensity.Q",    
                                                              "Norm.Intensity.C", "Norm.Intensity.H", "Norm.Intensity.M", "Norm.Intensity.R", 
                                                              "Norm.Intensity.D", "Norm.Intensity.I", "Norm.Intensity.N", "Norm.Intensity.S",    
                                                              "Norm.Intensity.E", "Norm.Intensity.J", "Norm.Intensity.O"                     )

glimpse(DF_for_correlation_plot)

#check if there are all NA rows # 19 samples
DF_for_correlation_plot$na.rows <- pbapply(DF_for_correlation_plot, 1, function(x) sum(is.na(x)))
glimpse(DF_for_correlation_plot)
table(DF_for_correlation_plot$na.rows)

DF_for_correlation_plot <- DF_for_correlation_plot %>% select(-na.rows )
glimpse(DF_for_correlation_plot)

colnames(DF_for_correlation_plot) <- c("GFP.1","GFP.2","GFP.3","GFP.4",  "FL.1","FL.2","FL.3","FL.4",  "dE18.1","dE18.2","dE18.3","dE18.4",  "FL-Bicc1.1","FL-Bicc1.2","FL-Bicc1.3","FL-Bicc1.4",   "dE18-Bicc1.1","dE18-Bicc1.2","dE18-Bicc1.3") #

#correlate
corr2 <- cor(DF_for_correlation_plot, method = "pearson", use = "na.or.complete") 
glimpse(corr2)
head(corr2)
min(corr2) 
max(corr2)
round(min(corr2), 2) #0.83

# prepare annotation data
annot_df_for_heatmap <- data.frame(samples = colnames(DF_for_correlation_plot),
                                   group = c(rep("GFP", 4), rep("FL", 4), rep("dE18", 4), rep("FL-Bicc1", 4), rep("dE18-Bicc1", 3))
                                   )

annot_df_for_heatmap 

# make shorter annotation data frame
annot_df_for_heatmap.short <- data.frame( 
  group  = annot_df_for_heatmap$group
)

glimpse(annot_df_for_heatmap.short)
annot_df_for_heatmap.short

# define colors for annotation bar
annot.colors = list(group = c("GFP"="yellow", "FL"="blue", "dE18"="red", "FL-Bicc1"="black", "dE18-Bicc1"="grey40")
)
annot.colors

# create the heatmap annotation
ha.fr <- HeatmapAnnotation(group=annot_df_for_heatmap.short$group,
                           treatment=annot_df_for_heatmap.short$treatment,
                           col = annot.colors,
                           annotation_legend_param = list(grid_height = unit(8, "mm"))
)
ha.fr

# define colors
FR.heatmap.colors.2 <-colorRamp2(c(min(corr2) , 1.0), c("grey90", "grey10"))

#pdf("cells.IMAC.global.phospho.correlation.no.Bicc1.sites.pdf", width=19/2.54, height=17/2.54, useDingbats=FALSE)
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
####################################################################################################################################
####################################################################################################################################
##############################################################################################################################################################
### cells IMAC single sample PTMSEA (ssPTMSEA) using mouse PTMSigDB

## cells IMAC: load normalized class 1 data
#load("DZ.IMAC.MBR.ON_expanded_log2_norm.class1.PS.Rdata")

glimpse(DZ.IMAC.MBR.ON_expanded_log2_norm.class1.PS) #16430

DZ.IMAC.MBR.ON_expanded_log2_norm.class1.PS_ssPTMSEA <- DZ.IMAC.MBR.ON_expanded_log2_norm.class1.PS

### sample overview
# A	1	  GFP
# F	2	  GFP
# K	3	  GFP
# P	4	  GFP
# B	5	  Fgfr2-FL
# G	6	  Fgfr2-FL
# L	7	  Fgfr2-FL
# Q	8	  Fgfr2-FL
# C	9	  Fgfr2-dE18
# H	10	Fgfr2-dE18
# M	11	Fgfr2-dE18
# R	12	Fgfr2-dE18
# D	13	Fgfr2-FL-Bicc1
# I	14	Fgfr2-FL-Bicc1
# N	15	Fgfr2-FL-Bicc1
# S	16	Fgfr2-FL-Bicc1
# E	17	Fgfr2-dE18-Bicc1
# J	18	Fgfr2-dE18-Bicc1
# O	19	Fgfr2-dE18-Bicc1
# T	20	Fgfr2-dE18-Bicc1


### start define function to make aa15.window
make.aa.15.window.from.31.sequence.window <- function(input){
  nchar(input)
  input <- unlist(str_split(string=input, pattern="")); input
  return(paste0(input[9:23], collapse="") )
}
make.aa.15.window.from.31.sequence.window(input="RLETSTSCFYQPQRRsVILDGRSGRQIE___")
make.aa.15.window.from.31.sequence.window(input="RLETSTSCFYQPQRRsVILDGRSGRQIE___")
### stop define function to make aa15.window



### prepare  aa15.window from 31 aa Sequence.window MaxQuant
# take first entry sequence window if there are two entries
DZ.IMAC.MBR.ON_expanded_log2_norm.class1.PS_ssPTMSEA$aa15.window <- pbsapply(DZ.IMAC.MBR.ON_expanded_log2_norm.class1.PS_ssPTMSEA$Sequence.window, function(x) unlist(str_split(string =x, pattern = ";"))[1]  )
DZ.IMAC.MBR.ON_expanded_log2_norm.class1.PS_ssPTMSEA$aa15.window <- pbsapply(DZ.IMAC.MBR.ON_expanded_log2_norm.class1.PS_ssPTMSEA$aa15.window, function(x) make.aa.15.window.from.31.sequence.window(input = x)  )

glimpse(DZ.IMAC.MBR.ON_expanded_log2_norm.class1.PS_ssPTMSEA)

table(duplicated(DZ.IMAC.MBR.ON_expanded_log2_norm.class1.PS_ssPTMSEA$aa15.window)) #2792 duplicated
which(duplicated(DZ.IMAC.MBR.ON_expanded_log2_norm.class1.PS_ssPTMSEA$aa15.window))

### prepare new data frame and remove aa15-window duplicates
DZ.IMAC.MBR.ON_expanded_log2_norm.class1.PS_ssPTMSEA_2 <- DZ.IMAC.MBR.ON_expanded_log2_norm.class1.PS_ssPTMSEA
glimpse(DZ.IMAC.MBR.ON_expanded_log2_norm.class1.PS_ssPTMSEA_2)

# find unique aa15.windows
unique.15aa.windows <- unique(DZ.IMAC.MBR.ON_expanded_log2_norm.class1.PS_ssPTMSEA_2$aa15.window)
glimpse(unique.15aa.windows) #13638

# PTMSEA does not accept duplicate aa15.windows
# if there are duplicate entries (e.g. because of _1_2_3) here, use the entry with the highest sum intensity over all samples, 
# this will favour enties with low NAs, good signal an can possibly distinguish duplicate situations with the same number of data points present
DZ.IMAC.MBR.ON_expanded_log2_norm.class1.PS_ssPTMSEA_3 <- c()
for (i in unique.15aa.windows){
  print(i)

  #i <- "RERSRTGSESSQTGA"; i #Eif4b_S422x
  
  temp <-DZ.IMAC.MBR.ON_expanded_log2_norm.class1.PS_ssPTMSEA_2 %>% filter(aa15.window %in% c(i))
  
  temp$temp.sum <- apply(temp[, c("Norm.Intensity.A", "Norm.Intensity.F", "Norm.Intensity.K", "Norm.Intensity.P", 
                                  "Norm.Intensity.B", "Norm.Intensity.G", "Norm.Intensity.L", "Norm.Intensity.Q",    
                                  "Norm.Intensity.C", "Norm.Intensity.H", "Norm.Intensity.M", "Norm.Intensity.R", 
                                  "Norm.Intensity.D", "Norm.Intensity.I", "Norm.Intensity.N", "Norm.Intensity.S",    
                                  "Norm.Intensity.E", "Norm.Intensity.J", "Norm.Intensity.O")], 1, function(x) sum(x, na.rm=T)) 
  
  temp <- temp %>% arrange(desc(temp.sum)) 
  temp <- temp[1,]
 
  DZ.IMAC.MBR.ON_expanded_log2_norm.class1.PS_ssPTMSEA_3 <- bind_rows(DZ.IMAC.MBR.ON_expanded_log2_norm.class1.PS_ssPTMSEA_3 , temp)
}

glimpse(DZ.IMAC.MBR.ON_expanded_log2_norm.class1.PS_ssPTMSEA_3)
table(duplicated(DZ.IMAC.MBR.ON_expanded_log2_norm.class1.PS_ssPTMSEA_3$aa15.window))

### add_p to end of 15aa window
DZ.IMAC.MBR.ON_expanded_log2_norm.class1.PS_ssPTMSEA_3$aa15.window.2 <- paste0(DZ.IMAC.MBR.ON_expanded_log2_norm.class1.PS_ssPTMSEA_3$aa15.window, "-p")
glimpse(DZ.IMAC.MBR.ON_expanded_log2_norm.class1.PS_ssPTMSEA_3) 


### prepare input for PTMSEA via https://cloud.genepattern.org
gene.pattern.PTM.SEA.input <- DZ.IMAC.MBR.ON_expanded_log2_norm.class1.PS_ssPTMSEA_3 %>% dplyr::select(aa15.window.2, 
                                                                                                       "Norm.Intensity.A", "Norm.Intensity.F", "Norm.Intensity.K", "Norm.Intensity.P", 
                                                                                                       "Norm.Intensity.B", "Norm.Intensity.G", "Norm.Intensity.L", "Norm.Intensity.Q",    
                                                                                                       "Norm.Intensity.C", "Norm.Intensity.H", "Norm.Intensity.M", "Norm.Intensity.R", 
                                                                                                       "Norm.Intensity.D", "Norm.Intensity.I", "Norm.Intensity.N", "Norm.Intensity.S",    
                                                                                                       "Norm.Intensity.E", "Norm.Intensity.J", "Norm.Intensity.O")
glimpse(gene.pattern.PTM.SEA.input)
nrow(gene.pattern.PTM.SEA.input) #13638

#impute a zero (e.g. for black and white regulation situations) 
gene.pattern.PTM.SEA.input[is.na(gene.pattern.PTM.SEA.input)] <- 0
glimpse(gene.pattern.PTM.SEA.input)

#prepare gct input file with cmapR package
gene.pattern.PTM.SEA.input_GCT <- as.matrix(gene.pattern.PTM.SEA.input[,c( "Norm.Intensity.A", "Norm.Intensity.F", "Norm.Intensity.K", "Norm.Intensity.P", 
                                                                           "Norm.Intensity.B", "Norm.Intensity.G", "Norm.Intensity.L", "Norm.Intensity.Q",    
                                                                           "Norm.Intensity.C", "Norm.Intensity.H", "Norm.Intensity.M", "Norm.Intensity.R", 
                                                                           "Norm.Intensity.D", "Norm.Intensity.I", "Norm.Intensity.N", "Norm.Intensity.S",    
                                                                           "Norm.Intensity.E", "Norm.Intensity.J", "Norm.Intensity.O")])
glimpse(gene.pattern.PTM.SEA.input_GCT)
colnames(gene.pattern.PTM.SEA.input_GCT) <- c( "Norm.Intensity.A", "Norm.Intensity.F", "Norm.Intensity.K", "Norm.Intensity.P", 
                                               "Norm.Intensity.B", "Norm.Intensity.G", "Norm.Intensity.L", "Norm.Intensity.Q",    
                                               "Norm.Intensity.C", "Norm.Intensity.H", "Norm.Intensity.M", "Norm.Intensity.R", 
                                               "Norm.Intensity.D", "Norm.Intensity.I", "Norm.Intensity.N", "Norm.Intensity.S",    
                                               "Norm.Intensity.E", "Norm.Intensity.J", "Norm.Intensity.O")
rownames(gene.pattern.PTM.SEA.input_GCT) <- gene.pattern.PTM.SEA.input$aa15.window.2
glimpse(gene.pattern.PTM.SEA.input_GCT)

gene.pattern.PTM.SEA.input_GCT.2 <- new("GCT", mat=gene.pattern.PTM.SEA.input_GCT)
glimpse(gene.pattern.PTM.SEA.input_GCT.2)

#save GCT file for PTMSEA input
write_gct(gene.pattern.PTM.SEA.input_GCT.2, "cells.IMAC.mouse.windows.zero.imputed.ssPTMSEA.input.gct")


### run PTMSEA via Gene Pattern at https://cloud.genepattern.org ################################################################
# pathwaydb: ptm.sig.db.all.flanking.mouse.v1.9.0.gmt see https://github.com/broadinstitute/ssGSEA2.0
# all settings default
# save resulsta and ...combined.gct file
# change ...combined.gct file type to .txt and remove first two rows


### analyze PTMSEA output from GenePattern ####################################################################################################################################

## load PTMSEA result cells IMAC mouse windows
#NES is named Norm.Intensity.[...] here (just like the input columns)
ss.PTMSEA.result <- fread("cells.IMAC.mouse.windows.zero.imputed.ssPTMSEA.combined.result.txt")

glimpse(ss.PTMSEA.result)
nrow(ss.PTMSEA.result) #37

#sample overview
# A	1	GFP
# F	2	GFP
# K	3	GFP
# P	4	GFP
# B	5	Fgfr2-FL
# G	6	Fgfr2-FL
# L	7	Fgfr2-FL
# Q	8	Fgfr2-FL
# C	9	  Fgfr2-dE18
# H	10	Fgfr2-dE18
# M	11	Fgfr2-dE18
# R	12	Fgfr2-dE18
# D	13	Fgfr2-FL-Bicc1
# I	14	Fgfr2-FL-Bicc1
# N	15	Fgfr2-FL-Bicc1
# S	16	Fgfr2-FL-Bicc1
# E	17	Fgfr2-dE18-Bicc1
# J	18	Fgfr2-dE18-Bicc1
# O	19	Fgfr2-dE18-Bicc1
# T	20	Fgfr2-dE18-Bicc1

#calculate some measures
for(i in 1:nrow(ss.PTMSEA.result)){
  print(i)
  
  #GFP vs dE18 t.test
  tempA <- t.test(as.numeric(ss.PTMSEA.result[i, c("Norm.Intensity.A", "Norm.Intensity.F", "Norm.Intensity.K", "Norm.Intensity.P")]),
                  as.numeric(ss.PTMSEA.result[i, c( "Norm.Intensity.C", "Norm.Intensity.H", "Norm.Intensity.M", "Norm.Intensity.R")]))
  ss.PTMSEA.result$GFP_dE18_t.test.pval[i] <- tempA$p.value
  
  
  #FL vs dE18 t.test
  tempB <- t.test(as.numeric(ss.PTMSEA.result[i, c("Norm.Intensity.B", "Norm.Intensity.G", "Norm.Intensity.L", "Norm.Intensity.Q")]),
                  as.numeric(ss.PTMSEA.result[i, c( "Norm.Intensity.C", "Norm.Intensity.H", "Norm.Intensity.M", "Norm.Intensity.R")]))
  ss.PTMSEA.result$FL_dE18_t.test.pval[i] <- tempB$p.value
  
  #GFP vs FL
  tempC <- t.test(as.numeric(ss.PTMSEA.result[i, c("Norm.Intensity.A", "Norm.Intensity.F", "Norm.Intensity.K", "Norm.Intensity.P")]),
                  as.numeric(ss.PTMSEA.result[i, c("Norm.Intensity.B", "Norm.Intensity.G", "Norm.Intensity.L", "Norm.Intensity.Q")]))
  ss.PTMSEA.result$GFP_FL_t.test.pval[i] <- tempC$p.value
  
  #add FC FL vs dE18
  ss.PTMSEA.result$FC_FL_dE18[i] <- foldchange(mean(as.numeric(ss.PTMSEA.result[i, c("Norm.Intensity.C", "Norm.Intensity.H", "Norm.Intensity.M", "Norm.Intensity.R")])), mean(as.numeric(ss.PTMSEA.result[i, c("Norm.Intensity.B", "Norm.Intensity.G", "Norm.Intensity.L", "Norm.Intensity.Q")])))
}

glimpse(ss.PTMSEA.result)

sign.ss.PTMSEA.result <- ss.PTMSEA.result %>% filter(GFP_dE18_t.test.pval < 0.05 | FL_dE18_t.test.pval < 0.05 | GFP_FL_t.test.pval < 0.05)
nrow(sign.ss.PTMSEA.result) #22

sign.ss.PTMSEA.result <- sign.ss.PTMSEA.result %>% arrange(FC_FL_dE18)
glimpse((sign.ss.PTMSEA.result))

order.plot.18.01.22 <- sign.ss.PTMSEA.result$id
order.plot.18.01.22 

custom.order.mouse.ssPTMSEA.plot.DZ.170122 <- c("PERT-PSP_SERUM",            "KINASE-PSP_ERK2/Mapk1",                                  "KINASE-PSP_JNK1/Mapk8",     "KINASE-PSP_Akt1",          
                                                 "KINASE-PSP_mTOR/Mtor",      "KINASE-PSP_p70S6K/Rps6kb1", "KINASE-PSP_MSK1/Rps6ka5",   "PERT-PSP_RAPAMYCIN",        
                                                 "PERT-PSP_EGF",                                                                        "PERT-PSP_INSULIN",          "PERT-PSP_LIF",             
                                                 "PERT-PSP_PDGF",                                                                       "PERT-PSP_UV",               "KINASE-PSP_AurB/Aurkb",    
                                                 "PERT-PSP_H2O2",             "KINASE-PSP_PDK1/Pdpk1",      
                                                 "KINASE-PSP_AMPKA1/Prkaa1",  "KINASE-PSP_CDK1/Cdk1",                                   
                                                 "KINASE-PSP_PKCA/Prkca",     "KINASE-PSP_PKACA/Prkaca",                                "PERT-PSP_BOMBESIN",        
                                                 "PERT-PSP_LPS")   
custom.order.mouse.ssPTMSEA.plot.DZ.170122


### z-score result

#start define functions
zscorescalebaseR<- function(x){
  #arg <- c(1, 2, 3, 4, 5)
  #temp <- (x-mean(x))/sd(x)
  temp <- (x-mean(x))/sd(x)
  return(temp)
}
x <- c(1, 2, 3, 4, 5); x
zscorescalebaseR(x = x)
mean(zscorescalebaseR(x = x))
#stop define functions

D <- sign.ss.PTMSEA.result %>% select(Norm.Intensity.A : Norm.Intensity.O)
D <- t(apply(D,1, function(x) zscorescalebaseR(x = x)))

sd(D[1,]) #should be 1

DD <- as_tibble(D)

sign.ss.PTMSEA.result <- bind_cols(sign.ss.PTMSEA.result %>% select(id), DD) #with z score

#reshape
#NES is named with Norm.Intensity.[...] here (the input column names)
melted.ss.PTMSEA.result.18.01.22 <- reshape2::melt(sign.ss.PTMSEA.result %>% select("id",  "Norm.Intensity.A", "Norm.Intensity.F", "Norm.Intensity.K", "Norm.Intensity.P", 
                                                                                    "Norm.Intensity.B", "Norm.Intensity.G", "Norm.Intensity.L", "Norm.Intensity.Q",    
                                                                                    "Norm.Intensity.C", "Norm.Intensity.H", "Norm.Intensity.M", "Norm.Intensity.R", 
                                                                                    "Norm.Intensity.D", "Norm.Intensity.I", "Norm.Intensity.N", "Norm.Intensity.S",    
                                                                                    "Norm.Intensity.E", "Norm.Intensity.J", "Norm.Intensity.O"))
glimpse(melted.ss.PTMSEA.result.18.01.22)

#plot
ggplot(melted.ss.PTMSEA.result.18.01.22, aes(x=variable, y=id, fill=value)) +  
  geom_tile() + 
  #scale_fill_gradient2(low = "blue", mid = "white",high = "red") + #color A
  #scale_fill_gradient2(low = "blue", mid = "white",high = "darkred") + #color B
  #scale_fill_scico(palette = 'vikO',name="NES", na.value = "grey50")+ #color C 
  scale_fill_paletteer_c("pals::ocean.balance")+ #color D
  scale_y_discrete(limits=  rev(custom.order.mouse.ssPTMSEA.plot.DZ.170122))+ 
  theme(legend.position="bottom", legend.justification = "center")+
  theme(axis.line = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  theme(legend.title.align=0.5)+
  labs(fill = "row z-scored\nNES")+ #legend title
  ggtitle("ssPTM-SEA  \n mouse windows \n cells IMAC") +
  theme(plot.title = element_text(hjust = 0.5, face="bold", size=16)) + 
  theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1))+ 
  theme(axis.text.y= element_text(size=14))+  
  xlab(NULL) + 
  ylab(NULL) +
  theme(axis.text.y=element_text(face= "plain", colour="black", size=8) )+
  theme(axis.text.x = element_blank()) +
  theme(axis.ticks.x = element_blank())+
  
  geom_vline(xintercept = c(4.5, 8.5, 12.5, 16.5) , size=1.0, linetype = "solid", color="white")+
  annotate("text", x = 2.5, y = nrow(sign.ss.PTMSEA.result)+2, label = "GFP", fontface = "bold", size=3)+
  annotate("text", x = 6.5, y = nrow(sign.ss.PTMSEA.result)+2, label = "FL", fontface = "bold", size=3)+
  annotate("text", x = 10.5, y = nrow(sign.ss.PTMSEA.result)+2, label = "dE18", fontface = "bold", size=3)+
  annotate("text", x = 14.5, y = nrow(sign.ss.PTMSEA.result)+2, label = "FL\nBicc1", fontface = "bold", size=3)+
  annotate("text", x = 18.0, y = nrow(sign.ss.PTMSEA.result)+2, label = "dE18\nBicc1", fontface = "bold", size=3)+
  annotate("text", x = 3.5, y = nrow(sign.ss.PTMSEA.result)+3.5, label = "",  color = "transparent")

#ggsave("ss.PTMSEA.cells.IMAC.mouse.windows.custom.order.pdf", useDingbats=FALSE,  width = 10, height =15, units = "cm") 


################################################################################################################################
################################################################################################################################
################################################################################################################################
################################################################################################################################
################################################################################################################################

### cells IMAC: select and heatmap visualize candidate phosphosites and their upstream kinases 
### perform several different two group analyses and subject results to RoKAI tool (https://rokai.io/)
### combine results and visualize selection

## perform separate two group analyses
# GFP vs. FL
# GFP vs. dE18
# GFP vs. FL.Bicc1
# GFP vs. dE18.Bicc1
# FL vs. dE18
# FL vs. FL.Bicc1
# FL vs dE18.Bicc1

## sample overview
# A	1	GFP
# F	2	GFP
# K	3	GFP
# P	4	GFP

# B	5	Fgfr2-FL
# G	6	Fgfr2-FL
# L	7	Fgfr2-FL
# Q	8	Fgfr2-FL

# C	9	Fgfr2-dE18
# H	10	Fgfr2-dE18
# M	11	Fgfr2-dE18
# R	12	Fgfr2-dE18

# D	13	Fgfr2-FL-Bicc1
# I	14	Fgfr2-FL-Bicc1
# N	15	Fgfr2-FL-Bicc1
# S	16	Fgfr2-FL-Bicc1

# E	17	Fgfr2-dE18-Bicc1
# J	18	Fgfr2-dE18-Bicc1
# O	19	Fgfr2-dE18-Bicc1
# T	20	Fgfr2-dE18-Bicc1 (excluded)

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



#load("DZ.IMAC.MBR.ON_expanded_log2_norm.class1.PS.Rdata")
#glimpse(DZ.IMAC.MBR.ON_expanded_log2_norm.class1.PS) #16430


### GFP vs FL & RoKAI ##########################################################################################################################################
### GFP vs FL & RoKAI ##########################################################################################################################################

## select columns of interest
GFP_FL <- dplyr::select(DZ.IMAC.MBR.ON_expanded_log2_norm.class1.PS,
                        "ID.FR.all.C1.PS",
                        "id",                          "Proteins",                    "Positions.within.proteins",   "Leading.proteins",            "Protein",                    
                        "Protein.names",               "Gene.names",                  "Fasta.headers",               "Localization.prob",           "Number.of.Phospho.(STY)",    
                        "Amino.acid",                  "Sequence.window",             "Phospho.(STY).Probabilities", "Position.in.peptide",         "Reverse",                    
                        "Potential.contaminant",       "Positions",                   "Position",                    "Peptide.IDs",                 "Mod..peptide.IDs",   
                        "PS_Multiplicity", "Prot.Name_AApos", "Gene.Name_AApos",
                        
                        "Norm.Intensity.A", "Norm.Intensity.F", "Norm.Intensity.K", "Norm.Intensity.P",
                        "Norm.Intensity.B", "Norm.Intensity.G", "Norm.Intensity.L", "Norm.Intensity.Q"
                        
)
glimpse(GFP_FL)
colnames(GFP_FL)

##check data completeness
GFP_FL$data.presence.all <- apply(GFP_FL[, c("Norm.Intensity.A", "Norm.Intensity.F", "Norm.Intensity.K", "Norm.Intensity.P",    
                                             "Norm.Intensity.B", "Norm.Intensity.G", "Norm.Intensity.L", "Norm.Intensity.Q")], 1 , function(x) sum(!is.na(x))) 
glimpse(GFP_FL) #16430

##remove complete zero rows
table(GFP_FL$data.presence.all == 0) #F14254 T 2176

GFP_FL <- GFP_FL %>% filter(data.presence.all > 0)
glimpse(GFP_FL) #14254

##X out of Y filter 
#count data presence and filter
GFP_FL$data.presence.left <- apply(GFP_FL[, c("Norm.Intensity.A", "Norm.Intensity.F", "Norm.Intensity.K", "Norm.Intensity.P")],1, function(x) sum(!is.na(x)))
GFP_FL$data.presence.right <- apply(GFP_FL[, c( "Norm.Intensity.B", "Norm.Intensity.G", "Norm.Intensity.L", "Norm.Intensity.Q")],1, function(x) sum(!is.na(x)))
glimpse(GFP_FL)

temp.left <- filter(GFP_FL, data.presence.left >= 3) ; nrow(temp.left)    #8462
temp.right <- filter(GFP_FL, data.presence.right >= 3) ; nrow(temp.right) #8544

temp.left.and.right <- bind_rows(temp.left , temp.right) %>% distinct() ; nrow(temp.left.and.right) #9714

GFP_FL.XYfilter <- temp.left.and.right
nrow(GFP_FL.XYfilter) #9714
length(unique(GFP_FL.XYfilter$Prot.Name_AApos)) #8572

#add FRID
GFP_FL.XYfilter$FRID <- 1:nrow(GFP_FL.XYfilter) #here FRID is also row number
glimpse(GFP_FL.XYfilter)

#find likely Limma NA pval situations 
possible.NApvals.A <- GFP_FL.XYfilter %>% filter((data.presence.left >= 4 & data.presence.right == 0) | (data.presence.left == 0 & data.presence.right >= 4)  )
glimpse(possible.NApvals.A)
length(unique(possible.NApvals.A$FRID))
FRIDS.possible.NApvals <- possible.NApvals.A$FRID
FRIDS.possible.NApvals 

# for these filtered  black and white situations add zero so that p.val can be calculated
W <- GFP_FL.XYfilter[FRIDS.possible.NApvals, ]
W <- W %>% mutate_at(vars("Norm.Intensity.A", "Norm.Intensity.F", "Norm.Intensity.K", "Norm.Intensity.P",    
                          "Norm.Intensity.B", "Norm.Intensity.G", "Norm.Intensity.L", "Norm.Intensity.Q"), ~replace(., is.na(.), 0))
glimpse(W)

#combine again
GFP_FL.XYfilter.2 <- GFP_FL.XYfilter %>% filter(! FRID %in% FRIDS.possible.NApvals ) 
GFP_FL.XYfilter.2 <- bind_rows(GFP_FL.XYfilter.2, W)
GFP_FL.XYfilter.2 <- GFP_FL.XYfilter.2 %>% arrange(FRID)
glimpse(GFP_FL.XYfilter.2)


### LIMMA
# get data
DataforLimma <- GFP_FL.XYfilter.2[, c("Norm.Intensity.A", "Norm.Intensity.F", "Norm.Intensity.K", "Norm.Intensity.P",    
                                      "Norm.Intensity.B", "Norm.Intensity.G", "Norm.Intensity.L", "Norm.Intensity.Q")]
glimpse(DataforLimma)
colnames(DataforLimma)

design <- cbind(left=c(rep(1, 4), rep(0, 4)), right=c(rep(0, 4), rep(1, 4))); design

fit <- lmFit(DataforLimma, design) # Warning: Partial NA coefficients for 255 probe(s) 
cont.matrix <- makeContrasts( left.vs.right =  right - left, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

Limma.p.value <- as.numeric(fit2$p.value)
Limma.adj.p.valueBH <- p.adjust(Limma.p.value, "BH")
Limma.t <- as.numeric(fit2$t)

GFP_FL.XYfilter.2$Limma.p.value <- as.numeric(Limma.p.value)
GFP_FL.XYfilter.2$Limma.adj.p.valueBH <- as.numeric(Limma.adj.p.valueBH)
GFP_FL.XYfilter.2$Limma.t <- as.numeric(Limma.t)

GFP_FL.XYfilter.B <- GFP_FL.XYfilter.2
glimpse(GFP_FL.XYfilter.B)


### calculate FC
for(i in 1:nrow(GFP_FL.XYfilter.B)){
  print(i)
  left.vals <- as.numeric(GFP_FL.XYfilter.B[i,c("Norm.Intensity.A", "Norm.Intensity.F", "Norm.Intensity.K", "Norm.Intensity.P")])
  right.vals <- as.numeric(GFP_FL.XYfilter.B[i,c("Norm.Intensity.B", "Norm.Intensity.G", "Norm.Intensity.L", "Norm.Intensity.Q")])
  FC <- FR_FC_LOGData_2.3.17(control = left.vals, treatment = right.vals, LOGtype = 2); FC
  GFP_FL.XYfilter.B$FC_FR[i] <- FC
}
glimpse(GFP_FL.XYfilter.B)

#check numbers unique Sequence.window
GFP_FL.XYfilter.B %>% filter(Limma.p.value < 0.05) %>% select(Sequence.window) %>% distinct() %>% nrow() #1864 total sign
GFP_FL.XYfilter.B %>% filter(Limma.p.value < 0.05) %>% filter(FC_FR > 1) %>% select(Sequence.window) %>% distinct() %>% nrow() #800 up
GFP_FL.XYfilter.B %>% filter(Limma.p.value < 0.05) %>% filter(FC_FR < 1) %>% select(Sequence.window) %>% distinct() %>% nrow() #1068 down


# filter na pvalues
GFP_FL.XYfilter.B <- GFP_FL.XYfilter.B %>% filter(!is.na(Limma.p.value))
glimpse(GFP_FL.XYfilter.B) #9459

# check +/- 10000 FCs
table(GFP_FL.XYfilter.B$FC_FR >= 10000) #102
table(GFP_FL.XYfilter.B$FC_FR <= -10000) #58

#find second min or max and add +/-5 for +/- 10000 FCs
temp <- GFP_FL.XYfilter.B %>% filter(!FC_FR %in% c(10000, -10000)) %>% arrange(desc(FC_FR))
head(temp$FC_FR)
tail(temp$FC_FR)

GFP_FL.XYfilter.B$FC_FR <- replace(GFP_FL.XYfilter.B$FC_FR, GFP_FL.XYfilter.B$FC_FR==10000, max(temp$FC_FR)+5)
GFP_FL.XYfilter.B$FC_FR <- replace(GFP_FL.XYfilter.B$FC_FR, GFP_FL.XYfilter.B$FC_FR==-10000, min(temp$FC_FR)-5)
glimpse(GFP_FL.XYfilter.B) #9459

### GFP_FL reshape for RoKAI
GFP_FL.RoKAI <- GFP_FL.XYfilter.B 

# make aa15.window
# take first entry sequence window if there are two entries
GFP_FL.RoKAI$aa15.window <- pbsapply(GFP_FL.RoKAI$Sequence.window, function(x) unlist(str_split(string =x, pattern = ";"))[1]  )
GFP_FL.RoKAI$aa15.window <- pbsapply(GFP_FL.RoKAI$aa15.window, function(x) make.aa.15.window.from.31.sequence.window(input = x)  )

glimpse(GFP_FL.RoKAI)

table(duplicated(GFP_FL.RoKAI$aa15.window))
which(duplicated(GFP_FL.RoKAI$aa15.window))

### if there are duplicate entries (e.g. because of _1_2_3) take the entry with the lowest p value
unique.15aa.windows <- unique(GFP_FL.RoKAI$aa15.window)
length(unique.15aa.windows) #8333
glimpse(GFP_FL.RoKAI) #9459


GFP_FL.RoKAI.2 <- c()
for (i in unique.15aa.windows){
  print(i)
  temp <- GFP_FL.RoKAI %>% filter(aa15.window %in% c(i))
  temp <- temp %>% arrange(Limma.p.value)
  temp <- temp[1,]
  GFP_FL.RoKAI.2 <- bind_rows(GFP_FL.RoKAI.2, temp)
}

glimpse(GFP_FL.RoKAI.2)
table(duplicated(GFP_FL.RoKAI.2$aa15.window))

#calculate log2(FC) and order
GFP_FL.RoKAI.2$log2FC <- log2(abs(GFP_FL.RoKAI.2$FC_FR ))
GFP_FL.RoKAI.2$sign.FC_FR <- sign(GFP_FL.RoKAI.2$FC_FR)
GFP_FL.RoKAI.2$log2FC <- sign(GFP_FL.RoKAI.2$FC_FR) * log2(abs(GFP_FL.RoKAI.2$FC_FR ))
glimpse(GFP_FL.RoKAI.2) 

GFP_FL.RoKAI.2 <- GFP_FL.RoKAI.2 %>% arrange(desc(log2FC))
head(GFP_FL.RoKAI.2$FC_FR)
tail(GFP_FL.RoKAI.2$FC_FR)

### prepare RoKAI input and save data
GFP_FL.RoKAI.2.23.11.21 <- GFP_FL.RoKAI.2 %>% select(Protein, Position, log2FC) 
table(is.na(GFP_FL.RoKAI.2.23.11.21$Protein))
GFP_FL.RoKAI.2.23.11.21 <- GFP_FL.RoKAI.2.23.11.21 %>% filter(!is.na(Protein))
glimpse(GFP_FL.RoKAI.2.23.11.21)
table(is.na(GFP_FL.RoKAI.2.23.11.21$Protein))
table(is.na(GFP_FL.RoKAI.2.23.11.21$Position))
table(is.na(GFP_FL.RoKAI.2.23.11.21$log2FC))
colnames(GFP_FL.RoKAI.2.23.11.21) <- c("Protein", "Position", "Quantification")


#write_csv(GFP_FL.RoKAI.2.23.11.21,"GFP_FL.RoKAI.2_23.11.21.csv", col_names=TRUE)
#save(GFP_FL.RoKAI, file="GFP_FL.RoKAI.Rdata")

# perform analysis via https://rokai.io/
# tool settings was default: Reference Proteome: Uniprot Mouse / Fold changes: normalized / Kinase substrate dataset: PSP + Signor / RoKAI network: KS+PPI+SD+CoEv
# save RoKAI output kinase_table.csv and kinase_targets.csv and load these later




### GFP vs dE18 & RoKAI ##########################################################################################################################################
### GFP vs dE18 & RoKAI ##########################################################################################################################################
glimpse(DZ.IMAC.MBR.ON_expanded_log2_norm.class1.PS)

## select columns of interest
GFP_dE18 <- dplyr::select(DZ.IMAC.MBR.ON_expanded_log2_norm.class1.PS,
                          "ID.FR.all.C1.PS",
                          "id",                          "Proteins",                    "Positions.within.proteins",   "Leading.proteins",            "Protein",                    
                          "Protein.names",               "Gene.names",                  "Fasta.headers",               "Localization.prob",           "Number.of.Phospho.(STY)",    
                          "Amino.acid",                  "Sequence.window",             "Phospho.(STY).Probabilities", "Position.in.peptide",         "Reverse",                    
                          "Potential.contaminant",       "Positions",                   "Position",                    "Peptide.IDs",                 "Mod..peptide.IDs",   
                          "PS_Multiplicity", "Prot.Name_AApos", "Gene.Name_AApos",
                          
                          "Norm.Intensity.A", "Norm.Intensity.F", "Norm.Intensity.K", "Norm.Intensity.P",
                          "Norm.Intensity.C", "Norm.Intensity.H", "Norm.Intensity.M", "Norm.Intensity.R"
                          
)
glimpse(GFP_dE18)
colnames(GFP_dE18)

##check data completeness
GFP_dE18$data.presence.all <- apply(GFP_dE18[, c("Norm.Intensity.A", "Norm.Intensity.F", "Norm.Intensity.K", "Norm.Intensity.P",    
                                                 "Norm.Intensity.C", "Norm.Intensity.H", "Norm.Intensity.M", "Norm.Intensity.R")], 1 , function(x) sum(!is.na(x))) #count data completeness
glimpse(GFP_dE18) #16430


#remove complete zero rows
table(GFP_dE18$data.presence.all == 0) #F14914 T1516

GFP_dE18 <- GFP_dE18 %>% filter(data.presence.all > 0)
glimpse(GFP_dE18) #14914

## X out of Y filter
#count data presence and filter
GFP_dE18$data.presence.left <- apply(GFP_dE18[, c("Norm.Intensity.A", "Norm.Intensity.F", "Norm.Intensity.K", "Norm.Intensity.P")],1, function(x) sum(!is.na(x)))
GFP_dE18$data.presence.right <- apply(GFP_dE18[, c( "Norm.Intensity.C", "Norm.Intensity.H", "Norm.Intensity.M", "Norm.Intensity.R")],1, function(x) sum(!is.na(x)))
glimpse(GFP_dE18)

temp.left <- filter(GFP_dE18, data.presence.left >= 3) ; nrow(temp.left)    #8462
temp.right <- filter(GFP_dE18, data.presence.right >= 3) ; nrow(temp.right) #8438

temp.left.and.right <- bind_rows(temp.left , temp.right) ; nrow(temp.left.and.right)
temp.left.and.right <- distinct(temp.left.and.right) ; nrow(temp.left.and.right) #9974

GFP_dE18.XYfilter <- temp.left.and.right
nrow(GFP_dE18.XYfilter) #9974
length(unique(GFP_dE18.XYfilter$Prot.Name_AApos)) #8801

#add FRID
GFP_dE18.XYfilter$FRID <- 1:nrow(GFP_dE18.XYfilter) #here FRID is also row number
glimpse(GFP_dE18.XYfilter)

#find likely Limma NA pval situations 
possible.NApvals.A <- GFP_dE18.XYfilter %>% filter((data.presence.left >= 4 & data.presence.right == 0) | (data.presence.left == 0 & data.presence.right >= 4)  )
glimpse(possible.NApvals.A)
length(unique(possible.NApvals.A$FRID))
FRIDS.possible.NApvals <- possible.NApvals.A$FRID
FRIDS.possible.NApvals 

# for these filtered  black and white situation add zero so that p.val can be calculated
W <- GFP_dE18.XYfilter[FRIDS.possible.NApvals, ]
glimpse(W)
W <- W %>% mutate_at(vars("Norm.Intensity.A", "Norm.Intensity.F", "Norm.Intensity.K", "Norm.Intensity.P",    
                          "Norm.Intensity.C", "Norm.Intensity.H", "Norm.Intensity.M", "Norm.Intensity.R"), ~replace(., is.na(.), 0))
glimpse(W)

#combine again
GFP_dE18.XYfilter.2 <- GFP_dE18.XYfilter %>% filter(! FRID %in% FRIDS.possible.NApvals ) 
GFP_dE18.XYfilter.2 <- bind_rows(GFP_dE18.XYfilter.2, W)
GFP_dE18.XYfilter.2 <- GFP_dE18.XYfilter.2 %>% arrange(FRID)
glimpse(GFP_dE18.XYfilter.2)


### LIMMA
# get data
DataforLimma <- GFP_dE18.XYfilter.2[, c("Norm.Intensity.A", "Norm.Intensity.F", "Norm.Intensity.K", "Norm.Intensity.P",    
                                        "Norm.Intensity.C", "Norm.Intensity.H", "Norm.Intensity.M", "Norm.Intensity.R")]
glimpse(DataforLimma)
colnames(DataforLimma)

design <- cbind(left=c(rep(1, 4), rep(0, 4)), right=c(rep(0, 4), rep(1, 4))); design

fit <- lmFit(DataforLimma, design) # Warning: Partial NA coefficients for 481 probe(s) 
cont.matrix <- makeContrasts( left.vs.right =  right - left, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

Limma.p.value <- as.numeric(fit2$p.value)
glimpse(Limma.p.value )
Limma.adj.p.valueBH <- p.adjust(Limma.p.value, "BH") #
glimpse(Limma.adj.p.valueBH)
Limma.t <- as.numeric(fit2$t)
glimpse(Limma.t)

GFP_dE18.XYfilter.2$Limma.p.value <- as.numeric(Limma.p.value)
GFP_dE18.XYfilter.2$Limma.adj.p.valueBH <- as.numeric(Limma.adj.p.valueBH)
GFP_dE18.XYfilter.2$Limma.t <- as.numeric(Limma.t)

GFP_dE18.XYfilter.B <- GFP_dE18.XYfilter.2
glimpse(GFP_dE18.XYfilter.B)


### calculate FC
for(i in 1:nrow(GFP_dE18.XYfilter.B)){
  print(i)
  left.vals <- as.numeric(GFP_dE18.XYfilter.B[i,c("Norm.Intensity.A", "Norm.Intensity.F", "Norm.Intensity.K", "Norm.Intensity.P")])
  right.vals <- as.numeric(GFP_dE18.XYfilter.B[i,c("Norm.Intensity.C", "Norm.Intensity.H", "Norm.Intensity.M", "Norm.Intensity.R")]) #
  FC <- FR_FC_LOGData_2.3.17(control = left.vals, treatment = right.vals, LOGtype = 2); FC
  GFP_dE18.XYfilter.B$FC_FR[i] <- FC
}
glimpse(GFP_dE18.XYfilter.B)


#check numbers unique Sequence.window
GFP_dE18.XYfilter.B %>% filter(Limma.p.value < 0.05) %>% select(Sequence.window) %>% distinct() %>% nrow() #3304 total sign
GFP_dE18.XYfilter.B %>% filter(Limma.p.value < 0.05) %>% filter(FC_FR > 1) %>% select(Sequence.window) %>% distinct() %>% nrow() #1480 up
GFP_dE18.XYfilter.B %>% filter(Limma.p.value < 0.05) %>% filter(FC_FR < 1) %>% select(Sequence.window) %>% distinct() %>% nrow() #1855 down

# filter na pvalues
table(is.na(GFP_dE18.XYfilter.B$Limma.p.value))
GFP_dE18.XYfilter.B <- GFP_dE18.XYfilter.B %>% filter(!is.na(Limma.p.value))
glimpse(GFP_dE18.XYfilter.B) #9493
table(is.na(GFP_dE18.XYfilter.B$Limma.p.value))


# check +/- 10000 FCs
table(GFP_dE18.XYfilter.B$FC_FR >= 10000) #219
table(GFP_dE18.XYfilter.B$FC_FR <= -10000) #165

#find second min or max and add +/-5 for +/- 10000 FCs
temp <- GFP_dE18.XYfilter.B %>% filter(!FC_FR %in% c(10000, -10000)) %>% arrange(desc(FC_FR))
head(temp$FC_FR)
tail(temp$FC_FR)
max(temp$FC_FR)
min(temp$FC_FR)

GFP_dE18.XYfilter.B$FC_FR <- replace(GFP_dE18.XYfilter.B$FC_FR, GFP_dE18.XYfilter.B$FC_FR==10000, max(temp$FC_FR)+5)
GFP_dE18.XYfilter.B$FC_FR <- replace(GFP_dE18.XYfilter.B$FC_FR, GFP_dE18.XYfilter.B$FC_FR==-10000, min(temp$FC_FR)-5)
summary(GFP_dE18.XYfilter.B$FC_FR)
glimpse(GFP_dE18.XYfilter.B) #9493

### GFP_dE18 reshape for RoKAI
GFP_dE18.RoKAI <- GFP_dE18.XYfilter.B 

# make aa15.window
# take first entry sequence window if there are two entries
GFP_dE18.RoKAI$aa15.window <- pbsapply(GFP_dE18.RoKAI$Sequence.window, function(x) unlist(str_split(string =x, pattern = ";"))[1]  )
GFP_dE18.RoKAI$aa15.window <- pbsapply(GFP_dE18.RoKAI$aa15.window, function(x) make.aa.15.window.from.31.sequence.window(input = x)  )

glimpse(GFP_dE18.RoKAI)

table(duplicated(GFP_dE18.RoKAI$aa15.window))
which(duplicated(GFP_dE18.RoKAI$aa15.window))


## if there are duplicate entries (e.g. because of _1_2_3) take the entry with the best p value
unique.15aa.windows <- unique(GFP_dE18.RoKAI$aa15.window)
length(unique.15aa.windows) #8355
glimpse(GFP_dE18.RoKAI) #9493

GFP_dE18.RoKAI.2 <- c()
for (i in unique.15aa.windows){
  print(i)
  temp <- GFP_dE18.RoKAI %>% filter(aa15.window %in% c(i))
  temp <- temp %>% arrange(Limma.p.value)
  temp <- temp[1,]
  GFP_dE18.RoKAI.2 <- bind_rows(GFP_dE18.RoKAI.2, temp)
}

glimpse(GFP_dE18.RoKAI.2)
table(duplicated(GFP_dE18.RoKAI.2$aa15.window))


#calculate log2(FC) and order
GFP_dE18.RoKAI.2$log2FC <- log2(abs(GFP_dE18.RoKAI.2$FC_FR ))
GFP_dE18.RoKAI.2$sign.FC_FR <- sign(GFP_dE18.RoKAI.2$FC_FR)
GFP_dE18.RoKAI.2$log2FC <- sign(GFP_dE18.RoKAI.2$FC_FR) * log2(abs(GFP_dE18.RoKAI.2$FC_FR ))
glimpse(GFP_dE18.RoKAI.2) 

GFP_dE18.RoKAI.2 <- GFP_dE18.RoKAI.2 %>% arrange(desc(log2FC))
head(GFP_dE18.RoKAI.2$FC_FR)
tail(GFP_dE18.RoKAI.2$FC_FR)

### prepare RoKAI input and save data
glimpse(GFP_dE18.RoKAI.2) 
GFP_dE18.RoKAI.2.23.11.21 <- GFP_dE18.RoKAI.2 %>% select(Protein, Position, log2FC)
table(is.na(GFP_dE18.RoKAI.2.23.11.21$Protein))
GFP_dE18.RoKAI.2.23.11.21 <- GFP_dE18.RoKAI.2.23.11.21 %>% filter(!is.na(Protein))
glimpse(GFP_dE18.RoKAI.2.23.11.21)
table(is.na(GFP_dE18.RoKAI.2.23.11.21$Protein))
table(is.na(GFP_dE18.RoKAI.2.23.11.21$Position))
table(is.na(GFP_dE18.RoKAI.2.23.11.21$log2FC))
colnames(GFP_dE18.RoKAI.2.23.11.21) <- c("Protein", "Position", "Quantification")

#write_csv(GFP_dE18.RoKAI.2.23.11.21,"GFP_dE18.RoKAI.2_23.11.21.csv", col_names=TRUE)
#save(GFP_dE18.RoKAI, file="GFP_dE18.RoKAI.Rdata")

# perform analysis via https://rokai.io/
# tool settings was default: Reference Proteome: Uniprot Mouse / Fold changes: normalized / Kinase substrate dataset: PSP + Signor / RoKAI network: KS+PPI+SD+CoEv
# save RoKAI output kinase_table.csv and kinase_targets.csv and load these later


### GFP vs FL.Bicc1 & RoKAI ##########################################################################################################################################
### GFP vs FL.Bicc1 & RoKAI ##########################################################################################################################################
glimpse(DZ.IMAC.MBR.ON_expanded_log2_norm.class1.PS)

# select columns of interest
GFP_FL.Bicc1 <- dplyr::select(DZ.IMAC.MBR.ON_expanded_log2_norm.class1.PS,
                              "ID.FR.all.C1.PS",
                              "id",                          "Proteins",                    "Positions.within.proteins",   "Leading.proteins",            "Protein",                    
                              "Protein.names",               "Gene.names",                  "Fasta.headers",               "Localization.prob",           "Number.of.Phospho.(STY)",    
                              "Amino.acid",                  "Sequence.window",             "Phospho.(STY).Probabilities", "Position.in.peptide",         "Reverse",                    
                              "Potential.contaminant",       "Positions",                   "Position",                    "Peptide.IDs",                 "Mod..peptide.IDs",   
                              "PS_Multiplicity", "Prot.Name_AApos", "Gene.Name_AApos",
                              
                              "Norm.Intensity.A", "Norm.Intensity.F", "Norm.Intensity.K", "Norm.Intensity.P",
                              "Norm.Intensity.D", "Norm.Intensity.I", "Norm.Intensity.N", "Norm.Intensity.S"
                              
)
glimpse(GFP_FL.Bicc1)
colnames(GFP_FL.Bicc1)

GFP_FL.Bicc1$data.presence.all <- apply(GFP_FL.Bicc1[, c("Norm.Intensity.A", "Norm.Intensity.F", "Norm.Intensity.K", "Norm.Intensity.P",    
                                                         "Norm.Intensity.D", "Norm.Intensity.I", "Norm.Intensity.N", "Norm.Intensity.S")], 1 , function(x) sum(!is.na(x))) #count data completeness
glimpse(GFP_FL.Bicc1) #16430


#remove complete zero rows
table(GFP_FL.Bicc1$data.presence.all == 0) #F14584 T1846

GFP_FL.Bicc1 <- GFP_FL.Bicc1 %>% filter(data.presence.all > 0)
glimpse(GFP_FL.Bicc1) #14584

#X out of Y filter
#count data presence and filter
GFP_FL.Bicc1$data.presence.left <- apply(GFP_FL.Bicc1[, c("Norm.Intensity.A", "Norm.Intensity.F", "Norm.Intensity.K", "Norm.Intensity.P")],1, function(x) sum(!is.na(x)))
GFP_FL.Bicc1$data.presence.right <- apply(GFP_FL.Bicc1[, c( "Norm.Intensity.D", "Norm.Intensity.I", "Norm.Intensity.N", "Norm.Intensity.S")],1, function(x) sum(!is.na(x)))
glimpse(GFP_FL.Bicc1)


temp.left <- filter(GFP_FL.Bicc1, data.presence.left >= 3) ; nrow(temp.left)    #8462
temp.right <- filter(GFP_FL.Bicc1, data.presence.right >= 3) ; nrow(temp.right) #8866

temp.left.and.right <- bind_rows(temp.left , temp.right) ; nrow(temp.left.and.right)
temp.left.and.right <- distinct(temp.left.and.right) ; nrow(temp.left.and.right) #10059

GFP_FL.Bicc1.XYfilter <- temp.left.and.right
nrow(GFP_FL.Bicc1.XYfilter) #10059
length(unique(GFP_FL.Bicc1.XYfilter$Prot.Name_AApos)) #8853

glimpse(GFP_FL.Bicc1.XYfilter)

#add FRID
GFP_FL.Bicc1.XYfilter$FRID <- 1:nrow(GFP_FL.Bicc1.XYfilter) #here FRID is also row number
glimpse(GFP_FL.Bicc1.XYfilter)

#find likely Limma NA pval situations 
possible.NApvals.A <- GFP_FL.Bicc1.XYfilter %>% filter((data.presence.left >= 4 & data.presence.right == 0) | (data.presence.left == 0 & data.presence.right >= 4)  )
glimpse(possible.NApvals.A)
length(unique(possible.NApvals.A$FRID))
FRIDS.possible.NApvals <- possible.NApvals.A$FRID
FRIDS.possible.NApvals 

# for these filtered  black and white situation add zero so that p.val can be calculated
W <- GFP_FL.Bicc1.XYfilter[FRIDS.possible.NApvals, ]
glimpse(W)
W <- W %>% mutate_at(vars("Norm.Intensity.A", "Norm.Intensity.F", "Norm.Intensity.K", "Norm.Intensity.P",    
                          "Norm.Intensity.D", "Norm.Intensity.I", "Norm.Intensity.N", "Norm.Intensity.S"), ~replace(., is.na(.), 0))
glimpse(W)

#combine again
GFP_FL.Bicc1.XYfilter.2 <- GFP_FL.Bicc1.XYfilter %>% filter(! FRID %in% FRIDS.possible.NApvals ) 
GFP_FL.Bicc1.XYfilter.2 <- bind_rows(GFP_FL.Bicc1.XYfilter.2, W)
GFP_FL.Bicc1.XYfilter.2 <- GFP_FL.Bicc1.XYfilter.2 %>% arrange(FRID)
glimpse(GFP_FL.Bicc1.XYfilter.2)

### LIMMA
# get data
DataforLimma <- GFP_FL.Bicc1.XYfilter.2[, c("Norm.Intensity.A", "Norm.Intensity.F", "Norm.Intensity.K", "Norm.Intensity.P",    
                                            "Norm.Intensity.D", "Norm.Intensity.I", "Norm.Intensity.N", "Norm.Intensity.S")]
glimpse(DataforLimma)
colnames(DataforLimma)

design <- cbind(left=c(rep(1, 4), rep(0, 4)), right=c(rep(0, 4), rep(1, 4))); design

fit <- lmFit(DataforLimma, design) # Warning: Partial NA coefficients for 302 probe(s) 
cont.matrix <- makeContrasts( left.vs.right =  right - left, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

Limma.p.value <- as.numeric(fit2$p.value)
glimpse(Limma.p.value )
Limma.adj.p.valueBH <- p.adjust(Limma.p.value, "BH")
glimpse(Limma.adj.p.valueBH)
Limma.t <- as.numeric(fit2$t)
glimpse(Limma.t)

GFP_FL.Bicc1.XYfilter.2$Limma.p.value <- as.numeric(Limma.p.value)
GFP_FL.Bicc1.XYfilter.2$Limma.adj.p.valueBH <- as.numeric(Limma.adj.p.valueBH)
GFP_FL.Bicc1.XYfilter.2$Limma.t <- as.numeric(Limma.t)

GFP_FL.Bicc1.XYfilter.B <- GFP_FL.Bicc1.XYfilter.2
glimpse(GFP_FL.Bicc1.XYfilter.B)


### calculate FC
for(i in 1:nrow(GFP_FL.Bicc1.XYfilter.B)){
  print(i)
  left.vals <- as.numeric(GFP_FL.Bicc1.XYfilter.B[i,c("Norm.Intensity.A", "Norm.Intensity.F", "Norm.Intensity.K", "Norm.Intensity.P")])
  right.vals <- as.numeric(GFP_FL.Bicc1.XYfilter.B[i,c("Norm.Intensity.D", "Norm.Intensity.I", "Norm.Intensity.N", "Norm.Intensity.S")])
  FC <- FR_FC_LOGData_2.3.17(control = left.vals, treatment = right.vals, LOGtype = 2); FC
  GFP_FL.Bicc1.XYfilter.B$FC_FR[i] <- FC
}
glimpse(GFP_FL.Bicc1.XYfilter.B)


#check numbers unique Sequence.window
GFP_FL.Bicc1.XYfilter.B %>% filter(Limma.p.value < 0.05) %>% select(Sequence.window) %>% distinct() %>% nrow() #2424 total sign
GFP_FL.Bicc1.XYfilter.B %>% filter(Limma.p.value < 0.05) %>% filter(FC_FR > 1) %>% select(Sequence.window) %>% distinct() %>% nrow() #1349 up
GFP_FL.Bicc1.XYfilter.B %>% filter(Limma.p.value < 0.05) %>% filter(FC_FR < 1) %>% select(Sequence.window) %>% distinct() %>% nrow() #1083 down

# filter na pvalues
table(is.na(GFP_FL.Bicc1.XYfilter.B$Limma.p.value))
GFP_FL.Bicc1.XYfilter.B <- GFP_FL.Bicc1.XYfilter.B %>% filter(!is.na(Limma.p.value))
glimpse(GFP_FL.Bicc1.XYfilter.B) #9757
table(is.na(GFP_FL.Bicc1.XYfilter.B$Limma.p.value))


# check +/- 10000 FCs
table(GFP_FL.Bicc1.XYfilter.B$FC_FR >= 10000) #193
table(GFP_FL.Bicc1.XYfilter.B$FC_FR <= -10000) #46
summary(GFP_FL.Bicc1.XYfilter.B$FC_FR)

#find second min or max and add +/-5 for +/- 10000 FCs
temp <- GFP_FL.Bicc1.XYfilter.B %>% filter(!FC_FR %in% c(10000, -10000)) %>% arrange(desc(FC_FR))
head(temp$FC_FR)
tail(temp$FC_FR)
max(temp$FC_FR)
min(temp$FC_FR)

GFP_FL.Bicc1.XYfilter.B$FC_FR <- replace(GFP_FL.Bicc1.XYfilter.B$FC_FR, GFP_FL.Bicc1.XYfilter.B$FC_FR==10000, max(temp$FC_FR)+5)
GFP_FL.Bicc1.XYfilter.B$FC_FR <- replace(GFP_FL.Bicc1.XYfilter.B$FC_FR, GFP_FL.Bicc1.XYfilter.B$FC_FR==-10000, min(temp$FC_FR)-5)
summary(GFP_FL.Bicc1.XYfilter.B$FC_FR)
glimpse(GFP_FL.Bicc1.XYfilter.B) #9757

## GFP_FL.Bicc1 reshape for RoKAI
GFP_FL.Bicc1.RoKAI <- GFP_FL.Bicc1.XYfilter.B 

# make aa15.window
# take first entry sequence window if there are two entries
GFP_FL.Bicc1.RoKAI$aa15.window <- pbsapply(GFP_FL.Bicc1.RoKAI$Sequence.window, function(x) unlist(str_split(string =x, pattern = ";"))[1]  )
GFP_FL.Bicc1.RoKAI$aa15.window <- pbsapply(GFP_FL.Bicc1.RoKAI$aa15.window, function(x) make.aa.15.window.from.31.sequence.window(input = x)  )


glimpse(GFP_FL.Bicc1.RoKAI)

table(duplicated(GFP_FL.Bicc1.RoKAI$aa15.window))
which(duplicated(GFP_FL.Bicc1.RoKAI$aa15.window))


## if there are duplicate entries (e.g. because of _1_2_3) take the entry with the best p value
unique.15aa.windows <- unique(GFP_FL.Bicc1.RoKAI$aa15.window)
length(unique.15aa.windows) #8571
glimpse(GFP_FL.Bicc1.RoKAI) #9757


GFP_FL.Bicc1.RoKAI.2 <- c()
for (i in unique.15aa.windows){
  print(i)
  temp <- GFP_FL.Bicc1.RoKAI %>% filter(aa15.window %in% c(i))
  temp <- temp %>% arrange(Limma.p.value)
  temp <- temp[1,]
  GFP_FL.Bicc1.RoKAI.2 <- bind_rows(GFP_FL.Bicc1.RoKAI.2, temp)
}

glimpse(GFP_FL.Bicc1.RoKAI.2)
table(duplicated(GFP_FL.Bicc1.RoKAI.2$aa15.window))



#calculate log2(FC) and order
GFP_FL.Bicc1.RoKAI.2$log2FC <- log2(abs(GFP_FL.Bicc1.RoKAI.2$FC_FR ))
GFP_FL.Bicc1.RoKAI.2$sign.FC_FR <- sign(GFP_FL.Bicc1.RoKAI.2$FC_FR)
GFP_FL.Bicc1.RoKAI.2$log2FC <- sign(GFP_FL.Bicc1.RoKAI.2$FC_FR) * log2(abs(GFP_FL.Bicc1.RoKAI.2$FC_FR ))
glimpse(GFP_FL.Bicc1.RoKAI.2) 

GFP_FL.Bicc1.RoKAI.2 <- GFP_FL.Bicc1.RoKAI.2 %>% arrange(desc(log2FC))
head(GFP_FL.Bicc1.RoKAI.2$FC_FR)
tail(GFP_FL.Bicc1.RoKAI.2$FC_FR)

### prepare RoKAI input and save data
glimpse(GFP_FL.Bicc1.RoKAI.2) 
GFP_FL.Bicc1.RoKAI.2.23.11.21 <- GFP_FL.Bicc1.RoKAI.2 %>% select(Protein, Position, log2FC)
table(is.na(GFP_FL.Bicc1.RoKAI.2.23.11.21$Protein))
GFP_FL.Bicc1.RoKAI.2.23.11.21 <- GFP_FL.Bicc1.RoKAI.2.23.11.21 %>% filter(!is.na(Protein))
glimpse(GFP_FL.Bicc1.RoKAI.2.23.11.21)
table(is.na(GFP_FL.Bicc1.RoKAI.2.23.11.21$Protein))
table(is.na(GFP_FL.Bicc1.RoKAI.2.23.11.21$Position))
table(is.na(GFP_FL.Bicc1.RoKAI.2.23.11.21$log2FC))
colnames(GFP_FL.Bicc1.RoKAI.2.23.11.21) <- c("Protein", "Position", "Quantification")

#write_csv(GFP_FL.Bicc1.RoKAI.2.23.11.21,"GFP_FL.Bicc1.RoKAI.2_23.11.21.csv", col_names=TRUE)
#save(GFP_FL.Bicc1.RoKAI , file="GFP_FL.Bicc1.RoKAI.Rdata")


# perform analysis via https://rokai.io/
# tool settings was default: Reference Proteome: Uniprot Mouse / Fold changes: normalized / Kinase substrate dataset: PSP + Signor / RoKAI network: KS+PPI+SD+CoEv
# save RoKAI output kinase_table.csv and kinase_targets.csv and load these later




### GFP vs dE18.Bicc1 & RoKAI ##########################################################################################################################################
### GFP vs dE18.Bicc1 & RoKAI ##########################################################################################################################################
glimpse(DZ.IMAC.MBR.ON_expanded_log2_norm.class1.PS)

# select columns of interest
GFP_dE18.Bicc1 <- dplyr::select(DZ.IMAC.MBR.ON_expanded_log2_norm.class1.PS,
                                "ID.FR.all.C1.PS",
                                "id",                          "Proteins",                    "Positions.within.proteins",   "Leading.proteins",            "Protein",                    
                                "Protein.names",               "Gene.names",                  "Fasta.headers",               "Localization.prob",           "Number.of.Phospho.(STY)",    
                                "Amino.acid",                  "Sequence.window",             "Phospho.(STY).Probabilities", "Position.in.peptide",         "Reverse",                    
                                "Potential.contaminant",       "Positions",                   "Position",                    "Peptide.IDs",                 "Mod..peptide.IDs",   
                                "PS_Multiplicity", "Prot.Name_AApos", "Gene.Name_AApos",
                                
                                "Norm.Intensity.A", "Norm.Intensity.F", "Norm.Intensity.K", "Norm.Intensity.P",
                                "Norm.Intensity.E", "Norm.Intensity.J", "Norm.Intensity.O"
                                
)
glimpse(GFP_dE18.Bicc1)
colnames(GFP_dE18.Bicc1)

##check data completeness
GFP_dE18.Bicc1$data.presence.all <- apply(GFP_dE18.Bicc1[, c("Norm.Intensity.A", "Norm.Intensity.F", "Norm.Intensity.K", "Norm.Intensity.P",    
                                                             "Norm.Intensity.E", "Norm.Intensity.J", "Norm.Intensity.O")], 1 , function(x) sum(!is.na(x))) #count data completeness
glimpse(GFP_dE18.Bicc1) #16430


#remove complete zero rows
table(GFP_dE18.Bicc1$data.presence.all == 0) #F15035 T1395

GFP_dE18.Bicc1 <- GFP_dE18.Bicc1 %>% filter(data.presence.all > 0)
glimpse(GFP_dE18.Bicc1) #15035

##X out of Y filter
#count data presence and filter
GFP_dE18.Bicc1$data.presence.left <- apply(GFP_dE18.Bicc1[, c("Norm.Intensity.A", "Norm.Intensity.F", "Norm.Intensity.K", "Norm.Intensity.P")],1, function(x) sum(!is.na(x)))
GFP_dE18.Bicc1$data.presence.right <- apply(GFP_dE18.Bicc1[, c( "Norm.Intensity.E", "Norm.Intensity.J", "Norm.Intensity.O")],1, function(x) sum(!is.na(x)))
glimpse(GFP_dE18.Bicc1)


temp.left <- filter(GFP_dE18.Bicc1, data.presence.left >= 3) ; nrow(temp.left)    #8462
temp.right <- filter(GFP_dE18.Bicc1, data.presence.right >= 3) ; nrow(temp.right) #7116

temp.left.and.right <- bind_rows(temp.left , temp.right) ; nrow(temp.left.and.right)
temp.left.and.right <- distinct(temp.left.and.right) ; nrow(temp.left.and.right) #9656

GFP_dE18.Bicc1.XYfilter <- temp.left.and.right
nrow(GFP_dE18.Bicc1.XYfilter) #9656
length(unique(GFP_dE18.Bicc1.XYfilter$Prot.Name_AApos)) #8551

glimpse(GFP_dE18.Bicc1.XYfilter)

#add FRID
GFP_dE18.Bicc1.XYfilter$FRID <- 1:nrow(GFP_dE18.Bicc1.XYfilter) #here FRID is also row number
glimpse(GFP_dE18.Bicc1.XYfilter)

#find likely Limma NA pval situations 
possible.NApvals.A <- GFP_dE18.Bicc1.XYfilter %>% filter((data.presence.left >= 4 & data.presence.right == 0) | (data.presence.left == 0 & data.presence.right >= 3)  )
glimpse(possible.NApvals.A)
length(unique(possible.NApvals.A$FRID))
FRIDS.possible.NApvals <- possible.NApvals.A$FRID
FRIDS.possible.NApvals 

# for these filtered  black and white situation add zero so that p.val can be calculated
W <- GFP_dE18.Bicc1.XYfilter[FRIDS.possible.NApvals, ]
glimpse(W)
W <- W %>% mutate_at(vars("Norm.Intensity.A", "Norm.Intensity.F", "Norm.Intensity.K", "Norm.Intensity.P",    
                          "Norm.Intensity.E", "Norm.Intensity.J", "Norm.Intensity.O"), ~replace(., is.na(.), 0))#https://stackoverflow.com/questions/8161836/how-do-i-replace-na-values-with-zeros-in-an-r-dataframe
glimpse(W)

#combine again
GFP_dE18.Bicc1.XYfilter.2 <- GFP_dE18.Bicc1.XYfilter %>% filter(! FRID %in% FRIDS.possible.NApvals ) 
GFP_dE18.Bicc1.XYfilter.2 <- bind_rows(GFP_dE18.Bicc1.XYfilter.2, W)
GFP_dE18.Bicc1.XYfilter.2 <- GFP_dE18.Bicc1.XYfilter.2 %>% arrange(FRID)
glimpse(GFP_dE18.Bicc1.XYfilter.2)


### LIMMA
# get data
DataforLimma <- GFP_dE18.Bicc1.XYfilter.2[, c("Norm.Intensity.A", "Norm.Intensity.F", "Norm.Intensity.K", "Norm.Intensity.P",    
                                              "Norm.Intensity.E", "Norm.Intensity.J", "Norm.Intensity.O")]
glimpse(DataforLimma)
colnames(DataforLimma)

design <- cbind(left=c(rep(1, 4), rep(0, 3)), right=c(rep(0, 4), rep(1, 3))); design

fit <- lmFit(DataforLimma, design) # Warning: Partial NA coefficients for 326  probe(s) 
cont.matrix <- makeContrasts( left.vs.right =  right - left, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)


Limma.p.value <- as.numeric(fit2$p.value)
glimpse(Limma.p.value )
Limma.adj.p.valueBH <- p.adjust(Limma.p.value, "BH") 
glimpse(Limma.adj.p.valueBH)
Limma.t <- as.numeric(fit2$t)
glimpse(Limma.t)

GFP_dE18.Bicc1.XYfilter.2$Limma.p.value <- as.numeric(Limma.p.value)
GFP_dE18.Bicc1.XYfilter.2$Limma.adj.p.valueBH <- as.numeric(Limma.adj.p.valueBH)
GFP_dE18.Bicc1.XYfilter.2$Limma.t <- as.numeric(Limma.t)

GFP_dE18.Bicc1.XYfilter.B <- GFP_dE18.Bicc1.XYfilter.2
glimpse(GFP_dE18.Bicc1.XYfilter.B)

### calculate FC ###
for(i in 1:nrow(GFP_dE18.Bicc1.XYfilter.B)){
  print(i)
  left.vals <- as.numeric(GFP_dE18.Bicc1.XYfilter.B[i,c("Norm.Intensity.A", "Norm.Intensity.F", "Norm.Intensity.K", "Norm.Intensity.P")])
  right.vals <- as.numeric(GFP_dE18.Bicc1.XYfilter.B[i,c("Norm.Intensity.E", "Norm.Intensity.J", "Norm.Intensity.O")])
  FC <- FR_FC_LOGData_2.3.17(control = left.vals, treatment = right.vals, LOGtype = 2); FC
  GFP_dE18.Bicc1.XYfilter.B$FC_FR[i] <- FC
}
glimpse(GFP_dE18.Bicc1.XYfilter.B)


#check numbers unique Sequence.window
GFP_dE18.Bicc1.XYfilter.B %>% filter(Limma.p.value < 0.05) %>% select(Sequence.window) %>% distinct() %>% nrow() #3935 total sign
GFP_dE18.Bicc1.XYfilter.B %>% filter(Limma.p.value < 0.05) %>% filter(FC_FR > 1) %>% select(Sequence.window) %>% distinct() %>% nrow() #1987 up
GFP_dE18.Bicc1.XYfilter.B %>% filter(Limma.p.value < 0.05) %>% filter(FC_FR < 1) %>% select(Sequence.window) %>% distinct() %>% nrow() #1977 down

# filter na pvalues
table(is.na(GFP_dE18.Bicc1.XYfilter.B$Limma.p.value))
GFP_dE18.Bicc1.XYfilter.B <- GFP_dE18.Bicc1.XYfilter.B %>% filter(!is.na(Limma.p.value))
glimpse(GFP_dE18.Bicc1.XYfilter.B) #9330
table(is.na(GFP_dE18.Bicc1.XYfilter.B$Limma.p.value))


# check +/- 10000 FCs
table(GFP_dE18.Bicc1.XYfilter.B$FC_FR >= 10000) #420
table(GFP_dE18.Bicc1.XYfilter.B$FC_FR <= -10000) #251
summary(GFP_dE18.Bicc1.XYfilter.B$FC_FR)

#find second min or max and add +/-5 for +/- 10000 FCs
temp <- GFP_dE18.Bicc1.XYfilter.B %>% filter(!FC_FR %in% c(10000, -10000)) %>% arrange(desc(FC_FR))
head(temp$FC_FR)
tail(temp$FC_FR)
max(temp$FC_FR)
min(temp$FC_FR)

GFP_dE18.Bicc1.XYfilter.B$FC_FR <- replace(GFP_dE18.Bicc1.XYfilter.B$FC_FR, GFP_dE18.Bicc1.XYfilter.B$FC_FR==10000, max(temp$FC_FR)+5)
GFP_dE18.Bicc1.XYfilter.B$FC_FR <- replace(GFP_dE18.Bicc1.XYfilter.B$FC_FR, GFP_dE18.Bicc1.XYfilter.B$FC_FR==-10000, min(temp$FC_FR)-5)
summary(GFP_dE18.Bicc1.XYfilter.B$FC_FR)
glimpse(GFP_dE18.Bicc1.XYfilter.B) #9330

## GFP_dE18.Bicc1 reshape for RoKAI
GFP_dE18.Bicc1.RoKAI <- GFP_dE18.Bicc1.XYfilter.B 

# make aa15.window
# take first entry sequence window if there are two entries
GFP_dE18.Bicc1.RoKAI$aa15.window <- pbsapply(GFP_dE18.Bicc1.RoKAI$Sequence.window, function(x) unlist(str_split(string =x, pattern = ";"))[1]  )
GFP_dE18.Bicc1.RoKAI$aa15.window <- pbsapply(GFP_dE18.Bicc1.RoKAI$aa15.window, function(x) make.aa.15.window.from.31.sequence.window(input = x)  )

glimpse(GFP_dE18.Bicc1.RoKAI)

table(duplicated(GFP_dE18.Bicc1.RoKAI$aa15.window))
which(duplicated(GFP_dE18.Bicc1.RoKAI$aa15.window))


## if there are duplicate entries (e.g. because of _1_2_3) take the entry with the best p value
unique.15aa.windows <- unique(GFP_dE18.Bicc1.RoKAI$aa15.window)
length(unique.15aa.windows) #8247
glimpse(GFP_dE18.Bicc1.RoKAI) #9330


GFP_dE18.Bicc1.RoKAI.2 <- c()
for (i in unique.15aa.windows){
  print(i)
  temp <- GFP_dE18.Bicc1.RoKAI %>% filter(aa15.window %in% c(i))
  temp <- temp %>% arrange(Limma.p.value)
  temp <- temp[1,]
  GFP_dE18.Bicc1.RoKAI.2 <- bind_rows(GFP_dE18.Bicc1.RoKAI.2, temp)
}

glimpse(GFP_dE18.Bicc1.RoKAI.2)
table(duplicated(GFP_dE18.Bicc1.RoKAI.2$aa15.window))


#calculate log2(FC) and order
GFP_dE18.Bicc1.RoKAI.2$log2FC <- log2(abs(GFP_dE18.Bicc1.RoKAI.2$FC_FR ))
GFP_dE18.Bicc1.RoKAI.2$sign.FC_FR <- sign(GFP_dE18.Bicc1.RoKAI.2$FC_FR)
GFP_dE18.Bicc1.RoKAI.2$log2FC <- sign(GFP_dE18.Bicc1.RoKAI.2$FC_FR) * log2(abs(GFP_dE18.Bicc1.RoKAI.2$FC_FR ))
glimpse(GFP_dE18.Bicc1.RoKAI.2) 

GFP_dE18.Bicc1.RoKAI.2 <- GFP_dE18.Bicc1.RoKAI.2 %>% arrange(desc(log2FC))
head(GFP_dE18.Bicc1.RoKAI.2$FC_FR)
tail(GFP_dE18.Bicc1.RoKAI.2$FC_FR)

### prepare RoKAI input and save data
glimpse(GFP_dE18.Bicc1.RoKAI.2) 
GFP_dE18.Bicc1.RoKAI.2.23.11.21 <- GFP_dE18.Bicc1.RoKAI.2 %>% select(Protein, Position, log2FC)
table(is.na(GFP_dE18.Bicc1.RoKAI.2.23.11.21$Protein))
GFP_dE18.Bicc1.RoKAI.2.23.11.21 <- GFP_dE18.Bicc1.RoKAI.2.23.11.21 %>% filter(!is.na(Protein))
glimpse(GFP_dE18.Bicc1.RoKAI.2.23.11.21)
table(is.na(GFP_dE18.Bicc1.RoKAI.2.23.11.21$Protein))
table(is.na(GFP_dE18.Bicc1.RoKAI.2.23.11.21$Position))
table(is.na(GFP_dE18.Bicc1.RoKAI.2.23.11.21$log2FC))
colnames(GFP_dE18.Bicc1.RoKAI.2.23.11.21) <- c("Protein", "Position", "Quantification")

#write_csv(GFP_dE18.Bicc1.RoKAI.2.23.11.21,"GFP_dE18.Bicc1.RoKAI.2_23.11.21.csv", col_names=TRUE)
#save(GFP_dE18.Bicc1.RoKAI, file="GFP_dE18.Bicc1.RoKAI.Rdata")


# perform analysis via https://rokai.io/
# tool settings was default: Reference Proteome: Uniprot Mouse / Fold changes: normalized / Kinase substrate dataset: PSP + Signor / RoKAI network: KS+PPI+SD+CoEv
# save RoKAI output kinase_table.csv and kinase_targets.csv and load these later


### FL vs dE18 & RoKAI ##########################################################################################################################################
### FL vs dE18 & RoKAI ##########################################################################################################################################
glimpse(DZ.IMAC.MBR.ON_expanded_log2_norm.class1.PS)

# select columns of interest
FL_dE18 <- dplyr::select(DZ.IMAC.MBR.ON_expanded_log2_norm.class1.PS,
                         "ID.FR.all.C1.PS",
                         "id",                          "Proteins",                    "Positions.within.proteins",   "Leading.proteins",            "Protein",                    
                         "Protein.names",               "Gene.names",                  "Fasta.headers",               "Localization.prob",           "Number.of.Phospho.(STY)",    
                         "Amino.acid",                  "Sequence.window",             "Phospho.(STY).Probabilities", "Position.in.peptide",         "Reverse",                    
                         "Potential.contaminant",       "Positions",                   "Position",                    "Peptide.IDs",                 "Mod..peptide.IDs",   
                         "PS_Multiplicity", "Prot.Name_AApos", "Gene.Name_AApos",
                         
                         "Norm.Intensity.B", "Norm.Intensity.G", "Norm.Intensity.L", "Norm.Intensity.Q",
                         "Norm.Intensity.C", "Norm.Intensity.H", "Norm.Intensity.M", "Norm.Intensity.R"
                         
)
glimpse(FL_dE18)
colnames(FL_dE18)

##check data completeness
FL_dE18$data.presence.all <- apply(FL_dE18[, c("Norm.Intensity.B", "Norm.Intensity.G", "Norm.Intensity.L", "Norm.Intensity.Q",    
                                               "Norm.Intensity.C", "Norm.Intensity.H", "Norm.Intensity.M", "Norm.Intensity.R")], 1 , function(x) sum(!is.na(x))) #count data completeness
glimpse(FL_dE18) #16430


#remove complete zero rows
table(FL_dE18$data.presence.all == 0) #F14416 T2014

FL_dE18 <- FL_dE18 %>% filter(data.presence.all > 0)
glimpse(FL_dE18) #14416

##X out of Y filter 
#count data presence and filter
FL_dE18$data.presence.left <- apply(FL_dE18[, c("Norm.Intensity.B", "Norm.Intensity.G", "Norm.Intensity.L", "Norm.Intensity.Q")],1, function(x) sum(!is.na(x)))
FL_dE18$data.presence.right <- apply(FL_dE18[, c( "Norm.Intensity.C", "Norm.Intensity.H", "Norm.Intensity.M", "Norm.Intensity.R")],1, function(x) sum(!is.na(x)))
glimpse(FL_dE18)


temp.left <- filter(FL_dE18, data.presence.left >= 3) ; nrow(temp.left)    #8544
temp.right <- filter(FL_dE18, data.presence.right >= 3) ; nrow(temp.right) #8438

temp.left.and.right <- bind_rows(temp.left , temp.right) ; nrow(temp.left.and.right)
temp.left.and.right <- distinct(temp.left.and.right) ; nrow(temp.left.and.right) #9783

FL_dE18.XYfilter <- temp.left.and.right
nrow(FL_dE18.XYfilter) #9783
length(unique(FL_dE18.XYfilter$Prot.Name_AApos)) #8630


#add FRID
FL_dE18.XYfilter$FRID <- 1:nrow(FL_dE18.XYfilter) #here FRID is also row number
glimpse(FL_dE18.XYfilter)

#find likely Limma NA pval situations 
possible.NApvals.A <- FL_dE18.XYfilter %>% filter((data.presence.left >= 4 & data.presence.right == 0) | (data.presence.left == 0 & data.presence.right >= 4)  )
glimpse(possible.NApvals.A)
length(unique(possible.NApvals.A$FRID))
FRIDS.possible.NApvals <- possible.NApvals.A$FRID
FRIDS.possible.NApvals 

# for these filtered  black and white situation add zero so that p.val can be calculated
W <- FL_dE18.XYfilter[FRIDS.possible.NApvals, ]
glimpse(W)
W <- W %>% mutate_at(vars("Norm.Intensity.B", "Norm.Intensity.G", "Norm.Intensity.L", "Norm.Intensity.Q",    
                          "Norm.Intensity.C", "Norm.Intensity.H", "Norm.Intensity.M", "Norm.Intensity.R"), ~replace(., is.na(.), 0))
glimpse(W)


#combine again
FL_dE18.XYfilter.2 <- FL_dE18.XYfilter %>% filter(! FRID %in% FRIDS.possible.NApvals ) 
FL_dE18.XYfilter.2 <- bind_rows(FL_dE18.XYfilter.2, W)
FL_dE18.XYfilter.2 <- FL_dE18.XYfilter.2 %>% arrange(FRID)
glimpse(FL_dE18.XYfilter.2)


### LIMMA
# get data
DataforLimma <- FL_dE18.XYfilter.2[, c("Norm.Intensity.B", "Norm.Intensity.G", "Norm.Intensity.L", "Norm.Intensity.Q",    
                                       "Norm.Intensity.C", "Norm.Intensity.H", "Norm.Intensity.M", "Norm.Intensity.R")]
glimpse(DataforLimma)
colnames(DataforLimma)

design <- cbind(left=c(rep(1, 4), rep(0, 4)), right=c(rep(0, 4), rep(1, 4))); design

fit <- lmFit(DataforLimma, design) # Warning: Partial NA coefficients for 407 probe(s) 
cont.matrix <- makeContrasts( left.vs.right =  right - left, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

Limma.p.value <- as.numeric(fit2$p.value)
glimpse(Limma.p.value )
Limma.adj.p.valueBH <- p.adjust(Limma.p.value, "BH") 
glimpse(Limma.adj.p.valueBH)
Limma.t <- as.numeric(fit2$t)
glimpse(Limma.t)

FL_dE18.XYfilter.2$Limma.p.value <- as.numeric(Limma.p.value)
FL_dE18.XYfilter.2$Limma.adj.p.valueBH <- as.numeric(Limma.adj.p.valueBH)
FL_dE18.XYfilter.2$Limma.t <- as.numeric(Limma.t)

FL_dE18.XYfilter.B <- FL_dE18.XYfilter.2
glimpse(FL_dE18.XYfilter.B)

### calculate FC
for(i in 1:nrow(FL_dE18.XYfilter.B)){
  print(i)
  left.vals <- as.numeric(FL_dE18.XYfilter.B[i,c("Norm.Intensity.B", "Norm.Intensity.G", "Norm.Intensity.L", "Norm.Intensity.Q")])
  right.vals <- as.numeric(FL_dE18.XYfilter.B[i,c("Norm.Intensity.C", "Norm.Intensity.H", "Norm.Intensity.M", "Norm.Intensity.R")])
  FC <- FR_FC_LOGData_2.3.17(control = left.vals, treatment = right.vals, LOGtype = 2); FC
  FL_dE18.XYfilter.B$FC_FR[i] <- FC
}
glimpse(FL_dE18.XYfilter.B)

#check numbers unique Sequence.window
FL_dE18.XYfilter.B %>% filter(Limma.p.value < 0.05) %>% select(Sequence.window) %>% distinct() %>% nrow() #2852 total sign
FL_dE18.XYfilter.B %>% filter(Limma.p.value < 0.05) %>% filter(FC_FR > 1) %>% select(Sequence.window) %>% distinct() %>% nrow() #1373 up
FL_dE18.XYfilter.B %>% filter(Limma.p.value < 0.05) %>% filter(FC_FR < 1) %>% select(Sequence.window) %>% distinct() %>% nrow() #1506 down

# filter na pvalues
table(is.na(FL_dE18.XYfilter.B$Limma.p.value))
FL_dE18.XYfilter.B <- FL_dE18.XYfilter.B %>% filter(!is.na(Limma.p.value))
glimpse(FL_dE18.XYfilter.B) #9376
table(is.na(FL_dE18.XYfilter.B$Limma.p.value))


# check +/- 10000 FCs
table(FL_dE18.XYfilter.B$FC_FR >= 10000) #140
table(FL_dE18.XYfilter.B$FC_FR <= -10000) #97
summary(FL_dE18.XYfilter.B$FC_FR)

#find second min or max and add +/-5 for +/- 10000 FCs
temp <- FL_dE18.XYfilter.B %>% filter(!FC_FR %in% c(10000, -10000)) %>% arrange(desc(FC_FR))
head(temp$FC_FR)
tail(temp$FC_FR)
max(temp$FC_FR)
min(temp$FC_FR)

FL_dE18.XYfilter.B$FC_FR <- replace(FL_dE18.XYfilter.B$FC_FR, FL_dE18.XYfilter.B$FC_FR==10000, max(temp$FC_FR)+5)
FL_dE18.XYfilter.B$FC_FR <- replace(FL_dE18.XYfilter.B$FC_FR, FL_dE18.XYfilter.B$FC_FR==-10000, min(temp$FC_FR)-5)
summary(FL_dE18.XYfilter.B$FC_FR)
glimpse(FL_dE18.XYfilter.B) #9376

## FL_dE18 reshape for RoKAI 
FL_dE18.RoKAI <- FL_dE18.XYfilter.B 

# make aa15.window
# take first entry sequence window if there are two entries
FL_dE18.RoKAI$aa15.window <- pbsapply(FL_dE18.RoKAI$Sequence.window, function(x) unlist(str_split(string =x, pattern = ";"))[1]  )
FL_dE18.RoKAI$aa15.window <- pbsapply(FL_dE18.RoKAI$aa15.window, function(x) make.aa.15.window.from.31.sequence.window(input = x)  )


glimpse(FL_dE18.RoKAI)

table(duplicated(FL_dE18.RoKAI$aa15.window))
which(duplicated(FL_dE18.RoKAI$aa15.window))


## if there are duplicate entries (e.g. because of _1_2_3) take the entry with the best p value
unique.15aa.windows <- unique(FL_dE18.RoKAI$aa15.window)
length(unique.15aa.windows) #8264
glimpse(FL_dE18.RoKAI) #9376


FL_dE18.RoKAI.2 <- c()
for (i in unique.15aa.windows){
  print(i)
  temp <- FL_dE18.RoKAI %>% filter(aa15.window %in% c(i))
  temp <- temp %>% arrange(Limma.p.value)
  temp <- temp[1,]
  FL_dE18.RoKAI.2 <- bind_rows(FL_dE18.RoKAI.2, temp)
}

glimpse(FL_dE18.RoKAI.2)
table(duplicated(FL_dE18.RoKAI.2$aa15.window))

#calculate log2(FC) and order
FL_dE18.RoKAI.2$log2FC <- log2(abs(FL_dE18.RoKAI.2$FC_FR ))
FL_dE18.RoKAI.2$sign.FC_FR <- sign(FL_dE18.RoKAI.2$FC_FR)
FL_dE18.RoKAI.2$log2FC <- sign(FL_dE18.RoKAI.2$FC_FR) * log2(abs(FL_dE18.RoKAI.2$FC_FR ))
glimpse(FL_dE18.RoKAI.2) 

FL_dE18.RoKAI.2 <- FL_dE18.RoKAI.2 %>% arrange(desc(log2FC))
head(FL_dE18.RoKAI.2$FC_FR)
tail(FL_dE18.RoKAI.2$FC_FR)

###prepare RoKAI input and save data
glimpse(FL_dE18.RoKAI.2) 
FL_dE18.RoKAI.2.23.11.21 <- FL_dE18.RoKAI.2 %>% select(Protein, Position, log2FC) 
table(is.na(FL_dE18.RoKAI.2.23.11.21$Protein))
FL_dE18.RoKAI.2.23.11.21 <- FL_dE18.RoKAI.2.23.11.21 %>% filter(!is.na(Protein))
glimpse(FL_dE18.RoKAI.2.23.11.21)
table(is.na(FL_dE18.RoKAI.2.23.11.21$Protein))
table(is.na(FL_dE18.RoKAI.2.23.11.21$Position))
table(is.na(FL_dE18.RoKAI.2.23.11.21$log2FC))
colnames(FL_dE18.RoKAI.2.23.11.21) <- c("Protein", "Position", "Quantification")

#write_csv(FL_dE18.RoKAI.2.23.11.21,"FL_dE18.RoKAI.2_23.11.21.csv", col_names=TRUE)
#save(FL_dE18.RoKAI, file="FL_dE18.RoKAI.Rdata")

# perform analysis via https://rokai.io/
# tool settings was default: Reference Proteome: Uniprot Mouse / Fold changes: normalized / Kinase substrate dataset: PSP + Signor / RoKAI network: KS+PPI+SD+CoEv
# save RoKAI output kinase_table.csv and kinase_targets.csv and load these later


### FL vs FL.Bicc1 & RoKAI ##########################################################################################################################################
### FL vs FL.Bicc1 &  RoKAI ##########################################################################################################################################
glimpse(DZ.IMAC.MBR.ON_expanded_log2_norm.class1.PS)

# select columns of interest
FL_FL.Bicc1 <- dplyr::select(DZ.IMAC.MBR.ON_expanded_log2_norm.class1.PS,
                             "ID.FR.all.C1.PS",
                             "id",                          "Proteins",                    "Positions.within.proteins",   "Leading.proteins",            "Protein",                    
                             "Protein.names",               "Gene.names",                  "Fasta.headers",               "Localization.prob",           "Number.of.Phospho.(STY)",    
                             "Amino.acid",                  "Sequence.window",             "Phospho.(STY).Probabilities", "Position.in.peptide",         "Reverse",                    
                             "Potential.contaminant",       "Positions",                   "Position",                    "Peptide.IDs",                 "Mod..peptide.IDs",   
                             "PS_Multiplicity", "Prot.Name_AApos", "Gene.Name_AApos",
                             
                             "Norm.Intensity.B", "Norm.Intensity.G", "Norm.Intensity.L", "Norm.Intensity.Q",
                             "Norm.Intensity.D", "Norm.Intensity.I", "Norm.Intensity.N", "Norm.Intensity.S"
                             
)
glimpse(FL_FL.Bicc1)
colnames(FL_FL.Bicc1)

##check data completeness
FL_FL.Bicc1$data.presence.all <- apply(FL_FL.Bicc1[, c("Norm.Intensity.B", "Norm.Intensity.G", "Norm.Intensity.L", "Norm.Intensity.Q",    
                                                       "Norm.Intensity.D", "Norm.Intensity.I", "Norm.Intensity.N", "Norm.Intensity.S")], 1 , function(x) sum(!is.na(x))) #count data completeness
glimpse(FL_FL.Bicc1) #16430


#remove complete zero rows
table(FL_FL.Bicc1$data.presence.all == 0) #F14394 T2036

FL_FL.Bicc1 <- FL_FL.Bicc1 %>% filter(data.presence.all > 0)
glimpse(FL_FL.Bicc1) #14394

##X out of Y filter
#count data presence and filter
FL_FL.Bicc1$data.presence.left <- apply(FL_FL.Bicc1[, c("Norm.Intensity.B", "Norm.Intensity.G", "Norm.Intensity.L", "Norm.Intensity.Q")],1, function(x) sum(!is.na(x)))
FL_FL.Bicc1$data.presence.right <- apply(FL_FL.Bicc1[, c( "Norm.Intensity.D", "Norm.Intensity.I", "Norm.Intensity.N", "Norm.Intensity.S")],1, function(x) sum(!is.na(x)))
glimpse(FL_FL.Bicc1)


temp.left <- filter(FL_FL.Bicc1, data.presence.left >= 3) ; nrow(temp.left)    #8544
temp.right <- filter(FL_FL.Bicc1, data.presence.right >= 3) ; nrow(temp.right) #8866

temp.left.and.right <- bind_rows(temp.left , temp.right) ; nrow(temp.left.and.right)
temp.left.and.right <- distinct(temp.left.and.right) ; nrow(temp.left.and.right) #9897

FL_FL.Bicc1.XYfilter <- temp.left.and.right
nrow(FL_FL.Bicc1.XYfilter) #9897
length(unique(FL_FL.Bicc1.XYfilter$Prot.Name_AApos)) #8714

#add FRID
FL_FL.Bicc1.XYfilter$FRID <- 1:nrow(FL_FL.Bicc1.XYfilter) #here FRID is also row number
glimpse(FL_FL.Bicc1.XYfilter)

#find likely Limma NA pval situations 
possible.NApvals.A <- FL_FL.Bicc1.XYfilter %>% filter((data.presence.left >= 4 & data.presence.right == 0) | (data.presence.left == 0 & data.presence.right >= 4)  )
glimpse(possible.NApvals.A)
length(unique(possible.NApvals.A$FRID))
FRIDS.possible.NApvals <- possible.NApvals.A$FRID
FRIDS.possible.NApvals 

# for these filtered  black and white situation add zero so that p.val can be calculated
W <- FL_FL.Bicc1.XYfilter[FRIDS.possible.NApvals, ]
glimpse(W)
W <- W %>% mutate_at(vars("Norm.Intensity.B", "Norm.Intensity.G", "Norm.Intensity.L", "Norm.Intensity.Q",    
                          "Norm.Intensity.D", "Norm.Intensity.I", "Norm.Intensity.N", "Norm.Intensity.S"), ~replace(., is.na(.), 0))#https://stackoverflow.com/questions/8161836/how-do-i-replace-na-values-with-zeros-in-an-r-dataframe
glimpse(W)

#combine again
FL_FL.Bicc1.XYfilter.2 <- FL_FL.Bicc1.XYfilter %>% filter(! FRID %in% FRIDS.possible.NApvals ) 
FL_FL.Bicc1.XYfilter.2 <- bind_rows(FL_FL.Bicc1.XYfilter.2, W)
FL_FL.Bicc1.XYfilter.2 <- FL_FL.Bicc1.XYfilter.2 %>% arrange(FRID)
glimpse(FL_FL.Bicc1.XYfilter.2)

### LIMMA
# get data
DataforLimma <- FL_FL.Bicc1.XYfilter.2[, c("Norm.Intensity.B", "Norm.Intensity.G", "Norm.Intensity.L", "Norm.Intensity.Q",    
                                           "Norm.Intensity.D", "Norm.Intensity.I", "Norm.Intensity.N", "Norm.Intensity.S")]
glimpse(DataforLimma)
colnames(DataforLimma)

design <- cbind(left=c(rep(1, 4), rep(0, 4)), right=c(rep(0, 4), rep(1, 4))); design

fit <- lmFit(DataforLimma, design) # Warning: Partial NA coefficients for 315 probe(s) 
cont.matrix <- makeContrasts( left.vs.right =  right - left, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

Limma.p.value <- as.numeric(fit2$p.value)
glimpse(Limma.p.value )
Limma.adj.p.valueBH <- p.adjust(Limma.p.value, "BH") 
glimpse(Limma.adj.p.valueBH)
Limma.t <- as.numeric(fit2$t)
glimpse(Limma.t)

FL_FL.Bicc1.XYfilter.2$Limma.p.value <- as.numeric(Limma.p.value)
FL_FL.Bicc1.XYfilter.2$Limma.adj.p.valueBH <- as.numeric(Limma.adj.p.valueBH)
FL_FL.Bicc1.XYfilter.2$Limma.t <- as.numeric(Limma.t)

FL_FL.Bicc1.XYfilter.B <- FL_FL.Bicc1.XYfilter.2
glimpse(FL_FL.Bicc1.XYfilter.B)



### calculate FC
for(i in 1:nrow(FL_FL.Bicc1.XYfilter.B)){
  print(i)
  left.vals <- as.numeric(FL_FL.Bicc1.XYfilter.B[i,c("Norm.Intensity.B", "Norm.Intensity.G", "Norm.Intensity.L", "Norm.Intensity.Q")])
  right.vals <- as.numeric(FL_FL.Bicc1.XYfilter.B[i,c("Norm.Intensity.D", "Norm.Intensity.I", "Norm.Intensity.N", "Norm.Intensity.S")])
  FC <- FR_FC_LOGData_2.3.17(control = left.vals, treatment = right.vals, LOGtype = 2); FC
  FL_FL.Bicc1.XYfilter.B$FC_FR[i] <- FC
}
glimpse(FL_FL.Bicc1.XYfilter.B)


#check numbers unique Sequence.window
FL_FL.Bicc1.XYfilter.B %>% filter(Limma.p.value < 0.05) %>% select(Sequence.window) %>% distinct() %>% nrow() #1744 total sign
FL_FL.Bicc1.XYfilter.B %>% filter(Limma.p.value < 0.05) %>% filter(FC_FR > 1) %>% select(Sequence.window) %>% distinct() %>% nrow() #1119 up
FL_FL.Bicc1.XYfilter.B %>% filter(Limma.p.value < 0.05) %>% filter(FC_FR < 1) %>% select(Sequence.window) %>% distinct() %>% nrow() #631 down


table(is.na(FL_FL.Bicc1.XYfilter.B$Limma.p.value))
FL_FL.Bicc1.XYfilter.B <- FL_FL.Bicc1.XYfilter.B %>% filter(!is.na(Limma.p.value))
glimpse(FL_FL.Bicc1.XYfilter.B) #9582
table(is.na(FL_FL.Bicc1.XYfilter.B$Limma.p.value))


# check +/- 10000 FCs
table(FL_FL.Bicc1.XYfilter.B$FC_FR >= 10000) #122
table(FL_FL.Bicc1.XYfilter.B$FC_FR <= -10000) #31
summary(FL_FL.Bicc1.XYfilter.B$FC_FR)

#find second min or max and add +/-5 for +/- 10000 FCs
temp <- FL_FL.Bicc1.XYfilter.B %>% filter(!FC_FR %in% c(10000, -10000)) %>% arrange(desc(FC_FR))
head(temp$FC_FR)
tail(temp$FC_FR)
max(temp$FC_FR)
min(temp$FC_FR)

FL_FL.Bicc1.XYfilter.B$FC_FR <- replace(FL_FL.Bicc1.XYfilter.B$FC_FR, FL_FL.Bicc1.XYfilter.B$FC_FR==10000, max(temp$FC_FR)+5)
FL_FL.Bicc1.XYfilter.B$FC_FR <- replace(FL_FL.Bicc1.XYfilter.B$FC_FR, FL_FL.Bicc1.XYfilter.B$FC_FR==-10000, min(temp$FC_FR)-5)
summary(FL_FL.Bicc1.XYfilter.B$FC_FR)
glimpse(FL_FL.Bicc1.XYfilter.B) #9582

## FL_FL.Bicc1 reshape for RoKAI 

FL_FL.Bicc1.RoKAI <- FL_FL.Bicc1.XYfilter.B 

# make aa15.window
# take first entry sequence window if there are two entries
FL_FL.Bicc1.RoKAI$aa15.window <- pbsapply(FL_FL.Bicc1.RoKAI$Sequence.window, function(x) unlist(str_split(string =x, pattern = ";"))[1]  )
FL_FL.Bicc1.RoKAI$aa15.window <- pbsapply(FL_FL.Bicc1.RoKAI$aa15.window, function(x) make.aa.15.window.from.31.sequence.window(input = x)  )


glimpse(FL_FL.Bicc1.RoKAI)

table(duplicated(FL_FL.Bicc1.RoKAI$aa15.window))
which(duplicated(FL_FL.Bicc1.RoKAI$aa15.window))


## if there are duplicate entries (e.g. because of _1_2_3) take the entry with the best p value
unique.15aa.windows <- unique(FL_FL.Bicc1.RoKAI$aa15.window)
length(unique.15aa.windows) #8434
glimpse(FL_FL.Bicc1.RoKAI) #9582


FL_FL.Bicc1.RoKAI.2 <- c()
for (i in unique.15aa.windows){
  print(i)
  temp <- FL_FL.Bicc1.RoKAI %>% filter(aa15.window %in% c(i))
  temp <- temp %>% arrange(Limma.p.value)
  temp <- temp[1,]
  FL_FL.Bicc1.RoKAI.2 <- bind_rows(FL_FL.Bicc1.RoKAI.2, temp)
}

glimpse(FL_FL.Bicc1.RoKAI.2)
table(duplicated(FL_FL.Bicc1.RoKAI.2$aa15.window))



#calculate log2(FC) and order
FL_FL.Bicc1.RoKAI.2$log2FC <- log2(abs(FL_FL.Bicc1.RoKAI.2$FC_FR ))
FL_FL.Bicc1.RoKAI.2$sign.FC_FR <- sign(FL_FL.Bicc1.RoKAI.2$FC_FR)
FL_FL.Bicc1.RoKAI.2$log2FC <- sign(FL_FL.Bicc1.RoKAI.2$FC_FR) * log2(abs(FL_FL.Bicc1.RoKAI.2$FC_FR ))
glimpse(FL_FL.Bicc1.RoKAI.2) 

FL_FL.Bicc1.RoKAI.2 <- FL_FL.Bicc1.RoKAI.2 %>% arrange(desc(log2FC))
head(FL_FL.Bicc1.RoKAI.2$FC_FR)
tail(FL_FL.Bicc1.RoKAI.2$FC_FR)

### prepare RoKAI input and save data
glimpse(FL_FL.Bicc1.RoKAI.2) 
FL_FL.Bicc1.RoKAI.2.23.11.21 <- FL_FL.Bicc1.RoKAI.2 %>% select(Protein, Position, log2FC) 
table(is.na(FL_FL.Bicc1.RoKAI.2.23.11.21$Protein))
FL_FL.Bicc1.RoKAI.2.23.11.21 <- FL_FL.Bicc1.RoKAI.2.23.11.21 %>% filter(!is.na(Protein))
glimpse(FL_FL.Bicc1.RoKAI.2.23.11.21)
table(is.na(FL_FL.Bicc1.RoKAI.2.23.11.21$Protein))
table(is.na(FL_FL.Bicc1.RoKAI.2.23.11.21$Position))
table(is.na(FL_FL.Bicc1.RoKAI.2.23.11.21$log2FC))
colnames(FL_FL.Bicc1.RoKAI.2.23.11.21) <- c("Protein", "Position", "Quantification")

#write_csv(FL_FL.Bicc1.RoKAI.2.23.11.21,"FL_FL.Bicc1.RoKAI.2_23.11.21.csv", col_names=TRUE)
#save(FL_FL.Bicc1.RoKAI, file="FL_FL.Bicc1.RoKAI.Rdata")

# perform analysis via https://rokai.io/
# tool settings was default: Reference Proteome: Uniprot Mouse / Fold changes: normalized / Kinase substrate dataset: PSP + Signor / RoKAI network: KS+PPI+SD+CoEv
# save RoKAI output kinase_table.csv and kinase_targets.csv and load these later

### FL vs dE18.Bicc1  RoKAI ##########################################################################################################################################
### FL vs dE18.Bicc1  RoKAI ##########################################################################################################################################
glimpse(DZ.IMAC.MBR.ON_expanded_log2_norm.class1.PS)

# select columns of interest
FL_dE18.Bicc1 <- dplyr::select(DZ.IMAC.MBR.ON_expanded_log2_norm.class1.PS,
                               "ID.FR.all.C1.PS",
                               "id",                          "Proteins",                    "Positions.within.proteins",   "Leading.proteins",            "Protein",                    
                               "Protein.names",               "Gene.names",                  "Fasta.headers",               "Localization.prob",           "Number.of.Phospho.(STY)",    
                               "Amino.acid",                  "Sequence.window",             "Phospho.(STY).Probabilities", "Position.in.peptide",         "Reverse",                    
                               "Potential.contaminant",       "Positions",                   "Position",                    "Peptide.IDs",                 "Mod..peptide.IDs",   
                               "PS_Multiplicity", "Prot.Name_AApos", "Gene.Name_AApos",
                               
                               "Norm.Intensity.B", "Norm.Intensity.G", "Norm.Intensity.L", "Norm.Intensity.Q",
                               "Norm.Intensity.E", "Norm.Intensity.J", "Norm.Intensity.O"
                               
)
glimpse(FL_dE18.Bicc1)
colnames(FL_dE18.Bicc1)

##check data completeness
FL_dE18.Bicc1$data.presence.all <- apply(FL_dE18.Bicc1[, c("Norm.Intensity.B", "Norm.Intensity.G", "Norm.Intensity.L", "Norm.Intensity.Q",    
                                                           "Norm.Intensity.E", "Norm.Intensity.J", "Norm.Intensity.O")], 1 , function(x) sum(!is.na(x))) #count data completeness
glimpse(FL_dE18.Bicc1) #16430


#remove complete zero rows
table(FL_dE18.Bicc1$data.presence.all == 0) #F14701 T1729

FL_dE18.Bicc1 <- FL_dE18.Bicc1 %>% filter(data.presence.all > 0)
glimpse(FL_dE18.Bicc1) #14701

#X out of Y filter 
#count data presence and filter
FL_dE18.Bicc1$data.presence.left <- apply(FL_dE18.Bicc1[, c("Norm.Intensity.B", "Norm.Intensity.G", "Norm.Intensity.L", "Norm.Intensity.Q")],1, function(x) sum(!is.na(x)))
FL_dE18.Bicc1$data.presence.right <- apply(FL_dE18.Bicc1[, c( "Norm.Intensity.E", "Norm.Intensity.J", "Norm.Intensity.O")],1, function(x) sum(!is.na(x)))
glimpse(FL_dE18.Bicc1)


temp.left <- filter(FL_dE18.Bicc1, data.presence.left >= 3) ; nrow(temp.left)    #8544
temp.right <- filter(FL_dE18.Bicc1, data.presence.right >= 3) ; nrow(temp.right) #7116

temp.left.and.right <- bind_rows(temp.left , temp.right) ; nrow(temp.left.and.right)
temp.left.and.right <- distinct(temp.left.and.right) ; nrow(temp.left.and.right) #9622

FL_dE18.Bicc1.XYfilter <- temp.left.and.right
nrow(FL_dE18.Bicc1.XYfilter) #9622
length(unique(FL_dE18.Bicc1.XYfilter$Prot.Name_AApos)) #8519

glimpse(FL_dE18.Bicc1.XYfilter)

#add FRID
FL_dE18.Bicc1.XYfilter$FRID <- 1:nrow(FL_dE18.Bicc1.XYfilter) #here FRID is also row number
glimpse(FL_dE18.Bicc1.XYfilter)

#find likely Limma NA pval situations 
possible.NApvals.A <- FL_dE18.Bicc1.XYfilter %>% filter((data.presence.left >= 4 & data.presence.right == 0) | (data.presence.left == 0 & data.presence.right >= 3)  )
glimpse(possible.NApvals.A)
length(unique(possible.NApvals.A$FRID))
FRIDS.possible.NApvals <- possible.NApvals.A$FRID
FRIDS.possible.NApvals 

# for these filtered  black and white situation add zero so that p.val can be calculated
W <- FL_dE18.Bicc1.XYfilter[FRIDS.possible.NApvals, ]
glimpse(W)
W <- W %>% mutate_at(vars("Norm.Intensity.B", "Norm.Intensity.G", "Norm.Intensity.L", "Norm.Intensity.Q",    
                          "Norm.Intensity.E", "Norm.Intensity.J", "Norm.Intensity.O"), ~replace(., is.na(.), 0))
glimpse(W)

#combine again
FL_dE18.Bicc1.XYfilter.2 <- FL_dE18.Bicc1.XYfilter %>% filter(! FRID %in% FRIDS.possible.NApvals ) 
FL_dE18.Bicc1.XYfilter.2 <- bind_rows(FL_dE18.Bicc1.XYfilter.2, W)
FL_dE18.Bicc1.XYfilter.2 <- FL_dE18.Bicc1.XYfilter.2 %>% arrange(FRID)
glimpse(FL_dE18.Bicc1.XYfilter.2)


### LIMMA
# get data
DataforLimma <- FL_dE18.Bicc1.XYfilter.2[, c("Norm.Intensity.B", "Norm.Intensity.G", "Norm.Intensity.L", "Norm.Intensity.Q",    
                                             "Norm.Intensity.E", "Norm.Intensity.J", "Norm.Intensity.O")]
glimpse(DataforLimma)
colnames(DataforLimma)

design <- cbind(left=c(rep(1, 4), rep(0, 3)), right=c(rep(0, 4), rep(1, 3))); design

fit <- lmFit(DataforLimma, design) # Warning: Partial NA coefficients for 258 probe(s) 
cont.matrix <- makeContrasts( left.vs.right =  right - left, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)


Limma.p.value <- as.numeric(fit2$p.value)
glimpse(Limma.p.value )
Limma.adj.p.valueBH <- p.adjust(Limma.p.value, "BH") 
glimpse(Limma.adj.p.valueBH)
Limma.t <- as.numeric(fit2$t)
glimpse(Limma.t)

FL_dE18.Bicc1.XYfilter.2$Limma.p.value <- as.numeric(Limma.p.value)
FL_dE18.Bicc1.XYfilter.2$Limma.adj.p.valueBH <- as.numeric(Limma.adj.p.valueBH)
FL_dE18.Bicc1.XYfilter.2$Limma.t <- as.numeric(Limma.t)

FL_dE18.Bicc1.XYfilter.B <- FL_dE18.Bicc1.XYfilter.2
glimpse(FL_dE18.Bicc1.XYfilter.B)


### calculate FC
for(i in 1:nrow(FL_dE18.Bicc1.XYfilter.B)){
  print(i)
  left.vals <- as.numeric(FL_dE18.Bicc1.XYfilter.B[i,c("Norm.Intensity.B", "Norm.Intensity.G", "Norm.Intensity.L", "Norm.Intensity.Q")])
  right.vals <- as.numeric(FL_dE18.Bicc1.XYfilter.B[i,c("Norm.Intensity.E", "Norm.Intensity.J", "Norm.Intensity.O")])
  FC <- FR_FC_LOGData_2.3.17(control = left.vals, treatment = right.vals, LOGtype = 2)
  FL_dE18.Bicc1.XYfilter.B$FC_FR[i] <- FC
}
glimpse(FL_dE18.Bicc1.XYfilter.B)


#check numbers unique Sequence.window
FL_dE18.Bicc1.XYfilter.B %>% filter(Limma.p.value < 0.05) %>% select(Sequence.window) %>% distinct() %>% nrow() #3957 total sign
FL_dE18.Bicc1.XYfilter.B %>% filter(Limma.p.value < 0.05) %>% filter(FC_FR > 1) %>% select(Sequence.window) %>% distinct() %>% nrow() #2162 up
FL_dE18.Bicc1.XYfilter.B %>% filter(Limma.p.value < 0.05) %>% filter(FC_FR < 1) %>% select(Sequence.window) %>% distinct() %>% nrow() #1836 down

# filter na pvalues
table(is.na(FL_dE18.Bicc1.XYfilter.B$Limma.p.value))
FL_dE18.Bicc1.XYfilter.B <- FL_dE18.Bicc1.XYfilter.B %>% filter(!is.na(Limma.p.value))
glimpse(FL_dE18.Bicc1.XYfilter.B) #9364
table(is.na(FL_dE18.Bicc1.XYfilter.B$Limma.p.value))


# check +/- 10000 FCs
table(FL_dE18.Bicc1.XYfilter.B$FC_FR >= 10000) #428
table(FL_dE18.Bicc1.XYfilter.B$FC_FR <= -10000) #235
summary(FL_dE18.Bicc1.XYfilter.B$FC_FR)

#find second min or max and add +/-5 for +/- 10000 FCs
temp <- FL_dE18.Bicc1.XYfilter.B %>% filter(!FC_FR %in% c(10000, -10000)) %>% arrange(desc(FC_FR))
head(temp$FC_FR)
tail(temp$FC_FR)
max(temp$FC_FR)
min(temp$FC_FR)

FL_dE18.Bicc1.XYfilter.B$FC_FR <- replace(FL_dE18.Bicc1.XYfilter.B$FC_FR, FL_dE18.Bicc1.XYfilter.B$FC_FR==10000, max(temp$FC_FR)+5)
FL_dE18.Bicc1.XYfilter.B$FC_FR <- replace(FL_dE18.Bicc1.XYfilter.B$FC_FR, FL_dE18.Bicc1.XYfilter.B$FC_FR==-10000, min(temp$FC_FR)-5)
summary(FL_dE18.Bicc1.XYfilter.B$FC_FR)
glimpse(FL_dE18.Bicc1.XYfilter.B) #9364

## FL_dE18.Bicc1 reshape for RoKAI 
FL_dE18.Bicc1.RoKAI <- FL_dE18.Bicc1.XYfilter.B 

# make aa15.window
# take first entry sequence window if there are two entries
FL_dE18.Bicc1.RoKAI$aa15.window <- pbsapply(FL_dE18.Bicc1.RoKAI$Sequence.window, function(x) unlist(str_split(string =x, pattern = ";"))[1]  )
FL_dE18.Bicc1.RoKAI$aa15.window <- pbsapply(FL_dE18.Bicc1.RoKAI$aa15.window, function(x) make.aa.15.window.from.31.sequence.window(input = x)  )


glimpse(FL_dE18.Bicc1.RoKAI)

table(duplicated(FL_dE18.Bicc1.RoKAI$aa15.window))
which(duplicated(FL_dE18.Bicc1.RoKAI$aa15.window))


## if there are duplicate entries (e.g. because of _1_2_3) take the entry with the best p value
unique.15aa.windows <- unique(FL_dE18.Bicc1.RoKAI$aa15.window)
length(unique.15aa.windows) #8277
glimpse(FL_dE18.Bicc1.RoKAI) #9364


FL_dE18.Bicc1.RoKAI.2 <- c()
for (i in unique.15aa.windows){
  print(i)
  temp <- FL_dE18.Bicc1.RoKAI %>% filter(aa15.window %in% c(i))
  temp <- temp %>% arrange(Limma.p.value)
  temp <- temp[1,]
  FL_dE18.Bicc1.RoKAI.2 <- bind_rows(FL_dE18.Bicc1.RoKAI.2, temp)
}

glimpse(FL_dE18.Bicc1.RoKAI.2)
table(duplicated(FL_dE18.Bicc1.RoKAI.2$aa15.window))

#calculate log2(FC) and order
FL_dE18.Bicc1.RoKAI.2$log2FC <- log2(abs(FL_dE18.Bicc1.RoKAI.2$FC_FR ))
FL_dE18.Bicc1.RoKAI.2$sign.FC_FR <- sign(FL_dE18.Bicc1.RoKAI.2$FC_FR)
FL_dE18.Bicc1.RoKAI.2$log2FC <- sign(FL_dE18.Bicc1.RoKAI.2$FC_FR) * log2(abs(FL_dE18.Bicc1.RoKAI.2$FC_FR ))
glimpse(FL_dE18.Bicc1.RoKAI.2) 

FL_dE18.Bicc1.RoKAI.2 <- FL_dE18.Bicc1.RoKAI.2 %>% arrange(desc(log2FC))
head(FL_dE18.Bicc1.RoKAI.2$FC_FR)
tail(FL_dE18.Bicc1.RoKAI.2$FC_FR)

### prepare RoKAI input and save data
glimpse(FL_dE18.Bicc1.RoKAI.2) 
FL_dE18.Bicc1.RoKAI.2.23.11.21 <- FL_dE18.Bicc1.RoKAI.2 %>% select(Protein, Position, log2FC) 
table(is.na(FL_dE18.Bicc1.RoKAI.2.23.11.21$Protein))
FL_dE18.Bicc1.RoKAI.2.23.11.21 <- FL_dE18.Bicc1.RoKAI.2.23.11.21 %>% filter(!is.na(Protein))
glimpse(FL_dE18.Bicc1.RoKAI.2.23.11.21)
table(is.na(FL_dE18.Bicc1.RoKAI.2.23.11.21$Protein))
table(is.na(FL_dE18.Bicc1.RoKAI.2.23.11.21$Position))
table(is.na(FL_dE18.Bicc1.RoKAI.2.23.11.21$log2FC))
colnames(FL_dE18.Bicc1.RoKAI.2.23.11.21) <- c("Protein", "Position", "Quantification")

#write_csv(FL_dE18.Bicc1.RoKAI.2.23.11.21,"FL_dE18.Bicc1.RoKAI.2_23.11.21.csv", col_names=TRUE)
#save(FL_dE18.Bicc1.RoKAI  , file="FL_dE18.Bicc1.RoKAI.Rdata")

# perform analysis via https://rokai.io/
# tool settings was default: Reference Proteome: Uniprot Mouse / Fold changes: normalized / Kinase substrate dataset: PSP + Signor / RoKAI network: KS+PPI+SD+CoEv
# save RoKAI output kinase_table.csv and kinase_targets.csv and load these later

#################################################################################################################################
#################################################################################################################################

### combine RoKAI results kinase table
### filter for number of substrates >= 3 and FDR < 0.05 for each comparison

GFP_FL_rokai_kinase.table <- fread("kinase_table_ROKAI_DZ.cells.IMAC.GFP_FL_23.11.21.csv")
GFP_FL_rokai_kinase.table <- GFP_FL_rokai_kinase.table  %>% filter(NumSubs >= 3)
GFP_FL_rokai_kinase.table <- GFP_FL_rokai_kinase.table  %>% filter(FDR < 0.05)
colnames(GFP_FL_rokai_kinase.table) <- c("KinaseID","Name","Gene",  paste0("GFP_FL_", colnames(GFP_FL_rokai_kinase.table[, c(4:9)])))
glimpse(GFP_FL_rokai_kinase.table)


GFP_dE18_rokai_kinase.table <- fread("kinase_table_ROKAI_DZ.cells.IMAC.GFP_dE18_23.11.21.csv")
GFP_dE18_rokai_kinase.table <- GFP_dE18_rokai_kinase.table  %>% filter(NumSubs >= 3)
GFP_dE18_rokai_kinase.table <- GFP_dE18_rokai_kinase.table  %>% filter(FDR < 0.05)
colnames(GFP_dE18_rokai_kinase.table) <- c("KinaseID","Name","Gene",  paste0("GFP_dE18_", colnames(GFP_dE18_rokai_kinase.table[, c(4:9)])))
glimpse(GFP_dE18_rokai_kinase.table)


GFP_FL.Bicc1_rokai_kinase.table <- fread("kinase_table_ROKAI_DZ.cells.IMAC.GFP_FL.Bicc1_23.11.21.csv")
GFP_FL.Bicc1_rokai_kinase.table <- GFP_FL.Bicc1_rokai_kinase.table  %>% filter(NumSubs >= 3)
GFP_FL.Bicc1_rokai_kinase.table <- GFP_FL.Bicc1_rokai_kinase.table  %>% filter(FDR < 0.05)
colnames(GFP_FL.Bicc1_rokai_kinase.table) <- c("KinaseID","Name","Gene",  paste0("GFP_FL.Bicc1_", colnames(GFP_FL.Bicc1_rokai_kinase.table[, c(4:9)])))
glimpse(GFP_FL.Bicc1_rokai_kinase.table)


GFP_dE18.Bicc1_rokai_kinase.table <- fread("kinase_table_ROKAI_DZ.cells.IMAC.GFP_dE18.Bicc1_23.11.21.csv")
GFP_dE18.Bicc1_rokai_kinase.table <- GFP_dE18.Bicc1_rokai_kinase.table  %>% filter(NumSubs >= 3)
GFP_dE18.Bicc1_rokai_kinase.table <- GFP_dE18.Bicc1_rokai_kinase.table  %>% filter(FDR < 0.05)
colnames(GFP_dE18.Bicc1_rokai_kinase.table) <- c("KinaseID","Name","Gene",  paste0("GFP_dE18.Bicc1_", colnames(GFP_dE18.Bicc1_rokai_kinase.table[, c(4:9)])))
glimpse(GFP_dE18.Bicc1_rokai_kinase.table)


FL_dE18_rokai_kinase.table <- fread("kinase_table_ROKAI_DZ.cells.IMAC.FL_dE18_23.11.21.csv")
FL_dE18_rokai_kinase.table <- FL_dE18_rokai_kinase.table  %>% filter(NumSubs >= 3)
FL_dE18_rokai_kinase.table <- FL_dE18_rokai_kinase.table  %>% filter(FDR < 0.05)
colnames(FL_dE18_rokai_kinase.table) <- c("KinaseID","Name","Gene",  paste0("FL_dE18_", colnames(FL_dE18_rokai_kinase.table[, c(4:9)])))
glimpse(FL_dE18_rokai_kinase.table)


FL_FL.Bicc1_rokai_kinase.table <- fread("kinase_table_ROKAI_DZ.cells.IMAC.FL_FL.Bicc1_23.11.21.csv")
FL_FL.Bicc1_rokai_kinase.table <- FL_FL.Bicc1_rokai_kinase.table  %>% filter(NumSubs >= 3)
FL_FL.Bicc1_rokai_kinase.table <- FL_FL.Bicc1_rokai_kinase.table  %>% filter(FDR < 0.05)
colnames(FL_FL.Bicc1_rokai_kinase.table) <- c("KinaseID","Name","Gene",  paste0("FL_FL.Bicc1_", colnames(FL_FL.Bicc1_rokai_kinase.table[, c(4:9)])))
glimpse(FL_FL.Bicc1_rokai_kinase.table )


FL_dE18.Bicc1_rokai_kinase.table <- fread("kinase_table_ROKAI_DZ.cells.IMAC.FL_dE18.Bicc1_23.11.21.csv")
FL_dE18.Bicc1_rokai_kinase.table <- FL_dE18.Bicc1_rokai_kinase.table  %>% filter(NumSubs >= 3)
FL_dE18.Bicc1_rokai_kinase.table <- FL_dE18.Bicc1_rokai_kinase.table  %>% filter(FDR< 0.05)
colnames(FL_dE18.Bicc1_rokai_kinase.table) <- c("KinaseID","Name","Gene",  paste0("FL_dE18.Bicc1_", colnames(FL_dE18.Bicc1_rokai_kinase.table[, c(4:9)])))
glimpse(FL_dE18.Bicc1_rokai_kinase.table )



###get all KinaseIDs corresponding Names and Genes and remove duplicates
comb.ROKAI.KinaseID <- bind_rows(GFP_FL_rokai_kinase.table %>% select(KinaseID, Name , Gene),
                                 GFP_dE18_rokai_kinase.table %>% select(KinaseID, Name , Gene),
                                 GFP_FL.Bicc1_rokai_kinase.table %>% select(KinaseID, Name , Gene),
                                 GFP_dE18.Bicc1_rokai_kinase.table %>% select(KinaseID, Name , Gene),
                                 
                                 FL_dE18_rokai_kinase.table %>% select(KinaseID, Name , Gene),
                                 FL_FL.Bicc1_rokai_kinase.table %>% select(KinaseID, Name , Gene),
                                 FL_dE18.Bicc1_rokai_kinase.table %>% select(KinaseID, Name , Gene)) %>% distinct()
glimpse(comb.ROKAI.KinaseID)
comb.ROKAI.KinaseID$Gene

#merge data 
comb.ROKAI.result <- merge(comb.ROKAI.KinaseID, GFP_FL_rokai_kinase.table[, c(1, 4:9)], by.x = "KinaseID", by.y = "KinaseID", all.x = T, all.y = T, sort = F)
comb.ROKAI.result <- merge(comb.ROKAI.result, GFP_dE18_rokai_kinase.table[, c(1, 4:9)], by.x = "KinaseID", by.y = "KinaseID", all.x = T, all.y = T, sort = F)
comb.ROKAI.result <- merge(comb.ROKAI.result, GFP_FL.Bicc1_rokai_kinase.table[, c(1, 4:9)], by.x = "KinaseID", by.y = "KinaseID", all.x = T, all.y = T, sort = F)
comb.ROKAI.result <- merge(comb.ROKAI.result, GFP_dE18.Bicc1_rokai_kinase.table[, c(1, 4:9)], by.x = "KinaseID", by.y = "KinaseID", all.x = T, all.y = T, sort = F)
comb.ROKAI.result <- merge(comb.ROKAI.result, FL_dE18_rokai_kinase.table[, c(1, 4:9)], by.x = "KinaseID", by.y = "KinaseID", all.x = T, all.y = T, sort = F)
comb.ROKAI.result <- merge(comb.ROKAI.result, FL_FL.Bicc1_rokai_kinase.table[, c(1, 4:9)], by.x = "KinaseID", by.y = "KinaseID", all.x = T, all.y = T, sort = F)
comb.ROKAI.result <- merge(comb.ROKAI.result,  FL_dE18.Bicc1_rokai_kinase.table[, c(1, 4:9)], by.x = "KinaseID", by.y = "KinaseID", all.x = T, all.y = T, sort = F)

glimpse(comb.ROKAI.result)

#plot kinase activity scores
plot.comb.ROKAI.result <- comb.ROKAI.result %>% arrange(desc(FL_dE18_ZScore))  %>% select(Gene,  GFP_FL_ZScore, GFP_dE18_ZScore, GFP_FL.Bicc1_ZScore, GFP_dE18.Bicc1_ZScore, FL_dE18_ZScore,FL_FL.Bicc1_ZScore,FL_dE18.Bicc1_ZScore   )
plot.comb.ROKAI.result$sum.z.scores <- apply(plot.comb.ROKAI.result[, c("GFP_FL_ZScore", "GFP_dE18_ZScore", "GFP_FL.Bicc1_ZScore", "GFP_dE18.Bicc1_ZScore", "FL_dE18_ZScore", "FL_FL.Bicc1_ZScore","FL_dE18.Bicc1_ZScore")], 1, function(x) sum(x, na.rm=T))
plot.comb.ROKAI.result <- plot.comb.ROKAI.result %>%  arrange(desc(sum.z.scores)) %>% select(-sum.z.scores)
glimpse(plot.comb.ROKAI.result)

order.plot.comb.ROKAI.result <- plot.comb.ROKAI.result$Gene
order.plot.comb.ROKAI.result

melted.plot.comb.ROKAI.result <- reshape2::melt(plot.comb.ROKAI.result) 
glimpse(melted.plot.comb.ROKAI.result)

#plot overview heatmap
ggplot(melted.plot.comb.ROKAI.result  , aes(x=variable, y=Gene, fill=value)) + 
  geom_tile(color = "white") + 
  scale_fill_viridis( option = "turbo", na.value = "grey60", name="z.score") + 
  coord_equal()+
  scale_y_discrete(limits= rev(order.plot.comb.ROKAI.result), expand=c(0,0))+
  theme(legend.position="right", legend.justification = "center")+
  ggtitle("IMAC cells - RoKAI Kinase activity z-score ") +
  theme(plot.title = element_text(hjust = 0.5, face="bold", size=16)) + #center title
  theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1))+  ## Vertical text on x axis
  xlab(NULL) + 
  ylab(NULL) +
  theme(axis.text.x = element_text(face= "plain", colour="black", size=10)) 


################################################################################################################################################################
################################################################################################################################################################

### combine and visualize kinase targets

glimpse(DZ.IMAC.MBR.ON_expanded_log2_norm.class1.PS)

#define kinases of interest
kinases.of.interest <- order.plot.comb.ROKAI.result
length(kinases.of.interest)
kinases.of.interest

### load Rokai Kinase target results
GFP_FL_rokai_kinase.targets <- fread("kinase_targets_ROKAI_DZ.cells.IMAC.GFP_FL_23.11.21.csv")
colnames(GFP_FL_rokai_kinase.targets) <- c("KinID","KinName","KinGene","SubsProtein","SubsGene","Position","DataSource","Flanking" ,  paste0("GFP_FL_", colnames(GFP_FL_rokai_kinase.targets[, c(9:12)])))
glimpse(GFP_FL_rokai_kinase.targets)


GFP_dE18_rokai_kinase.targets <- fread("kinase_targetsROKAI_DZ.cells.IMAC.GFP_dE18_23.11.21.csv")
colnames(GFP_dE18_rokai_kinase.targets) <- c("KinID","KinName","KinGene","SubsProtein","SubsGene","Position","DataSource","Flanking" ,  paste0("GFP_dE18_", colnames(GFP_dE18_rokai_kinase.targets[, c(9:12)])))
glimpse(GFP_dE18_rokai_kinase.targets)


GFP_FL.Bicc1_rokai_kinase.targets <- fread("kinase_targets_ROKAI_DZ.cells.IMAC.GFP_FL.Bicc1_23.11.21.csv")
colnames(GFP_FL.Bicc1_rokai_kinase.targets) <- c("KinID","KinName","KinGene","SubsProtein","SubsGene","Position","DataSource","Flanking" ,  paste0("GFP_FL.Bicc1_", colnames(GFP_FL.Bicc1_rokai_kinase.targets[, c(9:12)])))
glimpse(GFP_FL.Bicc1_rokai_kinase.targets)


GFP_dE18.Bicc1_rokai_kinase.targets <- fread("kinase_targets_ROKAI_DZ.cells.IMAC.GFP_dE18.Bicc1_23.11.21.csv")
colnames(GFP_dE18.Bicc1_rokai_kinase.targets) <- c("KinID","KinName","KinGene","SubsProtein","SubsGene","Position","DataSource","Flanking" ,  paste0("GFP_dE18.Bicc1_", colnames(GFP_dE18.Bicc1_rokai_kinase.targets[, c(9:12)])))
glimpse(GFP_dE18.Bicc1_rokai_kinase.targets)


FL_dE18_rokai_kinase.targets <- fread("kinase_targets_ROKAI_DZ.cells.IMAC.FL_dE18_23.11.21.csv")
colnames(FL_dE18_rokai_kinase.targets) <- c("KinID","KinName","KinGene","SubsProtein","SubsGene","Position","DataSource","Flanking" ,  paste0("FL_dE18_", colnames(FL_dE18_rokai_kinase.targets[, c(9:12)])))
glimpse(FL_dE18_rokai_kinase.targets)


FL_FL.Bicc1_rokai_kinase.targets <- fread("kinase_targets_ROKAI_DZ.cells.IMAC.FL_FL.Bicc1_23.11.21.csv")
colnames(FL_FL.Bicc1_rokai_kinase.targets) <- c("KinID","KinName","KinGene","SubsProtein","SubsGene","Position","DataSource","Flanking" ,  paste0("FL_FL.Bicc1_", colnames(FL_FL.Bicc1_rokai_kinase.targets[, c(9:12)])))
glimpse(FL_FL.Bicc1_rokai_kinase.targets )


FL_dE18.Bicc1_rokai_kinase.targets <- fread("kinase_targets_ROKAI_DZ.cells.IMAC.FL_dE18.Bicc1_23.11.21.csv")
colnames(FL_dE18.Bicc1_rokai_kinase.targets) <- c("KinID","KinName","KinGene","SubsProtein","SubsGene","Position","DataSource","Flanking" ,  paste0("FL_dE18.Bicc1_", colnames(FL_dE18.Bicc1_rokai_kinase.targets[, c(9:12)])))
glimpse(FL_dE18.Bicc1_rokai_kinase.targets )



###get all KinaseIDs, targtes corresponding Names and Genes and remove duplicates
comb.ROKAI.targets <- bind_rows(GFP_FL_rokai_kinase.targets %>% select(KinID, KinName, KinGene, SubsProtein, SubsGene, Position, DataSource, Flanking),
                                GFP_dE18_rokai_kinase.targets %>% select(KinID, KinName, KinGene, SubsProtein, SubsGene, Position, DataSource, Flanking),
                                GFP_FL.Bicc1_rokai_kinase.targets %>% select(KinID, KinName, KinGene, SubsProtein, SubsGene, Position, DataSource, Flanking),
                                GFP_dE18.Bicc1_rokai_kinase.targets %>% select(KinID, KinName, KinGene, SubsProtein, SubsGene, Position, DataSource, Flanking),
                                
                                FL_dE18_rokai_kinase.targets %>% select(KinID, KinName, KinGene, SubsProtein, SubsGene, Position, DataSource, Flanking),
                                FL_FL.Bicc1_rokai_kinase.targets %>% select(KinID, KinName, KinGene, SubsProtein, SubsGene, Position, DataSource, Flanking),
                                FL_dE18.Bicc1_rokai_kinase.targets %>% select(KinID, KinName, KinGene, SubsProtein, SubsGene, Position, DataSource, Flanking)) %>% distinct()
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

#load("GFP_FL.RoKAI.Rdata")
glimpse(GFP_FL.RoKAI) 
GFP_FL_2group.filter.sign.FC <- GFP_FL.RoKAI %>% filter(Limma.p.value < 0.05) %>% filter(FC_FR >=1.5 | FC_FR <= -1.5) %>% select(ID.FR.all.C1.PS) %>% pull()
glimpse(GFP_FL_2group.filter.sign.FC)

#load("GFP_dE18.RoKAI.Rdata")
glimpse(GFP_dE18.RoKAI)
GFP_dE18_2group.filter.sign.FC <- GFP_dE18.RoKAI  %>% filter(Limma.p.value < 0.05) %>% filter(FC_FR >=1.5 | FC_FR <= -1.5) %>% select(ID.FR.all.C1.PS) %>% pull()
glimpse(GFP_dE18_2group.filter.sign.FC)

#load("GFP_FL.Bicc1.RoKAI.Rdata")
glimpse(GFP_FL.Bicc1.RoKAI)
GFP_FL.Bicc1_2group.filter.sign.FC <- GFP_FL.Bicc1.RoKAI %>% filter(Limma.p.value < 0.05) %>% filter(FC_FR >=1.5 | FC_FR <= -1.5) %>% select(ID.FR.all.C1.PS) %>% pull()
glimpse(GFP_FL.Bicc1_2group.filter.sign.FC)

#load("GFP_dE18.Bicc1.RoKAI.Rdata")
glimpse(GFP_dE18.Bicc1.RoKAI)
GFP_dE18.Bicc1_2group.filter.sign.FC <- GFP_dE18.Bicc1.RoKAI %>% filter(Limma.p.value < 0.05) %>% filter(FC_FR >=1.5 | FC_FR <= -1.5) %>% select(ID.FR.all.C1.PS) %>% pull()
glimpse(GFP_dE18.Bicc1_2group.filter.sign.FC)

#load("FL_dE18.RoKAI.Rdata")
glimpse(FL_dE18.RoKAI)
FL_dE18_2group.filter.sign.FC <- FL_dE18.RoKAI %>% filter(Limma.p.value < 0.05) %>% filter(FC_FR >=1.5 | FC_FR <= -1.5) %>% select(ID.FR.all.C1.PS) %>% pull()
glimpse(FL_dE18_2group.filter.sign.FC)

#load("FL_FL.Bicc1.RoKAI.Rdata")
glimpse(FL_FL.Bicc1.RoKAI)
FL_FL.Bicc1_2group.filter.sign.FC <- FL_FL.Bicc1.RoKAI %>% filter(Limma.p.value < 0.05) %>% filter(FC_FR >=1.5 | FC_FR <= -1.5) %>% select(ID.FR.all.C1.PS) %>% pull()
glimpse(FL_FL.Bicc1_2group.filter.sign.FC)

#load("FL_dE18.Bicc1.RoKAI.Rdata")
glimpse(FL_dE18.Bicc1.RoKAI)
FL_dE18.Bicc1_2group.filter.sign.FC <- FL_dE18.Bicc1.RoKAI  %>% filter(Limma.p.value < 0.05) %>% filter(FC_FR >=1.5 | FC_FR <= -1.5) %>% select(ID.FR.all.C1.PS) %>% pull()
glimpse(FL_dE18.Bicc1_2group.filter.sign.FC)


### get unique candidate IDs 
all.2group.filter.sign.FC.C1.IDs <- unique( c(GFP_FL_2group.filter.sign.FC, GFP_dE18_2group.filter.sign.FC, GFP_FL.Bicc1_2group.filter.sign.FC, GFP_dE18.Bicc1_2group.filter.sign.FC,
                                              FL_dE18_2group.filter.sign.FC, FL_FL.Bicc1_2group.filter.sign.FC, FL_dE18.Bicc1_2group.filter.sign.FC) )
glimpse(all.2group.filter.sign.FC.C1.IDs)




###get all the intensity data for all samples and the selected sites (all.2group.filter.sign.FC.C1.IDs)
glimpse(DZ.IMAC.MBR.ON_expanded_log2_norm.class1.PS)

DZ.IMAC.MBR.ON_expanded_log2_norm.class1.PS_filter.sign.FC <- DZ.IMAC.MBR.ON_expanded_log2_norm.class1.PS %>% filter(ID.FR.all.C1.PS %in% all.2group.filter.sign.FC.C1.IDs)
glimpse(DZ.IMAC.MBR.ON_expanded_log2_norm.class1.PS_filter.sign.FC)

# make aa15.window
# take first entry sequence window if there are two entries
DZ.IMAC.MBR.ON_expanded_log2_norm.class1.PS_filter.sign.FC$aa15.window <- pbsapply(DZ.IMAC.MBR.ON_expanded_log2_norm.class1.PS_filter.sign.FC$Sequence.window, function(x) unlist(str_split(string =x, pattern = ";"))[1]  )
DZ.IMAC.MBR.ON_expanded_log2_norm.class1.PS_filter.sign.FC$aa15.window <- pbsapply(DZ.IMAC.MBR.ON_expanded_log2_norm.class1.PS_filter.sign.FC$aa15.window, function(x) make.aa.15.window.from.31.sequence.window(input = x)  )
glimpse(DZ.IMAC.MBR.ON_expanded_log2_norm.class1.PS_filter.sign.FC)


# filter for the combined RoKAI target sites
glimpse(comb.ROKAI.targets_2_unique.flanking)

DZ.IMAC.MBR.ON_expanded_log2_norm.class1.PS_filter.sign.FC_rokai.filter <- DZ.IMAC.MBR.ON_expanded_log2_norm.class1.PS_filter.sign.FC %>% filter(aa15.window %in% comb.ROKAI.targets_2_unique.flanking)
glimpse(DZ.IMAC.MBR.ON_expanded_log2_norm.class1.PS_filter.sign.FC_rokai.filter)



#merge information on upstream kinase
glimpse(comb.ROKAI.targets_2)

DZ.IMAC.MBR.ON_expanded_log2_norm.class1.PS_filter.sign.FC_rokai.filter_merge.rokai.upstream.Kins <- merge(DZ.IMAC.MBR.ON_expanded_log2_norm.class1.PS_filter.sign.FC_rokai.filter, comb.ROKAI.targets_2, by.x = "aa15.window", by.y = "Flanking", all.x = T, all.y = F, sort = F)
glimpse(DZ.IMAC.MBR.ON_expanded_log2_norm.class1.PS_filter.sign.FC_rokai.filter_merge.rokai.upstream.Kins)

#prepare a unique rowname to be able to plot everything
DZ.IMAC.MBR.ON_expanded_log2_norm.class1.PS_filter.sign.FC_rokai.filter_merge.rokai.upstream.Kins$row_Gene.Name_AApos <- paste0(paste0(rownames(DZ.IMAC.MBR.ON_expanded_log2_norm.class1.PS_filter.sign.FC_rokai.filter_merge.rokai.upstream.Kins), "_", DZ.IMAC.MBR.ON_expanded_log2_norm.class1.PS_filter.sign.FC_rokai.filter_merge.rokai.upstream.Kins$Gene.Name_AApos))
glimpse(DZ.IMAC.MBR.ON_expanded_log2_norm.class1.PS_filter.sign.FC_rokai.filter_merge.rokai.upstream.Kins)


### prepare heatmap facet plot for results combined RoKAI #################################################################################


facet.hm.plot <- DZ.IMAC.MBR.ON_expanded_log2_norm.class1.PS_filter.sign.FC_rokai.filter_merge.rokai.upstream.Kins %>% select(
  "Gene.Name_AApos", 
  "KinGene",
  "Norm.Intensity.A", "Norm.Intensity.F", "Norm.Intensity.K", "Norm.Intensity.P", 
  "Norm.Intensity.B", "Norm.Intensity.G", "Norm.Intensity.L", "Norm.Intensity.Q",    
  "Norm.Intensity.C", "Norm.Intensity.H", "Norm.Intensity.M", "Norm.Intensity.R", 
  "Norm.Intensity.D", "Norm.Intensity.I", "Norm.Intensity.N", "Norm.Intensity.S",    
  "Norm.Intensity.E", "Norm.Intensity.J", "Norm.Intensity.O"
)

glimpse(facet.hm.plot)

#get overview how many observations per upstream Kinase  (some sites duplicated per kinase term)
print(facet.hm.plot %>% group_by(KinGene) %>% summarize(count = n()) %>% arrange(desc(count)), n=38)

#define order via number of substrates per upstream kinase
order.upstreamKin.via.num.subs <- facet.hm.plot %>% group_by(KinGene) %>% summarize(count = n()) %>% arrange(desc(count)) %>% select(KinGene) %>% pull()
order.upstreamKin.via.num.subs
order.upstreamKin.via.num.subs[1] #the kinase with the highest substrates here: Mapk1

# remove duplicate PS per upstream Kinase start with Kinase that has highest substrate number
glimpse(facet.hm.plot)
facet.hm.plot_2 <- facet.hm.plot %>% filter(KinGene %in% c(order.upstreamKin.via.num.subs[1]))
glimpse(facet.hm.plot_2)
temp.substrates.start <- facet.hm.plot %>% filter(KinGene %in% c(order.upstreamKin.via.num.subs[1])) %>% select(Gene.Name_AApos) %>% pull()
temp.substrates.start

for(i in order.upstreamKin.via.num.subs[2:length(order.upstreamKin.via.num.subs)]){
  print(i)
  temp.substrates.next.kin <- facet.hm.plot %>% filter(KinGene %in% c(i)) %>% select(Gene.Name_AApos) %>% pull()
  temp.substrates.next.kin
  temp.substrates.next.kin.include <- setdiff(temp.substrates.next.kin, temp.substrates.start)
  temp.substrates.next.kin.include
  
  temp.df.include <- facet.hm.plot %>% filter(KinGene %in% c(i)) %>% filter(Gene.Name_AApos %in% temp.substrates.next.kin.include)
  
  facet.hm.plot_2 <- bind_rows(facet.hm.plot_2, temp.df.include)
  
  temp.substrates.next.kin.exclude <- setdiff( temp.substrates.start, temp.substrates.next.kin)
  temp.substrates.next.kin.exclude
  
  temp.substrates.start <- c(temp.substrates.start, temp.substrates.next.kin.include)
  temp.substrates.start
}
glimpse(facet.hm.plot_2)


# change column names and check
colnames(facet.hm.plot_2) <- c("Gene.Name_AApos", "KinGene" ,"GFP.1","GFP.2","GFP.3","GFP.4",  "FL.1","FL.2","FL.3","FL.4",  "dE18.1","dE18.2","dE18.3","dE18.4",  "FL-Bicc1.1","FL-Bicc1.2","FL-Bicc1.3","FL-Bicc1.4",   "dE18-Bicc1.1","dE18-Bicc1.2","dE18-Bicc1.3")
glimpse(facet.hm.plot_2)
table(duplicated(facet.hm.plot_2$Gene.Name_AApos))


###for better visualization collapse _1/_2/_3 sites
###collapsing _1/_2/_3 via using sum
# note: sum is based on _1/_2/_3 intensities in the dataframe at hand => unfiltered data might contain more _1/_2/_3 PS that are not included in the sum here
# unlog, collapse, relog

facet.hm.plot_2$Gene.Name_AApos.collapsed <- pbsapply(facet.hm.plot_2$Gene.Name_AApos, function(x) str_remove_all(string=x, pattern="x_\\d") )
glimpse(facet.hm.plot_2)
# unlog, do collapsing
facet.hm.plot_2[, c("GFP.1","GFP.2","GFP.3","GFP.4","FL.1","FL.2","FL.3","FL.4","dE18.1","dE18.2","dE18.3","dE18.4","FL-Bicc1.1","FL-Bicc1.2","FL-Bicc1.3","FL-Bicc1.4","dE18-Bicc1.1","dE18-Bicc1.2","dE18-Bicc1.3")] <- 2^facet.hm.plot_2[, c("GFP.1","GFP.2","GFP.3","GFP.4","FL.1","FL.2","FL.3","FL.4","dE18.1","dE18.2","dE18.3","dE18.4","FL-Bicc1.1","FL-Bicc1.2","FL-Bicc1.3","FL-Bicc1.4","dE18-Bicc1.1","dE18-Bicc1.2","dE18-Bicc1.3")] 
glimpse(facet.hm.plot_2)
#collapse
facet.hm.plot_2_B <-facet.hm.plot_2 %>% group_by(Gene.Name_AApos.collapsed) %>% summarise_at(c("GFP.1","GFP.2","GFP.3","GFP.4","FL.1","FL.2","FL.3","FL.4","dE18.1","dE18.2","dE18.3","dE18.4","FL-Bicc1.1","FL-Bicc1.2","FL-Bicc1.3","FL-Bicc1.4","dE18-Bicc1.1","dE18-Bicc1.2","dE18-Bicc1.3"), sum, na.rm = TRUE)
glimpse(facet.hm.plot_2_B)
facet.hm.plot_2_B[facet.hm.plot_2_B == 0] <- NA
#relog
facet.hm.plot_2_B[c("GFP.1","GFP.2","GFP.3","GFP.4","FL.1","FL.2","FL.3","FL.4","dE18.1","dE18.2","dE18.3","dE18.4","FL-Bicc1.1","FL-Bicc1.2","FL-Bicc1.3","FL-Bicc1.4","dE18-Bicc1.1","dE18-Bicc1.2","dE18-Bicc1.3")] <- log2(facet.hm.plot_2_B[c("GFP.1","GFP.2","GFP.3","GFP.4","FL.1","FL.2","FL.3","FL.4","dE18.1","dE18.2","dE18.3","dE18.4","FL-Bicc1.1","FL-Bicc1.2","FL-Bicc1.3","FL-Bicc1.4","dE18-Bicc1.1","dE18-Bicc1.2","dE18-Bicc1.3")] )
glimpse(facet.hm.plot_2_B)
#get KinGene column back
facet.hm.plot_2_B <- merge(facet.hm.plot_2_B, facet.hm.plot_2[, c("Gene.Name_AApos.collapsed", "KinGene")], by.x="Gene.Name_AApos.collapsed", by.y="Gene.Name_AApos.collapsed", all.x=T, all.y=F, sort=F)
glimpse(facet.hm.plot_2_B)

#make new data frame
facet.hm.plot_3 <- facet.hm.plot_2_B

### reorder data frame (collapsed sites) and calculate FC for FL vs. dE18
second.order.upstreamKin.via.num.subs <- facet.hm.plot_3 %>% group_by(KinGene) %>% summarize(count = n()) %>% arrange(desc(count)) %>% select(KinGene) %>% pull()
second.order.upstreamKin.via.num.subs

#calculate FC for FL vs. dE18
glimpse(facet.hm.plot_3)
for(i in 1:nrow(facet.hm.plot_3)){
  print(i)
  leftvals <- as.numeric(facet.hm.plot_3[i,c("FL.1", "FL.2", "FL.3", "FL.4")])
  rightvals <- as.numeric(facet.hm.plot_3[i,c("dE18.1", "dE18.2", "dE18.3", "dE18.4" )])
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
facet.hm.plot_3 <- facet.hm.plot_3 %>% select(Gene.Name_AApos.collapsed, KinGene, "GFP.1","GFP.2","GFP.3","GFP.4","FL.1","FL.2","FL.3","FL.4","dE18.1","dE18.2", "dE18.3","dE18.4","FL-Bicc1.1","FL-Bicc1.2","FL-Bicc1.3","FL-Bicc1.4","dE18-Bicc1.1","dE18-Bicc1.2","dE18-Bicc1.3")
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
#stop define function

facet.hm.plot_3_Z <- facet.hm.plot_3 %>% select("GFP.1","GFP.2","GFP.3","GFP.4","FL.1","FL.2","FL.3","FL.4","dE18.1","dE18.2", "dE18.3","dE18.4","FL-Bicc1.1","FL-Bicc1.2","FL-Bicc1.3","FL-Bicc1.4","dE18-Bicc1.1","dE18-Bicc1.2","dE18-Bicc1.3")
glimpse(facet.hm.plot_3_Z)

facet.hm.plot_3_Z <- t(apply(facet.hm.plot_3_Z, 1, function(x) zscorescalebaseR(x = x)))
glimpse(facet.hm.plot_3_Z)
sd(facet.hm.plot_3_Z[1,], na.rm = T) #should be 1

facet.hm.plot_3_Z <- as.data.frame(facet.hm.plot_3_Z)
glimpse(facet.hm.plot_3_Z)
facet.hm.plot_3_Z$Gene.Name_AApos.collapsed <- facet.hm.plot_3$Gene.Name_AApos.collapsed
facet.hm.plot_3_Z$KinGene <- facet.hm.plot_3$KinGene
facet.hm.plot_3_Z <- facet.hm.plot_3_Z %>% select(Gene.Name_AApos.collapsed, KinGene, "GFP.1","GFP.2","GFP.3","GFP.4","FL.1","FL.2","FL.3","FL.4","dE18.1","dE18.2", "dE18.3","dE18.4","FL-Bicc1.1","FL-Bicc1.2","FL-Bicc1.3","FL-Bicc1.4","dE18-Bicc1.1","dE18-Bicc1.2","dE18-Bicc1.3")
glimpse(facet.hm.plot_3_Z)



# melt data frame for plot and change column names
melted.facet.hm.plot <- reshape2::melt(facet.hm.plot_3_Z)

colnames(melted.facet.hm.plot ) <- c("PG.Genes.3", "pathway",    "variable",   "log2.int")    
glimpse(melted.facet.hm.plot)

#plot all data found
ggplot()+
  geom_tile(data=melted.facet.hm.plot, aes(y = fct_inorder(PG.Genes.3), x = variable, fill = log2.int))+ 
  scale_x_discrete(limits=c("GFP.1","GFP.2","GFP.3","GFP.4",  "FL.1","FL.2","FL.3","FL.4",  "dE18.1","dE18.2","dE18.3","dE18.4",  "FL-Bicc1.1","FL-Bicc1.2","FL-Bicc1.3","FL-Bicc1.4",   "dE18-Bicc1.1","dE18-Bicc1.2","dE18-Bicc1.3"))+
  
  scale_fill_scico(palette = 'vikO',name="row z-score \n log2(int.)", na.value = "grey50")+ #
  
  xlab(NULL) + 
  ylab(NULL)+
  theme(strip.background=element_rect(colour="transparent", fill="transparent"))+ # facet box color
  theme(plot.title = element_text(hjust = 0.5, face="bold"))+  #center title
  theme(legend.position = "right")+
  theme(axis.text.x = element_text(size = 10)) +
  theme(axis.text.y = element_text(size = 12))+
  facet_grid(fct_inorder(pathway) ~ ., scales = "free", space = "free")+
  theme(strip.text.y = element_text(angle=0, face="bold", colour="black", size=10, lineheight = 5))+ # facet box text size
  theme(axis.text.x  = element_text(angle=90, hjust = 1.0, vjust = 0.5, face="bold"))+
  theme(axis.text.y = element_text(size = 7))+
  geom_vline(xintercept = c(4.5, 8.5, 12.5, 16.5) , size=1.5, linetype = "solid", color="white")

#ggsave("DZ.cells.IMAC.global.phospho.active.kinase.sites.no.duplicates.FC1.5.pdf", useDingbats=FALSE,  width = 15, height = 60, units = "cm") #, width = 15, height = 25, units = "cm"



### visualize only selected site candidates, group/adjust names upstream kinases

## define selected sites
DZ.candidates.29.11.21_B <- tibble(Gene.Name_AApos.collapsed = c("Pak1_T212",                        "Pxn_S83",                          "Pxn_S126",                         "Pxn_S130",                        
                                                                 "Map2k1;Map2k2_S222;226",           "Map2k1_T292",                      "Map2k1_T386",                      "Mapk1_T183",                      
                                                                 "Mapk3_T203",                       "Mapk3_Y205",                       "Erf_T529",                         "Jun_S63",                         
                                                                 "Junb_S256",                        "Gsk3b_S9",                         "Map1b_S1260",                      "Map1b_T1784",                     
                                                                 "Akt1s1_T247",                      "Acly_S455",                        "Bad_S112",                         "Tbc1d1_S621",                     
                                                                 "Tbc1d1_T590",                      "Foxo3_S252",                       "Rps6ka1_S352",                     "Rps6ka1;Rps6ka3;Rps6ka2_T562;577",
                                                                 "Rps6ka3_S386",                     "Rps6_S235",                        "Rps6_S236",                        "Rps6_S240",                       
                                                                 "Rps6_S244",                        "Eif4b_S422",                       "Eif4ebp1_T36",                     "Eif4ebp1_S64",                    
                                                                 "Eif4ebp1_T69",                     "Eif4ebp2_T37",                     "Cad_S1859",                       "Larp1_S743",                      
                                                                 "Larp1_S751"),
                                   Group.DZ = c("MAPK signaling", "MAPK signaling", "MAPK signaling", "MAPK signaling", "MAPK signaling", "MAPK signaling", "MAPK signaling", "MAPK signaling",
                                                "MAPK signaling", "MAPK signaling", "MAPK signaling", "MAPK signaling", "MAPK signaling", "MAPK signaling", "MAPK signaling", "MAPK signaling",
                                                "AKT signaling",  "AKT signaling",  "AKT signaling",  "AKT signaling",  "AKT signaling",  "AKT signaling",  "mTOR signaling", "mTOR signaling",
                                                "mTOR signaling", "mTOR signaling", "mTOR signaling", "mTOR signaling", "mTOR signaling", "mTOR signaling", "mTOR signaling", "mTOR signaling",
                                                "mTOR signaling", "mTOR signaling", "mTOR signaling", "mTOR signaling", "mTOR signaling"),
                                   Order.DZ = c(1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16,  1,  2,  3,  4,  5,  6,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15),
                                   Order.DZ.rev = c(16, 15, 14, 13, 12, 11, 10,  9,  8,  7,  6,  5,  4,  3,  2,  1,    6, 5, 4, 3, 2, 1,       15, 14, 13, 12, 11, 10,  9,  8,  7,  6,  5,  4,  3,  2,  1),
                                   Order.FR = c(1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37),
                                   Order.FR.rev = c(37, 36, 35, 34, 33, 32, 31, 30, 29, 28, 27, 26, 25, 24, 23, 22, 21, 20, 19, 18, 17, 16, 15, 14, 13, 12, 11, 10,  9,  8,  7,  6,  5,  4,  3,  2,  1)
)

DZ.candidates.29.11.21_B

#filter for selected candidates
glimpse(facet.hm.plot_3_Z)

facet.hm.plot_filter.custom.DZ <- facet.hm.plot_3_Z %>% filter(Gene.Name_AApos.collapsed %in% c(DZ.candidates.29.11.21_B$Gene.Name_AApos.collapsed)) %>% unique()
glimpse(facet.hm.plot_filter.custom.DZ)
table(duplicated(facet.hm.plot_filter.custom.DZ$Gene.Name_AApos.collapsed))

#merge order and adjust names upstream kinases arrange to custom order
facet.hm.plot_filter.custom.DZ <- merge(facet.hm.plot_filter.custom.DZ, DZ.candidates.29.11.21_B, by.x = "Gene.Name_AApos.collapsed", by.y = "Gene.Name_AApos.collapsed", all.x=T, all.y = T)
facet.hm.plot_filter.custom.DZ <- facet.hm.plot_filter.custom.DZ  %>% arrange(match(Group.DZ, c("MAPK signaling", "AKT signaling", "mTOR signaling")), Order.DZ.rev) %>% unique()


glimpse(facet.hm.plot_filter.custom.DZ)

## reshape for plot
melted.facet.hm.plot_filter.custom.DZ <- reshape2::melt(facet.hm.plot_filter.custom.DZ  %>% select(Gene.Name_AApos.collapsed, Group.DZ, "GFP.1","GFP.2","GFP.3","GFP.4",  "FL.1","FL.2","FL.3","FL.4",  "dE18.1","dE18.2","dE18.3","dE18.4",  "FL-Bicc1.1","FL-Bicc1.2","FL-Bicc1.3","FL-Bicc1.4",   "dE18-Bicc1.1","dE18-Bicc1.2","dE18-Bicc1.3"))
glimpse(melted.facet.hm.plot_filter.custom.DZ)

### facet plot selected candidates DZ
ggplot()+
  geom_tile(data=melted.facet.hm.plot_filter.custom.DZ, aes(y = fct_inorder(Gene.Name_AApos.collapsed), x = variable, fill = value))+ #, color = "white"
  scale_x_discrete(limits=c("GFP.1","GFP.2","GFP.3","GFP.4",  "FL.1","FL.2","FL.3","FL.4",  "dE18.1","dE18.2","dE18.3","dE18.4",  "FL-Bicc1.1","FL-Bicc1.2","FL-Bicc1.3","FL-Bicc1.4",   "dE18-Bicc1.1","dE18-Bicc1.2","dE18-Bicc1.3"))+

  scale_fill_scico(palette = 'vikO',name="row z-score \n log2(int.)", na.value = "grey50")+ #
 
  xlab(NULL) + 
  ylab(NULL)+
  theme(strip.background=element_rect(colour="transparent", fill="transparent"))+ # facet box color
  theme(plot.title = element_text(hjust = 0.5, face="bold"))+  #center title
  theme(legend.position="bottom", legend.justification = "center")+
  theme(axis.text.x = element_text(size = 10)) +
  theme(axis.text.y = element_text(size = 12))+
  facet_grid(fct_inorder(Group.DZ) ~ ., scales = "free", space = "free")+
  theme(strip.text.y = element_text(angle=0, face="bold", colour="black", size=10, lineheight = 5))+ #theme(strip.text.x = element_text(size =12))+ # facet box text size
  theme(axis.text.x  = element_text(angle=90, hjust = 1.0, vjust = 0.5, face="bold"))+
  theme(axis.text.y = element_text(size = 7))+
  geom_vline(xintercept = c(4.5, 8.5, 12.5, 16.5) , size=1.5, linetype = "solid", color="white")

#ggsave("DZ.cells.IMAC.global.phospho.selected.candidates.facet.plot.pdf", useDingbats=FALSE,  width = 14, height =14, units = "cm") #, width = 15, height = 30, units = "cm"


##################################################################################################################################
# > sessionInfo()
# R version 4.1.2 (2021-11-01)
# Platform: aarch64-apple-darwin20 (64-bit)
# Running under: macOS Monterey 12.2.1
# 
# Matrix products: default
# LAPACK: /Library/Frameworks/R.framework/Versions/4.1-arm64/Resources/lib/libRlapack.dylib
# 
# locale:
#   [1] de_DE.UTF-8/de_DE.UTF-8/de_DE.UTF-8/C/de_DE.UTF-8/de_DE.UTF-8
# 
# attached base packages:
#   [1] grid      stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] paletteer_1.4.0       gtools_3.9.2          msigdbr_7.4.1         cmapR_1.6.0           viridis_0.6.2        
# [6] viridisLite_0.4.0     circlize_0.4.14       ComplexHeatmap_2.10.0 pbapply_1.5-0         data.table_1.14.2    
# [11] forcats_0.5.1         stringr_1.4.0         dplyr_1.0.8           purrr_0.3.4           readr_2.1.2          
# [16] tidyr_1.2.0           tibble_3.1.6          ggplot2_3.3.5         tidyverse_1.3.1      
# 
# loaded via a namespace (and not attached):
#   [1] bitops_1.0-7                matrixStats_0.61.0          fs_1.5.2                    lubridate_1.8.0            
# [5] doParallel_1.0.17           RColorBrewer_1.1-2          httr_1.4.2                  GenomeInfoDb_1.30.1        
# [9] tools_4.1.2                 backports_1.4.1             utf8_1.2.2                  R6_2.5.1                   
# [13] DBI_1.1.2                   BiocGenerics_0.40.0         colorspace_2.0-3            GetoptLong_1.0.5           
# [17] withr_2.5.0                 tidyselect_1.1.2            gridExtra_2.3               compiler_4.1.2             
# [21] cli_3.2.0                   rvest_1.0.2                 Biobase_2.54.0              xml2_1.3.3                 
# [25] DelayedArray_0.20.0         prismatic_1.1.0             labeling_0.4.2              scales_1.1.1               
# [29] digest_0.6.29               XVector_0.34.0              dichromat_2.0-0             pkgconfig_2.0.3            
# [33] MatrixGenerics_1.6.0        maps_3.4.0                  dbplyr_2.1.1                rlang_1.0.2                
# [37] GlobalOptions_0.1.2         readxl_1.3.1                pals_1.7                    flowCore_2.6.0             
# [41] rstudioapi_0.13             farver_2.1.0                shape_1.4.6                 generics_0.1.2             
# [45] jsonlite_1.8.0              RCurl_1.98-1.6              magrittr_2.0.2              GenomeInfoDbData_1.2.7     
# [49] RProtoBufLib_2.6.0          Matrix_1.4-0                Rcpp_1.0.8.2                munsell_0.5.0              
# [53] S4Vectors_0.32.3            fansi_1.0.2                 babelgene_21.4              lifecycle_1.0.1            
# [57] stringi_1.7.6               SummarizedExperiment_1.24.0 zlibbioc_1.40.0             plyr_1.8.6                 
# [61] parallel_4.1.2              crayon_1.5.0                lattice_0.20-45             haven_2.4.3                
# [65] mapproj_1.2.8               hms_1.1.1                   pillar_1.7.0                GenomicRanges_1.46.1       
# [69] rjson_0.2.21                reshape2_1.4.4              codetools_0.2-18            stats4_4.1.2               
# [73] reprex_2.0.1                glue_1.6.2                  RcppParallel_5.1.5          BiocManager_1.30.16        
# [77] modelr_0.1.8                png_0.1-7                   vctrs_0.3.8                 tzdb_0.2.0                 
# [81] foreach_1.5.2               cellranger_1.1.0            gtable_0.3.0                rematch2_2.1.2             
# [85] clue_0.3-60                 assertthat_0.2.1            broom_0.7.12                iterators_1.0.14           
# [89] cytolib_2.6.2               IRanges_2.28.0              cluster_2.1.2               ellipsis_0.3.2 





























