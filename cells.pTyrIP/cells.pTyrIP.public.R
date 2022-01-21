#_______________________________________________________________________________________________________________________
# 19.01.2022
# 
# Project: proteomics.truncated.FGFR2.is.oncogene.cancer.Zingg.et.al 
# 
# Script name: cells.pTyrIP.public
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
### load data from MQ Phospho (STY)Sites.txt human search   ##################################################################################################


### sample overview labels FR
#1	GFP	1	A
#2	GFP	2	F
#3	GFP	3	K
#4	FL	1	B
#5	FL	2	G
#6	FL	3	L
#7	dE18	1	C
#8	dE18	2	H
#9	dE18	3	M
#10	FL.Bicc1	1	D
#11	FL.Bicc1	2	I
#12	FL.Bicc1	3	N
#13	dE18.Bicc1	1	E
#14	dE18.Bicc1	2	J
#15	dE18.Bicc1	3	O


### sample overview labels RdH
Lables.RdH <- tibble(label.RdH =  c(1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15),
                     number.FR = c(1,  4,  7, 10, 13,  2,  5,  8, 11, 14,  3,  6,  9, 12, 15),
                     groupFR = c("GFP","FL","dE18","FL.Bicc1","dE18.Bicc1","GFP","FL","dE18","FL.Bicc1","dE18.Bicc1","GFP","FL","dE18","FL.Bicc1","dE18.Bicc1"),
                     label.RdH.extended = c("pTyr_1","pTyr_2","pTyr_3","pTyr_4","pTyr_5","pTyr_6","pTyr_7","pTyr_8","pTyr_9","pTyr_10","pTyr_11","pTyr_12","pTyr_13","pTyr_14","pTyr_15")
                     )
Lables.RdH <- Lables.RdH %>% arrange(number.FR)
Lables.RdH 


#####################################################################################################################################################################
### load data cells IMAC DDA from MaxQuant Phospho(STY)Sites.txt

#DZ.pTyr.IP.MBR.ON <- fread("/Users/frankrolfs/NLPostDR/VUmc OPL/Phosphoproteomics with Daniel Zingg NKI/cells harvest 1 IP pTyr/txt-pTyr-QE3_210621_OPL1025_FR_pTyrIP_MEC/Phospho (STY)Sites.txt",integer64 = "numeric") #

DZ.pTyr.IP.MBR.ON <- fread("cells.pTyrIP_Phospho(STY)Sites.txt",integer64 = "numeric") #

glimpse(DZ.pTyr.IP.MBR.ON)
nrow(DZ.pTyr.IP.MBR.ON) #8477
ncol(DZ.pTyr.IP.MBR.ON) #246
colnames(DZ.pTyr.IP.MBR.ON)
length(unique(DZ.pTyr.IP.MBR.ON$id)) #8477



### select columns of interest
DZ.pTyr.IP.MBR.ON.2 <- DZ.pTyr.IP.MBR.ON %>% select(
  
  "id",                          "Proteins",                    "Positions within proteins",   "Leading proteins",            "Protein",                    
  "Protein names",               "Gene names",                  "Fasta headers",               "Localization prob",           "Number of Phospho (STY)",    
  "Amino acid",                  "Sequence window",             "Phospho (STY) Probabilities", "Position in peptide",         "Reverse",                    
  "Potential contaminant",       "Positions",                   "Position",                    "Peptide IDs",                 "Mod. peptide IDs",      
  
  matches("___[0-9]"), -Intensity___1, -Intensity___2, -Intensity___3
  
)

glimpse(DZ.pTyr.IP.MBR.ON.2)


### remove spaces from column names
column.names <- colnames(DZ.pTyr.IP.MBR.ON.2)
column.names
column.names.changed <- unname(sapply(column.names , function(x) str_replace_all(string=x, pattern=" ", replacement=".")))
column.names.changed
colnames(DZ.pTyr.IP.MBR.ON.2) <- column.names.changed

glimpse(DZ.pTyr.IP.MBR.ON.2)
colnames(DZ.pTyr.IP.MBR.ON.2)


#change column classes
DZ.pTyr.IP.MBR.ON.2 <- mutate_at(DZ.pTyr.IP.MBR.ON.2,  21:ncol(DZ.pTyr.IP.MBR.ON.2), list(as.numeric) )
glimpse(DZ.pTyr.IP.MBR.ON.2)


### remove reverse hits
table(DZ.pTyr.IP.MBR.ON.2$Reverse) #83+
DZ.pTyr.IP.MBR.ON.3 <- DZ.pTyr.IP.MBR.ON.2[DZ.pTyr.IP.MBR.ON.2$Reverse =="",]
nrow(DZ.pTyr.IP.MBR.ON.3) #8394


###remove contaminants
table(DZ.pTyr.IP.MBR.ON.3$Potential.contaminant) #41+
DZ.pTyr.IP.MBR.ON.3 <- DZ.pTyr.IP.MBR.ON.3[DZ.pTyr.IP.MBR.ON.3$Potential.contaminant == "",]
nrow(DZ.pTyr.IP.MBR.ON.3) #8353




### remove rows completely zero
glimpse(DZ.pTyr.IP.MBR.ON.3)
colnames(DZ.pTyr.IP.MBR.ON.3)
ncol(DZ.pTyr.IP.MBR.ON.3) #65
DZ.pTyr.IP.MBR.ON.3$zero.check <- pbapply(DZ.pTyr.IP.MBR.ON.3[, c(21:65)], 1, function(x) sum(x) ) #
table(DZ.pTyr.IP.MBR.ON.3$zero.check >0) #F1039 T7314

DZ.pTyr.IP.MBR.ON.3.complete.zero.rows <- filter(DZ.pTyr.IP.MBR.ON.3, zero.check <= 0)
glimpse(DZ.pTyr.IP.MBR.ON.3.complete.zero.rows) #

DZ.pTyr.IP.MBR.ON.3 <- filter(DZ.pTyr.IP.MBR.ON.3, zero.check > 0)
nrow(DZ.pTyr.IP.MBR.ON.3) #7314
glimpse(DZ.pTyr.IP.MBR.ON.3)
table(DZ.pTyr.IP.MBR.ON.3$zero.check > 0)


length(unique(DZ.pTyr.IP.MBR.ON.3$id))#7314 = total number PS
table(DZ.pTyr.IP.MBR.ON.3$Amino.acid) #S960 T484 Y5870
length(unique(DZ.pTyr.IP.MBR.ON.3$Sequence.window)) #7311

### remove column zero.check
DZ.pTyr.IP.MBR.ON.3 <- DZ.pTyr.IP.MBR.ON.3 %>% select(-zero.check)
glimpse(DZ.pTyr.IP.MBR.ON.3)


#create variable with short sample names with order as in Phospho(STY)Sites.txt
samples.in.experiment <- colnames(DZ.pTyr.IP.MBR.ON.3)
samples.in.experiment <- samples.in.experiment[21:ncol(DZ.pTyr.IP.MBR.ON.3)]
samples.in.experiment <- str_remove(string=samples.in.experiment, pattern="Intensity.")
samples.in.experiment <- str_remove(string=samples.in.experiment, pattern="___[0-9]")
samples.in.experiment <- unique(samples.in.experiment)
samples.in.experiment

### loop over data to count number of PS per sample ##############################################################################
colnames(DZ.pTyr.IP.MBR.ON.3)
ncol(DZ.pTyr.IP.MBR.ON.3) #65

result.PS.count.DZ.pTyr.IP.MBR.ON <- c()
sample.name.count <- 1

for(i.1 in seq(from=21, to=65, by=3)){ 
  print(i.1)
 
  i.2 <- i.1+1
  i.3 <- i.2+1
  
  sample.name <- samples.in.experiment[sample.name.count]
  print(sample.name)
  col.names <- colnames(DZ.pTyr.IP.MBR.ON.3)
  col.names.select <- col.names[c(i.1, i.2, i.3)]
  
  temp <- select(DZ.pTyr.IP.MBR.ON.3,Amino.acid, col.names.select)
  temp$sum_1_2_3 <- pbapply(temp[, 2:4], 1, function(x) sum(x))
  
  numberPS.A <- filter(temp, sum_1_2_3 > 0)
  numberPS <- nrow(numberPS.A); numberPS
  
  nMVs <- filter(temp, sum_1_2_3 <= 0)
  nMVs <- nrow(nMVs); nMVs
  
  percentageSTY <- table(numberPS.A$Amino.acid)
  
  temp.result <- tibble(sample=sample.name, numberPS=numberPS, nMVs=nMVs, nS=percentageSTY[1], nT=percentageSTY[2], nY=percentageSTY[3])
  
  result.PS.count.DZ.pTyr.IP.MBR.ON  <- bind_rows(result.PS.count.DZ.pTyr.IP.MBR.ON, temp.result)
  
  sample.name.count <- sample.name.count +1
  
}
result.PS.count.DZ.pTyr.IP.MBR.ON

glimpse(result.PS.count.DZ.pTyr.IP.MBR.ON)

mean(result.PS.count.DZ.pTyr.IP.MBR.ON$numberPS) #4204.333

##note: sample 8 suffered from a bad desalting tip 

Lables.RdH <- Lables.RdH %>% arrange(number.FR)
Lables.RdH 

#define order samples for e.g. plot
DZ.pTyr.IP.MBR.ON.samples.order <- c(
  "pTyr_1", "pTyr_6", "pTyr_11", #GFP
  
  
  "pTyr_2", "pTyr_7", "pTyr_12",  #FL
  
  "pTyr_3", "pTyr_8", "pTyr_13",  #dE18
  
  "pTyr_4", "pTyr_9", "pTyr_14", #FL-Bicc1
  
  "pTyr_5", "pTyr_10", "pTyr_15" #dE18-Bicc1
  
  
)
DZ.pTyr.IP.MBR.ON.samples.order

# calculate percentage MV PS all classes level
nrow(DZ.pTyr.IP.MBR.ON.3)  #7314 = total number PS

result.PS.count.DZ.pTyr.IP.MBR.ON$percentageMV <- pbsapply(result.PS.count.DZ.pTyr.IP.MBR.ON$nMVs, function(x) (100*x)/7314) #
result.PS.count.DZ.pTyr.IP.MBR.ON


# calculate percentage pS/pT/pY PS all classes level
for(i in 1:nrow(result.PS.count.DZ.pTyr.IP.MBR.ON)){
  print(i)
  
  total.PS.sample <- result.PS.count.DZ.pTyr.IP.MBR.ON$numberPS[i]
  n.pS <- result.PS.count.DZ.pTyr.IP.MBR.ON$nS[i]
  n.pT <- result.PS.count.DZ.pTyr.IP.MBR.ON$nT[i]
  n.pY <- result.PS.count.DZ.pTyr.IP.MBR.ON$nY[i]
  
  result.PS.count.DZ.pTyr.IP.MBR.ON$percentage.pS[i] <- (100*n.pS)/total.PS.sample
  result.PS.count.DZ.pTyr.IP.MBR.ON$percentage.pT[i] <- (100*n.pT)/total.PS.sample
  result.PS.count.DZ.pTyr.IP.MBR.ON$percentage.pY[i] <- (100*n.pY)/total.PS.sample
}
result.PS.count.DZ.pTyr.IP.MBR.ON
glimpse(result.PS.count.DZ.pTyr.IP.MBR.ON)

#reorder and add alternative saple name
result.PS.count.DZ.pTyr.IP.MBR.ON <- result.PS.count.DZ.pTyr.IP.MBR.ON[match(DZ.pTyr.IP.MBR.ON.samples.order, result.PS.count.DZ.pTyr.IP.MBR.ON$sample),]
result.PS.count.DZ.pTyr.IP.MBR.ON

result.PS.count.DZ.pTyr.IP.MBR.ON$sample2 <- c("GFP.1",
                                               "GFP.2",
                                               "GFP.3",
                                               "FL.1",
                                               "FL.2",
                                               "FL.3",
                                               "dE18.1",
                                               "dE18.2",
                                               "dE18.3",
                                               "FL.Bicc1.1",
                                               "FL.Bicc1.2",
                                               "FL.Bicc1.3",
                                               "dE18.Bicc1.1",
                                               "dE18.Bicc1.2",
                                               "dE18.Bicc1.3"
)

result.PS.count.DZ.pTyr.IP.MBR.ON

# calcuate means
mean.PS.allclassPS.per.sample <- mean(result.PS.count.DZ.pTyr.IP.MBR.ON$numberPS); mean.PS.allclassPS.per.sample #4204.333

mean.PS.allclassPS.per.sample.pS <- mean(result.PS.count.DZ.pTyr.IP.MBR.ON$percentage.pS);round(mean.PS.allclassPS.per.sample.pS, 1) #%pS 8
mean.PS.allclassPS.per.sample.pT <- mean(result.PS.count.DZ.pTyr.IP.MBR.ON$percentage.pT);round(mean.PS.allclassPS.per.sample.pT, 1) #%pT 3.9
mean.PS.allclassPS.per.sample.pY <- mean(result.PS.count.DZ.pTyr.IP.MBR.ON$percentage.pY);round(mean.PS.allclassPS.per.sample.pY, 1) #%pY 88.1

#mean count all class PS  GFP
result.PS.count.DZ.pTyr.IP.MBR.ON %>% filter(sample %in% c( "pTyr_1", "pTyr_6", "pTyr_11")) %>% select(numberPS) %>% pull() %>% mean() #2697.333
#mean count all class PS  FL
result.PS.count.DZ.pTyr.IP.MBR.ON %>% filter(sample %in% c("pTyr_2", "pTyr_7", "pTyr_12")) %>% select(numberPS) %>% pull() %>% mean() #4704.667
#mean count all class PS  dE18
result.PS.count.DZ.pTyr.IP.MBR.ON %>% filter(sample %in% c("pTyr_3", "pTyr_8", "pTyr_13")) %>% select(numberPS) %>% pull() %>% mean() #4701.667
#mean count all class PS  FL-Bicc1
result.PS.count.DZ.pTyr.IP.MBR.ON %>% filter(sample %in% c("pTyr_4", "pTyr_9", "pTyr_14")) %>% select(numberPS) %>% pull() %>% mean() #3891
#mean count all class PS  dE18-Bicc1
result.PS.count.DZ.pTyr.IP.MBR.ON %>% filter(sample %in% c("pTyr_5", "pTyr_10", "pTyr_15")) %>% select(numberPS) %>% pull() %>% mean() #5027




### loop over data to count number of PS per sample class 1 PS ###################################################################################################################################

glimpse(DZ.pTyr.IP.MBR.ON.3)

# Class1 PS filter for PS with localization probability > 0.75
table(DZ.pTyr.IP.MBR.ON.3$Localization.prob >= 0.75) 
DZ.pTyr.IP.MBR.ON.3.class1.COUNT <- DZ.pTyr.IP.MBR.ON.3[which(DZ.pTyr.IP.MBR.ON.3$Localization.prob >= 0.75),]
nrow(DZ.pTyr.IP.MBR.ON.3.class1.COUNT) #6069
table(DZ.pTyr.IP.MBR.ON.3.class1.COUNT$Localization.prob >= 0.75) 

### calculate number of class 1 phosphosite per sample
glimpse(DZ.pTyr.IP.MBR.ON.3.class1.COUNT)
colnames(DZ.pTyr.IP.MBR.ON.3.class1.COUNT)
ncol(DZ.pTyr.IP.MBR.ON.3.class1.COUNT) #65

#loop over data to count number of class 1 PS per sample
result.class1.PS.count.DZ.pTyr.IP.MBR.ON <- c()
sample.name.count <- 1

for(i.1 in seq(from=21, to=65, by=3)){ #53
  print(i.1)
 
  i.2 <- i.1+1
  i.3 <- i.2+1
  
  sample.name <- samples.in.experiment[sample.name.count]
  col.names <- colnames(DZ.pTyr.IP.MBR.ON.3.class1.COUNT)
  col.names.select <- col.names[c(i.1, i.2, i.3)]
  
  temp <- select(DZ.pTyr.IP.MBR.ON.3.class1.COUNT, Amino.acid, col.names.select)
  temp$sum_1_2_3 <- pbapply(temp[, 2:4], 1, function(x) sum(x))
  
  numberPS.A <- filter(temp, sum_1_2_3 > 0)
  numberPS <- nrow(numberPS.A)
  
  nMVs <- filter(temp, sum_1_2_3 <= 0)
  nMVs <- nrow(nMVs); nMVs
  
  percentageSTY <- table(numberPS.A$Amino.acid)
  
  temp.result <- tibble(class1.sample=sample.name, class1.numberPS=numberPS, class1.nMVs=nMVs, class1.nS=percentageSTY[1], class1.nT=percentageSTY[2], class1.nY=percentageSTY[3])
  
  result.class1.PS.count.DZ.pTyr.IP.MBR.ON  <- bind_rows(result.class1.PS.count.DZ.pTyr.IP.MBR.ON , temp.result)
  
  sample.name.count <- sample.name.count +1
}

result.class1.PS.count.DZ.pTyr.IP.MBR.ON 


# calculate percentage MV PS class1
result.class1.PS.count.DZ.pTyr.IP.MBR.ON$percentageMV <- pbsapply(result.class1.PS.count.DZ.pTyr.IP.MBR.ON $class1.nMVs, function(x) (100*x)/6069) 
result.class1.PS.count.DZ.pTyr.IP.MBR.ON 


# calculate percentage pS/pT/pY PS class1
for(i in 1:nrow(result.class1.PS.count.DZ.pTyr.IP.MBR.ON )){
  print(i)
  
  total.PS.sample <- result.class1.PS.count.DZ.pTyr.IP.MBR.ON $class1.numberPS[i]
  n.pS <- result.class1.PS.count.DZ.pTyr.IP.MBR.ON $class1.nS[i]
  n.pT <- result.class1.PS.count.DZ.pTyr.IP.MBR.ON $class1.nT[i]
  n.pY <- result.class1.PS.count.DZ.pTyr.IP.MBR.ON $class1.nY[i]
  
  result.class1.PS.count.DZ.pTyr.IP.MBR.ON $class1.percentage.pS[i] <- (100*n.pS)/total.PS.sample
  result.class1.PS.count.DZ.pTyr.IP.MBR.ON $class1.percentage.pT[i] <- (100*n.pT)/total.PS.sample
  result.class1.PS.count.DZ.pTyr.IP.MBR.ON $class1.percentage.pY[i] <- (100*n.pY)/total.PS.sample
}
glimpse(result.class1.PS.count.DZ.pTyr.IP.MBR.ON )

#reorder and add alternative sample name
result.class1.PS.count.DZ.pTyr.IP.MBR.ON <- result.class1.PS.count.DZ.pTyr.IP.MBR.ON[match(DZ.pTyr.IP.MBR.ON.samples.order, result.class1.PS.count.DZ.pTyr.IP.MBR.ON$class1.sample),]
result.class1.PS.count.DZ.pTyr.IP.MBR.ON


result.class1.PS.count.DZ.pTyr.IP.MBR.ON$sample2 <- c("GFP.1",
                                                      "GFP.2",
                                                      "GFP.3",
                                                      "FL.1",
                                                      "FL.2",
                                                      "FL.3",
                                                      "dE18.1",
                                                      "dE18.2",
                                                      "dE18.3",
                                                      "FL.Bicc1.1",
                                                      "FL.Bicc1.2",
                                                      "FL.Bicc1.3",
                                                      "dE18.Bicc1.1",
                                                      "dE18.Bicc1.2",
                                                      "dE18.Bicc1.3")


result.class1.PS.count.DZ.pTyr.IP.MBR.ON

# calcuate means
mean.result.class1.PS.count.DZ.pTyr.IP.MBR.ON  <- mean(result.class1.PS.count.DZ.pTyr.IP.MBR.ON $class1.numberPS); mean.result.class1.PS.count.DZ.pTyr.IP.MBR.ON   #3868.4

#mean count all class PS  GFP
result.class1.PS.count.DZ.pTyr.IP.MBR.ON %>% filter(class1.sample %in% c("pTyr_1", "pTyr_3", "pTyr_11")) %>% select(class1.numberPS) %>% pull() %>% mean() #3258.33
#mean count all class PS  FL
result.class1.PS.count.DZ.pTyr.IP.MBR.ON %>% filter(class1.sample %in% c("pTyr_2", "pTyr_7", "pTyr_12")) %>% select(class1.numberPS) %>% pull() %>% mean() #4305.667
#mean count all class PS  dE18
result.class1.PS.count.DZ.pTyr.IP.MBR.ON %>% filter(class1.sample %in% c("pTyr_3", "pTyr_8", "pTyr_13")) %>% select(class1.numberPS) %>% pull() %>% mean() #4296.33 with sample dE18.2
result.class1.PS.count.DZ.pTyr.IP.MBR.ON %>% filter(class1.sample %in% c("pTyr_3",           "pTyr_13")) %>% select(class1.numberPS) %>% pull() %>% mean() #4296.33 without sample (8) dE18.2
#mean count all class PS  FL-Bicc1
result.class1.PS.count.DZ.pTyr.IP.MBR.ON %>% filter(class1.sample %in% c("pTyr_4", "pTyr_9", "pTyr_14")) %>% select(class1.numberPS) %>% pull() %>% mean() #3597
#mean count all class PS  dE18-Bicc1
result.class1.PS.count.DZ.pTyr.IP.MBR.ON %>% filter(class1.sample %in% c("pTyr_5", "pTyr_10", "pTyr_15")) %>% select(class1.numberPS) %>% pull() %>% mean() #4646



### plot number of PS CLASS1 AND number of PS (all classes)
ggplot()+
  geom_col(data = result.PS.count.DZ.pTyr.IP.MBR.ON, aes(x=sample2, y=numberPS))+
  geom_col(data= result.class1.PS.count.DZ.pTyr.IP.MBR.ON, aes(x=sample2, y=class1.numberPS), color="black", fill="yellow")+
  scale_x_discrete(limits= c("GFP.1",
                             "GFP.2",
                             "GFP.3",
                             "FL.1",
                             "FL.2",
                             "FL.3",
                             "dE18.1",
                             "dE18.2",
                             "dE18.3",
                             "FL.Bicc1.1",
                             "FL.Bicc1.2",
                             "FL.Bicc1.3",
                             "dE18.Bicc1.1",
                             "dE18.Bicc1.2",
                             "dE18.Bicc1.3"))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1.0, vjust = 0.5))+
  ylab("number PS")+
  xlab(NULL)+
  scale_y_continuous(breaks=seq(from=0, to=5500, by=500), limits = c(0,5500))+
  ggtitle("# PS (all classes) per sample - grey \n # class 1PS per sample - yellow") +
  theme(plot.title = element_text(hjust = 0.5, face= "bold"))+  #center title
  theme(axis.text.x = element_text(face= "plain", colour="black", size=14))+
  theme(axis.text.y = element_text(face= "plain", colour="black", size=14))+
  geom_vline(xintercept = c(3.5, 6.5, 9.5, 12.5) , size=0.25, linetype = 2)


### plot number class 1 pY
ggplot()+
  geom_col(data= result.class1.PS.count.DZ.pTyr.IP.MBR.ON, aes(x=sample2, y=class1.nY), color="black", fill="yellow")+
  geom_col(data= result.class1.PS.count.DZ.pTyr.IP.MBR.ON, aes(x=sample2, y=class1.nS), color="black", fill="blue")+
  geom_col(data= result.class1.PS.count.DZ.pTyr.IP.MBR.ON, aes(x=sample2, y=class1.nT), color="black", fill="black")+
  #scale_x_discrete(limits= c(DZ.pTyr.IP.MBR.ON.samples.order))+
  scale_x_discrete(limits= c("GFP.1",
                             "GFP.2",
                             "GFP.3",
                             "FL.1",
                             "FL.2",
                             "FL.3",
                             "dE18.1",
                             "dE18.2",
                             "dE18.3",
                             "FL.Bicc1.1",
                             "FL.Bicc1.2",
                             "FL.Bicc1.3",
                             "dE18.Bicc1.1",
                             "dE18.Bicc1.2",
                             "dE18.Bicc1.3"))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1.0, vjust = 0.5))+
  ylab("number class 1 pY")+
  xlab(NULL)+
  scale_y_continuous(breaks=seq(from=0, to=5500, by=500), limits = c(0,5500))+
  ggtitle(" # PS class 1 Tyr - yellow \n # PS class 1 Ser - blue  \n # PS class 1 Thr - black") +
  theme(plot.title = element_text(hjust = 0.5, face= "bold"))+  #center title
  theme(axis.text.x = element_text(face= "plain", colour="black", size=14))+
  theme(axis.text.y = element_text(face= "plain", colour="black", size=14))+
  geom_vline(xintercept = c(3.5, 6.5, 9.5, 12.5) , size=0.25, linetype = 2)


## means class1.percentage.pY
#mean count all class PS  GFP
result.class1.PS.count.DZ.pTyr.IP.MBR.ON %>% filter(class1.sample %in% c("pTyr_1", "pTyr_6", "pTyr_11")) %>% select(class1.percentage.pY) %>% pull() %>% mean() #93.04
#mean count all class PS  FL
result.class1.PS.count.DZ.pTyr.IP.MBR.ON %>% filter(class1.sample %in% c("pTyr_2", "pTyr_7", "pTyr_12")) %>% select(class1.percentage.pY) %>% pull() %>% mean() #93.14
#mean count all class PS  dE18
result.class1.PS.count.DZ.pTyr.IP.MBR.ON %>% filter(class1.sample %in% c("pTyr_3", "pTyr_8", "pTyr_13")) %>% select(class1.percentage.pY) %>% pull() %>% mean() #93.35
#mean count all class PS  FL-Bicc1
result.class1.PS.count.DZ.pTyr.IP.MBR.ON %>% filter(class1.sample %in% c("pTyr_4", "pTyr_9", "pTyr_14")) %>% select(class1.percentage.pY) %>% pull() %>% mean() #93.18
#mean count all class PS  dE18-Bicc1
result.class1.PS.count.DZ.pTyr.IP.MBR.ON %>% filter(class1.sample %in% c("pTyr_5", "pTyr_10", "pTyr_15")) %>% select(class1.percentage.pY) %>% pull() %>% mean() #94.41

## means class1 number pY
#mean count all class PS  GFP
result.class1.PS.count.DZ.pTyr.IP.MBR.ON %>% filter(class1.sample %in% c("pTyr_1", "pTyr_6", "pTyr_11")) %>% select(class1.nY) %>% pull() %>% mean() #2323
#mean count all class PS  FL
result.class1.PS.count.DZ.pTyr.IP.MBR.ON %>% filter(class1.sample %in% c("pTyr_2", "pTyr_7", "pTyr_12")) %>% select(class1.nY) %>% pull() %>% mean() #4010.333
#mean count all class PS  dE18
result.class1.PS.count.DZ.pTyr.IP.MBR.ON %>% filter(class1.sample %in% c("pTyr_3", "pTyr_8", "pTyr_13")) %>% select(class1.nY) %>% pull() %>% mean() #4009.33 with dE18.2
result.class1.PS.count.DZ.pTyr.IP.MBR.ON %>% filter(class1.sample %in% c("pTyr_3",           "pTyr_13")) %>% select(class1.nY) %>% pull() %>% mean() #4359 without dE18.2
#mean count all class PS  FL-Bicc1
result.class1.PS.count.DZ.pTyr.IP.MBR.ON %>% filter(class1.sample %in% c("pTyr_4", "pTyr_9", "pTyr_14")) %>% select(class1.nY) %>% pull() %>% mean() #3351.667
#mean count all class PS  dE18-Bicc1
result.class1.PS.count.DZ.pTyr.IP.MBR.ON %>% filter(class1.sample %in% c("pTyr_5", "pTyr_10", "pTyr_15")) %>% select(class1.nY) %>% pull() %>% mean() #4386.667




#### cells pTyrIP expand with PerseusPSexpand function (see above) ################################################################################################################
#### cells pTyrIP expand with PerseusPSexpand function (see above) ################################################################################################################

glimpse(DZ.pTyr.IP.MBR.ON.3)

DZ.pTyr.IP.MBR.ON_expanded <- PerseusPSexpand.17.02.21(DZ.pTyr.IP.MBR.ON.3)


glimpse(DZ.pTyr.IP.MBR.ON_expanded)



#remove rows that have all zero intensities
glimpse(DZ.pTyr.IP.MBR.ON_expanded)
colnames(DZ.pTyr.IP.MBR.ON_expanded)

DZ.pTyr.IP.MBR.ON_expanded$ZEROcheckINT <- pbapply(DZ.pTyr.IP.MBR.ON_expanded[, c(21:35)], MARGIN=1, function(x) sum(x)) #
table(DZ.pTyr.IP.MBR.ON_expanded$ZEROcheckINT > 0)
nrow(DZ.pTyr.IP.MBR.ON_expanded) #21942

DZ.pTyr.IP.MBR.ON_expanded <- filter(DZ.pTyr.IP.MBR.ON_expanded, ZEROcheckINT > 0)
nrow(DZ.pTyr.IP.MBR.ON_expanded) #8105


#remove column ZEROcheckINT
DZ.pTyr.IP.MBR.ON_expanded <- select(DZ.pTyr.IP.MBR.ON_expanded,  -ZEROcheckINT)
glimpse(DZ.pTyr.IP.MBR.ON_expanded)



# cells pTyrIP normalization & histogram standard median shift version #########################################################################################################################################
# cells pTyrIP normalization & histogram standard median shift version #########################################################################################################################################

glimpse(DZ.pTyr.IP.MBR.ON_expanded)
colnames(DZ.pTyr.IP.MBR.ON_expanded)

#log2 transform & replace -Inf with NA
colnames(DZ.pTyr.IP.MBR.ON_expanded)
Log2Trafo <- DZ.pTyr.IP.MBR.ON_expanded[, c(21:35)] #
Log2Trafo <- log2(Log2Trafo)
Log2Trafo[Log2Trafo  == -Inf] <- NA

glimpse(Log2Trafo)

DZ.pTyr.IP.MBR.ON_expanded_log2 <- DZ.pTyr.IP.MBR.ON_expanded
DZ.pTyr.IP.MBR.ON_expanded_log2[, c(21:35)] <- Log2Trafo #

glimpse(DZ.pTyr.IP.MBR.ON_expanded_log2, list.len=900)


## reorder and calculate median per column RawINT
colnames(DZ.pTyr.IP.MBR.ON_expanded_log2)


DZ.pTyr.IP.MBR.ON_expanded_log2 <- DZ.pTyr.IP.MBR.ON_expanded_log2 %>% select("id",                          "Proteins",                    "Positions.within.proteins",   "Leading.proteins",            "Protein",                    
                                                                              "Protein.names",               "Gene.names",                  "Fasta.headers",               "Localization.prob",           "Number.of.Phospho.(STY)",    
                                                                              "Amino.acid",                  "Sequence.window",             "Phospho.(STY).Probabilities", "Position.in.peptide",         "Reverse",                    
                                                                              "Potential.contaminant",       "Positions",                   "Position",                    "Peptide.IDs",                 "Mod..peptide.IDs",   
                                                                              "Intensity.pTyr_1", "Intensity.pTyr_6" , "Intensity.pTyr_11", 
                                                                              "Intensity.pTyr_2" , "Intensity.pTyr_7" , "Intensity.pTyr_12", 
                                                                              "Intensity.pTyr_3" , "Intensity.pTyr_8" , "Intensity.pTyr_13",
                                                                              "Intensity.pTyr_4" , "Intensity.pTyr_9" , "Intensity.pTyr_14",
                                                                              "Intensity.pTyr_5" , "Intensity.pTyr_10" , "Intensity.pTyr_15" ,
                                                                              
                                                                              "PS_Multiplicity"
)

glimpse(DZ.pTyr.IP.MBR.ON_expanded_log2)


ColMedianfindDF <- glimpse(DZ.pTyr.IP.MBR.ON_expanded_log2[, c(21:35)]) #
vectormediancolumnsINT <- c()
for(i in 1:ncol(ColMedianfindDF)){
  print(i)
  temp <- median(ColMedianfindDF[, i], na.rm = TRUE); temp
  vectormediancolumnsINT <- append(vectormediancolumnsINT, temp)
}
vectormediancolumnsINT

median(vectormediancolumnsINT)


# data frame to plot medians as lines
Mdf <- data.frame(variable=colnames(DZ.pTyr.IP.MBR.ON_expanded_log2[, c(21:35)]), median=vectormediancolumnsINT); Mdf
Mdf$variable2 <- c("GFP.1",
                   "GFP.2",
                   "GFP.3",
                   "FL.1",
                   "FL.2",
                   "FL.3",
                   "dE18.1",
                   "dE18.2",
                   "dE18.3",
                   "FL.Bicc1.1",
                   "FL.Bicc1.2",
                   "FL.Bicc1.3",
                   "dE18.Bicc1.1",
                   "dE18.Bicc1.2",
                   "dE18.Bicc1.3")

Mdf$group <- c(rep("GFP", 3), rep("FL", 3), rep("dE18", 3), rep("FL.Bicc1", 3), rep("dE18.Bicc1", 3))
Mdf

## Plot
DensityPlotData <- DZ.pTyr.IP.MBR.ON_expanded_log2[, c(21:35)]

glimpse(DensityPlotData)

DensityPlotDataMelted <- reshape2::melt(DensityPlotData)
glimpse(DensityPlotDataMelted)

#density plot
ggplot(DensityPlotDataMelted, aes(x=value, colour = variable) ) + 
  geom_density() +
  theme(legend.position="bottom") +
  geom_vline(data=Mdf, aes(xintercept=median, color=variable))






## Start define function to calculate median shift factors 
MedianShiftFactorFinder <- function(arg){
  tempfactor <- 23.2535/arg #log2(1e07) = 23.2535
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
temp_for_median_shift <-  DZ.pTyr.IP.MBR.ON_expanded_log2[, c(21:35)]
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

#plot density plot check
ggplot(DensityPlotDataCHECKMelted, aes(x=value, colour = variable) ) + 
  geom_density() +
  theme(legend.position="bottom") +
  geom_vline(data=MdfCHECK, aes(xintercept=median, color=variable))



### rename columns after normalization
newCOLNAMESnorm <- colnames(temp_for_median_shift); newCOLNAMESnorm

newCOLNAMESnorm <- paste0("Norm.",newCOLNAMESnorm); newCOLNAMESnorm

colnames(temp_for_median_shift) <- newCOLNAMESnorm
glimpse(temp_for_median_shift)


### add norm INT to original Dataframe 

DZ.pTyr.IP.MBR.ON_expanded_log2_norm <- DZ.pTyr.IP.MBR.ON_expanded_log2

DZ.pTyr.IP.MBR.ON_expanded_log2_norm  <- cbind(DZ.pTyr.IP.MBR.ON_expanded_log2_norm , temp_for_median_shift)
glimpse(DZ.pTyr.IP.MBR.ON_expanded_log2_norm, list.len=999)


## box plot after normalization
data.for.boxplot.after.norm.DZ.pTyr.IP.MBR.ON <- DZ.pTyr.IP.MBR.ON_expanded_log2_norm %>% select("Norm.Intensity.pTyr_1", "Norm.Intensity.pTyr_6" , "Norm.Intensity.pTyr_11", 
                                                                                                 "Norm.Intensity.pTyr_2" , "Norm.Intensity.pTyr_7" , "Norm.Intensity.pTyr_12", 
                                                                                                 "Norm.Intensity.pTyr_3" , "Norm.Intensity.pTyr_8" , "Norm.Intensity.pTyr_13",
                                                                                                 "Norm.Intensity.pTyr_4" , "Norm.Intensity.pTyr_9" , "Norm.Intensity.pTyr_14",
                                                                                                 "Norm.Intensity.pTyr_5" , "Norm.Intensity.pTyr_10" , "Norm.Intensity.pTyr_15")
glimpse(data.for.boxplot.after.norm.DZ.pTyr.IP.MBR.ON)
colnames(data.for.boxplot.after.norm.DZ.pTyr.IP.MBR.ON) <- c("GFP.1",
                                                             "GFP.2",
                                                             "GFP.3",
                                                             "FL.1",
                                                             "FL.2",
                                                             "FL.3",
                                                             "dE18.1",
                                                             "dE18.2",
                                                             "dE18.3",
                                                             "FL.Bicc1.1",
                                                             "FL.Bicc1.2",
                                                             "FL.Bicc1.3",
                                                             "dE18.Bicc1.1",
                                                             "dE18.Bicc1.2",
                                                             "dE18.Bicc1.3")
glimpse(data.for.boxplot.after.norm.DZ.pTyr.IP.MBR.ON)

melted.data.for.boxplot.after.norm.DZ.pTyr.IP.MBR.ON <- reshape2::melt(data.for.boxplot.after.norm.DZ.pTyr.IP.MBR.ON)
glimpse(melted.data.for.boxplot.after.norm.DZ.pTyr.IP.MBR.ON )
unique(melted.data.for.boxplot.after.norm.DZ.pTyr.IP.MBR.ON$variable)

ggplot(data = melted.data.for.boxplot.after.norm.DZ.pTyr.IP.MBR.ON)+
  geom_boxplot(mapping = aes(x= variable, y = value, fill=variable), outlier.size = 0.5, notch = TRUE, notchwidth = 0.1)+
  scale_x_discrete(limits= c("GFP.1",
                             "GFP.2",
                             "GFP.3",
                             "FL.1",
                             "FL.2",
                             "FL.3",
                             "dE18.1",
                             "dE18.2",
                             "dE18.3",
                             "FL.Bicc1.1",
                             "FL.Bicc1.2",
                             "FL.Bicc1.3",
                             "dE18.Bicc1.1",
                             "dE18.Bicc1.2",
                             "dE18.Bicc1.3"))+
  scale_fill_manual(values = c("yellow", "yellow", "yellow",  "blue", "blue", "blue",  "red", "red", "red", "black","black","black", "grey40", "grey40", "grey40") )+ 
  theme(axis.text.x = element_text(angle = 90, hjust = 1.0, vjust = 0.5))+
  theme(legend.position="none", legend.justification = "center")+
  ggtitle("cells pTyrIP - normalized data - PS level") +
  theme(plot.title = element_text(hjust = 0.5, face= "bold"))+ 
  theme(axis.text.x = element_text(face= "plain", colour="black", size=14))+
  theme(axis.text.y = element_text(face= "plain", colour="black", size=14))+
  xlab(NULL)+ 
  ylab("log2(Int.)") 


## box plot WITHOUT normalization
glimpse(DZ.pTyr.IP.MBR.ON_expanded_log2_norm)

data.for.boxplot.after.norm.DZ.pTyr.IP.MBR.ON <- DZ.pTyr.IP.MBR.ON_expanded_log2_norm %>% select("Intensity.pTyr_1", "Intensity.pTyr_6" , "Intensity.pTyr_11", 
                                                                                                 "Intensity.pTyr_2" , "Intensity.pTyr_7" , "Intensity.pTyr_12", 
                                                                                                 "Intensity.pTyr_3" , "Intensity.pTyr_8" , "Intensity.pTyr_13",
                                                                                                 "Intensity.pTyr_4" , "Intensity.pTyr_9" , "Intensity.pTyr_14",
                                                                                                 "Intensity.pTyr_5" , "Intensity.pTyr_10" , "Intensity.pTyr_15")
glimpse(data.for.boxplot.after.norm.DZ.pTyr.IP.MBR.ON)
colnames(data.for.boxplot.after.norm.DZ.pTyr.IP.MBR.ON) <- c("GFP.1",
                                                             "GFP.2",
                                                             "GFP.3",
                                                             "FL.1",
                                                             "FL.2",
                                                             "FL.3",
                                                             "dE18.1",
                                                             "dE18.2",
                                                             "dE18.3",
                                                             "FL.Bicc1.1",
                                                             "FL.Bicc1.2",
                                                             "FL.Bicc1.3",
                                                             "dE18.Bicc1.1",
                                                             "dE18.Bicc1.2",
                                                             "dE18.Bicc1.3")
glimpse(data.for.boxplot.after.norm.DZ.pTyr.IP.MBR.ON)

melted.data.for.boxplot.after.norm.DZ.pTyr.IP.MBR.ON <- reshape2::melt(data.for.boxplot.after.norm.DZ.pTyr.IP.MBR.ON)
glimpse(melted.data.for.boxplot.after.norm.DZ.pTyr.IP.MBR.ON )
unique(melted.data.for.boxplot.after.norm.DZ.pTyr.IP.MBR.ON$variable)

ggplot(data = melted.data.for.boxplot.after.norm.DZ.pTyr.IP.MBR.ON)+
  geom_boxplot(mapping = aes(x= variable, y = value, fill=variable), outlier.size = 0.5, notch = TRUE, notchwidth = 0.1)+
  scale_x_discrete(limits= c("GFP.1",
                             "GFP.2",
                             "GFP.3",
                             "FL.1",
                             "FL.2",
                             "FL.3",
                             "dE18.1",
                             "dE18.2",
                             "dE18.3",
                             "FL.Bicc1.1",
                             "FL.Bicc1.2",
                             "FL.Bicc1.3",
                             "dE18.Bicc1.1",
                             "dE18.Bicc1.2",
                             "dE18.Bicc1.3"))+
  scale_fill_manual(values = c("yellow", "yellow", "yellow",  "blue", "blue", "blue",  "red", "red", "red", "black","black","black", "grey40", "grey40", "grey40") )+ 
  theme(axis.text.x = element_text(angle = 90, hjust = 1.0, vjust = 0.5))+
  theme(legend.position="none", legend.justification = "center")+
  ggtitle("cells pTyrIP - NOT normalized data - PS level") +
  theme(plot.title = element_text(hjust = 0.5, face= "bold"))+  #center title
  theme(axis.text.x = element_text(face= "plain", colour="black", size=14))+
  theme(axis.text.y = element_text(face= "plain", colour="black", size=14))+
  xlab(NULL)+ 
  ylab("log2(Int.)") 


#####################################################################################################################################################################
#### analysis via class 1 PS ##################################################################################################################################
glimpse(DZ.pTyr.IP.MBR.ON_expanded_log2_norm) #

DZ.pTyr.IP.MBR.ON_expanded_log2_norm.class1.PS <- DZ.pTyr.IP.MBR.ON_expanded_log2_norm %>% filter(Localization.prob>= 0.75)
glimpse(DZ.pTyr.IP.MBR.ON_expanded_log2_norm.class1.PS) # 6798


# add additional column to better calculate number PS 
# take Proteins and add Psite without multiplicity info ###
for(i in 1:nrow(DZ.pTyr.IP.MBR.ON_expanded_log2_norm.class1.PS)){
  print(i)
  temp_prot_name <- DZ.pTyr.IP.MBR.ON_expanded_log2_norm.class1.PS$Proteins[i]
  tempPSite <- DZ.pTyr.IP.MBR.ON_expanded_log2_norm.class1.PS$Positions[i]
  tempAA <- DZ.pTyr.IP.MBR.ON_expanded_log2_norm.class1.PS$Amino.acid[i]
  DZ.pTyr.IP.MBR.ON_expanded_log2_norm.class1.PS$Prot.Name_AApos[i] <- paste0(temp_prot_name, "_", tempAA, tempPSite)
}
length(unique(DZ.pTyr.IP.MBR.ON_expanded_log2_norm.class1.PS$Prot.Name_AApos)) #6069

glimpse(DZ.pTyr.IP.MBR.ON_expanded_log2_norm.class1.PS)

#  prepare new column
# take genename and add Psite ###
for(i in 1:nrow(DZ.pTyr.IP.MBR.ON_expanded_log2_norm.class1.PS)){
  print(i)
  temp_gene_name <- DZ.pTyr.IP.MBR.ON_expanded_log2_norm.class1.PS$Gene.names[i]
  tempPSite <- DZ.pTyr.IP.MBR.ON_expanded_log2_norm.class1.PS$Positions[i]
  tempAA <- DZ.pTyr.IP.MBR.ON_expanded_log2_norm.class1.PS$Amino.acid[i]
  tempPSmultiplicity <- DZ.pTyr.IP.MBR.ON_expanded_log2_norm.class1.PS$PS_Multiplicity[i]
  
  DZ.pTyr.IP.MBR.ON_expanded_log2_norm.class1.PS$Gene.Name_AApos[i] <- paste0(temp_gene_name, "_", tempAA, tempPSite, "x", tempPSmultiplicity)
}
glimpse(DZ.pTyr.IP.MBR.ON_expanded_log2_norm.class1.PS)


#add identifier for plotting reason PTMSEA later
glimpse(DZ.pTyr.IP.MBR.ON_expanded_log2_norm.class1.PS)

DZ.pTyr.IP.MBR.ON_expanded_log2_norm.class1.PS$ID.FR.all.C1.PS <- 1:nrow(DZ.pTyr.IP.MBR.ON_expanded_log2_norm.class1.PS)

#save data
#save(DZ.pTyr.IP.MBR.ON_expanded_log2_norm.class1.PS, file = "cells.pTyr.IP.MBR.ON_expanded_log2_norm.class1.PS.Rdata")
#load("cells.pTyr.IP.MBR.ON_expanded_log2_norm.class1.PS.Rdata")
#glimpse(DZ.pTyr.IP.MBR.ON_expanded_log2_norm.class1.PS)




############################################################################################################################################################
############################################################################################################################################################

### cells pTyrIP
### heatmap for all Fgfr2 sites in dataset for all groups with collapsing _1/_2/_3 via sum 
### PS numbers also for NP_963895.2

glimpse(DZ.pTyr.IP.MBR.ON_expanded_log2_norm.class1.PS) 

HM.plot.data.pTyr <- DZ.pTyr.IP.MBR.ON_expanded_log2_norm.class1.PS%>% filter(stringr::str_detect(Proteins, "P21803(?!-)"))

glimpse(HM.plot.data.pTyr)

## pick columns of interest and rename
HM.plotdata.pTyr <- HM.plot.data.pTyr %>% select("Gene.Name_AApos", 
                                                 "Norm.Intensity.pTyr_1", "Norm.Intensity.pTyr_6" , "Norm.Intensity.pTyr_11", 
                                                 "Norm.Intensity.pTyr_2" , "Norm.Intensity.pTyr_7" , "Norm.Intensity.pTyr_12", 
                                                 "Norm.Intensity.pTyr_3" , "Norm.Intensity.pTyr_8" , "Norm.Intensity.pTyr_13",
                                                 "Norm.Intensity.pTyr_4" , "Norm.Intensity.pTyr_9" , "Norm.Intensity.pTyr_14",
                                                 "Norm.Intensity.pTyr_5" , "Norm.Intensity.pTyr_10" , "Norm.Intensity.pTyr_15" ,
                                                 "Proteins",
                                                 "Positions.within.proteins",
                                                 "Amino.acid",
                                                 "Gene.names" ,
                                                 "PS_Multiplicity")
glimpse(HM.plotdata.pTyr)

colnames(HM.plotdata.pTyr) <- c("Gene.Name_AApos", 
                                "GFP.1","GFP.2","GFP.3",  
                                "FL.1","FL.2","FL.3",  
                                "dE18.1","dE18.2","dE18.3", 
                                "FL-Bicc1.1","FL-Bicc1.2","FL-Bicc1.3", 
                                "dE18-Bicc1.1","dE18-Bicc1.2","dE18-Bicc1.3",
                                "Proteins",
                                "Positions.within.proteins",
                                "Amino.acid",
                                "Gene.names" ,
                                "PS_Multiplicity")
glimpse(HM.plotdata.pTyr)


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
for(i in 1:nrow(HM.plotdata.pTyr)){
  print(i)
  HM.plotdata.pTyr$isoform.1.corrected.PS.position[i] <- position.in.mq.string(prot.of.interest ="P21803", #isoform 1
                                                                               Gene.names=HM.plotdata.pTyr$Gene.names[i],
                                                                               Amino.acid=HM.plotdata.pTyr$Amino.acid[i], 
                                                                               mq.position.string=HM.plotdata.pTyr$Positions.within.proteins[i],
                                                                               mq.string = HM.plotdata.pTyr$Proteins[i],
                                                                               PS_Multiplicity=HM.plotdata.pTyr$PS_Multiplicity[i])
}
HM.plotdata.pTyr$isoform.1.corrected.Gene.Name_AApos = paste0(HM.plotdata.pTyr$Gene.names, "_", HM.plotdata.pTyr$Amino.acid, HM.plotdata.pTyr$isoform.1.corrected.PS.position, "x", HM.plotdata.pTyr$PS_Multiplicity)
glimpse(HM.plotdata.pTyr)

# isoform 2 ="P21803-2"
for(i in 1:nrow(HM.plotdata.pTyr)){
  print(i)
  HM.plotdata.pTyr$isoform.2.corrected.PS.position[i] <- position.in.mq.string(prot.of.interest ="P21803-2", #isoform 2
                                                                               Gene.names=HM.plotdata.pTyr$Gene.names[i],
                                                                               Amino.acid=HM.plotdata.pTyr$Amino.acid[i], 
                                                                               mq.position.string=HM.plotdata.pTyr$Positions.within.proteins[i],
                                                                               mq.string = HM.plotdata.pTyr$Proteins[i],
                                                                               PS_Multiplicity=HM.plotdata.pTyr$PS_Multiplicity[i])
}
HM.plotdata.pTyr$isoform.2.corrected.Gene.Name_AApos = paste0(HM.plotdata.pTyr$Gene.names, "_", HM.plotdata.pTyr$Amino.acid, HM.plotdata.pTyr$isoform.2.corrected.PS.position, "x", HM.plotdata.pTyr$PS_Multiplicity)

glimpse(HM.plotdata.pTyr)

#prepare helper data frame for later order
HM.plotdata.pTyr.melted <- reshape2::melt(HM.plotdata.pTyr)
glimpse(HM.plotdata.pTyr.melted)
colnames(HM.plotdata.pTyr.melted) <- c("Gene.Name_AApos", "Proteins","Positions.within.proteins", "Amino.acid","Gene.names", "PS_Multiplicity", "isoform.1.corrected.PS.position", "isoform.1.corrected.Gene.Name_AApos","isoform.2.corrected.PS.position", "isoform.2.corrected.Gene.Name_AApos",  "variable","log2.int")
glimpse(HM.plotdata.pTyr.melted)


order.df.HM.plotdata.pTyr.melted <- unique(tibble(Gene.Name_AApos = HM.plotdata.pTyr.melted$Gene.Name_AApos,
                                                  Proteins = HM.plotdata.pTyr.melted$Proteins,
                                                  Positions.within.proteins = HM.plotdata.pTyr.melted$Positions.within.proteins,
                                                  Amino.acid = HM.plotdata.pTyr.melted$Amino.acid,
                                                  Gene.names = HM.plotdata.pTyr.melted$Gene.names,
                                                  PS_Multiplicity =HM.plotdata.pTyr.melted$PS_Multiplicity,
                                                  isoform.1.corrected.PS.position =HM.plotdata.pTyr.melted$isoform.1.corrected.PS.position,
                                                  isoform.1.corrected.Gene.Name_AApos = paste0(HM.plotdata.pTyr.melted$Gene.names, "_", HM.plotdata.pTyr.melted$Amino.acid, HM.plotdata.pTyr.melted$isoform.1.corrected.PS.position, "x", HM.plotdata.pTyr.melted$PS_Multiplicity),
                                                  isoform.2.corrected.PS.position =HM.plotdata.pTyr.melted$isoform.2.corrected.PS.position,
                                                  isoform.2.corrected.Gene.Name_AApos = paste0(HM.plotdata.pTyr.melted$Gene.names, "_", HM.plotdata.pTyr.melted$Amino.acid, HM.plotdata.pTyr.melted$isoform.2.corrected.PS.position, "x", HM.plotdata.pTyr.melted$PS_Multiplicity)
))
glimpse(order.df.HM.plotdata.pTyr.melted)




### for better visualization/plotting => collapse _1/_2/_3 via  sum
### note: sum is based on _1/_2/_3 intensities in dataframe at hand => unfiltered data might contain more _1/_2/_3 PS that are not included in the sum here

HM.plotdata.pTyr.collapsed <- HM.plotdata.pTyr
glimpse(HM.plotdata.pTyr.collapsed )

HM.plotdata.pTyr.collapsed$isoform.1.corrected.Gene.Name_AApos.collapsed <- pbsapply(HM.plotdata.pTyr.collapsed$isoform.1.corrected.Gene.Name_AApos, function(x) str_remove_all(string=x, pattern="x_\\d") )
HM.plotdata.pTyr.collapsed$isoform.2.corrected.Gene.Name_AApos.collapsed <- pbsapply(HM.plotdata.pTyr.collapsed$isoform.2.corrected.Gene.Name_AApos, function(x) str_remove_all(string=x, pattern="x_\\d") )
glimpse(HM.plotdata.pTyr.collapsed)
## unlog, do collapsing, relog
HM.plotdata.pTyr.collapsed [, c("GFP.1","GFP.2","GFP.3","FL.1","FL.2","FL.3","dE18.1","dE18.2","dE18.3","FL-Bicc1.1","FL-Bicc1.2","FL-Bicc1.3","dE18-Bicc1.1","dE18-Bicc1.2","dE18-Bicc1.3")] <- 2^HM.plotdata.pTyr.collapsed [, c("GFP.1","GFP.2","GFP.3","FL.1","FL.2","FL.3","dE18.1","dE18.2","dE18.3","FL-Bicc1.1","FL-Bicc1.2","FL-Bicc1.3","dE18-Bicc1.1","dE18-Bicc1.2","dE18-Bicc1.3")] 
glimpse(HM.plotdata.pTyr.collapsed )
##collapse
HM.plotdata.pTyr.collapsed.2  <- HM.plotdata.pTyr.collapsed %>% group_by(isoform.1.corrected.Gene.Name_AApos.collapsed, isoform.2.corrected.Gene.Name_AApos.collapsed) %>% summarise_at(c("GFP.1","GFP.2","GFP.3","FL.1","FL.2","FL.3","dE18.1","dE18.2","dE18.3","FL-Bicc1.1","FL-Bicc1.2","FL-Bicc1.3","dE18-Bicc1.1","dE18-Bicc1.2","dE18-Bicc1.3"), sum, na.rm = TRUE)
glimpse(HM.plotdata.pTyr.collapsed.2)
HM.plotdata.pTyr.collapsed.2[HM.plotdata.pTyr.collapsed.2  == 0] <- NA
##relog
HM.plotdata.pTyr.collapsed.2[c("GFP.1","GFP.2","GFP.3","FL.1","FL.2","FL.3","dE18.1","dE18.2","dE18.3","FL-Bicc1.1","FL-Bicc1.2","FL-Bicc1.3","dE18-Bicc1.1","dE18-Bicc1.2","dE18-Bicc1.3")] <- log2(HM.plotdata.pTyr.collapsed.2[c("GFP.1","GFP.2","GFP.3","FL.1","FL.2","FL.3","dE18.1","dE18.2","dE18.3","FL-Bicc1.1","FL-Bicc1.2","FL-Bicc1.3","dE18-Bicc1.1","dE18-Bicc1.2","dE18-Bicc1.3")] )
glimpse(HM.plotdata.pTyr.collapsed.2)


### reshape and reorder, define order for plot
HM.plotdata.pTyr.collapsed.melted <- reshape2::melt(HM.plotdata.pTyr.collapsed.2)
glimpse(HM.plotdata.pTyr.collapsed.melted)
colnames(HM.plotdata.pTyr.collapsed.melted ) <- c("isoform.1.corrected.Gene.Name_AApos.collapsed","isoform.2.corrected.Gene.Name_AApos.collapsed", "variable","log2.int")
glimpse(HM.plotdata.pTyr.collapsed.melted )


#define order for plot, get gene names, corrected positions, aminco acid etc. and order according to number
order.df.HM.plotdata.pTyr.collapsed.melted <- unique(tibble(isoform.1.corrected.Gene.Name_AApos.collapsed = HM.plotdata.pTyr.collapsed.melted$isoform.1.corrected.Gene.Name_AApos.collapsed,
                                                            isoform.2.corrected.Gene.Name_AApos.collapsed = HM.plotdata.pTyr.collapsed.melted$isoform.2.corrected.Gene.Name_AApos.collapsed))
order.df.HM.plotdata.pTyr.collapsed.melted 

order.df.HM.plotdata.pTyr.collapsed.melted$Gene <- pbsapply(order.df.HM.plotdata.pTyr.collapsed.melted$isoform.1.corrected.Gene.Name_AApos.collapsed , function(x) unlist(str_split(string=x,  pattern="_") )[1])

order.df.HM.plotdata.pTyr.collapsed.melted$isoform.1.corrected.AAPos <- pbsapply(order.df.HM.plotdata.pTyr.collapsed.melted$isoform.1.corrected.Gene.Name_AApos.collapsed, function(x) unlist(str_split(string=x,  pattern="_") )[2])
order.df.HM.plotdata.pTyr.collapsed.melted$isoform.2.corrected.AAPos <- pbsapply(order.df.HM.plotdata.pTyr.collapsed.melted$isoform.2.corrected.Gene.Name_AApos.collapsed, function(x) unlist(str_split(string=x,  pattern="_") )[2])

order.df.HM.plotdata.pTyr.collapsed.melted$isoform.corrected.Amino.acid <- pbsapply(order.df.HM.plotdata.pTyr.collapsed.melted$isoform.1.corrected.AAPos, function(x) unlist((str_extract_all(string=x,  pattern="[A-Z]") )))

order.df.HM.plotdata.pTyr.collapsed.melted$isoform.1.corrected.Pos <- pbsapply(order.df.HM.plotdata.pTyr.collapsed.melted$isoform.1.corrected.AAPos, function(x) unlist((str_extract_all(string=x,  pattern="[0-9]+") )))
order.df.HM.plotdata.pTyr.collapsed.melted$isoform.2.corrected.Pos <- pbsapply(order.df.HM.plotdata.pTyr.collapsed.melted$isoform.2.corrected.AAPos, function(x) unlist((str_extract_all(string=x,  pattern="[0-9]+") )))

order.df.HM.plotdata.pTyr.collapsed.melted <- order.df.HM.plotdata.pTyr.collapsed.melted %>% arrange(desc(isoform.2.corrected.Pos, isoform.2.corrected.Amino.acid))
order.df.HM.plotdata.pTyr.collapsed.melted


isoformcorrected.collapes.order.29.11.21 <- order.df.HM.plotdata.pTyr.collapsed.melted$isoform.2.corrected.Gene.Name_AApos.collapsed
isoformcorrected.collapes.order.29.11.21


glimpse(HM.plotdata.pTyr.collapsed.melted)


#plot the heatmap with site numbers for P21803-2
ggplot(
  HM.plotdata.pTyr.collapsed.melted , aes(x=variable, y=isoform.2.corrected.Gene.Name_AApos.collapsed, fill=log2.int)) + 
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
  geom_vline(xintercept = c(3.5, 6.5, 9.5, 12.5) , size=1.0, linetype = "solid", color="white")+
  annotate("text", x = 2.0, y = nrow(HM.plotdata.pTyr.collapsed.2)+2, label = "GFP", fontface = "bold", size=2.5)+
  annotate("text", x = 5.0, y = nrow(HM.plotdata.pTyr.collapsed.2)+2, label = "FL", fontface = "bold", size=2.5)+
  annotate("text", x = 8.0, y = nrow(HM.plotdata.pTyr.collapsed.2)+2, label = "dE18", fontface = "bold", size=2.5)+
  annotate("text", x = 11.0, y = nrow(HM.plotdata.pTyr.collapsed.2)+2, label = "FL\nBicc1", fontface = "bold", size=2.5)+
  annotate("text", x = 14.0, y = nrow(HM.plotdata.pTyr.collapsed.2)+2, label = "dE18\nBicc1", fontface = "bold", size=2.5)+
  annotate("text", x = 14.0, y = nrow(HM.plotdata.pTyr.collapsed.2)+3.5, label = "",  color = "transparent")

### add corresponding human FGFR2 phosphosite numbers to the heatmap plot ######################################################

glimpse(HM.plotdata.pTyr.collapsed.2)
order.df.HM.plotdata.pTyr.collapsed.melted

# Online Clustal Omega alignment showed that P21803-2 (707 aa) and NP_963895.2 (726 aa) are very similar. The difference is N-terminal: 
# 19aa  (MGLPSTWRYGRGPGIGTVT) are added for  NP_963895.2. Everything else is identical.
# Thus P21803-2 position +19 should be the positions for NP_963895.2

#prepare additional data frame and calculate the new site numbers
HM.plotdata.pTyr.collapsed.3 <- HM.plotdata.pTyr.collapsed.2

HM.plotdata.pTyr.collapsed.3$isoform.1.corrected.AApos <- pbsapply(HM.plotdata.pTyr.collapsed.3$isoform.1.corrected.Gene.Name_AApos.collapsed, function(x) unlist(str_split(string=x, pattern="_"))[2])
HM.plotdata.pTyr.collapsed.3$isoform.1.corrected.pos <- pbsapply(HM.plotdata.pTyr.collapsed.3$isoform.1.corrected.AApos, function(x) unlist(str_extract_all(string=x, pattern="[0-9]+")))

HM.plotdata.pTyr.collapsed.3$isoform.2.corrected.AApos <- pbsapply(HM.plotdata.pTyr.collapsed.3$isoform.2.corrected.Gene.Name_AApos.collapsed, function(x) unlist(str_split(string=x, pattern="_"))[2])
HM.plotdata.pTyr.collapsed.3$isoform.2.corrected.pos <- pbsapply(HM.plotdata.pTyr.collapsed.3$isoform.2.corrected.AApos, function(x) unlist(str_extract_all(string=x, pattern="[0-9]+")))
glimpse(HM.plotdata.pTyr.collapsed.3)


HM.plotdata.pTyr.collapsed.3$pos.NP_963895.2 <- as.numeric(HM.plotdata.pTyr.collapsed.3$isoform.2.corrected.pos) +19
glimpse(HM.plotdata.pTyr.collapsed.3)



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
HM.plotdata.pTyr.collapsed.4 <- merge(HM.plotdata.pTyr.collapsed.3, human.PSP.site.number.lookup.2[, c("site.Fgfr2.mouse.iso1", "site.FGFR2.human.iso1" )], by.x = "isoform.1.corrected.AApos", by.y = "site.Fgfr2.mouse.iso1", all.x = T, all.y = F, sort = F)
HM.plotdata.pTyr.collapsed.4$Gene <- pbsapply(HM.plotdata.pTyr.collapsed.4$isoform.1.corrected.Gene.Name_AApos.collapsed, function(x) unlist(str_split(string=x,  pattern="_") )[1])
HM.plotdata.pTyr.collapsed.4 <- HM.plotdata.pTyr.collapsed.4  %>% arrange(HM.plotdata.pTyr.collapsed.4$isoform.1.corrected.pos)

glimpse(HM.plotdata.pTyr.collapsed.4)

#some human sites have to be entered manually since not all sites are in PSP yet (looked up in online Clustal Omega alignment)
#they correspond mostly to numbers from P21803-1

unname(HM.plotdata.pTyr.collapsed.4$isoform.1.corrected.AApos)
HM.plotdata.pTyr.collapsed.4$site.FGFR2.human.iso1

HM.plotdata.pTyr.collapsed.4$AA.site.FGFR2.human.iso1_manual.change <- c("Y229", "Y237", "T243", "Y244", "S453", "S464", "Y466","Y561", "Y566", "Y575", "Y586", "S587", "Y588", "T607", "Y608", "Y616", "Y656", "Y657", "T660", "Y733", "Y805", "Y812", "S818")
HM.plotdata.pTyr.collapsed.4$site.FGFR2.human.iso1_manual.change <- pbsapply(HM.plotdata.pTyr.collapsed.4$AA.site.FGFR2.human.iso1_manual.change, function(x) str_remove_all(string=x, pattern="[A-Z]"))
glimpse(HM.plotdata.pTyr.collapsed.4)

#make.new.plot.label
#e.g. Fgfr2_Y229|114|133|229 => 1st number: P21803-1; 2nd number: P21803-2; 3rd number: NP_963895.2; 4th number: P21802-1
HM.plotdata.pTyr.collapsed.4$combined.plot.label <- paste0(HM.plotdata.pTyr.collapsed.4$isoform.1.corrected.Gene.Name_AApos.collapsed, "|", HM.plotdata.pTyr.collapsed.4$isoform.2.corrected.pos, "|",  HM.plotdata.pTyr.collapsed.4$pos.NP_963895.2, "|", HM.plotdata.pTyr.collapsed.4$site.FGFR2.human.iso1_manual.change)
glimpse(HM.plotdata.pTyr.collapsed.4)

#order for plot
isoform.corrected.collapsed.order.02.12.2021 <- rev(HM.plotdata.pTyr.collapsed.4$combined.plot.label)

# reshape for plot
melted.HM.plotdata.pTyr.collapsed.4 <- reshape2::melt(HM.plotdata.pTyr.collapsed.4 %>% select(combined.plot.label, GFP.1    : `dE18-Bicc1.3`))
glimpse(melted.HM.plotdata.pTyr.collapsed.4)

#plot
ggplot(
  melted.HM.plotdata.pTyr.collapsed.4 , aes(x=variable, y=combined.plot.label, fill=value)) + 
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
  theme(plot.title = element_text(hjust = 0.5, face="bold", size=16)) +
  theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1))+
  theme(axis.text.y= element_text(size=8))+ 
  xlab(NULL) + 
  xlab(NULL) + 
  ylab(NULL) +
  theme(axis.text.x = element_blank()) +
  theme(axis.ticks.x = element_blank())+
  geom_vline(xintercept = c(3.5, 6.5, 9.5, 12.5) , size=1.0, linetype = "solid", color="white")+
  annotate("text", x = 2.0, y = nrow(HM.plotdata.pTyr.collapsed.4)+2, label = "GFP", fontface = "bold", size=2.5)+
  annotate("text", x = 5.0, y = nrow(HM.plotdata.pTyr.collapsed.4)+2, label = "FL", fontface = "bold", size=2.5)+
  annotate("text", x = 8.0, y = nrow(HM.plotdata.pTyr.collapsed.4)+2, label = "dE18", fontface = "bold", size=2.5)+
  annotate("text", x = 11.0, y = nrow(HM.plotdata.pTyr.collapsed.4)+2, label = "FL\nBicc1", fontface = "bold", size=2.5)+
  annotate("text", x = 14.0, y = nrow(HM.plotdata.pTyr.collapsed.4)+2, label = "dE18\nBicc1", fontface = "bold", size=2.5)+
  annotate("text", x = 14.0, y = nrow(HM.plotdata.pTyr.collapsed.4)+3.5, label = "",  color = "transparent")

#ggsave("cells.pTyrIP.Fgfr2.sites.collapsed.with.isoform.numbers.pdf", useDingbats=FALSE,  width = 12, height =15, units = "cm") 



#############################################################################################################################################################
### cells pTyrIP correlation #################################################################################################################################
### cells pTyrIP correlation #################################################################################################################################

### correlation without Bicc1 sites 

glimpse(DZ.pTyr.IP.MBR.ON_expanded_log2_norm.class1.PS)
length(unique(DZ.pTyr.IP.MBR.ON_expanded_log2_norm.class1.PS$Sequence.window))

#define sites to exclude
sites_excluded_for_correlation_plot.cells.ip <- DZ.pTyr.IP.MBR.ON_expanded_log2_norm.class1.PS %>% filter(stringr::str_detect(Gene.Name_AApos, "Bicc1"))
glimpse(sites_excluded_for_correlation_plot.cells.ip) #47

#exclude Bicc1 sites
DF_for_correlation_plot.pTyr.IP.class1PS <- DZ.pTyr.IP.MBR.ON_expanded_log2_norm.class1.PS %>% filter(!stringr::str_detect(Gene.Name_AApos, "Bicc1")  )
glimpse(DF_for_correlation_plot.pTyr.IP.class1PS) #6751

#select columns and order
DF_for_correlation_plot.pTyr.IP.class1PS <-DF_for_correlation_plot.pTyr.IP.class1PS %>% select(contains("Norm.Intensity.pTyr")               )

glimpse(DF_for_correlation_plot.pTyr.IP.class1PS)



##check if there are all NA rows # 15 samples
DF_for_correlation_plot.pTyr.IP.class1PS$na.rows <- pbapply(DF_for_correlation_plot.pTyr.IP.class1PS, 1, function(x) sum(is.na(x)))
glimpse(DF_for_correlation_plot.pTyr.IP.class1PS)
table(DF_for_correlation_plot.pTyr.IP.class1PS$na.rows)


DF_for_correlation_plot.pTyr.IP.class1PS <- DF_for_correlation_plot.pTyr.IP.class1PS %>% select(-na.rows )
glimpse(DF_for_correlation_plot.pTyr.IP.class1PS)

colnames(DF_for_correlation_plot.pTyr.IP.class1PS) <- c("GFP.1",
                                                        "GFP.2",
                                                        "GFP.3",
                                                        "FL.1",
                                                        "FL.2",
                                                        "FL.3",
                                                        "dE18.1",
                                                        "dE18.2",
                                                        "dE18.3",
                                                        "FL.Bicc1.1",
                                                        "FL.Bicc1.2",
                                                        "FL.Bicc1.3",
                                                        "dE18.Bicc1.1",
                                                        "dE18.Bicc1.2",
                                                        "dE18.Bicc1.3")


#correlate
corr2 <- cor(DF_for_correlation_plot.pTyr.IP.class1PS, method = "pearson", use = "na.or.complete")
glimpse(corr2)
head(corr2)
min(corr2) #0.498
max(corr2)
round(min(corr2), 2)

# prepare annotation data
annot_df_for_heatmap <- data.frame(samples = colnames(DF_for_correlation_plot.pTyr.IP.class1PS),
                                   group = c(rep("GFP", 3), rep("FL", 3), rep("dE18", 3), rep("FL-Bicc1", 3), rep("dE18-Bicc1", 3)) 
                                   )

annot_df_for_heatmap 


# make shorter annotation data frame
annot_df_for_heatmap.short <- data.frame( 
  group  = annot_df_for_heatmap$group
)

glimpse(annot_df_for_heatmap.short)
annot_df_for_heatmap.short

# define colors for annotation bar
annot.colors = list(group = c("GFP"="yellow", "FL"="blue", "dE18"="red", "FL-Bicc1"="black", "dE18-Bicc1"="grey40"))
annot.colors


# create the heatmap annotation
ha.fr <- HeatmapAnnotation(group=annot_df_for_heatmap.short$group,
                           treatment=annot_df_for_heatmap.short$treatment,
                           col = annot.colors,
                           annotation_legend_param = list(grid_height = unit(8, "mm"))
)
ha.fr

# define colors
FR.heatmap.colors.2 <-colorRamp2(c(0.5, 1.0), c("grey90", "grey10"))



#pdf("cells.pTyrIP.correlation.pdf", width=19/2.54, height=17/2.54, useDingbats=FALSE) 
Heatmap(corr2, 
        name = "corr.coeff", 
        
        top_annotation = ha.fr,
        col = FR.heatmap.colors.2,
        
        clustering_distance_rows = function(m) as.dist((1-m)/2), #"euclidean"
        clustering_method_rows = "ward.D2",
        
        clustering_distance_columns = function(m) as.dist((1-m)/2), #"euclidean",
        clustering_method_columns = "ward.D2",
        
        column_dend_height = unit(40, "mm"),
        row_dend_width = unit(40, "mm"),
        
        heatmap_legend_param = list(ncol = 1, nrow = 1, legend_height = unit(60, "mm"))
        
)
#dev.off()





























