#_______________________________________________________________________________________________________________________
# 21.01.2022
# 
# Project: proteomics.truncated.FGFR2.is.oncogene.cancer.Zingg.et.al 
# 
# Script name: tumors.KB1P.global.phospho.public
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


###### define Perseus expand function ###############################################################################################################################


### define Perseus-like expand function
PerseusPSexpand.17.02.21 <- function(DF1){
  print("Wait a moment... ...")
  
  #get all column names data
  allDFcolNAMES <- colnames(DF1) #
  
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
### load data KB1P tumors global phospho DDA from MaxQuant Phospho(STY)Sites.txt

#/Users/frankrolfs/NLPostDR/NKI/EWA Phosphoproteomics with OPL/EWA Raw Data Thang version 10.2015/QE2_150602_OPL1005_EG_TiOx_BRCA_mouse Phospho (STY)Sites.txt
#/Users/frankrolfs/NLPostDR/NKI/EWA Phosphoproteomics with OPL/Ewa Tiox Evidence txt file/txt-tiox-28Sep2015/Phospho (STY)Sites.txt

##Tiox.Ewa.ms.tums <- fread("/Users/frankrolfs/NLPostDR/NKI/EWA Phosphoproteomics with OPL/Ewa Tiox Evidence txt file/txt-tiox-28Sep2015/Phospho (STY)Sites.txt",integer64 = "numeric") #
##glimpse(Tiox.Ewa.ms.tums)

#Tiox.Ewa.ms.tums <- fread("tumors.KB1P.Phospho(STY)Sites.txt",integer64 = "numeric") #
#glimpse(Tiox.Ewa.ms.tums)
#save(Tiox.Ewa.ms.tums, file = "Tiox.Ewa.ms.tums.Rdata")

load("Tiox.Ewa.ms.tums.Rdata")
glimpse(Tiox.Ewa.ms.tums)

nrow(Tiox.Ewa.ms.tums) #18681
ncol(Tiox.Ewa.ms.tums) #1011
length(unique(Tiox.Ewa.ms.tums$id)) #18681
glimpse(Tiox.Ewa.ms.tums)

### select columns of interest
Tiox.Ewa.ms.tums.2 <- Tiox.Ewa.ms.tums %>% select(
  
  "id",                          "Proteins",                    "Positions within proteins",   "Leading proteins",            "Protein",                    
  "Protein names",               "Gene names",                  "Fasta headers",               "Localization prob",           "Number of Phospho (STY)",    
  "Amino acid",                  "Sequence window",             "Phospho (STY) Probabilities", "Position in peptide",         "Reverse",                    
  "Potential contaminant",       "Positions",                   "Position",                    "Peptide IDs",                 "Mod. peptide IDs",      
  
  matches("___[0-9]"), -Intensity___1, -Intensity___2, -Intensity___3
  
)

glimpse(Tiox.Ewa.ms.tums.2)


### remove spaces from column names
column.names <- colnames(Tiox.Ewa.ms.tums.2)
column.names.changed <- unname(sapply(column.names , function(x) str_replace_all(string=x, pattern=" ", replacement=".")))
column.names.changed
colnames(Tiox.Ewa.ms.tums.2) <- column.names.changed

glimpse(Tiox.Ewa.ms.tums.2)


#change column classes
Tiox.Ewa.ms.tums.2 <- mutate_at(Tiox.Ewa.ms.tums.2,  21:ncol(Tiox.Ewa.ms.tums.2), list(as.numeric) )
glimpse(Tiox.Ewa.ms.tums.2)


### remove reverse hits
table(Tiox.Ewa.ms.tums.2$Reverse) #56
Tiox.Ewa.ms.tums.3 <- Tiox.Ewa.ms.tums.2[Tiox.Ewa.ms.tums.2$Reverse =="",]
nrow(Tiox.Ewa.ms.tums.3) #18625


### remove contaminants
table(Tiox.Ewa.ms.tums.3$Potential.contaminant) #145+
Tiox.Ewa.ms.tums.3 <- Tiox.Ewa.ms.tums.3[Tiox.Ewa.ms.tums.3$Potential.contaminant == "",]
nrow(Tiox.Ewa.ms.tums.3) #18480


### remove rows completely zero
colnames(Tiox.Ewa.ms.tums.3)
ncol(Tiox.Ewa.ms.tums.3) #308
Tiox.Ewa.ms.tums.3$zero.check <- pbapply(Tiox.Ewa.ms.tums.3[, c(21:308)], 1, function(x) sum(x) ) #
table(Tiox.Ewa.ms.tums.3$zero.check >0) #F3819 T14661

Tiox.Ewa.ms.tums.3.complete.zero.rows <- filter(Tiox.Ewa.ms.tums.3, zero.check <=0)
glimpse(Tiox.Ewa.ms.tums.3.complete.zero.rows) #

Tiox.Ewa.ms.tums.3 <- filter(Tiox.Ewa.ms.tums.3, zero.check > 0)
nrow(Tiox.Ewa.ms.tums.3) #14661
glimpse(Tiox.Ewa.ms.tums.3)
table(Tiox.Ewa.ms.tums.3$zero.check >0)


length(unique(Tiox.Ewa.ms.tums.3$id))#14661 = total number PS
table(Tiox.Ewa.ms.tums.3$Amino.acid) #S12644 T1903 Y114
length(unique(Tiox.Ewa.ms.tums.3$Sequence.window)) #14656

### remove column zero.check
Tiox.Ewa.ms.tums.3 <- Tiox.Ewa.ms.tums.3 %>% select(-zero.check)
glimpse(Tiox.Ewa.ms.tums.3)


### create variable with short sample names with order as in Phospho(STY)Sites.txt
samples.in.experiment <- colnames(Tiox.Ewa.ms.tums.3)
samples.in.experiment <- samples.in.experiment[21:ncol(Tiox.Ewa.ms.tums.3)]
samples.in.experiment <- str_remove(string=samples.in.experiment, pattern="Intensity.")
samples.in.experiment <- str_remove(string=samples.in.experiment, pattern="___[0-9]")
samples.in.experiment <- unique(samples.in.experiment)
samples.in.experiment

################################################################################################################################################################
################################################################################################################################################################

#### expand with PerseusPSexpand function (see above)################################################################################################################

glimpse(Tiox.Ewa.ms.tums.3)

# exclude sample 31 => bad quality
Tiox.Ewa.ms.tums.3 <- Tiox.Ewa.ms.tums.3 %>% select(-Intensity.TiOx_31___1, -Intensity.TiOx_31___2, -Intensity.TiOx_31___3)
glimpse(Tiox.Ewa.ms.tums.3)

#expand
Tiox.Ewa.ms.tums_expanded <- PerseusPSexpand.17.02.21(Tiox.Ewa.ms.tums.3)

glimpse(Tiox.Ewa.ms.tums_expanded )

#remove rows that have all zero intensities
glimpse(Tiox.Ewa.ms.tums_expanded)
colnames(Tiox.Ewa.ms.tums_expanded)

Tiox.Ewa.ms.tums_expanded$ZEROcheckINT <- pbapply(Tiox.Ewa.ms.tums_expanded[, c(21:115)], MARGIN=1, function(x) sum(x)) #
table(Tiox.Ewa.ms.tums_expanded$ZEROcheckINT > 0)
nrow(Tiox.Ewa.ms.tums_expanded) #43983

Tiox.Ewa.ms.tums_expanded <- filter(Tiox.Ewa.ms.tums_expanded, ZEROcheckINT > 0)
nrow(Tiox.Ewa.ms.tums_expanded) #17587

#remove column ZEROcheckINT
Tiox.Ewa.ms.tums_expanded <- select(Tiox.Ewa.ms.tums_expanded,  -ZEROcheckINT)
glimpse(Tiox.Ewa.ms.tums_expanded)



# normalization & histogram standard median shift version #########################################################################################################################################
# normalization & histogram standard median shift version #########################################################################################################################################

glimpse(Tiox.Ewa.ms.tums_expanded)
colnames(Tiox.Ewa.ms.tums_expanded)

#log2 transform & replace -Inf with NA
colnames(Tiox.Ewa.ms.tums_expanded)
Log2Trafo <- Tiox.Ewa.ms.tums_expanded[, c(21:115)] #
Log2Trafo <- log2(Log2Trafo)
Log2Trafo[Log2Trafo  == -Inf] <- NA

glimpse(Log2Trafo)

Tiox.Ewa.ms.tums_expanded_log2 <- Tiox.Ewa.ms.tums_expanded
Tiox.Ewa.ms.tums_expanded_log2[, c(21:115)] <- Log2Trafo #

glimpse(Tiox.Ewa.ms.tums_expanded_log2, list.len=900)


# calculate median per column RawINT
glimpse(Tiox.Ewa.ms.tums_expanded_log2)
colnames(Tiox.Ewa.ms.tums_expanded_log2)

ColMedianfindDF <- glimpse(Tiox.Ewa.ms.tums_expanded_log2[, c(21:115)]) #
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
Mdf <- data.frame(variable=colnames(Tiox.Ewa.ms.tums_expanded_log2[, c(21:115)]), median=vectormediancolumnsINT); Mdf
Mdf

## Plot
DensityPlotData <- Tiox.Ewa.ms.tums_expanded_log2[, c(21:115)]
glimpse(DensityPlotData)

DensityPlotDataMelted <- reshape2::melt(DensityPlotData)

ggplot(DensityPlotDataMelted, aes(x=value, colour = variable) ) + 
  geom_density() +
  theme(legend.position="bottom") +
  geom_vline(data=Mdf, aes(xintercept=median, color=variable))


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
temp_for_median_shift <-  Tiox.Ewa.ms.tums_expanded_log2[, c(21:115)]
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
Tiox.Ewa.ms.tums_expanded_log2_norm <- Tiox.Ewa.ms.tums_expanded_log2
Tiox.Ewa.ms.tums_expanded_log2_norm  <- cbind(Tiox.Ewa.ms.tums_expanded_log2_norm , temp_for_median_shift)
glimpse(Tiox.Ewa.ms.tums_expanded_log2_norm, list.len=999)


### box plot after normalization
data.for.boxplot.after.norm.DZ.pTyr.IP.MBR.ON <- Tiox.Ewa.ms.tums_expanded_log2_norm %>% select(contains("Norm.Intensity.TiOx_") )
glimpse(data.for.boxplot.after.norm.DZ.pTyr.IP.MBR.ON)

melted.data.for.boxplot.after.norm.DZ.pTyr.IP.MBR.ON <- reshape2::melt(data.for.boxplot.after.norm.DZ.pTyr.IP.MBR.ON)
glimpse(melted.data.for.boxplot.after.norm.DZ.pTyr.IP.MBR.ON )

ggplot(data = melted.data.for.boxplot.after.norm.DZ.pTyr.IP.MBR.ON)+
  geom_boxplot(mapping = aes(x= variable, y = value), outlier.size = 0.5, notch = TRUE, notchwidth = 0.1)+
  scale_fill_viridis(discrete=T, option="inferno")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1.0, vjust = 0.5))+
  theme(legend.position="none", legend.justification = "center")+
  ggtitle("tumors KB1P global phospho - normalized data - PS level") +
  theme(plot.title = element_text(hjust = 0.5, face= "bold"))+  #center title
  theme(axis.text.x = element_text(face= "plain", colour="black", size=14))+
  theme(axis.text.y = element_text(face= "plain", colour="black", size=14))+
  xlab(NULL)+ 
  ylab("log2(Int.)") 



### box plot WITHOUT normalization
glimpse(Tiox.Ewa.ms.tums_expanded_log2_norm)

data.for.boxplot.after.norm.DZ.pTyr.IP.MBR.ON <- Tiox.Ewa.ms.tums_expanded_log2_norm %>% select( Intensity.TiOx_01 : Intensity.TiOx_96)
glimpse(data.for.boxplot.after.norm.DZ.pTyr.IP.MBR.ON)

melted.data.for.boxplot.after.norm.DZ.pTyr.IP.MBR.ON <- reshape2::melt(data.for.boxplot.after.norm.DZ.pTyr.IP.MBR.ON)

ggplot(data = melted.data.for.boxplot.after.norm.DZ.pTyr.IP.MBR.ON)+
  geom_boxplot(mapping = aes(x= variable, y = value), outlier.size = 0.5, notch = TRUE, notchwidth = 0.1)+ 
  scale_x_discrete(limits= c())+
  theme(axis.text.x = element_text(angle = 90, hjust = 1.0, vjust = 0.5))+
  theme(legend.position="none", legend.justification = "center")+
  ggtitle("tumors KB1P global phospho - NOT normalized data - PS level") +
  theme(plot.title = element_text(hjust = 0.5, face= "bold"))+  #center title
  theme(axis.text.x = element_text(face= "plain", colour="black", size=14))+
  theme(axis.text.y = element_text(face= "plain", colour="black", size=14))+
  xlab(NULL)+ 
  ylab("log2(Int.)") 


#####################################################################################################################################################################
#### analysis via class 1 PS ##################################################################################################################################

#filter class 1
Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS <- Tiox.Ewa.ms.tums_expanded_log2_norm %>% filter(Localization.prob>= 0.75)
glimpse(Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS) #13771
length(unique(Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS$Sequence.window)) #11132

# add additional column to better calculate number PS 
# take Proteins and add Psite without multiplicity info ###
for(i in 1:nrow(Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS)){
  print(i)
  temp_prot_name <- Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS$Proteins[i]
  tempPSite <- Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS$Positions[i]
  tempAA <- Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS$Amino.acid[i]
  Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS$Prot.Name_AApos[i] <- paste0(temp_prot_name, "_", tempAA, tempPSite)
}
length(unique(Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS$Prot.Name_AApos)) #11136 == the number of class 1 PS in exp

glimpse(Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS)

#  prepare new column
# take genename and add Psite ###
for(i in 1:nrow(Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS)){
  print(i)
  temp_gene_name <- Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS$Gene.names[i]
  tempPSite <- Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS$Positions[i]
  tempAA <- Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS$Amino.acid[i]
  tempPSmultiplicity <- Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS$PS_Multiplicity[i]
  Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS$Gene.Name_AApos[i] <- paste0(temp_gene_name, "_", tempAA, tempPSite, "x", tempPSmultiplicity)
}
glimpse(Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS)


#add identifier for plotting reason e.g. PTMSEA later
glimpse(Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS)

Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS$ID.FR.all.C1.PS <- 1:nrow(Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS)
glimpse(Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS)


#save data
#save(Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS, file="Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS.Rdata")
#load("Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS")
#glimpse(Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS)


### average replicates outside log space #######################################################################################################################

###START define function to average replicate values version for NA####
average.replicatesEWATiox.withNAs.04.11.16_A <- function(arg){
  arg <- as.vector(arg, mode = "numeric")
  
  if( !is.na(arg[1]) & !is.na(arg[2]) ) {
    output <- (arg[1] + arg[2])/2
  } else if( is.na(arg[1]) & is.na(arg[2]) ){ 
    output <- (arg[1] + arg[2])/2
  } else if( is.na(arg[1])& arg[2] > 0 )  {
    output <- arg[2]
  } else if( is.na(arg[2]) & arg[1] > 0) {
    output <- arg[1]
  } else {
    output <- "error"
  }
  return(output)  
}
###END define function to average replicate values version for NA####

glimpse(Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS)
colnames(Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS)

#unlog for averaging
temp.data.only <- Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS[, c(117 : 211)]
temp.data.only <- 2^temp.data.only
glimpse(temp.data.only)

Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS[, c(117 : 211)] <- temp.data.only

glimpse(Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS)
colnames(Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS)

#average
Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS$Int.1_2 <- apply(Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS[, c(117,118)], 1, function(x) average.replicatesEWATiox.withNAs.04.11.16_A(x) ) 
Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS$Int.3_4 <- apply(Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS[, c(119,120)], 1, function(x) average.replicatesEWATiox.withNAs.04.11.16_A(x) ) 
Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS$Int.5_6 <- apply(Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS[, c(121,122)], 1, function(x) average.replicatesEWATiox.withNAs.04.11.16_A(x) ) 
Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS$Int.7_8 <- apply(Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS[, c(123,124)], 1, function(x) average.replicatesEWATiox.withNAs.04.11.16_A(x) ) 
Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS$Int.9_10 <- apply(Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS[, c(125,126)], 1, function(x) average.replicatesEWATiox.withNAs.04.11.16_A(x) ) 

Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS$Int.11_12 <- apply(Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS[, c(127,128)], 1, function(x) average.replicatesEWATiox.withNAs.04.11.16_A(x) ) 
Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS$Int.13_14 <- apply(Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS[, c(129,130)], 1, function(x) average.replicatesEWATiox.withNAs.04.11.16_A(x) ) 
Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS$Int.15_16 <- apply(Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS[, c(131,132)], 1, function(x) average.replicatesEWATiox.withNAs.04.11.16_A(x) ) 
Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS$Int.17_18 <- apply(Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS[, c(133,134)], 1, function(x) average.replicatesEWATiox.withNAs.04.11.16_A(x) ) 

Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS$Int.19_20 <- apply(Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS[, c(135,136)], 1, function(x) average.replicatesEWATiox.withNAs.04.11.16_A(x) ) 
Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS$Int.21_22 <- apply(Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS[, c(137,138)], 1, function(x) average.replicatesEWATiox.withNAs.04.11.16_A(x) ) 
Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS$Int.23_24 <- apply(Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS[, c(139,140)], 1, function(x) average.replicatesEWATiox.withNAs.04.11.16_A(x) ) 
Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS$Int.25_26 <- apply(Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS[, c(141,142)], 1, function(x) average.replicatesEWATiox.withNAs.04.11.16_A(x) ) 
Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS$Int.27_28 <- apply(Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS[, c(143,144)], 1, function(x) average.replicatesEWATiox.withNAs.04.11.16_A(x) ) 

Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS$Int.29_30 <- apply(Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS[, c(145,146)], 1, function(x) average.replicatesEWATiox.withNAs.04.11.16_A(x) ) 
Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS$Int.32_32 <- apply(Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS[, c(147,147)], 1, function(x) average.replicatesEWATiox.withNAs.04.11.16_A(x) ) 
Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS$Int.33_34 <- apply(Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS[, c(148,149)], 1, function(x) average.replicatesEWATiox.withNAs.04.11.16_A(x) ) 
Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS$Int.35_36 <- apply(Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS[, c(150,151)], 1, function(x) average.replicatesEWATiox.withNAs.04.11.16_A(x) ) 
Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS$Int.37_38 <- apply(Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS[, c(152,153)], 1, function(x) average.replicatesEWATiox.withNAs.04.11.16_A(x) ) 

Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS$Int.39_40 <- apply(Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS[, c(154,155)], 1, function(x) average.replicatesEWATiox.withNAs.04.11.16_A(x) ) 
Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS$Int.41_42 <- apply(Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS[, c(156,157)], 1, function(x) average.replicatesEWATiox.withNAs.04.11.16_A(x) ) 
Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS$Int.43_44 <- apply(Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS[, c(158,159)], 1, function(x) average.replicatesEWATiox.withNAs.04.11.16_A(x) ) 
Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS$Int.45_46 <- apply(Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS[, c(160,161)], 1, function(x) average.replicatesEWATiox.withNAs.04.11.16_A(x) ) 
Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS$Int.47_48 <- apply(Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS[, c(162,163)], 1, function(x) average.replicatesEWATiox.withNAs.04.11.16_A(x) ) 

Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS$Int.49_50 <- apply(Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS[, c(164,165)], 1, function(x) average.replicatesEWATiox.withNAs.04.11.16_A(x) ) 
Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS$Int.51_52 <- apply(Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS[, c(166,167)], 1, function(x) average.replicatesEWATiox.withNAs.04.11.16_A(x) ) 
Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS$Int.53_54 <- apply(Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS[, c(168,169)], 1, function(x) average.replicatesEWATiox.withNAs.04.11.16_A(x) ) 
Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS$Int.55_56 <- apply(Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS[, c(170,171)], 1, function(x) average.replicatesEWATiox.withNAs.04.11.16_A(x) ) 
Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS$Int.57_58 <- apply(Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS[, c(172,173)], 1, function(x) average.replicatesEWATiox.withNAs.04.11.16_A(x) ) 

Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS$Int.59_60 <- apply(Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS[, c(174,175)], 1, function(x) average.replicatesEWATiox.withNAs.04.11.16_A(x) ) 
Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS$Int.61_62 <- apply(Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS[, c(176,177)], 1, function(x) average.replicatesEWATiox.withNAs.04.11.16_A(x) ) 
Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS$Int.63_64 <- apply(Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS[, c(178,179)], 1, function(x) average.replicatesEWATiox.withNAs.04.11.16_A(x) ) 
Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS$Int.65_66 <- apply(Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS[, c(180,181)], 1, function(x) average.replicatesEWATiox.withNAs.04.11.16_A(x) ) 
Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS$Int.67_68 <- apply(Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS[, c(182,183)], 1, function(x) average.replicatesEWATiox.withNAs.04.11.16_A(x) ) 

Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS$Int.69_70 <- apply(Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS[, c(184,185)], 1, function(x) average.replicatesEWATiox.withNAs.04.11.16_A(x) )  
Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS$Int.71_72 <- apply(Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS[, c(186,187)], 1, function(x) average.replicatesEWATiox.withNAs.04.11.16_A(x) ) 
Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS$Int.73_74 <- apply(Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS[, c(188,189)], 1, function(x) average.replicatesEWATiox.withNAs.04.11.16_A(x) ) 
Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS$Int.75_76 <- apply(Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS[, c(190,191)], 1, function(x) average.replicatesEWATiox.withNAs.04.11.16_A(x) ) 
Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS$Int.77_78 <- apply(Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS[, c(192,193)], 1, function(x) average.replicatesEWATiox.withNAs.04.11.16_A(x) ) 

Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS$Int.79_80 <- apply(Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS[, c(194,195)], 1, function(x) average.replicatesEWATiox.withNAs.04.11.16_A(x) ) 
Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS$Int.81_82 <- apply(Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS[, c(196,197)], 1, function(x) average.replicatesEWATiox.withNAs.04.11.16_A(x) ) 
Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS$Int.83_84 <- apply(Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS[, c(198,199)], 1, function(x) average.replicatesEWATiox.withNAs.04.11.16_A(x) ) 
Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS$Int.85_86 <- apply(Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS[, c(200,201)], 1, function(x) average.replicatesEWATiox.withNAs.04.11.16_A(x) ) 
Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS$Int.87_88 <- apply(Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS[, c(202,203)], 1, function(x) average.replicatesEWATiox.withNAs.04.11.16_A(x) ) 

Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS$Int.89_90 <- apply(Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS[, c(204,205)], 1, function(x) average.replicatesEWATiox.withNAs.04.11.16_A(x) ) 
Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS$Int.91_92 <- apply(Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS[, c(206,207)], 1, function(x) average.replicatesEWATiox.withNAs.04.11.16_A(x) ) 
Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS$Int.93_94 <- apply(Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS[, c(208,209)], 1, function(x) average.replicatesEWATiox.withNAs.04.11.16_A(x) ) 
Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS$Int.95_96 <- apply(Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS[, c(210,211)], 1, function(x) average.replicatesEWATiox.withNAs.04.11.16_A(x) ) 

glimpse(Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS) # data averaged outside log space
colnames(Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS) # data averaged outside log space

# relog
temp.data.only <- Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS[, c(215 : 262)]
temp.data.only <- log2(temp.data.only)
glimpse(temp.data.only)

Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS[, c(215 : 262)] <- temp.data.only

glimpse(Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS) # 
colnames(Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS) #




##############################################################################################################################################################
##############################################################################################################################################################
##############################################################################################################################################################
#tumors KB1P global phospho:  single sample PTMSEA (ssPTMSEA) using mouse PTMSigDB

glimpse(Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS)


#select samples of interest (24)
Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS_ssPTMSEA <- Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS  %>% select("id": "Mod..peptide.IDs","PS_Multiplicity","Prot.Name_AApos","Gene.Name_AApos","ID.FR.all.C1.PS",
                                                                                                                  "Int.1_2", "Int.3_4", "Int.9_10","Int.11_12","Int.17_18","Int.19_20","Int.25_26","Int.27_28", "Int.33_34","Int.35_36", "Int.41_42","Int.43_44",
                                                                                                                  "Int.65_66","Int.67_68","Int.73_74","Int.75_76", "Int.81_82","Int.83_84","Int.89_90","Int.91_92", "Int.57_58","Int.59_60",
                                                                                                                  "Int.49_50","Int.51_52")
glimpse(Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS_ssPTMSEA )



##sample selection
# "Int.1_2", "Int.3_4", "Int.9_10","Int.11_12","Int.17_18","Int.19_20","Int.25_26","Int.27_28", "Int.33_34","Int.35_36", "Int.41_42","Int.43_44"
# "Int.65_66","Int.67_68","Int.73_74","Int.75_76", "Int.81_82","Int.83_84","Int.89_90","Int.91_92", 
# "Int.49_50","Int.51_52", "Int.57_58","Int.59_60"


### prepare data for ssPTMSEA 

## start define function to make aa15.window
make.aa.15.window.from.31.sequence.window <- function(input){
  nchar(input)
  input <- unlist(str_split(string=input, pattern="")); input
  return(paste0(input[9:23], collapse="") )
}
make.aa.15.window.from.31.sequence.window(input="RLETSTSCFYQPQRRsVILDGRSGRQIE___")
make.aa.15.window.from.31.sequence.window(input="RLETSTSCFYQPQRRsVILDGRSGRQIE___")
## stop define function to make aa15.window


# make aa15.window
# take first entry sequence window if there are two entries
Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS_ssPTMSEA$aa15.window <- pbsapply(Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS_ssPTMSEA$Sequence.window, function(x) unlist(str_split(string =x, pattern = ";"))[1]  )
Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS_ssPTMSEA$aa15.window <- pbsapply(Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS_ssPTMSEA$aa15.window, function(x) make.aa.15.window.from.31.sequence.window(input = x)  )

glimpse(Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS_ssPTMSEA)

table(duplicated(Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS_ssPTMSEA$aa15.window))#2672
which(duplicated(Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS_ssPTMSEA$aa15.window))


Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS_ssPTMSEA_2 <- Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS_ssPTMSEA
glimpse(Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS_ssPTMSEA_2  )


# find unique aa15.windows
unique.15aa.windows <- unique(Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS_ssPTMSEA_2 $aa15.window)
glimpse(unique.15aa.windows) #11099

# PTMSEA does not accept duplicate aa15.windows
# if there are duplicate entries (e.g. because of _1_2_3) here, use the entry with the highest sum intensity over all samples, 
# this will favour entries with low NAs, good signal an can possibly distinguish duplicate situations with the same number of data points present
Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS_ssPTMSEA_3 <- c()
for (i in unique.15aa.windows){
  print(i)
  temp <-Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS_ssPTMSEA_2 %>% filter(aa15.window %in% c(i))
  temp$temp.sum<-apply(temp[c(
    "Int.1_2", "Int.3_4", "Int.9_10","Int.11_12","Int.17_18","Int.19_20","Int.25_26","Int.27_28", "Int.33_34","Int.35_36", "Int.41_42","Int.43_44",
    "Int.65_66","Int.67_68","Int.73_74","Int.75_76", "Int.81_82","Int.83_84","Int.89_90","Int.91_92", "Int.57_58","Int.59_60",
    "Int.49_50","Int.51_52")], 1, function(x) sum(x, na.rm=T))
  temp <- temp %>% arrange(desc(temp.sum))
  temp <- temp[1,]
  Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS_ssPTMSEA_3 <- bind_rows(Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS_ssPTMSEA_3 , temp)
}

glimpse(Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS_ssPTMSEA_3)
table(duplicated(Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS_ssPTMSEA_3$aa15.window))

# add_p to end of 15aa window
Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS_ssPTMSEA_3$aa15.window.2 <- paste0(Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS_ssPTMSEA_3$aa15.window, "-p")
glimpse(Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS_ssPTMSEA_3) 



### prepare input for PTMSEA via https://cloud.genepattern.org
gene.pattern.PTM.SEA.input <- Tiox.Ewa.ms.tums_expanded_log2_normclass1.PS_ssPTMSEA_3 %>% dplyr::select(aa15.window.2, 
                                                                                                        "Int.1_2", "Int.3_4", "Int.9_10","Int.11_12","Int.17_18","Int.19_20","Int.25_26","Int.27_28", "Int.33_34","Int.35_36", "Int.41_42","Int.43_44",
                                                                                                        "Int.65_66","Int.67_68","Int.73_74","Int.75_76", "Int.81_82","Int.83_84","Int.89_90","Int.91_92", "Int.57_58","Int.59_60",
                                                                                                        "Int.49_50","Int.51_52")
glimpse(gene.pattern.PTM.SEA.input)
nrow(gene.pattern.PTM.SEA.input) #11099


#impute a zero (e.g. for black and white regulation situations) 
gene.pattern.PTM.SEA.input[is.na(gene.pattern.PTM.SEA.input)] <- 0
glimpse(gene.pattern.PTM.SEA.input)


#prepare gct input file with cmapR package
gene.pattern.PTM.SEA.input_GCT <- as.matrix(gene.pattern.PTM.SEA.input[,c(
  "Int.1_2", "Int.3_4", "Int.9_10","Int.11_12","Int.17_18","Int.19_20","Int.25_26","Int.27_28", "Int.33_34","Int.35_36", "Int.41_42","Int.43_44",
  "Int.65_66","Int.67_68","Int.73_74","Int.75_76", "Int.81_82","Int.83_84","Int.89_90","Int.91_92", "Int.57_58","Int.59_60",
  "Int.49_50","Int.51_52")])
glimpse(gene.pattern.PTM.SEA.input_GCT)
colnames(gene.pattern.PTM.SEA.input_GCT) <- c(
  "Int.1_2", "Int.3_4", "Int.9_10","Int.11_12","Int.17_18","Int.19_20","Int.25_26","Int.27_28", "Int.33_34","Int.35_36", "Int.41_42","Int.43_44",
  "Int.65_66","Int.67_68","Int.73_74","Int.75_76", "Int.81_82","Int.83_84","Int.89_90","Int.91_92", "Int.57_58","Int.59_60",
  "Int.49_50","Int.51_52")
rownames(gene.pattern.PTM.SEA.input_GCT) <- gene.pattern.PTM.SEA.input$aa15.window.2

glimpse(gene.pattern.PTM.SEA.input_GCT)


gene.pattern.PTM.SEA.input_GCT.2 <- new("GCT", mat=gene.pattern.PTM.SEA.input_GCT)
glimpse(gene.pattern.PTM.SEA.input_GCT.2)

#save GCT file for PTMSEA input
write_gct(gene.pattern.PTM.SEA.input_GCT.2, "KB1P.tums.class1.PS.ssPTMSEA.input.gct")


### run PTMSEA via Gene Pattern at https://cloud.genepattern.org ################################################################
# pathwaydb: ptm.sig.db.all.flanking.mouse.v1.9.0.gmt see https://github.com/broadinstitute/ssGSEA2.0
# all settings default
# save results and ...combined.gct file
# change ...combined.gct file type to .txt and remove first two rows


# result ssPTMSEA stored under "KB1P.tumors.ssPTMSEA.result.combined.txt" 
# compare script "tumors.global.phospho.IMAC.public.R" for the combination with FGFR2 IMAC tumors and the corresponding plot


##################################################################################################################################
##################################################################################################################################

(.packages())
#[1] "ggprism"         "gtools"          "openxlsx"        "ggsci"           "paletteer"       "readxl"          "scico"          
#[8] "reshape2"        "splitstackshape" "viridis"         "viridisLite"     "circlize"        "ComplexHeatmap"  "grid"           
#[15] "pbapply"         "data.table"      "forcats"         "stringr"         "dplyr"           "purrr"           "readr"          
#[22] "tidyr"           "tibble"          "ggplot2"         "tidyverse"       "cmapR"           "limma"           "stats"          
#[29] "graphics"        "grDevices"       "utils"           "datasets"        "methods"         "base"           

sessionInfo()

# R version 4.1.2 (2021-11-01)
# Platform: aarch64-apple-darwin20 (64-bit)
# Running under: macOS Monterey 12.1
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
#   [1] ggprism_1.0.3         gtools_3.9.2          openxlsx_4.2.5        ggsci_2.9             paletteer_1.4.0      
# [6] readxl_1.3.1          scico_1.3.0           reshape2_1.4.4        splitstackshape_1.4.8 viridis_0.6.2        
# [11] viridisLite_0.4.0     circlize_0.4.13       ComplexHeatmap_2.10.0 pbapply_1.5-0         data.table_1.14.2    
# [16] forcats_0.5.1         stringr_1.4.0         dplyr_1.0.7           purrr_0.3.4           readr_2.1.1          
# [21] tidyr_1.1.4           tibble_3.1.6          ggplot2_3.3.5         tidyverse_1.3.1       cmapR_1.6.0          
# [26] limma_3.50.0         
# 
# loaded via a namespace (and not attached):
#   [1] colorspace_2.0-2            rjson_0.2.20                ellipsis_0.3.2              cytolib_2.6.1              
# [5] XVector_0.34.0              GenomicRanges_1.46.1        GlobalOptions_0.1.2         fs_1.5.2                   
# [9] dichromat_2.0-0             clue_0.3-60                 rstudioapi_0.13             farver_2.1.0               
# [13] bit64_4.0.5                 fansi_0.5.0                 lubridate_1.8.0             xml2_1.3.3                 
# [17] codetools_0.2-18            doParallel_1.0.16           jsonlite_1.7.2              broom_0.7.11               
# [21] cluster_2.1.2               dbplyr_2.1.1                png_0.1-7                   mapproj_1.2.7              
# [25] compiler_4.1.2              httr_1.4.2                  backports_1.4.1             assertthat_0.2.1           
# [29] Matrix_1.4-0                cli_3.1.0                   tools_4.1.2                 gtable_0.3.0               
# [33] glue_1.6.0                  GenomeInfoDbData_1.2.7      maps_3.4.0                  Rcpp_1.0.7                 
# [37] Biobase_2.54.0              cellranger_1.1.0            vctrs_0.3.8                 iterators_1.0.13           
# [41] rvest_1.0.2                 lifecycle_1.0.1             zlibbioc_1.40.0             scales_1.1.1               
# [45] vroom_1.5.7                 RProtoBufLib_2.6.0          hms_1.1.1                   MatrixGenerics_1.6.0       
# [49] parallel_4.1.2              SummarizedExperiment_1.24.0 rematch2_2.1.2              RColorBrewer_1.1-2         
# [53] prismatic_1.1.0             gridExtra_2.3               stringi_1.7.6               S4Vectors_0.32.3           
# [57] foreach_1.5.1               flowCore_2.6.0              BiocGenerics_0.40.0         zip_2.2.0                  
# [61] shape_1.4.6                 GenomeInfoDb_1.30.0         pals_1.7                    rlang_0.4.12               
# [65] pkgconfig_2.0.3             matrixStats_0.61.0          bitops_1.0-7                lattice_0.20-45            
# [69] labeling_0.4.2              bit_4.0.4                   tidyselect_1.1.1            plyr_1.8.6                 
# [73] magrittr_2.0.1              R6_2.5.1                    IRanges_2.28.0              generics_0.1.1             
# [77] DelayedArray_0.20.0         DBI_1.1.2                   pillar_1.6.4                haven_2.4.3                
# [81] withr_2.4.3                 RCurl_1.98-1.5              modelr_0.1.8                crayon_1.4.2               
# [85] utf8_1.2.2                  tzdb_0.2.0                  GetoptLong_1.0.5            reprex_2.0.1               
# [89] digest_0.6.29               RcppParallel_5.1.5          stats4_4.1.2                munsell_0.5.0  


