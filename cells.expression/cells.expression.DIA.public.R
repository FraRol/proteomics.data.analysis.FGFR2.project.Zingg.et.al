#_______________________________________________________________________________________________________________________
# 12.01.2022
# 
# Project: proteomics.truncated.FGFR2.is.oncogene.cancer.Zingg.et.al 
# 
# Script name: cells.expression.DIA.public
#
# Purpose of script: analysis protein expression data cells truncated FGFR2 is oncogene cancer project Zingg et al 
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
library(cmapR)
library(msigdbr)
library(gtools)
library(paletteer)

###################################################################################################################################

### load expresion data cells DIA Zingg
load("Zingg.cells.expression.DIA.Spectronaut.report.long.Rdata")
glimpse(DZ.cells.expr)

### changes for sample names
DZ.cells.expr <- mutate(DZ.cells.expr, Condition=R.FileName) #
glimpse(DZ.cells.expr)

### select less columns: PG level only
DZ.cells.expr <- select(DZ.cells.expr, Condition, contains("PG")) %>% distinct()
glimpse(DZ.cells.expr)
nrow(DZ.cells.expr) #358711

### detect PG.ProteinGroups with "CON__" and remove them
# CON entries do not have data like PG.Quantity, or Q.value ...
DZ.cells.expr$contamination <- pbsapply(DZ.cells.expr$PG.ProteinGroups, function(x) if(str_detect(string= x, pattern="CON__")){"CON"} else {"ok"})
table(DZ.cells.expr$contamination) #1791 CON

DZ.cells.expr.CON__ <- filter(DZ.cells.expr, contamination == "CON" )
nrow(DZ.cells.expr.CON__) #1791
glimpse(DZ.cells.expr.CON__)


DZ.cells.expr.2 <- filter(DZ.cells.expr, contamination == "ok" )
nrow(DZ.cells.expr.2) #356920
glimpse(DZ.cells.expr.2) 


### get unique PG.ProteinAccessions => needed for long format to wide format transformation step
unique.PG.ProteinGroups.after.quality.filters.V2 <- unique(DZ.cells.expr.2$PG.ProteinGroups)
glimpse(unique.PG.ProteinGroups.after.quality.filters.V2) #6160


#### long format to wide format spread
DZ.cells.expr.2.PG.wide <- spread(DZ.cells.expr.2, key=Condition ,  value = PG.Quantity)
glimpse(DZ.cells.expr.2.PG.wide)
length(unique(DZ.cells.expr.2.PG.wide$PG.ProteinGroups)) #6160


### change sample names 

# Proteomics Sample Label #Investigator's Sample Label	#Source Of Material
# A	1	GFP
# F	2	GFP
# K	3	GFP
# P	4	GFP
# B	5	Fgfr2-FL
# G	6	Fgfr2-FL
# L	7	Fgfr2-FL
# Q	8	Fgfr2-FL
# C	9	Fgfr2-dE18.
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


colnames(DZ.cells.expr.2.PG.wide)
colnames(DZ.cells.expr.2.PG.wide)<-c("PG.Qvalue",                                       
                                     "PG.IsSingleHit",                                  
                                     "PG.FastaHeaders",                                
                                     "PG.Genes",                                        
                                     "PG.Organisms",  
                                     "PG.ProteinAccessions",      
                                     "PG.FastaFiles",                                  
                                     "PG.ProteinDescriptions",                          
                                     "PG.NrOfStrippedSequencesMeasured",                
                                     "PG.NrOfModifiedSequencesMeasured",               
                                     "PG.NrOfPrecursorsMeasured",                       
                                     "PG.ProteinGroups",         
                                     "PG.CellularComponent",                            
                                     "PG.BiologicalProcess",                           
                                     "PG.MolecularFunction",                            
                                     "PG.ProteinNames",                                 
                                     "PG.UniProtIds",                                  
                                     "contamination",                                   
                                     "GFP.1_A",                 "GFP.1_B", 
                                     "GFP.1_C",  "Fgfr2-dE18.10_A",                "Fgfr2-dE18.10_B",
                                     "Fgfr2-dE18.10_C", "Fgfr2-dE18.11_A",                "Fgfr2-dE18.11_B",
                                     "Fgfr2-dE18.11_C", "Fgfr2-dE18.12_A",                "Fgfr2-dE18.12_B",
                                     "Fgfr2-dE18.12_C", "Fgfr2-FL-Bicc1.13_A",                "Fgfr2-FL-Bicc1.13_B",
                                     "Fgfr2-FL-Bicc1.13_C", "Fgfr2-FL-Bicc1.14_A",                "Fgfr2-FL-Bicc1.14_B",
                                     "Fgfr2-FL-Bicc1.14_C", "Fgfr2-FL-Bicc1.15_A",                "Fgfr2-FL-Bicc1.15_B",
                                     "Fgfr2-FL-Bicc1.15_C", "Fgfr2-FL-Bicc1.16_A",                "Fgfr2-FL-Bicc1.16_B",
                                     "Fgfr2-FL-Bicc1.16_C", "Fgfr2-dE18-Bicc1.17_A",                "Fgfr2-dE18-Bicc1.17_B",
                                     "Fgfr2-dE18-Bicc1.17_C", "Fgfr2-dE18-Bicc1.18_A",                "Fgfr2-dE18-Bicc1.18_B",
                                     "Fgfr2-dE18-Bicc1.18_C", "Fgfr2-dE18-Bicc1.19_A",                "Fgfr2-dE18-Bicc1.19_B",
                                     "Fgfr2-dE18-Bicc1.19_C", "GFP.2_A",                 "GFP.2_B", 
                                     "GFP.2_C",  "Fgfr2-dE18-Bicc1.20_A",                "Fgfr2-dE18-Bicc1.20_B",
                                     "Fgfr2-dE18-Bicc1.20_C", "GFP.3_A",                 "GFP.3_B", 
                                     "GFP.3_C",  "GFP.4_A",                 "GFP.4_B", 
                                     "GFP.4_C",  "Fgfr2-FL.5_A",                 "Fgfr2-FL.5_B", 
                                     "Fgfr2-FL.5_C",  "Fgfr2-FL.6_A",                 "Fgfr2-FL.6_B", 
                                     "Fgfr2-FL.6_C",  "Fgfr2-FL.7_A",                 "Fgfr2-FL.7_B", 
                                     "Fgfr2-FL.7_C",  "Fgfr2-FL.8_A",                 "Fgfr2-FL.8_B", 
                                     "Fgfr2-FL.8_C",  "Fgfr2-dE18.9_A",                 "Fgfr2-dE18.9_B", 
                                     "Fgfr2-dE18.9_C" )

glimpse(DZ.cells.expr.2.PG.wide)

### reorder samples
DZ.cells.expr.2.PG.wide <- DZ.cells.expr.2.PG.wide %>% select(PG.Qvalue :contamination, 
                                                              contains("GFP"), 
                                                              "Fgfr2-FL.5_A", "Fgfr2-FL.5_B", "Fgfr2-FL.5_C",
                                                              "Fgfr2-FL.6_A", "Fgfr2-FL.6_B", "Fgfr2-FL.6_C",
                                                              "Fgfr2-FL.7_A", "Fgfr2-FL.7_B", "Fgfr2-FL.7_C",
                                                              "Fgfr2-FL.8_A", "Fgfr2-FL.8_B", "Fgfr2-FL.8_C",
                                                              
                                                              "Fgfr2-dE18.9_A","Fgfr2-dE18.9_B","Fgfr2-dE18.9_C",
                                                              "Fgfr2-dE18.10_A","Fgfr2-dE18.10_B","Fgfr2-dE18.10_C",
                                                              "Fgfr2-dE18.11_A","Fgfr2-dE18.11_B","Fgfr2-dE18.11_C",
                                                              "Fgfr2-dE18.12_A","Fgfr2-dE18.12_B","Fgfr2-dE18.12_C",
                                                              
                                                              "Fgfr2-FL-Bicc1.13_A","Fgfr2-FL-Bicc1.13_B","Fgfr2-FL-Bicc1.13_C",
                                                              "Fgfr2-FL-Bicc1.14_A","Fgfr2-FL-Bicc1.14_B","Fgfr2-FL-Bicc1.14_C",
                                                              "Fgfr2-FL-Bicc1.15_A","Fgfr2-FL-Bicc1.15_B","Fgfr2-FL-Bicc1.15_C",
                                                              "Fgfr2-FL-Bicc1.16_A","Fgfr2-FL-Bicc1.16_B","Fgfr2-FL-Bicc1.16_C",
                                                              
                                                              "Fgfr2-dE18-Bicc1.17_A","Fgfr2-dE18-Bicc1.17_B","Fgfr2-dE18-Bicc1.17_C",
                                                              "Fgfr2-dE18-Bicc1.18_A","Fgfr2-dE18-Bicc1.18_B","Fgfr2-dE18-Bicc1.18_C",
                                                              "Fgfr2-dE18-Bicc1.19_A","Fgfr2-dE18-Bicc1.19_B","Fgfr2-dE18-Bicc1.19_C",
                                                              "Fgfr2-dE18-Bicc1.20_A","Fgfr2-dE18-Bicc1.20_B","Fgfr2-dE18-Bicc1.20_C")

glimpse(DZ.cells.expr.2.PG.wide)

### boxplot and density plot samples to check normalization #####################################################################################

#get data
data.for.boxplot <- DZ.cells.expr.2.PG.wide %>% select(GFP.1_A : `Fgfr2-dE18-Bicc1.20_A` )
data.for.boxplot <- log2(data.for.boxplot)
glimpse(data.for.boxplot)

# reshape for plot
melted.data.for.boxplot <- reshape2::melt(data.for.boxplot)
glimpse(melted.data.for.boxplot)
unique(melted.data.for.boxplot$variable)

#box plot
ggplot(data = melted.data.for.boxplot)+
  geom_boxplot(mapping = aes(x= variable, y = value, fill=variable), outlier.size = 0.5, notch = TRUE, notchwidth = 0.1)+
  scale_fill_manual(values = c(rep("yellow", 12), rep("blue", 12), rep("red", 12), rep("black", 12), rep("grey", 12) ) )+ 
  theme(axis.text.x = element_text(angle = 90, hjust = 1.0, vjust = 0.5))+
  theme(legend.position="none", legend.justification = "center")+
  ggtitle("cells expression DIA") +
  theme(plot.title = element_text(hjust = 0.5, face= "bold"))+
  theme(axis.text.x = element_text(face= "plain", colour="black", size=14))+
  theme(axis.text.y = element_text(face= "plain", colour="black", size=14))+
  xlab(NULL)+ 
  ylab("log2(Int.)") 

# density plot
glimpse(melted.data.for.boxplot)
ggplot(melted.data.for.boxplot, aes(x=value, colour = variable) ) + 
  geom_density() +
  scale_color_manual(values = c(rep("yellow", 12), rep("blue", 12), rep("red", 12), rep("black", 12), rep("grey", 12) ))+
  theme(legend.position="bottom")+
  ggtitle("cells 1 DZ expression DIA \n data distribution") +
  theme(plot.title = element_text(hjust = 0.5, face= "bold"))+  #center title
  theme(axis.text.x = element_text(face= "plain", colour="black", size=14))+
  theme(axis.text.y = element_text(face= "plain", colour="black", size=14))


### count indentifications per sample ###################################################################################################################################
glimpse(DZ.cells.expr.2.PG.wide)
colnames(DZ.cells.expr.2.PG.wide)

result.EXPR.count.LFQ <- c()
sample.name.count <- 1

for(i in 19:78 ){ #23:82 
  print(i)
  #i=22; i
  
  sample.name <- colnames(DZ.cells.expr.2.PG.wide)
  sample.name <- sample.name[i]
  sample.name
  print(sample.name)
  
  temp <- DZ.cells.expr.2.PG.wide %>% select(all_of(sample.name)) %>% pull()
  temp
  temp.data <- sum(!is.na(temp))
  temp.na<- sum(is.na(temp))
  glimpse(temp)
  temp.data
  
  temp.result <- tibble(sample=sample.name, identifications=temp.data, MVs=temp.na); temp.result
  
  result.EXPR.count.LFQ  <- bind_rows(result.EXPR.count.LFQ, temp.result)
  
  
}
print(result.EXPR.count.LFQ, n=60)
glimpse(result.EXPR.count.LFQ)

mean(result.EXPR.count.LFQ$identifications) #5948.667


### plot numbers expression

ggplot()+
  geom_col(data= result.EXPR.count.LFQ, aes(x=fct_inorder(sample) , y=identifications), color="black", fill="blue")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1.0, vjust = 0.5))+
  ylab("identificatins")+
  xlab(NULL)+
  #scale_y_continuous(breaks=seq(from=0, to=4000, by=250), limits = c(0,4000))+
  ggtitle("cells expr. DIA - identifications") +
  theme(plot.title = element_text(hjust = 0.5, face= "bold"))+  #center title
  theme(axis.text.x = element_text(face= "plain", colour="black", size=14))+
  theme(axis.text.y = element_text(face= "plain", colour="black", size=14))






#############################################################################################################################################################
### correlation EXPRESSION proteins with averaged technical replicates

glimpse(DZ.cells.expr.2.PG.wide)


### select columns and order and log2 transform
DF_for_correlation_plot.EXPR.Prot <- DZ.cells.expr.2.PG.wide %>% select(GFP.1_A :  `Fgfr2-dE18-Bicc1.20_C`  )
glimpse(DF_for_correlation_plot.EXPR.Prot)

DF_for_correlation_plot.EXPR.Prot$GFP.1 <- pbapply(DF_for_correlation_plot.EXPR.Prot[, c("GFP.1_A", "GFP.1_B", "GFP.1_C")], 1, function(x) mean(x))
DF_for_correlation_plot.EXPR.Prot$GFP.2 <- pbapply(DF_for_correlation_plot.EXPR.Prot[, c("GFP.2_A", "GFP.2_B", "GFP.2_C")], 1, function(x) mean(x))
DF_for_correlation_plot.EXPR.Prot$GFP.3 <- pbapply(DF_for_correlation_plot.EXPR.Prot[, c("GFP.3_A", "GFP.3_B", "GFP.3_C")], 1, function(x) mean(x))
DF_for_correlation_plot.EXPR.Prot$GFP.4 <- pbapply(DF_for_correlation_plot.EXPR.Prot[, c("GFP.4_A", "GFP.4_B", "GFP.4_C")], 1, function(x) mean(x))

DF_for_correlation_plot.EXPR.Prot$FL.1 <- pbapply(DF_for_correlation_plot.EXPR.Prot[, c("Fgfr2-FL.5_A", "Fgfr2-FL.5_B", "Fgfr2-FL.5_C")], 1, function(x) mean(x))
DF_for_correlation_plot.EXPR.Prot$FL.2 <- pbapply(DF_for_correlation_plot.EXPR.Prot[, c("Fgfr2-FL.6_A", "Fgfr2-FL.6_B", "Fgfr2-FL.6_C")], 1, function(x) mean(x))
DF_for_correlation_plot.EXPR.Prot$FL.3 <- pbapply(DF_for_correlation_plot.EXPR.Prot[, c("Fgfr2-FL.7_A", "Fgfr2-FL.7_B", "Fgfr2-FL.7_C")], 1, function(x) mean(x))
DF_for_correlation_plot.EXPR.Prot$FL.4 <- pbapply(DF_for_correlation_plot.EXPR.Prot[, c("Fgfr2-FL.8_A", "Fgfr2-FL.8_B", "Fgfr2-FL.8_C")], 1, function(x) mean(x))

DF_for_correlation_plot.EXPR.Prot$dE18.1 <- pbapply(DF_for_correlation_plot.EXPR.Prot[, c("Fgfr2-dE18.9_A", "Fgfr2-dE18.9_B", "Fgfr2-dE18.9_C")], 1, function(x) mean(x))
DF_for_correlation_plot.EXPR.Prot$dE18.2 <- pbapply(DF_for_correlation_plot.EXPR.Prot[, c("Fgfr2-dE18.10_A", "Fgfr2-dE18.10_B", "Fgfr2-dE18.10_C")], 1, function(x) mean(x))
DF_for_correlation_plot.EXPR.Prot$dE18.3 <- pbapply(DF_for_correlation_plot.EXPR.Prot[, c("Fgfr2-dE18.11_A", "Fgfr2-dE18.11_B", "Fgfr2-dE18.11_C")], 1, function(x) mean(x))
DF_for_correlation_plot.EXPR.Prot$dE18.4 <- pbapply(DF_for_correlation_plot.EXPR.Prot[, c("Fgfr2-dE18.12_A", "Fgfr2-dE18.12_B", "Fgfr2-dE18.12_C")], 1, function(x) mean(x))

DF_for_correlation_plot.EXPR.Prot$FL.Bicc1.1 <- pbapply(DF_for_correlation_plot.EXPR.Prot[, c("Fgfr2-FL-Bicc1.13_A", "Fgfr2-FL-Bicc1.13_B", "Fgfr2-FL-Bicc1.13_C")], 1, function(x) mean(x))
DF_for_correlation_plot.EXPR.Prot$FL.Bicc1.2 <- pbapply(DF_for_correlation_plot.EXPR.Prot[, c("Fgfr2-FL-Bicc1.14_A", "Fgfr2-FL-Bicc1.14_B", "Fgfr2-FL-Bicc1.14_C")], 1, function(x) mean(x))
DF_for_correlation_plot.EXPR.Prot$FL.Bicc1.3 <- pbapply(DF_for_correlation_plot.EXPR.Prot[, c("Fgfr2-FL-Bicc1.15_A", "Fgfr2-FL-Bicc1.15_B", "Fgfr2-FL-Bicc1.15_C")], 1, function(x) mean(x))
DF_for_correlation_plot.EXPR.Prot$FL.Bicc1.4 <- pbapply(DF_for_correlation_plot.EXPR.Prot[, c("Fgfr2-FL-Bicc1.16_A", "Fgfr2-FL-Bicc1.16_B", "Fgfr2-FL-Bicc1.16_C")], 1, function(x) mean(x))

DF_for_correlation_plot.EXPR.Prot$dE18.Bicc1.1 <- pbapply(DF_for_correlation_plot.EXPR.Prot[, c("Fgfr2-dE18-Bicc1.17_A", "Fgfr2-dE18-Bicc1.17_B", "Fgfr2-dE18-Bicc1.17_C")], 1, function(x) mean(x))
DF_for_correlation_plot.EXPR.Prot$dE18.Bicc1.2 <- pbapply(DF_for_correlation_plot.EXPR.Prot[, c("Fgfr2-dE18-Bicc1.18_A", "Fgfr2-dE18-Bicc1.18_B", "Fgfr2-dE18-Bicc1.18_C")], 1, function(x) mean(x))
DF_for_correlation_plot.EXPR.Prot$dE18.Bicc1.3 <- pbapply(DF_for_correlation_plot.EXPR.Prot[, c("Fgfr2-dE18-Bicc1.19_A", "Fgfr2-dE18-Bicc1.19_B", "Fgfr2-dE18-Bicc1.19_C")], 1, function(x) mean(x))
DF_for_correlation_plot.EXPR.Prot$dE18.Bicc1.4 <- pbapply(DF_for_correlation_plot.EXPR.Prot[, c("Fgfr2-dE18-Bicc1.20_A", "Fgfr2-dE18-Bicc1.20_B", "Fgfr2-dE18-Bicc1.20_C")], 1, function(x) mean(x))

glimpse(DF_for_correlation_plot.EXPR.Prot)

DF_for_correlation_plot.EXPR.Prot <- DF_for_correlation_plot.EXPR.Prot %>% select("GFP.1","GFP.2","GFP.3","GFP.4",  "FL.1","FL.2","FL.3","FL.4",  "dE18.1","dE18.2","dE18.3","dE18.4",  "FL.Bicc1.1","FL.Bicc1.2","FL.Bicc1.3","FL.Bicc1.4",   "dE18.Bicc1.1","dE18.Bicc1.2","dE18.Bicc1.3","dE18.Bicc1.4")
glimpse(DF_for_correlation_plot.EXPR.Prot)


DF_for_correlation_plot.EXPR.Prot <- log2(DF_for_correlation_plot.EXPR.Prot)
glimpse(DF_for_correlation_plot.EXPR.Prot)



### check NA rows 
DF_for_correlation_plot.EXPR.Prot$na.rows <- pbapply(DF_for_correlation_plot.EXPR.Prot, 1, function(x) sum(is.na(x)))
glimpse(DF_for_correlation_plot.EXPR.Prot)
table(DF_for_correlation_plot.EXPR.Prot$na.rows)

DF_for_correlation_plot.EXPR.Prot <- DF_for_correlation_plot.EXPR.Prot %>% select(-na.rows )
glimpse(DF_for_correlation_plot.EXPR.Prot)
colnames(DF_for_correlation_plot.EXPR.Prot)


### correlation
corr2 <- cor(DF_for_correlation_plot.EXPR.Prot, method = "pearson", use = "na.or.complete") 
glimpse(corr2)
min(corr2) #0.8786569

### prepare heatmap annotation
annot_df_for_heatmap <- data.frame(samples = colnames(DF_for_correlation_plot.EXPR.Prot),
                                   group = c(rep("GFP", 4), rep("FL", 4), rep("dE18", 4), rep("FL-Bicc1", 4), rep("dE18-Bicc1", 4)) #
                                   #group = c(                rep("FL", 4), rep("dE18", 4), rep("FL-Bicc1", 4), rep("dE18-Bicc1", 3)) #
                                   
)

annot_df_for_heatmap 


### make shorter annotation data frame
annot_df_for_heatmap.short <- data.frame( 
  group  = annot_df_for_heatmap$group
)

glimpse(annot_df_for_heatmap.short)
annot_df_for_heatmap.short

### define colors for annotation bar
annot.colors = list(group = c("GFP"="yellow", "FL"="blue", "dE18"="red", "FL-Bicc1"="black", "dE18-Bicc1"="grey40")
)
annot.colors


### Create the heatmap annotation
ha.fr <- HeatmapAnnotation(group=annot_df_for_heatmap.short$group,
                           treatment=annot_df_for_heatmap.short$treatment,
                           col = annot.colors,
                           annotation_legend_param = list(grid_height = unit(8, "mm"))
)
ha.fr

#prepare heatmap colors
FR.heatmap.colors.2 <-colorRamp2(c(min(corr2), 1.0), c("grey90", "grey10"))


#pdf("DZ.cells.expr.DIA.correlation.heatmap.averaged.pdf", width=19/2.54, height=17/2.54, useDingbats=FALSE)
Heatmap(corr2, 
        name = "corr.coeff", #title of legend
        
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

##################################################################################################################################
### multiple candidates plot to check expression
### average technical replicates

mtpl.candidates <- c("FGFR2_MOUSE",
                     #"BICC1_MOUSE",
                     "PAK1_MOUSE",  "PAXI_MOUSE",  "MP2K1_MOUSE", "MP2K2_MOUSE", "MK01_MOUSE",  "MK03_MOUSE",  "ERF_MOUSE",   "JUN_MOUSE",   "JUNB_MOUSE",  "GSK3B_MOUSE", "MAP1B_MOUSE",
                     "AKTS1_MOUSE", "ACLY_MOUSE",  "BAD_MOUSE",   "TBCD1_MOUSE", "FOXO3_MOUSE", 
                     "KS6A3_MOUSE", "RS6_MOUSE",   "IF4B_MOUSE",  "4EBP1_MOUSE", "PYR1_MOUSE",  "LARP1_MOUSE"
                     )

EXPR.mtplplot.data <- DZ.cells.expr.2.PG.wide %>% filter(PG.ProteinNames %in% mtpl.candidates)
glimpse(EXPR.mtplplot.data)



### Reorder data according to order candidates, average technical replicates, log2 transform, reshape for plot
EXPR.mtplplot.data <- left_join(data.frame(PG.ProteinNames = mtpl.candidates),
                                EXPR.mtplplot.data,
                                by = "PG.ProteinNames")

EXPR.mtplplot.data.2A <- EXPR.mtplplot.data  %>% select(PG.Genes, PG.ProteinNames,PG.UniProtIds) 
glimpse(EXPR.mtplplot.data.2A)

EXPR.mtplplot.data.2B <- EXPR.mtplplot.data  %>% select(
  contains("GFP"), 
  "Fgfr2-FL.5_A", "Fgfr2-FL.5_B", "Fgfr2-FL.5_C",
  "Fgfr2-FL.6_A", "Fgfr2-FL.6_B", "Fgfr2-FL.6_C",
  "Fgfr2-FL.7_A", "Fgfr2-FL.7_B", "Fgfr2-FL.7_C",
  "Fgfr2-FL.8_A", "Fgfr2-FL.8_B", "Fgfr2-FL.8_C",
  
  "Fgfr2-dE18.9_A","Fgfr2-dE18.9_B","Fgfr2-dE18.9_C",
  "Fgfr2-dE18.10_A","Fgfr2-dE18.10_B","Fgfr2-dE18.10_C",
  "Fgfr2-dE18.11_A","Fgfr2-dE18.11_B","Fgfr2-dE18.11_C",
  "Fgfr2-dE18.12_A","Fgfr2-dE18.12_B","Fgfr2-dE18.12_C",
  
  "Fgfr2-FL-Bicc1.13_A","Fgfr2-FL-Bicc1.13_B","Fgfr2-FL-Bicc1.13_C",
  "Fgfr2-FL-Bicc1.14_A","Fgfr2-FL-Bicc1.14_B","Fgfr2-FL-Bicc1.14_C",
  "Fgfr2-FL-Bicc1.15_A","Fgfr2-FL-Bicc1.15_B","Fgfr2-FL-Bicc1.15_C",
  "Fgfr2-FL-Bicc1.16_A","Fgfr2-FL-Bicc1.16_B","Fgfr2-FL-Bicc1.16_C",
  
  "Fgfr2-dE18-Bicc1.17_A","Fgfr2-dE18-Bicc1.17_B","Fgfr2-dE18-Bicc1.17_C",
  "Fgfr2-dE18-Bicc1.18_A","Fgfr2-dE18-Bicc1.18_B","Fgfr2-dE18-Bicc1.18_C",
  "Fgfr2-dE18-Bicc1.19_A","Fgfr2-dE18-Bicc1.19_B","Fgfr2-dE18-Bicc1.19_C",
  "Fgfr2-dE18-Bicc1.20_A","Fgfr2-dE18-Bicc1.20_B","Fgfr2-dE18-Bicc1.20_C") 
glimpse(EXPR.mtplplot.data.2B)

# average technical replicates
EXPR.mtplplot.data.2C <- EXPR.mtplplot.data.2B

EXPR.mtplplot.data.2C$GFP.1 <- pbapply(EXPR.mtplplot.data.2C[, c("GFP.1_A", "GFP.1_B", "GFP.1_C")], 1, function(x) mean(x))
EXPR.mtplplot.data.2C$GFP.2 <- pbapply(EXPR.mtplplot.data.2C[, c("GFP.2_A", "GFP.2_B", "GFP.2_C")], 1, function(x) mean(x))
EXPR.mtplplot.data.2C$GFP.3 <- pbapply(EXPR.mtplplot.data.2C[, c("GFP.3_A", "GFP.3_B", "GFP.3_C")], 1, function(x) mean(x))
EXPR.mtplplot.data.2C$GFP.4 <- pbapply(EXPR.mtplplot.data.2C[, c("GFP.4_A", "GFP.4_B", "GFP.4_C")], 1, function(x) mean(x))

EXPR.mtplplot.data.2C$FL.1 <- pbapply(EXPR.mtplplot.data.2C[, c("Fgfr2-FL.5_A", "Fgfr2-FL.5_B", "Fgfr2-FL.5_C")], 1, function(x) mean(x))
EXPR.mtplplot.data.2C$FL.2 <- pbapply(EXPR.mtplplot.data.2C[, c("Fgfr2-FL.6_A", "Fgfr2-FL.6_B", "Fgfr2-FL.6_C")], 1, function(x) mean(x))
EXPR.mtplplot.data.2C$FL.3 <- pbapply(EXPR.mtplplot.data.2C[, c("Fgfr2-FL.7_A", "Fgfr2-FL.7_B", "Fgfr2-FL.7_C")], 1, function(x) mean(x))
EXPR.mtplplot.data.2C$FL.4 <- pbapply(EXPR.mtplplot.data.2C[, c("Fgfr2-FL.8_A", "Fgfr2-FL.8_B", "Fgfr2-FL.8_C")], 1, function(x) mean(x))

EXPR.mtplplot.data.2C$dE18.1 <- pbapply(EXPR.mtplplot.data.2C[, c("Fgfr2-dE18.9_A", "Fgfr2-dE18.9_B", "Fgfr2-dE18.9_C")], 1, function(x) mean(x))
EXPR.mtplplot.data.2C$dE18.2 <- pbapply(EXPR.mtplplot.data.2C[, c("Fgfr2-dE18.10_A", "Fgfr2-dE18.10_B", "Fgfr2-dE18.10_C")], 1, function(x) mean(x))
EXPR.mtplplot.data.2C$dE18.3 <- pbapply(EXPR.mtplplot.data.2C[, c("Fgfr2-dE18.11_A", "Fgfr2-dE18.11_B", "Fgfr2-dE18.11_C")], 1, function(x) mean(x))
EXPR.mtplplot.data.2C$dE18.4 <- pbapply(EXPR.mtplplot.data.2C[, c("Fgfr2-dE18.12_A", "Fgfr2-dE18.12_B", "Fgfr2-dE18.12_C")], 1, function(x) mean(x))

EXPR.mtplplot.data.2C$FL.Bicc1.1 <- pbapply(EXPR.mtplplot.data.2C[, c("Fgfr2-FL-Bicc1.13_A", "Fgfr2-FL-Bicc1.13_B", "Fgfr2-FL-Bicc1.13_C")], 1, function(x) mean(x))
EXPR.mtplplot.data.2C$FL.Bicc1.2 <- pbapply(EXPR.mtplplot.data.2C[, c("Fgfr2-FL-Bicc1.14_A", "Fgfr2-FL-Bicc1.14_B", "Fgfr2-FL-Bicc1.14_C")], 1, function(x) mean(x))
EXPR.mtplplot.data.2C$FL.Bicc1.3 <- pbapply(EXPR.mtplplot.data.2C[, c("Fgfr2-FL-Bicc1.15_A", "Fgfr2-FL-Bicc1.15_B", "Fgfr2-FL-Bicc1.15_C")], 1, function(x) mean(x))
EXPR.mtplplot.data.2C$FL.Bicc1.4 <- pbapply(EXPR.mtplplot.data.2C[, c("Fgfr2-FL-Bicc1.16_A", "Fgfr2-FL-Bicc1.16_B", "Fgfr2-FL-Bicc1.16_C")], 1, function(x) mean(x))

EXPR.mtplplot.data.2C$dE18.Bicc1.1 <- pbapply(EXPR.mtplplot.data.2C[, c("Fgfr2-dE18-Bicc1.17_A", "Fgfr2-dE18-Bicc1.17_B", "Fgfr2-dE18-Bicc1.17_C")], 1, function(x) mean(x))
EXPR.mtplplot.data.2C$dE18.Bicc1.2 <- pbapply(EXPR.mtplplot.data.2C[, c("Fgfr2-dE18-Bicc1.18_A", "Fgfr2-dE18-Bicc1.18_B", "Fgfr2-dE18-Bicc1.18_C")], 1, function(x) mean(x))
EXPR.mtplplot.data.2C$dE18.Bicc1.3 <- pbapply(EXPR.mtplplot.data.2C[, c("Fgfr2-dE18-Bicc1.19_A", "Fgfr2-dE18-Bicc1.19_B", "Fgfr2-dE18-Bicc1.19_C")], 1, function(x) mean(x))
EXPR.mtplplot.data.2C$dE18.Bicc1.4 <- pbapply(EXPR.mtplplot.data.2C[, c("Fgfr2-dE18-Bicc1.20_A", "Fgfr2-dE18-Bicc1.20_B", "Fgfr2-dE18-Bicc1.20_C")], 1, function(x) mean(x))

EXPR.mtplplot.data.2C <- EXPR.mtplplot.data.2C %>% select("GFP.1","GFP.2","GFP.3","GFP.4",  "FL.1","FL.2","FL.3","FL.4",  "dE18.1","dE18.2","dE18.3","dE18.4",  "FL.Bicc1.1","FL.Bicc1.2","FL.Bicc1.3","FL.Bicc1.4",   "dE18.Bicc1.1","dE18.Bicc1.2","dE18.Bicc1.3","dE18.Bicc1.4")
glimpse(EXPR.mtplplot.data.2C)

# log2 data transformation
EXPR.mtplplot.data.2D <- log2(EXPR.mtplplot.data.2C)
glimpse(EXPR.mtplplot.data.2D)

# combine data frames
EXPR.mtplplot.data.3 <- bind_cols(EXPR.mtplplot.data.2A, EXPR.mtplplot.data.2D) #
glimpse(EXPR.mtplplot.data.3)
nrow(EXPR.mtplplot.data.3)

# prepare additional candidate names
EXPR.mtplplot.data.3$ProteinNames_Genes_UniprotID <- paste0(EXPR.mtplplot.data.3$PG.ProteinNames, "_", EXPR.mtplplot.data.3$PG.Genes, "_", EXPR.mtplplot.data.3$PG.UniProtIds)
glimpse(EXPR.mtplplot.data.3)

EXPR.mtplplot.data.3$PG.UniProtIds_2 <- pbsapply(EXPR.mtplplot.data.3$PG.UniProtIds, function(x)  unlist(str_split(string=x, pattern=";"))[1])
EXPR.mtplplot.data.3$PG.Genes_2 <- paste0(EXPR.mtplplot.data.3$PG.Genes, "_", EXPR.mtplplot.data.3$PG.UniProtIds_2)
glimpse(EXPR.mtplplot.data.3)

# reshape
melted.EXPR.mtplplot.data.3 <- reshape2::melt(EXPR.mtplplot.data.3)
glimpse(melted.EXPR.mtplplot.data.3)

# define order in plot
order.in.plot2 <- EXPR.mtplplot.data.3$PG.Genes 
order.in.plot2



### ggplot
ggplot( melted.EXPR.mtplplot.data.3 , aes(x=variable, y=PG.Genes, fill=value)) +  
  geom_tile() + 
  scale_fill_viridis( option = "inferno", na.value = "grey30",name="log2(intensity)") + #
  coord_equal()+
  scale_y_discrete(limits=  rev(order.in.plot2))+
  
  theme(legend.position="bottom", legend.justification = "center")+
  
  theme(axis.line = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  
  theme(legend.title.align=0.5)+
  
  #ggtitle("cells expression candidates") +
  
  theme(plot.title = element_text(hjust = 0.5, face="bold", size=16)) +
  theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1))+ 
  theme(axis.text.y= element_text(size=9))+  
  xlab(NULL) + 
  xlab(NULL) + 
  ylab(NULL) +
  
  theme(axis.text.x = element_blank()) +
  theme(axis.ticks.x = element_blank())+
  
  geom_vline(xintercept = c(4.5, 8.5, 12.5, 16.5) , size=1.0, linetype = "solid", color="white")+
  annotate("text", x = 2.5, y = nrow(EXPR.mtplplot.data.3)+2, label = "GFP", fontface = "bold",size=4.0)+
  annotate("text", x = 6.5, y = nrow(EXPR.mtplplot.data.3)+2, label = "FL", fontface = "bold",size=4.0)+
  annotate("text", x = 10.5, y = nrow(EXPR.mtplplot.data.3)+2, label = "dE18", fontface = "bold",size=4.0)+
  annotate("text", x = 14.5, y = nrow(EXPR.mtplplot.data.3)+2, label = "FL\nBicc1", fontface = "bold",size=4.0)+
  annotate("text", x = 18.5, y = nrow(EXPR.mtplplot.data.3)+2, label = "dE18\nBicc1", fontface = "bold" ,size=4.0)+
  annotate("text", x = 22.5, y = nrow(EXPR.mtplplot.data.3)+3.0, label = "",  color = "transparent")
 
#ggsave("DZ.cells.expr.DIA.candidates.averaged.pdf", useDingbats=FALSE,  width = 14, height =14, units = "cm")



#################################################################################################################################
### ssGSEA

#get data and reshape for ssGSEA
glimpse(DZ.cells.expr.2.PG.wide)

input.ssGSEA <- DZ.cells.expr.2.PG.wide
glimpse(input.ssGSEA)

input.ssGSEA$GFP.1 <- pbapply(input.ssGSEA[, c("GFP.1_A", "GFP.1_B", "GFP.1_C")], 1, function(x) mean(x))
input.ssGSEA$GFP.2 <- pbapply(input.ssGSEA[, c("GFP.2_A", "GFP.2_B", "GFP.2_C")], 1, function(x) mean(x))
input.ssGSEA$GFP.3 <- pbapply(input.ssGSEA[, c("GFP.3_A", "GFP.3_B", "GFP.3_C")], 1, function(x) mean(x))
input.ssGSEA$GFP.4 <- pbapply(input.ssGSEA[, c("GFP.4_A", "GFP.4_B", "GFP.4_C")], 1, function(x) mean(x))

input.ssGSEA$FL.1 <- pbapply(input.ssGSEA[, c("Fgfr2-FL.5_A", "Fgfr2-FL.5_B", "Fgfr2-FL.5_C")], 1, function(x) mean(x))
input.ssGSEA$FL.2 <- pbapply(input.ssGSEA[, c("Fgfr2-FL.6_A", "Fgfr2-FL.6_B", "Fgfr2-FL.6_C")], 1, function(x) mean(x))
input.ssGSEA$FL.3 <- pbapply(input.ssGSEA[, c("Fgfr2-FL.7_A", "Fgfr2-FL.7_B", "Fgfr2-FL.7_C")], 1, function(x) mean(x))
input.ssGSEA$FL.4 <- pbapply(input.ssGSEA[, c("Fgfr2-FL.8_A", "Fgfr2-FL.8_B", "Fgfr2-FL.8_C")], 1, function(x) mean(x))

input.ssGSEA$dE18.1 <- pbapply(input.ssGSEA[, c("Fgfr2-dE18.9_A", "Fgfr2-dE18.9_B", "Fgfr2-dE18.9_C")], 1, function(x) mean(x))
input.ssGSEA$dE18.2 <- pbapply(input.ssGSEA[, c("Fgfr2-dE18.10_A", "Fgfr2-dE18.10_B", "Fgfr2-dE18.10_C")], 1, function(x) mean(x))
input.ssGSEA$dE18.3 <- pbapply(input.ssGSEA[, c("Fgfr2-dE18.11_A", "Fgfr2-dE18.11_B", "Fgfr2-dE18.11_C")], 1, function(x) mean(x))
input.ssGSEA$dE18.4 <- pbapply(input.ssGSEA[, c("Fgfr2-dE18.12_A", "Fgfr2-dE18.12_B", "Fgfr2-dE18.12_C")], 1, function(x) mean(x))

input.ssGSEA$FL.Bicc1.1 <- pbapply(input.ssGSEA[, c("Fgfr2-FL-Bicc1.13_A", "Fgfr2-FL-Bicc1.13_B", "Fgfr2-FL-Bicc1.13_C")], 1, function(x) mean(x))
input.ssGSEA$FL.Bicc1.2 <- pbapply(input.ssGSEA[, c("Fgfr2-FL-Bicc1.14_A", "Fgfr2-FL-Bicc1.14_B", "Fgfr2-FL-Bicc1.14_C")], 1, function(x) mean(x))
input.ssGSEA$FL.Bicc1.3 <- pbapply(input.ssGSEA[, c("Fgfr2-FL-Bicc1.15_A", "Fgfr2-FL-Bicc1.15_B", "Fgfr2-FL-Bicc1.15_C")], 1, function(x) mean(x))
input.ssGSEA$FL.Bicc1.4 <- pbapply(input.ssGSEA[, c("Fgfr2-FL-Bicc1.16_A", "Fgfr2-FL-Bicc1.16_B", "Fgfr2-FL-Bicc1.16_C")], 1, function(x) mean(x))

input.ssGSEA$dE18.Bicc1.1 <- pbapply(input.ssGSEA[, c("Fgfr2-dE18-Bicc1.17_A", "Fgfr2-dE18-Bicc1.17_B", "Fgfr2-dE18-Bicc1.17_C")], 1, function(x) mean(x))
input.ssGSEA$dE18.Bicc1.2 <- pbapply(input.ssGSEA[, c("Fgfr2-dE18-Bicc1.18_A", "Fgfr2-dE18-Bicc1.18_B", "Fgfr2-dE18-Bicc1.18_C")], 1, function(x) mean(x))
input.ssGSEA$dE18.Bicc1.3 <- pbapply(input.ssGSEA[, c("Fgfr2-dE18-Bicc1.19_A", "Fgfr2-dE18-Bicc1.19_B", "Fgfr2-dE18-Bicc1.19_C")], 1, function(x) mean(x))
input.ssGSEA$dE18.Bicc1.4 <- pbapply(input.ssGSEA[, c("Fgfr2-dE18-Bicc1.20_A", "Fgfr2-dE18-Bicc1.20_B", "Fgfr2-dE18-Bicc1.20_C")], 1, function(x) mean(x))

glimpse(input.ssGSEA)

# pick first gene name
input.ssGSEA$PG.Genes.2 <- pbsapply(input.ssGSEA$PG.Genes, function(x)  unlist(str_split(string = x, pattern = ";"))[1])
glimpse(input.ssGSEA)

# remove possible NAs in gene names and if gene is duplicated use min rowsum entry
table(input.ssGSEA$PG.Genes.2 =="")
input.ssGSEA <- input.ssGSEA %>% filter(PG.Genes.2 !="")

#inspect duplicated observations
duplicated.observations <- input.ssGSEA  %>% count(PG.Genes.2) %>% filter(n > 1) 
glimpse(duplicated.observations)
duplicated.observations

unique.genes.in.list <- unique(input.ssGSEA$PG.Genes.2)
unique.genes.in.list


input.ssGSEA.2 <- c()
for(i in unique.genes.in.list){
  print(i)
  temp.obs <- as.data.frame(filter(input.ssGSEA, PG.Genes.2 %in% c(i)))
  temp.obs$sum.int.samples <-apply(temp.obs[,c("GFP.1","GFP.2","GFP.3","GFP.4",  "FL.1","FL.2","FL.3","FL.4",  "dE18.1","dE18.2","dE18.3","dE18.4",  "FL.Bicc1.1","FL.Bicc1.2","FL.Bicc1.3","FL.Bicc1.4",   "dE18.Bicc1.1","dE18.Bicc1.2","dE18.Bicc1.3","dE18.Bicc1.4")], 1, function(x) sum(x))
  temp.obs <- temp.obs[order(desc(temp.obs$sum.int.samples)),]
  temp.obs <- temp.obs[1,]
  input.ssGSEA.2 <- bind_rows(input.ssGSEA.2, temp.obs) 
}
glimpse(input.ssGSEA.2) #6038

#select columns of interest and log2 transform
input.ssGSEA.3a <- input.ssGSEA.2 %>% select(PG.Genes.2)
input.ssGSEA.3b <- log2(input.ssGSEA.2 %>% select("GFP.1","GFP.2","GFP.3","GFP.4",  "FL.1","FL.2","FL.3","FL.4",  "dE18.1","dE18.2","dE18.3","dE18.4",  "FL.Bicc1.1","FL.Bicc1.2","FL.Bicc1.3","FL.Bicc1.4",   "dE18.Bicc1.1","dE18.Bicc1.2","dE18.Bicc1.3","dE18.Bicc1.4"))
input.ssGSEA.3 <- bind_cols(input.ssGSEA.3a, input.ssGSEA.3b)
glimpse(input.ssGSEA.3)


# prepare gct file 
gene.pattern.ssGSEA.input_GCT <- input.ssGSEA.3[,c("GFP.1","GFP.2","GFP.3","GFP.4",  
                                                       "FL.1","FL.2","FL.3","FL.4",  
                                                       "dE18.1","dE18.2","dE18.3","dE18.4",  
                                                       "FL.Bicc1.1","FL.Bicc1.2","FL.Bicc1.3","FL.Bicc1.4",   
                                                       "dE18.Bicc1.1","dE18.Bicc1.2","dE18.Bicc1.3","dE18.Bicc1.4")]
glimpse(gene.pattern.ssGSEA.input_GCT)

gene.pattern.ssGSEA.input_GCT <- as.matrix(gene.pattern.ssGSEA.input_GCT[,c("GFP.1","GFP.2","GFP.3","GFP.4",  
                                                                            "FL.1","FL.2","FL.3","FL.4",  
                                                                            "dE18.1","dE18.2","dE18.3","dE18.4",  
                                                                            "FL.Bicc1.1","FL.Bicc1.2","FL.Bicc1.3","FL.Bicc1.4",   
                                                                            "dE18.Bicc1.1","dE18.Bicc1.2","dE18.Bicc1.3","dE18.Bicc1.4")])
  

table(is.na(gene.pattern.ssGSEA.input_GCT))
gene.pattern.ssGSEA.input_GCT[is.na(gene.pattern.ssGSEA.input_GCT)] <- 0 #impute missing values with zero


rownames(gene.pattern.ssGSEA.input_GCT) <- input.ssGSEA.3$PG.Genes.2 #add row names
glimpse(gene.pattern.ssGSEA.input_GCT)


gene.pattern.ssGSEA.input_GCT.2 <- new("GCT", mat=gene.pattern.ssGSEA.input_GCT)
glimpse(gene.pattern.ssGSEA.input_GCT.2)

#save GCT file for ssGSEA input
write_gct(gene.pattern.ssGSEA.input_GCT.2, "cells.expr.ssGSEA_input.a",ver = 2) # use file, add column "Description" in e.g. Excel and adjust column names accordingly to make file compatible => see file  "cells.expr.ssGSEA_input.b_n20x6038.gct"


## prepare gene sets database file (gmt) for Hallmarks from Msigdb via msigdbr package
packageVersion("msigdbr")
msigdbr_species() #list the species available in the msigdbr package
print(msigdbr_collections(), n=23) #show the available collections

#use Hallmark "H" gene sets Mus musculus
MsigDBofinterest <- msigdbr(species = "mouse", category = "H")
glimpse(MsigDBofinterest)
table(MsigDBofinterest$gs_cat)
table(MsigDBofinterest$gs_subcat)
unique(MsigDBofinterest$gs_name)


all.geneset.names <- unique(MsigDBofinterest$gs_name)
glimpse(all.geneset.names)

# reshape pathway definitions to list
pathway.list.with.ids <- c()
for(i in all.geneset.names){
  print(i)
  temp <- MsigDBofinterest %>% filter(gs_name %in% i)
  out <- list("pathway.name" = temp$gene_symbol)
  names(out) <- i
  pathway.list.with.ids <- append(pathway.list.with.ids, out)
}
glimpse(pathway.list.with.ids)


##write pathways gmt file see https://github.com/lwaldron/LeviRmisc/blob/master/R/writeGMT.R

# start define function
writeGMT <- function #Create a gmt (gene matrix transposed) file
### Createss a gmt (gene matrix transposed) file such as those
### provided by mSigDB or geneSigDB, from an R list object.
### Function by Levi Waldron.
(object,
 ### R list object that will be converted to GMT file.  Each element
 ### should contain a vector of gene names, and the names of the
 ### elements will used for the gene set names
 fname
 ### Output file name for .gmt file
){
  if (class(object) != "list") stop("object should be of class 'list'")
  if(file.exists(fname)) unlink(fname)
  for (iElement in 1:length(object)){
    write.table(t(c(make.names(rep(names(object)[iElement],2)),object[[iElement]])),
                sep="\t",quote=FALSE,
                file=fname,append=TRUE,col.names=FALSE,row.names=FALSE)
  }
  ### Called for the effect of writing a .gmt file
}
# stop define function


writeGMT(pathway.list.with.ids, "msigdbr_mouse_hallmarks.gmt")



###perform analysis via cloud.genepattern.org => ssGSEA module v10.0.11, standard settings



### load ssGSEA result from gene pattern and visualize

ssGSEA_result.df <- fread("cells.expr.ssGSEA_result.txt") 
glimpse(ssGSEA_result.df)


#calculate some measures
for(i in 1:nrow(ssGSEA_result.df)){
  print(i)
  
  #FL vs dE18 t.test
  tempB <- t.test(as.numeric(ssGSEA_result.df[i, c("FL.1","FL.2","FL.3","FL.4")]),
                  as.numeric(ssGSEA_result.df[i, c("dE18.1","dE18.2","dE18.3","dE18.4")]))
  ssGSEA_result.df$FL_dE18_t.test.pval[i] <- tempB$p.value
  
}
glimpse(ssGSEA_result.df)

#p.adjust
ssGSEA_result.df$FL_dE18_t.test.pval.adj.BH <- p.adjust(ssGSEA_result.df$FL_dE18_t.test.pval, "BH")

table(ssGSEA_result.df$FL_dE18_t.test.pval < 0.05)

# change column names, filter sign. observations, 
colnames(ssGSEA_result.df) <- c("pathways",                "Description",         "GFP.1",                "GFP.2" ,               "GFP.3",               
                          "GFP.4"     ,           "FL.1"   ,              "FL.2"  ,               "FL.3"   ,              "FL.4"   ,             
                          "dE18.1"   ,            "dE18.2"     ,          "dE18.3"      ,         "dE18.4"   ,            "FL.Bicc1.1"  ,        
                          "FL.Bicc1.2"   ,        "FL.Bicc1.3"        ,   "FL.Bicc1.4" ,          "dE18.Bicc1.1"  ,       "dE18.Bicc1.2"  ,      
                          "dE18.Bicc1.3"   ,      "dE18.Bicc1.4"    ,    "FL_dE18_t.test.pval", "FL_dE18_t.test.pval.adj.BH")

### ggplot heatmap
sign.ssGSEA.result <- ssGSEA_result.df %>% filter(FL_dE18_t.test.pval < 0.05 ) 
nrow(sign.ssGSEA.result) #18

custom.order.mouse.ssGSEA <- c("HALLMARK_KRAS_SIGNALING",
                               "HALLMARK_KRAS_SIGNALING_UP",
                               "HALLMARK_PI3K_AKT_MTOR_SIGNALING",
                               "HALLMARK_MTORC1_SIGNALING",
                               "HALLMARK_MYC_TARGETS_V1",
                               "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
                               "HALLMARK_TGF_BETA_SIGNALING",
                               "HALLMARK_WNT_BETA_CATENIN_SIGNALING",
                               "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
                               "HALLMARK_ANGIOGENESIS",
                               "HALLMARK_ESTROGEN_RESPONSE_EARLY",
                               "HALLMARK_ANDROGEN_RESPONSE",
                               "HALLMARK_UNFOLDED_PROTEIN_RESPONSE",
                               "HALLMARK_CHOLESTEROL_HOMEOSTASIS",
                               "HALLMARK_APICAL_SURFACE",
                               "HALLMARK_INTERFERON_GAMMA_RESPONSE",
                               "HALLMARK_INTERFERON_ALPHA_RESPONSE",
                               "HALLMARK_HEDGEHOG_SIGNALING")   
custom.order.mouse.ssGSEA

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

A <- sign.ssGSEA.result %>% select(GFP.1    : dE18.Bicc1.4 )
A <- t(apply(A,1, function(x) zscorescalebaseR(x = x)))

sd(A[1,]) #should be 1

B <- as_tibble(A)

sign.ssGSEA.result <- bind_cols(sign.ssGSEA.result %>% select(pathways), B) #with z score
glimpse(sign.ssGSEA.result)


#reshape for plot
melted.ssGSEA.result <- reshape2::melt(sign.ssGSEA.result %>% select("pathways",  
                                                                     "GFP.1","GFP.2","GFP.3","GFP.4",  
                                                                     "FL.1","FL.2","FL.3","FL.4",  
                                                                     "dE18.1","dE18.2","dE18.3","dE18.4",  
                                                                     "FL.Bicc1.1","FL.Bicc1.2","FL.Bicc1.3","FL.Bicc1.4",   
                                                                     "dE18.Bicc1.1","dE18.Bicc1.2","dE18.Bicc1.3","dE18.Bicc1.4"))
glimpse(melted.ssGSEA.result)

###plot
ggplot(melted.ssGSEA.result, aes(x=variable, y=pathways, fill=value)) +  
  geom_tile() + 
  #scale_fill_gradient2(low = "blue", mid = "white",high = "red") + #color A
  #scale_fill_gradient2(low = "blue", mid = "white",high = "darkred") + #color B
  #scale_fill_scico(palette = 'vikO',name="score", na.value = "grey50")+ #color C 
  scale_fill_paletteer_c("pals::ocean.balance")+ #color D
  scale_y_discrete(limits=  rev(custom.order.mouse.ssGSEA))+ 
  theme(legend.position="bottom", legend.justification = "center")+
  theme(axis.line = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  theme(legend.title.align=0.5)+
  labs(fill = "row z-scored\nNES")+ #legend title
  ggtitle("ssGSEA\ncells expression") +
  theme(plot.title = element_text(hjust = 0.5, face="bold", size=16)) + 
  theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1))+ 
  theme(axis.text.y= element_text(size=14))+  
  xlab(NULL) + 
  ylab(NULL) +
  theme(axis.text.y=element_text(face= "plain", colour="black", size=8) )+
  theme(axis.text.x = element_blank()) +
  theme(axis.ticks.x = element_blank())+
  
  geom_vline(xintercept = c(4.5, 8.5, 12.5, 16.5) , size=1.0, linetype = "solid", color="white")+
  annotate("text", x = 2.5, y = nrow(sign.ssGSEA.result)+2, label = "GFP", fontface = "bold", size=3)+
  annotate("text", x = 6.5, y = nrow(sign.ssGSEA.result)+2, label = "FL", fontface = "bold", size=3)+
  annotate("text", x = 10.5, y = nrow(sign.ssGSEA.result)+2, label = "dE18", fontface = "bold", size=3)+
  annotate("text", x = 14.5, y = nrow(sign.ssGSEA.result)+2, label = "FL\nBicc1", fontface = "bold", size=3)+
  annotate("text", x = 18.5, y = nrow(sign.ssGSEA.result)+2, label = "dE18\nBicc1", fontface = "bold", size=3)+
  annotate("text", x = 3.5, y = nrow(sign.ssGSEA.result)+3.5, label = "",  color = "transparent")

#ggsave("cells.expr.ssGSEA.result.heatmap.pdf", useDingbats=FALSE,  width = 14, height =13, units = "cm") 













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









