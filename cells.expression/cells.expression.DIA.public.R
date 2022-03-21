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
                     "BICC1_MOUSE",
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




















