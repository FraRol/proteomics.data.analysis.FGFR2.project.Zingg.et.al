#_______________________________________________________________________________________________________________________
# 20.01.2022
# 
# Project: proteomics.truncated.FGFR2.is.oncogene.cancer.Zingg.et.al 
# 
# Script name: tumors.expression.DIA.public
#
# Purpose of script: analysis protein expression data tumors truncated FGFR2 is oncogene cancer project Zingg et al 
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

### long report (FR PG level) no pivot style, WITH filter analysis review checked (yes, on) 8.11.2021 the one to go for
#DZ.tumors.expr_V2 <- fread("/Users/frankrolfs/NLPostDR/VUmc OPL/Phosphoproteomics with Daniel Zingg NKI/mouse tumors/Zing.Fgfr2.tumors.expression.DIA/19.10.21_FR_Zingg_Fgfr2_tumors_DIA_analysis_with.MQ.direct.DIA.library_Report_PG_long_filter.analysis.review.xls",
#                           integer64="numeric",
#                           header=TRUE,
#                           dec=","
#)

### load expresion data cells DIA Zingg
#save(DZ.tumors.expr_V2, file="Zingg.tumors.expression.Spectronaut.report.long.Rdata")
load("Zingg.tumors.expression.Spectronaut.report.long.Rdata")

glimpse(DZ.tumors.expr_V2) #188802

length(unique(DZ.tumors.expr_V2$PG.Genes)) #6757
length(unique(DZ.tumors.expr_V2$PG.ProteinGroups)) #6956


# changes for sample names
DZ.tumors.expr_V2 <- mutate(DZ.tumors.expr_V2, Condition=R.FileName) #
glimpse(DZ.tumors.expr_V2)


# select less columns: PG level only
DZ.tumors.expr_V2 <- select(DZ.tumors.expr_V2, Condition, contains("PG")) %>% distinct()
glimpse(DZ.tumors.expr_V2)
nrow(DZ.tumors.expr_V2) #188802


#detect PG.ProteinGroups with "CON__" and remove them
# CON entries do not have data like PG.Quantity, or Q.value ...
DZ.tumors.expr_V2$contamination <- pbsapply(DZ.tumors.expr_V2$PG.ProteinGroups, function(x) if(str_detect(string= x, pattern="CON__")){"CON"} else {"ok"})
table(DZ.tumors.expr_V2$contamination) #712 CON

DZ.tumors.expr_V2.CON__ <- filter(DZ.tumors.expr_V2, contamination == "CON" )
nrow(DZ.tumors.expr_V2.CON__) #712
glimpse(DZ.tumors.expr_V2.CON__)


DZ.tumors.expr_V2.2 <- filter(DZ.tumors.expr_V2, contamination == "ok" )
nrow(DZ.tumors.expr_V2.2) #188090
glimpse(DZ.tumors.expr_V2.2) 

#quick overview number identifications
print(DZ.tumors.expr_V2.2 %>% group_by(Condition) %>% summarize(number.identifications=n()), n=32)

# get unique PG.ProteinAccessions => needed for long format to wide format transformation step
unique.PG.ProteinGroups.after.quality.filters.V2 <- unique(DZ.tumors.expr_V2.2$PG.ProteinGroups)
glimpse(unique.PG.ProteinGroups.after.quality.filters.V2) #6916



#### long format to wide format spread
DZ.tumors.expr_V2.2.PG.wide <- spread(DZ.tumors.expr_V2.2, key=Condition ,  value = PG.Quantity)
glimpse(DZ.tumors.expr_V2.2.PG.wide) #6916 obs
length(unique(DZ.tumors.expr_V2.2.PG.wide$PG.ProteinGroups)) #6916

### change sample names 

#sample overview
# A	Fgfr2	
# B	Fgfr2-E18-C2	
# C	Fgfr2-dE18	
# D	Fgfr2-Bicc1	
# E	Fgfr2-dE18-Bicc1	
# F	Fgfr2-Ate1	
# G	Fgfr2-dE18-Ate1	
# H	Fgfr2-Tacc2	
# I	Fgfr2-dE18-Tacc2	
# J	Fgfr2-dE18-IGR1	
# K	Fgfr2-dE18-IGR2	
# L	Fgfr2-E18-C3	
# M	Fgfr2-E18-C4	
# N	Fgfr2	2	1.79	
# O	Fgfr2-dE18	
# P	Fgfr2-Bicc1	
# Q	Fgfr2-dE18-Bicc1	
# R	Fgfr2-Ate1	
# S	Fgfr2-dE18-Ate1	
# T	Fgfr2-Tacc2	
# U	Fgfr2-dE18-Tacc2	
# V	Fgfr2-dE18-IGR1	
# W	Fgfr2-dE18-IGR2	
# X	Fgfr2-E18-C2	
# Y	Fgfr2-E18-C3	
# Z	Fgfr2-E18-C4	
# AA	Fgfr2-dE18	
# AB	Fgfr2-Bicc1	
# AC	Fgfr2-dE18-Bicc1
# AD	Fgfr2-Tacc2	
# AE	Fgfr2-dE18-Tacc2
# AF	Fgfr2-E18-C2


# load Zingg tumor expression DIA sample information
load("Zingg.Fgfr2.tumors.sample.info.Rdata")
glimpse(Zingg.Fgfr2.tumors.sample.info)

colnames(DZ.tumors.expr_V2.2.PG.wide)
colnames(DZ.tumors.expr_V2.2.PG.wide)<-c("PG.Qvalue",                                "PG.IsSingleHit",                           "PG.FastaHeaders",                         
                                         "PG.Genes",                                 "PG.Organisms",                            "PG.ProteinAccessions",                    
                                         "PG.FastaFiles",                            "PG.ProteinDescriptions",                   "PG.NrOfStrippedSequencesMeasured",        
                                         "PG.NrOfModifiedSequencesMeasured",         "PG.NrOfPrecursorsMeasured",                "PG.ProteinGroups",                        
                                         "PG.CellularComponent",                     "PG.BiologicalProcess",                     "PG.MolecularFunction",                    
                                         "PG.ProteinNames",                          "PG.UniProtIds",                            "contamination",          
                                         
                                         "A_1_Fgfr2",              "AA_6_Fgfr2-dE18",        "AB_9_Fgfr2-Bicc1",       "AC_12_Fgfr2-dE18-Bicc1", "AD_19_Fgfr2-Tacc2",      "AE_22_Fgfr2-dE18-Tacc2",
                                         "AF_28_Fgfr2-E18-C2",     "B_3_Fgfr2-E18-C2",       "C_4_Fgfr2-dE18",         "D_7_Fgfr2-Bicc1",        "E_10_Fgfr2-dE18-Bicc1",  "F_13_Fgfr2-Ate1",       
                                         "G_15_Fgfr2-dE18-Ate1",   "H_17_Fgfr2-Tacc2",       "I_20_Fgfr2-dE18-Tacc2",  "J_23_Fgfr2-dE18-IGR1",   "K_25_Fgfr2-dE18-IGR2",   "L_29_Fgfr2-E18-C3",     
                                         "M_31_Fgfr2-E18-C4",      "N_2_Fgfr2",              "O_5_Fgfr2-dE18",         "P_8_Fgfr2-Bicc1",        "Q_11_Fgfr2-dE18-Bicc1",  "R_14_Fgfr2-Ate1",       
                                         "S_16_Fgfr2-dE18-Ate1",   "T_18_Fgfr2-Tacc2",       "U_21_Fgfr2-dE18-Tacc2",  "V_24_Fgfr2-dE18-IGR1",   "W_26_Fgfr2-dE18-IGR2",   "X_27_Fgfr2-E18-C2",     
                                         "Y_30_Fgfr2-E18-C3",      "Z_32_Fgfr2-E18-C4"     
)

glimpse(DZ.tumors.expr_V2.2.PG.wide)

### reorder samples
DZ.tumors.expr_V2.2.PG.wide <- DZ.tumors.expr_V2.2.PG.wide %>% select(PG.Qvalue :contamination, 
                                                                      "A_1_Fgfr2",
                                                                      "N_2_Fgfr2",  
                                                                      
                                                                      "F_13_Fgfr2-Ate1",
                                                                      "R_14_Fgfr2-Ate1",  
                                                                      
                                                                      "D_7_Fgfr2-Bicc1",
                                                                      "P_8_Fgfr2-Bicc1", 
                                                                      "AB_9_Fgfr2-Bicc1",  
                                                                      
                                                                      "H_17_Fgfr2-Tacc2",
                                                                      "T_18_Fgfr2-Tacc2", 
                                                                      "AD_19_Fgfr2-Tacc2", 
                                                                      
                                                                      "C_4_Fgfr2-dE18",
                                                                      "O_5_Fgfr2-dE18", 
                                                                      "AA_6_Fgfr2-dE18", 
                                                                      
                                                                      "G_15_Fgfr2-dE18-Ate1", 
                                                                      "S_16_Fgfr2-dE18-Ate1",
                                                                      
                                                                      "E_10_Fgfr2-dE18-Bicc1", 
                                                                      "Q_11_Fgfr2-dE18-Bicc1", 
                                                                      "AC_12_Fgfr2-dE18-Bicc1",
                                                                      
                                                                      "I_20_Fgfr2-dE18-Tacc2",
                                                                      "U_21_Fgfr2-dE18-Tacc2",
                                                                      "AE_22_Fgfr2-dE18-Tacc2",
                                                                      
                                                                      "J_23_Fgfr2-dE18-IGR1",
                                                                      "V_24_Fgfr2-dE18-IGR1",
                                                                      
                                                                      "K_25_Fgfr2-dE18-IGR2",
                                                                      "W_26_Fgfr2-dE18-IGR2",
                                                                      
                                                                      "B_3_Fgfr2-E18-C2",      
                                                                      "X_27_Fgfr2-E18-C2",  
                                                                      "AF_28_Fgfr2-E18-C2",
                                                                      
                                                                      "L_29_Fgfr2-E18-C3",
                                                                      "Y_30_Fgfr2-E18-C3",
                                                                      
                                                                      "M_31_Fgfr2-E18-C4", 
                                                                      "Z_32_Fgfr2-E18-C4" )

glimpse(DZ.tumors.expr_V2.2.PG.wide)


### boxplot and density plot samples to check normalization #####################################################################################

data.for.boxplot <- DZ.tumors.expr_V2.2.PG.wide %>% select(A_1_Fgfr2   : `Z_32_Fgfr2-E18-C4`)
data.for.boxplot <- log2(data.for.boxplot)
glimpse(data.for.boxplot)



# reshape for plot
melted.data.for.boxplot <- reshape2::melt(data.for.boxplot)
glimpse(melted.data.for.boxplot)
unique(melted.data.for.boxplot$variable)

#box plot
ggplot(data = melted.data.for.boxplot)+
  geom_boxplot(mapping = aes(x= variable, y = value, fill=variable), outlier.size = 0.5, notch = TRUE, notchwidth = 0.1)+
  scale_fill_viridis(discrete=T, option="inferno")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1.0, vjust = 0.5))+
  theme(legend.position="none", legend.justification = "center")+
  ggtitle("mouse tumors expression DIA") +
  theme(plot.title = element_text(hjust = 0.5, face= "bold"))+ 
  theme(axis.text.x = element_text(face= "plain", colour="black", size=14))+
  theme(axis.text.y = element_text(face= "plain", colour="black", size=14))+
  xlab(NULL)+ 
  ylab("log2(Int.)") 


# density plot
glimpse(melted.data.for.boxplot)
ggplot(melted.data.for.boxplot, aes(x=value, colour = variable) ) + 
  geom_density() +
  scale_color_viridis(discrete=T, option="inferno")+
  theme(legend.position="bottom")+
  ggtitle("mouse tumors expression DIA \n data distribution") +
  theme(plot.title = element_text(hjust = 0.5, face= "bold"))+
  theme(axis.text.x = element_text(face= "plain", colour="black", size=14))+
  theme(axis.text.y = element_text(face= "plain", colour="black", size=14))


### count indentifications per sample ###################################################################################################################################
glimpse(DZ.tumors.expr_V2.2.PG.wide)
colnames(DZ.tumors.expr_V2.2.PG.wide)

result.EXPR.count.LFQ <- c()
sample.name.count <- 1
for(i in 19:50){ 
  print(i)
  sample.name <- colnames(DZ.tumors.expr_V2.2.PG.wide)
  sample.name <- sample.name[i]
  print(sample.name)
  
  temp <- DZ.tumors.expr_V2.2.PG.wide %>% select(sample.name) %>% pull()
  temp.data <- sum(!is.na(temp))
  temp.na<- sum(is.na(temp))
  temp.result <- tibble(sample=sample.name, identifications=temp.data, MVs=temp.na)
  
  result.EXPR.count.LFQ  <- bind_rows(result.EXPR.count.LFQ, temp.result)
}
print(result.EXPR.count.LFQ, n=32)
mean(result.EXPR.count.LFQ$identifications) #5877.406


#add group information to count results
result.EXPR.count.LFQ$variant.group <- c("Fgfr2",            "Fgfr2",            "Fgfr2-Ate1",       "Fgfr2-Ate1",       "Fgfr2-Bicc1",      "Fgfr2-Bicc1",      "Fgfr2-Bicc1",      "Fgfr2-Tacc2",     
                                         "Fgfr2-Tacc2",      "Fgfr2-Tacc2",      "Fgfr2-dE18",       "Fgfr2-dE18",       "Fgfr2-dE18",       "Fgfr2-dE18-Ate1",  "Fgfr2-dE18-Ate1",  "Fgfr2-dE18-Bicc1",
                                         "Fgfr2-dE18-Bicc1", "Fgfr2-dE18-Bicc1", "Fgfr2-dE18-Tacc2", "Fgfr2-dE18-Tacc2", "Fgfr2-dE18-Tacc2", "Fgfr2-dE18-IGR1",  "Fgfr2-dE18-IGR1",  "Fgfr2-dE18-IGR2", 
                                         "Fgfr2-dE18-IGR2",  "Fgfr2-E18-C2",     "Fgfr2-E18-C2",     "Fgfr2-E18-C2",     "Fgfr2-E18-C3",     "Fgfr2-E18-C3",     "Fgfr2-E18-C4",     "Fgfr2-E18-C4")
result.EXPR.count.LFQ

### plot numbers expression
ggplot(data = result.EXPR.count.LFQ )+
  geom_col(aes(x=sample, y=identifications, fill=variant.group), color="black")+
  theme(legend.position="none") +
  scale_x_discrete(limits= c( "A_1_Fgfr2",
                              "N_2_Fgfr2",  
                              
                              "F_13_Fgfr2-Ate1",
                              "R_14_Fgfr2-Ate1",  
                              
                              "D_7_Fgfr2-Bicc1",
                              "P_8_Fgfr2-Bicc1", 
                              "AB_9_Fgfr2-Bicc1",  
                              
                              "H_17_Fgfr2-Tacc2",
                              "T_18_Fgfr2-Tacc2", 
                              "AD_19_Fgfr2-Tacc2", 
                              
                              "C_4_Fgfr2-dE18",
                              "O_5_Fgfr2-dE18", 
                              "AA_6_Fgfr2-dE18", 
                              
                              "G_15_Fgfr2-dE18-Ate1", 
                              "S_16_Fgfr2-dE18-Ate1",
                              
                              "E_10_Fgfr2-dE18-Bicc1", 
                              "Q_11_Fgfr2-dE18-Bicc1", 
                              "AC_12_Fgfr2-dE18-Bicc1",
                              
                              "I_20_Fgfr2-dE18-Tacc2",
                              "U_21_Fgfr2-dE18-Tacc2",
                              "AE_22_Fgfr2-dE18-Tacc2",
                              
                              "J_23_Fgfr2-dE18-IGR1",
                              "V_24_Fgfr2-dE18-IGR1",
                              
                              "K_25_Fgfr2-dE18-IGR2",
                              "W_26_Fgfr2-dE18-IGR2",
                              
                              "B_3_Fgfr2-E18-C2",      
                              "X_27_Fgfr2-E18-C2",  
                              "AF_28_Fgfr2-E18-C2",
                              
                              "L_29_Fgfr2-E18-C3",
                              "Y_30_Fgfr2-E18-C3",
                              
                              "M_31_Fgfr2-E18-C4", 
                              "Z_32_Fgfr2-E18-C4"))+
  scale_fill_viridis(discrete=T, option="inferno")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1.0, vjust = 0.5))+
  ylab("number PS")+
  xlab(NULL)+
  scale_y_continuous(breaks=seq(from=0, to=6600, by=500), limits = c(0,6600))+
  ggtitle("mouse tumors - identifications per sample" )+
  theme(plot.title = element_text(hjust = 0.5, face= "bold"))+  #center title
  theme(axis.text.x = element_text(face= "plain", colour="black", size=14))+
  theme(axis.text.y = element_text(face= "plain", colour="black", size=14))+
  geom_vline(xintercept = c(2.5, 4.5, 7.5, 10.5, 13.5, 15.5, 18.5, 21.5, 23.5, 25.5, 28.5, 30.5) , size=0.25, linetype = 2)


#############################################################################################################################################################
#############################################################################################################################################################
### correlation EXPRESSION proteins

glimpse(DZ.tumors.expr_V2.2.PG.wide)

#select columns and log2 transform
DF_for_correlation_plot.EXPR.Prot <- DZ.tumors.expr_V2.2.PG.wide %>% select(A_1_Fgfr2  :  `Z_32_Fgfr2-E18-C4`)
DF_for_correlation_plot.EXPR.Prot <- log2(DF_for_correlation_plot.EXPR.Prot)
glimpse(DF_for_correlation_plot.EXPR.Prot)

#prepare better sample names
tum.expr.name.change <- tibble(orig.col.name = colnames(DF_for_correlation_plot.EXPR.Prot))
tum.expr.name.change$new.col.name.1 <- pbsapply(tum.expr.name.change$orig.col.name, function(x) unlist(str_split(string=x, pattern="_"))[3]   )
tum.expr.name.change$new.col.name.2 <- paste0(tum.expr.name.change$new.col.name.1, ".", c(1, 2,   1, 2,   1, 2, 3,   1, 2, 3,    1, 2, 3,   1, 2,   1, 2, 3,    1, 2, 3,    1, 2,  1, 2,   1, 2, 3 ,   1, 2,   1, 2))
tum.expr.name.change

colnames(DF_for_correlation_plot.EXPR.Prot) <- tum.expr.name.change$new.col.name.2

##check NA rows 
DF_for_correlation_plot.EXPR.Prot$na.rows <- pbapply(DF_for_correlation_plot.EXPR.Prot, 1, function(x) sum(is.na(x)))
glimpse(DF_for_correlation_plot.EXPR.Prot)
table(DF_for_correlation_plot.EXPR.Prot$na.rows)

DF_for_correlation_plot.EXPR.Prot <- DF_for_correlation_plot.EXPR.Prot %>% select(-na.rows )
glimpse(DF_for_correlation_plot.EXPR.Prot)
colnames(DF_for_correlation_plot.EXPR.Prot)


### correlation
corr2 <- cor(DF_for_correlation_plot.EXPR.Prot, method = "pearson", use = "na.or.complete")
glimpse(corr2)
head(corr2)
min(corr2) #0.3725943
max(corr2)
round(min(corr2), 2)



### prepare heatmap annotation
annot_df_for_heatmap <- data.frame(samples = colnames(DF_for_correlation_plot.EXPR.Prot),
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

glimpse(annot_df_for_heatmap.short)
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

### Create the heatmap annotation
ha.fr <- HeatmapAnnotation(group=annot_df_for_heatmap.short$group,
                           treatment=annot_df_for_heatmap.short$treatment,
                           col = annot.colors,
                           annotation_legend_param = list(grid_height = unit(8, "mm"))
)
ha.fr


#prepare heatmap colors
FR.heatmap.colors.2 <-colorRamp2(c(min(corr2), 1.0), c("grey90", "grey10"))

#plot
#pdf("DZ.tumors.expr.DIA.correlation.pdf", width=25/2.54, height=20/2.54, useDingbats=FALSE)
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




##################################################################################################################################
##################################################################################################################################
### multiple candidates plot to check expression


# multiple heatmap plot: ###################################################################################################################################
# multiple heatmap plot: ###################################################################################################################################

### define expression candidates to plot A ###
mtpl.candidates <- c("FGFR2_MOUSE", "BICC1_MOUSE", "TACC2_MOUSE", "ATE1_MOUSE") 

###filter for candidates
EXPR.mtplplot.data <- DZ.tumors.expr_V2.2.PG.wide %>% filter(PG.ProteinNames %in% mtpl.candidates)
glimpse(EXPR.mtplplot.data)

### Reorder data according to order candidates, average technical replicates, log2 transform, reshape for plot
EXPR.mtplplot.data <- left_join(data.frame(PG.ProteinNames = mtpl.candidates),
                                EXPR.mtplplot.data,
                                by = "PG.ProteinNames")

EXPR.mtplplot.data.2A <- EXPR.mtplplot.data  %>% select(PG.Genes, PG.ProteinNames,PG.UniProtIds
) 
glimpse(EXPR.mtplplot.data.2A)

EXPR.mtplplot.data.2B <- EXPR.mtplplot.data  %>% select(
  "A_1_Fgfr2",
  "N_2_Fgfr2",  
  
  "F_13_Fgfr2-Ate1",
  "R_14_Fgfr2-Ate1",  
  
  "D_7_Fgfr2-Bicc1",
  "P_8_Fgfr2-Bicc1", 
  "AB_9_Fgfr2-Bicc1",  
  
  "H_17_Fgfr2-Tacc2",
  "T_18_Fgfr2-Tacc2", 
  "AD_19_Fgfr2-Tacc2", 
  
  "C_4_Fgfr2-dE18",
  "O_5_Fgfr2-dE18", 
  "AA_6_Fgfr2-dE18", 
  
  "G_15_Fgfr2-dE18-Ate1", 
  "S_16_Fgfr2-dE18-Ate1",
  
  "E_10_Fgfr2-dE18-Bicc1", 
  "Q_11_Fgfr2-dE18-Bicc1", 
  "AC_12_Fgfr2-dE18-Bicc1",
  
  "I_20_Fgfr2-dE18-Tacc2",
  "U_21_Fgfr2-dE18-Tacc2",
  "AE_22_Fgfr2-dE18-Tacc2",
  
  "J_23_Fgfr2-dE18-IGR1",
  "V_24_Fgfr2-dE18-IGR1",
  
  "K_25_Fgfr2-dE18-IGR2",
  "W_26_Fgfr2-dE18-IGR2",
  
  "B_3_Fgfr2-E18-C2",      
  "X_27_Fgfr2-E18-C2",  
  "AF_28_Fgfr2-E18-C2",
  
  "L_29_Fgfr2-E18-C3",
  "Y_30_Fgfr2-E18-C3",
  
  "M_31_Fgfr2-E18-C4", 
  "Z_32_Fgfr2-E18-C4") 
glimpse(EXPR.mtplplot.data.2B)


#log2 data transformation
EXPR.mtplplot.data.2D <- log2(EXPR.mtplplot.data.2B)
glimpse(EXPR.mtplplot.data.2D)

# combine data frames
EXPR.mtplplot.data.3 <- bind_cols(EXPR.mtplplot.data.2A, EXPR.mtplplot.data.2D) 
glimpse(EXPR.mtplplot.data.3)
nrow(EXPR.mtplplot.data.3)

# prepare additional candidate names
EXPR.mtplplot.data.3$ProteinNames_Genes_UniprotID <- paste0(EXPR.mtplplot.data.3$PG.Genes, "_", EXPR.mtplplot.data.3$PG.UniProtIds)
glimpse(EXPR.mtplplot.data.3)

# define order in plot
order.in.plot <- EXPR.mtplplot.data.3$ProteinNames_Genes_UniprotID 
order.in.plot

order.in.plot2 <- EXPR.mtplplot.data.3$PG.Genes 
order.in.plot2

# prepare additional sample names and order
name.change.DZ.tum.expr.dia <- tibble(orig.col.name = colnames(EXPR.mtplplot.data.3 %>% select(-PG.Genes, -PG.ProteinNames, -PG.UniProtIds, -ProteinNames_Genes_UniprotID)) )
name.change.DZ.tum.expr.dia$new.col.name.1 <- pbsapply(name.change.DZ.tum.expr.dia$orig.col.name, function(x) unlist(str_split(string=x, pattern="_"))[3]   )
name.change.DZ.tum.expr.dia$new.col.name.2 <- paste0(name.change.DZ.tum.expr.dia$new.col.name.1, ".", c(1, 2,   1, 2,   1, 2, 3,   1, 2, 3,    1, 2, 3,   1, 2,   1, 2, 3,    1, 2, 3,    1, 2,  1, 2,   1, 2, 3 ,   1, 2,   1, 2))
print(name.change.DZ.tum.expr.dia, n=32)


custom.sample.order.DZ.tum.expr.dia <- c("Fgfr2.1", "Fgfr2.2", 
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
custom.sample.order.DZ.tum.expr.dia 

### prepare and reshape data for plot
glimpse(EXPR.mtplplot.data.3)
EXPR.mtplplot.data.4 <- EXPR.mtplplot.data.3
colnames(EXPR.mtplplot.data.4) <-c("PG.Genes", "PG.ProteinNames", "PG.UniProtIds", name.change.DZ.tum.expr.dia$new.col.name.2,   "ProteinNames_Genes_UniprotID")

melted.EXPR.mtplplot.data.4 <- reshape2::melt(EXPR.mtplplot.data.4)
glimpse(melted.EXPR.mtplplot.data.4)


#heatmap plot DZ Expression Fgfr2 Bicc1 Tacc2 Ate 1 
ggplot(
  melted.EXPR.mtplplot.data.4 , aes(x=variable, y=PG.Genes, fill=value)) + 
  geom_tile() + 
  scale_fill_viridis(limits = c(6, max(melted.EXPR.mtplplot.data.4$value, na.rm=T)), option = "inferno", na.value = "grey30",name="log2(intensity)") +
  coord_equal()+
  scale_y_discrete(limits=  rev(order.in.plot2))+
  scale_x_discrete(limits=  custom.sample.order.DZ.tum.expr.dia)+
  theme(legend.position="bottom", legend.justification = "center")+
  theme(axis.line = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  theme(legend.title.align=0.5)+
  theme(plot.title = element_text(hjust = 0.5, face="bold", size=16)) + 
  theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1))+ 
  theme(axis.text.y= element_text(size=14))+
  xlab(NULL) + 
  xlab(NULL) + 
  ylab(NULL) +
  geom_vline(xintercept = c(2.5, 5.5, 8.5, 11.5, 13.5, 15.5, 18.5, 21.5, 23.5, 25.5, 28.5, 30.5) , size=1.0, linetype = "solid", color="white")+
  annotate("text", x = 1.5, y = nrow(EXPR.mtplplot.data.4)+1, label = "Fgfr2", fontface = "bold", angle='90', hjust = 0, size=5)+
  annotate("text", x = 4, y = nrow(EXPR.mtplplot.data.4)+1, label = "Fgfr2-dE18", fontface = "bold", angle='90', hjust = 0, size=5)+
  annotate("text", x = 7, y = nrow(EXPR.mtplplot.data.4)+1, label = "Fgfr2-Bicc1", fontface = "bold", angle='90', hjust = 0, size=5)+
  annotate("text", x = 10, y = nrow(EXPR.mtplplot.data.4)+1, label = "Fgfr2-dE18-Bicc1", fontface = "bold", angle='90', hjust = 0, size=5)+
  annotate("text", x = 12.5, y = nrow(EXPR.mtplplot.data.4)+1, label = "Fgfr2-Ate1", fontface = "bold", angle='90', hjust = 0, size=5)+
  annotate("text", x = 14.5, y = nrow(EXPR.mtplplot.data.4)+1, label = "Fgfr2-dE18-Ate1", fontface = "bold", angle='90', hjust = 0, size=5)+
  annotate("text", x = 17, y = nrow(EXPR.mtplplot.data.4)+1, label = "Fgfr2-Tacc2", fontface = "bold", angle='90', hjust = 0, size=5)+
  annotate("text", x = 20, y = nrow(EXPR.mtplplot.data.4)+1, label = "Fgfr2-dE18-Tacc2", fontface = "bold", angle='90', hjust = 0, size=5)+
  annotate("text", x = 22.5, y = nrow(EXPR.mtplplot.data.4)+1, label = "Fgfr2-dE18-IGR1", fontface = "bold", angle='90', hjust = 0, size=5)+
  annotate("text", x = 24.5, y = nrow(EXPR.mtplplot.data.4)+1, label = "Fgfr2-dE18-IGR2", fontface = "bold", angle='90', hjust = 0, size=5)+
  annotate("text", x = 27, y = nrow(EXPR.mtplplot.data.4)+1, label = "Fgfr2-E18-C2", fontface = "bold", angle='90', hjust = 0, size=5)+
  annotate("text", x = 29.5, y = nrow(EXPR.mtplplot.data.4)+1, label = "Fgfr2-E18-C3", fontface = "bold", angle='90', hjust = 0, size=5)+
  annotate("text", x = 31.5, y = nrow(EXPR.mtplplot.data.4)+1, label = "Fgfr2-E18-C4", fontface = "bold", angle='90', hjust = 0, size=5)+
  annotate("text", x = 32, y = nrow(EXPR.mtplplot.data.4)+10, label = "",  color = "transparent")

#ggsave("DZ.tumors.expression.DIA.candidatesA.pdf", useDingbats=FALSE,  width = 18, height =14, units = "cm")




### define expression candidates to plot B ###

mtpl.candidates <- c("AKTS1_MOUSE", "ACLY_MOUSE",  "BAD_MOUSE",   "F262_MOUSE",  "STX7_MOUSE",  "TBCD1_MOUSE", 
                     
                     "CCDC6_MOUSE", "CDK7_MOUSE",  "COR1A_MOUSE",
                     "RUNX1_MOUSE", "SMAD3_MOUSE", "TP53B_MOUSE", 
                     
                     "CTNB1_MOUSE", "CTND1_MOUSE", "HDAC2_MOUSE",
                     
                     "PAK2_MOUSE",  "FAK1_MOUSE",  "RIPK1_MOUSE",
                     "MP2K1_MOUSE", "MP2K2_MOUSE", "MP2K4_MOUSE", "MK01_MOUSE",  "MK03_MOUSE",  "KSR1_MOUSE",  "JUN_MOUSE",   "JUNB_MOUSE",  "EPS8_MOUSE", 
                     "SP1_MOUSE",   "SP3_MOUSE",   "GSK3B_MOUSE", "MAP1B_MOUSE", 	
                     "MAP2_MOUSE",
                     
                     "KS6A1_MOUSE", "KS6A3_MOUSE", "RS6_MOUSE",   "IF4B_MOUSE",  "4EBP1_MOUSE", "IF2B1_MOUSE",
                     "PDCD4_MOUSE", "TIF1B_MOUSE") #
mtpl.candidates 

###filter for candidates
EXPR.mtplplot.data <- DZ.tumors.expr_V2.2.PG.wide %>% filter(PG.ProteinNames %in% mtpl.candidates)
glimpse(EXPR.mtplplot.data)

### Reorder data according to order candidates, average technical replicates, log2 transform, reshape for plot
EXPR.mtplplot.data <- left_join(data.frame(PG.ProteinNames = mtpl.candidates),
                                EXPR.mtplplot.data,
                                by = "PG.ProteinNames")

EXPR.mtplplot.data.2A <- EXPR.mtplplot.data  %>% select(PG.Genes, PG.ProteinNames,PG.UniProtIds
) 
glimpse(EXPR.mtplplot.data.2A)

EXPR.mtplplot.data.2B <- EXPR.mtplplot.data  %>% select(
  "A_1_Fgfr2",
  "N_2_Fgfr2",  
  
  "F_13_Fgfr2-Ate1",
  "R_14_Fgfr2-Ate1",  
  
  "D_7_Fgfr2-Bicc1",
  "P_8_Fgfr2-Bicc1", 
  "AB_9_Fgfr2-Bicc1",  
  
  "H_17_Fgfr2-Tacc2",
  "T_18_Fgfr2-Tacc2", 
  "AD_19_Fgfr2-Tacc2", 
  
  "C_4_Fgfr2-dE18",
  "O_5_Fgfr2-dE18", 
  "AA_6_Fgfr2-dE18", 
  
  "G_15_Fgfr2-dE18-Ate1", 
  "S_16_Fgfr2-dE18-Ate1",
  
  "E_10_Fgfr2-dE18-Bicc1", 
  "Q_11_Fgfr2-dE18-Bicc1", 
  "AC_12_Fgfr2-dE18-Bicc1",
  
  "I_20_Fgfr2-dE18-Tacc2",
  "U_21_Fgfr2-dE18-Tacc2",
  "AE_22_Fgfr2-dE18-Tacc2",
  
  "J_23_Fgfr2-dE18-IGR1",
  "V_24_Fgfr2-dE18-IGR1",
  
  "K_25_Fgfr2-dE18-IGR2",
  "W_26_Fgfr2-dE18-IGR2",
  
  "B_3_Fgfr2-E18-C2",      
  "X_27_Fgfr2-E18-C2",  
  "AF_28_Fgfr2-E18-C2",
  
  "L_29_Fgfr2-E18-C3",
  "Y_30_Fgfr2-E18-C3",
  
  "M_31_Fgfr2-E18-C4", 
  "Z_32_Fgfr2-E18-C4") 
glimpse(EXPR.mtplplot.data.2B)


#log2 data transformation
EXPR.mtplplot.data.2D <- log2(EXPR.mtplplot.data.2B)
glimpse(EXPR.mtplplot.data.2D)

# combine data frames
EXPR.mtplplot.data.3 <- bind_cols(EXPR.mtplplot.data.2A, EXPR.mtplplot.data.2D) 
glimpse(EXPR.mtplplot.data.3)
nrow(EXPR.mtplplot.data.3)

# prepare additional candidate names
EXPR.mtplplot.data.3$ProteinNames_Genes_UniprotID <- paste0(EXPR.mtplplot.data.3$PG.Genes, "_", EXPR.mtplplot.data.3$PG.UniProtIds)
glimpse(EXPR.mtplplot.data.3)

# define order in plot
order.in.plot <- EXPR.mtplplot.data.3$ProteinNames_Genes_UniprotID 
order.in.plot

order.in.plot2 <- EXPR.mtplplot.data.3$PG.Genes 
order.in.plot2

# prepare additional sample names and order
name.change.DZ.tum.expr.dia <- tibble(orig.col.name = colnames(EXPR.mtplplot.data.3 %>% select(-PG.Genes, -PG.ProteinNames, -PG.UniProtIds, -ProteinNames_Genes_UniprotID)) )
name.change.DZ.tum.expr.dia$new.col.name.1 <- pbsapply(name.change.DZ.tum.expr.dia$orig.col.name, function(x) unlist(str_split(string=x, pattern="_"))[3]   )
name.change.DZ.tum.expr.dia$new.col.name.2 <- paste0(name.change.DZ.tum.expr.dia$new.col.name.1, ".", c(1, 2,   1, 2,   1, 2, 3,   1, 2, 3,    1, 2, 3,   1, 2,   1, 2, 3,    1, 2, 3,    1, 2,  1, 2,   1, 2, 3 ,   1, 2,   1, 2))
print(name.change.DZ.tum.expr.dia, n=32)

custom.sample.order.DZ.tum.expr.dia <- c("Fgfr2.1", "Fgfr2.2", 
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
custom.sample.order.DZ.tum.expr.dia 

### prepare and reshape data for plot
glimpse(EXPR.mtplplot.data.3)
EXPR.mtplplot.data.4 <- EXPR.mtplplot.data.3
colnames(EXPR.mtplplot.data.4) <-c("PG.Genes", "PG.ProteinNames", "PG.UniProtIds", name.change.DZ.tum.expr.dia$new.col.name.2,   "ProteinNames_Genes_UniprotID")

melted.EXPR.mtplplot.data.4 <- reshape2::melt(EXPR.mtplplot.data.4)
glimpse(melted.EXPR.mtplplot.data.4)

### heatmap plot expression candidates B 

#exlude Ctnnd1_P30999-2;P30999-3
melted.EXPR.mtplplot.data.4 <- melted.EXPR.mtplplot.data.4 %>% filter(!ProteinNames_Genes_UniprotID %in% c("Ctnnd1_P30999-2;P30999-3"))

ggplot(
  melted.EXPR.mtplplot.data.4 , aes(x=variable, y=PG.Genes , fill=value)) + 
  geom_tile() + 
  scale_fill_viridis( option = "inferno", na.value = "grey30",name="log2(intensity)") + 
  coord_equal()+
  scale_y_discrete(limits=  rev(c(
    "Akt1s1",   "Acly",     "Bad",      "Pfkfb2",   "Stx7",     "Tbc1d1",   "Ccdc6",    "Cdk7",     "Coro1a",   "Runx1",    "Smad3",   
    "Tp53bp1",  "Ctnnb1",   "Ctnnd1",               "Hdac2",    "Pak2",     "Ptk2",     "Ripk1",    "Map2k1",   "Map2k2",   "Map2k4",  
    "Mapk1",    "Mapk3",    "Ksr1",     "Jun",      "Junb",     "Eps8",     "Sp1",      "Sp3",      "Gsk3b",    "Map1b",    "Metap2",  
    "Rps6ka1",  "Rps6ka3",  "Rps6",     "Eif4b",    "Eif4ebp1", "Igf2bp1",  "Pdcd4",    "Trim28"  ))   )+
  scale_x_discrete(limits=  c(    
    #  #Group 1
    #  "1_Fgfr2","2_Fgfr2", "13_Fgfr2-Ate1","14_Fgfr2-Ate1",
    "Fgfr2.1", "Fgfr2.2", "Fgfr2-Ate1.1", "Fgfr2-Ate1.2",
    
    #  #Group 3
    #  "7_Fgfr2-Bicc1","8_Fgfr2-Bicc1","9_Fgfr2-Bicc1","17_Fgfr2-Tacc2","18_Fgfr2-Tacc2","19_Fgfr2-Tacc2",
    "Fgfr2-Bicc1.1", "Fgfr2-Bicc1.2", "Fgfr2-Bicc1.3", "Fgfr2-Tacc2.1", "Fgfr2-Tacc2.2", "Fgfr2-Tacc2.3",
    
    #  #Group 5
    #  "3_Fgfr2-E18-C2","27_Fgfr2-E18-C2","28_Fgfr2-E18-C2",
    "Fgfr2-E18-C2.1", "Fgfr2-E18-C2.2", "Fgfr2-E18-C2.3",
    
    #  #Group 2
    #  "4_Fgfr2-dE18","5_Fgfr2-dE18","6_Fgfr2-dE18","23_Fgfr2-dE18-IGR1","24_Fgfr2-dE18-IGR1","25_Fgfr2-dE18-IGR2","26_Fgfr2-dE18-IGR2","29_Fgfr2-E18-C3","30_Fgfr2-E18-C3","31_Fgfr2-E18-C4","32_Fgfr2-E18-C4",
    "Fgfr2-dE18.1", "Fgfr2-dE18.2", "Fgfr2-dE18.3", "Fgfr2-dE18-IGR1.1", "Fgfr2-dE18-IGR1.2", "Fgfr2-dE18-IGR2.1", "Fgfr2-dE18-IGR2.2", "Fgfr2-E18-C3.1", "Fgfr2-E18-C3.2", "Fgfr2-E18-C4.1", "Fgfr2-E18-C4.2",
    
    #  #Group 4
    #  "15_Fgfr2-dE18-Ate1","16_Fgfr2-dE18-Ate1","10_Fgfr2-dE18-Bicc1","11_Fgfr2-dE18-Bicc1","12_Fgfr2-dE18-Bicc1","20_Fgfr2-dE18-Tacc2","21_Fgfr2-dE18-Tacc2","22_Fgfr2-dE18-Tacc2"
    "Fgfr2-dE18-Ate1.1", "Fgfr2-dE18-Ate1.2", "Fgfr2-dE18-Bicc1.1", "Fgfr2-dE18-Bicc1.2", "Fgfr2-dE18-Bicc1.3", "Fgfr2-dE18-Tacc2.1", "Fgfr2-dE18-Tacc2.2", "Fgfr2-dE18-Tacc2.3"
  ) )+
  theme(legend.position="bottom", legend.justification = "center")+
  theme(axis.line = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  theme(legend.title.align=0.5)+
  theme(plot.title = element_text(hjust = 0.5, face="bold", size=16)) + #center title
  theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1))+  ## Vertical text on x axis
  theme(axis.text.y= element_text(size=14))+  ## Vertical text on x axis
  xlab(NULL) + 
  xlab(NULL) + 
  ylab(NULL) +
  geom_vline(xintercept = c(      10.5, 13.5, 24.5) , size=1.0, linetype = "solid", color="white")+ #G1+G3, G5, G2, G4
  annotate("text", x = 5, y = nrow(EXPR.mtplplot.data.4)+2, label = "G1&G3", fontface = "bold")+
  annotate("text", x = 12, y = nrow(EXPR.mtplplot.data.4)+2, label = "G5", fontface = "bold")+
  annotate("text", x = 19.0, y = nrow(EXPR.mtplplot.data.4)+2, label = "G2", fontface = "bold")+
  annotate("text", x = 28.5, y = nrow(EXPR.mtplplot.data.4)+2, label = "G4", fontface = "bold")+
  annotate("text", x = 3.5, y = nrow(EXPR.mtplplot.data.4)+3, label = "",  color = "transparent")

#ggsave("DZ.tumors.expression.DIA.candidatesB.pdf", useDingbats=FALSE,  width = 18, height =22, units = "cm")









