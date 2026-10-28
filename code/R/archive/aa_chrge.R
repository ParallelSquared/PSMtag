
# Packages used: 

list.of.packages <- c("reshape2", "ggridges","tidyverse","ggpubr","ggpattern",
                      "arrow","gridExtra","grid","eulerr","ggh4x","xml2",
                      "googledrive","readxl","dplyr","patchwork")
new.packages <- list.of.packages[!(list.of.packages %in%
                                     installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)


#remotes::install_github("coolbutuseless/ggpattern")
#install.packages("ggh4x")
library(reshape2)
library(ggridges)
library(patchwork)
library(dplyr)
library(tidyverse)
library(ggpubr)
library(ggpattern)
library(arrow)
library(gridExtra)
library(grid)
library(eulerr)
library(ggh4x)
library(xml2)
library(googledrive)
library(readxl)

# Analysis-wide paramters 

theme_tag<-function() {
  theme_pubr() + theme(text = element_text(size=18))
}

image_px_w<-500
image_px_h<-500

image_in_w<-5
image_in_h<-5


################################################################################################################
# remove.duplicates - remove duplicates across multiple columns
################################################################################################################

# Find unique entries by column1 and column2: 
# Ex remove.duplicates(TMT50,c("Sequence","Charge"))
remove.duplicates<-function(data,Cols){
  
  return(data[!duplicated(data[,Cols]),])
  
}







ev<-read.delim("/Volumes/Lab/MA/new_tags/235_sage/235-4_24/results.sage.tsv")

head(ev)

summary(ev$peptide_len)

# Correlating pKa of AA to labeling efficiency:

AA<-c("A","D","C","E","G","I","L","K","F","S","T","Y","V","R","N","Q","H","M","P","W" )
AA<-AA[order(AA)]
pKa_C<-c(2.35,2.09,1.71,2.19,2.34,2.36,2.36,2.18,1.83,2.21,2.63, 2.20,2.32,2.17,2.02,2.17, 5.0,2.28,1.99,2.38)
pKa_N<-c(9.69,9.82, 8.33,9.67,9.6, 9.68, 9.60, 8.95,9.13,9.15, 10.43, 9.11,9.62,9.04, 8.8,9.13,9.7,9.21,10.60,9.39)
hAA<-c("I",	"V",	"L","F",	"C",	"M",	"A",	"W",	"G",	"T",	"S",	"Y",	"P",	"H",	"N",	"D",	"Q",	"E",	"K",	"R")
hydropathy<-c(4.5,	4.2,	3.8,	2.8,	2.5,	1.9,	1.8,	-0.9,	-0.4,	-0.7,	-0.8,	-1.3,	-1.6,	-3.2,	-3.5,	-3.5,	-3.5,	-3.5,	-3.9,	-4.5)

library(stringr)

ev<-ev[!grepl("K",ev$peptide),]
#ev<-ev[grepl("[+657.35986]",ev$peptide, fixed=T),]
ev<-ev[!grepl("[+657.35986]",ev$peptide, fixed=T),]
#ev$peptide<-gsub("[+657.35986]-","",ev$peptide,fixed=T)
ev<-ev[ev$missed_cleavages==0,]
ev<-ev[!grepl("]",ev$peptide, fixed=T),]
#ev<-ev[ev$peptide_len>1 & ev$peptide_len<20,]
ev<-ev[ev$charge==4,]



ev2 <- ev%>%group_by(peptide_len) %>% summarise(mch = mean(charge))

ggplot(ev2,aes(x=peptide_len, y=mch)) + 
  geom_point() + 
  theme_tag() + 
  ylab("Mean Charge\n") + 
  xlab("\nPeptide length") 
  



# play
ev$peptide<-str_sub(ev$peptide,1,str_length(ev$peptide)-1)
#ev$peptide<-str_sub(ev$peptide,1,6)
#ev$peptide<-str_sub(ev$peptide,str_length(ev$peptide)-4,str_length(ev$peptide)-1)

#hist(ev$peptide_len)


ev[,AA]<-NA
for(X in AA){
  
  ev[,X]<-str_count(ev$peptide,X) 
  #ev[,X]<-str_count(ev$Nterm,X) 
  
}

#ev<-ev[ev$peptide_len>6,]

cors<-c()
for(X in AA){
  
  cors<-c(cors, cor(ev[,X]/ev$peptide_len, ev$peptide_len ,method="spearman"))
  #cors<-c(cors, cor(ev[,X], ev$peptide_len ,method="spearman"))
  
}

names(cors)<-AA

print(cors)

hdf<-data.frame(hAA, hydropathy)
adf<-data.frame(names(cors),cors)
df<-merge(adf,hdf, by.x="names.cors.", by.y="hAA")

plot(y=df$cors, df$hydropathy, ylim=c(-0.1,0.5))


h_order<-hydropathy[order(hAA)]
names(h_order)<-hAA[order(hAA)]

mat<-as.matrix(ev[,AA])

sum_h<-mat%*%h_order

cor(sum_h,ev$peptide_len)

gas<-read.delim("/Users/hs/Downloads/amino_acid_gb_tsv (1).txt")

head(gas)

sum_b<-mat%*%gas$Gas.Phase.Basicity..kcal.mol.[1:20] 

sumb2<-sum_b/rowSums(mat)
#sumb2<-sum_h/rowSums(mat)

cor(sumb2,ev$peptide_len, use = "complete.obs")
cor(sumb2,ev$charge, use = "complete.obs")

dfm<-data.frame(sumb2,ev$peptide_len)
colnames(dfm)<-c("gas","len")
ggplot(dfm,aes(x=len, y=gas)) + 
  geom_point() +
  #geom_density_2d() + 
  geom_smooth() + 
  theme_tag() + 
  ylab("Mean Gas-phase basicity\n") + 
  xlab("\nPeptide length") + 
  geom_line(aes(y=195),linetype="dashed", color="red") +  
  geom_line(aes(y=205),linetype="dashed", color="red")  + 
  #ggtitle("235-4_24, Charge = 2",subtitle="Omitting C-term AA") 
ggtitle("Label-free, Charge = 4",subtitle="Omitting C-term AA") 



# dfm<-data.frame(sumb2,ev$expmass)
# colnames(dfm)<-c("gas","len")
# ggplot(dfm,aes(x=len, y=gas)) +
#   geom_point() +
#   #geom_density_2d() +
#   geom_smooth() +
#   theme_tag() +
#   ylab("Mean Gas-phase basicity\n") +
#   xlab("\nPeptide length")

# ggplot(dfm,aes(x=len, y=gas)) + 
#   stat_ecdf() + 
#   #geom_density_2d() + 
#   #geom_smooth() + 
#   theme_tag() + 
#   ylab("Mean Gas-phase basicity\n") + 
#   xlab("\nPeptide length")








# proportion of r/k by length bin - do 2 or 1 tags make it easier to ID long 



