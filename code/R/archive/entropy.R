
ev<-read.csv("/Volumes/Lab/KMD/JModPreprint/PSMtag/JModSearches/2025-05-12_d0-mixD_200pg_20win_24nce_1_diann_tag6_Astral_MBR_anhy_May13_jmodUpdate130525_20ppm_3m_unmatchc_DECOYrev_libfrac0.5_RT_Dino_iso3.0_tag6__Anhy_rtfit/filtered_IDs.csv")
head(ev)
colnames(ev)

ev$obs_ent<-NA
ev$lib_ent<-NA

for(i in 1:nrow(ev)){
  
  vt<-as.numeric(unlist(str_split(ev$obs_int[i],";")))
  
  vt<-vt/sum(vt)
  
  ev$obs_ent[i]<- -sum(vt*log2(vt))
  
  vt<-as.numeric(unlist(str_split(ev$frag_int[i],";")))
  
  vt<-vt/sum(vt)
  
  ev$lib_ent[i]<- -sum(vt*log2(vt))
  
  
}

ev$aa<-"R"
ev$aa[grepl("K",ev$seq)]<-"K"

ggplot(ev,aes(x=obs_ent,color=aa)) + 
  geom_density() + 
  theme_tag() + 
  xlab("\n Entropy")

ev2<-ev[ev$pep_len>9 & ev$pep_len<15,]
p1<-ggplot(ev2,aes(x=obs_ent,color=aa)) + 
  geom_density() + 
  theme_tag() + 
  xlab("\n Entropy")

facet(p1, facet.by = "pep_len")

ggplot(ev,aes(x=obs_ent/lib_ent,color=aa)) + 
  geom_density() + 
  theme_tag() + 
  xlab("\n Entropy")




d1<-read.delim("/Users/hs/Library/CloudStorage/GoogleDrive-hspecht@parallelsq.org/Shared drives/PTI Shared Drive/Projects/TeamTag/tag6_publication/processed_data/astral_CE_DDA_final/T6_24/results.sage.tsv")
d2<-read.delim("/Users/hs/Library/CloudStorage/GoogleDrive-hspecht@parallelsq.org/Shared drives/PTI Shared Drive/Projects/TeamTag/tag6_publication/processed_data/astral_CE_DDA_final/T6_24/matched_fragments.sage.tsv")

head(d1)

d1$peptide<-gsub("[+308.1161]","",d1$peptide,fixed=T)
d1$peptide<-gsub("-","",d1$peptide,fixed=T)
d1$peptide<-gsub("[+57.021465]","",d1$peptide,fixed=T)

d1$obs_ent<-NA
for(i in 1:nrow(d1)){
  
  vt<-as.numeric(d2$fragment_intensity[d2$psm_id%in%d1$psm_id[i]])
  it<-(d2$fragment_type[d2$psm_id%in%d1$psm_id[i]])
  
  vt<-vt[it%in%"b"]
  
  vt<-vt/sum(vt)
  
  d1$obs_ent[i]<- -sum(vt*log2(vt))
  #d1$obs_ent[i]<- length(vt)
  
}

d1$aa<-"R"
d1$aa[grepl("K",d1$peptide)]<-"K"

dia_seq<-ev$seq

dia_seq<-gsub("(tag6-0)","",dia_seq,fixed=T)
dia_seq<-gsub("(UniMod:4)","",dia_seq,fixed=T)

d1$obs_dia<-d1$peptide%in%dia_seq

p1<-ggplot(d1,aes(x=obs_ent,color=aa)) + 
  geom_density() + 
  theme_tag() + 
  xlab("\n Entropy")

facet(p1, facet.by = "obs_dia")

p1<-ggplot(d1,aes(x=obs_ent,color=obs_dia)) + 
  geom_density() + 
  theme_tag() + 
  xlab("\n Entropy")

facet(p1, facet.by = "aa")

