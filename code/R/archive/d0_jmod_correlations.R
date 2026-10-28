# library -- where da K go

ev<-read.csv("/Volumes/Lab/KMD/JModPreprint/PSMtag/JModSearches/2025-05-12_d0-mixD_200pg_20win_24nce_1_diann_tag6_Astral_MBR_anhy_May13_jmodUpdate130525_20ppm_3m_unmatchc_DECOYrev_libfrac0.5_RT_Dino_iso3.0_tag6__Anhy_rtfit/filtered_IDs.csv")
head(ev)
colnames(ev)

ev<-ev[grepl("K",ev$seq),]

#ggplot(ev,aes(x=,y=)) + geom_point()


#### correlation of results and features


nums <- unlist(lapply(ev, is.numeric), use.names = FALSE)

comp_cor<-ev[,nums]

comp_mat<-as.matrix(comp_cor)
#rownames(comp_mat)<-ev$Name
colnames(comp_mat)<-colnames(ev[,nums])

cor_mat<-cor(comp_mat, method = "pearson", use="pairwise.complete.obs")

kc<-colnames(cor_mat)[colSums(is.na(cor_mat)) < 54]

#kc<-c("Mass.x","Charge.x","Retention.time.x","Retention.length..FWHM.","Intensity","DP.time.difference.x","DP.score.x","LogP","NumAmide","NumHAcceptors","NumHDonors","NumAromaticCarbocycles","NumRotatableBonds")


d_rows <- dist(cor_mat[kc,kc])
hc_rows <- hclust(d_rows, method = "average")

dmat<-cor_mat[kc,kc]

dfm<-melt(dmat[hc_rows$order, hc_rows$order])

# dfm$Var1<-gsub(".x","",dfm$Var1)
# dfm$Var1<-gsub("\\."," ",dfm$Var1)
#
# dfm$Var2<-gsub(".x","",dfm$Var2)
# dfm$Var2<-gsub("\\."," ",dfm$Var2)


ggplot(dfm,aes(x=Var1,y=Var2,fill=value)) +
  geom_tile() +
  theme_tag() +
  xlab("") +
  ylab("") +
  scale_fill_gradient2(high = "red3", mid = "white", low = "blue3", name="Correlation" ) +
  #rremove("x.text") +
  #rremove("x.ticks")
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(axis.text = element_text(size=14)) +
  ggtitle("K peptides")


# ggsave("/Users/hs/Library/CloudStorage/GoogleDrive-hspecht@parallelsq.org/Shared drives/PTI Shared Drive/Projects/TeamTag/tag6_publication/figures/jmod_cors_K.pdf",
#        device="pdf",
#        width = image_in_w*3,
#        height = image_in_h*3,
#        units="in")



# library -- where da K go

ev<-read.csv("/Volumes/Lab/KMD/JModPreprint/PSMtag/JModSearches/2025-05-12_d0-mixD_200pg_20win_24nce_1_diann_tag6_Astral_MBR_anhy_May13_jmodUpdate130525_20ppm_3m_unmatchc_DECOYrev_libfrac0.5_RT_Dino_iso3.0_tag6__Anhy_rtfit/filtered_IDs.csv")
head(ev)
colnames(ev)

ev<-ev[!grepl("K",ev$seq),]

#ggplot(ev,aes(x=,y=)) + geom_point()


#### correlation of results and features


nums <- unlist(lapply(ev, is.numeric), use.names = FALSE)

comp_cor<-ev[,nums]

comp_mat<-as.matrix(comp_cor)
#rownames(comp_mat)<-ev$Name
colnames(comp_mat)<-colnames(ev[,nums])

cor_mat<-cor(comp_mat, method = "pearson", use="pairwise.complete.obs")

kc<-colnames(cor_mat)[colSums(is.na(cor_mat)) < 54]

#kc<-c("Mass.x","Charge.x","Retention.time.x","Retention.length..FWHM.","Intensity","DP.time.difference.x","DP.score.x","LogP","NumAmide","NumHAcceptors","NumHDonors","NumAromaticCarbocycles","NumRotatableBonds")


d_rows <- dist(cor_mat[kc,kc])
hc_rows <- hclust(d_rows, method = "average")

dmat<-cor_mat[kc,kc]

dfm2<-melt(dmat[hc_rows$order, hc_rows$order])

# dfm$Var1<-gsub(".x","",dfm$Var1)
# dfm$Var1<-gsub("\\."," ",dfm$Var1)
#
# dfm$Var2<-gsub(".x","",dfm$Var2)
# dfm$Var2<-gsub("\\."," ",dfm$Var2)


ggplot(dfm2,aes(x=Var1,y=Var2,fill=value)) +
  geom_tile() +
  theme_tag() +
  xlab("") +
  ylab("") +
  scale_fill_gradient2(high = "red3", mid = "white", low = "blue3", name="Correlation" ) +
  #rremove("x.text") +
  #rremove("x.ticks")
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(axis.text = element_text(size=14)) +
  ggtitle("R peptides")


#dfm$pep<-"K"; dfm2$pep<-"R"
#dfm3<-rbind(dfm,dfm2)

dfm$comp<-paste0(dfm$Var1,"-vs-",dfm$Var2)
dfm2$comp<-paste0(dfm2$Var1,"-vs-",dfm2$Var2)

dfx<-merge(dfm[,c("comp","value")], dfm2[,c("comp","value")], by="comp")

dfx

dfx$diff<-dfx$value.x - dfx$value.y

dfx<-dfx[order(dfx$diff),]

(dfx[seq(1,30,2),])


dfx<-dfx[order(abs(dfx$diff),decreasing = T),]

(dfx[seq(1,40,2),])
