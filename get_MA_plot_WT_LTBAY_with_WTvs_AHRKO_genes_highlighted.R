setwd("~/Desktop/PhD_Project_related/COMPARING_WT_LT_BAY_AHRKO_DMSO_SAMP")



norm_counts=read.csv(file="Normalised_CPM_count.txt",header = T,stringsAsFactors = F,sep="\t")
ltbay=read.csv(file="Comparison_WT_vs_LTBAY_both_DMSO.txt",header = T,stringsAsFactors = F,sep="\t")
ahrko=read.csv(file="Comparison_WT_vs_AHRKO_both_DMSO.txt",header = T,stringsAsFactors = F,sep="\t")

adjp=0.01

ahrko=ahrko[complete.cases(ahrko),]
ahrko$baseMean_log=log2(ahrko$baseMean+1)

ltbay=ltbay[complete.cases(ltbay),]
ltbay$baseMean_log=log2(ltbay$baseMean+1)




library(dplyr)
library(biomaRt)

ensembl_symbol=read.table(file="GRCm38_gene_symbol.txt",header = T,sep="\t",stringsAsFactors = F)

ensembl_symbol=subset(ensembl_symbol,ensembl_symbol$Gene.name!="")
ensembl_symbol=ensembl_symbol[,c(1,3)]
ensembl_symbol=ensembl_symbol[complete.cases(ensembl_symbol),]


ahrko=left_join(ahrko,ensembl_symbol,by=c("Ensembl"="Gene.stable.ID"))
ltbay=left_join(ltbay,ensembl_symbol,by=c("Ensembl"="Gene.stable.ID"))

colnames(ahrko)[9]="Symbol"
colnames(ltbay)[9]="Symbol"

#all_combined=rbind(up,down,no_change)
library(ggrepel)
#ggplot(all_combined, aes(avg_exp_log, log2FoldChange)) +
#  geom_point(color = all_combined$Color) +
# theme_classic(base_size = 16)+
# geom_point(data=up_logfc_4, aes(x=avg_exp_log, y=log2FoldChange), colour="red", size=5)

ahrko$Neg_log_p_val=-log10(ahrko$padj)
ltbay$Neg_log_p_val=-log10(ltbay$padj)


ahrko_sig_up=subset(ahrko,ahrko$log2FoldChange>1 & ahrko$padj < 0.01)
ahrko_sig_dn=subset(ahrko,ahrko$log2FoldChange < -1 & ahrko$padj < 0.01)




ahrko_sig_up$ahrko_Direction="AHRko_Up"
ahrko_sig_dn$ahrko_Direction="AHRko_Down"

ahrko_sig=rbind(ahrko_sig_up,ahrko_sig_dn)

ahrko_sig=ahrko_sig[,c("Ensembl","ahrko_Direction")]


ahrko=left_join(ahrko,ahrko_sig,by=c("Ensembl"))
ahrko$ahrko_Direction[is.na(ahrko$ahrko_Direction)]="No sig change"




p=ggplot(ahrko, aes(baseMean_log, log2FoldChange)) +
  theme_classic(base_size = 16)+
  geom_point(data=ahrko, aes(x=baseMean_log, y=log2FoldChange), colour="grey", size=2)
p1 <- p +  geom_point(data = ahrko_sig_up, aes(x=baseMean_log, y=log2FoldChange) ,size=3,color="red")
p2 <- p1 +  geom_point(data = ahrko_sig_dn, aes(x=baseMean_log, y=log2FoldChange) ,size=3,color="blue")
p2+ggtitle("MA plot for AHRKO DMSO vs WT DMSO")+theme(plot.title = element_text(hjust = 0.5))+annotate(geom="text", x=15, y=-15, label="701 genes up in AHRKO vs WT",color="red",size=5)+annotate(geom="text", x=15, y=-16.5, label="1159 genes down in AHRKO vs WT",color="blue",size=5)+xlab("Log2(Mean+1)")+ylim(-18,10)+geom_text_repel(data=top_15_up_dn,aes(x=baseMean_log, y=log2FoldChange,label=Symbol),color="black",arrow=arrow(ends="last",type="open"),max.overlaps = 50)+xlim(0,18)


ahrko_inf=subset(ahrko,ahrko$Neg_log_p_val==Inf)
ahrko_not_inf=subset(ahrko,ahrko$Neg_log_p_val!=Inf)

ahrko_inf$Neg_log_p_val=400

ahrko_new=rbind(ahrko_inf,ahrko_not_inf)


ahrko_sig_up=subset(ahrko_new,ahrko_new$log2FoldChange>1 & ahrko_new$padj < 0.01)
ahrko_sig_dn=subset(ahrko_new,ahrko_new$log2FoldChange < -1 & ahrko_new$padj < 0.01)

ahrko_sig_up$ahrko_Direction=ifelse(ahrko_sig_up$log2FoldChange > 0,"Up")
ahrko_sig_dn$ahrko_Direction=ifelse(ahrko_sig_dn$log2FoldChange < 0,"Down")


ahrko_sig_up=ahrko_sig_up[order(ahrko_sig_up$padj),]
ahrko_sig_dn=ahrko_sig_dn[order(ahrko_sig_dn$padj),]

top_15_up=ahrko_sig_up[1:15,]
top_15_dn=ahrko_sig_dn[1:15,]

top_15_up_dn=rbind(top_15_up,top_15_dn)


p=ggplot(ahrko_new, aes(log2FoldChange, Neg_log_p_val)) +
  theme_classic(base_size = 16)+
  geom_point(data=ahrko_new, aes(x=log2FoldChange, y=Neg_log_p_val), colour="grey", size=2)
p1 <- p +  geom_point(data = ahrko_sig_up, aes(x=log2FoldChange, y=Neg_log_p_val) ,size=3,color="red")
p2 <- p1 +  geom_point(data = ahrko_sig_dn, aes(x=log2FoldChange, y=Neg_log_p_val) ,size=3,color="blue")
p2+ggtitle("Volcano plot for AHRKO DMSO vs WT DMSO")+theme(plot.title = element_text(hjust = 0.5))+annotate(geom="text", x=15, y=350, label="701 genes up in AHRKO vs WT",color="red",size=5)+annotate(geom="text", x=15, y=330, label="1159 genes down in AHRKO vs WT",color="blue",size=5)+xlab("Log2(Mean+1)")+ylim(0,400)+geom_text_repel(data=top_15_up_dn,aes(x=log2FoldChange, y=Neg_log_p_val,label=Symbol),color="black",arrow=arrow(ends="last",type="open"),max.overlaps = 50)




intersect(ahrko_sig_up$Ensembl,ltbay_sig_up$Ensembl)
intersect(ahrko_sig_dn$Ensembl,ltbay_sig_dn$Ensembl)


ltbay_sig_up=subset(ltbay,ltbay$log2FoldChange>1 & ltbay$padj < 0.01)
ltbay_sig_dn=subset(ltbay,ltbay$log2FoldChange < -1 & ltbay$padj < 0.01)




ltbay_and_ahrko_up=subset(ltbay_sig_up,ltbay_sig_up$Ensembl %in% ahrko_sig_up$Ensembl)
ltbay_and_ahrko_dn=subset(ltbay_sig_dn,ltbay_sig_dn$Ensembl %in% ahrko_sig_dn$Ensembl)

`%ni%`=Negate(`%in%`)

ltbay_and_not_ahrko_up=subset(ltbay_sig_up,ltbay_sig_up$Ensembl %ni% ahrko_sig_up$Ensembl)
ltbay_and_not_ahrko_dn=subset(ltbay_sig_dn,ltbay_sig_dn$Ensembl %ni% ahrko_sig_dn$Ensembl)




#table(ahrko=ahrko_sig_up$ahrko_Direction,`ahrko Up`=ahrko_sig_up$ahrko_Direction)


p=ggplot(ahrko, aes(baseMean_log, log2FoldChange)) +
  theme_classic(base_size = 16)+
  geom_point(data=ahrko, aes(x=baseMean_log, y=log2FoldChange), colour="grey", size=2)
p1 <- p +  geom_point(data = ltbay_and_ahrko_up, aes(x=baseMean_log, y=log2FoldChange) ,size=3,color="purple")
p2 <- p1 +  geom_point(data = ltbay_and_ahrko_dn, aes(x=baseMean_log, y=log2FoldChange) ,size=3,color="green")
p2+ggtitle("MA plot for WT DMSO vs LT-BAY DMSO with WT vs AHRKo DMSO genes highlighted")+theme(plot.title = element_text(hjust = 0.5))+annotate(geom="text", x=15, y=-15, label="295 genes up in LT-BAY vs WT and AHRKO vs WT",color="purple",size=3)+annotate(geom="text", x=15, y=-16.5, label="419 genes down in LT-BAY vs WT and AHRKO vs WT",color="green",size=3)+xlab("Log2(Mean+1)")+ylim(-18,10)+geom_text_repel(data=top_15_up_dn,aes(x=baseMean_log, y=log2FoldChange,label=Symbol),color="black",arrow=arrow(ends="last",type="open"))+xlim(0,18)



ltbay_and_ahrko_up=ltbay_and_ahrko_up[order(-ltbay_and_ahrko_up$Neg_log_p_val),]
ltbay_and_ahrko_dn=ltbay_and_ahrko_dn[order(-ltbay_and_ahrko_dn$Neg_log_p_val),]


top_5_up=ltbay_and_ahrko_up[1:5,]
top_5_dn=ltbay_and_ahrko_dn[1:5,]

top_5_up_dn=rbind(top_5_up,top_5_dn)



goi=c("Ahr", "Cyp1a1", "Cyp1b1", "Tiparp", "Ahrr")

goi_df=subset(ltbay, ltbay$Symbol %in% goi)

top_5_up_dn_goi=rbind(top_5_up_dn,goi_df)


library(org.Mm.eg.db)

ltbay$Entrez <- mapIds(org.Mm.eg.db, ltbay$Ensembl,keytype="ENSEMBL", column="ENTREZID")
ltbay$Symbol <- mapIds(org.Mm.eg.db, ltbay$Entrez,keytype="ENTREZID", column="SYMBOL")
ltbay$GeneName <- mapIds(org.Mm.eg.db, ltbay$Ensembl,keytype="ENSEMBL", column="GENENAME")

ltbay$Entrez=as.character(ltbay$Entrez)
ltbay$Symbol=as.character(ltbay$Symbol)
ltbay$GeneName=as.character(ltbay$GeneName)


ltbay_noGm_riken=ltbay[- grep("RIKEN",ltbay$GeneName),]
ltbay_noGm_riken=ltbay_noGm_riken[- grep("Riken",ltbay_noGm_riken$GeneName),]

ltbay_noGm_riken=ltbay_noGm_riken[- grep("predicted",ltbay_noGm_riken$GeneName),]

ltBay_inf=subset(ltbay_noGm_riken,ltbay_noGm_riken$Neg_log_p_val==Inf)
ltBay_not_inf=subset(ltbay_noGm_riken,ltbay_noGm_riken$Neg_log_p_val!=Inf)

ltBay_inf$Neg_log_p_val=300

ltBay_new=rbind(ltBay_inf,ltBay_not_inf)








ahrko$Entrez <- mapIds(org.Mm.eg.db, ahrko$Ensembl,keytype="ENSEMBL", column="ENTREZID")
ahrko$Symbol <- mapIds(org.Mm.eg.db, ahrko$Entrez,keytype="ENTREZID", column="SYMBOL")
ahrko$GeneName <- mapIds(org.Mm.eg.db, ahrko$Ensembl,keytype="ENSEMBL", column="GENENAME")

ahrko$Entrez=as.character(ahrko$Entrez)
ahrko$Symbol=as.character(ahrko$Symbol)
ahrko$GeneName=as.character(ahrko$GeneName)


ahrko_noGm_riken=ahrko[- grep("RIKEN",ahrko$GeneName),]
ahrko_noGm_riken=ahrko_noGm_riken[- grep("Riken",ahrko_noGm_riken$GeneName),]

ahrko_noGm_riken=ahrko_noGm_riken[- grep("predicted",ahrko_noGm_riken$GeneName),]


ltbay_sig_up=subset(ltBay_new,ltBay_new$log2FoldChange>1 & ltBay_new$padj < 0.01)
ltbay_sig_dn=subset(ltBay_new,ltBay_new$log2FoldChange < -1 & ltBay_new$padj < 0.01)



ahrko_sig_up=subset(ahrko_noGm_riken,ahrko_noGm_riken$log2FoldChange>1 & ahrko_noGm_riken$padj < 0.01)
ahrko_sig_dn=subset(ahrko_noGm_riken,ahrko_noGm_riken$log2FoldChange < -1 & ahrko_noGm_riken$padj < 0.01)







goi=c("Ahr", "Cyp1a1", "Cyp1b1", "Tiparp", "Ahrr")

goi_df=subset(ltbay_noGm_riken, ltbay_noGm_riken$Symbol %in% goi)

top_5_up_dn_goi=rbind(top_5_up_dn,goi_df)


ahrko_sig_up=subset(ahrko_new,ahrko_new$log2FoldChange>1 & ahrko_new$padj < 0.01)
ahrko_sig_dn=subset(ahrko_new,ahrko_new$log2FoldChange < -1 & ahrko_new$padj < 0.01)




ltbay_and_ahrko_up=subset(ltbay_sig_up,ltbay_sig_up$Ensembl %in% ahrko_sig_up$Ensembl)
ltbay_and_ahrko_dn=subset(ltbay_sig_dn,ltbay_sig_dn$Ensembl %in% ahrko_sig_dn$Ensembl)


ltbay_and_ahrko_up=ltbay_and_ahrko_up[order(-ltbay_and_ahrko_up$Neg_log_p_val),]
ltbay_and_ahrko_dn=ltbay_and_ahrko_dn[order(-ltbay_and_ahrko_dn$Neg_log_p_val),]


top_5_up=ltbay_and_ahrko_up[1:5,]
top_5_dn=ltbay_and_ahrko_dn[1:5,]

top_5_up_dn=rbind(top_5_up,top_5_dn)


top_5_up_dn_goi=rbind(top_5_up_dn,goi_df)


p=ggplot(ltbay_noGm_riken, aes(log2FoldChange, Neg_log_p_val)) +
  theme_classic(base_size = 16)+
  geom_point(data=ltbay_noGm_riken, aes(x=log2FoldChange, y=Neg_log_p_val), colour="grey", size=2)
p1 <- p +  geom_point(data = ltbay_and_ahrko_up, aes(x=log2FoldChange, y=Neg_log_p_val) ,size=3,color="grey")
p2 <- p1 +  geom_point(data = ltbay_and_ahrko_dn, aes(x=log2FoldChange, y=Neg_log_p_val) ,size=3,color="grey")
p3 <- p2 +  geom_point(data = top_5_up_dn_goi, aes(x=log2FoldChange, y=Neg_log_p_val) ,size=3,color="green")

p4=p3+ggtitle("Volcano plot for WT DMSO vs LT-BAY DMSO with only selected genes highlighted")+theme(plot.title = element_text(hjust = 0.5,size=10))+xlab("Log2 Fold Change")+xlim(-27,24) + ylab("-Log10(adj. p value)") + ylim(0,320) + geom_vline(xintercept = c(-1,1), linetype="dotted",  color = "black", size=1) + geom_hline(yintercept = c(2), linetype="dotted",color = "black", size=1)+geom_text_repel(data=top_5_up_dn_goi,aes(x=log2FoldChange, y=Neg_log_p_val,label=Symbol),color="purple",arrow=arrow(ends="last",type="open",length = unit(0.05,"inches")))
   p4                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             

tiff(file="WT vs LT-BAY Volcano plot with selected genes highlighted.tiff",res=300,height = 2000,width = 3000)
p4
dev.off()




p=ggplot(ltbay_noGm_riken, aes(log2FoldChange, Neg_log_p_val)) +
  theme_classic(base_size = 16)+
  geom_point(data=ltbay_noGm_riken, aes(x=log2FoldChange, y=Neg_log_p_val), colour="grey", size=2)
p1 <- p +  geom_point(data = ltbay_and_ahrko_up, aes(x=log2FoldChange, y=Neg_log_p_val) ,size=3,color="grey")
p2 <- p1 +  geom_point(data = ltbay_and_ahrko_dn, aes(x=log2FoldChange, y=Neg_log_p_val) ,size=3,color="grey")
p3 <- p2 +  geom_point(data = top_5_up_dn_goi, aes(x=log2FoldChange, y=Neg_log_p_val) ,size=3,color="green")
p4=p3+ggtitle("Volcano plot for WT DMSO vs LT-BAY DMSO with selected genes highlighted")+theme(plot.title = element_text(hjust = 0.5))+xlab("Log2 Fold Change")+ylim(-15,10)+geom_text_repel(data=top_5_up_dn_goi,aes(x=log2FoldChange, y=Neg_log_p_val,label=Symbol),color="purple",arrow=arrow(ends="last",type="open",length = unit(0.05,"inches")),max.overlaps = 10)+xlim(-27,24) + ylab("-Log10(adj. p value)") + ylim(0,320) + geom_vline(xintercept = c(-1,1), linetype="dotted",color = "black", size=1) + geom_hline(yintercept = c(2), linetype="dotted",color = "black", size=1) 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   
p4
tiff(file="Volcano plot with only genes of interest highlighted.tiff",res=300,height = 1500,width = 3000)
p4
dev.off()




