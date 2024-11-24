setwd("~/Desktop/PhD_Project_related/COMPARING_WT_LT_BAY_AHRKO_DMSO_SAMP")



norm_counts=read.csv(file="Normalised_CPM_count.txt",header = T,stringsAsFactors = F,sep="\t")
bay_vs_dmso_wt=read.csv(file="Comparison_WT_vs_LTBAY_both_DMSO.txt",header = T,stringsAsFactors = F,sep="\t")
ahrko_vs_dmso_wt=read.csv(file="Comparison_WT_vs_AHRKO_both_DMSO.txt",header = T,stringsAsFactors = F,sep="\t")


adjp=0.01

bay_vs_dmso_wt=bay_vs_dmso_wt[complete.cases(bay_vs_dmso_wt),]
bay_vs_dmso_wt$baseMean_log=log2(bay_vs_dmso_wt$baseMean+1)


library(dplyr)
library(biomaRt)

ensembl_symbol=read.table(file="GRCm38_gene_symbol.txt",header = T,sep="\t",stringsAsFactors = F)

ensembl_symbol=subset(ensembl_symbol,ensembl_symbol$Gene.name!="")
ensembl_symbol=ensembl_symbol[,c(1,3)]
ensembl_symbol=ensembl_symbol[complete.cases(ensembl_symbol),]


#bay_vs_dmso_wt=left_join(bay_vs_dmso_wt,ensembl_symbol,by=c("Ensembl"="Gene.stable.ID"))
#colnames(bay_vs_dmso_wt)[9]="Symbol"

#all_combined=rbind(up,down,no_change)
library(ggrepel)
#ggplot(all_combined, aes(avg_exp_log, log2FoldChange)) +
#  geom_point(color = all_combined$Color) +
# theme_classic(base_size = 16)+
# geom_point(data=up_logfc_4, aes(x=avg_exp_log, y=log2FoldChange), colour="red", size=5)

bay_vs_dmso_wt$Neg_log_p_val=-log10(bay_vs_dmso_wt$padj)

ahrko_sig_up=subset(ahrko_vs_dmso_wt,ahrko_vs_dmso_wt$log2FoldChange>1 & ahrko_vs_dmso_wt$padj < 0.01)
ahrko_sig_dn=subset(ahrko_vs_dmso_wt,ahrko_vs_dmso_wt$log2FoldChange < -1 & ahrko_vs_dmso_wt$padj < 0.01)

bay_vs_dmso_wt_sig_up=subset(bay_vs_dmso_wt,bay_vs_dmso_wt$log2FoldChange > 1 & bay_vs_dmso_wt$padj < 0.01)
bay_vs_dmso_wt_sig_dn=subset(bay_vs_dmso_wt,bay_vs_dmso_wt$log2FoldChange < -1 & bay_vs_dmso_wt$padj < 0.01)


ahrko_and_bay_vs_dmso_wt_sig_up=subset(bay_vs_dmso_wt_sig_up,bay_vs_dmso_wt_sig_up$Ensembl %in% ahrko_sig_up$Ensembl)
ahrko_and_bay_vs_dmso_wt_sig_dn=subset(bay_vs_dmso_wt_sig_dn,bay_vs_dmso_wt_sig_dn$Ensembl %in% ahrko_sig_dn$Ensembl)


ahrko_and_bay_vs_dmso_wt_sig_up$bay_vs_dmso_wt_Direction="bay_vs_dmso_wt_Up"
ahrko_and_bay_vs_dmso_wt_sig_dn$bay_vs_dmso_wt_Direction="bay_vs_dmso_wt_Down"

bay_vs_dmso_wt_sig=rbind(ahrko_and_bay_vs_dmso_wt_sig_up,ahrko_and_bay_vs_dmso_wt_sig_dn)

#bay_vs_dmso_wt_sig=bay_vs_dmso_wt_sig[,c("Ensembl","bay_vs_dmso_wt_Direction")]


#bay_vs_dmso_wt=left_join(bay_vs_dmso_wt,bay_vs_dmso_wt_sig,by=c("Ensembl"))
#bay_vs_dmso_wt$bay_vs_dmso_wt_Direction[is.na(bay_vs_dmso_wt$bay_vs_dmso_wt_Direction)]="No sig change"


bay_vs_dmso_wt_sig_up=subset(bay_vs_dmso_wt_sig,bay_vs_dmso_wt_sig$log2FoldChange>0)
bay_vs_dmso_wt_sig_dn=subset(bay_vs_dmso_wt_sig,bay_vs_dmso_wt_sig$log2FoldChange < 0 )

bay_vs_dmso_wt_sig_up$bay_vs_dmso_wt_Direction="Up"
bay_vs_dmso_wt_sig_dn$bay_vs_dmso_wt_Direction="Down"


bay_vs_dmso_wt_sig_up=bay_vs_dmso_wt_sig_up[order(bay_vs_dmso_wt_sig_up$padj),]
bay_vs_dmso_wt_sig_dn=bay_vs_dmso_wt_sig_dn[order(bay_vs_dmso_wt_sig_dn$padj),]


#table(bay_vs_dmso_wt=bay_vs_dmso_wt_sig_up$bay_vs_dmso_wt_Direction,`bay_vs_dmso_wt Up`=bay_vs_dmso_wt_sig_up$bay_vs_dmso_wt_Direction)


p=ggplot(bay_vs_dmso_wt_sig, aes(baseMean_log, log2FoldChange)) +
  theme_classic(base_size = 16)+
  geom_point(data=bay_vs_dmso_wt_sig, aes(x=baseMean_log, y=log2FoldChange), colour="grey", size=2)
p1 <- p +  geom_point(data = bay_vs_dmso_wt_sig_up, aes(x=baseMean_log, y=log2FoldChange) ,size=3,color="red")
p2 <- p1 +  geom_point(data = bay_vs_dmso_wt_sig_dn, aes(x=baseMean_log, y=log2FoldChange) ,size=3,color="blue")
p2+ggtitle("MA plot for BAY vs DMSO significant genes overlapping with AHRKO vs DMSO in WT")+theme(plot.title = element_text(hjust = 0.5))+annotate(geom="text", x=15, y=-15, label="295 genes up after BAY vs DMSO, WT",color="red",size=6)+xlab("Log2(Mean+1)")+ylim(-23,23)+geom_text_repel(data=top_15_up_dn,aes(x=baseMean_log, y=log2FoldChange,label=Symbol),color="black",arrow=arrow(ends="last",type="open"))+xlim(0,19)



#boxplot(bay_vs_dmso_wt_sig_up_bay_vs_dmso_wt_up$log2FoldChange,bay_vs_dmso_wt_sig_up_bay_vs_dmso_wt_no_change$log2FoldChange,outline=FALSE)

#dn_tab=bay_vs_dmso_wt_sig_dn

#top_15_dn=dn_tab[1:15,]

#top_15_up_dn=rbind(top_15_up,top_15_dn)
inf_values=subset(bay_vs_dmso_wt_sig,bay_vs_dmso_wt_sig$Neg_log_p_val==Inf)
not_inf_values=subset(bay_vs_dmso_wt_sig,bay_vs_dmso_wt_sig$Neg_log_p_val!=Inf)

inf_values$Neg_log_p_val=300

bay_vs_dmso_wt_sig_new=rbind(inf_values,not_inf_values)

bay_vs_dmso_wt_sig_up=subset(bay_vs_dmso_wt_sig_new,bay_vs_dmso_wt_sig_new$log2FoldChange > 0)
bay_vs_dmso_wt_sig_dn=subset(bay_vs_dmso_wt_sig_new,bay_vs_dmso_wt_sig_new$log2FoldChange < 0)


top_15_up=bay_vs_dmso_wt_sig_up[1:5,]
top_15_dn=bay_vs_dmso_wt_sig_dn[1:5,]

top_15_up_dn=rbind(top_15_up,top_15_dn)

  
p=ggplot(bay_vs_dmso_wt_sig_new, aes(log2FoldChange, Neg_log_p_val)) +
  theme_classic(base_size = 16)+
  geom_point(data=bay_vs_dmso_wt_sig_new, aes(x=log2FoldChange, y=Neg_log_p_val), colour="grey", size=2)
p1 <- p +  geom_point(data = bay_vs_dmso_wt_sig_up, aes(x=log2FoldChange, y=Neg_log_p_val) ,size=3,color="red")
p2 <- p1 +  geom_point(data = bay_vs_dmso_wt_sig_dn, aes(x=log2FoldChange, y=Neg_log_p_val) ,size=3,color="blue")
p3=p2+ggtitle("Volcano plot for LT-BAY vs DMSO")+theme(plot.title = element_text(hjust = 0.5))+annotate(geom="text", x=10, y=320, label="295 genes up",color="red",size=4)+annotate(geom="text", x=10, y=298, label="419 genes down",color="blue",size=4)+xlab("Log2FoldChange")+xlim(-24,24)+ylim(0,350)+ylab("-log10(adj. p value)")+geom_text_repel(data=top_15_up_dn,aes(x=log2FoldChange, y=Neg_log_p_val,label=Symbol),color="black",arrow=arrow(ends="last",type="open"),fontface="bold",max.overlaps = 50)

tiff(file="Volcano plot LT-BAY vs DMSO with AHRKO vs DMSO overlapping genes.tiff",res=300,height = 1500,width = 3000)
p3
dev.off()



#rownames(norm_counts)=norm_counts$X
#norm_counts=norm_counts[,-c(1)]
avg_wt=apply(norm_counts[,c(1:4)],1,mean)
avg_bay=apply(norm_counts[,c(5:8)],1,mean)

avg_wt=as.data.frame(avg_wt)
avg_bay=as.data.frame(avg_bay)



rownames(bay_vs_dmso_wt)=bay_vs_dmso_wt$Ensembl

bay_vs_dmso_wt_with_exp=cbind(bay_vs_dmso_wt,avg_wt[rownames(bay_vs_dmso_wt),])
bay_vs_dmso_wt_with_exp=cbind(bay_vs_dmso_wt_with_exp,avg_bay[rownames(bay_vs_dmso_wt_with_exp),])
colnames(bay_vs_dmso_wt_with_exp)[c(12:13)]=c("DMSO","BAY")

bay_vs_dmso_wt_with_exp$DMSO_log=log2(bay_vs_dmso_wt_with_exp$DMSO+1)
bay_vs_dmso_wt_with_exp$bay_log=log2(bay_vs_dmso_wt_with_exp$bay+1)



bay_vs_dmso_wt_sig_up=subset(bay_vs_dmso_wt_with_exp,bay_vs_dmso_wt_with_exp$log2FoldChange>1 & bay_vs_dmso_wt_with_exp$padj < 0.01)
bay_vs_dmso_wt_sig_dn=subset(bay_vs_dmso_wt_with_exp,bay_vs_dmso_wt_with_exp$log2FoldChange < -1 & bay_vs_dmso_wt_with_exp$padj < 0.01)


bay_vs_dmso_wt_sig_up$bay_vs_dmso_wt_Direction=ifelse(bay_vs_dmso_wt_sig_up$log2FoldChange > 0,"Up")
bay_vs_dmso_wt_sig_dn$bay_vs_dmso_wt_Direction=ifelse(bay_vs_dmso_wt_sig_dn$log2FoldChange < 0,"Down")


bay_vs_dmso_wt_sig_up=bay_vs_dmso_wt_sig_up[order(bay_vs_dmso_wt_sig_up$padj),]
bay_vs_dmso_wt_sig_dn=bay_vs_dmso_wt_sig_dn[order(bay_vs_dmso_wt_sig_dn$padj),]

top_15_up=bay_vs_dmso_wt_sig_up[1:15,]
top_15_dn=bay_vs_dmso_wt_sig_dn[1:15,]

top_15_up_dn=rbind(top_15_up,top_15_dn)



p=ggplot(bay_vs_dmso_wt_with_exp, aes(DMSO_log, bay_log)) +
  theme_classic(base_size = 16)+
  geom_point(data=bay_vs_dmso_wt_with_exp, aes(x=DMSO_log, y=bay_log), colour="grey", size=2)
p1 <- p +  geom_point(data = bay_vs_dmso_wt_sig_up, aes(x=DMSO_log, y=bay_log) ,size=3,color="red")
p2 <- p1 +  geom_point(data = bay_vs_dmso_wt_sig_dn, aes(x=DMSO_log, y=bay_log) ,size=3,color="blue")
p2+ggtitle("Correlation plot for DMSO vs BAY-2416964 WT")+theme(plot.title = element_text(hjust = 0.5))+xlab("Log2(Mean+1) WT DMSO")+ylab("Log2(Mean+1) WT BAY-2416964")+xlim(0,15)+ylim(0,15)+geom_abline(slope=1, intercept=0,linetype="dotted")+annotate(geom="text", x=4, y=11, label="41 genes up after BAY-2416964 vs DMSO WT",color="red",size=6)+annotate(geom="text", x=4, y=12, label="15 genes down after BAY-2416964 vs DMSO WT",color="blue",size=6)+geom_text_repel(data=top_15_up_dn,aes(x=DMSO_log, y=bay_log,label=Symbol),color="black",arrow=arrow(ends="last",type="open"),fontface="bold")



