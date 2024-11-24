setwd("~/Desktop/PhD_Project_related/COMPARING_WT_LT_BAY_AHRKO_DMSO_SAMP")



norm_counts=read.csv(file="Normalised_CPM_count.txt",header = T,stringsAsFactors = F,sep="\t")
ltbay_vs_wt=read.csv(file="Comparison_WT_vs_LTBAY_both_DMSO.txt",header = T,stringsAsFactors = F,sep="\t")
adjp=0.01

ltbay_vs_wt=ltbay_vs_wt[complete.cases(ltbay_vs_wt),]
ltbay_vs_wt$baseMean_log=log2(ltbay_vs_wt$baseMean+1)


library(dplyr)
library(biomaRt)

ensembl_symbol=read.table(file="GRCm38_gene_symbol.txt",header = T,sep="\t",stringsAsFactors = F)

ensembl_symbol=subset(ensembl_symbol,ensembl_symbol$Gene.name!="")
ensembl_symbol=ensembl_symbol[,c(1,3)]
ensembl_symbol=ensembl_symbol[complete.cases(ensembl_symbol),]


ltbay_vs_wt=left_join(ltbay_vs_wt,ensembl_symbol,by=c("Ensembl"="Gene.stable.ID"))
colnames(ltbay_vs_wt)[9]="Symbol"

#all_combined=rbind(up,down,no_change)
library(ggrepel)
#ggplot(all_combined, aes(avg_exp_log, log2FoldChange)) +
#  geom_point(color = all_combined$Color) +
# theme_classic(base_size = 16)+
# geom_point(data=up_logfc_4, aes(x=avg_exp_log, y=log2FoldChange), colour="red", size=5)

ltbay_vs_wt$Neg_log_p_val=-log10(ltbay_vs_wt$padj)

ltbay_vs_wt_sig_up=subset(ltbay_vs_wt,ltbay_vs_wt$log2FoldChange>1 & ltbay_vs_wt$padj < 0.01)
ltbay_vs_wt_sig_dn=subset(ltbay_vs_wt,ltbay_vs_wt$log2FoldChange < -1 & ltbay_vs_wt$padj < 0.01)

ltbay_vs_wt_sig_up$ltbay_vs_wt_Direction="ltbay_vs_wt_Up"
ltbay_vs_wt_sig_dn$ltbay_vs_wt_Direction="ltbay_vs_wt_Down"

ltbay_vs_wt_sig=rbind(ltbay_vs_wt_sig_up,ltbay_vs_wt_sig_dn)

ltbay_vs_wt_sig=ltbay_vs_wt_sig[,c("Ensembl","ltbay_vs_wt_Direction")]


ltbay_vs_wt=left_join(ltbay_vs_wt,ltbay_vs_wt_sig,by=c("Ensembl"))
ltbay_vs_wt$ltbay_vs_wt_Direction[is.na(ltbay_vs_wt$ltbay_vs_wt_Direction)]="No sig change"




#table(ltbay_vs_wt=ltbay_vs_wt_sig_up$ltbay_vs_wt_Direction,`ltbay_vs_wt Up`=ltbay_vs_wt_sig_up$ltbay_vs_wt_Direction)


p=ggplot(ltbay_vs_wt, aes(baseMean_log, log2FoldChange)) +
  theme_classic(base_size = 16)+
  geom_point(data=ltbay_vs_wt, aes(x=baseMean_log, y=log2FoldChange), colour="grey", size=2)
p1 <- p +  geom_point(data = ltbay_vs_wt_sig_up, aes(x=baseMean_log, y=log2FoldChange) ,size=3,color="red")
p2 <- p1 +  geom_point(data = ltbay_vs_wt_sig_dn, aes(x=baseMean_log, y=log2FoldChange) ,size=3,color="blue")
p2+ggtitle("MA plot for LT-BAY DMSO vs WT DMSO")+theme(plot.title = element_text(hjust = 0.5))+annotate(geom="text", x=15, y=-15, label="895 genes up after LT-BAY DMSO vs WT DMSO",color="red",size=3)+annotate(geom="text", x=15, y=-17, label="1187 genes down after LT-BAY DMSO vs WT DMSO",color="blue",size=3)+xlab("Log2(Mean+1)")+ylim(-27,27)+geom_text_repel(data=top_15_up_dn,aes(x=baseMean_log, y=log2FoldChange,label=Symbol),color="black",arrow=arrow(ends="last",type="open"))+xlim(0,17)

ltbay_vs_wt_inf=subset(ltbay_vs_wt,ltbay_vs_wt$Neg_log_p_val==Inf)
ltbay_vs_wt_not_inf=subset(ltbay_vs_wt,ltbay_vs_wt$Neg_log_p_val!=Inf)

ltbay_vs_wt_inf$Neg_log_p_val=300

ltbay_vs_wt_new=rbind(ltbay_vs_wt_inf,ltbay_vs_wt_not_inf)

ltbay_vs_wt_sig_up=subset(ltbay_vs_wt_new,ltbay_vs_wt_new$log2FoldChange>1 & ltbay_vs_wt_new$padj < 0.01)
ltbay_vs_wt_sig_dn=subset(ltbay_vs_wt_new,ltbay_vs_wt_new$log2FoldChange < -1 & ltbay_vs_wt_new$padj < 0.01)

ltbay_vs_wt_sig_up$ltbay_vs_wt_Direction=ifelse(ltbay_vs_wt_sig_up$log2FoldChange > 0,"Up")
ltbay_vs_wt_sig_dn$ltbay_vs_wt_Direction=ifelse(ltbay_vs_wt_sig_dn$log2FoldChange < 0,"Down")


ltbay_vs_wt_sig_up=ltbay_vs_wt_sig_up[order(ltbay_vs_wt_sig_up$padj),]
ltbay_vs_wt_sig_dn=ltbay_vs_wt_sig_dn[order(ltbay_vs_wt_sig_dn$padj),]

top_15_up=ltbay_vs_wt_sig_up[1:15,]
top_15_dn=ltbay_vs_wt_sig_dn[1:15,]

top_15_up_dn=rbind(top_15_up,top_15_dn)



p=ggplot(ltbay_vs_wt_new, aes(log2FoldChange, Neg_log_p_val)) +
  theme_classic(base_size = 16)+
  geom_point(data=ltbay_vs_wt_new, aes(x=log2FoldChange, y=Neg_log_p_val), colour="grey", size=2)
p1 <- p +  geom_point(data = ltbay_vs_wt_sig_up, aes(x=log2FoldChange, y=Neg_log_p_val) ,size=3,color="red")
p2 <- p1 +  geom_point(data = ltbay_vs_wt_sig_dn, aes(x=log2FoldChange, y=Neg_log_p_val) ,size=3,color="blue")
p2+ggtitle("Volcano plot for LT-BAY DMSO vs WT DMSO")+theme(plot.title = element_text(hjust = 0.5))+annotate(geom="text", x=15, y=290, label="895 genes up after LT-BAY DMSO vs WT DMSO",color="red",size=4)+annotate(geom="text", x=15, y=280, label="1187 genes down after LT-BAY DMSO vs WT DMSO",color="blue",size=4)+xlab("Log2(Mean+1)")+ylim(0,350)+geom_text_repel(data=top_15_up_dn,aes(x=log2FoldChange, y=Neg_log_p_val,label=Symbol),color="black",arrow=arrow(ends="last",type="open"))+xlim(-27,27)+ylab("-Log10(adj. pval)")



#boxplot(ltbay_vs_wt_sig_up_ltbay_vs_wt_up$log2FoldChange,ltbay_vs_wt_sig_up_ltbay_vs_wt_no_change$log2FoldChange,outline=FALSE)

#dn_tab=ltbay_vs_wt_sig_dn

#top_15_dn=dn_tab[1:15,]

#top_15_up_dn=rbind(top_15_up,top_15_dn)

p=ggplot(ltbay_vs_wt, aes(log2FoldChange, Neg_log_p_val)) +
  theme_classic(base_size = 16)+
  geom_point(data=ltbay_vs_wt, aes(x=log2FoldChange, y=Neg_log_p_val), colour="grey", size=2)
p1 <- p +  geom_point(data = ltbay_vs_wt_sig_up, aes(x=log2FoldChange, y=Neg_log_p_val) ,size=3,color="red")
p2 <- p1 +  geom_point(data = ltbay_vs_wt_sig_dn, aes(x=log2FoldChange, y=Neg_log_p_val) ,size=3,color="blue")
p2+ggtitle("Volcano plot for DMSO vs BAY-2416964 WT")+theme(plot.title = element_text(hjust = 0.5))+annotate(geom="text", x=-2, y=55, label="41 genes up after BAY-2416964 vs DMSO WT",color="red",size=6)+annotate(geom="text", x=-2, y=53, label="15 genes down after BAY-2416964 vs DMSO WT",color="blue",size=6)+xlab("Log2FoldChange")+xlim(-5,10)+ylim(0,2)+ylab("-log10(adj. p value)")


#rownames(norm_counts)=norm_counts$X
#norm_counts=norm_counts[,-c(1)]
avg_wt=apply(norm_counts[,c(1:4)],1,mean)
avg_bay=apply(norm_counts[,c(5:8)],1,mean)

avg_wt=as.data.frame(avg_wt)
avg_bay=as.data.frame(avg_bay)



rownames(ltbay_vs_wt)=ltbay_vs_wt$Ensembl

ltbay_vs_wt_with_exp=cbind(ltbay_vs_wt,avg_wt[rownames(ltbay_vs_wt),])
ltbay_vs_wt_with_exp=cbind(ltbay_vs_wt_with_exp,avg_bay[rownames(ltbay_vs_wt_with_exp),])
colnames(ltbay_vs_wt_with_exp)[c(12:13)]=c("DMSO","BAY")

ltbay_vs_wt_with_exp$DMSO_log=log2(ltbay_vs_wt_with_exp$DMSO+1)
ltbay_vs_wt_with_exp$bay_log=log2(ltbay_vs_wt_with_exp$bay+1)



ltbay_vs_wt_sig_up=subset(ltbay_vs_wt_with_exp,ltbay_vs_wt_with_exp$log2FoldChange>1 & ltbay_vs_wt_with_exp$padj < 0.01)
ltbay_vs_wt_sig_dn=subset(ltbay_vs_wt_with_exp,ltbay_vs_wt_with_exp$log2FoldChange < -1 & ltbay_vs_wt_with_exp$padj < 0.01)


ltbay_vs_wt_sig_up$ltbay_vs_wt_Direction=ifelse(ltbay_vs_wt_sig_up$log2FoldChange > 0,"Up")
ltbay_vs_wt_sig_dn$ltbay_vs_wt_Direction=ifelse(ltbay_vs_wt_sig_dn$log2FoldChange < 0,"Down")


ltbay_vs_wt_sig_up=ltbay_vs_wt_sig_up[order(ltbay_vs_wt_sig_up$padj),]
ltbay_vs_wt_sig_dn=ltbay_vs_wt_sig_dn[order(ltbay_vs_wt_sig_dn$padj),]

top_15_up=ltbay_vs_wt_sig_up[1:15,]
top_15_dn=ltbay_vs_wt_sig_dn[1:15,]

top_15_up_dn=rbind(top_15_up,top_15_dn)



p=ggplot(ltbay_vs_wt_with_exp, aes(DMSO_log, bay_log)) +
  theme_classic(base_size = 16)+
  geom_point(data=ltbay_vs_wt_with_exp, aes(x=DMSO_log, y=bay_log), colour="grey", size=2)
p1 <- p +  geom_point(data = ltbay_vs_wt_sig_up, aes(x=DMSO_log, y=bay_log) ,size=3,color="red")
p2 <- p1 +  geom_point(data = ltbay_vs_wt_sig_dn, aes(x=DMSO_log, y=bay_log) ,size=3,color="blue")
p2+ggtitle("Correlation plot for DMSO vs BAY-2416964 WT")+theme(plot.title = element_text(hjust = 0.5))+xlab("Log2(Mean+1) WT DMSO")+ylab("Log2(Mean+1) WT BAY-2416964")+xlim(0,15)+ylim(0,15)+geom_abline(slope=1, intercept=0,linetype="dotted")+annotate(geom="text", x=4, y=11, label="41 genes up after BAY-2416964 vs DMSO WT",color="red",size=6)+annotate(geom="text", x=4, y=12, label="15 genes down after BAY-2416964 vs DMSO WT",color="blue",size=6)+geom_text_repel(data=top_15_up_dn,aes(x=DMSO_log, y=bay_log,label=Symbol),color="black",arrow=arrow(ends="last",type="open"),fontface="bold")



