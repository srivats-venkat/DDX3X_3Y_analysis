library("ggplot2")
library("DESeq2")
library("SaTAnn")
library("DEXSeq")
library("riborex")
library("GenomicAlignments")
library("GenomicFiles")
library("GenomicFeatures")
library("rtracklayer")
library("gplots")
library("reshape2")
library("ggrepel")
library("topGO")
library("org.Hs.eg.db")
library("gridExtra")
library("cowplot")
library("ggridges")
library("mixtools")
library("corrgram")
library("Biostrings")
library("ggpubr")
library("patchwork")
library("RiboseQC")
library("coRdon")
library("msa")
library("seqinr")
library("ape")
library("biomaRt")
library("tidytree")
library("ggtree")
library("effsize")

#laods the annotation files
load_annotation("/Users/srivatsv/Box Sync/PostDoc Stephen Floor Lab/annotations/Homo.sapiens/gencode.v33.primary_assembly.annotation.gtf_Rannot")
load("/Users/srivatsv/Box Sync/PostDoc Stephen Floor Lab/annotations/Homo.sapiens/tx_regions_features.RData")

load("./counts/DMSO_HCT116_rep1_ribo_Aligned.sortedByCoord.out.bam_counts_summary_25min")
dmso_ribo_rep1_25min <- res_counts$counts_gene_ex_unq
load("./counts/DMSO_HCT116_rep2_ribo_Aligned.sortedByCoord.out.bam_counts_summary_25min")
dmso_ribo_rep2_25min <- res_counts$counts_gene_ex_unq
load("./counts/IAA_HCT116_rep1_ribo_Aligned.sortedByCoord.out.bam_counts_summary_25min")
iaa_ribo_rep1_25min <- res_counts$counts_gene_ex_unq
load("./counts/IAA_HCT116_rep2_ribo_Aligned.sortedByCoord.out.bam_counts_summary_25min")
iaa_ribo_rep2_25min <- res_counts$counts_gene_ex_unq
load("./counts/DMSO_HCT116_rep1_rna_Aligned.sortedByCoord.out.bam_counts_summary_25min")
dmso_rna_rep1_25min <- res_counts$counts_gene_ex_unq
load("./counts/DMSO_HCT116_rep2_rna_Aligned.sortedByCoord.out.bam_counts_summary_25min")
dmso_rna_rep2_25min <- res_counts$counts_gene_ex_unq
load("./counts/IAA_HCT116_rep1_rna_Aligned.sortedByCoord.out.bam_counts_summary_25min")
iaa_rna_rep1_25min <- res_counts$counts_gene_ex_unq
load("./counts/IAA_HCT116_rep2_rna_Aligned.sortedByCoord.out.bam_counts_summary_25min")
iaa_rna_rep2_25min <- res_counts$counts_gene_ex_unq
rm(res_counts)

gene_ex_unq_25min <- cbind(dmso_ribo_rep1_25min, 
                           dmso_ribo_rep2_25min, 
                           iaa_ribo_rep1_25min, 
                           iaa_ribo_rep2_25min, 
                           dmso_rna_rep1_25min, 
                           dmso_rna_rep2_25min, 
                           iaa_rna_rep1_25min, 
                           iaa_rna_rep2_25min)

colnames(gene_ex_unq_25min) <- c("dmso_ribo_rep1_25min", 
                                 "dmso_ribo_rep2_25min", 
                                 "iaa_ribo_rep1_25min", 
                                 "iaa_ribo_rep2_25min", 
                                 "dmso_rna_rep1_25min", 
                                 "dmso_rna_rep2_25min", 
                                 "iaa_rna_rep1_25min", 
                                 "iaa_rna_rep2_25min")

gene_ex_unq_25min <- as.data.frame(gene_ex_unq_25min)

rm(list = ls(pattern = "rep"))

#loads exp_design file and sets parameter values
exp_design <- read.delim("./sample_info.txt")
exp_design$Genotype <- factor(exp_design$Genotype, levels = c("dmso","iaa"))
exp_design$Condition <- factor(exp_design$Condition, levels = c( "rna","ribo"))

#DESeq for the riboRNA interaction
deseq_dataset_ex_25min <- DESeqDataSetFromMatrix(gene_ex_unq_25min, exp_design, design = ~ Condition + Genotype + Condition:Genotype)
deseq_dataset_ex_25min <- DESeq(deseq_dataset_ex_25min, test = "LRT", reduced = ~ Condition + Genotype)
deseq_results_ex_25min <- results(deseq_dataset_ex_25min, independentFiltering = FALSE)
deseq_results_ex_25min$pvalue[deseq_results_ex_25min$baseMean < 20] <- NA
deseq_results_ex_25min$padj <- p.adjust(deseq_results_ex_25min$pvalue, method = "BH")

#subsets the exp design and count files for rna only and ribo only DE
exp_design_rnaonly <- exp_design[c(5:8),]
exp_design_riboonly <- exp_design[c(1:4),]

gene_ex_unq_25min_rnaonly <- gene_ex_unq_25min[,c(5:8)]
gene_ex_unq_25min_riboonly <- gene_ex_unq_25min[,c(1:4)]

#DESeq for the individual rna and ribo terms
deseq_dataset_ex_25min_rnaonly <- DESeqDataSetFromMatrix(gene_ex_unq_25min_rnaonly, exp_design_rnaonly, design = ~ Genotype)
deseq_dataset_ex_25min_rnaonly <- DESeq(deseq_dataset_ex_25min_rnaonly)
deseq_results_ex_25min_rnaonly <- results(deseq_dataset_ex_25min_rnaonly)
deseq_results_ex_25min_rnaonly$pvalue[deseq_results_ex_25min_rnaonly$baseMean < 20] <- NA
deseq_results_ex_25min_rnaonly$padj <- p.adjust(deseq_results_ex_25min_rnaonly$pvalue, method = "BH")

deseq_dataset_ex_25min_riboonly <- DESeqDataSetFromMatrix(gene_ex_unq_25min_riboonly, exp_design_riboonly, design = ~ Genotype)
deseq_dataset_ex_25min_riboonly <- DESeq(deseq_dataset_ex_25min_riboonly)
deseq_results_ex_25min_riboonly <- results(deseq_dataset_ex_25min_riboonly)
deseq_results_ex_25min_riboonly$pvalue[deseq_results_ex_25min_riboonly$baseMean < 20] <- NA
deseq_results_ex_25min_riboonly$padj <- p.adjust(deseq_results_ex_25min_riboonly$pvalue, method = "BH")

#calls annotations from the GTF. Gene name, biotype 
deseq_results_ex_25min$symb<-GTF_annotation$trann$gene_name[match(rownames(deseq_results_ex_25min),GTF_annotation$trann$gene_id)]
deseq_results_ex_25min$biotype<-GTF_annotation$trann$gene_biotype[match(rownames(deseq_results_ex_25min),GTF_annotation$trann$gene_id)]
deseq_results_ex_25min_rnaonly$symb<-GTF_annotation$trann$gene_name[match(rownames(deseq_results_ex_25min_rnaonly),GTF_annotation$trann$gene_id)]
deseq_results_ex_25min_rnaonly$biotype<-GTF_annotation$trann$gene_biotype[match(rownames(deseq_results_ex_25min_rnaonly),GTF_annotation$trann$gene_id)]
deseq_results_ex_25min_riboonly$symb<-GTF_annotation$trann$gene_name[match(rownames(deseq_results_ex_25min_riboonly),GTF_annotation$trann$gene_id)]
deseq_results_ex_25min_riboonly$biotype<-GTF_annotation$trann$gene_biotype[match(rownames(deseq_results_ex_25min_riboonly),GTF_annotation$trann$gene_id)]

#builds a common results table with relevant parameters
results_ex_25min <- deseq_results_ex_25min[,c(1,2,6,7,8)]
colnames(results_ex_25min) <- c("baseMean_riborna", "log2FoldChange_riborna", "padj_riborna","symb", "biotype")

results_ex_25min_rnaonly <- deseq_results_ex_25min_rnaonly[,c(1,2,6)]
colnames(results_ex_25min_rnaonly) <- c("baseMean_rna", "log2FoldChange_rna", "padj_rna")

results_ex_25min_riboonly <- deseq_results_ex_25min_riboonly[,c(1,2,6)]
colnames(results_ex_25min_riboonly) <- c("baseMean_ribo", "log2FoldChange_ribo", "padj_ribo")

results_ex_25min <- cbind(results_ex_25min, results_ex_25min_rnaonly, results_ex_25min_riboonly)

results_ex_25min <- results_ex_25min[order(results_ex_25min$padj_riborna),]

#imposes read number cutoffs
cutoff_ex<-which(results_ex_25min$baseMean_rna>20 & results_ex_25min$baseMean_ribo>20)
results_ex_25min <- results_ex_25min[cutoff_ex,]

# classifies genes based on the effect of perturbation
results_ex_25min$tx_type<-"Not_significant"
results_ex_25min$tx_type[results_ex_25min$padj_rna<.01 & results_ex_25min$log2FoldChange_rna<0]<-"RNA_down"
results_ex_25min$tx_type[results_ex_25min$padj_rna<.01 & results_ex_25min$log2FoldChange_rna>0]<-"RNA_up"
results_ex_25min$tx_type[results_ex_25min$padj_riborna<.05 & results_ex_25min$log2FoldChange_riborna<0]<-"TE_down"
results_ex_25min$tx_type[results_ex_25min$padj_riborna<.05 & results_ex_25min$log2FoldChange_riborna>0]<-"TE_up"

#log of p-value
results_ex_25min$log10_padj_riborna<-log10(results_ex_25min$padj_riborna)
results_ex_25min$log10_padj_riborna[is.na(results_ex_25min$log10_padj_riborna)]<-0

results_ex_25min$tx_type<-factor(results_ex_25min$tx_type,levels=c("TE_up","TE_down","RNA_up","RNA_down","Not_significant"))
results_ex_25min$biotype[grep(results_ex_25min$symb,pattern = "HIST|H2B")]<-"Histone genes"
results_ex_25min$biotype[grep(results_ex_25min$symb,pattern = "HIST|H2B",invert =  T)]<-"all genes"
results_ex_25min$biotype<-factor(results_ex_25min$biotype,levels=c("all genes","Histone genes"))

results_ex_25min$symb_deplot<-as.character(results_ex_25min$symb)
results_ex_25min$symb_deplot<-NA
toprna_ex<-results_ex_25min[results_ex_25min$tx_type=="RNA_up",]
toprna_ex<-toprna_ex[order(toprna_ex$log2FoldChange_rna,toprna_ex$log2FoldChange_ribo,decreasing = T),"symb"]
if(length(toprna_ex)>4){toprna_ex=toprna_ex[1:3]}
topribo_ex<-results_ex_25min[results_ex_25min$tx_type=="TE_up",]
topribo_ex<-topribo_ex[order(topribo_ex$log2FoldChange_riborna,decreasing = T),"symb"]
if(length(topribo_ex)>5){topribo_ex=topribo_ex[1:4]}
bottomrna_ex<-results_ex_25min[results_ex_25min$tx_type=="RNA_down",]
bottomrna_ex<-bottomrna_ex[order(bottomrna_ex$log2FoldChange_rna,bottomrna_ex$log2FoldChange_ribo,decreasing = F),"symb"]
if(length(bottomrna_ex)>4){bottomrna_ex=bottomrna_ex[1:3]}
bottomribo_ex<-results_ex_25min[results_ex_25min$tx_type=="TE_down",]
bottomribo_ex<-bottomribo_ex[order(bottomribo_ex$log2FoldChange_riborna,decreasing = F),"symb"]
if(length(bottomribo_ex)>5){bottomribo_ex=c(bottomribo_ex[1:4],bottomribo_ex[bottomribo_ex%in%c("ODC1","DVL1","DVL2","HMBS","PRKRA")])}
results_ex_25min$symb_deplot[results_ex_25min$symb%in%c(toprna_ex,topribo_ex,bottomrna_ex,bottomribo_ex)]<-results_ex_25min$symb[results_ex_25min$symb%in%c(toprna_ex,topribo_ex,bottomrna_ex,bottomribo_ex)]

results_ex_25min$nudgey<-NA
results_ex_25min$nudgey[which(nchar(results_ex_25min$symb_deplot)>2)]<-ifelse(results_ex_25min$log2FoldChange_ribo[which(nchar(results_ex_25min$symb_deplot)>2)]>0, .4, -.5)
results_ex_25min$nudgex<-NA
results_ex_25min$nudgex[which(nchar(results_ex_25min$symb_deplot)>2)]<-ifelse(results_ex_25min$log2FoldChange_rna[which(nchar(results_ex_25min$symb_deplot)>2)]>0, .3, -.4)
results_ex_25min$nudgex[results_ex_25min$tx_type=="TE_up"]<-0
results_ex_25min$nudgey[results_ex_25min$tx_type=="RNA_up"]<-0
results_ex_25min$nudgey[results_ex_25min$tx_type=="RNA_down"]<-0
results_ex_25min$nudgex[results_ex_25min$tx_type=="RNA_up"]<-0
results_ex_25min$nudgex[results_ex_25min$tx_type=="RNA_down"]<-0

results_ex_25min_df <- as.data.frame(results_ex_25min)

write.csv(results_ex_25min_df, "data/results_ex_25min_df.csv")

DE_plot_ex_25min <- ggplot(results_ex_25min_df,aes(x=log2FoldChange_rna,y=log2FoldChange_ribo,color=tx_type,size=-log10_padj_riborna,label=symb_deplot)) + 
  geom_point() + 
  theme_bw() +
  ylab("Ribo fold change (IAA/DMSO)") +
  xlab("RNA fold change (IAA/DMSO)") +
  geom_vline(xintercept = 0, alpha = 0.5) +
  geom_hline(yintercept = 0, alpha = 0.5) +
  theme(axis.title.x = element_text(size=22),axis.text.x  = element_text(angle=45, vjust=0.5, size=18)) +
  theme(axis.title.y = element_text(size=24),axis.text.y  = element_text(angle=45, vjust=0.5, size=18))  +
  scale_color_manual(values = alpha(c("dark red","blue","firebrick1","cornflowerblue","gray24"),c(.8,.8,.8,.8,.3)),"Tx_class") +
  scale_size_continuous(name = "-log10 adj.pval\nRiboRNA") +
  scale_x_continuous(limits = c(-2,2)) + 
  scale_y_continuous(limits = c(-2,2)) +
  geom_text_repel(size=5,nudge_x = results_ex_25min_df$nudgex[which(nchar(results_ex_25min_df$symb_deplot)>2)],
                  nudge_y = results_ex_25min_df$nudgey[which(nchar(results_ex_25min_df$symb_deplot)>2)],force=20) +
  theme(legend.text = element_text(size = rel(1.5)), legend.title = element_text(size = rel(1.5)))

pdf(file = "figures/DE_plot_ex_25min_fewer.pdf", width = 10, height = 7)
print(DE_plot_ex_25min)
dev.off()

png(file = "figures/DE_plot_ex_25min_fewer.png", width = 10, height = 7, units = "in", res = 600)
print(DE_plot_ex_25min)
dev.off()

split_results_ex_25min <- results_ex_25min_df
split_results_ex_25min$gene_id <- rownames(split_results_ex_25min)
split_results_ex_25min <- split(split_results_ex_25min$gene_id, f = split_results_ex_25min$tx_type)
for(i in names(split_results_ex_25min)[1:5]){
  tx_regions_feats$tx_type[tx_regions_feats$gene_id%in%split_results_ex_25min[[i]]]<-i
}

features_df <- as.data.frame(mcols(tx_regions_feats))
features_df$tx_type<-factor(features_df$tx_type,levels=rev(c("Not_significant",
                                                             "RNA_up",
                                                             "RNA_down",
                                                             "TE_up",
                                                             "TE_down")))
features_df_completecases <- features_df[complete.cases(features_df),]

five_utr_length <- ggplot(features_df_completecases,aes(len_5ut+1,group=tx_type,y=tx_type,fill=tx_type))  +
  theme_ridges() +
  scale_fill_manual(values=c("blue","dark red","cornflowerblue","firebrick1","dark grey")) +
  xlab("5'UTR length") +
  stat_density_ridges(quantile_lines = TRUE, quantiles = 2,from = 0) +
  theme(legend.position="none") +
  theme(axis.title.x = element_text(size=22),axis.text.x  = element_text( vjust=0.5, size=18)) +
  theme(axis.title.y = element_text(size=24),axis.text.y  = element_text( vjust=0.5, size=18))  +
  scale_x_log10(breaks=c(1,11,101,1001,10001),labels=c(0,10,100,1000,10000))

five_utr_gc <- ggplot(features_df_completecases,aes(G.C,group=tx_type,y=tx_type,fill=tx_type))  +
  theme_ridges() +
  scale_fill_manual(values=c("blue","dark red","cornflowerblue","firebrick1","dark grey")) +
  xlab("%GC content in 5'UTR")  +
  stat_density_ridges(quantile_lines = TRUE, quantiles = 2) +
  theme(legend.position="none") +
  theme(axis.title.x = element_text(size=22),axis.text.x  = element_text( vjust=0.5, size=18)) +
  theme(axis.title.y = element_text(size=24),axis.text.y  = element_text( vjust=0.5, size=18))

three_utr_length <- ggplot(features_df_completecases,aes(len_3ut+1,group=tx_type,y=tx_type,fill=tx_type))  +
  theme_ridges() +
  scale_fill_manual(values=c("blue","dark red","cornflowerblue","firebrick1","dark grey")) +
  xlab("3'UTR length") +
  stat_density_ridges(quantile_lines = TRUE, quantiles = 2,from = 0) +
  theme(legend.position="none") +
  theme(axis.title.x = element_text(size=22),axis.text.x  = element_text( vjust=0.5, size=18)) +
  theme(axis.title.y = element_text(size=24),axis.text.y  = element_text( vjust=0.5, size=18))  +
  scale_x_log10(breaks=c(1,11,101,1001,10001),labels=c(0,10,100,1000,10000))

three_utr_gc <- ggplot(features_df_completecases,aes(G.C.2,group=tx_type,y=tx_type,fill=tx_type))  +
  theme_ridges() +
  scale_fill_manual(values=c("blue","dark red","cornflowerblue","firebrick1","dark grey")) +
  xlab("%GC content in 3'UTR")  +
  stat_density_ridges(quantile_lines = TRUE, quantiles = 2) +
  theme(legend.position="none") +
  theme(axis.title.x = element_text(size=22),axis.text.x  = element_text( vjust=0.5, size=18)) +
  theme(axis.title.y = element_text(size=24),axis.text.y  = element_text( vjust=0.5, size=18))

cds_length <- ggplot(features_df_completecases,aes(len_cds+1,group=tx_type,y=tx_type,fill=tx_type))  +
  theme_ridges() +
  scale_fill_manual(values=c("blue","dark red","cornflowerblue","firebrick1","dark grey")) +
  xlab("CDS length") +
  stat_density_ridges(quantile_lines = TRUE, quantiles = 2,from = 0) +
  theme(legend.position="none") +
  theme(axis.title.x = element_text(size=22),axis.text.x  = element_text( vjust=0.5, size=18)) +
  theme(axis.title.y = element_text(size=24),axis.text.y  = element_text( vjust=0.5, size=18))  +
  scale_x_log10(breaks=c(1,11,101,1001,10001),labels=c(0,10,100,1000,10000))

cds_gc <- ggplot(features_df_completecases,aes(G.C.1,group=tx_type,y=tx_type,fill=tx_type))  +
  theme_ridges() +
  scale_fill_manual(values=c("blue","dark red","cornflowerblue","firebrick1","dark grey")) +
  xlab("%GC content in CDS")  +
  stat_density_ridges(quantile_lines = TRUE, quantiles = 2) +
  theme(axis.title.x = element_text(size=22),axis.text.x  = element_text( vjust=0.5, size=18)) +
  theme(axis.title.y = element_text(size=24),axis.text.y  = element_text( vjust=0.5, size=18))  +
  theme(legend.position="none")

tx_length <- ggplot(features_df_completecases,aes(lentx+1,group=tx_type,y=tx_type,fill=tx_type))  +
  theme_ridges() +
  scale_fill_manual(values=c("blue","dark red","cornflowerblue","firebrick1","dark grey")) +
  xlab("Transcript length") +
  stat_density_ridges(quantile_lines = TRUE, quantiles = 2,from = 0) +
  theme(legend.position="none") +
  theme(axis.title.x = element_text(size=22),axis.text.x  = element_text( vjust=0.5, size=18)) +
  theme(axis.title.y = element_text(size=24),axis.text.y  = element_text( vjust=0.5, size=18))  +
  scale_x_log10(breaks=c(1,11,101,1001,10001),labels=c(0,10,100,1000,10000))

tx_gc <- ggplot(features_df_completecases,aes(G.C.3,group=tx_type,y=tx_type,fill=tx_type))  +
  theme_ridges() +
  scale_fill_manual(values=c("blue","dark red","cornflowerblue","firebrick1","dark grey")) +
  xlab("%GC content in transcript")  +
  stat_density_ridges(quantile_lines = TRUE, quantiles = 2) +
  theme(axis.title.x = element_text(size=22),axis.text.x  = element_text( vjust=0.5, size=18)) +
  theme(axis.title.y = element_text(size=24),axis.text.y  = element_text( vjust=0.5, size=18))  +
  theme(legend.position="none")

DE_5utr_feats <- five_utr_length | five_utr_gc
DE_3utr_feats <- three_utr_length | three_utr_gc
DE_cds_feats <- cds_length | cds_gc
DE_tx_feats <- tx_length | tx_gc

pdf(file = "figures/DE_5utr_feats.pdf",width = 12.5,height = 5)
print(DE_5utr_feats)
dev.off()

pdf(file = "figures/DE_3utr_feats.pdf",width = 12.5,height = 5)
print(DE_3utr_feats)
dev.off()

pdf(file = "figures/DE_cds_feats.pdf",width = 12.5,height = 5)
print(DE_cds_feats)
dev.off()

pdf(file = "figures/DE_tx_feats.pdf",width = 12.5,height = 5)
print(DE_tx_feats)
dev.off()

#loads and prepares a table for 5'UTR/CDS read ratios for DMSO and IAA
load("./counts/DMSO_HCT116_rep1_ribo_Aligned.sortedByCoord.out.bam_counts_summary_25min")
utr_counts <- res_counts$counts_5cds3_unq
utr_counts <- as.data.frame(utr_counts)
colnames(utr_counts) <- c("fiveutr", "cds" , "threeutr")
utr_counts$fiveutrbycds <- utr_counts$fiveutr/utr_counts$cds
utr_counts_ratio <- as.data.frame(utr_counts[,4])
rownames(utr_counts_ratio) <- rownames(utr_counts)
load("./counts/DMSO_HCT116_rep2_ribo_Aligned.sortedByCoord.out.bam_counts_summary_25min")
utr_counts <- res_counts$counts_5cds3_unq
utr_counts <- as.data.frame(utr_counts)
colnames(utr_counts) <- c("fiveutr", "cds" , "threeutr")
utr_counts$fiveutrbycds <- utr_counts$fiveutr/utr_counts$cds
utr_counts_ratio <- cbind(utr_counts_ratio, utr_counts$fiveutrbycds)
load("./counts/IAA_HCT116_rep1_ribo_Aligned.sortedByCoord.out.bam_counts_summary_25min")
utr_counts <- res_counts$counts_5cds3_unq
utr_counts <- as.data.frame(utr_counts)
colnames(utr_counts) <- c("fiveutr", "cds" , "threeutr")
utr_counts$fiveutrbycds <- utr_counts$fiveutr/utr_counts$cds
utr_counts_ratio <- cbind(utr_counts_ratio, utr_counts$fiveutrbycds)
load("./counts/IAA_HCT116_rep2_ribo_Aligned.sortedByCoord.out.bam_counts_summary_25min")
utr_counts <- res_counts$counts_5cds3_unq
utr_counts <- as.data.frame(utr_counts)
colnames(utr_counts) <- c("fiveutr", "cds" , "threeutr")
utr_counts$fiveutrbycds <- utr_counts$fiveutr/utr_counts$cds
utr_counts_ratio <- cbind(utr_counts_ratio, utr_counts$fiveutrbycds)

rm(res_counts)

colnames(utr_counts_ratio) <- c("dmso_ribo_rep1", "dmso_ribo_rep2", "iaa_ribo_rep1", "iaa_ribo_rep2")

utr_counts_ratio <- utr_counts_ratio[complete.cases(utr_counts_ratio), ]
utr_counts_ratio <- utr_counts_ratio[!is.infinite(rowSums(utr_counts_ratio)),]

utr_counts_ratio$dmso <- (utr_counts_ratio$dmso_ribo_rep1 + utr_counts_ratio$dmso_ribo_rep2)/2
utr_counts_ratio$iaa <- (utr_counts_ratio$iaa_ribo_rep1 + utr_counts_ratio$iaa_ribo_rep2)/2
utr_counts_ratio <- utr_counts_ratio[,c(5:6)]

utr_counts_ratio$iaa_dmso <- log2(utr_counts_ratio$iaa/utr_counts_ratio$dmso)
utr_counts_ratio <- utr_counts_ratio[complete.cases(utr_counts_ratio), ]
utr_counts_ratio <- utr_counts_ratio[!is.infinite(rowSums(utr_counts_ratio)),]

utr_counts_ratio$symb<-GTF_annotation$trann$gene_name[match(rownames(utr_counts_ratio),GTF_annotation$trann$gene_id)]
utr_counts_ratio$biotype<-GTF_annotation$trann$gene_biotype[match(rownames(utr_counts_ratio),GTF_annotation$trann$gene_id)]

TE_down_utrratio <- results_ex_25min_df[which(results_ex_25min_df$tx_type == "TE_down"),]
TE_up_utrratio <- results_ex_25min_df[which(results_ex_25min_df$tx_type == "TE_up"),]
RNA_down_utrratio <- results_ex_25min_df[which(results_ex_25min_df$tx_type == "RNA_down"),]
RNA_up_utrratio <- results_ex_25min_df[which(results_ex_25min_df$tx_type == "RNA_up"),]
not_significant_utrratio <- results_ex_25min_df[which(results_ex_25min_df$tx_type == "Not_significant"),]

TE_down_utrratio$five_cds_iaa_dmso<-utr_counts_ratio$iaa_dmso[match(rownames(TE_down_utrratio),rownames(utr_counts_ratio))]
TE_up_utrratio$five_cds_iaa_dmso<-utr_counts_ratio$iaa_dmso[match(rownames(TE_up_utrratio),rownames(utr_counts_ratio))]
RNA_down_utrratio$five_cds_iaa_dmso<-utr_counts_ratio$iaa_dmso[match(rownames(RNA_down_utrratio),rownames(utr_counts_ratio))]
RNA_up_utrratio$five_cds_iaa_dmso<-utr_counts_ratio$iaa_dmso[match(rownames(RNA_up_utrratio),rownames(utr_counts_ratio))]
not_significant_utrratio$five_cds_iaa_dmso<-utr_counts_ratio$iaa_dmso[match(rownames(not_significant_utrratio),rownames(utr_counts_ratio))]

all_utrratio <- results_ex_25min_df 
all_utrratio$five_cds_iaa_dmso<-utr_counts_ratio$iaa_dmso[match(rownames(all_utrratio),rownames(utr_counts_ratio))]

load("/Users/srivatsv/Box Sync/PostDoc Stephen Floor Lab/Data/DDX3X_DDX3Y_riboprofile_Nov2019/counts/IAA3XWT_HCT116_rep1_ribo_Aligned.sortedByCoord.out.bam_counts_summary_25min")
utr_counts <- res_counts$counts_5cds3_unq
utr_counts <- as.data.frame(utr_counts)
colnames(utr_counts) <- c("fiveutr", "cds" , "threeutr")
utr_counts$fiveutrbycds <- utr_counts$fiveutr/utr_counts$cds
utr_counts_ratio <- as.data.frame(utr_counts[,4])
rownames(utr_counts_ratio) <- rownames(utr_counts)
load("/Users/srivatsv/Box Sync/PostDoc Stephen Floor Lab/Data/DDX3X_DDX3Y_riboprofile_Nov2019/counts/IAA3XWT_HCT116_rep2_ribo_Aligned.sortedByCoord.out.bam_counts_summary_25min")
utr_counts <- res_counts$counts_5cds3_unq
utr_counts <- as.data.frame(utr_counts)
colnames(utr_counts) <- c("fiveutr", "cds" , "threeutr")
utr_counts$fiveutrbycds <- utr_counts$fiveutr/utr_counts$cds
utr_counts_ratio <- cbind(utr_counts_ratio, utr_counts$fiveutrbycds)
load("/Users/srivatsv/Box Sync/PostDoc Stephen Floor Lab/Data/DDX3X_DDX3Y_riboprofile_Nov2019/counts/IAA3Y_HCT116_rep1_ribo_Aligned.sortedByCoord.out.bam_counts_summary_25min")
utr_counts <- res_counts$counts_5cds3_unq
utr_counts <- as.data.frame(utr_counts)
colnames(utr_counts) <- c("fiveutr", "cds" , "threeutr")
utr_counts$fiveutrbycds <- utr_counts$fiveutr/utr_counts$cds
utr_counts_ratio <- cbind(utr_counts_ratio, utr_counts$fiveutrbycds)
load("/Users/srivatsv/Box Sync/PostDoc Stephen Floor Lab/Data/DDX3X_DDX3Y_riboprofile_Nov2019/counts/IAA3Y_HCT116_rep2_ribo_Aligned.sortedByCoord.out.bam_counts_summary_25min")
utr_counts <- res_counts$counts_5cds3_unq
utr_counts <- as.data.frame(utr_counts)
colnames(utr_counts) <- c("fiveutr", "cds" , "threeutr")
utr_counts$fiveutrbycds <- utr_counts$fiveutr/utr_counts$cds
utr_counts_ratio <- cbind(utr_counts_ratio, utr_counts$fiveutrbycds)

rm(res_counts)

colnames(utr_counts_ratio) <- c("IAA3XWT_ribo_rep1", "IAA3XWT_ribo_rep2", "IAA3Y_ribo_rep1", "IAA3Y_ribo_rep2")

utr_counts_ratio <- utr_counts_ratio[complete.cases(utr_counts_ratio), ]
utr_counts_ratio <- utr_counts_ratio[!is.infinite(rowSums(utr_counts_ratio)),]

utr_counts_ratio$IAA3XWT <- (utr_counts_ratio$IAA3XWT_ribo_rep1 + utr_counts_ratio$IAA3XWT_ribo_rep2)/2
utr_counts_ratio$IAA3Y <- (utr_counts_ratio$IAA3Y_ribo_rep1 + utr_counts_ratio$IAA3Y_ribo_rep2)/2
utr_counts_ratio <- utr_counts_ratio[,c(5:6)]

utr_counts_ratio$IAA3Y_IAA3XWT <- log2(utr_counts_ratio$IAA3Y/utr_counts_ratio$IAA3XWT)
utr_counts_ratio <- utr_counts_ratio[complete.cases(utr_counts_ratio), ]
utr_counts_ratio <- utr_counts_ratio[!is.infinite(rowSums(utr_counts_ratio)),]

utr_counts_ratio$symb<-GTF_annotation$trann$gene_name[match(rownames(utr_counts_ratio),GTF_annotation$trann$gene_id)]
utr_counts_ratio$biotype<-GTF_annotation$trann$gene_biotype[match(rownames(utr_counts_ratio),GTF_annotation$trann$gene_id)]

TE_down_utrratio$five_cds_IAA3Y_IAA3XWT<-utr_counts_ratio$IAA3Y_IAA3XWT[match(rownames(TE_down_utrratio),rownames(utr_counts_ratio))]
TE_up_utrratio$five_cds_IAA3Y_IAA3XWT<-utr_counts_ratio$IAA3Y_IAA3XWT[match(rownames(TE_up_utrratio),rownames(utr_counts_ratio))]
RNA_down_utrratio$five_cds_IAA3Y_IAA3XWT<-utr_counts_ratio$IAA3Y_IAA3XWT[match(rownames(RNA_down_utrratio),rownames(utr_counts_ratio))]
RNA_up_utrratio$five_cds_IAA3Y_IAA3XWT<-utr_counts_ratio$IAA3Y_IAA3XWT[match(rownames(RNA_up_utrratio),rownames(utr_counts_ratio))]
not_significant_utrratio$five_cds_IAA3Y_IAA3XWT<-utr_counts_ratio$IAA3Y_IAA3XWT[match(rownames(not_significant_utrratio),rownames(utr_counts_ratio))]

all_utrratio$five_cds_IAA3Y_IAA3XWT<-utr_counts_ratio$IAA3Y_IAA3XWT[match(rownames(all_utrratio),rownames(utr_counts_ratio))]
all_utrratio$Tx_class<-factor(all_utrratio$tx_type,levels=rev(c("Not_significant","RNA_up","RNA_down","TE_up","TE_down")))

medianWithoutNA<-function(x) {
  median(x[which(!is.na(x))])
}

write.csv(all_utrratio, "data/all_utrratio.csv")

iaavsdmso_5riboskew <- ggplot(data = all_utrratio, aes(y = five_cds_iaa_dmso, x = Tx_class)) +
  theme_minimal() +
  ylim(-2,2) +
  geom_violin(aes(fill = Tx_class)) +
  geom_boxplot(width = 0.1, alpha = 0) +
  scale_fill_manual(values=c("blue","dark red","cornflowerblue","firebrick1","dark grey"), name = "Transcript Class") +
  xlab("Transcript Class") +
  ylab("Log2FC 5' Riboskew (IAA/DMSO)") +
  theme(plot.title = element_text(size = 24)) +
  theme(axis.title.x = element_text(size = 16)) +
  theme(axis.title.y = element_text(size = 16)) +
  theme(axis.text.x = element_text(size = 12)) +
  theme(axis.text.y = element_text(size = 12)) +
  theme(legend.title = element_text(size = 16)) +
  theme(legend.text = element_text(size = 12)) +
  theme(legend.spacing.y = unit(20, 'points')) +
  theme(legend.key.size = unit(20, 'points')) +
  geom_hline(yintercept= medianWithoutNA(TE_down_utrratio$five_cds_iaa_dmso), linetype="dashed", color = "darkgreen", size = 1) +
  ggtitle("DDX3X degron vs. control")

YvsX_5riboskew <- ggplot(data = all_utrratio, aes(y = five_cds_IAA3Y_IAA3XWT, x = Tx_class)) +
  theme_minimal() +
  ylim(-2,2) +
  geom_violin(aes(fill = Tx_class)) +
  geom_boxplot(width = 0.1, alpha = 0) +
  scale_fill_manual(values=c("blue","dark red","cornflowerblue","firebrick1","dark grey"), name = "Transcript Class") +
  xlab("Transcript Class (as in IAA vs. DMSO)") +
  ylab("Log2FC 5' Riboskew (DDX3Y/DDX3X)") +
  theme(plot.title = element_text(size = 24)) +
  theme(axis.title.x = element_text(size = 16)) +
  theme(axis.title.y = element_text(size = 16)) +
  theme(axis.text.x = element_text(size = 12)) +
  theme(axis.text.y = element_text(size = 12)) +
  theme(legend.title = element_text(size = 16)) +
  theme(legend.text = element_text(size = 12)) +
  theme(legend.spacing.y = unit(20, 'points')) +
  theme(legend.key.size = unit(20, 'points')) +
  geom_hline(yintercept= medianWithoutNA(TE_down_utrratio$five_cds_IAA3Y_IAA3XWT), linetype="dashed", color = "darkgreen", size = 1) +
  ggtitle("DDX3Y vs. DDX3X")

riboskew_compare <- ggarrange(iaavsdmso_5riboskew,YvsX_5riboskew,nrow = 1)
ggsave(filename = "5riboskew_IAAvsDMSO_3Yvs3X.pdf", plot = riboskew_compare, device = "pdf", units = "in", width = 14, height = 10.8, path = "/Users/srivatsv/Box Sync/PostDoc Stephen Floor Lab/Data/DDX3X_degron_riboprofile_Nov2019/figures")

iaavsdmso_5riboskew_ridge <- ggplot(all_utrratio,aes(five_cds_iaa_dmso,group=Tx_class,y=Tx_class,fill=Tx_class))  +
  theme_ridges() +
  scale_fill_manual(values=c("blue","dark red","cornflowerblue","firebrick1","dark grey")) +
  xlab("5' Riboskew IAA/DMSO")  +
  stat_density_ridges(quantile_lines = TRUE, quantiles = 2) +
  theme(axis.title.x = element_text(size=22),axis.text.x  = element_text( vjust=0.5, size=18)) +
  theme(axis.title.y = element_text(size=24),axis.text.y  = element_text( vjust=0.5, size=18))  +
  scale_x_continuous(limits = c(-2,2)) +
  theme(legend.position="none")

YvsX_5riboskew_ridge <- ggplot(all_utrratio,aes(five_cds_IAA3Y_IAA3XWT,group=Tx_class,y=Tx_class,fill=Tx_class))  +
  theme_ridges() +
  scale_fill_manual(values=c("blue","dark red","cornflowerblue","firebrick1","dark grey")) +
  xlab("5' Riboskew DDX3Y/DDX3X")  +
  ylab("Tx_class \n (as in IAA/DMSO)") +
  stat_density_ridges(quantile_lines = TRUE, quantiles = 2) +
  theme(axis.title.x = element_text(size=22),axis.text.x  = element_text( vjust=0.5, size=18)) +
  theme(axis.title.y = element_text(size=24),axis.text.y  = element_text( vjust=0.5, size=18))  +
  scale_x_continuous(limits = c(-2,2)) +
  theme(legend.position="none")

riboskew_compare_ridges <- iaavsdmso_5riboskew_ridge | YvsX_5riboskew_ridge

temp2 <- all_utrratio[which(all_utrratio$tx_type == "Not_significant"),]
temp2 <- temp2$five_cds_iaa_dmso
temp <- all_utrratio[which(all_utrratio$tx_type == "TE_down"),]
temp <- temp$five_cds_iaa_dmso
riboskew_iaadmso_clifd <- cliff.delta(temp, temp2, return.dm = TRUE)

temp2 <- all_utrratio[which(all_utrratio$tx_type == "Not_significant"),]
temp2 <- temp2$five_cds_IAA3Y_IAA3XWT
temp <- all_utrratio[which(all_utrratio$tx_type == "TE_down"),]
temp <- temp$five_cds_IAA3Y_IAA3XWT
riboskew_3Yvs3X_clifd <- cliff.delta(temp, temp2, return.dm = TRUE)

pdf(file = "figures/5'riboskew_compare_ridges.pdf",width = 10,height = 5)
riboskew_compare_ridges
dev.off()

png(file = "figures/5'riboskew_compare_ridges.png", width = 10, height = 5, units = "in", res = 600)
print(riboskew_compare_ridges)
dev.off()

deltaTE_vs_baseTE_degron_df <- results_ex_25min_df
gene_ex_unq_25min$TE_base <- log2((gene_ex_unq_25min$dmso_ribo_rep1_25min + gene_ex_unq_25min$dmso_ribo_rep2_25min) / (gene_ex_unq_25min$dmso_rna_rep1_25min + gene_ex_unq_25min$dmso_rna_rep2_25min))
deltaTE_vs_baseTE_degron_df$TE_base <- gene_ex_unq_25min$TE_base[match(rownames(deltaTE_vs_baseTE_degron_df), rownames(gene_ex_unq_25min))]

deltaTE_vs_baseTE_degron <- ggplot(deltaTE_vs_baseTE_degron_df,aes(x=TE_base,y=log2FoldChange_riborna,color=tx_type,size=-log10_padj_riborna)) + 
  geom_point() + 
  theme_bw() +
  ylab("Log2 FoldChange TE \n(IAA/DMSO) ") +
  xlab("Baseline TE in DMSO") +
  geom_vline(xintercept = 0, alpha = 0.5) +
  geom_hline(yintercept = 0, alpha = 0.5) +
  theme(axis.title.x = element_text(size=22),axis.text.x  = element_text(angle=45, vjust=0.5, size=18)) +
  theme(axis.title.y = element_text(size=24),axis.text.y  = element_text(angle=45, vjust=0.5, size=18))  +
  scale_color_manual(values = alpha(c("dark red","blue","firebrick1","cornflowerblue","gray24"),c(.8,.8,.8,.8,.3)),"Tx_class") +
  scale_size_continuous(name = "-log10 adj.pval\nRiboRNA") +
  scale_x_continuous(limits = c(-5,2.5)) +
  scale_y_continuous(limits = c(-2,2)) +
  theme(legend.text = element_text(size = rel(1.5)), legend.title = element_text(size = rel(1.5)))

pdf(file = "figures/deltaTE_vs_baseTE_degron.pdf", height = 7, width = 10)
print(deltaTE_vs_baseTE_degron)
dev.off()

png(file = "figures/deltaTE_vs_baseTE_degron.png", width = 10, height = 7, units = "in", res = 600)
print(deltaTE_vs_baseTE_degron)
dev.off()


deltaTE_vs_baseTE_degron_df$Tx_class<-factor(deltaTE_vs_baseTE_degron_df$tx_type,
                                             levels=rev(c("Not_significant","RNA_up","RNA_down","TE_up","TE_down")))

txclass_vs_baseTE_degron_violin <- ggplot(data = deltaTE_vs_baseTE_degron_df, aes(y = TE_base, x = Tx_class)) +
  theme_minimal() +
  ylim(-5,2.5) +
  geom_violin(aes(fill = Tx_class)) +
  geom_boxplot(width = 0.1, alpha = 0) +
  scale_fill_manual(values=c("blue","dark red","cornflowerblue","firebrick1","dark grey"), name = "Transcript Class") +
  xlab("Transcript Class") +
  ylab("Baseline TE in DMSO") + 
  theme(plot.title = element_text(size = 24)) +
  theme(axis.title.x = element_text(size = 16)) +
  theme(axis.title.y = element_text(size = 16)) +
  theme(axis.text.x = element_text(size = 12)) +
  theme(axis.text.y = element_text(size = 12)) +
  theme(legend.title = element_text(size = 16)) +
  theme(legend.text = element_text(size = 12)) +
  theme(legend.spacing.y = unit(20, 'points')) +
  theme(legend.key.size = unit(20, 'points'))

pdf(file = "figures/txclass_vs_baseTE_degron_violin.pdf", height = 8, width = 8)
print(txclass_vs_baseTE_degron_violin)
dev.off()

deltaTE_vs_baseTE_degron_df <- deltaTE_vs_baseTE_degron_df[,-c(14:16)]
deltaTE_vs_baseTE_degron_df <- deltaTE_vs_baseTE_degron_df[complete.cases(deltaTE_vs_baseTE_degron_df),]
deltaTE_vs_baseTE_degron_df <- deltaTE_vs_baseTE_degron_df[is.finite(deltaTE_vs_baseTE_degron_df$TE_base),]

txclass_vs_baseTE_degron_ridges <- ggplot(deltaTE_vs_baseTE_degron_df,aes(TE_base,group=Tx_class,y=Tx_class,fill=Tx_class))  +
  theme_ridges() +
  scale_fill_manual(values=c("blue","dark red","cornflowerblue","firebrick1","dark grey")) +
  xlab("Base translation efficiency (DMSO)") +
  stat_density_ridges(quantile_lines = TRUE, quantiles = 2) +
  theme(axis.title.x = element_text(size=22),axis.text.x  = element_text( vjust=0.5, size=18)) +
  theme(axis.title.y = element_text(size=24),axis.text.y  = element_text( vjust=0.5, size=18))  +
  theme(legend.position="none")

pdf(file = "figures/txclass_vs_baseTE_degron_ridges.pdf", height = 5, width = 5)
print(txclass_vs_baseTE_degron_ridges)
dev.off()

png(file = "figures/txclass_vs_baseTE_degron_ridges.png", height = 5, width = 5, units = "in", res = 600)
print(txclass_vs_baseTE_degron_ridges)
dev.off()

txclass_vs_tAI_degron_ridges <- ggplot(features_df_completecases,aes(tai_fiveutr,group=tx_type,y=tx_type,fill=tx_type))  +
  theme_ridges() +
  scale_fill_manual(values=c("blue","dark red","cornflowerblue","firebrick1","dark grey")) +
  xlab("tRNA adaptation index") +
  stat_density_ridges(quantile_lines = TRUE, quantiles = 2,from = 0) +
  theme(legend.position="none") + 
  theme(axis.title.x = element_text(size=22),axis.text.x  = element_text( vjust=0.5, size=18)) +
  theme(axis.title.y = element_text(size=24),axis.text.y  = element_text( vjust=0.5, size=18))  +
  scale_x_continuous(breaks = c(0.2,0.3,0.4,0.5,0.6), limits = c(0.29,0.51))

pdf(file = "figures/txclass_vs_tAI_degron_ridges.pdf", width = 5, height = 5)
print(txclass_vs_tAI_degron_ridges)
dev.off()

png(file = "figures/txclass_vs_tAI_degron_ridges.png", width = 5, height = 5, units = "in", res = 600)
print(txclass_vs_tAI_degron_ridges)
dev.off()

txclass_vs_MFE_ridges <- ggplot(features_df_completecases,aes(MFE_fiveutr,group=tx_type,y=tx_type,fill=tx_type))  +
  theme_ridges() +
  scale_fill_manual(values=c("blue","dark red","cornflowerblue","firebrick1","dark grey")) +
  xlab("5'UTR MFE") +
  stat_density_ridges(quantile_lines = TRUE, quantiles = 2) +
  theme(axis.title.x = element_text(size=22),axis.text.x  = element_text( vjust=0.5, size=18)) +
  theme(axis.title.y = element_text(size=24),axis.text.y  = element_text( vjust=0.5, size=18))  +
  theme(legend.position="none")

txclass_vs_normMFE_ridges <- ggplot(features_df_completecases,aes(normMFE,group=tx_type,y=tx_type,fill=tx_type))  +
  theme_ridges() +
  scale_fill_manual(values=c("blue","dark red","cornflowerblue","firebrick1","dark grey")) +
  xlab("Length-Normalized 5'UTR MFE") +
  stat_density_ridges(quantile_lines = TRUE, quantiles = 2) +
  theme(axis.title.x = element_text(size=22),axis.text.x  = element_text( vjust=0.5, size=18)) +
  theme(axis.title.y = element_text(size=24),axis.text.y  = element_text( vjust=0.5, size=18))  +
  theme(legend.position="none")

pdf(file = "figures/txclass_vs_MFE_ridges.pdf", width = 5, height = 5)
print(txclass_vs_MFE_ridges)
dev.off()

pdf(file = "figures/txclass_vs_normMFE_ridges.pdf", width = 5, height = 5)
print(txclass_vs_normMFE_ridges)
dev.off()

DE_5utr_feats_withMFE <- five_utr_length | five_utr_gc | txclass_vs_normMFE_ridges

temp2 <- features_df_completecases[which(features_df_completecases$tx_type == "Not_significant"),]
temp2 <- temp2$len_5ut
temp <- features_df_completecases[which(features_df_completecases$tx_type == "TE_down"),]
temp <- temp$len_5ut
len_5ut_clifd <- cliff.delta(temp, temp2, return.dm = TRUE)

temp2 <- features_df_completecases[which(features_df_completecases$tx_type == "Not_significant"),]
temp2 <- temp2$G.C
temp <- features_df_completecases[which(features_df_completecases$tx_type == "TE_down"),]
temp <- temp$G.C
G.C_clifd <- cliff.delta(temp, temp2, return.dm = TRUE)

temp2 <- features_df_completecases[which(features_df_completecases$tx_type == "Not_significant"),]
temp2 <- temp2$normMFE
temp <- features_df_completecases[which(features_df_completecases$tx_type == "TE_down"),]
temp <- temp$normMFE
normMFE_clifd <- cliff.delta(temp, temp2, return.dm = TRUE)

pdf(file = "figures/DE_5utr_feats_withMFE.pdf",width = 17.5,height = 5)
print(DE_5utr_feats_withMFE)
dev.off()

uAUG_number <- ggplot(features_df_completecases,aes(ATG_number,group=tx_type,y=tx_type,fill=tx_type))  +
  theme_ridges() +
  scale_fill_manual(values=c("blue","dark red","cornflowerblue","firebrick1","dark grey")) +
  xlab("Number of uAUGs") +
  stat_density_ridges(quantile_lines = TRUE, quantiles = 2,from = 0) +
  theme(legend.position="none") +
  theme(axis.title.x = element_text(size=22),axis.text.x  = element_text( vjust=0.5, size=18)) +
  theme(axis.title.y = element_text(size=24),axis.text.y  = element_text( vjust=0.5, size=18)) +
  scale_x_continuous(limits = c(0,10))

uAUG_density <- ggplot(features_df_completecases,aes(ATG_density,group=tx_type,y=tx_type,fill=tx_type))  +
  theme_ridges() +
  scale_fill_manual(values=c("blue","dark red","cornflowerblue","firebrick1","dark grey")) +
  xlab("Density of uAUGs") +
  stat_density_ridges(quantile_lines = TRUE, quantiles = 2,from = 0) +
  theme(legend.position="none") +
  theme(axis.title.x = element_text(size=22),axis.text.x  = element_text( vjust=0.5, size=18)) +
  theme(axis.title.y = element_text(size=24),axis.text.y  = element_text( vjust=0.5, size=18)) +
  scale_x_continuous(limits = c(0,0.05))

longest_cognate_uORF <- ggplot(features_df_completecases,aes(log(len_longest_uORF+1,2),group=tx_type,y=tx_type,fill=tx_type))  +
  theme_ridges() +
  scale_fill_manual(values=c("blue","dark red","cornflowerblue","firebrick1","dark grey")) +
  xlab("Log2 longest cognate-uORF") +
  stat_density_ridges(quantile_lines = TRUE, quantiles = 2,from = 0) +
  theme(legend.position="none") +
  theme(axis.title.x = element_text(size=22),axis.text.x  = element_text( vjust=0.5, size=18)) +
  theme(axis.title.y = element_text(size=24),axis.text.y  = element_text( vjust=0.5, size=18)) 

G4_number <- ggplot(features_df_completecases,aes(number_GQuad,group=tx_type,y=tx_type,fill=tx_type))  +
  theme_ridges() +
  scale_fill_manual(values=c("blue","dark red","cornflowerblue","firebrick1","dark grey")) +
  xlab("Number of G4s") +
  stat_density_ridges(quantile_lines = TRUE, quantiles = 2,from = 0) +
  theme(legend.position="none") +
  theme(axis.title.x = element_text(size=22),axis.text.x  = element_text( vjust=0.5, size=18)) +
  theme(axis.title.y = element_text(size=24),axis.text.y  = element_text( vjust=0.5, size=18)) +
  scale_x_continuous(limits = c(0,20), breaks = c(0,5,10,20), labels = c(0,5,10,20))

G4_score_cumu <- ggplot(features_df_completecases,aes(score_GQuad,group=tx_type,y=tx_type,fill=tx_type))  +
  theme_ridges() +
  scale_fill_manual(values=c("blue","dark red","cornflowerblue","firebrick1","dark grey")) +
  xlab("Cumulative G4-propensity") +
  stat_density_ridges(quantile_lines = TRUE, quantiles = 2,from = 0) +
  theme(legend.position="none") +
  theme(axis.title.x = element_text(size=22),axis.text.x  = element_text( vjust=0.5, size=18)) +
  theme(axis.title.y = element_text(size=24),axis.text.y  = element_text( vjust=0.5, size=18)) +
  scale_x_continuous(limits = c(0,1000), breaks = c(0,250,500,750,1000), labels = c(0,250,500,750,1000))

G4_density <- ggplot(features_df_completecases,aes((GQuad_density),group=tx_type,y=tx_type,fill=tx_type))  +
  theme_ridges() +
  scale_fill_manual(values=c("blue","dark red","cornflowerblue","firebrick1","dark grey")) +
  xlab("G4-density") +
  stat_density_ridges(quantile_lines = TRUE, quantiles = 2,from = 0) +
  theme(legend.position="none") +
  theme(axis.title.x = element_text(size=22),axis.text.x  = element_text( vjust=0.5, size=18)) +
  theme(axis.title.y = element_text(size=24),axis.text.y  = element_text( vjust=0.5, size=18)) +
  scale_x_continuous(limits = c(0,0.06), breaks = c(0,0.02,0.04,0.06), labels = c(0,0.02,0.04,0.06))

features_df_completecases$GQuad_density_corr <- features_df_completecases$GQuad_density / features_df_completecases$G.C

G4_density_GCadj <- ggplot(features_df_completecases,aes((GQuad_density_corr*100),group=tx_type,y=tx_type,fill=tx_type))  +
  theme_ridges() +
  scale_fill_manual(values=c("blue","dark red","cornflowerblue","firebrick1","dark grey")) +
  xlab("G4-density, GC-content adjusted") +
  stat_density_ridges(quantile_lines = TRUE, quantiles = 2,from = 0) +
  theme(legend.position="none") +
  theme(axis.title.x = element_text(size=22),axis.text.x  = element_text( vjust=0.5, size=18)) +
  theme(axis.title.y = element_text(size=24),axis.text.y  = element_text( vjust=0.5, size=18)) +
  scale_x_continuous(limits = c(0,0.04), breaks = c(0,0.01,0.02,0.03,0.04), labels = c(0,0.01,0.02,0.03,0.04))

G4_score_average <- ggplot(features_df_completecases,aes(GQuad_av_strength,group=tx_type,y=tx_type,fill=tx_type))  +
  theme_ridges() +
  scale_fill_manual(values=c("blue","dark red","cornflowerblue","firebrick1","dark grey")) +
  xlab("G4-propensity per instance") +
  stat_density_ridges(quantile_lines = TRUE, quantiles = 2,from = 0) +
  theme(legend.position="none") +
  theme(axis.title.x = element_text(size=22),axis.text.x  = element_text( vjust=0.5, size=18)) +
  theme(axis.title.y = element_text(size=24),axis.text.y  = element_text( vjust=0.5, size=18)) +
  scale_x_continuous(limits = c(0,125), breaks = c(0,25,50,75,100,125), labels = c(0,25,50,75,100,125))

uORF_plots <- uAUG_number | uAUG_density | longest_cognate_uORF
GQuad_fiveUTR_plots <- (G4_number | G4_density | G4_density_GCadj) /
  (G4_score_cumu | G4_score_average)
GQuad_fiveUTR_plots_withMFE <- (txclass_vs_MFE_ridges | G4_number | G4_density ) /
  (G4_density_GCadj | G4_score_cumu | G4_score_average)

pdf(file = "figures/GQuad_fiveUTR_plots_withMFE.pdf",width = 17.5,height = 10)
print(GQuad_fiveUTR_plots_withMFE)
dev.off()

pdf(file = "figures/GQuad_fiveUTR_plots.pdf",width = 17.5,height = 10)
print(GQuad_fiveUTR_plots)
dev.off()

pdf(file = "figures/uORF_plots.pdf",width = 17.5,height = 5)
print(uORF_plots)
dev.off()

pdf(file = "figures/uAUG_number.pdf", width = 5, height = 5)
print(uAUG_number)
dev.off()

pdf(file = "figures/uAUG_density.pdf", width = 5, height = 5)
print(uAUG_density)
dev.off()

pdf(file = "figures/longest_cognate_uORF.pdf", width = 5, height = 5)
print(longest_cognate_uORF)
dev.off()

pdf(file = "figures/G4_number.pdf", width = 5, height = 5)
print(G4_number)
dev.off()

pdf(file = "figures/G4_density.pdf", width = 5, height = 5)
print(G4_density)
dev.off()

pdf(file = "figures/G4_density_GCadj.pdf", width = 5, height = 5)
print(G4_density_GCadj)
dev.off()

pdf(file = "figures/G4_score_cumu.pdf", width = 5, height = 5)
print(G4_score_cumu)
dev.off()

pdf(file = "figures/G4_score_average.pdf", width = 5, height = 5)
print(G4_score_average)
dev.off()

###################################################################################################

load("./counts/RP_ctrl_rep1_Aligned.sortedByCoord.out.bam_counts_summary_25min")
control_ribo_rep1 <- res_counts$counts_gene_ex_unq
load("./counts/RP_ctrl_rep2_Aligned.sortedByCoord.out.bam_counts_summary_25min")
control_ribo_rep2 <- res_counts$counts_gene_ex_unq
load("./counts/RP_siDDX3_rep1_Aligned.sortedByCoord.out.bam_counts_summary_25min")
siDDX3X_ribo_rep1 <- res_counts$counts_gene_ex_unq
load("./counts/RP_siDDX3_rep2_Aligned.sortedByCoord.out.bam_counts_summary_25min")
siDDX3X_ribo_rep2 <- res_counts$counts_gene_ex_unq
load("./counts/RNA_ctrl_rep1_Aligned.sortedByCoord.out.bam_counts_summary_25min")
control_rna_rep1 <- res_counts$counts_gene_ex_unq
load("./counts/RNA_ctrl_rep2_Aligned.sortedByCoord.out.bam_counts_summary_25min")
control_rna_rep2 <- res_counts$counts_gene_ex_unq
load("./counts/RNA_siDDX3_rep1_Aligned.sortedByCoord.out.bam_counts_summary_25min")
siDDX3X_rna_rep1 <- res_counts$counts_gene_ex_unq
load("./counts/RNA_siDDX3_rep2_Aligned.sortedByCoord.out.bam_counts_summary_25min")
siDDX3X_rna_rep2 <- res_counts$counts_gene_ex_unq
rm(res_counts)

siDDX3X_markus_gene_ex_unq <- cbind(control_ribo_rep1, 
                                    control_ribo_rep2, 
                                    siDDX3X_ribo_rep1, 
                                    siDDX3X_ribo_rep2, 
                                    control_rna_rep1, 
                                    control_rna_rep2, 
                                    siDDX3X_rna_rep1, 
                                    siDDX3X_rna_rep2)

colnames(siDDX3X_markus_gene_ex_unq) <- c("control_ribo_rep1", 
                                          "control_ribo_rep2", 
                                          "siDDX3X_ribo_rep1", 
                                          "siDDX3X_ribo_rep2", 
                                          "control_rna_rep1", 
                                          "control_rna_rep2", 
                                          "siDDX3X_rna_rep1", 
                                          "siDDX3X_rna_rep2")

siDDX3X_markus_gene_ex_unq <- as.data.frame(siDDX3X_markus_gene_ex_unq)

rm(list = ls(pattern = "rep"))

exp_design <- read.delim("./sample_info_siDDX3X_markus.txt")
exp_design$Genotype <- factor(exp_design$Genotype, levels = c("control","siDDX3X"))
exp_design$Condition <- factor(exp_design$Condition, levels = c( "rna","ribo"))

siDDX3X_markus_deseq_dataset_ex <- DESeqDataSetFromMatrix(siDDX3X_markus_gene_ex_unq, exp_design, design = ~ Condition + Genotype + Condition:Genotype)
siDDX3X_markus_deseq_dataset_ex <- DESeq(siDDX3X_markus_deseq_dataset_ex, test = "LRT", reduced = ~ Condition + Genotype)
siDDX3X_markus_deseq_results_ex <- results(siDDX3X_markus_deseq_dataset_ex, independentFiltering = FALSE)
siDDX3X_markus_deseq_results_ex$pvalue[siDDX3X_markus_deseq_results_ex$baseMean < 20] <- NA
siDDX3X_markus_deseq_results_ex$padj <- p.adjust(siDDX3X_markus_deseq_results_ex$pvalue, method = "BH")

exp_design_rnaonly <- exp_design[c(5:8),]
exp_design_riboonly <- exp_design[c(1:4),]
siDDX3X_markus_gene_ex_unq_rnaonly <- siDDX3X_markus_gene_ex_unq[,c(5:8)]
siDDX3X_markus_gene_ex_unq_riboonly <- siDDX3X_markus_gene_ex_unq[,c(1:4)]

siDDX3X_markus_deseq_dataset_ex_rnaonly <- DESeqDataSetFromMatrix(siDDX3X_markus_gene_ex_unq_rnaonly, exp_design_rnaonly, design = ~ Genotype)
siDDX3X_markus_deseq_dataset_ex_rnaonly <- DESeq(siDDX3X_markus_deseq_dataset_ex_rnaonly)
siDDX3X_markus_deseq_results_ex_rnaonly <- results(siDDX3X_markus_deseq_dataset_ex_rnaonly, independentFiltering = FALSE)
siDDX3X_markus_deseq_results_ex_rnaonly$pvalue[siDDX3X_markus_deseq_results_ex_rnaonly$baseMean < 20] <- NA
siDDX3X_markus_deseq_results_ex_rnaonly$padj <- p.adjust(siDDX3X_markus_deseq_results_ex_rnaonly$pvalue, method = "BH")

siDDX3X_markus_deseq_dataset_ex_riboonly <- DESeqDataSetFromMatrix(siDDX3X_markus_gene_ex_unq_riboonly, exp_design_riboonly, design = ~ Genotype)
siDDX3X_markus_deseq_dataset_ex_riboonly <- DESeq(siDDX3X_markus_deseq_dataset_ex_riboonly)
siDDX3X_markus_deseq_results_ex_riboonly <- results(siDDX3X_markus_deseq_dataset_ex_riboonly, independentFiltering = FALSE)
siDDX3X_markus_deseq_results_ex_riboonly$pvalue[siDDX3X_markus_deseq_results_ex_riboonly$baseMean < 20] <- NA
siDDX3X_markus_deseq_results_ex_riboonly$padj <- p.adjust(siDDX3X_markus_deseq_results_ex_riboonly$pvalue, method = "BH")

siDDX3X_markus_deseq_results_ex$symb<-GTF_annotation$trann$gene_name[match(rownames(siDDX3X_markus_deseq_results_ex),GTF_annotation$trann$gene_id)]
siDDX3X_markus_deseq_results_ex$biotype<-GTF_annotation$trann$gene_biotype[match(rownames(siDDX3X_markus_deseq_results_ex),GTF_annotation$trann$gene_id)]
siDDX3X_markus_deseq_results_ex_rnaonly$symb<-GTF_annotation$trann$gene_name[match(rownames(siDDX3X_markus_deseq_results_ex_rnaonly),GTF_annotation$trann$gene_id)]
siDDX3X_markus_deseq_results_ex_rnaonly$biotype<-GTF_annotation$trann$gene_biotype[match(rownames(siDDX3X_markus_deseq_results_ex_rnaonly),GTF_annotation$trann$gene_id)]
siDDX3X_markus_deseq_results_ex_riboonly$symb<-GTF_annotation$trann$gene_name[match(rownames(siDDX3X_markus_deseq_results_ex_riboonly),GTF_annotation$trann$gene_id)]
siDDX3X_markus_deseq_results_ex_riboonly$biotype<-GTF_annotation$trann$gene_biotype[match(rownames(siDDX3X_markus_deseq_results_ex_riboonly),GTF_annotation$trann$gene_id)]

results_ex_siDDX3X_markus <- siDDX3X_markus_deseq_results_ex[,c(1,2,6,7,8)]
colnames(results_ex_siDDX3X_markus) <- c("baseMean_riborna", "log2FoldChange_riborna", "padj_riborna", "symb", "biotype")

results_ex_rnaonly_siDDX3X_markus <- siDDX3X_markus_deseq_results_ex_rnaonly[,c(1,2,6)]
colnames(results_ex_rnaonly_siDDX3X_markus) <- c("baseMean_rna", "log2FoldChange_rna", "padj_rna")

results_ex_riboonly_siDDX3X_markus <- siDDX3X_markus_deseq_results_ex_riboonly[,c(1,2,6)]
colnames(results_ex_riboonly_siDDX3X_markus) <- c("baseMean_ribo", "log2FoldChange_ribo", "padj_ribo")

results_ex_siDDX3X_markus <- cbind(results_ex_siDDX3X_markus, results_ex_rnaonly_siDDX3X_markus, results_ex_riboonly_siDDX3X_markus)
results_ex_siDDX3X_markus <- results_ex_siDDX3X_markus[order(results_ex_siDDX3X_markus$padj_riborna),]

cutoff_ex_siDDX3X_markus<-which(results_ex_siDDX3X_markus$baseMean_rna>20 & results_ex_siDDX3X_markus$baseMean_ribo>20)
results_ex_siDDX3X_markus <- results_ex_siDDX3X_markus[cutoff_ex_siDDX3X_markus,]

results_ex_siDDX3X_markus$tx_type<-"Not_significant"
results_ex_siDDX3X_markus$tx_type[results_ex_siDDX3X_markus$padj_rna<.01 & results_ex_siDDX3X_markus$log2FoldChange_rna<0]<-"RNA_down"
results_ex_siDDX3X_markus$tx_type[results_ex_siDDX3X_markus$padj_rna<.01 & results_ex_siDDX3X_markus$log2FoldChange_rna>0]<-"RNA_up"
results_ex_siDDX3X_markus$tx_type[results_ex_siDDX3X_markus$padj_riborna<.05 & results_ex_siDDX3X_markus$log2FoldChange_riborna<0]<-"TE_down"
results_ex_siDDX3X_markus$tx_type[results_ex_siDDX3X_markus$padj_riborna<.05 & results_ex_siDDX3X_markus$log2FoldChange_riborna>0]<-"TE_up"

results_ex_siDDX3X_markus$log10_padj_riborna<-log10(results_ex_siDDX3X_markus$padj_riborna)
results_ex_siDDX3X_markus$log10_padj_riborna[is.na(results_ex_siDDX3X_markus$log10_padj_riborna)]<-0

results_ex_siDDX3X_markus$tx_type<-factor(results_ex_siDDX3X_markus$tx_type,levels=c("TE_up","TE_down","RNA_up","RNA_down","Not_significant"))
results_ex_siDDX3X_markus$biotype[grep(results_ex_siDDX3X_markus$symb,pattern = "HIST|H2B")]<-"Histone genes"
results_ex_siDDX3X_markus$biotype[grep(results_ex_siDDX3X_markus$symb,pattern = "HIST|H2B",invert =  T)]<-"all genes"
results_ex_siDDX3X_markus$biotype<-factor(results_ex_siDDX3X_markus$biotype,levels=c("all genes","Histone genes"))

results_ex_siDDX3X_markus$symb_deplot<-as.character(results_ex_siDDX3X_markus$symb)
results_ex_siDDX3X_markus$symb_deplot<-NA
toprna_ex_siDDX3X_markus<-results_ex_siDDX3X_markus[results_ex_siDDX3X_markus$tx_type=="RNA_up",]
toprna_ex_siDDX3X_markus<-toprna_ex_siDDX3X_markus[order(toprna_ex_siDDX3X_markus$log2FoldChange_rna,toprna_ex_siDDX3X_markus$log2FoldChange_ribo,decreasing = T),"symb"]
if(length(toprna_ex_siDDX3X_markus)>6){toprna_ex_siDDX3X_markus=toprna_ex_siDDX3X_markus[1:5]}
topribo_ex_siDDX3X_markus<-results_ex_siDDX3X_markus[results_ex_siDDX3X_markus$tx_type=="TE_up",]
topribo_ex_siDDX3X_markus<-topribo_ex_siDDX3X_markus[order(topribo_ex_siDDX3X_markus$log2FoldChange_riborna,decreasing = T),"symb"]
if(length(topribo_ex_siDDX3X_markus)>6){topribo_ex_siDDX3X_markus=topribo_ex_siDDX3X_markus[1:4]}
bottomrna_ex_siDDX3X_markus<-results_ex_siDDX3X_markus[results_ex_siDDX3X_markus$tx_type=="RNA_down",]
bottomrna_ex_siDDX3X_markus<-bottomrna_ex_siDDX3X_markus[order(bottomrna_ex_siDDX3X_markus$log2FoldChange_rna,bottomrna_ex_siDDX3X_markus$log2FoldChange_ribo,decreasing = F),"symb"]
if(length(bottomrna_ex_siDDX3X_markus)>6){bottomrna_ex_siDDX3X_markus=bottomrna_ex_siDDX3X_markus[1:4]}
bottomribo_ex_siDDX3X_markus<-results_ex_siDDX3X_markus[results_ex_siDDX3X_markus$tx_type=="TE_down",]
bottomribo_ex_siDDX3X_markus<-bottomribo_ex_siDDX3X_markus[order(bottomribo_ex_siDDX3X_markus$log2FoldChange_riborna,decreasing = F),"symb"]
if(length(bottomribo_ex_siDDX3X_markus)>10){bottomribo_ex_siDDX3X_markus=c(bottomribo_ex_siDDX3X_markus[1:4],bottomribo_ex_siDDX3X_markus[bottomribo_ex_siDDX3X_markus%in%c("ODC1","DVL1","DVL2","HMBS","PRKRA")])}
results_ex_siDDX3X_markus$symb_deplot[results_ex_siDDX3X_markus$symb%in%c(toprna_ex_siDDX3X_markus,topribo_ex_siDDX3X_markus,bottomrna_ex_siDDX3X_markus,bottomribo_ex_siDDX3X_markus)]<-results_ex_siDDX3X_markus$symb[results_ex_siDDX3X_markus$symb%in%c(toprna_ex_siDDX3X_markus,topribo_ex_siDDX3X_markus,bottomrna_ex_siDDX3X_markus,bottomribo_ex_siDDX3X_markus)]

results_ex_siDDX3X_markus$nudgey<-NA
results_ex_siDDX3X_markus$nudgey[which(nchar(results_ex_siDDX3X_markus$symb_deplot)>2)]<-ifelse(results_ex_siDDX3X_markus$log2FoldChange_ribo[which(nchar(results_ex_siDDX3X_markus$symb_deplot)>2)]>0, .4, -.5)
results_ex_siDDX3X_markus$nudgex<-NA
results_ex_siDDX3X_markus$nudgex[which(nchar(results_ex_siDDX3X_markus$symb_deplot)>2)]<-ifelse(results_ex_siDDX3X_markus$log2FoldChange_rna[which(nchar(results_ex_siDDX3X_markus$symb_deplot)>2)]>0, .3, -.4)
results_ex_siDDX3X_markus$nudgex[results_ex_siDDX3X_markus$tx_type=="TE_up"]<-0
results_ex_siDDX3X_markus$nudgey[results_ex_siDDX3X_markus$tx_type=="RNA_up"]<-0
results_ex_siDDX3X_markus$nudgey[results_ex_siDDX3X_markus$tx_type=="RNA_down"]<-0
results_ex_siDDX3X_markus$nudgex[results_ex_siDDX3X_markus$tx_type=="RNA_up"]<-0
results_ex_siDDX3X_markus$nudgex[results_ex_siDDX3X_markus$tx_type=="RNA_down"]<-0

results_ex_siDDX3X_markus_df <- as.data.frame(results_ex_siDDX3X_markus)

write.csv(results_ex_siDDX3X_markus_df, "data/results_ex_siDDX3X_markus_df.csv")

results_ex_siDDX3X_markus_df$nudgey[results_ex_siDDX3X_markus_df$symb == "ODC1"] <- 0.5
DE_plot_ex_siDDX3X_markus <- ggplot(results_ex_siDDX3X_markus_df,aes(x=log2FoldChange_rna,y=log2FoldChange_ribo,color=tx_type,size=-log10_padj_riborna,shape=biotype,label=symb_deplot)) + 
  geom_point() + 
  theme_bw() +
  ylab("Ribo fold change (siDDX3X/control)") +
  xlab("RNA fold change (siDDX3X/control)") +
  geom_vline(xintercept = 0, alpha = 0.5) +
  geom_hline(yintercept = 0, alpha = 0.5) +
  theme(axis.title.x = element_text(size=22),axis.text.x  = element_text(angle=45, vjust=0.5, size=18)) +
  theme(axis.title.y = element_text(size=24),axis.text.y  = element_text(angle=45, vjust=0.5, size=18))  +
  scale_color_manual(values = alpha(c("dark red","blue","firebrick1","cornflowerblue","gray24"),c(.8,.8,.8,.8,.3)),"Tx_class") +
  scale_size_continuous(name = "-log10 adj.pval\nRiboRNA") +
  scale_x_continuous(limits = c(-3,3)) + 
  scale_y_continuous(limits = c(-3,3)) +
  geom_text_repel(size=5,nudge_x = results_ex_siDDX3X_markus_df$nudgex[which(nchar(results_ex_siDDX3X_markus_df$symb_deplot)>2)],
                  nudge_y = results_ex_siDDX3X_markus_df$nudgey[which(nchar(results_ex_siDDX3X_markus_df$symb_deplot)>2)],force=20)
pdf(file = "figures/DE_plot_ex_siDDX3X_markus.pdf")
print(DE_plot_ex_siDDX3X_markus)
dev.off()

results_ex_25min_df$log2FCTE <- results_ex_25min_df$log2FoldChange_ribo - results_ex_25min_df$log2FoldChange_rna
results_ex_siDDX3X_markus_df$log2FCTE <- results_ex_siDDX3X_markus_df$log2FoldChange_ribo - results_ex_siDDX3X_markus_df$log2FoldChange_rna

merge_ddx3xdegron_siDDX3X <- merge(results_ex_25min_df,results_ex_siDDX3X_markus_df,by.x="symb",by.y="symb",all.x=TRUE, all.y=TRUE)
cor<-cor.test(merge_ddx3xdegron_siDDX3X$log2FCTE.x, merge_ddx3xdegron_siDDX3X$log2FCTE.y, method = "pearson")
cor2<-cor.test(merge_ddx3xdegron_siDDX3X$log2FCTE.x, merge_ddx3xdegron_siDDX3X$log2FCTE.y, method = "spearman")
grob = grobTree(textGrob(paste("Pearson Correlation : ", round(cor$estimate, 2) ), x = 0.5, y = 0.1, hjust = 0, gp = gpar(col = "red", fontsize = 15)))
grob2 = grobTree(textGrob(paste("Spearman Correlation : ", round(cor2$estimate, 2) ), x = 0.5, y = 0.05, hjust = 0, gp = gpar(col = "red", fontsize = 15)))
TEcorr_degron_siDDX3X <- ggplot(merge_ddx3xdegron_siDDX3X,aes(x=log2FCTE.x,y=log2FCTE.y)) + 
  geom_point(alpha=0.3) + 
  stat_smooth(method =lm, geom='line', alpha=0.5, se=FALSE, color = "darkred", size = 2) +
  theme_bw() +
  xlab("TE fold change DDX3X (IAA/DMSO)") +
  ylab("TE fold change siDDX3X/siNTC") +
  theme(axis.title.x = element_text(size=24),axis.text.x  = element_text(angle=45, vjust=0.5, size=16)) +
  theme(axis.title.y = element_text(size=24),axis.text.y  = element_text(angle=45, vjust=0.5, size=16))  +
  scale_x_continuous(limits = c(-2,2)) + 
  scale_y_continuous(limits = c(-2,2)) +
  annotation_custom(grob) +
  annotation_custom(grob2)

pdf(file = "figures/TE_compare_siDDX3XvsDegron.pdf", width = 8, height = 8)
print(TEcorr_degron_siDDX3X)
dev.off()

png(file = "figures/TE_compare_siDDX3XvsDegron.png", width = 8, height = 8, units = "in", res = 600)
print(TEcorr_degron_siDDX3X)
dev.off()

write.csv(merge_ddx3xdegron_siDDX3X, "data/merge_ddx3xdegron_siDDX3X_cds.csv")


############################################################################################################################################

load("./counts/IAA3XWT_HCT116_rep1_ribo_Aligned.sortedByCoord.out.bam_counts_summary_25min")
ddx3x_ribo_rep1_25min <- res_counts$counts_gene_ex_unq
load("./counts/IAA3XWT_HCT116_rep2_ribo_Aligned.sortedByCoord.out.bam_counts_summary_25min")
ddx3x_ribo_rep2_25min <- res_counts$counts_gene_ex_unq
load("./counts/IAA3Y_HCT116_rep1_ribo_Aligned.sortedByCoord.out.bam_counts_summary_25min")
ddx3y_ribo_rep1_25min <- res_counts$counts_gene_ex_unq
load("./counts/IAA3Y_HCT116_rep2_ribo_Aligned.sortedByCoord.out.bam_counts_summary_25min")
ddx3y_ribo_rep2_25min <- res_counts$counts_gene_ex_unq
load("./counts/IAA3XWT_HCT116_rep1_rna_Aligned.sortedByCoord.out.bam_counts_summary_25min")
ddx3x_rna_rep1_25min <- res_counts$counts_gene_ex_unq
load("./counts/IAA3XWT_HCT116_rep2_rna_Aligned.sortedByCoord.out.bam_counts_summary_25min")
ddx3x_rna_rep2_25min <- res_counts$counts_gene_ex_unq
load("./counts/IAA3Y_HCT116_rep1_rna_Aligned.sortedByCoord.out.bam_counts_summary_25min")
ddx3y_rna_rep1_25min <- res_counts$counts_gene_ex_unq
load("./counts/IAA3Y_HCT116_rep2_rna_Aligned.sortedByCoord.out.bam_counts_summary_25min")
ddx3y_rna_rep2_25min <- res_counts$counts_gene_ex_unq
rm(res_counts)

gene_ex_unq_25min <- cbind(ddx3x_ribo_rep1_25min, 
                           ddx3x_ribo_rep2_25min, 
                           ddx3y_ribo_rep1_25min, 
                           ddx3y_ribo_rep2_25min, 
                           ddx3x_rna_rep1_25min, 
                           ddx3x_rna_rep2_25min, 
                           ddx3y_rna_rep1_25min, 
                           ddx3y_rna_rep2_25min)
colnames(gene_ex_unq_25min) <- c("ddx3x_ribo_rep1_25min", 
                                 "ddx3x_ribo_rep2_25min", 
                                 "ddx3y_ribo_rep1_25min", 
                                 "ddx3y_ribo_rep2_25min", 
                                 "ddx3x_rna_rep1_25min", 
                                 "ddx3x_rna_rep2_25min", 
                                 "ddx3y_rna_rep1_25min", 
                                 "ddx3y_rna_rep2_25min")
gene_ex_unq_25min <- as.data.frame(gene_ex_unq_25min)

rm(list = ls(pattern = "rep"))

#loads exp_design file and sets parameter values
exp_design <- read.delim("./sample_info.txt")
exp_design$Genotype <- factor(exp_design$Genotype, levels = c("ddx3x","ddx3y"))
exp_design$Condition <- factor(exp_design$Condition, levels = c( "rna","ribo"))

#DESeq for the riboRNA interaction
deseq_dataset_ex_25min <- DESeqDataSetFromMatrix(gene_ex_unq_25min, exp_design, design = ~ Condition + Genotype + Condition:Genotype)
deseq_dataset_ex_25min <- DESeq(deseq_dataset_ex_25min, test = "LRT", reduced = ~ Condition + Genotype)
deseq_results_ex_25min <- results(deseq_dataset_ex_25min)
deseq_results_ex_25min$pvalue[deseq_results_ex_25min$baseMean < 20] <- NA
deseq_results_ex_25min$padj <- p.adjust(deseq_results_ex_25min$pvalue, method = "BH")

#subsets the exp design and count files for rna only and ribo only DE
exp_design_rnaonly <- exp_design[c(5:8),]
exp_design_riboonly <- exp_design[c(1:4),]
gene_ex_unq_25min_rnaonly <- gene_ex_unq_25min[,c(5:8)]
gene_ex_unq_25min_riboonly <- gene_ex_unq_25min[,c(1:4)]

#DESeq for the individual rna and ribo terms
deseq_dataset_ex_25min_rnaonly <- DESeqDataSetFromMatrix(gene_ex_unq_25min_rnaonly, exp_design_rnaonly, design = ~ Genotype)
deseq_dataset_ex_25min_rnaonly <- DESeq(deseq_dataset_ex_25min_rnaonly)
deseq_results_ex_25min_rnaonly <- results(deseq_dataset_ex_25min_rnaonly)
deseq_results_ex_25min_rnaonly$pvalue[deseq_results_ex_25min_rnaonly$baseMean < 20] <- NA
deseq_results_ex_25min_rnaonly$padj <- p.adjust(deseq_results_ex_25min_rnaonly$pvalue, method = "BH")

deseq_dataset_ex_25min_riboonly <- DESeqDataSetFromMatrix(gene_ex_unq_25min_riboonly, exp_design_riboonly, design = ~ Genotype)
deseq_dataset_ex_25min_riboonly <- DESeq(deseq_dataset_ex_25min_riboonly)
deseq_results_ex_25min_riboonly <- results(deseq_dataset_ex_25min_riboonly)
deseq_results_ex_25min_riboonly$pvalue[deseq_results_ex_25min_riboonly$baseMean < 20] <- NA
deseq_results_ex_25min_riboonly$padj <- p.adjust(deseq_results_ex_25min_riboonly$pvalue, method = "BH")

#calls annotations from the GTF. Gene name, biotype 
deseq_results_ex_25min$symb<-GTF_annotation$trann$gene_name[match(rownames(deseq_results_ex_25min),GTF_annotation$trann$gene_id)]
deseq_results_ex_25min$biotype<-GTF_annotation$trann$gene_biotype[match(rownames(deseq_results_ex_25min),GTF_annotation$trann$gene_id)]
deseq_results_ex_25min_rnaonly$symb<-GTF_annotation$trann$gene_name[match(rownames(deseq_results_ex_25min_rnaonly),GTF_annotation$trann$gene_id)]
deseq_results_ex_25min_rnaonly$biotype<-GTF_annotation$trann$gene_biotype[match(rownames(deseq_results_ex_25min_rnaonly),GTF_annotation$trann$gene_id)]
deseq_results_ex_25min_riboonly$symb<-GTF_annotation$trann$gene_name[match(rownames(deseq_results_ex_25min_riboonly),GTF_annotation$trann$gene_id)]
deseq_results_ex_25min_riboonly$biotype<-GTF_annotation$trann$gene_biotype[match(rownames(deseq_results_ex_25min_riboonly),GTF_annotation$trann$gene_id)]

#builds a common results table with relevant parameters
results_ex_25min <- deseq_results_ex_25min[,c(1,2,6,7,8)]
colnames(results_ex_25min) <- c("baseMean_riborna", "log2FoldChange_riborna", "padj_riborna","symb", "biotype")

results_ex_25min_rnaonly <- deseq_results_ex_25min_rnaonly[,c(1,2,6)]
colnames(results_ex_25min_rnaonly) <- c("baseMean_rna", "log2FoldChange_rna", "padj_rna")

results_ex_25min_riboonly <- deseq_results_ex_25min_riboonly[,c(1,2,6)]
colnames(results_ex_25min_riboonly) <- c("baseMean_ribo", "log2FoldChange_ribo", "padj_ribo")

results_ex_25min <- cbind(results_ex_25min, results_ex_25min_rnaonly, results_ex_25min_riboonly)
results_ex_25min <- results_ex_25min[order(results_ex_25min$padj_riborna),]

#imposes read number cutoffs
cutoff_ex<-which(results_ex_25min$baseMean_rna>20 & results_ex_25min$baseMean_ribo>20)
results_ex_25min <- results_ex_25min[cutoff_ex,]

# classifies genes based on the effect of perturbation
results_ex_25min$tx_type<-"Not_significant"
results_ex_25min$tx_type[results_ex_25min$padj_rna<.01 & results_ex_25min$log2FoldChange_rna<0]<-"RNA_down"
results_ex_25min$tx_type[results_ex_25min$padj_rna<.01 & results_ex_25min$log2FoldChange_rna>0]<-"RNA_up"
results_ex_25min$tx_type[results_ex_25min$padj_riborna<.05 & results_ex_25min$log2FoldChange_riborna<0]<-"TE_down"
results_ex_25min$tx_type[results_ex_25min$padj_riborna<.05 & results_ex_25min$log2FoldChange_riborna>0]<-"TE_up"

#log of p-value
results_ex_25min$log10_padj_riborna<-log10(results_ex_25min$padj_riborna)
results_ex_25min$log10_padj_riborna[is.na(results_ex_25min$log10_padj_riborna)]<-0

results_ex_25min$tx_type<-factor(results_ex_25min$tx_type,levels=c("TE_up","TE_down","RNA_up","RNA_down","Not_significant"))
results_ex_25min$biotype[grep(results_ex_25min$symb,pattern = "HIST|H2B")]<-"Histone genes"
results_ex_25min$biotype[grep(results_ex_25min$symb,pattern = "HIST|H2B",invert =  T)]<-"all genes"
results_ex_25min$biotype<-factor(results_ex_25min$biotype,levels=c("all genes","Histone genes"))

results_ex_25min$symb_deplot<-as.character(results_ex_25min$symb)
results_ex_25min$symb_deplot<-NA
toprna_ex<-results_ex_25min[results_ex_25min$tx_type=="RNA_up",]
toprna_ex<-toprna_ex[order(toprna_ex$log2FoldChange_rna,toprna_ex$log2FoldChange_ribo,decreasing = T),"symb"]
if(length(toprna_ex)>5){toprna_ex=toprna_ex[1:5]}
topribo_ex<-results_ex_25min[results_ex_25min$tx_type=="TE_up",]
topribo_ex<-topribo_ex[order(topribo_ex$log2FoldChange_riborna,decreasing = T),"symb"]
if(length(topribo_ex)>3){topribo_ex=topribo_ex[1:2]}
bottomrna_ex<-results_ex_25min[results_ex_25min$tx_type=="RNA_down",]
bottomrna_ex<-bottomrna_ex[order(bottomrna_ex$log2FoldChange_rna,bottomrna_ex$log2FoldChange_ribo,decreasing = F),"symb"]
if(length(bottomrna_ex)>3){bottomrna_ex=bottomrna_ex[1:2]}
bottomribo_ex<-results_ex_25min[results_ex_25min$tx_type=="TE_down",]
bottomribo_ex<-bottomribo_ex[order(bottomribo_ex$log2FoldChange_riborna,decreasing = F),"symb"]
if(length(bottomribo_ex)>3){bottomribo_ex=c(bottomribo_ex[1:3],bottomribo_ex[bottomribo_ex%in%c("ODC1","DVL1","DVL2","HMBS","PRKRA")])}
results_ex_25min$symb_deplot[results_ex_25min$symb%in%c(toprna_ex,topribo_ex,bottomrna_ex,bottomribo_ex)]<-results_ex_25min$symb[results_ex_25min$symb%in%c(toprna_ex,topribo_ex,bottomrna_ex,bottomribo_ex)]

results_ex_25min$nudgey<-NA
results_ex_25min$nudgey[which(nchar(results_ex_25min$symb_deplot)>2)]<-ifelse(results_ex_25min$log2FoldChange_ribo[which(nchar(results_ex_25min$symb_deplot)>2)]>0, .4, -.5)
results_ex_25min$nudgex<-NA
results_ex_25min$nudgex[which(nchar(results_ex_25min$symb_deplot)>2)]<-ifelse(results_ex_25min$log2FoldChange_rna[which(nchar(results_ex_25min$symb_deplot)>2)]>0, .3, -.4)
results_ex_25min$nudgex[results_ex_25min$tx_type=="TE_up"]<-0
results_ex_25min$nudgey[results_ex_25min$tx_type=="RNA_up"]<-0
results_ex_25min$nudgey[results_ex_25min$tx_type=="RNA_down"]<-0
results_ex_25min$nudgex[results_ex_25min$tx_type=="RNA_up"]<-0
results_ex_25min$nudgex[results_ex_25min$tx_type=="RNA_down"]<-0

results_ex_25min_df <- as.data.frame(results_ex_25min)

#flattens log10(padj) to a maximum value of 15. Aesthetic purposes only, the huge p-value for DDX3Y was skewing the point sizes in the plot
results_ex_25min_df$log10_padj_riborna_flat <- results_ex_25min_df$log10_padj_riborna
results_ex_25min_df$log10_padj_riborna_flat[results_ex_25min_df$log10_padj_riborna < -15]<- -15

write.csv(results_ex_25min_df, "data/DE_table_3X3Y_ex_3X3Y.csv")

DE_plot_ex_25min <- ggplot(results_ex_25min_df,aes(x=log2FoldChange_rna,y=log2FoldChange_ribo,color=tx_type,size=-log10_padj_riborna_flat,label=symb_deplot)) + 
  geom_point() + 
  theme_bw() +
  ylab("Ribo fold change (DDX3Y/DDX3X)") +
  xlab("RNA fold change (DDX3Y/DDX3X)") +
  geom_vline(xintercept = 0, alpha = 0.5) +
  geom_hline(yintercept = 0, alpha = 0.5) +
  theme(axis.title.x = element_text(size=22),axis.text.x  = element_text(angle=45, vjust=0.5, size=18)) +
  theme(axis.title.y = element_text(size=24),axis.text.y  = element_text(angle=45, vjust=0.5, size=18))  +
  scale_color_manual(values = alpha(c("dark red","blue","firebrick1","cornflowerblue","gray24"),c(.8,.8,.8,.8,.3)),"Tx_class") +
  scale_size_continuous(name = "-log10 adj.pval\nRiboRNA") +
  scale_x_continuous(limits = c(-15,15)) + 
  scale_y_continuous(limits = c(-10,10)) +
  geom_text_repel(size=5,nudge_x = results_ex_25min_df$nudgex[which(nchar(results_ex_25min_df$symb_deplot)>2)],
                  nudge_y = results_ex_25min_df$nudgey[which(nchar(results_ex_25min_df$symb_deplot)>2)],force=20) +
  theme(legend.text = element_text(size = rel(1.5)), legend.title = element_text(size = rel(1.5)))
pdf(file = "figures/DE_plot_ex_25min.pdf", width = 10, height = 7)
print(DE_plot_ex_25min)
dev.off()

DE_plot_ex_25min_trunc <- ggplot(results_ex_25min_df,aes(x=log2FoldChange_rna,y=log2FoldChange_ribo,color=tx_type,size=-log10_padj_riborna_flat,label=symb_deplot)) + 
  geom_point() + 
  theme_bw() +
  ylab("Ribo fold change (DDX3Y/DDX3X)") +
  xlab("RNA fold change (DDX3Y/DDX3X)") +
  geom_vline(xintercept = 0, alpha = 0.5) +
  geom_hline(yintercept = 0, alpha = 0.5) +
  theme(axis.title.x = element_text(size=22),axis.text.x  = element_text(angle=45, vjust=0.5, size=18)) +
  theme(axis.title.y = element_text(size=24),axis.text.y  = element_text(angle=45, vjust=0.5, size=18))  +
  scale_color_manual(values = alpha(c("dark red","blue","firebrick1","cornflowerblue","gray24"),c(.8,.8,.8,.8,.3)),"Tx_class") +
  scale_size_continuous(name = "-log10 adj.pval\nRiboRNA") +
  scale_x_continuous(limits = c(-2,2)) + 
  scale_y_continuous(limits = c(-2,2)) +
  geom_text_repel(size=5,nudge_x = results_ex_25min_df$nudgex[which(nchar(results_ex_25min_df$symb_deplot)>2)],
                  nudge_y = results_ex_25min_df$nudgey[which(nchar(results_ex_25min_df$symb_deplot)>2)],force=20) +
  theme(legend.text = element_text(size = rel(1.5)), legend.title = element_text(size = rel(1.5)))
pdf(file = "figures/DE_plot_ex_25min_trunc.pdf", width = 10, height = 7)
print(DE_plot_ex_25min_trunc)
dev.off()

png(file = "figures/DE_plot_ex_25min_fewer.png", width = 10, height = 7, units = "in", res = 600)
print(DE_plot_ex_25min)
dev.off()
png(file = "figures/DE_plot_ex_25min_fewer_trunc.png", width = 10, height = 7, units = "in", res = 600)
print(DE_plot_ex_25min_trunc)
dev.off()

