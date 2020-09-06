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

#working directory requires the following sub-directories: ./counts (with .bam.count.summary files), ./riboseQC_annot (with .gtf.Rannot file) and ./figures
setwd(".")

#laods the annotation files
load_annotation("./annotations/Homo.sapiens/gencode.v33.primary_assembly.annotation.gtf_Rannot")
load("./annotations/Homo.sapiens/tx_regions_features.RData")

#loads and concatenates the counts (cds and exonic)
load("./counts/DMSO_HCT116_rep1_ribo_Aligned.sortedByCoord.out.bam_counts_summary_25min")
dmso_ribo_rep1_25min <- res_counts$counts_gene_cds_unq
load("./counts/DMSO_HCT116_rep2_ribo_Aligned.sortedByCoord.out.bam_counts_summary_25min")
dmso_ribo_rep2_25min <- res_counts$counts_gene_cds_unq
load("./counts/IAA_HCT116_rep1_ribo_Aligned.sortedByCoord.out.bam_counts_summary_25min")
iaa_ribo_rep1_25min <- res_counts$counts_gene_cds_unq
load("./counts/IAA_HCT116_rep2_ribo_Aligned.sortedByCoord.out.bam_counts_summary_25min")
iaa_ribo_rep2_25min <- res_counts$counts_gene_cds_unq
load("./counts/DMSO_HCT116_rep1_rna_Aligned.sortedByCoord.out.bam_counts_summary_25min")
dmso_rna_rep1_25min <- res_counts$counts_gene_cds_unq
load("./counts/DMSO_HCT116_rep2_rna_Aligned.sortedByCoord.out.bam_counts_summary_25min")
dmso_rna_rep2_25min <- res_counts$counts_gene_cds_unq
load("./counts/IAA_HCT116_rep1_rna_Aligned.sortedByCoord.out.bam_counts_summary_25min")
iaa_rna_rep1_25min <- res_counts$counts_gene_cds_unq
load("./counts/IAA_HCT116_rep2_rna_Aligned.sortedByCoord.out.bam_counts_summary_25min")
iaa_rna_rep2_25min <- res_counts$counts_gene_cds_unq
rm(res_counts)

gene_cds_unq_25min <- cbind(dmso_ribo_rep1_25min, 
                            dmso_ribo_rep2_25min, 
                            iaa_ribo_rep1_25min, 
                            iaa_ribo_rep2_25min, 
                            dmso_rna_rep1_25min, 
                            dmso_rna_rep2_25min, 
                            iaa_rna_rep1_25min, 
                            iaa_rna_rep2_25min)

colnames(gene_cds_unq_25min) <- c("dmso_ribo_rep1_25min", 
                                  "dmso_ribo_rep2_25min", 
                                  "iaa_ribo_rep1_25min", 
                                  "iaa_ribo_rep2_25min", 
                                  "dmso_rna_rep1_25min", 
                                  "dmso_rna_rep2_25min", 
                                  "iaa_rna_rep1_25min", 
                                  "iaa_rna_rep2_25min")

gene_cds_unq_25min <- as.data.frame(gene_cds_unq_25min)

rm(list = ls(pattern = "rep"))

#loads exp_design file and sets parameter values
exp_design <- read.delim("./sample_info.txt")
exp_design$Genotype <- factor(exp_design$Genotype, levels = c("dmso","iaa"))
exp_design$Condition <- factor(exp_design$Condition, levels = c( "rna","ribo"))

#DESeq for the riboRNA interaction
deseq_dataset_cds_25min <- DESeqDataSetFromMatrix(gene_cds_unq_25min, exp_design, design = ~ Condition + Genotype + Condition:Genotype)
deseq_dataset_cds_25min <- DESeq(deseq_dataset_cds_25min, test = "LRT", reduced = ~ Condition + Genotype)
deseq_results_cds_25min <- results(deseq_dataset_cds_25min, independentFiltering = FALSE)
deseq_results_cds_25min$pvalue[deseq_results_cds_25min$baseMean < 20] <- NA
deseq_results_cds_25min$padj <- p.adjust(deseq_results_cds_25min$pvalue, method = "BH")

#subsets the exp design and count files for rna only and ribo only DE
exp_design_rnaonly <- exp_design[c(5:8),]
exp_design_riboonly <- exp_design[c(1:4),]

gene_cds_unq_25min_rnaonly <- gene_cds_unq_25min[,c(5:8)]
gene_cds_unq_25min_riboonly <- gene_cds_unq_25min[,c(1:4)]

#DESeq for the individual rna and ribo terms
deseq_dataset_cds_25min_rnaonly <- DESeqDataSetFromMatrix(gene_cds_unq_25min_rnaonly, exp_design_rnaonly, design = ~ Genotype)
deseq_dataset_cds_25min_rnaonly <- DESeq(deseq_dataset_cds_25min_rnaonly)
deseq_results_cds_25min_rnaonly <- results(deseq_dataset_cds_25min_rnaonly, independentFiltering = FALSE)
deseq_results_cds_25min_rnaonly$pvalue[deseq_results_cds_25min_rnaonly$baseMean < 20] <- NA
deseq_results_cds_25min_rnaonly$padj <- p.adjust(deseq_results_cds_25min_rnaonly$pvalue, method = "BH")

deseq_dataset_cds_25min_riboonly <- DESeqDataSetFromMatrix(gene_cds_unq_25min_riboonly, exp_design_riboonly, design = ~ Genotype)
deseq_dataset_cds_25min_riboonly <- DESeq(deseq_dataset_cds_25min_riboonly)
deseq_results_cds_25min_riboonly <- results(deseq_dataset_cds_25min_riboonly, independentFiltering = FALSE)
deseq_results_cds_25min_riboonly$pvalue[deseq_results_cds_25min_riboonly$baseMean < 20] <- NA
deseq_results_cds_25min_riboonly$padj <- p.adjust(deseq_results_cds_25min_riboonly$pvalue, method = "BH")

#calls annotations from the GTF. Gene name, biotype 
deseq_results_cds_25min$symb<-GTF_annotation$trann$gene_name[match(rownames(deseq_results_cds_25min),GTF_annotation$trann$gene_id)]
deseq_results_cds_25min$biotype<-GTF_annotation$trann$gene_biotype[match(rownames(deseq_results_cds_25min),GTF_annotation$trann$gene_id)]
deseq_results_cds_25min_rnaonly$symb<-GTF_annotation$trann$gene_name[match(rownames(deseq_results_cds_25min_rnaonly),GTF_annotation$trann$gene_id)]
deseq_results_cds_25min_rnaonly$biotype<-GTF_annotation$trann$gene_biotype[match(rownames(deseq_results_cds_25min_rnaonly),GTF_annotation$trann$gene_id)]
deseq_results_cds_25min_riboonly$symb<-GTF_annotation$trann$gene_name[match(rownames(deseq_results_cds_25min_riboonly),GTF_annotation$trann$gene_id)]
deseq_results_cds_25min_riboonly$biotype<-GTF_annotation$trann$gene_biotype[match(rownames(deseq_results_cds_25min_riboonly),GTF_annotation$trann$gene_id)]

#builds a common results table with relevant parameters
results_cds_25min <- deseq_results_cds_25min[,c(1,2,6,7,8)]
colnames(results_cds_25min) <- c("baseMean_riborna", "log2FoldChange_riborna", "padj_riborna", "symb", "biotype")

results_cds_25min_rnaonly <- deseq_results_cds_25min_rnaonly[,c(1,2,6)]
colnames(results_cds_25min_rnaonly) <- c("baseMean_rna", "log2FoldChange_rna", "padj_rna")

results_cds_25min_riboonly <- deseq_results_cds_25min_riboonly[,c(1,2,6)]
colnames(results_cds_25min_riboonly) <- c("baseMean_ribo", "log2FoldChange_ribo", "padj_ribo")

results_cds_25min <- cbind(results_cds_25min, results_cds_25min_rnaonly, results_cds_25min_riboonly)
results_cds_25min <- results_cds_25min[order(results_cds_25min$padj_riborna),]

#imposes read number cutoffs
cutoff_cds<-which(results_cds_25min$baseMean_rna>20 & results_cds_25min$baseMean_ribo>20)

results_cds_25min <- results_cds_25min[cutoff_cds,]

# classifies genes based on the effect of perturbation
results_cds_25min$tx_type<-"Not_significant"
results_cds_25min$tx_type[results_cds_25min$padj_rna<.01 & results_cds_25min$log2FoldChange_rna<0]<-"RNA_down"
results_cds_25min$tx_type[results_cds_25min$padj_rna<.01 & results_cds_25min$log2FoldChange_rna>0]<-"RNA_up"
results_cds_25min$tx_type[results_cds_25min$padj_riborna<.05 & results_cds_25min$log2FoldChange_riborna<0]<-"TE_down"
results_cds_25min$tx_type[results_cds_25min$padj_riborna<.05 & results_cds_25min$log2FoldChange_riborna>0]<-"TE_up"

#log of p-value
results_cds_25min$log10_padj_riborna<-log10(results_cds_25min$padj_riborna)
results_cds_25min$log10_padj_riborna[is.na(results_cds_25min$log10_padj_riborna)]<-0

results_cds_25min$tx_type<-factor(results_cds_25min$tx_type,levels=c("TE_up","TE_down","RNA_up","RNA_down","Not_significant"))
results_cds_25min$biotype[grep(results_cds_25min$symb,pattern = "HIST|H2B")]<-"Histone genes"
results_cds_25min$biotype[grep(results_cds_25min$symb,pattern = "HIST|H2B",invert =  T)]<-"all genes"
results_cds_25min$biotype<-factor(results_cds_25min$biotype,levels=c("all genes","Histone genes"))

results_cds_25min$symb_deplot<-as.character(results_cds_25min$symb)
results_cds_25min$symb_deplot<-NA
toprna_cds<-results_cds_25min[results_cds_25min$tx_type=="RNA_up",]
toprna_cds<-toprna_cds[order(toprna_cds$log2FoldChange_rna,toprna_cds$log2FoldChange_ribo,decreasing = T),"symb"]
if(length(toprna_cds)>9){toprna_cds=toprna_cds[1:2]}
topribo_cds<-results_cds_25min[results_cds_25min$tx_type=="TE_up",]
topribo_cds<-topribo_cds[order(topribo_cds$log2FoldChange_riborna,decreasing = T),"symb"]
if(length(topribo_cds)>9){topribo_cds=topribo_cds[1:4]}
bottomrna_cds<-results_cds_25min[results_cds_25min$tx_type=="RNA_down",]
bottomrna_cds<-bottomrna_cds[order(bottomrna_cds$log2FoldChange_rna,bottomrna_cds$log2FoldChange_ribo,decreasing = F),"symb"]
if(length(bottomrna_cds)>9){bottomrna_cds=bottomrna_cds[1:2]}
bottomribo_cds<-results_cds_25min[results_cds_25min$tx_type=="TE_down",]
bottomribo_cds<-bottomribo_cds[order(bottomribo_cds$log2FoldChange_riborna,decreasing = F),"symb"]
if(length(bottomribo_cds)>9){bottomribo_cds=c(bottomribo_cds[1:3],bottomribo_cds[bottomribo_cds%in%c("ODC1","DVL1","DVL2","HMBS","PRKRA")])}
results_cds_25min$symb_deplot[results_cds_25min$symb%in%c(toprna_cds,topribo_cds,bottomrna_cds,bottomribo_cds)]<-results_cds_25min$symb[results_cds_25min$symb%in%c(toprna_cds,topribo_cds,bottomrna_cds,bottomribo_cds)]

results_cds_25min$nudgey<-NA
results_cds_25min$nudgey[which(nchar(results_cds_25min$symb_deplot)>2)]<-ifelse(results_cds_25min$log2FoldChange_ribo[which(nchar(results_cds_25min$symb_deplot)>2)]>0, .4, -.5)
results_cds_25min$nudgex<-NA
results_cds_25min$nudgex[which(nchar(results_cds_25min$symb_deplot)>2)]<-ifelse(results_cds_25min$log2FoldChange_rna[which(nchar(results_cds_25min$symb_deplot)>2)]>0, .3, -.4)
results_cds_25min$nudgex[results_cds_25min$tx_type=="TE_up"]<-0
results_cds_25min$nudgey[results_cds_25min$tx_type=="RNA_up"]<-0
results_cds_25min$nudgey[results_cds_25min$tx_type=="RNA_down"]<-0
results_cds_25min$nudgex[results_cds_25min$tx_type=="RNA_up"]<-0
results_cds_25min$nudgex[results_cds_25min$tx_type=="RNA_down"]<-0

results_cds_25min_df <- as.data.frame(results_cds_25min)
results_cds_25min_df_iaadmso <- results_cds_25min_df
write.csv(results_cds_25min_df, "data/results_cds_25min_df.csv")

DE_plot_cds_25min <- ggplot(results_cds_25min_df,aes(x=log2FoldChange_rna,y=log2FoldChange_ribo,color=tx_type,size=-log10_padj_riborna,label=symb_deplot)) + 
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
  geom_text_repel(size=5,nudge_x = results_cds_25min_df$nudgex[which(nchar(results_cds_25min_df$symb_deplot)>2)],
                  nudge_y = results_cds_25min_df$nudgey[which(nchar(results_cds_25min_df$symb_deplot)>2)],force=20) +
  theme(legend.text = element_text(size = rel(1.5)), legend.title = element_text(size = rel(1.5)))
pdf(file = "figures/DE_plot_cds_25min.pdf", width = 10, height = 7)
print(DE_plot_cds_25min)
dev.off()
png(file = "figures/DE_plot_cds_25min.png", width = 10, height = 7, units = "in", res = 600)
print(DE_plot_cds_25min)
dev.off()

split_results_cds_25min <- results_cds_25min_df
split_results_cds_25min$gene_id <- rownames(split_results_cds_25min)
split_results_cds_25min <- split(split_results_cds_25min$gene_id, f = split_results_cds_25min$tx_type)
for(i in names(split_results_cds_25min)[1:5]){
  tx_regions_feats$tx_type[tx_regions_feats$gene_id%in%split_results_cds_25min[[i]]]<-i
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

TE_down_utrratio <- results_cds_25min_df[which(results_cds_25min_df$tx_type == "TE_down"),]
TE_up_utrratio <- results_cds_25min_df[which(results_cds_25min_df$tx_type == "TE_up"),]
RNA_down_utrratio <- results_cds_25min_df[which(results_cds_25min_df$tx_type == "RNA_down"),]
RNA_up_utrratio <- results_cds_25min_df[which(results_cds_25min_df$tx_type == "RNA_up"),]
not_significant_utrratio <- results_cds_25min_df[which(results_cds_25min_df$tx_type == "Not_significant"),]

TE_down_utrratio$five_cds_iaa_dmso<-utr_counts_ratio$iaa_dmso[match(rownames(TE_down_utrratio),rownames(utr_counts_ratio))]
TE_up_utrratio$five_cds_iaa_dmso<-utr_counts_ratio$iaa_dmso[match(rownames(TE_up_utrratio),rownames(utr_counts_ratio))]
RNA_down_utrratio$five_cds_iaa_dmso<-utr_counts_ratio$iaa_dmso[match(rownames(RNA_down_utrratio),rownames(utr_counts_ratio))]
RNA_up_utrratio$five_cds_iaa_dmso<-utr_counts_ratio$iaa_dmso[match(rownames(RNA_up_utrratio),rownames(utr_counts_ratio))]
not_significant_utrratio$five_cds_iaa_dmso<-utr_counts_ratio$iaa_dmso[match(rownames(not_significant_utrratio),rownames(utr_counts_ratio))]

all_utrratio <- results_cds_25min_df 
all_utrratio$five_cds_iaa_dmso<-utr_counts_ratio$iaa_dmso[match(rownames(all_utrratio),rownames(utr_counts_ratio))]

load("./counts/IAA3XWT_HCT116_rep1_ribo_Aligned.sortedByCoord.out.bam_counts_summary_25min")
utr_counts <- res_counts$counts_5cds3_unq
utr_counts <- as.data.frame(utr_counts)
colnames(utr_counts) <- c("fiveutr", "cds" , "threeutr")
utr_counts$fiveutrbycds <- utr_counts$fiveutr/utr_counts$cds
utr_counts_ratio <- as.data.frame(utr_counts[,4])
rownames(utr_counts_ratio) <- rownames(utr_counts)
load("./counts/IAA3XWT_HCT116_rep2_ribo_Aligned.sortedByCoord.out.bam_counts_summary_25min")
utr_counts <- res_counts$counts_5cds3_unq
utr_counts <- as.data.frame(utr_counts)
colnames(utr_counts) <- c("fiveutr", "cds" , "threeutr")
utr_counts$fiveutrbycds <- utr_counts$fiveutr/utr_counts$cds
utr_counts_ratio <- cbind(utr_counts_ratio, utr_counts$fiveutrbycds)
load("./counts/IAA3Y_HCT116_rep1_ribo_Aligned.sortedByCoord.out.bam_counts_summary_25min")
utr_counts <- res_counts$counts_5cds3_unq
utr_counts <- as.data.frame(utr_counts)
colnames(utr_counts) <- c("fiveutr", "cds" , "threeutr")
utr_counts$fiveutrbycds <- utr_counts$fiveutr/utr_counts$cds
utr_counts_ratio <- cbind(utr_counts_ratio, utr_counts$fiveutrbycds)
load("./counts/IAA3Y_HCT116_rep2_ribo_Aligned.sortedByCoord.out.bam_counts_summary_25min")
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
ggsave(filename = "5riboskew_IAAvsDMSO_3Yvs3X.pdf", plot = riboskew_compare, device = "pdf", units = "in", width = 14, height = 10.8, path = "./figures")

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

deltaTE_vs_baseTE_degron_df <- results_cds_25min_df
gene_cds_unq_25min$TE_base <- log2((gene_cds_unq_25min$dmso_ribo_rep1_25min + gene_cds_unq_25min$dmso_ribo_rep2_25min) / (gene_cds_unq_25min$dmso_rna_rep1_25min + gene_cds_unq_25min$dmso_rna_rep2_25min))
deltaTE_vs_baseTE_degron_df$TE_base <- gene_cds_unq_25min$TE_base[match(rownames(deltaTE_vs_baseTE_degron_df), rownames(gene_cds_unq_25min))]

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




deltaTE_vs_baseTE_degron_df$Tx_class<-factor(deltaTE_vs_baseTE_degron_df$tx_type,levels=rev(c("Not_significant","RNA_up","RNA_down","TE_up","TE_down")))

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

############################################################################################################################################


#working directory requires the following sub-directories: ./counts (with .bam.count.summary files), ./riboseQC_annot (with .gtf.Rannot file) and ./figures
setwd(".")

#loads and concatenates the counts (cds and exonic)
load("./counts/IAA3XWT_HCT116_rep1_ribo_Aligned.sortedByCoord.out.bam_counts_summary_25min")
ddx3x_ribo_rep1_25min <- res_counts$counts_gene_cds_unq
load("./counts/IAA3XWT_HCT116_rep2_ribo_Aligned.sortedByCoord.out.bam_counts_summary_25min")
ddx3x_ribo_rep2_25min <- res_counts$counts_gene_cds_unq
load("./counts/IAA3Y_HCT116_rep1_ribo_Aligned.sortedByCoord.out.bam_counts_summary_25min")
ddx3y_ribo_rep1_25min <- res_counts$counts_gene_cds_unq
load("./counts/IAA3Y_HCT116_rep2_ribo_Aligned.sortedByCoord.out.bam_counts_summary_25min")
ddx3y_ribo_rep2_25min <- res_counts$counts_gene_cds_unq
load("./counts/IAA3XWT_HCT116_rep1_rna_Aligned.sortedByCoord.out.bam_counts_summary_25min")
ddx3x_rna_rep1_25min <- res_counts$counts_gene_cds_unq
load("./counts/IAA3XWT_HCT116_rep2_rna_Aligned.sortedByCoord.out.bam_counts_summary_25min")
ddx3x_rna_rep2_25min <- res_counts$counts_gene_cds_unq
load("./counts/IAA3Y_HCT116_rep1_rna_Aligned.sortedByCoord.out.bam_counts_summary_25min")
ddx3y_rna_rep1_25min <- res_counts$counts_gene_cds_unq
load("./counts/IAA3Y_HCT116_rep2_rna_Aligned.sortedByCoord.out.bam_counts_summary_25min")
ddx3y_rna_rep2_25min <- res_counts$counts_gene_cds_unq
rm(res_counts)

gene_cds_unq_25min <- cbind(ddx3x_ribo_rep1_25min, 
                            ddx3x_ribo_rep2_25min, 
                            ddx3y_ribo_rep1_25min, 
                            ddx3y_ribo_rep2_25min, 
                            ddx3x_rna_rep1_25min, 
                            ddx3x_rna_rep2_25min, 
                            ddx3y_rna_rep1_25min, 
                            ddx3y_rna_rep2_25min)
colnames(gene_cds_unq) <- c("ddx3x_ribo_rep1_25min", 
                            "ddx3x_ribo_rep2_25min", 
                            "ddx3y_ribo_rep1_25min", 
                            "ddx3y_ribo_rep2_25min", 
                            "ddx3x_rna_rep1_25min", 
                            "ddx3x_rna_rep2_25min", 
                            "ddx3y_rna_rep1_25min", 
                            "ddx3y_rna_rep2_25min")
gene_cds_unq_25min <- as.data.frame(gene_cds_unq_25min)


rm(list = ls(pattern = "rep"))

#loads exp_design file and sets parameter values
exp_design <- read.delim("./sample_info.txt")
exp_design$Genotype <- factor(exp_design$Genotype, levels = c("ddx3x","ddx3y"))
exp_design$Condition <- factor(exp_design$Condition, levels = c( "rna","ribo"))

#DESeq for the riboRNA interaction
deseq_dataset_cds_25min <- DESeqDataSetFromMatrix(gene_cds_unq_25min, exp_design, design = ~ Condition + Genotype + Condition:Genotype)
deseq_dataset_cds_25min <- DESeq(deseq_dataset_cds_25min, test = "LRT", reduced = ~ Condition + Genotype)
deseq_results_cds_25min <- results(deseq_dataset_cds_25min, independentFiltering = FALSE)
deseq_results_cds_25min$pvalue[deseq_results_cds_25min$baseMean < 20] <- NA
deseq_results_cds_25min$padj <- p.adjust(deseq_results_cds_25min$pvalue, method = "BH")

#subsets the exp design and count files for rna only and ribo only DE
exp_design_rnaonly <- exp_design[c(5:8),]
exp_design_riboonly <- exp_design[c(1:4),]
gene_cds_unq_25min_rnaonly <- gene_cds_unq_25min[,c(5:8)]
gene_cds_unq_25min_riboonly <- gene_cds_unq_25min[,c(1:4)]

#DESeq for the individual rna and ribo terms
deseq_dataset_cds_25min_rnaonly <- DESeqDataSetFromMatrix(gene_cds_unq_25min_rnaonly, exp_design_rnaonly, design = ~ Genotype)
deseq_dataset_cds_25min_rnaonly <- DESeq(deseq_dataset_cds_25min_rnaonly)
deseq_results_cds_25min_rnaonly <- results(deseq_dataset_cds_25min_rnaonly, independentFiltering = FALSE)
deseq_results_cds_25min_rnaonly$pvalue[deseq_results_cds_25min_rnaonly$baseMean < 20] <- NA
deseq_results_cds_25min_rnaonly$padj <- p.adjust(deseq_results_cds_25min_rnaonly$pvalue, method = "BH")

deseq_dataset_cds_25min_riboonly <- DESeqDataSetFromMatrix(gene_cds_unq_25min_riboonly, exp_design_riboonly, design = ~ Genotype)
deseq_dataset_cds_25min_riboonly <- DESeq(deseq_dataset_cds_25min_riboonly)
deseq_results_cds_25min_riboonly <- results(deseq_dataset_cds_25min_riboonly, independentFiltering = FALSE)
deseq_results_cds_25min_riboonly$pvalue[deseq_results_cds_25min_riboonly$baseMean < 20] <- NA
deseq_results_cds_25min_riboonly$padj <- p.adjust(deseq_results_cds_25min_riboonly$pvalue, method = "BH")

#calls annotations from the GTF. Gene name, biotype 
deseq_results_cds_25min$symb<-GTF_annotation$trann$gene_name[match(rownames(deseq_results_cds_25min),GTF_annotation$trann$gene_id)]
deseq_results_cds_25min$biotype<-GTF_annotation$trann$gene_biotype[match(rownames(deseq_results_cds_25min),GTF_annotation$trann$gene_id)]
deseq_results_cds_25min_rnaonly$symb<-GTF_annotation$trann$gene_name[match(rownames(deseq_results_cds_25min_rnaonly),GTF_annotation$trann$gene_id)]
deseq_results_cds_25min_rnaonly$biotype<-GTF_annotation$trann$gene_biotype[match(rownames(deseq_results_cds_25min_rnaonly),GTF_annotation$trann$gene_id)]
deseq_results_cds_25min_riboonly$symb<-GTF_annotation$trann$gene_name[match(rownames(deseq_results_cds_25min_riboonly),GTF_annotation$trann$gene_id)]
deseq_results_cds_25min_riboonly$biotype<-GTF_annotation$trann$gene_biotype[match(rownames(deseq_results_cds_25min_riboonly),GTF_annotation$trann$gene_id)]

#builds a common results table with relevant parameters
results_cds_25min <- deseq_results_cds_25min[,c(1,2,6,7,8)]
colnames(results_cds_25min) <- c("baseMean_riborna", "log2FoldChange_riborna", "padj_riborna", "symb", "biotype")

results_cds_25min_rnaonly <- deseq_results_cds_25min_rnaonly[,c(1,2,6)]
colnames(results_cds_25min_rnaonly) <- c("baseMean_rna", "log2FoldChange_rna", "padj_rna")

results_cds_25min_riboonly <- deseq_results_cds_25min_riboonly[,c(1,2,6)]
colnames(results_cds_25min_riboonly) <- c("baseMean_ribo", "log2FoldChange_ribo", "padj_ribo")

results_cds_25min <- cbind(results_cds_25min, results_cds_25min_rnaonly, results_cds_25min_riboonly)
results_cds_25min <- results_cds_25min[order(results_cds_25min$padj_riborna),]

#imposes read number cutoffs
cutoff_cds<-which(results_cds_25min$baseMean_rna>20 & results_cds_25min$baseMean_ribo>20)

results_cds_25min <- results_cds_25min[cutoff_cds,]

# classifies genes based on the effect of perturbation
results_cds_25min$tx_type<-"Not_significant"
results_cds_25min$tx_type[results_cds_25min$padj_rna<.01 & results_cds_25min$log2FoldChange_rna<0]<-"RNA_down"
results_cds_25min$tx_type[results_cds_25min$padj_rna<.01 & results_cds_25min$log2FoldChange_rna>0]<-"RNA_up"
results_cds_25min$tx_type[results_cds_25min$padj_riborna<.05 & results_cds_25min$log2FoldChange_riborna<0]<-"TE_down"
results_cds_25min$tx_type[results_cds_25min$padj_riborna<.05 & results_cds_25min$log2FoldChange_riborna>0]<-"TE_up"

#log of p-value
results_cds_25min$log10_padj_riborna<-log10(results_cds_25min$padj_riborna)
results_cds_25min$log10_padj_riborna[is.na(results_cds_25min$log10_padj_riborna)]<-0

results_cds_25min$tx_type<-factor(results_cds_25min$tx_type,levels=c("TE_up","TE_down","RNA_up","RNA_down","Not_significant"))
results_cds_25min$biotype[grep(results_cds_25min$symb,pattern = "HIST|H2B")]<-"Histone genes"
results_cds_25min$biotype[grep(results_cds_25min$symb,pattern = "HIST|H2B",invert =  T)]<-"all genes"
results_cds_25min$biotype<-factor(results_cds_25min$biotype,levels=c("all genes","Histone genes"))

results_cds_25min$symb_deplot<-as.character(results_cds_25min$symb)
results_cds_25min$symb_deplot<-NA
toprna_cds<-results_cds_25min[results_cds_25min$tx_type=="RNA_up",]
toprna_cds<-toprna_cds[order(toprna_cds$log2FoldChange_rna,toprna_cds$log2FoldChange_ribo,decreasing = T),"symb"]
if(length(toprna_cds)>6){toprna_cds=toprna_cds[1:5]}
topribo_cds<-results_cds_25min[results_cds_25min$tx_type=="TE_up",]
topribo_cds<-topribo_cds[order(topribo_cds$log2FoldChange_riborna,decreasing = T),"symb"]
if(length(topribo_cds)>1){topribo_cds=topribo_cds[1:1]}
bottomrna_cds<-results_cds_25min[results_cds_25min$tx_type=="RNA_down",]
bottomrna_cds<-bottomrna_cds[order(bottomrna_cds$log2FoldChange_rna,bottomrna_cds$log2FoldChange_ribo,decreasing = F),"symb"]
if(length(bottomrna_cds)>2){bottomrna_cds=bottomrna_cds[1:2]}
bottomribo_cds<-results_cds_25min[results_cds_25min$tx_type=="TE_down",]
bottomribo_cds<-bottomribo_cds[order(bottomribo_cds$log2FoldChange_riborna,decreasing = F),"symb"]
if(length(bottomribo_cds)>1){bottomribo_cds=c(bottomribo_cds[1:1],bottomribo_cds[bottomribo_cds%in%c("ODC1","DVL1","DVL2","HMBS","PRKRA")])}
results_cds_25min$symb_deplot[results_cds_25min$symb%in%c(toprna_cds,topribo_cds,bottomrna_cds,bottomribo_cds)]<-results_cds_25min$symb[results_cds_25min$symb%in%c(toprna_cds,topribo_cds,bottomrna_cds,bottomribo_cds)]
results_cds_25min$symb_deplot[which(results_cds_25min$symb=="ODC1")] <- "ODC1"

results_cds_25min$nudgey<-NA
results_cds_25min$nudgey[which(nchar(results_cds_25min$symb_deplot)>2)]<-ifelse(results_cds_25min$log2FoldChange_ribo[which(nchar(results_cds_25min$symb_deplot)>2)]>0, .4, -.5)
results_cds_25min$nudgex<-NA
results_cds_25min$nudgex[which(nchar(results_cds_25min$symb_deplot)>2)]<-ifelse(results_cds_25min$log2FoldChange_rna[which(nchar(results_cds_25min$symb_deplot)>2)]>0, .3, -.4)
results_cds_25min$nudgex[results_cds_25min$tx_type=="TE_up"]<-0
results_cds_25min$nudgey[results_cds_25min$tx_type=="RNA_up"]<-0
results_cds_25min$nudgey[results_cds_25min$tx_type=="RNA_down"]<-0
results_cds_25min$nudgex[results_cds_25min$tx_type=="RNA_up"]<-0
results_cds_25min$nudgex[results_cds_25min$tx_type=="RNA_down"]<-0

results_cds_25min_df <- as.data.frame(results_cds_25min)

#flattens log10(padj) to a maximum value of 15. Aesthetic purposes only, the huge p-value for DDX3Y was skewing the point sizes in the plot
results_cds_25min_df$log10_padj_riborna_flat <- results_cds_25min_df$log10_padj_riborna
results_cds_25min_df$log10_padj_riborna_flat[results_cds_25min_df$log10_padj_riborna < -15]<- -15

results_cds_25min_df_3X3Y <-results_cds_25min_df
results_cds_25min_df_3X3Y$tx_class_IAAvsDMSO <- results_cds_25min_df_iaadmso$tx_type[match(rownames(results_cds_25min_df_3X3Y), rownames(results_cds_25min_df_iaadmso))]
write.csv(results_cds_25min_df_3X3Y, "data/DE_table_3X3Y_cds_3X3Y.csv")


DE_plot_cds_25min <- ggplot(results_cds_25min_df,aes(x=log2FoldChange_rna,y=log2FoldChange_ribo,color=tx_type,size=-log10_padj_riborna_flat,label=symb_deplot)) + 
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
  geom_text_repel(size=5,nudge_x = results_cds_25min_df$nudgex[which(nchar(results_cds_25min_df$symb_deplot)>2)],
                  nudge_y = results_cds_25min_df$nudgey[which(nchar(results_cds_25min_df$symb_deplot)>2)],force=20)
pdf(file = "figures/DE_plot_cds_25min.pdf", width = 10, height = 7) +
  theme(legend.text = element_text(size = rel(1.5)), legend.title = element_text(size = rel(1.5)))
print(DE_plot_cds_25min)
dev.off()

DE_plot_cds_25min_trunc <- ggplot(results_cds_25min_df,aes(x=log2FoldChange_rna,y=log2FoldChange_ribo,color=tx_type,size=-log10_padj_riborna_flat,label=symb_deplot)) + 
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
  geom_text_repel(size=5,nudge_x = results_cds_25min_df$nudgex[which(nchar(results_cds_25min_df$symb_deplot)>2)],
                  nudge_y = results_cds_25min_df$nudgey[which(nchar(results_cds_25min_df$symb_deplot)>2)],force=50) +
  theme(legend.text = element_text(size = rel(1.5)), legend.title = element_text(size = rel(1.5)))
pdf(file = "figures/DE_plot_cds_25min_trunc.pdf", width = 10, height = 7)
print(DE_plot_cds_25min_trunc)
dev.off()

png(file = "figures/DE_plot_cds_25min_fewer.png", width = 10, height = 7, units = "in", res = 600)
print(DE_plot_cds_25min)
dev.off()

png(file = "figures/DE_plot_cds_25min_trunc_fewer.png", width = 10, height = 7, units = "in", res = 600)
print(DE_plot_cds_25min_trunc)
dev.off()

results_cds_25min_df_iaadmso$deltaTE_IAA3Y_IAA3XWT<-results_cds_25min_df$log2FoldChange_riborna[match(rownames(results_cds_25min_df_iaadmso),rownames(results_cds_25min_df))]
results_cds_25min_df_iaadmso$tx_type<-factor(results_cds_25min_df_iaadmso$tx_type,levels=rev(c("Not_significant",
                                                                                               "RNA_up",
                                                                                               "RNA_down",
                                                                                               "TE_up",
                                                                                               "TE_down")))
deltaTE3Yvs3X<- ggplot(results_cds_25min_df_iaadmso,aes(deltaTE_IAA3Y_IAA3XWT,group=tx_type,y=tx_type,fill=tx_type))  +
  theme_ridges() +
  scale_fill_manual(values=c("blue","dark red","cornflowerblue","firebrick1","dark grey")) +
  xlab("log2FC TE DDX3Y/DDX3X")  +
  ylab("tx_type as in IAA/DMSO") +
  stat_density_ridges(quantile_lines = TRUE, quantiles = 2) +
  theme(axis.title.x = element_text(size=22),axis.text.x  = element_text( vjust=0.5, size=18)) +
  theme(axis.title.y = element_text(size=24),axis.text.y  = element_text( vjust=0.5, size=18))  +
  scale_x_continuous(limits = c(-2,2)) +
  theme(legend.position="none")

png(file = "figures/deltaTE3Yvs3X.png", width = 10, height = 7, units = "in", res = 600)
print(deltaTE3Yvs3X)
dev.off()

pdf(file = "figures/deltaTE3Yvs3X.pdf", width = 5, height = 5)
print(deltaTE3Yvs3X)
dev.off()

temp <- results_cds_25min_df
temp <- temp[order(temp$log2FoldChange_riborna),]
temp$tx_type[c(1:100)]<- "TE_down"
temp <- temp[order(temp$log2FoldChange_riborna, decreasing = T),]
temp$tx_type[c(1:100)]<- "TE_up"

load("./annotations/Homo.sapiens/tx_regions_features.RData")
split_results_cds_25min <- temp
split_results_cds_25min$gene_id <- rownames(split_results_cds_25min)
split_results_cds_25min <- split(split_results_cds_25min$gene_id, f = split_results_cds_25min$tx_type)
for(i in names(split_results_cds_25min)[1:5]){
  tx_regions_feats$tx_type[tx_regions_feats$gene_id%in%split_results_cds_25min[[i]]]<-i
}

features_df <- as.data.frame(mcols(tx_regions_feats))
features_df$tx_type<-factor(features_df$tx_type,levels=rev(c("TE_up",
                                                             "TE_down")))
features_df_completecases <- features_df[complete.cases(features_df),]

fiveUTR_GC_top100XvsY <- ggplot(features_df_completecases,aes(G.C,group=tx_type,y=tx_type,fill=tx_type))  +
  theme_ridges() +
  scale_fill_manual(values=c("blue","dark red","cornflowerblue","firebrick1","dark grey")) +
  xlab("%GC content in 5'UTR")  +
  ylab("Top 100 genes changes in TE \n DDX3Y vs. DDX3X") +
  stat_density_ridges(quantile_lines = TRUE, quantiles = 2) +
  theme(legend.position="none") +
  theme(axis.title.x = element_text(size=22),axis.text.x  = element_text( vjust=0.5, size=18)) +
  theme(axis.title.y = element_text(size=24),axis.text.y  = element_text( vjust=0.5, size=18))

temp <- results_cds_25min_df_iaadmso
temp <- temp[order(temp$log2FoldChange_riborna),]
temp$tx_type[c(1:100)]<- "TE_down"
temp <- temp[order(temp$log2FoldChange_riborna, decreasing = T),]
temp$tx_type[c(1:100)]<- "TE_up"

load("./annotations/Homo.sapiens/tx_regions_features.RData")
split_results_cds_25min <- temp
split_results_cds_25min$gene_id <- rownames(split_results_cds_25min)
split_results_cds_25min <- split(split_results_cds_25min$gene_id, f = split_results_cds_25min$tx_type)
for(i in names(split_results_cds_25min)[1:5]){
  tx_regions_feats$tx_type[tx_regions_feats$gene_id%in%split_results_cds_25min[[i]]]<-i
}

features_df <- as.data.frame(mcols(tx_regions_feats))
features_df$tx_type<-factor(features_df$tx_type,levels=rev(c("TE_up",
                                                             "TE_down")))
features_df_completecases <- features_df[complete.cases(features_df),]

fiveUTR_GC_top100IAAvsDMSO <- ggplot(features_df_completecases,aes(G.C,group=tx_type,y=tx_type,fill=tx_type))  +
  theme_ridges() +
  scale_fill_manual(values=c("blue","dark red","cornflowerblue","firebrick1","dark grey")) +
  xlab("%GC content in 5'UTR")  +
  ylab("Top 100 genes changes in TE \n IAA vs. DMSO") +
  stat_density_ridges(quantile_lines = TRUE, quantiles = 2) +
  theme(legend.position="none") +
  theme(axis.title.x = element_text(size=22),axis.text.x  = element_text( vjust=0.5, size=18)) +
  theme(axis.title.y = element_text(size=24),axis.text.y  = element_text( vjust=0.5, size=18))

fiveUTR_GC_top100 <- fiveUTR_GC_top100XvsY | fiveUTR_GC_top100IAAvsDMSO

png(file = "figures/fiveUTR_GC_top100.png", width = 20, height = 7, units = "in", res = 600)
print(fiveUTR_GC_top100)
dev.off()

pdf(file = "figures/fiveUTR_GC_top100.pdf", width = 10, height = 5)
print(fiveUTR_GC_top100)
dev.off()


effsize <- as.data.frame(c(100,200,300,400,500,600,700,800,900,1000))
colnames(effsize) <- "N_top"
effsize$YvsX <- c(0,0,0,0,0,0,0,0,0,0)
load("./annotations/Homo.sapiens/tx_regions_features.RData")
k=0
for (j in c(100,200,300,400,500,600,700,800,900,1000)) {
  k=k+1
  temp <- results_cds_25min_df
  temp$tx_type <- "Not_significant"
  temp <- temp[order(temp$log2FoldChange_riborna),]
  temp$tx_type[c(1:j)]<- "TE_down"
  temp <- temp[order(temp$log2FoldChange_riborna, decreasing = T),]
  temp$tx_type[c(1:j)]<- "TE_up"
  split_results_cds_25min <- temp
  split_results_cds_25min$gene_id <- rownames(split_results_cds_25min)
  split_results_cds_25min <- split(split_results_cds_25min$gene_id, f = split_results_cds_25min$tx_type)
  for(i in names(split_results_cds_25min)[1:5]){
    tx_regions_feats$tx_type[tx_regions_feats$gene_id%in%split_results_cds_25min[[i]]]<-i
  }
  features_df <- as.data.frame(mcols(tx_regions_feats))
  temp <-features_df$G.C[which(features_df$tx_type=="TE_down")]
  temp2 <-features_df$G.C[which(features_df$tx_type=="TE_up")]
  effsize$YvsX[k] <- cliff.delta(temp, temp2, return.dm = T)$estimate
}

effsize$IAAvsDMSO <- c(0,0,0,0,0,0,0,0,0,0)
load("./annotations/Homo.sapiens/tx_regions_features.RData")
k=0
for (j in c(100,200,300,400,500,600,700,800,900,1000)) {
  k=k+1
  temp <- results_cds_25min_df_iaadmso
  temp$tx_type <- "Not_significant"
  temp <- temp[order(temp$log2FoldChange_riborna),]
  temp$tx_type[c(1:j)]<- "TE_down"
  temp <- temp[order(temp$log2FoldChange_riborna, decreasing = T),]
  temp$tx_type[c(1:j)]<- "TE_up"
  split_results_cds_25min <- temp
  split_results_cds_25min$gene_id <- rownames(split_results_cds_25min)
  split_results_cds_25min <- split(split_results_cds_25min$gene_id, f = split_results_cds_25min$tx_type)
  for(i in names(split_results_cds_25min)[1:5]){
    tx_regions_feats$tx_type[tx_regions_feats$gene_id%in%split_results_cds_25min[[i]]]<-i
  }
  features_df <- as.data.frame(mcols(tx_regions_feats))
  temp <-features_df$G.C[which(features_df$tx_type=="TE_down")]
  temp2 <-features_df$G.C[which(features_df$tx_type=="TE_up")]
  effsize$IAAvsDMSO[k] <- cliff.delta(temp, temp2, return.dm = T)$estimate
}

effsize_melt <- melt(effsize, id.vars = c("N_top"))
colnames(effsize_melt) <- c("N_top","comparison","clifsD_estimate")

fiveUTR_GC_top_effsize <- ggplot(effsize_melt, aes(x=N_top, y=clifsD_estimate, color = comparison)) +
  geom_point(size =4) +
  geom_line(linetype = "dashed", alpha =0.5) +
  theme_pubclean() +
  geom_hline(yintercept = 0, alpha = 0.5) +
  scale_color_manual(values=c("black","dark grey"))+
  scale_y_continuous(limits = c(-0.5,0.5), breaks = c(-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5)) +
  scale_x_continuous(limits = c(0,1100), breaks = c(100,200,300,400,500,600,700,800,900,1000)) +
  ggtitle(label = "Effect size - top and bottom ∆TE\n 5' UTR GC content") +
  ylab("Cliffs Delta estimate") +
  xlab("Number of genes considered") +
  theme(axis.title.x = element_text(size=22),axis.text.x  = element_text(angle=45, vjust=0.5, size=18)) +
  theme(axis.title.y = element_text(size=22),axis.text.y  = element_text(angle=0, vjust=0.5, size=18)) +
  theme(plot.title = element_text(size = 30)) +
  theme(legend.text = element_text(size = rel(1.5)), legend.title = element_text(size = rel(1.5)))

fiveUTR_GC_top_effsize

png(file = "figures/fiveUTR_GC_top_effsize.png", width = 10, height = 7, units = "in", res = 600)
print(fiveUTR_GC_top_effsize)
dev.off()

pdf(file = "figures/fiveUTR_GC_top_effsize.pdf", width = 10, height = 7)
print(fiveUTR_GC_top_effsize)
dev.off()


