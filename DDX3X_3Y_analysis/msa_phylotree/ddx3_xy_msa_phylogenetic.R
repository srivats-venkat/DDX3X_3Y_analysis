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

setwd("/Users/srivatsv/Box Sync/PostDoc Stephen Floor Lab/Venkataramanan et al. DDX3X_3Y paper 2020/DDX3X_DDX3Y_manuscript_FloorLab_forGitHub/msa_phylotree/")
load_annotation("/Users/srivatsv/Box Sync/PostDoc Stephen Floor Lab/annotations/Homo.sapiens/gencode.v33.primary_assembly.annotation.gtf_Rannot")


DDX3_XY_humanonly_refseq_protein <- readAAStringSet(filepath = "./DDX3_XY_humanonly_refseq_protein.fasta", format = "fasta")
DDX3_XY_humanonly_refseq_protein_msa <- msa(DDX3_XY_humanonly_refseq_protein)
msaPrettyPrint(DDX3_XY_humanonly_refseq_protein_msa, 
               output="pdf", 
               showNames="left",
               showLogo = "none",
               askForOverwrite=FALSE, 
               verbose=FALSE,
               showConsensus = "none",
               paperWidth = 9)
                          

DDX3_XY_refseq_protein_subset_v3 <- readAAStringSet(filepath = "./DDX3_XY_refseq_protein_new_subset_v3.fasta", format = "fasta")
DDX3_XY_refseq_protein_msa_subset_v3 <- msa(DDX3_XY_refseq_protein_subset_v3)

DDX3_XY_refseq_protein_msa2_subset_v3 <- msaConvert(DDX3_XY_refseq_protein_msa_subset_v3, type="seqinr::alignment")
DDX3_XY_refseq_protein_msa2_subset_v3_dist <- dist.alignment(DDX3_XY_refseq_protein_msa2_subset_v3, "identity")
DDX3_XY_refseq_protein_subset_v3_tree <- nj(DDX3_XY_refseq_protein_msa2_subset_v3_dist)

DDX3_XY_refseq_protein_subset_v3_tree_plot <-  ggtree(DDX3_XY_refseq_protein_subset_v3_tree) + 
  geom_tiplab(size = 7) + 
  xlim(0,0.35)

pdf(file = "./DDX3_XY_refseq_protein_subset_v3_tree_plot.pdf", width = 18, height = 10)
print(DDX3_XY_refseq_protein_subset_v3_tree_plot)
dev.off()

DDX3_XY_refseq_protein_v4 <- readAAStringSet(filepath = "./DDX3_XY_v4.fasta", format = "fasta")
DDX3_XY_refseq_protein_v4_msa <- msa(DDX3_XY_refseq_protein_v4)

DDX3_XY_refseq_protein_v4_msa2 <- msaConvert(DDX3_XY_refseq_protein_v4_msa, type="seqinr::alignment")
DDX3_XY_refseq_protein_v4_msa2_dist <- dist.alignment(DDX3_XY_refseq_protein_v4_msa2, "identity")
DDX3_XY_refseq_protein_v4_tree <- nj(DDX3_XY_refseq_protein_v4_msa2_dist)

DDX3_XY_refseq_protein_v4_tree_plot <-  ggtree(DDX3_XY_refseq_protein_v4_tree) + 
  geom_tiplab(size = 1.5) + xlim(0,0.3)

pdf(file = "./DDX3_XY_refseq_protein_v4_tree_plot.pdf", width = 10, height = 10)
print(DDX3_XY_refseq_protein_v4_tree_plot)
dev.off()
