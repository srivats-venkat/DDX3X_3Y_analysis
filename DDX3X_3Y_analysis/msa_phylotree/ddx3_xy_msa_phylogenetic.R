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

setwd(".")
load_annotation("./gencode.v33.primary_assembly.annotation.gtf_Rannot")


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