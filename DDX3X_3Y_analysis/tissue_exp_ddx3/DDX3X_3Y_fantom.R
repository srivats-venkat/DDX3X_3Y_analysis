library(ggplot2)
library(patchwork)
library(dplyr)
library(wordspace)
library(matrixStats)
library(scales)
library(reshape2)

setwd(".")
rna_tissue_fantom <- as.data.frame(read.table(file = 'rna_tissue_fantom.tsv', sep = '\t', header = TRUE))


rna_tissue_fantom$sex_spec <- " Other Common Tissues"
rna_tissue_fantom$sex_spec[which(rna_tissue_fantom$Tissue == "testis" | rna_tissue_fantom$Tissue == "prostate" | rna_tissue_fantom$Tissue == "seminal vesicle" | rna_tissue_fantom$Tissue == "epididymis" | rna_tissue_fantom$Tissue =="ductus deferens")] <- "Male Tissues"
rna_tissue_fantom$sex_spec[which(rna_tissue_fantom$Tissue == "breast" | rna_tissue_fantom$Tissue == "ovary" | rna_tissue_fantom$Tissue == "endometrium" | rna_tissue_fantom$Tissue == "cervix, uterine" | rna_tissue_fantom$Tissue == "fallopian tube" | rna_tissue_fantom$Tissue == "vagina")] <- "Female Tissues"
rna_tissue_fantom$sex_spec[which(rna_tissue_fantom$Tissue == "amygdala" | rna_tissue_fantom$Tissue == "basal ganglia" | rna_tissue_fantom$Tissue == "cerebellum" | rna_tissue_fantom$Tissue == "cerebral cortex" | rna_tissue_fantom$Tissue == "corpus callosum" | rna_tissue_fantom$Tissue == "corpus callosum" | rna_tissue_fantom$Tissue == "hippocampal formation" | rna_tissue_fantom$Tissue == "midbrain" | rna_tissue_fantom$Tissue == "olfactory region" | rna_tissue_fantom$Tissue == "pons and medulla" | rna_tissue_fantom$Tissue == "thalamus")] <- "Brain Tissues"

subset_fantom <- rna_tissue_fantom[which(rna_tissue_fantom$Gene.name == "DDX3X" | rna_tissue_fantom$Gene.name == "DDX3Y"),]
subset_fantom$Gene.name <- factor(subset_fantom$Gene.name)

DDX3X_3Y_fantom <- ggplot(data = subset_fantom, mapping = aes(x = Tissue, 
                                                              y = Gene.name, 
                                                              fill = Scaled.tags.per.million)) +
  geom_tile() +
  xlab(label = element_blank()) +
  ylab(label = element_blank()) +
  facet_grid(~ subset_fantom$sex_spec, switch = "x", scales = "free_x", space = "free_x") +
  scale_fill_gradient(name = "Abundance",
                      low = "#FFFFFF",
                      high = "#31a354") +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = "#EEEEEE", color = "#FFFFFF")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, face = "bold", size = 12),
        axis.title.x = element_text(face = "bold", size = 12),
        strip.text.x = element_text(face = "bold", size = 12)) +
  theme(axis.text.y = element_text(angle = 90, hjust = 0.5, face = "bold", size = 12),
        axis.title.y = element_text(face = "bold", size = 12)) +
  theme(legend.text = element_text(face = "bold", size = 8),
        legend.title = element_text(face = "bold", size = 12))

DDX3X_3Y_fantom

DDX3X_3Y_fantom_vert <- ggplot(data = subset_fantom, mapping = aes(x = Gene.name, 
                                                                   y = Tissue, 
                                                                   fill = Scaled.tags.per.million)) +
  geom_tile() +
  xlab(label = element_blank()) +
  ylab(label = element_blank()) +
  facet_grid(subset_fantom$sex_spec ~ . , scales = "free_y", switch = "y", space = "free_y") +
  scale_fill_gradient(name = "Abundance",
                      low = "#FFFFFF",
                      high = "#31a354") +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = "#EEEEEE", color = "#FFFFFF")) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, face = "bold", size = 14),
        axis.title.x = element_text(face = "bold", size = 14),
        strip.text.x = element_text(face = "bold", size = 12)) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, face = "bold", size = 14),
        axis.title.y = element_text(face = "bold", size = 12),
        strip.text.y  = element_text(face = "bold", size = 9)) +
  theme(legend.text = element_text(face = "bold", size = 12),
        legend.title = element_text(face = "bold", size = 14))

DDX3X_3Y_fantom_vert


pdf(file = "DDX3X_3Y_fantom.pdf", width = 20, height = 4)
print(DDX3X_3Y_fantom)
dev.off()


pdf(file = "DDX3X_3Y_fantom_vert.pdf", width = 6, height = 10)
print(DDX3X_3Y_fantom_vert)
dev.off()


##############Quantile normalization begins###############

temp_X <- subset_fantom[which(subset_fantom$Gene.name=="DDX3X"),]
temp_Y <- subset_fantom[which(subset_fantom$Gene.name=="DDX3Y"),]
rownames(temp_X)<- temp_X$Tissue
rownames(temp_Y)<- temp_Y$Tissue
temp_XY <- cbind(temp_X$Scaled.tags.per.million, temp_Y$Scaled.tags.per.million)
rownames(temp_XY) <- row.names(temp_X)
colnames(temp_XY) <- c("DDX3X", "DDX3Y")
temp_XY_rank <- apply(temp_XY,2,rank,ties.method="min")
temp_XY_sorted <- data.frame(apply(temp_XY, 2, sort))
temp_XY_mean <- apply(temp_XY_sorted, 1, mean)

index_to_mean <- function(my_index, my_mean){
  return(my_mean[my_index])
}

temp_XY_final <- apply(temp_XY_rank, 2, index_to_mean, my_mean=temp_XY_mean)
rownames(temp_XY_final)<- temp_Y$Tissue
colnames(temp_XY_final) <- c("DDX3X", "DDX3Y")
temp_XY_final$Tissue <- row.names()
temp_XY_final_melt <- melt((temp_XY_final))

DDX3X_3Y_fantom_quantnorm <- ggplot(data = temp_XY_final_melt, mapping = aes(x = Var1, 
                                                             y = Var2, 
                                                             fill = value)) +
  geom_tile() +
  xlab(label = element_blank()) +
  ylab(label = element_blank()) +
  facet_grid(~ subset_fantom$sex_spec, switch = "x", scales = "free_x", space = "free_x") +
  scale_fill_gradient(name = "Abundance",
                      low = "#FFFFFF",
                      high = "#31a354") +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = "#EEEEEE", color = "#FFFFFF")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, face = "bold", size = 12),
        axis.title.x = element_text(face = "bold", size = 12),
        strip.text.x = element_text(face = "bold", size = 12)) +
  theme(axis.text.y = element_text(angle = 90, hjust = 0.5, face = "bold", size = 12),
        axis.title.y = element_text(face = "bold", size = 12)) +
  theme(legend.text = element_text(face = "bold", size = 8),
        legend.title = element_text(face = "bold", size = 12))


pdf(file = "DDX3X_3Y_fantom_quantnorm.pdf", width = 20, height = 4)
print(DDX3X_3Y_fantom_quantnorm)
dev.off()

df2 <- temp_XY_final_melt[c(46:90),]
df3 <- temp_XY_final_melt[c(1:45),]
temp_XY_final_melt <- rbind(df2,df3)

subset_fantom$quantile_norm <- temp_XY_final_melt$value

DDX3X_3Y_fantom_vert_quantilenorm <- ggplot(data = subset_fantom, mapping = aes(x = Gene.name, 
                                                                   y = Tissue, 
                                                                   fill = quantile_norm)) +
  geom_tile() +
  xlab(label = element_blank()) +
  ylab(label = element_blank()) +
  facet_grid(subset_fantom$sex_spec ~ . , scales = "free_y", switch = "y", space = "free_y") +
  scale_fill_gradient(name = "Abundance",
                      low = "#FFFFFF",
                      high = "#31a354") +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = "#EEEEEE", color = "#FFFFFF")) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, face = "bold", size = 14),
        axis.title.x = element_text(face = "bold", size = 14),
        strip.text.x = element_text(face = "bold", size = 12)) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, face = "bold", size = 14),
        axis.title.y = element_text(face = "bold", size = 12),
        strip.text.y  = element_text(face = "bold", size = 9)) +
  theme(legend.text = element_text(face = "bold", size = 12),
        legend.title = element_text(face = "bold", size = 14))

pdf(file = "DDX3X_3Y_fantom_vert_quantilenorm.pdf", width = 6, height = 10)
print(DDX3X_3Y_fantom_vert_quantilenorm)
dev.off()

