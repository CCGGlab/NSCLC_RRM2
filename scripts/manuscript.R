library(ggplot2)
library(cowplot)
library(ggrepel)
library(ggvenn)
library(ggpubr)

# Load data and functions
source("scripts/functions/helpers.R")
load("data/NSCLCL_data.RData")

res_proc<- as.data.frame(res_diff_expr$CUTO9.24h_Lorlatinib)
res_proc<- res_proc[!is.na(res_proc$padj),]

# Decide which genes to label:
genes_to_label<- c(
  "ALK",
  "EML4",
  "RRM2",
  "DUSP5",
  "DUSP6",
  "ETV1",
  "ETV4",
  "ETV5",
  "FOSL1"
)

# Plot
p_volc<- plot_volcano(DE_results = res_proc, gene = genes_to_label, labelAll = T, labelCol = "black", useHyperbolicTH = T,logFC_cu = 1,curve = 0.1, plotTH=F, p_cu = 0.05, plot_nominal_p = F)
p_volc<- p_volc +
  theme(
    plot.title = element_text(hjust = 0.5, size=8), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.line = element_line(colour = "black", size=0.2),
    axis.ticks = element_line(colour = "black", size = 0.2),
    axis.text = element_text(size=6),  
    axis.title = element_text(size=8)
  ) +
  scale_x_continuous(name = "log2(Fold Change)", limits = c(min(res_proc$log2FoldChange),max(res_proc$log2FoldChange))) +
  scale_y_continuous(name = "-log10(Padj)", limits = c(0,max(-log10(res_proc$padj))))
 
# p_volc

# Venn diagram
genes_DE_ls <- list()
for(cl in c("CUTO8", "CUTO9", "CUTO29","YU1077")){
  res_proc<- as.data.frame(res_diff_expr[[paste0(cl,".24h_Lorlatinib")]])
  res_proc<- res_proc[!is.na(res_proc$log2FoldChange),]
  DE_genes <- get_DE(logFC = res_proc$log2FoldChange, P = res_proc$padj, genes = res_proc$HGNC, th_logFC = 1, th_logP = -log10(0.05), curve = 0.1)
  genes_DE_ls[[cl]] <- DE_genes$DE
}
genes_DE<- sapply(genes_DE_ls, function(x) x)

# Plot
p_venn <- ggvenn(
  genes_DE, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  stroke_size = 0.2, set_name_size = 8*0.36,
  text_size = 6*0.36,
  show_percentage = F
)
# p_venn

# Enrichment overlap venn
DE_lorla <- Reduce(intersect, genes_DE_ls)
GSEA<- do_GSEA2(genes_retrieved = DE_lorla, genes_all = res_diff_expr$CUTO29.24h_Alectinib$HGNC, GSEA_db = geneset_ls$Ha_ls, min_genes = 2, isList = T)
rownames(GSEA) <- gsub("HALLMARK_","",rownames(GSEA))
rownames(GSEA) <- gsub("EPITHELIAL_MESENCHYMAL_TRANSITION","EMT",rownames(GSEA))
GSEA$pathway<- factor(rownames(GSEA),levels=rev(rownames(GSEA)))
GSEA$q<- as.numeric(GSEA$q)
GSEA_sel<- GSEA[GSEA$q<0.01,]

# Barplot
p_gsea <- ggplot(GSEA_sel, aes(x = pathway, y = -log10(q))) +
  geom_bar(stat="identity", fill="darkblue") +
  coord_flip(expand = F) +
  ylab("-log10(Padj)") +
  xlab("") +
  theme(
    plot.title = element_text(hjust = 0.5, size=8), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.line = element_line(colour = "black", size=0.2),
    axis.ticks = element_line(colour = "black", size = 0.2),
    axis.text = element_text(size=6),  
    axis.title = element_text(size=8),
  )

GSEA_sel[,"q",drop=F]

# E2F targers in top, this includes RRM2
# E2F_TARGETS         7.736995e-35
# G2M_CHECKPOINT      6.677628e-21
# MYC_TARGETS_V1      4.530066e-05
# DNA_REPAIR          2.528127e-04
# MYC_TARGETS_V2      8.035498e-04
# KRAS_SIGNALING_UP   1.238968e-03
# EMT                 5.352815e-03
# IL2_STAT5_SIGNALING 5.352815e-03

# RRM2 boxplot
data <- normalized_counts[rownames(normalized_counts)=="RRM2",grep("24h",colnames(normalized_counts))]
data$Gene <- rownames(data)
data <- reshape2::melt(data)
data <- merge(data,sample_info,by.x = 2, by.y = 0)
data$Condition <- paste0(gsub("C29","CUTO29",gsub("YO","YU1077",gsub("C8","CUTO8",gsub("C9","CUTO9",data$CL)))),"+",data$treatment)
data$Condition <- gsub("\\+Control"," Control",data$Condition)
data$value <- log(data$value)

# Order samples correctly
sample_order <- c("CUTO8 Control","CUTO8+Lorlatinib","CUTO8+Alectinib","CUTO8+Brigatinib",
                  "CUTO9 Control","CUTO9+Lorlatinib","CUTO9+Alectinib","CUTO9+Brigatinib",
                  "CUTO29 Control","CUTO29+Lorlatinib","CUTO29+Alectinib","CUTO29+Brigatinib",
                  "YU1077 Control","YU1077+Lorlatinib","YU1077+Brigatinib")

level_order <- factor(data$Condition, level = sample_order)

# T-test comparisons
my_comparisons <- list( c("CUTO29+Lorlatinib", "CUTO29 Control"),
                        c("CUTO29+Alectinib", "CUTO29 Control"),
                        c("CUTO29+Brigatinib", "CUTO29 Control"),
                        
                        c("CUTO8+Lorlatinib", "CUTO8 Control"),
                        c("CUTO8+Alectinib", "CUTO8 Control"),
                        c("CUTO8+Brigatinib", "CUTO8 Control"),
                        
                        c("CUTO9+Lorlatinib", "CUTO9 Control"),
                        c("CUTO9+Alectinib", "CUTO9 Control"),
                        c("CUTO9+Brigatinib", "CUTO9 Control"),
                        
                        c("YU1077+Lorlatinib", "YU1077 Control"),
                        c("YU1077+Brigatinib", "YU1077 Control") )

# Correct treatment name
colnames(data)[6] <- "Treatment"

# Plot
p_boxplot <- ggplot(data, aes(x = factor(Condition, level = sample_order), y=value,fill = Treatment)) +
  geom_boxplot(lwd=0.3) +
  stat_compare_means(comparisons = my_comparisons,method = "t.test",size = 3,label = "p.signif",label.y = c(8,8.5,9,
                                                                                                            7.5,8,8.5,
                                                                                                            9,9.5,10,
                                                                                                            8.2,8.7)) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
    axis.text.y = element_text(size=6),
    axis.title.y = element_text(size = 8),
    axis.line = element_line(colour = "black", size=0.2),
    axis.ticks = element_line(colour = "black", size = 0.2),
    legend.key = element_rect(fill = "transparent", colour = "transparent"),
    legend.position = "top",
    axis.ticks.x = element_blank(),
    legend.title = element_blank(),
    legend.key.size = unit(5, 'mm')
  ) +
  labs(x = "", y = "RRM2 exp (log(normalized counts))") +
  scale_fill_manual(values=c("#C4961A","#D16103","#52854C","#4E84C4")) +
  scale_x_discrete(labels= c("CUTO8","","","","CUTO9","","","","CUTO29","","","","YU1077","","","")) 



p1 <- plot_grid(
  NULL,p_venn,p_gsea,NULL,
  nrow = 1,
  rel_widths = c(0.3,0.7,1,0.3)
)

p2 <- plot_grid(
  p_volc,p_boxplot,NULL,
  nrow = 1,
  rel_widths = c(1,1,0.15)
)

p <- plot_grid(
  p1,p2,NULL,
  nrow = 3,
  rel_heights = c(0.6,1.3,2)
)


# Save main plot
ggsave("results/figs/manuscript.pdf",  p, width = 178, height = 265, units = "mm")

