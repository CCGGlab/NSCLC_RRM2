library(ggplot2)
library(cowplot)
library(ggdendro)
library(ggrepel)
# library(ggvenn)
library(WriteXLS)

# Load data and functions
source("scripts/functions/helpers.R")
load("data/NSCLCL_data.RData")

# Venn diagrams for alectinib and brigatinib
# p_venns <- list()
# for(tr in c("Brigatinib","Alectinib")){
#   DE_genes <- list()
#   for (cl in c("CUTO8", "CUTO9", "CUTO29","YU1077")) {
#     if(cl=="YU1077" & tr=="Alectinib") next
#     res_proc<- as.data.frame(res_diff_expr[[paste0(cl,".24h_",tr)]])
#     res_proc<- res_proc[!is.na(res_proc$log2FoldChange),]
#     DE_genes_temp <- get_DE(logFC = res_proc$log2FoldChange, P = res_proc$padj, genes = res_proc$HGNC, th_logFC = 1, th_logP = -log10(0.05), curve = 0.1)
#     DE_genes[[tr]][[cl]] <- DE_genes_temp$DE
#   }
#   
#   p_venns[[tr]] <- ggvenn(
#     DE_genes[[tr]],
#     fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
#     stroke_size = 0.2, set_name_size = 8*0.36,
#     text_size = 6*0.36,
#     show_percentage = F
#   ) + theme(plot.title = element_text(hjust = 0.5,size = 9),plot.margin=grid::unit(c(0,0,0,0), "mm"),)
#   print(tr)
# }

# GSEA for alectinib
# GSEA_genes <- Reduce(intersect, DE_genes$Alectinib)
# GSEA<- do_GSEA2(genes_retrieved = GSEA_genes, genes_all = res_diff_expr$CUTO29.24h_Alectinib$HGNC, GSEA_db = geneset_ls$Ha_ls, min_genes = 2, isList = T)
# rownames(GSEA) <- gsub("HALLMARK_","",rownames(GSEA))
# GSEA$pathway<- factor(rownames(GSEA),levels=rev(rownames(GSEA)))
# GSEA$q<- as.numeric(GSEA$q)
# 
# GSEA_sel<- GSEA[GSEA$q<0.01,]
# 
# # Barplot
# p_GSEA_alectinib <- ggplot(GSEA_sel, aes(x = pathway, y = -log10(q))) +
#   geom_bar(stat="identity", fill="darkblue") +
#   ggtitle("Alectinib") +
#   coord_flip(expand = F) +
#   ylab("-log10(Padj)") +
#   xlab("") +
#   theme(
#     plot.title = element_text(hjust = 0.5, size=8),
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     panel.background = element_blank(),
#     axis.line = element_line(colour = "black", size=0.2),
#     axis.ticks = element_line(colour = "black", size = 0.2),
#     axis.text = element_text(size=6),
#     axis.title = element_text(size=8),
#     plot.margin=grid::unit(c(0,0,0,0), "mm"),
#   )

# Enrichment heatmap for all significant Hallmark gene sets
res_diff_expr <- res_diff_expr[grep("24h",names(res_diff_expr))]
all_enr_tabs_gsea <- list()
for (i in 1:length(res_diff_expr)) {
  res_proc <- res_diff_expr[[i]]
  res_proc<- res_proc[!is.na(res_proc$log2FoldChange),]
  DE_genes <- get_DE(logFC = res_proc$log2FoldChange, P = res_proc$padj, genes = res_proc$HGNC, th_logFC = 1, th_logP = -log10(0.05), curve = 0.1)
  
  GSEA<- do_GSEA2(genes_retrieved = DE_genes$DE, genes_all = res_diff_expr$CUTO29.24h_Alectinib$HGNC, GSEA_db = geneset_ls$Ha_ls, min_genes = 2, isList = T)
  GSEA$q<- as.numeric(GSEA$q)
  GSEA <- GSEA[order(rownames(GSEA)),]
  all_enr_tabs_gsea[[names(res_diff_expr)[i]]] <- data.frame(Pathway = rownames(GSEA),p_val = GSEA$q)
  print(i)
}

all_enr_tabs_gsea <- as.data.frame(do.call(cbind,all_enr_tabs_gsea))
rownames(all_enr_tabs_gsea) <- gsub("EPITHELIAL_MESENCHYMAL_TRANSITION","EMT",gsub("HALLMARK_","",all_enr_tabs_gsea$CUTO29.24h_Alectinib.Pathway))

all_enr_tabs_gsea <- all_enr_tabs_gsea[,grep("p_val",colnames(all_enr_tabs_gsea))]
colnames(all_enr_tabs_gsea) <- gsub("\\.p_val","",gsub("\\.24h_","+",colnames(all_enr_tabs_gsea)))

# Removing rows where there is no significant hits
all_enr_tabs_gsea <- all_enr_tabs_gsea[all_enr_tabs_gsea< 0.01,]
all_enr_tabs_gsea <- all_enr_tabs_gsea[!is.na(rowSums(all_enr_tabs_gsea)),]

df1 <- all_enr_tabs_gsea
df1$Pathway <- rownames(all_enr_tabs_gsea)
df2 <- reshape2::melt(df1, id = "Pathway")

colnames(df2) <- c("Pathway","Fusion","p_val")
df2$`-log10(Padj)` <- -log10(df2$p_val)

dendro1 <- as.dendrogram(hclust(d = dist(x = -log10(all_enr_tabs_gsea))))
dendro2 <- as.dendrogram(hclust(d = dist(x = t(-log10(all_enr_tabs_gsea)))))
order1 <- order.dendrogram(dendro1)
order2 <- order.dendrogram(dendro2)

dendro.plot1 <- ggdendrogram(data = dendro1, rotate = TRUE,labels = FALSE,theme_dendro = T)
dendro.plot2 <- ggdendrogram(data = dendro2, rotate = FALSE,labels = FALSE,theme_dendro = T)

# Order rows by clustering from dendro
df2$Pathway <- factor(x = df2$Pathway,
                      levels = df1$Pathway[order1],
                      ordered = TRUE)

df2$text_col <- ifelse(df2$`-log10(Padj)`>20,"white","black")

p_heatmap <- ggplot(data = df2, aes(x = factor(Fusion, levels = colnames(all_enr_tabs_gsea)[order2]), y = Pathway)) +
  geom_tile(aes(fill = `-log10(Padj)`)) +
  geom_text(show.legend = FALSE,aes(label = signif(p_val, 2),color = text_col), size=2) + 
  scale_fill_gradient2(high = "darkblue",low = "red") +
  theme(axis.ticks.y = element_blank(),
        panel.background = element_blank(),
        legend.position = "bottom",
        legend.key.height = unit(0.4, 'cm'),
        plot.margin=grid::unit(c(0,0,0,0), "mm"),
        axis.text.x = element_text(angle = 45,  hjust=1),
        axis.text = element_text(size=8),  
        axis.title = element_text(size=8)) +
  labs(x="", y="") +
  scale_colour_manual(values=c("black", "white")) 


# Volcano plots
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

res_diff_expr <- res_diff_expr[grep("24h",names(res_diff_expr))]

p_volc <- list()
for (cl in c("CUTO8","CUTO29","YU1077")) {
  res_proc<- as.data.frame(res_diff_expr[[paste0(cl,".24h_Lorlatinib")]])
  res_proc<- res_proc[!is.na(res_proc$padj),]
  
  # Plot
  p_volc_tmp <- plot_volcano(DE_results = res_proc, gene = genes_to_label, labelAll = T, labelCol = "black", useHyperbolicTH = T,logFC_cu = 1,curve = 0.1, plotTH=F, p_cu = 0.05, plot_nominal_p = F)
  p_volc[[cl]] <- p_volc_tmp +
    theme(
      plot.title = element_text(hjust = 0.5, size=8), 
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      panel.background = element_blank(), 
      axis.line = element_line(colour = "black", size=0.2),
      axis.ticks = element_line(colour = "black", size = 0.2),
      axis.text = element_text(size=6),  
      axis.title = element_text(size=8),
      plot.margin=grid::unit(c(0,0,0,0), "mm"),
    ) +
    ggtitle(cl) +
    scale_x_continuous(name = "log2(Fold Change)", limits = c(min(res_proc$log2FoldChange),max(res_proc$log2FoldChange))) +
    scale_y_continuous(name = "-log10(Padj)", limits = c(0,max(-log10(res_proc$padj))))
}

# Merge plots
# p1 <- plot_grid(
#   p_GSEA_alectinib,p_venns$Brigatinib,
#   nrow = 2,
#   rel_heights = c(0.4,1)
# )
# 
# p2 <- plot_grid(
#   p_venns$Alectinib,p1,NULL,
#   nrow = 1
# )
# 
# p3 <- plot_grid(
#   p_volc$CUTO8,p_volc$CUTO29,p_volc$YU1077,
#   nrow = 1
# )
# 
# p <- plot_grid(
#   p2,p_heatmap,p3,NULL,
#   nrow = 4,
#   rel_heights = c(1,2,1.2,0.5)
# )

p0 <- plot_grid(
  p_volc$CUTO8,p_volc$CUTO29,p_volc$YU1077,
  nrow = 1
)

p <- plot_grid(
  p_heatmap,p0,NULL,
  nrow = 3,
  rel_heights = c(1.1,0.7,1)
)

ggsave("results/figs/manuscript_suppl.pdf",  p, width = 178, height = 265, units = "mm")

# Supplementary RNA DESeq result excel file

# Get only significant genes
cell_lines <- names(res_diff_expr)
res_diff_expr_sig <- list()
for (cl in cell_lines) {
  res_proc <- res_diff_expr[[cl]]
  res_proc<- res_proc[!is.na(res_proc$log2FoldChange),]
  DE_genes <- get_DE(logFC = res_proc$log2FoldChange, P = res_proc$padj, genes = res_proc$HGNC, th_logFC = 1, th_logP = -log10(0.05), curve = 0.1)
  res_proc <- res_proc[res_proc$HGNC %in% DE_genes$DE,]
  res_diff_expr_sig[[cl]] <- res_proc
}

# Save to excel sheet
WriteXLS(x = res_diff_expr_sig,
         ExcelFileName = "results/tables/manuscript_RNA_res_suppl.xls",
         SheetNames = names(res_diff_expr_sig),
         row.names = TRUE)
