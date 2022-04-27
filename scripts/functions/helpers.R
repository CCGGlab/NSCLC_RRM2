library(plotly)
library(fgsea)
library(dplyr)

do_GSEA2<- function(genes_retrieved,genes_all,GSEA_db,min_genes,isList=FALSE){
  if(isList) GSEA_names<-names(GSEA_db)
  else GSEA_names<-rownames(GSEA_db)
  GSEA_table<-matrix(NA, length(GSEA_names), 9, dimnames=list(GSEA_names,c("n_genes_pw","n_genes_pw_neg","n_genes_pw_pos","prop_neg","prop_pos","OR","p","q","genes")))
  for(i in 1:nrow(GSEA_table)){
    # cat(i," ")
    if(isList) genes_temp<-unlist(GSEA_db[[i]])
    else{
      genes_temp<-unique(as.character(GSEA_db[i,-1]))
      genes_temp<-genes_temp[genes_temp!=""]
    }
    # genes_overlap<- genes_temp[genes_temp%in%genes_all]
    genes_overlap<- intersect(genes_temp,genes_all)
    GSEA_table[i,"n_genes_pw"]<-length(genes_overlap)
    if(length(genes_overlap)<min_genes) next #No use analysing 
    else{
      isRetrieved<- genes_all%in%genes_retrieved  
      inDB=factor(genes_all%in%genes_overlap,levels=c(FALSE,TRUE))
      compare_pw_temp_t<-table(isRetrieved,inDB)
      GSEA_table[i,c("n_genes_pw_neg","n_genes_pw_pos")]<-compare_pw_temp_t[,"TRUE"]
      GSEA_table[i,c("prop_neg","prop_pos")]<-round(100*prop.table(compare_pw_temp_t,1)[,"TRUE"],1)
      GSEA_table[i,"genes"]<-paste(genes_all[isRetrieved&inDB==TRUE],collapse=",")
      GSEA_table[i,"OR"]<-fisher.test(compare_pw_temp_t,alternative = "greater")$estimate   
      GSEA_table[i,"p"]<-fisher.test(compare_pw_temp_t,alternative = "greater")$p.value   
    }
  }
  GSEA_table<-GSEA_table[!is.na(GSEA_table[,"p"])&!duplicated(rownames(GSEA_table)),]
  # GSEA_table[,"q"]<- qvalue(as.numeric(GSEA_table[,"p"]))$qvalues # To use Storey method (library qvalue)
  GSEA_table[,"q"]<-p.adjust(as.numeric(GSEA_table[,"p"]),"fdr")
  GSEA_table<-as.data.frame(GSEA_table[order(as.numeric(GSEA_table[,"p"])),])
  return(GSEA_table)
}

get_DE<- function(logFC,P,genes,th_logFC=1,th_logP= -log10(0.05),curve=1){
  idx_DE<- which(-log10(P)>sapply(logFC,function(x) get_hyperbolic_th_P(x,th_logFC=th_logFC,th_logP= th_logP,curve=curve)))
  idx_up<- intersect(idx_DE,which(logFC>0))
  idx_down<- intersect(idx_DE,which(logFC<0))
  genes_up<- genes[idx_up]
  genes_down<- genes[idx_down]
  genes_DE<- genes[idx_DE]
  res<- list(up=genes_up,down=genes_down,DE=genes_DE)
  return(res)
}

get_hyperbolic_th_P<- function(logFC,th_logFC=1,th_logP= -log10(0.05),curve=1){ 
  if(abs(logFC)<th_logFC) res<- Inf
  else res<- th_logP+(curve/sqrt(logFC^2-th_logFC^2))
  # if(logFC>=0 ) res<- th_logP+(1/(logFC-th_logFC))
  # else res<- th_logP+(1/(logFC+th_logFC))
  return(res)
}

plot_volcano<- function(DE_results, gene, useHyperbolicTH=F, labelAll=F, labelCol="red", p_cu=0.01, logFC_cu=2, curve=0.5, plotTH=T, plot_nominal_p=F, th_nominal_p=F){
  
  # Get selected data
  res_proc<- DE_results
  res_proc<- res_proc[!is.na(res_proc[,"padj"]),] # remove NA (i.e. no expression before/after treatment)
  res_proc[res_proc[,"padj"]==0,"padj"]<- min(res_proc[res_proc[,"padj"]!=0,"padj"])/100 # if q=0, put at lowest possible value
  res_proc[res_proc[,"pvalue"]==0,"pvalue"]<- min(res_proc[res_proc[,"pvalue"]!=0,"pvalue"])/100 # if q=0, put at lowest possible value
  
  # DE thresholds
  if(th_nominal_p){
    if(useHyperbolicTH) res_proc$isDE<- rownames(res_proc)%in%get_DE(res_proc$log2FoldChange,res_proc$pvalue,rownames(res_proc),logFC_cu,-log10(p_cu),curve)$DE
    else res_proc$isDE<- abs(res_proc$log2FoldChange)>=logFC_cu & -log10(res_proc$pvalue)>= -log10(p_cu)
  }
  else{
    if(useHyperbolicTH) res_proc$isDE<- rownames(res_proc)%in%get_DE(res_proc$log2FoldChange,res_proc$padj,rownames(res_proc),logFC_cu,-log10(p_cu),curve)$DE
    else res_proc$isDE<- abs(res_proc$log2FoldChange)>=logFC_cu & -log10(res_proc$padj)>= -log10(p_cu)
  }
  
  # Gene ids
  res_proc$gene_id<- rownames(res_proc)
  res_proc$isLabel<- res_proc$gene_id%in%gene
  
  # Plot
  res_proc<- res_proc[order(res_proc$isDE),] # Blue before grey
  if(plot_nominal_p) p<- ggplot(res_proc, aes(x = log2FoldChange, y = -log10(pvalue), color= isDE, key=gene_id))
  else p<- ggplot(res_proc, aes(x = log2FoldChange, y = -log10(padj), color= isDE, key=gene_id))
  p<- p +
    geom_point() +
    scale_color_manual(values = c("#d0d3d4", "#3498db")) +
    theme(legend.position = "none")
  
  if(length(gene)>=1) p<- p+ geom_point(data=res_proc[gene,], colour=labelCol) # this adds a red point
  
  if(labelAll==T){
    p<- p +
      geom_text_repel(data = subset(res_proc, res_proc$isLabel),
                      aes(label = gene_id, fontface=3),
                      size = 2,
                      colour="black",
                      box.padding   = 0.5,
                      point.padding = 0.5,
                      max.overlaps = 25,
                      segment.color = 'grey50')
    # p<- p+ geom_text(data=res_proc[gene,], label=gene, vjust=0, hjust=0, colour=labelCol) # this adds a label for the red point
  }
  
  # Hyperbolic treshold?
  if(plotTH==T){
    if(useHyperbolicTH){
      x_lim = c(round(min(res_proc$log2FoldChange)-1),round(max(res_proc$log2FoldChange)+1))
      if(plot_nominal_p) y_lim = c(0,round(max(-log10(res_proc$pvalue))+1))
      else y_lim = c(0,round(max(-log10(res_proc$padj))+1))
      th_data<- data.frame(
        x=seq(x_lim[1],x_lim[2],0.01),
        y=get_hyperbolic_th_P(seq(x_lim[1],x_lim[2],0.01),th_logFC = logFC_cu,th_logP = -log10(p_cu),curve = curve),
        gene_id="x",
        isDE=F
      )
      th_data<- th_data[th_data$y<y_lim[2],]
      p<- p +
        geom_line(data=th_data[th_data$x<0,],aes(x=x,y=y),linetype = "dashed",colour="grey20") +
        geom_line(data=th_data[th_data$x>0,],aes(x=x,y=y),linetype = "dashed",colour="grey20")
    }
    else{
      p<- p +
        geom_hline(yintercept = c(-log10(p_cu)), linetype = "dashed",colour="grey20") +
        geom_vline(xintercept = c(-logFC_cu,logFC_cu), linetype = "dashed",colour="grey20")
    }
  }
  
  # return
  return(p)
}


# Script to customize survival plots (ggsurvplot)
customize_labels <- function (p, font.title = NULL,
                              font.subtitle = NULL, font.caption = NULL,
                              font.x = NULL, font.y = NULL, font.xtickslab = NULL, font.ytickslab = NULL)
{
  original.p <- p
  if(is.ggplot(original.p)) list.plots <- list(original.p)
  else if(is.list(original.p)) list.plots <- original.p
  else stop("Can't handle an object of class ", class (original.p))
  .set_font <- function(font){
    font <- ggpubr:::.parse_font(font)
    ggtext::element_markdown (size = font$size, face = font$face, colour = font$color)
  }
  for(i in 1:length(list.plots)){
    p <- list.plots[[i]]
    if(is.ggplot(p)){
      if (!is.null(font.title)) p <- p + theme(plot.title = .set_font(font.title))
      if (!is.null(font.subtitle)) p <- p + theme(plot.subtitle = .set_font(font.subtitle))
      if (!is.null(font.caption)) p <- p + theme(plot.caption = .set_font(font.caption))
      if (!is.null(font.x)) p <- p + theme(axis.title.x = .set_font(font.x))
      if (!is.null(font.y)) p <- p + theme(axis.title.y = .set_font(font.y))
      if (!is.null(font.xtickslab)) p <- p + theme(axis.text.x = .set_font(font.xtickslab))
      if (!is.null(font.ytickslab)) p <- p + theme(axis.text.y = .set_font(font.ytickslab))
      list.plots[[i]] <- p
    }
  }
  if(is.ggplot(original.p)) list.plots[[1]]
  else list.plots
}

get_GSEA_stat <- function(DE_results, isMice=F, MGI_to_HGNC=NULL){
  res_proc<- DE_results
  if(!"stat" %in% colnames(res_proc)) res_proc$stat<- res_proc$log2FoldChange # Use logFC to rank
  stat<- res_proc$stat + 10e-10*res_proc$log2FoldChange # Addition ro prevent ties
  names(stat)<- rownames(res_proc)
  
  if(isMice){
    MGI_to_HGNC<- MGI_to_HGNC[!duplicated(MGI_to_HGNC$MGI.symbol),]
    rownames(MGI_to_HGNC)<- MGI_to_HGNC$MGI.symbol
    names(stat)<- MGI_to_HGNC[names(stat),"HGNC.symbol"]
  }
  
  stat<- stat[!duplicated(names(stat))]
  stat<- sort(stat,decreasing = T)
  stat
}

do_fGSEA <- function(db, stat, minSize=15,maxSize=500){
  Res<- fgsea::fgseaMultilevel(pathways = db, 
                               stats = stat,
                               minSize=minSize,
                               maxSize=maxSize,
                               eps=0)    
  Res<- Res[order(Res$pval),]
  Res
}