library(DESeq2)
library(biomaRt)
# Get ENSG_HGNC map

# Get ENSG genenames fom ht-seq-count output
sampleFiles <- grep("counts",list.files("raw/gene_counts",full.names = T), value=TRUE)
genenames_ENSG<- read.table(sampleFiles[1],stringsAsFactors = F)[,1]

ensembl_hs <- useEnsembl(biomart="ensembl",dataset="hsapiens_gene_ensembl") # according to https://support.bioconductor.org/p/74906/#74957 
ENSG_HGNC_map_table<- getBM(attributes = c("ensembl_gene_id","ensembl_gene_id_version","external_gene_name","gene_biotype"),filters = "ensembl_gene_id_version",values=genenames_ENSG,mart=ensembl_hs,uniqueRows = TRUE)
rownames(ENSG_HGNC_map_table)<- ENSG_HGNC_map_table[,"ensembl_gene_id_version"]  

# Load sample information
sample_information<- readRDS(file="data/sample_info.rds")

# Create column for different groups
sample_information$group <- paste(sample_information$CL,sample_information$time,sample_information$treatment, sep = ".")

# Get quantified files
sampleNames<- rownames(sample_information)
sampleFiles<- paste0("raw/gene_counts/",sample_information$fileNames1,".gene_counts")

# Add condition from sample info and create sample table
sampleCondition <- as.factor(sample_information$group)
sampleTable<-data.frame(sampleName=sampleNames, fileName=sampleFiles, condition=sampleCondition)

# From this phase combine with sample information!
ddsHTSeq<- DESeqDataSetFromHTSeqCount(sampleTable=sampleTable, design=~condition)
colData(ddsHTSeq)$condition<- factor(colData(ddsHTSeq)$condition, levels=levels(sampleTable$condition))
dds<- DESeq(ddsHTSeq)
normalized_counts <- counts(dds, normalized=T)
rownames(normalized_counts) <- gsub("\\..*","",rownames(normalized_counts))
normalized_counts <- as.data.frame(normalized_counts)

normalized_counts <- merge(ENSG_HGNC_map_table,normalized_counts,by.x = 1, by.y = 0)
normalized_counts <- normalized_counts[normalized_counts$gene_biotype=="protein_coding",]
normalized_counts <- normalized_counts[!duplicated(normalized_counts$external_gene_name),]
rownames(normalized_counts) <- normalized_counts$external_gene_name
normalized_counts <- normalized_counts[,-c(1:4)]

cl_time <- unique(gsub("\\.A.*","",gsub("\\.L.*","",gsub("\\.B.*","",gsub("\\.C.*","",sample_information$group)))))
treatment <- unique(gsub(".*\\.","",sample_information$group))[c(1,2,4)]

res_proc_all <- list()
for (i in 1:length(cl_time)) {
  for (j in 1:length(treatment)) {
    # YO fusion missing Alectinib treatment, so skipping this one
    if (any(i==7,i==8) & j==1) next
    print(c("condition",paste(cl_time[i],treatment[j],sep = "."),paste(cl_time[i],"Control",sep = ".")))
    # Do test for all treatments against control at all time points
    res_proc <- results(dds, contrast = c("condition",paste(cl_time[i],treatment[j],sep = "."),paste(cl_time[i],"Control",sep = ".")))
    
    # Get only protein coding genes
    res_proc <- merge(as.data.frame(res_proc),ENSG_HGNC_map_table,by = "row.names")
    res_proc <- res_proc[res_proc$gene_biotype == "protein_coding",]
    
    # Convert to HGNC names
    # rownames(res_proc) <- res_proc$external_gene_name
    res_proc$HGNC <- res_proc$external_gene_name
    res_proc$ENSG <- res_proc$ensembl_gene_id_version
    res_proc <- res_proc[,colnames(res_proc) %in% c("baseMean","log2FoldChange","lfcSE","stat","pvalue","padj","HGNC","ENSG")]
    res_proc <- res_proc[!duplicated(res_proc$HGNC),]
    rownames(res_proc) <- res_proc$HGNC
    
    # Sort on genenames
    res_proc<- res_proc[order(res_proc[,"HGNC"]),]
    
    # Redo padj for selected (coding) genes
    res_proc$padj<- p.adjust(res_proc$pvalue,"fdr")
    res_proc_all[[paste(cl_time[i],treatment[j],sep = "_")]] <- res_proc
  }
}

# Order list of data frames
res_proc_all <- res_proc_all[order(names(res_proc_all))]
# It's not YO, it's YU
names(res_proc_all) <- gsub("YO","YU1077",gsub("C9","CUTO9",gsub("C8","CUTO8",gsub("C29","CUTO29",names(res_proc_all)))))

# Save
saveRDS(normalized_counts,file = "data/normalized_counts.rds")
saveRDS(res_proc_all,file = "data/res_diff_expr.rds")


