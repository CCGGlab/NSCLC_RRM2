# Create sample information table based on available information in folder names

# Sample names
###############

# File names
fileNames <- list.files("raw/NSCLCL/",recursive = T)
fileNames1 <- gsub(".*\\/","",fileNames[grep("_1\\.fq.gz",fileNames)])
fileNames2 <- gsub(".*\\/","",fileNames[grep("_2\\.fq.gz",fileNames)])

fileNames <- fileNames[grep("_1\\.fq.gz",fileNames)]
fileNames <- gsub(".*\\/","",fileNames)

# CU8 and C8 are the same
sample_name<- gsub("_1.*","",fileNames)
CL <- rep(NA,length(sample_name))
CL[grep("CU8|C8",sample_name)] <- "C8"
CL[grep("CU9|C9",sample_name)] <- "C9"
CL[grep("CU29|C29",sample_name)] <- "C29"
CL[grep("CU29|C29",sample_name)] <- "C29"
CL[grep("YO",sample_name)] <- "YO"

# Time points
time <- rep(0,length(sample_name))
time[grep("6h",sample_name)] <- "6h"
time[grep("24h",sample_name)] <- "24h"

# Treatment
treatment <- rep(NA,length(sample_name))
treatment[grep("A",sample_name,ignore.case = F)] <- "Alectinib"
treatment[grep("L",sample_name,ignore.case = F)] <- "Lorlatinib"
treatment[grep("B",sample_name,ignore.case = F)] <- "Brigatinib"
treatment[grep("c",sample_name,ignore.case = F)] <- "Control"

# Fuse
sample_matrix<- data.frame(CL,time,treatment,fileNames1,fileNames2,row.names = sample_name)

# Save
saveRDS(sample_matrix,file="data/sample_info.rds")
