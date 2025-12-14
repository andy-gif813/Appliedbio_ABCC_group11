# DESeq2_basic_analysis.R

# 1. Load packages
library(BiocManager)
#BiocManager::install("DESeq2")
library(DESeq2)
library(pheatmap)
library(RColorBrewer)

# 2. Read data
counts <- read.csv("~/Desktop/count_matrix_raw.csv", header = TRUE, row.names = 1)
sample_info <- read.table("~/Desktop/sample_info.txt", header = TRUE, row.names = 1)
orf<-read.csv("~/Desktop/ORF.csv",header = FALSE)
# Round expression counts to integers
counts_int <- round(counts)
# Set wild_type as the reference level
sample_info$condition <- factor(sample_info$condition, 
                                levels = c("engineered","wild_type"))
# Extract Verified ORFs from SGD
f_orf<-orf[grepl("Verified", orf[, 2]), ]

# 3. Create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = counts_int,
                              colData = sample_info,
                              design = ~ condition)

# 4. Filter low-expression genes
keep <- rowSums(counts(dds)) >0
dds <- dds[keep, ]

# 5. Run differential analysis
dds <- DESeq(dds)

# 6. Get results
res <- results(dds, contrast = c("condition","engineered","wild_type"),alpha = 0.001)


# 7. Extract significant genes
sig_genes <- subset(res, padj < 0.001 & abs(log2FoldChange) >=1)
# Extract Verified ORF genes
sig_genes<-data.frame(sig_genes)
a<-sig_genes[order(sig_genes$pvalue),]
sig_genes<- cbind(ID = rownames(sig_genes), sig_genes)
o_sig_genes <- merge(sig_genes,f_orf, by.x = colnames(sig_genes)[1], 
                     by.y = colnames(f_orf)[1])

# 8. Save results
#write.csv(as.data.frame(sig_genes), "~/Desktop/DESeq2_results.csv")
#write.csv(o_sig_genes, "~/Desktop/Verified_ORF_results.csv")

# 9. Generate heatmap
# Normalize expression data
dds_rlog <- rlog(dds, blind = FALSE)
rlog_counts <- assay(dds_rlog)

# Create color gradient
my_col <- colorRampPalette(c("blue","cyan","green","yellow", "orange","red"))(100)

# Plot gene expression heatmap
pheatmap(rlog_counts,
         color = my_col,  
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         display_numbers = FALSE,
         show_rownames = FALSE, 
         fontsize = 15,
         breaks = seq(min(rlog_counts), max(rlog_counts), length.out = 100)
         )


#10 Calculate Pearson correlation between samples
cor_matrix <- cor(counts, method = "pearson")
#write.csv(cor_matrix, "~/Desktop/pearson.csv")

cat("Total number of genes:", nrow(res), "\n")
cat("Number of significant differentially expressed genes:", nrow(sig_genes), "\n")
cat("Number of Verified ORF genes:", nrow(o_sig_genes), "\n")

