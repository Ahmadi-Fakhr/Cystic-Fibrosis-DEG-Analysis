library(DESeq2)
library(ggplot2)

counts_data <- read.delim('counts_data.csv', sep = ',', header = TRUE, row.names = 1)
condition <- factor(c("Control", "Control", "Control", "Control", "Disease", "Disease", "Disease", "Disease"))
coldata <- data.frame(row.names = colnames(counts_data), condition)
Groups <- make.names(c("Control", "Disease"))

dds <- DESeqDataSetFromMatrix(countData = counts_data,
                              colData = coldata,
                              design = ~ condition)

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]
dds <- DESeq(dds)
r <- results(dds)
r <- subset(r, padj < 0.05 & abs(log2FoldChange) > 1)

# top
top_ones <- r[order(r$padj)[1:10], ]
top_ones

### Plots
r <- results(dds)
old.pal <- palette(c("#00BFFF", "#FF3030"))
par(mar=c(4,4,2,1), cex.main=1.5)
plot(r$log2FoldChange, -log10(r$padj), main=paste(Groups[1], "vs", Groups[2]),
     xlab="log2FC", ylab="-log10(Padj)", pch=20, cex=0.5)
with(subset(r, padj<0.05 & abs(log2FoldChange) >= 0),
     points(log2FoldChange, -log10(padj), pch=20, col=(sign(log2FoldChange) + 3)/2, cex=1))
legend("topright", title=paste("Padj<", 0.05, sep=""), legend=c("down", "up"), pch=20,col=1:2)

par(mar=c(4,4,2,1), cex.main=1.5)
plot(log10(r$baseMean), r$log2FoldChange, main=paste(Groups[1], "vs", Groups[2]),
     xlab="log10(mean of normalized counts)", ylab="log2FoldChange", pch=20, cex=0.5)
with(subset(r, padj<0.05 & abs(log2FoldChange) >= 0),
     points(log10(baseMean), log2FoldChange, pch=20, col=(sign(log2FoldChange) + 3)/2, cex=1))
legend("topright", title=paste("Padj<", 0.05, sep=""), legend=c("down", "up"), pch=20,col=1:2)
abline(h=0)
palette(old.pal)
