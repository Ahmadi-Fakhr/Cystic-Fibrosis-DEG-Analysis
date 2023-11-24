# Report

First, I obtained samples for Cystic Fibrosis from the GEO database (GSE141535). Then, I performed some data preprocessing by removing rows with a sum of counts less than 10. After that, I created the DESeq2 object and ran the function. Finally, I filtered the end results to include only the genes that were either **upregulated** or **downregulated** with a LFC (log2foldchange) of more than 1, and had a p-adjusted value of less than 0.05. This was done to identify genes that were differentially expressed in the disease group compared to the control group. The process was carried out to compare the disease group versus the control group. Eventually, I found 855 differentially expressed genes between the healthy and disease groups. These findings can be used for further analyses to generate biologically meaningful insights into this rare genetic disorder, Cystic Fibrosis.

## Tools

- ggplot2
- DESeq2

> More information can be found in the 'Report' file.

# File Guide
- `all_results`: Final raw/unfiltered output of the DESeq2 function.
- `counts_data`: Original dataset used.
- `filtered_results`: Filtered results based on p.adj < 0.05 & LFC > 1.
- `main.R`: Main R script.
- `top_DEGs`: Top DEG hits.
