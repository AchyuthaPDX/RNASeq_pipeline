# Differential Gene Expression Analysis with DESeq2

This project demonstrates a differential gene expression analysis using the DESeq2 package in R, along with data visualization using ggplot2 and an interactive Shiny application.

## Prerequisites

Before you begin, ensure you have the following packages installed:

- DESeq2
- tidyverse
- airway
- pheatmap
- ggplot2
- shiny

You can install these packages using the following commands:

```R
install.packages("tidyverse")
install.packages("pheatmap")
install.packages("ggplot2")
install.packages("shiny")
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("DESeq2")
BiocManager::install("airway")
```

Three main components include:
- Data Extraction
- Differential Gene Expression
- Interactive Visualization with Shiny

Acknowledgments
The airway package for providing sample data.
The DESeq2, tidyverse, pheatmap, ggplot2, and shiny packages for their functionalities in data analysis and visualization.
