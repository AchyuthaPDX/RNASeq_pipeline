library(DESeq2)
library(tidyverse)
library(airway)


counts_data <- read.csv('counts_data.csv')
head(counts_data)
colData <- read.csv('sample_info.csv')
all(colnames(counts_data) %in% rownames(colData))

all(colnames(counts_data) == rownames(colData))




dds <- DESeqDataSetFromMatrix(countData = counts_data,
                              colData = colData,
                              design = ~ dexamethasone)

dds

# pre-filtering: removing rows with low gene counts
# keeping rows that have at least 10 reads total
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

dds

# set the factor level
dds$dexamethasone <- relevel(dds$dexamethasone, ref = "untreated")

dds <- DESeq(dds)
res <- results(dds)

res


summary(res)

res0.01 <- results(dds, alpha = 0.01)
summary(res0.01)

# contrasts
resultsNames(dds)


results(dds, contrast = c("dexamethasone", "treated", "untreated"))

# MA plot
plotMA(res)

library(DESeq2)
library(pheatmap)

sigGenes <- res[which(res$padj < 0.05), ]
normCounts <- assay(vst(dds[rownames(sigGenes), ]))
pheatmap(normCounts, cluster_rows = TRUE, cluster_cols = TRUE, scale = "row")


library(ggplot2)

res_df <- as.data.frame(res)


ggplot(res_df, aes(x=log2FoldChange, y=-log10(pvalue))) + 
  geom_point(alpha=0.3) + 
  theme_minimal() + 
  labs(x = "Log2 Fold Change", y = "-Log10 p-value", title = "Volcano Plot")


library(shiny)
library(DESeq2)
library(ggplot2)
library(pheatmap)


ui <- fluidPage(
  titlePanel("DESeq2 Visualization"),
  sidebarLayout(
    sidebarPanel(
      sliderInput("pvalSlider",
                  "Adjusted p-value threshold:",
                  min = 0,
                  max = 0.05,
                  value = 0.05,
                  step = 0.001),
      actionButton("update", "Update Charts")
    ),
    mainPanel(
      plotOutput("maPlot"),
      plotOutput("volcanoPlot"),
      uiOutput("heatmapPlot")
    )
  )
)

server <- function(input, output, session) {
  threshold <- reactiveVal(0.05)
  
  observeEvent(input$update, {
    threshold(input$pvalSlider)
  })
  
  # MA Plot
  output$maPlot <- renderPlot({
    res_adj <- results(dds, alpha = threshold())
    plotMA(res_adj, main = "MA Plot")
  })
  
  # Volcano Plot
  output$volcanoPlot <- renderPlot({
    res_df <- as.data.frame(results(dds, alpha = threshold()))
    ggplot(res_df, aes(x = log2FoldChange, y = -log10(pvalue))) + 
      geom_point(alpha = 0.3) + 
      theme_minimal() + 
      labs(x = "Log2 Fold Change", y = "-Log10 p-value", title = "Volcano Plot")
  })
  
  # Heatmap Plot (as an image, due to the complexity of the plot)
  output$heatmapPlot <- renderUI({
    sigGenes <- results(dds, alpha = threshold())[which(res$padj < threshold()), ]
    # Normalize counts
    normCounts <- assay(vst(dds[rownames(sigGenes), ]))
    
    imgPath <- paste0(tempfile(), ".png")
    
    # Save the heatmap to a temporary file
    png(imgPath, width = 800, height = 600)
    pheatmap(normCounts, cluster_rows = TRUE, cluster_cols = TRUE, scale = "row")
    dev.off()
    
    # Render the saved image
    tags$img(src = imgPath, height = 600, width = 800, style = "display:block; margin:auto;")
  })
}
# Run the application
shinyApp(ui = ui, server = server)
