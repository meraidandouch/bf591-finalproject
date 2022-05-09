#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(ggplot2)
library(colourpicker) 
library(shinyjs)
library(tidyverse)
library(gridExtra)
library(limma)
library(edgeR)
library(DT)
library(RColorBrewer)
library(pheatmap)
library(gplots)
library(SummarizedExperiment)
library(DESeq2)
library(biomaRt)
library(testthat)
library(fgsea)
# library(STRINGdb)

# Define UI for application that draws a histogram
ui <- fluidPage(
      navbarPage(
        "BF591 Final Project: Bioinformatics Processes ",
        tabPanel("Samples",
            sidebarPanel(
                   fileInput("sample_data", "Load sample results", accept = ".csv"),
                   markdown("After uploading a samples.csv file, you can explore the data further by viewing the summary table, entire sortable data-table and density plot.")),
                 mainPanel(
                   tabsetPanel(
                     tabPanel("Summary",
                              h4("This tab provides a summary of the metadata"),
                              tableOutput("sample_tableInfo")),
                     tabPanel("Data",
                              h4("Below shows the entire datatable for metadata file."),
                              DT::dataTableOutput("sampleTable")),
                     tabPanel("Plots", 
                              h4("Density plot with multiple lines for total clean paired reads for each genotype"),
                              plotOutput("histogram")))
            )
        ),
        tabPanel("Counts", 
            sidebarPanel(
                   fileInput("count_data", "Load count data results", accept = ".csv"),
                   markdown("Enter your raw read counts data, and filter the data using variance and non-zero variables."), 
                   br(),
                   markdown("For example, you might want to see genes that have more variance than **20%** of the data. 
                       And genes that have at least **20** nonzero samples. Try it out!"),
                     sliderInput("slider1",
                                 "Select to include only genes with at least X percentile of variance:",
                                 min = 0,
                                 max = 1,
                                 value = .20),
                     sliderInput("slider2",
                               "Select to include genes with at least X samples that are non-zero:",
                               min = 0,
                               max = 79,
                               value = 20),
                    submitButton("Submit")),
            mainPanel(
                    tabsetPanel(
                      tabPanel("Summary",
                               h4("This tab provides a summary of genes that are passing/not passing after the filters are applied."),
                               tableOutput("countsTable")),
                      tabPanel("Scatterplots",
                               h4("Diagnostic scatter plots are shown below to provide a visual for variance and number of zeros in the counts data."),
                               plotOutput("plotvarmean")),
                      tabPanel("Clustered Heatmap",
                               h4("Below shows a heatmap figure for the top 30 highly variable genes using log counts of normalized cpm data."),
                               plotOutput("plothm")),
                      tabPanel("PCA",
                               h4("Two-dimensional projection of raw counts data using PCA"),
                               selectInput("dd1",
                                           label = "Choose 1st PCA:",
                                           choice = c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "P7")),
                               selectInput("dd2",
                                           label = "Choose 2nd PCA:",
                                           choice = c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7")),
                               submitButton("submit"),
                               plotOutput("pcaplot"))))
        ),
        tabPanel("DE",
            sidebarPanel(
                    fileInput("deseq_data", "Load differential expression results", accept = ".csv"),
                    markdown ("A volcano plot can be generated with **log2 fold-change** on the x-axis and **p-adjusted** on the y-axis."),
                    radioButtons("button1", "Choose the column for the x-axis",
                                 c("baseMean", "log2FoldChange",'lfcSE', 'stat', 
                                   'pvalue', 'padj')),
                    radioButtons("button2", "Choose the column for the y-axis",
                                 c("baseMean", "log2FoldChange",'lfcSE', 'stat', 
                                   'pvalue', 'padj')),
                    colourpicker::colourInput("col1", "Base point color", "#AAC8E6"),
                    colourpicker::colourInput("col2", "Highlight point color", "#000000"),
                    sliderInput("slider3",
                                "Select the magnitude of the p adjusted coloring:",
                                min = -300,
                                max = 0,
                                value = -150),
                    submitButton("submit")),
            mainPanel(
                  tabsetPanel(
                    tabPanel("Data",
                              h4("This tab shows all rows from DE results"),
                              DT::dataTableOutput("de_summ_table")),
                    tabPanel("Table",
                             h4("This tab summarizes DEG after filtering"),
                             tableOutput("table")),
                    tabPanel("Plots",
                              h4("This tab shows a plot using two variables of your choice!"),
                             plotOutput("volcano"))))
        ),
        tabPanel("GSEA", 
            sidebarPanel(
                   fileInput("fgsea_data", "Load fgsea pathway enrichment results", accept = ".csv"),
                   markdown("Enter GSEA data, and filter the data using p-adjusted values."),
                   sliderInput("slider4",
                               "Select the magnitude of the p adjusted coloring:",
                               min = 0,
                               max = 1,
                               value = 0.001),
                   submitButton("Submit")),
            mainPanel(
                   tabsetPanel(
                     tabPanel("Barplot",
                              h4("Barplot of fgsea NES for top pathways selected by slider"),
                              plotOutput("fgsea_barplot")),
                     tabPanel("Table",
                              h4("Filtered table by adjusted p-value"),
                              markdown("Please upload your file first before applying changes or pressing the download button!"),
                              radioButtons("NES_button", "Choose to select all postive or negtive NES pathways",
                                           c("postive", "negative")),
                              submitButton("Apply Changes"),
                              downloadButton('download_fgsea', "Download csv"),
                              tableOutput("fgsea_table")),
                     tabPanel("Scatterplot",
                              h4("Scatter plot of NES on x-axis and -log10 adj-pval"),
                              plotOutput("fgsea_scatterplot"))))
        ),
    )
)

# Define server logic required to draw a histogram
options(shiny.maxRequestSize = 30*1024^2)
server <- function(input, output, session) {
  #' load_Data
  load_sampledata <- reactive({
    req(input$sample_data$datapath)
    data <- read_csv(input$sample_data$datapath)
    names(data) <- make.names(names(data))
    return(data)
  })
  
  load_countdata <- reactive({
    req(input$count_data$datapath)
    data <- read_csv(input$count_data$datapath)
    data <- subset(data, select = -c(1, 3)) #remove extra index column and gene_symbols column
    data <- data[!duplicated(data), ] #remove duplicates
    return(data)
  })
  
  load_deseqdata <- reactive({
    req(input$deseq_data$datapath)
    data <- read_csv(input$deseq_data$datapath)
    return(data)
  })
  
  load_fgseadata <- reactive({
    req(input$fgsea_data$datapath)
    data <- read_csv(input$fgsea_data$datapath)
    return(data)
  })
  
  fgsea_filtered <- reactive({
    draw_fgseatable(load_fgseadata(), input$NES_button, input$slider4)
  })
  
  
  #' Volcano plot
  #'
  #' @param dataf The loaded data frame.
  #' @param x_name The column name to plot on the x-axis
  #' @param y_name The column name to plot on the y-axis
  #' @param slider A negative integer value representing the magnitude of
  #' p-adjusted values to color. Most of our data will be between -1 and -300.
  #' @param color1 One of the colors for the points.
  #' @param color2 The other colors for the points. Hexadecimal strings: "#CDC4B5"
  #'
  #' @return A ggplot object of a volcano plot
  #' @details I bet you're tired of these plots by now. Me too, don't worry.
  #' This is _just_ a normal function. No reactivity, no bells, no whistles. 
  #' Write a normal volcano plot using geom_point, and integrate all the above 
  #' values into it as shown in the example app. The testing script will treat 
  #' this as a normal function.
  #' 
  #' !!sym() may be required to access column names in ggplot aes().
  #'
  #' @examples volcano_plot(df, "log2fc", "padj", -100, "blue", "taupe")
  volcano_plot <- function(dataf, x_name, y_name, slider, color1, color2) {
    vp <- dataf %>% as_tibble() %>% 
      ggplot(mapping=aes(x=!!sym(x_name), y=-log10(!!sym(y_name)))) + 
      geom_point(aes(color = padj< 1*10^(slider)))+
      labs( color = str_glue('{y_name} < 1x10^{slider}'))+ 
      theme_minimal() +
      scale_color_manual(values=c(color1, color2))  +
      theme(legend.position="bottom") + 
      ggtitle("Volcano plot of DESeq2 differential expression results")
    return(vp)
  }
  
  #' Density/Histogram plot
  #'
  #' @param dataf The loaded data frame.
  #' 
  #' @return A ggplot object of a density plot 
  #'
  #' @examples histo_plot(data)
  histo_plot <- function(data) {
    hp <- data %>% as_tibble() %>% 
      ggplot(aes(x=clean.paired.reads..CP., color=Genotype)) + 
      geom_density() + 
      ggtitle("Density plot of cleaned raw count data against Genotype") + 
      theme(axis.text.x = element_text(angle = 90, hjust = 1))
    return(hp)
  }
  #' plot_variance_vs_mean
  #'
  #' @param data tibble: a (n x _S_) data set
  #' @param slider_cp Slider to include genes with at least X percentile of variance
  #' @param slider_cs Slider to include genes with at least X samples that are non-zero
  #'
  #' @return ggplot: A scatter plot where the x-axis is the rank of gene ordered by mean
  #' count over all samples, and the y-axis is the observed variance of the
  #' given gene. Each dot should have their transparency increased. The scatter
  #' plot should also be accompanied by a line representing the average median and
  #' variance values.
  #'
  plot_variance_vs_mean <- function(data, slider_cp, slider_cs) {
    #validate(need(load_countdata))
    d <- data[,2:ncol(data)]
    threshold <- quantile(apply(d, 1, var), probs = c(slider_cp)) %>% round(3)
    median <- apply(d, 1, median)
    variances <- apply(d, 1, var) #apply to rows
    plot_data <- tibble::tibble(median=median, variance=variances)
    plot_data$rank <- rank(plot_data$median)
    mv_plot <- plot_data %>% 
      ggplot(aes(x=rank, y=variances, color=variances>threshold)) +
      labs(color = str_glue('variances > {slider_cp}'))+
      ggplot2::geom_point(alpha=0.5) +
      ggplot2::geom_smooth(method = 'gam', formula = y ~ s(x, bs = "cs")) + 
      ggplot2::xlab("Rank(Median)") + 
      ggplot2::ylab("Variance") + 
      ggplot2::ggtitle("Variance vs. Median (Log10(Raw Counts))") + 
      scale_color_manual(values=c("gray", "black"))
    #scale y axis
    mv_plot <- mv_plot + ggplot2::scale_y_log10()
    mv_plot
    
    nonzero <- rowSums(d!=0) #number of zeros
    plot_data2 <- tibble::tibble(median=median, nonzeros=nonzero)
    plot_data2$rank <- rank(plot_data2$median)
    mv_plot2 <- plot_data2 %>% 
      ggplot(aes(x=rank, y=nonzeros, color=nonzeros > slider_cs)) +
      labs(color = str_glue('Nonzeros > {slider_cs}'))+
      ggplot2::geom_point(alpha=0.5) +
      ggplot2::xlab("Rank(Median)") + 
      ggplot2::ylab("NonZeros") + 
      ggplot2::ggtitle("NonZeros vs. Median(Raw Counts)") + 
      scale_color_manual(values=c("gray", "black"))
    mv_plot2
    
    mv<-grid.arrange(mv_plot, mv_plot2)
    return(mv)
  }
  
  #' plot_heatmap: a function that takes the intensity values for significant probes and
  #' creates a color-blind friendly heatmap
  #'
  #' @param data (matrix): The matrix of raw read count values 
  #' @param slider_cp Slider to include genes with at least X percentile of variance
  #' @param slider_cs Slider to include genes with at least X samples that are non-zero
  #' 
  #' @return A heatmap displaying the intensity values for the highly variable counts
  #'
  plot_heatmap <- function(data, slider_cp, slider_cs) {
    data_copy <- data
    x <- data[,2:ncol(data)]
    samples <- colnames(x)
    group <- substr(samples, 1, 5) %>% factor()
    genes <- as.data.frame(data_copy[,1:2])
    y <- DGEList(counts=x,group=group, genes=genes)
    keep <- filterByExpr(y)
    y <- y[keep,,keep.lib.sizes=FALSE]
    y <- calcNormFactors(y)
    dgeList <- calcNormFactors(y, method="TMM") #Normalize your raw data using the TMM method
    dgeList <- cbind(dgeList$genes, dgeList$counts) %>% as.data.frame()
    rownames(dgeList) <- dgeList$Gene_id
    dgeList <- subset(dgeList, select = -c(1) )
    logcounts <- cpm(dgeList, log=TRUE) #Obtain the TMM values of your data
    var_genes <- apply(logcounts, 1, var)
    select_var <- names(sort(var_genes, decreasing=TRUE))[1:30]
    highly_variable_lcpm <- logcounts[select_var, ]
    col.pal <- RColorBrewer::brewer.pal(11, 'Set3')
    hm <- pheatmap(highly_variable_lcpm,trace="none", main="Top 30 variable genes across conditions",scale="row",margins=c(10,10))
    return(hm)
  }
  
  #' pca_results
  #' Return PCA results
  #' 
  #' PCA is performed over genes, and samples should be colored by diet.
  #' @param data tibble: a (n x _S_) data set
  #'
  #' @return dataframe of PCA results
  #'
  pca_results <- function(data) {
    data_copy <- data
    data_copy <- as.data.frame(data_copy)
    rownames(data) <- data$Gene_id
    data <- subset(data, select = -c(1) )
    t_data <- t(data) #transpose data so samples are on rows and genes are columns 
    pca <- prcomp(t_data)
    return(pca)
  }
  
  #' Perform and plot PCA using processed data.
  #' 
  #' PCA is performed over genes, and samples should be colored by Diet.
  #' Both `y` and `x` axis should have percent of explained variance included.
  #'
  #'
  #' @param data tibble: a (n x _S_) data set
  #' @param x user chosen principle component to plot 
  #' @param y user chosen principle component to plot 
  #'
  #' @return ggplot: scatter plot showing each sample in the first two PCs. 
  plot_pca <- function(data, x, y) {
    pca <- pca_results(data)
    percent_var <- pca$sdev^2 / sum(pca$sdev^2)
    
    plot_data <- pca$x %>% as.data.frame()
    group <- substr(rownames(plot_data), 1, 5) %>% factor()
    
    x_num = substr(x, 3, 3) %>% as.numeric()
    y_num = substr(y, 3, 3) %>% as.numeric()
    
    biplot <- ggplot2::ggplot(plot_data, ggplot2::aes(x=!!sym(x), y=!!sym(y), color=group)) + 
      ggplot2::geom_point() + 
      ggplot2:: xlab(str_glue("{x}  ", round(percent_var[x_num] * 100, 3),"% variance"))+ 
      ggplot2:: ylab(str_glue("{y}  ", round(percent_var[y_num] * 100, 3),"% variance"))+
      ggplot2::ggtitle("Principal Component Analysis Projections")
    
    return(biplot)
  }
  
  #' plot_fgseablot: Barplot of fgsea NES for top pathways selected by slider
  #'
  #' @param data tibble: adata set
  #' @param slider_padj user chosen padj threshold
  #'
  #' @return ggplot: barplot showing padj significant pathways that pass the threshold. 
  plot_fgseablot <- function(data, slider_padj){
    top_pos <- data %>% filter(padj < slider_padj) %>% pull(pathway)
    top_neg <- data %>% filter(padj < slider_padj) %>% pull(pathway)
    
    subset <- data %>% 
      filter(pathway %in% c(top_pos, top_neg)) %>%
      mutate(pathway = factor(pathway)) %>%
      mutate(plot_name = str_replace_all(pathway, '_', ' '))
    
    plot <- subset %>% 
      mutate(plot_name = forcats::fct_reorder(factor(plot_name), padj)) %>%
      ggplot() +
      geom_bar(aes(x=plot_name, y=NES, fill = NES > 0), stat='identity', show.legend = FALSE) +
      scale_fill_manual(values = c('TRUE' = 'red', 'FALSE' = 'blue')) + 
      theme_minimal(base_size = 8) +
      ggtitle('fgsea results for MSigDB gene sets') +
      ylab('Normalized Enrichment Score (NES)') +
      xlab('') +
      scale_x_discrete(labels = function(x) str_wrap(x, width = 80)) +
      coord_flip()
    return(plot)
  }
  
  #' plot_fgseasplot: Scatterplot of fgsea NES for top pathways selected by slider
  #' 
  #' @param data tibble: a data set
  #' @param slider_padj user chosen padj threshold
  #'
  #' @return ggplot: barplot showing padj significant pathways that pass the threshold. 
  plot_fgseasplot <- function (data, slider_padj){
    fgsea_splot <- data %>% 
      ggplot(aes(x=NES, y=-log10(padj), color=padj < slider_padj)) + 
      geom_point()+
      labs(color = str_glue('padj < {slider_padj}')) + 
      scale_color_manual(values=c("gray", 'black')) +
      ggtitle('Scatter plot of NES against -log10(padj)')
    return(fgsea_splot)
    }
  
  
  #' Draw and filter table
  #'
  #' @param dataf Data frame loaded by load_data()
  #' @param slider Negative number, typically from the slider input.
  #'
  #' @return Data frame filtered to p-adjusted values that are less than 
  #' 1 * 10^slider, columns for p-value and p-adjusted value have more digits 
  #' displayed.
  #' @details Same as above, this function is a standard R function. Tests will 
  #' evaluate it normally. Not only does this function filter the data frame to 
  #' rows that are above the slider magnitude, it should also change the format 
  #' of the p-value columns to display more digits. This is so that it looks 
  #' better when displayed on the web page. I would suggest the function 
  #' `formatC()`
  #'
  #' @examples draw_table(deseq_df, -210)
  #'    X  baseMean     log2FC     lfcSE      stat       pvalue         padj
  #'gene1 11690.780   9.852926 0.2644650  37.25607 8.45125e-304 1.54472e-299
  #'gene2  3550.435  -6.183714 0.1792708 -34.49369 9.97262e-261 9.11398e-257
  draw_table <- function(dataf, slider) {
    dataf$padj <- as.numeric(dataf$padj)
    dataf <- dplyr::filter(dataf, padj< (1*10^slider)) 
    dataf <- dataf %>% summarise(
      gene = Gene_id, 
      baseMean = formatC(baseMean),
      log2FC = formatC(log2FoldChange),
      lfcSE = formatC(lfcSE),
      stat = formatC(stat), 
      pvalue = formatC(pvalue),
      padj = formatC(padj)
    )
    return(dataf)
  }
  
  #' draw_countsTable: Draw and filter read counts table
  #'
  #' @param data Data frame loaded by load_countdata()
  #' @param slider_cp Slider to include genes with at least X percentile of variance
  #' @param slider_cs Slider to include genes with at least X samples that are non-zero
  #'
  #' @return Data frame summarizing filtered table after adjusted according to percentile 
  #' of variance and samples that are non-zeros
  draw_countsTable <- function(data, slider_cp, slider_cs){
    #get total number of genes and samples 
    d2 <- data
    gene_no <- nrow(d2)
    sample_no <- ncol(d2) - 1
    
    #exclude first two columns which are gene ids and symbols
    d1 <- d2[,2:ncol(d2)]#reactive value.. 
    #filter genes using variance 
    threshold <- quantile(apply(d1, 1, var), probs = c(slider_cp))
    new_data <- d1[(apply(d1, 1, var)) > threshold, ]
    
    # filter genes using samples using nonzeros
    new_data <- new_data[rowSums(new_data!=0)>slider_cs, ]
    #get number and percentage of genes passing and not passing after filtering
    num_of_genes <- nrow(new_data)
    per_of_genes <- ((nrow(new_data))/gene_no) * 100 
    diffgene <- gene_no - num_of_genes
    diffper <- 100-per_of_genes
    
    count_table <- tibble(TotalSamples=sample_no, TotalGenes=gene_no, Passing=paste(formatC(num_of_genes), ', ', formatC(per_of_genes), '%'), NotPassing=paste(formatC(diffgene), ', ', formatC(diffper), '%'))
    return(count_table)
  }
  
  #' draw_sampleInfo
  draw_sampleInfo<- function(data){
    get_coldata <- function(col_x) {
      # second find out if type is a factor or num 
      if (is.factor(col_x)){
        l <- levels(col_x)
        return(toString(l))
      }
      # otherwise it must be num, so check it for
      else {
        m <- mean(col_x)
        std <- sd(col_x)
        return(paste(formatC(m), " +/-",formatC(std)))
      }
    }
    metadata <- data
    metadata[sapply(metadata, is.character)] <- lapply(metadata[sapply(metadata, is.character)], as.factor)
    metadata$Tag <- factor(metadata$Tag)
    
    Type <- c("Factor", "Factor", "Factor", "Factor", "Factor", "Factor", "num", "num", "num")
    
    Information <- sapply(metadata, get_coldata)
    
    samples <- tibble(Columns = colnames(metadata))
    samples <- cbind(samples, Type)
    samples <- cbind(samples, Information)
    return(samples)
  }
  
  #draw_fgseatable
  draw_fgseatable <- function(data, choice, slider_padj){
    if (choice == "postive"){
      rs <- data %>% filter(padj < slider_padj) %>% slice_max(NES, n=10)
      return(rs)
    }
    if (choice == "negative"){
      rs <- data %>% filter(padj < slider_padj) %>% slice_min(NES, n=10)
      return(rs)
    }
  }
  
  #' These outputs aren't really functions, so they don't get a full skeleton, 
  #' but use the renderPlot() and renderTabel() functions to return() a plot 
  #' or table object, and those will be displayed in your application.
  output$volcano <- renderPlot({
    volcano_plot(load_deseqdata(), input$button1, input$button2, input$slider3, input$col1, input$col2)})
  
  #' These outputs aren't really functions, so they don't get a full skeleton, 
  #' but use the renderPlot() and renderTabel() functions to return() a plot 
  #' or table object, and those will be displayed in your application.
  output$histogram <- renderPlot({
    histo_plot(load_sampledata())
  })
  
  #' These outputs aren't really functions, so they don't get a full skeleton, 
  #' but use the renderPlot() and renderTabel() functions to return() a plot 
  #' or table object, and those will be displayed in your application.
  output$plotvarmean <- renderPlot({
    plot_variance_vs_mean(load_countdata(), input$slider1, input$slider2)
  })
  
  #' These outputs aren't really functions, so they don't get a full skeleton, 
  #' but use the renderPlot() and renderTabel() functions to return() a plot 
  #' or table object, and those will be displayed in your application.
  output$plothm <- renderPlot({
    plot_heatmap(load_countdata(), input$slider1, input$slider2)
  })
  
  #' These outputs aren't really functions, so they don't get a full skeleton, 
  #' but use the renderPlot() and renderTabel() functions to return() a plot 
  #' or table object, and those will be displayed in your application.
  output$pcaplot <- renderPlot({
    plot_pca(load_countdata(), input$dd1, input$dd2)
  })
  
  output$fgsea_barplot <- renderPlot({
    plot_fgseablot(load_fgseadata(), input$slider4)
  })
  
  output$fgsea_scatterplot <- renderPlot({
    plot_fgseasplot(load_fgseadata(), input$slider4)
  })
  
  #RenderTable Outputs here 
  # Same here, just return the table as you want to see it in the web page
  output$sample_tableInfo <-renderTable({
    draw_sampleInfo(load_sampledata())
  }) 
  
  # Same here, just return the table as you want to see it in the web page
  output$sampleTable <-DT::renderDataTable({load_sampledata()}) 
  
  # Same here, just return the table as you want to see it in the web page
  output$countsTable <-renderTable({
    draw_countsTable(load_countdata(), input$slider1, input$slider2)
  }) 
  
  # Same here, just return the table as you want to see it in the web page
  output$table <-renderTable({
    draw_table(load_deseqdata(), input$slider3)
  }) 
  
  output$de_summ_table <- DT::renderDataTable({load_deseqdata()})
  
  output$fgsea_table <- renderTable({
    draw_fgseatable(load_fgseadata(), input$NES_button, input$slider4)
  })
  
  output$download_fgsea <- downloadHandler(
    filename = function() {
      paste('filtered-data-', Sys.Date(), '.csv', sep = '')
    },
    content = function(file){
      write.csv(fgsea_filtered(),file)
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
