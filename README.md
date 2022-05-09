# bf591-finalproject

This RShiny app attempts to summarize RNASeq bioinformatics studies and processes into multiple tabs. This app only works specifically for data generated by Hérault et al., 2019 study. 

The app can be broken into four parts: 
1. **Samples** - This tab allows users to explore metadata through a series of functionalities such as viewing the summary table, entire sortable data-table and density plot of sample groups. 
2. **Counts** - After uploading raw read counts data here, users will have access to a plethora of plots (diagnostic scatter plot, clustered heatmap, PCA) that they can interact with by adjusting percentile of gene count variance and number of nonzero gene counts. 
3. **DeSeq** - Differential expression identifies which genes, if any, are implicated in a specific biological comparison. This component allows the user to load and explore a differential expression dataset.
4. **FGSEA** - A table of fgsea results from the differential expression data is displayed in this tab, along with a barplot of fgsea NES for top pathways selected by slider, scatter plot of NES on x-axis and -log10 adjusted p-value on y-axis, with gene sets below threshold in grey color, and a download-able current filtered table results.

Things I still need to work on: 
- Heatmap shows only top highly variable 30 gene counts and does not interact with user input 
- Error display in samples and Counts tabe when first loading up the app 
- Making input choices more dynamic and not showing up before user has uploaded their files!
- Column data types are manually entered and future work should include automating this process 

**Reference**: Hérault, F., Houée-Bigot, M., Baéza, E., Bouchez, O., Esquerré, D., Klopp, C. and Diot, C., 2019. RNA-seq analysis of hepatic gene expression of common Pekin, Muscovy, mule and hinny ducks fed ad libitum or overfed. BMC genomics, 20(1), pp.1-14.
