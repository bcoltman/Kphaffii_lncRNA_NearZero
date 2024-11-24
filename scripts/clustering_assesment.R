
#!/usr/bin/Rscript
## -----------------------------------------------------------------------------
##
## Script name: clustering_assesment.r
##
## Purpose of script: To assess and analyse the heirarchical clustering results 
## from the MPLN clustering package.
##
##

# input arguments
args <- commandArgs(TRUE)
results_dir <- args[1] 
gtf <- args[2]
lncrna_list <- args[3]

#results_dir <- "functional_predictions/clustering/mpln"
#gtf <- "transcriptome_assembly/stringtie.all.transcripts.gtf"
#lncrna_list <- "lncrna_annotation/monoexonic_filter/first_lncrna_list.txt"

suppressMessages(library(GenomicFeatures))
library("MPLNClust")
library("dplyr")
library("gprofiler2")
library("openxlsx")
library("ggplot2")
library(RColorBrewer)

VisualizeLine <- function(dataset,
                              clusterMembershipVector = NA,
                              fileName = paste0('Plot_', Sys.Date()),
                              output_dir = "") {

  if(max(clusterMembershipVector) > 17) {
    qualColPals <- RColorBrewer::brewer.pal.info[brewer.pal.info$category == 'qual', ]
    coloursBarPlot <- unlist(mapply(RColorBrewer::brewer.pal,
                                    qualColPals$maxcolors,
                                    rownames(qualColPals)))
  } else {
    coloursBarPlot <- c('#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6',
                        '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324',
                        '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1',
                        '#000075', '#808080')
  }


  # Log transform the dataset
  dataset <- log(dataset + 1)

  # Combine log-transformed dataset with cluster membership vector
  DataPlusLabs <- cbind(dataset, cluster=clusterMembershipVector)

  # Convert DataPlusLabs to a dataframe
  df <- as.data.frame(DataPlusLabs)
  df$genes <- rownames(df)
  
  df_long <- reshape2::melt(df, id.vars = c("genes","cluster"))
  cluster_means <- df_long %>% group_by(variable, cluster) %>% summarise(mean_value = mean(value, na.rm = TRUE))
  cluster_means <- as.data.frame(cluster_means)


  # Plot using ggplot2 with facet_grid

  png(file = paste0(results_dir, "/mpln/",fileName, ".png"),width=800, height=600)
  p <- ggplot(data = df_long, aes(x = variable, y = value, color = as.factor(cluster), group = genes)) +
  geom_line() +
  geom_line(data = df_long %>% group_by(variable, cluster) %>% summarise(mean_value = mean(value, na.rm = TRUE)),
            aes(group = cluster, x = variable, y = mean_value),color = "black") +
  scale_color_manual(values = coloursBarPlot) +
  labs(x = "Samples", y = "Expression (log counts)")+
  theme_minimal() + facet_wrap(~ cluster, ncol = 5) + 
  theme(legend.position="none",axis.text = element_text(size = 20),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),strip.text = element_text(size=20),
        axis.title=element_text(size=20,face="bold"))
   
  print(p)
  dev.off()
}



load(file=paste0(results_dir, "/mpln/sig_mplnResults.RData"))


# lncrna transcript list
lncrna_transcripts <- read.table(lncrna_list,  
                                 header = F, stringsAsFactors = F)$V1

#lncrna_transcripts_conservative <- read.table(
 #                                header = F, stringsAsFactors = F)$V1

# load gene GTF
stringtie_txdb <- makeTxDbFromGFF(gtf)

# collapse lncRNA transcripts to genes 
# will return a warning: 'select()' returned 1:1 mapping between keys and columns. Therefore, suppress


suppressMessages(lncrna_genes_all <- AnnotationDbi::select(stringtie_txdb,
                                            keys = lncrna_transcripts,
                                            columns = c("TXCHROM", 
                                                        "GENEID", 
                                                        "TXSTRAND", 
                                                        "TXSTART", 
                                                        "TXEND"), 
                                            keytype = "TXNAME"))

# lncrna genes
lncrna_genes <- lncrna_genes_all[!is.na(lncrna_genes_all$GENEID), ]$GENEID


opt_clusters <- 4


clusters <- data.frame(ID=rownames(mplnResults$dataset),Cluster=mplnResults$allResults[[opt_clusters]]$clusterlabels)

# Initialize an empty data frame to store the processed data
processed_clusters <- data.frame()

# Process each row in the input data
for (i in 1:nrow(clusters)) {
  # Split the first column by "|"
  ids <- strsplit(as.character(clusters[i, 1]), "\\|")[[1]]
  n <- length(ids)
  
  if (n == 1) {
    # If only one element, keep it as is
    processed_clusters <- rbind(processed_clusters, clusters[i, ])
  } else {
    # If multiple elements, create multiple rows with the same value in the second column
    for (j in 2:n) {
      new_row <- data.frame(ids[j], clusters[i, 2])
      colnames(new_row) <- colnames(clusters)
      processed_clusters <- rbind(processed_clusters, new_row)
    }
  }
}





# Assuming your data is in a data frame named processed_data
# Column 1: Gene IDs
# Column 2: Cluster Labels

##########
# Get unique cluster labels
unique_clusters <- sort(unique(processed_clusters[, 2]))

# Create an Excel workbook to store the data
output_wb <- createWorkbook()

# Loop through unique cluster labels
for (cluster_id in unique_clusters) {

  cluster_label <- paste("Cluster", cluster_id)  # Create a unique sheet name
  
  # Filter the data for the current cluster label
  cluster_data <- subset(processed_clusters, processed_clusters[, 2] == cluster_id)
  
  # Extract gene IDs as a vector
  cluster_genes <- cluster_data[, 1]
  
  # Perform gene ontology analysis using gprofiler2's gost function for GO enrichment
# upload_GMT_file(gmtfile="functional_predictions/GO_analysis/gprofiler_converted.gmt")
  GO_name <- "gp__Uku7_mqlG_9fI"  
  GO_result <- gost(cluster_genes, organism = GO_name)
  
  # Perform gene ontology analysis using gprofiler2's gost function for the kegg ontology
# upload_GMT_file(gmtfile="functional_predictions/GO_analysis/Kegg_Ontology.gmt")
  kegg_name <- "gp__PEAO_RSQ8_tq0" 
  kegg_result <- gost(cluster_genes, organism = kegg_name)

  # Create a list to store the results for this cluster
  cluster_result <- list(
    cluster_label = cluster_label,
    genes = cluster_genes,
    GO_enrichment = NA,
    kegg_enrichment = NA
  )
 # Check if both results are empty, if not, fill in the data frame
  if (length(GO_result$result) > 0) {
    cluster_result$GO_enrichment <- GO_result$result
    num_empty_rows <- nrow(cluster_result$GO_enrichment)
  } else { num_empty_rows <- 0 }

  
  if (length(kegg_result$result) > 0) {
    cluster_result$kegg_enrichment <- kegg_result$result
  }

# Create a new sheet with the cluster label as the sheet name
  addWorksheet(output_wb, sheetName = cluster_label)

# Write the genes to the first column
  writeData(output_wb, sheet = cluster_label, cluster_result$genes, startCol = 1, startRow = 1)

  writeData(output_wb, sheet = cluster_label, "GO enrichment", 
              startCol = 5, startRow = 1)

  # Write the GO_enrichment data to the sheet
  writeData(output_wb, sheet = cluster_label, cluster_result$GO_enrichment, startCol = 5, startRow = 2)

  # Add two empty rows after the GO_enrichment data
  for (i in 1:2) {
    writeData(output_wb, sheet = cluster_label, "", 
              startCol = 5, startRow = num_empty_rows + i + 2)
  }

  writeData(output_wb, sheet = cluster_label, "Kegg enrichment", 
              startCol = 5, startRow = num_empty_rows + 5)

  # Write the kegg_enrichment data after the empty rows
  writeData(output_wb, sheet = cluster_label, cluster_result$kegg_enrichment, 
            startCol = 5, startRow = num_empty_rows + 6)


}

############################ NEED TO ADD A SHEET WITH THE PROBABilitites of each lnCRNA for each cluster

optimal <- data.frame(mplnResults$allResults[[opt_clusters]]$clusterlabels)
colnames(optimal) <- "Optimal"
rownames(optimal) <- rownames(mplnResults$dataset)

membership <- as.data.frame(mplnResults$allResults[[opt_clusters]]$probaPost)
colnames(membership) <- paste("Cluster", 1:ncol(membership), sep = "")
rownames(membership) <- rownames(mplnResults$dataset)

membership <- cbind(optimal, membership)

addWorksheet(output_wb, sheetName = "lncRNA membership")
#membership <- cbind(cluster_label=cl$cluster, cl$membership)
lncrna_membership <- membership[rownames(membership) %in% unique(lncrna_genes),]
writeData(output_wb, sheet = "lncRNA membership", lncrna_membership, 
              startCol = 1, startRow = 1, rowNames=TRUE)


# Save the workbook to the specified file
saveWorkbook(output_wb, file =paste0(results_dir, "/mpln/mpln_cluster_data.xlsx"), overwrite=TRUE)


###############
####Plots

#for (i in 2:20) {
#    MPLNVisuals <- VisualizeLine(dataset = mplnResults$dataset,
#	clusterMembershipVector =mplnResults$allResults[[i]]$clusterlabels,
#	fileName = paste0('/Line',i),
#	output_dir=results_dir)
#  }

MPLNVisuals <- VisualizeLine(dataset = mplnResults$dataset,
	clusterMembershipVector =mplnResults$allResults[[opt_clusters]]$clusterlabels,
	fileName = '/Line',
	output_dir=paste0(results_dir, "/mpln/"))

setwd(paste0(results_dir,"/mpln/"))

mplnViz <- mplnVisualizeBar(vectorObservations=1:nrow(mplnResults$dataset), 
			probabilities=mplnResults$allResults[[opt_clusters]]$probaPost,
			clusterMembershipVector=mplnResults$allResults[[opt_clusters]]$clusterlabels,
			format = 'png',fileName = 'PosteriorProbabilities',
			printPlot=TRUE)

mplnViz <- mplnVisualizeHeatmap(dataset=mplnResults$dataset, 
			clusterMembershipVector=mplnResults$allResults[[opt_clusters]]$clusterlabels,
			printPlot=TRUE,fileName = 'Heatmap',
			format = 'png')

png(file="Alluvial_AlternativeClustering.png", width=800, height=600)
MPLNVisualsAlluvial <- mplnVisualizeAlluvial(nObservations = nrow(mplnResults$dataset),
			firstGrouping = mplnResults$BICresults$BICmodelSelectedLabels,
			secondGrouping = mplnResults$ICLresults$ICLmodelSelectedLabels,
			thirdGrouping =  mplnResults$AICresults$AICmodelSelectedLabels,
			fourthGrouping = mplnResults$AIC3results$AIC3modelSelectedLabels, 
			fileName = "Alluvial",
			printPlot = FALSE) 
dev.off()

