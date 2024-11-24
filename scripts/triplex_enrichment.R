#!/usr/bin/env Rscript
## -----------------------------------------------------------------------------
##
## Script name: triplex_enrichment.R
##
## Purpose: This script performs Gene Ontology (GO) and KEGG enrichment analysis on lncRNA interacting pairs and outputs the results into an Excel workbook.
##
## Inputs:
##   1. Interacting pairs file (text file with "Source" and "Target" columns)
##   2. Output file name for the Excel workbook
##
##
## -----------------------------------------------------------------------------
##
## Description: 
##   - Reads interacting lncRNA-target gene pairs and identifies unique lncRNAs (sources).
##   - For each unique lncRNA, performs GO and KEGG enrichment analysis on its target genes.
##   - Writes results to an Excel file, with a sheet dedicated to each lncRNA, including separate sections for GO and KEGG results.
##
## -----------------------------------------------------------------------------


suppressMessages(library("gprofiler2"))
suppressMessages(library("openxlsx"))

args <- commandArgs(TRUE)
interacting_pairs <- args[1]
out_file <- args[2]

triplex_pairs <- read.table(interacting_pairs)
colnames(triplex_pairs) <- c("Source", "Target")
uniq_triplex_source <- unique(triplex_pairs$Source)


# Create a new Excel workbook
output_wb <- createWorkbook()

# Create a new sheet with the cluster label as the sheet name
addWorksheet(output_wb, sheetName = "lncRNAs")

row_counter = 1


for (i in uniq_triplex_source) {
	pairs <- triplex_pairs[triplex_pairs$Source == i,"Target"]

 # Perform gene ontology analysis using gprofiler2's gost function for GO enrichment
# upload_GMT_file(gmtfile="functional_predictions/GO_analysis/gprofiler_converted.gmt")
	GO_name <- "gp__Uku7_mqlG_9fI"  
	suppressMessages(GO_result <- gost(pairs, organism = GO_name))
  
  # Perform gene ontology analysis using gprofiler2's gost function for the kegg ontology
# upload_GMT_file(gmtfile="functional_predictions/GO_analysis/Kegg_Ontology.gmt")
	kegg_name <- "gp__PEAO_RSQ8_tq0" 
	suppressMessages(kegg_result <- gost(pairs, organism = kegg_name))

 # Check if both results are empty, if not, fill in the data frame
	if (length(GO_result$result) > 0) {
		GO_rows <- nrow(GO_result$result)
		} else { GO_rows <- 0 }

  
	if (length(kegg_result$result) > 0) {
		kegg_rows <- nrow(kegg_result$result)
		} else { kegg_rows <- 0 }
		
	
	num_rows = max(GO_rows, kegg_rows)

	if (num_rows != 0)  {

	writeData(output_wb, sheet = "lncRNAs", i, 
              startCol = 1, startRow = row_counter)

	writeData(output_wb, sheet = "lncRNAs", "GO enrichment", 
              startCol = 3, startRow = row_counter)

  # Write the GO_enrichment data to the sheet
	writeData(output_wb, sheet = "lncRNAs", GO_result$result, startCol = 3, startRow = row_counter + 1)

	writeData(output_wb, sheet = "lncRNAs", "Kegg enrichment", 
              startCol = 18, startRow = row_counter)

  # Write the kegg_enrichment data after the empty rows
	writeData(output_wb, sheet = "lncRNAs", kegg_result$result, 
            startCol = 18, startRow = row_counter + 1)

	row_counter = row_counter + num_rows + 3}
  
}


# Save the workbook to the specified file
saveWorkbook(output_wb, file =out_file, overwrite=TRUE)

