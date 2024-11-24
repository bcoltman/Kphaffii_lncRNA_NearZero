#!/usr/bin/Rscript
## -----------------------------------------------------------------------------
##
## Script name: simplify_lncrna_class.R
##
## Purpose: Simplify FEELnc lncRNA classifications for easier analysis
##
##
## -----------------------------------------------------------------------------
##
## Description: This script simplifies and organizes lncRNA classifications produced by FEELnc.
##
## -----------------------------------------------------------------------------

args <- commandArgs(TRUE)
class_file <- args[1]
out_class_file <- args[2]

feelnc_class <- read.table(class_file,
                           header = T,
                           stringsAsFactors = F,
                           sep = "\t")

feelnc_class <- feelnc_class[feelnc_class$isBest == 1, ] # kept only the best

mod_class <- feelnc_class[, c(1:5, 11,8)]
position <- feelnc_class$location
position[position == "exonic"] <- "overlapping"
position[position == "intronic"] <- "overlapping"
mod_class$strand <- feelnc_class$direction
mod_class$position <- position

summary_class <- matrix(, dim(mod_class)[1])
for (i in 1:dim(mod_class)[1]) {
  if ((mod_class$strand[i] == "antisense") &&
      (mod_class$position[i] == "upstream") &&
      (mod_class$distance[i] <= 2000)) {
    summary_class[i] <- "divergent"

  } else if ((mod_class$strand[i] == "antisense") &&
             (mod_class$position[i] == "upstream") &&
             (mod_class$distance[i] > 2000)) {
    summary_class[i] <- "upstream antisense"

  } else if ((mod_class$strand[i] == "antisense") &&
             (mod_class$position[i] == "downstream")) {
    summary_class[i] <- "downstream antisense"

  } else if ((mod_class$strand[i] == "antisense") &&
             (mod_class$position[i] == "overlapping")) {
    summary_class[i] <- "antisense"

  } else if ((mod_class$strand[i] == "sense") &&
             (mod_class$position[i] == "overlapping")) {
    summary_class[i] <- "overlap sense"

  } else if ((mod_class$strand[i] == "sense") &&
             (mod_class$position[i] == "downstream")) {
    summary_class[i] <- "downstream sense"
  } else if ((mod_class$strand[i] == "sense") &&
             (mod_class$position[i] == "upstream")) {
    summary_class[i] <- "upstream sense"
  }
}

mod_class$summary_class <- summary_class
write.csv(mod_class, file = out_class_file)
