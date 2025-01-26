#!/usr/bin/Rscript

suppressMessages(library(rtracklayer))
suppressMessages(library(karyoploteR))
suppressMessages(library(GenomicFeatures))

plotCoverageKaryoploteR <- function(zoom_region, chromosome, output_file, custom.genome, txdb) {
   
  #---------------------------------------------------------------------
  # 1) Define replicate-specific parameters for 0.025 h^-1 replicates
  #    (R1, R2, R3) in Data Panel 1
  #---------------------------------------------------------------------
  replicates_025 <- data.frame(
    replicate = c("R1", "R2", "R3"),
    # Overall coverage BAM files (used only to get max coverage)
    overall = c(
      "data/retentostat/mapped/SRR25915888Aligned.sortedByCoord.out.bam",  # R1
      "data/retentostat/mapped/SRR25915870Aligned.sortedByCoord.out.bam",  # R2
      "data/retentostat/mapped/SRR25915879Aligned.sortedByCoord.out.bam"   # R3
    ),
    # Split (forward) coverage BAM files
    plus = c(
      "data/retentostat/mapped_split_by_strand/SRR25915888_fwd.bam",
      "data/retentostat/mapped_split_by_strand/SRR25915870_fwd.bam",
      "data/retentostat/mapped_split_by_strand/SRR25915879_fwd.bam"
    ),
    # Split (reverse) coverage BAM files
    minus = c(
      "data/retentostat/mapped_split_by_strand/SRR25915888_rev.bam",
      "data/retentostat/mapped_split_by_strand/SRR25915870_rev.bam",
      "data/retentostat/mapped_split_by_strand/SRR25915879_rev.bam"
    ),
    # Track positions for the plus and minus tracks
    plus_r0  = c(0.375, 0.625, 0.875),
    plus_r1  = c(0.5, 0.75, 1.0),
    minus_r0 = c(0.375, 0.625, 0.875),
    minus_r1 = c(0.25, 0.5, 0.75),
    stringsAsFactors = FALSE
  )
  
  #---------------------------------------------------------------------
  # 2) Define replicate-specific parameters for R10 replicates
  #    (R1, R2, R3) in Data Panel 2
  #---------------------------------------------------------------------
  replicates_R10 <- data.frame(
    replicate = c("R1", "R2", "R3"),
    overall = c(
      "data/retentostat/mapped/SRR25915871Aligned.sortedByCoord.out.bam",  # R1
      "data/retentostat/mapped/SRR25915867Aligned.sortedByCoord.out.bam",  # R2
      "data/retentostat/mapped/SRR25915875Aligned.sortedByCoord.out.bam"   # R3
    ),
    plus = c(
      "data/retentostat/mapped_split_by_strand/SRR25915871_fwd.bam",
      "data/retentostat/mapped_split_by_strand/SRR25915867_fwd.bam",
      "data/retentostat/mapped_split_by_strand/SRR25915875_fwd.bam"
    ),
    minus = c(
      "data/retentostat/mapped_split_by_strand/SRR25915871_rev.bam",
      "data/retentostat/mapped_split_by_strand/SRR25915867_rev.bam",
      "data/retentostat/mapped_split_by_strand/SRR25915875_rev.bam"
    ),
    # Note that for some tracks r0 > r1, matching your original code
    plus_r0  = c(0.15, 0.45, 0.75),
    plus_r1  = c(0.0,  0.3,  0.6),
    minus_r0 = c(0.15, 0.45, 0.75),
    minus_r1 = c(0.3,  0.6,  0.9),
    stringsAsFactors = FALSE
  )
  
  #---------------------------------------------------------------------
  # 3) Open the PNG device
  #---------------------------------------------------------------------
  png(filename = output_file, width = 465, height = 225, units = 'mm', res = 300)
  
  #---------------------------------------------------------------------
  # 4) Set up karyoploteR parameters and create the main plot
  #---------------------------------------------------------------------
  pp <- getDefaultPlotParams(plot.type = 2)
  pp$leftmargin     <- 0.15
  pp$data1height    <- 300
  pp$data2height    <- 250
  pp$ideogramheight <- 10
  
  kp <- plotKaryotype(
    genome       = custom.genome,
    plot.params  = pp,
    plot.type    = 2,
    chromosomes  = chromosome,  # adjust if needed
    cex          = 1,
    zoom         = zoom_region
  )
  
  # Add base numbers (tick marks)
  kpAddBaseNumbers(kp, tick.dist = 1000, add.units = TRUE, digits=3)
  
  # Backgrounds on data panels
  kpDataBackground(kp, data.panel = 1, r0 = 0, r1 = 1, color = "#ffffff")
  kpDataBackground(kp, data.panel = 2, r0 = 0, r1 = 1, color = "#ffffff")
  
  #---------------------------------------------------------------------
  # 5) Plot Genes
  #---------------------------------------------------------------------
  kpPlotGenes(
    kp,
    data                       = txdb,
    plot.transcripts           = TRUE,
    r0                         = 0,
    r1                         = 0.25,
    data.panel                 = 1,
    col                        = "#46dc37",
    marks.col                  = "black",
    plot.transcripts.structure = FALSE,
    add.strand.marks           = TRUE,
    non.coding.exons.col       = "#dc4837",
    coding.exons.col           = "#46dc37",
    gene.name.position         = "top",
    mark.height                = 0.5,
    gene.name.col              = "black",
    transcript.margin          = 1.1
  )
  
  #---------------------------------------------------------------------
  # 6) Helper function to plot one replicate's coverage
  #---------------------------------------------------------------------
  plotOneRep <- function(kp, overall_bam, plus_bam, minus_bam,
                         plus_r0, plus_r1, minus_r0, minus_r1,
                         plus_label, minus_label,
                         data.panel.plus = 1,
                         data.panel.minus = 1) {
    
    # 1) Plot overall coverage (invisible track just to get max coverage)
    kp <- kpPlotBAMCoverage(
      kp, data = overall_bam,
      ymin = 0, r0 = 0, r1 = 0, col = "#0e87eb", border = NA
    )
    ymax <- kp$latest.plot$computed.values$max.coverage
    
    # 2) Plot the plus strand coverage
    kp <- kpPlotBAMCoverage(
      kp, data = plus_bam,
      ymin = 0, r0 = plus_r0, r1 = plus_r1,
      col = "#0e87eb", border = NA,
      data.panel = data.panel.plus,
      ymax = ymax
    )
    kpAxis(kp, r0 = plus_r0, r1 = plus_r1, data.panel = data.panel.plus)
    kpAddLabels(kp, labels = plus_label, pos = 1,
                r0 = plus_r0, r1 = plus_r1, label.margin = 0.0725,
                data.panel = data.panel.plus)
    
    # 3) Plot the minus strand coverage
    kp <- kpPlotBAMCoverage(
      kp, data = minus_bam,
      ymin = 0, r0 = minus_r0, r1 = minus_r1,
      col = "#dc4837", border = NA,
      data.panel = data.panel.minus,
      ymax = ymax
    )
    kpAxis(kp, r0 = minus_r0, r1 = minus_r1, data.panel = data.panel.minus)
    kpAddLabels(kp, labels = minus_label, pos = 1,
                r0 = minus_r0, r1 = minus_r1, label.margin = 0.0725,
                data.panel = data.panel.minus)
    
    return(kp)
  }
  
  #---------------------------------------------------------------------
  # 7) Loop over 0.025 h^-1 replicates, plotting in data.panel=1
  #---------------------------------------------------------------------
  for(i in seq_len(nrow(replicates_025))) {
    rinfo <- replicates_025[i, ]
    kp <- plotOneRep(
      kp,
      overall_bam = rinfo$overall,
      plus_bam    = rinfo$plus,
      minus_bam   = rinfo$minus,
      plus_r0     = rinfo$plus_r0,
      plus_r1     = rinfo$plus_r1 - 0.025,
      minus_r0    = rinfo$minus_r0,
      minus_r1    = rinfo$minus_r1 + 0.025,
      plus_label  = paste0("0.025 h^-1 ", rinfo$replicate, " +ve"),
      minus_label = paste0("0.025 h^-1 ", rinfo$replicate, " -ve"),
      data.panel.plus = 1,
      data.panel.minus = 1
    )
  }
  
  #---------------------------------------------------------------------
  # 8) Loop over R10 replicates, splitted coverage in data.panel=2
  #---------------------------------------------------------------------
  for(i in seq_len(nrow(replicates_R10))) {
    rinfo <- replicates_R10[i, ]
    kp <- plotOneRep(
      kp,
      overall_bam = rinfo$overall,
      plus_bam    = rinfo$plus,
      minus_bam   = rinfo$minus,
      plus_r0     = rinfo$plus_r0,
      plus_r1     = rinfo$plus_r1 + 0.03,
      minus_r0    = rinfo$minus_r0,
      minus_r1    = rinfo$minus_r1 - 0.03,
      plus_label  = paste0("R10 ", rinfo$replicate, " +ve"),
      minus_label = paste0("R10 ", rinfo$replicate, " -ve"),
      data.panel.plus = 2,
      data.panel.minus = 2
    )
  }
  
  #---------------------------------------------------------------------
  # 9) Close the device
  #---------------------------------------------------------------------
  dev.off()
}



################################################################
################################################################
################################################################

chr_sizes <- read.delim("reference_genome/cbs7435chromSizes.txt",header=FALSE)
colnames(chr_sizes) <- c("chr","end")
chr_sizes$start <- 1
chr_sizes <- chr_sizes[,c("chr","start","end")]
custom.genome <- toGRanges(chr_sizes)

txdb <- makeTxDbFromGFF("reference_genome/cbs7435_updated.gtf", format = "gtf")

################################################################
################################################################

# MSTRG.830 
mstrg.830.region <- toGRanges(data.frame("cbs7435_chr1",  1383000, 1390000))

plotCoverageKaryoploteR(
  zoom_region  = mstrg.830.region,
  chromosome   = "cbs7435_chr1",
  output_file  = "lncrna_annotation/plots/MSTRG.830_Scaled_All_Densities.png",
  custom.genome = custom.genome,
  txdb         = txdb
)

################################################################

# MSTRG.1680
mstrg.1680.region <- toGRanges(data.frame("cbs7435_chr1",  2826000, 2832000))

plotCoverageKaryoploteR(
  zoom_region  = mstrg.1680.region,
  chromosome   = "cbs7435_chr1",
  output_file  = "lncrna_annotation/plots/MSTRG.1680_Scaled_All_Densities.png",
  custom.genome = custom.genome,
  txdb         = txdb
)


################################################################

# MSTRG.3443 

mstrg.3443.region <- toGRanges(data.frame("cbs7435_chr3",  598000,  607000))

plotCoverageKaryoploteR(
  zoom_region  = mstrg.3443.region,
  chromosome   = "cbs7435_chr3",
  output_file  = "lncrna_annotation/plots/MSTRG.3443_Scaled_All_Densities.png",
  custom.genome = custom.genome,
  txdb         = txdb
)
