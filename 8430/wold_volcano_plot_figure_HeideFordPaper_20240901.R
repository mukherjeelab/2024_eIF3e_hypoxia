#' ---
#' title: "Analysis of K.Webb's thermal profiling analysis for Heide Ford Paper 9/1/2024"
#' author: "WMO"
#' date: "Sept 1 2024"
#' ---
#  
#  

library(tidyverse)
library(limma)
library(lattice)
library(grid)
library(gridExtra)
library(EnhancedVolcano)
library(magrittr)
# define functions
# 
source('./fdrfunctions.R')

stats_df <- readr::read_csv(file = "df_wide.csv", col_names = T)

stats_df <- stats_df %>%
  dplyr::relocate(Gene.names, ID,
                  D_vs_C_CI.L, D_vs_C_CI.R, D_vs_C_diff, D_vs_C_p.val, D_vs_C_p.adj, D_vs_C_significant,
                  Protein.names, .after = "name")



####################

genes2lab <- stats_df$name[stats_df$D_vs_C_p.adj < 0.01]

cpd8430_volcano <- EnhancedVolcano(stats_df,
                                      lab = stats_df$name,
                                     # selectLab = selected_sites,
                                      x = 'D_vs_C_diff',
                                      y = 'D_vs_C_p.val',
                                      title = NULL,
                                      subtitle = NULL,
                                      caption = NULL,
                                      selectLab = genes2lab,
                                      ylim = c(0, 12),
                                      xlim = c(-1,1)*0.8,
                                      xlab = bquote(~Log[2]~ '(Cmpd8430/DMSO)'),
                                      pCutoff = 0.01,
                                      pCutoffCol = 'D_vs_C_p.adj',
                                      FCcutoff = 0.2,
                                      pointSize = 3.0,
                                      labSize = 4.0,
                                      boxedLabels = TRUE,
                                      colAlpha = 4/5,
                                      gridlines.major = F,
                                      gridlines.minor = FALSE,
                                      legendPosition = 'right',
                                      legendLabSize = 10,
                                      legendIconSize = 4.0,
                                      drawConnectors = TRUE,
                                      widthConnectors = 1.5,
                                      lengthConnectors = unit(0.5, "npc"),
                                      max.overlaps = 20,
                                      min.segment.length = 0,
                                      arrowheads = F,
                                      #  col=c('black', 'black', 'black', 'red3'),
                                      legendLabels=c('Not sig.','|Log2FC| > 0.2','q-value < 0.01',
                                                     '|Log2FC| > 0.2 & q-value < 0.01'))

windows(width = 20, height = 10)
cpd8430_volcano


ggsave(cpd8430_volcano, 
       width = 20, height = 10, units = "in",
       device = "pdf", 
       filename = paste0('Cmpd8430_HfordPaper_volcanoplot_',
                         Sys.Date(), '.pdf'))




