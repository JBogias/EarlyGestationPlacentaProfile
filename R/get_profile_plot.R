#' @title Visualise gene profiles
#' 
#' @description Visualise the gene expression, transcript expression, and transcript proportions of a gene. This function also visualises the intron-exon structures of every transcript within the gene as well.
#' 
#' @return A four panel figure with the gene and transcript expression profiles along with the transcript proportions and intron-exon structures.
#' 
#' 


require(edgeR)
require(ggplot2)
require(tidyverse)
require(tximport)
require(EnsDb.Hsapiens.v75)
require(cowplot)
require(here)

source(here("R/plot_expression.R"))
source(here("R/plot_proportions.R"))
source(here("R/plot_structures.R"))

txprofiler <- function(dgelist,
                       dtelist,
                       target_gene,
                       variable,
                       transcript_variable = NULL,
                       dtu_variable = NULL,
                       abundance_path = NULL,
                       tpm_method = "none",
                       gene_x_label = NULL,
                       transcript_x_label = NULL,
                       dtu_x_label = NULL) {
 
   gene_plot <- plot_txp(
    dgelist = dgelist,
    target_gene = target_gene,
    variable = variable,
    level = "gene",
    path = NULL,
    keep_legend = FALSE
  ) +
    labs(x = if_else(is.null(gene_x_label), variable, gene_x_label)) +
    ggtitle(NULL)
  
  transcript_plot <- plot_txp(
    dgelist = dtelist,
    target_gene = target_gene,
    variable = if_else(
      is.null(transcript_variable), variable, transcript_variable
    ),
    level = "transcript",
    tpm_method = tpm_method,
    path = abundance_path,
    keep_legend = FALSE
  ) +
    labs(x = if_else(
      is.null(transcript_x_label), variable, transcript_x_label)
      ) +
    ggtitle(NULL)
  
  proportion_plot <- plot_proportions(
    dtelist = dtelist,
    target_gene = target_gene,
    variable = if_else(is.null(dtu_variable), variable, dtu_variable),
    tpm_method = tpm_method,
    path = abundance_path,
    keep_legend = FALSE
  ) +
    labs(x = if_else(is.null(dtu_x_label), variable, dtu_x_label)) +
    ggtitle(NULL)
  
  structure_plot <- plot_structures(
    dtelist = dtelist,
    target_gene = target_gene,
    annotation = EnsDb.Hsapiens.v75
  )
  
  left_col <- cowplot::plot_grid(structure_plot,
                                 proportion_plot,
                                 labels = c("A", "C"),
                                 nrow = 2)
  
  right_col <- cowplot::plot_grid(gene_plot,
                                  transcript_plot,
                                  labels = c("B", "D"),
                                  nrow = 2)
  
  legend_plot <- cowplot::get_legend(
    plot_txp(
      dgelist = dtelist,
      target_gene = target_gene,
      variable = if_else(
        is.null(transcript_variable), variable, transcript_variable
      ),
      level = "transcript",
      path = abundance_path,
      keep_legend = TRUE
    )
  )
  
  profile_plot <- cowplot::plot_grid(
    left_col,
    right_col,
    legend_plot,
    ncol = 3,
    rel_widths = c(1, 1, 0.4)
  )
  
  title_theme <- cowplot::ggdraw() +
    draw_label(
      paste0(target_gene, " profile"),
      colour = "black",
      fontface = "bold",
      size = 30
    )
  
  cowplot::plot_grid(title_theme,
                     profile_plot,
                     nrow = 2,
                     rel_heights = c(0.05, 1))
}



