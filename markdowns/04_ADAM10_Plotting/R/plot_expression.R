#' @title Visualise expression of individual transcripts
#'
#' @description Transcript expression will be visualised from a kallisto or salmon output. Gene-level expression may be incorporated to through a gene counts file typically sourced from featureCounts
#'
#' @return Expression plotted at the individual transcript level, with gene-level expression if desired
#'
#' @param dgelist A DGEList -or SummarizedExperiment- object containing transcript counts, gene/txp information, and metadata
#' @param target_gene Character vector of gene of interest


# Add functionality to work with summarized experiments
# Need to adjust how the text size is edited, put all of the changeable
# elements in each plot individually, instead of in the cowplot

tpm <- function(dgelist) {
    rpk <- dgelist$counts/(dgelist$genes$length/1000)
    scaling <- colSums(rpk)/1e6
    tpm <- rpk/scaling
    tpm
}



plot_txp <- function(
    dgelist,
    target_gene,
    variable,
    level = "transcript",
    tpm_method = "standard",
    path = NULL,
    keep_legend = TRUE
) {
  
  require(edgeR)
  require(ggplot2)
  require(magrittr)
  require(dplyr)
  require(readr)
  require(stringr)
  require(tidyr)
  require(tibble)
  require(magrittr)
  require(tximport)
  require(here)
  
  if (level == "transcript") {
    if (tpm_method == "tximport") {
      
      abundance_salmon <- read_csv(path)
      
      tpm <- abundance_salmon %>%
        as.data.frame() %>%
        set_rownames(abundance_salmon$feature_id) %>%
        dplyr::select(-feature_id) %>%
        as.matrix()
      
      lengthMat <- sapply(
        seq_len(ncol(tpm)),
        function(x) {dtelist$genes$length}
      ) %>%
        set_rownames(rownames(dtelist$counts)) %>%
        set_colnames(colnames(dtelist$counts))
      
      tx_counts <- tximport::makeCountsFromAbundance(
        countsMat = dtelist$counts,
        abundanceMat = tpm,
        lengthMat = lengthMat,
        countsFromAbundance = "lengthScaledTPM"
      )
      
      rna_counts <- tx_counts[rownames(dtelist$counts), ]
      
    } else if (tpm_method == "standard") {
      
      rna_counts <- tpm(dgelist)
      
    } else if (tpm_method == "none") {
      
      message("No tpm selected, assuming gene mode")
      
    } else {
      
      print("Either specify a path or just use standard")
      
    }
  } else if (level == "gene") {
    message("Running in gene mode, TPM not calculated...")
    rna_counts <- cpm(dgelist,
                      log = TRUE)
  }
  
  version_name_change <- dtelist$genes %>%
    dplyr::select("grch38_name" = gene_name,
                  "grch37_name" = transcript_name) %>%
    mutate(grch37_name = str_replace(grch37_name, "\\-.*", ""),
           same_name = grch38_name == grch37_name) %>%
    distinct()
  
  if(isFALSE(version_name_change %>%
             dplyr::filter(
               grch37_name == target_gene
             ) %>% .[["same_name"]]
  )) {
    target_gene_38 <- version_name_change %>%
      dplyr::filter(
        grch37_name == target_gene
      ) %>% .[["grch38_name"]]
  } else {
    target_gene_38 <- target_gene
  }
  
  metadata <- dgelist$samples
  gene_anno <- dgelist$genes
  
  gene_info <- gene_anno[gene_anno$gene_name %in% target_gene_38, ]
  
  rna_counts <- rownames_to_column(
    as.data.frame(rna_counts),
    if (level == "transcript") {
      var = "transcript_id"
    } else if (level == "gene") {
      var = "gene_id"
    }
  )
  
  rna_counts <- if (level == "transcript") {
    rna_counts[
      rna_counts$transcript_id %in% gene_info$transcript_id,
    ]
  } else if (level == "gene") {
    rna_counts[
      rna_counts$gene_id %in% gene_info$gene_id,
    ]
  }
  
  rna_ggplot <- pivot_longer(
    if (level == "transcript") {
      left_join(rna_counts,
                dplyr::select(gene_info,
                              "transcript_id",
                              "transcript_name"),
                by = "transcript_id")
    } else if (level == "gene") {
      left_join(rna_counts,
                dplyr::select(gene_info,
                              "gene_id",
                              "gene_name"),
                by = "gene_id")
    },
    cols = metadata$ID,
    names_to = "ID",
    values_to = "counts"
  )
  
  rna_ggplot <- left_join(
    rna_ggplot,
    dplyr::select(metadata,
                  "ID",
                  all_of(variable)),
    by = "ID"
  )
  
  if (level == "gene") {
    rna_ggplot <- rna_ggplot %>%
      left_join(version_name_change,
                by = c("gene_name" = "grch38_name")) %>%
      dplyr::rename("grch38_name" = "gene_name") %>%
      dplyr::rename("gene_name" = "grch37_name")
  } else {
    message("Transcript names are in GRCh37")
  }
  
  if (level == "transcript") {
    if (isTRUE(is.numeric(metadata[[variable]]))) {
      ggplot(rna_ggplot) +
        geom_point(
          aes(x = .data[[variable]],
              y = log2(counts),
              colour = transcript_name)
        ) +
        geom_smooth(
          aes(x = .data[[variable]],
              y = log2(counts),
              colour = transcript_name,
              fill = transcript_name)
        ) +
        scale_colour_viridis_d() +
        scale_fill_viridis_d() +
        expand_limits(y = 0) +
        xlab(variable) +
        ylab("Transcript Expression (log2 TPM)") +
        ggtitle(
          paste0(target_gene,
                 " transcript expression")
        ) +
        theme_bw() +
        theme(
          legend.title = element_blank()#,
          # legend.text = element_text(size = 15,
          #                            colour = "black"),
          # legend.key.size = unit(2, "lines"),
          # axis.title.x = element_text(size = 15,
          #                             colour = "black"),
          # axis.text.x = element_text(size = 15,
          #                            colour = "black"),
          # axis.title.y = element_text(size = 15,
          #                             colour = "black"),
          # axis.text.y = element_text(size = 15,
          #                            colour = "black"),
          # plot.title = element_text(size = 15,
          #                           face = "bold")
        ) +
        if (isFALSE(keep_legend)) {
          theme(legend.position = "none")
        }
    } else {
      ggplot(rna_ggplot) +
        geom_boxplot(
          aes(x = .data[[variable]],
              y = log2(counts),
              fill = transcript_name),
          colour = "black"
        ) +
        scale_fill_viridis_d() +
        expand_limits(y = 0) +
        xlab(variable) +
        ylab("Transcript Expression (log2 TPM)") +
        ggtitle(
          paste0(target_gene, " transcript expression")
        ) +
        facet_grid(~ .data[[variable]],
                   scales = "free") +
        theme_bw() +
        theme(
          legend.title = element_blank()#,
          # legend.text = element_text(size = 15,
          #                            colour = "black"),
          # legend.key.size = unit(2, "lines"),
          # axis.title.x = element_text(size = 15,
          #                             colour = "black"),
          # axis.text.x = element_text(size = 15,
          #                            colour = "black"),
          # axis.title.y = element_text(size = 15,
          #                             colour = "black"),
          # axis.text.y = element_text(size = 15,
          #                            colour = "black"),
          # plot.title = element_text(size = 15,
          #                           face = "bold")
        ) +
        if (isFALSE(keep_legend)) {
          theme(legend.position = "none")
        }
    }
    
  } else if (level == "gene") {
    if (isTRUE(is.numeric(metadata[[variable]]))) {
      ggplot(rna_ggplot) +
        geom_point(
          aes(x = .data[[variable]],
              y = counts,
              colour = gene_name)
        ) +
        geom_smooth(
          aes(x = .data[[variable]],
              y = counts,
              colour = gene_name,
              fill = gene_name)
        ) +
        scale_colour_viridis_d() +
        scale_fill_viridis_d() +
        expand_limits(y = 0) +
        xlab(variable) +
        ylab("Gene Expression (log2 CPM)") +
        ggtitle(
          paste0(target_gene, " gene expression")
        ) +
        theme_bw() +
        theme(
          legend.title = element_blank()#,
          # legend.text = element_text(size = 15,
          #                            colour = "black"),
          # legend.key.size = unit(2, "lines"),
          # axis.title.x = element_text(size = 15,
          #                             colour = "black"),
          # axis.text.x = element_text(size = 15,
          #                            colour = "black"),
          # axis.title.y = element_text(size = 15,
          #                             colour = "black"),
          # axis.text.y = element_text(size = 15,
          #                            colour = "black"),
          # plot.title = element_text(size = 15,
          #                           face = "bold")
        ) +
        if (isFALSE(keep_legend)) {
          theme(legend.position = "none")
        }
    } else {
      ggplot(rna_ggplot) +
        geom_boxplot(
          aes(x = .data[[variable]],
              y = log2(counts),
              fill = gene_name),
          colour = "black"
        ) +
        scale_fill_viridis_d() +
        expand_limits(y = 0) +
        xlab(variable) +
        ylab("Gene Expression (log2 CPM)") +
        ggtitle(
          paste0(target_gene, " gene expression")
        ) +
        theme_bw() +
        theme(
          legend.title = element_blank()#,
          # legend.text = element_text(size = 15,
          #                            colour = "black"),
          # legend.key.size = unit(2, "lines"),
          # axis.title.x = element_text(size = 15,
          #                             colour = "black"),
          # axis.text.x = element_text(size = 15,
          #                            colour = "black"),
          # axis.title.y = element_text(size = 15,
          #                             colour = "black"),
          # axis.text.y = element_text(size = 15,
          #                            colour = "black"),
          # plot.title = element_text(size = 15,
          #                           face = "bold")
        ) +
        if (isFALSE(keep_legend)) {
          theme(legend.position = "none")
        }
    }
  }
}

