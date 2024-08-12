#' @title Plot transcript structures
#'
#' @description Plot the structures of each transcript
#'
#' @param dgelist A dgelist with counts, gene annotations and samples
#' @param target_gene gene of interest
#' @param annotation an ensembl annotation object
#' @param rescale_introns
#' @param new_intron_length
#' @param flanking_length
#' @param connect_exons
#' @param transcript_label Should transcripts be labelled?
#' @param region_coords

plot_structures <- function(
  dtelist,
  target_gene,
  annotation,
  rescale_introns = TRUE,
  new_intron_length = 50,
  flanking_length = c(50, 50),
  connect_exons = TRUE,
  transcript_label = TRUE,
  region_coords = NULL
) {
  
  require(ggplot2)
  require(dplyr)
  require(readr)
  require(stringr)
  require(tidyr)
  require(tibble)
  require(magrittr)
  require(wiggleplotr)
  require(assertthat)
  require(tximport)
  require(here)
  
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
  
  target_gene_37 <- target_gene

  tx_names <- dtelist$genes %>%
    dplyr::filter(gene_name == target_gene_38) %>%
    dplyr::select(transcript_id) %>%
    as.matrix() %>%
    as.character()

  tx_annot <- wiggleplotr:::extractTranscriptAnnotationsFromEnsembldb(
    annotation,
    target_gene_37,
    tx_names
  )

  exons <- tx_annot$exons
  cdss <- tx_annot$cdss
  transcript_annotations <- tx_annot$transcript_annotations

  if(is.null(cdss)){
    cdss = exons
  }

  # Check exons and cdss
  assertthat::assert_that(
    is.list(exons)|| is(exons, "GRangesList")
  )

  # Check that exons and cdss objects are lists
  assertthat::assert_that(
    is.list(cdss) || is(exons, "GRangesList")
  )

  # Join exons together
  joint_exons = wiggleplotr:::joinExons(exons)

  # Extract chromosome name
  chromosome_name = as.vector(
    GenomicRanges::seqnames(joint_exons)[1]
  )

  # If region_coords is specificed, then ignore the flanking_length attrbute and compute
  # flanking_length form region_coords
  if (!is.null(region_coords)) {
    gene_range = wiggleplotr:::constructGeneRange(
      joint_exons,
      c(0,0)
    )
    min_start = min(
      GenomicRanges::start(gene_range)
    )
    max_end = max(
      GenomicRanges::end(gene_range)
    )
    flanking_length = c(
      min_start - region_coords[1],
      region_coords[2] - max_end
    )
  }
  # Make sure that flanking_length is a vector of two elements
  assertthat::assert_that(
    length(flanking_length) == 2
    )

  # Rescale introns
  if (rescale_introns) {
    tx_annotations = wiggleplotr:::rescaleIntrons(
      exons,
      cdss,
      joint_exons,
      new_intron_length = new_intron_length,
      flanking_length
    )
    xlabel = "Distance from region start (bp)"
  } else {
    old_introns = wiggleplotr:::intronsFromJointExonRanges(
      GenomicRanges::ranges(
        joint_exons
      ),
      flanking_length = flanking_length
    )
    tx_annotations = list(
      exon_ranges = lapply(
        exons,
        GenomicRanges::ranges
      ),
      cds_ranges = lapply(
        cdss,
        GenomicRanges::ranges
      ),
      old_introns = old_introns,
      new_introns = old_introns
    )

    xlabel = paste(
      "Chromosome",
      chromosome_name,
      "position (bp)"
    )
  }

  # If transcript annotations are not supplied then construct them manually from the GRanges list
  if (is.null(transcript_annotations)) {
    plotting_annotations = dplyr::tibble(
      transcript_id = names(exons),
      strand = extractStrandsFromGrangesList(exons)
    ) %>%
      prepareTranscriptAnnotations()
  } else {
    plotting_annotations = wiggleplotr:::prepareTranscriptAnnotations(
      transcript_annotations
      )
  }

  # Plot transcript structures
  limits = c(
    min(
      IRanges::start(tx_annotations$new_introns)
    ),
    max(
      IRanges::end(tx_annotations$new_introns)
    )
  )
  structure = wiggleplotr:::prepareTranscriptStructureForPlotting(
    tx_annotations$exon_ranges,
    tx_annotations$cds_ranges,
    plotting_annotations
  )

  structure$transcript_label <- structure$transcript_label %>%
    gsub(
      pattern = ".*:",
      replacement = ""
    ) %>%
    gsub(
      pattern = " >",
      replacement = ""
    ) %>%
    as.data.frame() %>%
    set_colnames("IDs") %>%
    left_join(
      dtelist$genes,
      by = c("IDs" = "transcript_id")
    ) %>%
    dplyr::select("transcript_name") %>%
    as.matrix() %>%
    as.character()

  exons_df <- structure
  limits <- limits
  connect_exons <- connect_exons
  transcript_label <- transcript_label

  # Extract the position for plotting transcript name
  transcript_annot = dplyr::group_by(
    exons_df,
    transcript_id
  ) %>%
    dplyr::filter(
      feature_type == "exon"
    ) %>%
    dplyr::arrange(
      'transcript_id',
      'start'
    ) %>%
    dplyr::filter(
      row_number() == 1
    ) %>%
    mutate(
      colour_info = paste0(
        transcript_label,
        feature_type
      )
    )

  exons_df <- exons_df %>%
    mutate(
      colour_info = paste0(
        transcript_label,
        feature_type
      )
    )

  exons_df <- structure %>%
    dplyr::filter(
      feature_type == "exon"
    )

  cds_df <- structure %>%
    dplyr::filter(
      feature_type == "cds"
    )

  exons_df <- exons_df %>%
    mutate(
      transcript_label = replace_na(
        transcript_label,
        "LRG_151t1"
      )
    )

  transcript_annot <- transcript_annot %>%
    mutate(
      transcript_label = replace_na(
        transcript_label,
        "LRG_151t1"
      )
    )

  exons_df %>%
    ggplot() +
    geom_rect(
      aes(
        xmin = start,
        xmax = end,
        ymax = transcript_rank+0.25,
        ymin = transcript_rank-0.25,
        fill = transcript_label
      ),
      alpha = 0.5
    ) +
    geom_rect(
      aes(
        xmin = start,
        xmax = end,
        ymax = transcript_rank+0.25,
        ymin = transcript_rank-0.25,
        colour = transcript_label,
        fill = transcript_label
      ),
      data = cds_df
    ) +
    geom_line(
      aes(
        x = start,
        y = transcript_rank,
        group = transcript_rank,
        colour = transcript_label
      )
    ) +
    geom_line(
      aes(
        x = end,
        y = transcript_rank,
        group = transcript_rank,
        colour = transcript_label
      )
    ) +
    xlab(xlabel) +
    theme_bw() +
    theme(
      plot.margin = unit(
        c(0, 0, 0, 0),
        "line"
      ),
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      legend.position="none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.text.y = element_text(
        colour = "black"
      ),
      strip.background = element_rect(
        fill = "grey85"
      )
    ) +
    scale_fill_viridis_d() +
    scale_colour_viridis_d() +
    geom_text(
      aes(
        x = start,
        y = transcript_rank + 0.30,
        label = transcript_label
      ),
      data = transcript_annot,
      hjust = 0,
      vjust = 0,
      size = 4
    )
}
