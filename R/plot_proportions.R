#' @title Plot transcript proportions
#'
#' @description Calculate proportions of transcripts per gene and visualise global DTU profile
#'
#' @param dtelist either a DGEList or SummarizedExperiment
#' @param variable variable of interest to test the data
#' @param target_gene gene of interest to visualise proportions


require(ggplot2)
require(tidyverse)
require(tximport)
require(here)

plot_proportions <- function (
    dtelist,
    target_gene,
    variable,
    tpm_method = "standard",
    path = NULL,
    keep_legend = TRUE
) {

    if (tpm_method == "standard") {

        rpk <- dtelist$counts/(dtelist$genes$length/1000)
        scaling <- colSums(rpk)/1e6
        tpm <- rpk/scaling

    } else if (tpm_method == "tximport") {
        abundance_salmon <- read_csv(path)

        tpm <- abundance_salmon %>%
            as.data.frame() %>%
            set_rownames(abundance_salmon$feature_id) %>%
            dplyr::select(-feature_id) %>%
            as.matrix()
    } else {
        stop(
            "Need either abundance from 'tximport' or 'standard' method"
            )
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

    gene_anno <- dtelist$genes

    metadata <- dtelist$samples

    metadata$variable <- metadata[[variable]]

    lengthMat <- sapply(seq_len(ncol(tpm)), function(x) {
        dtelist$genes$length
    })
    rownames(lengthMat) <- rownames(dtelist$counts)
    colnames(lengthMat) <- colnames(dtelist$counts)

    dtu_counts <- makeCountsFromAbundance(
        dtelist$counts,
        tpm,
        lengthMat,
        "lengthScaledTPM"
    )

    if (isTRUE(is.numeric(metadata[[variable]]))) {

        joined_data <- left_join(
            rownames_to_column(
                as.data.frame(
                    prop.table(
                        dtu_counts,
                        margin = 2
                    )
                ),
                var = "transcript_id"
            ),
            gene_anno,
            by = "transcript_id"
        )
        rownames(joined_data) <- joined_data$transcript_id

        prop2plot <- left_join(
            pivot_longer(
                rownames_to_column(
                    as.data.frame(
                        prop.table(
                            as.matrix(
                                joined_data[
                                    joined_data$gene_name %in% target_gene_38,
                                    colnames(joined_data) %in% metadata$ID
                                ]
                            ),
                            margin = 2
                        )
                    ), var = "transcript_id"
                ),
                cols = starts_with("PAC"),
                names_to = "ID",
                values_to = "Proportions"
            ),
            metadata,
            by = "ID"
        )

        prop2plot <- left_join(
            prop2plot,
            gene_anno[ c("transcript_id", "transcript_name") ],
            by = "transcript_id"
        )

        ggplot(
            prop2plot,
            aes(
                x = .data[[variable]],
                y = Proportions,
                fill = transcript_name
            )
        ) +
            geom_bar(
                position = "fill",
                stat = "identity",
                width = 0.95
            ) +
            ggtitle(
                paste0(
                    target_gene,
                    " transcript proportions"
                )
            ) +
            xlab(paste0(variable)) +
            ylab("Proportions") +
            scale_fill_viridis_d() +
            # facet_grid(
            #     ~factor(),
            #     scales = "free_x",
            #     space = "free",
            #     shrink = FALSE
            # ) +
            theme(
                # axis.text.x = element_blank(),
                # axis.ticks.x = element_blank(),
                # legend.title = element_text(
                #     size = 15,
                #     colour = "black"
                # ),
                # legend.text = element_text(
                #     size = 15,
                #     colour = "black"
                # ),
                # axis.title = element_text(
                #     size = 15,
                #     colour = "black"
                # ),
                # axis.text = element_text(
                #     size = 15,
                #     colour = "black"
                # ),
                # strip.text.x = element_text(
                #     size = 15
                # ),
                # plot.title = element_text(
                #     size = 15,
                #     face = "bold"
                # ),
                panel.background = element_rect(
                    fill = "white",
                    colour = NA
                ),
                panel.grid = element_line(
                    colour = "white"
                )
            ) +
            if (isFALSE(keep_legend)) {
                theme(
                    legend.position = "none"
                )
            }

    } else {

        # Here we get the proportions of transcripts making up a gene
        box_props <- left_join(
            rownames_to_column(
                as.data.frame(
                    prop.table(
                        dtu_counts,
                        margin = 2
                    )
                ), var = "transcript_id"
            ),
            gene_anno[c("gene_id", "transcript_id", "gene_name")],
            by = "transcript_id"
        )

        box_props <- box_props[
            box_props$gene_name %in% target_gene_38,
        ]
        rownames(box_props) <- box_props$transcript_id

        # Here we calculate the proportions for the gene of interest
        # There is something in the way I plot these which cause the graphic to
        # change, I've check and found that the dataframe says the same thing
        # as the original one in the markdown
        box_props <- left_join(
            pivot_longer(
                rownames_to_column(
                    as.data.frame(
                        prop.table(
                            as.matrix(
                                box_props[metadata$ID]),
                            margin = 2)
                    ),
                    var = "transcript_id"),
                cols = starts_with("PAC"),
                names_to = "ID",
                values_to = "Proportions"),
            metadata,
            by = "ID"
        )

        box_props <- left_join(
            box_props,
            gene_anno[c("transcript_id", "transcript_name")],
            by = "transcript_id"
        )


        ggplot(box_props) +
            geom_boxplot(
                aes(
                    x = variable,
                    y = Proportions,
                    fill = transcript_name
                ),
                colour = "black"
            ) +
            ggtitle(
                paste0(
                    target_gene,
                    " transcript proportions"
                )
            ) +
            xlab(
                paste0(variable)
            ) +
            ylab(
                "Median transcript proportions"
            ) +
            scale_fill_viridis_d() +
            facet_grid(
                ~factor(transcript_name),
                scales = "free_x",
                space = "free",
                shrink = FALSE
            ) +
            theme_bw() +
            theme(
                axis.ticks.x = element_blank(),
                # legend.title = element_text(
                #     size = 15,
                #     colour = "black"
                # ),
                # legend.text = element_text(
                #     size = 15,
                #     colour = "black"
                # ),
                # axis.title.x = element_text(
                #     size = 15,
                #     colour = "black"
                # ),
                # axis.text.x = element_text(
                #     size = 14,
                #     colour = "black"
                # ),
                # axis.text.y = element_text(
                #     size = 15,
                #     colour = "black"
                # ),
                # axis.title.y = element_text(
                #     size = 15,
                #     colour = "black"
                # ),
                # strip.text.x = element_text(
                #     size = 10
                # ),
                # plot.title = element_text(
                #     size = 15,
                #     face = "bold"
                # ),
                panel.background = element_rect(
                    fill = "white",
                    colour = NA
                ),
                panel.grid = element_line(
                    colour = "white"
                )
            ) +
            if (isFALSE(keep_legend)) {
                theme(
                    legend.position = "none"
                )
            }
    }
}




