
source(here("R/plot_proportions.R"))

dtelist = dtelist
target_gene = "ADAM10"
variable = "GestationGroup"
dtu_variable = "GestationGroup"
abundance_path = here("data/counts/salmon_abundance.csv.gz")
tpm_method = "tximport"
dtu_x_label = "Gestation Group (weeks)"

annotation_df <- data.frame(
  Transcript = c("ADAM10-001", "ADAM10-015"),
  start = "6-10",
  end = "11-23",
  y = c(0.82),
  label = "*",
  panel = "a"
)

white_rect <- data.frame(
  xmin = 0,
  xmax = 1,
  ymin = 0.7,
  ymax = 0.9,
  panel = "b"
)

# Note that transcript_name here had to be the same as the facet code
anno <- data.frame(x1 = c(1, 0, 0, 1),
                   x2 = c(2, 0, 0, 2), 
                   y1 = c(0.78, 0, 0, 0.78),
                   y2 = c(0.84, 0, 0, 0.84), 
                   xstar = c(1.5, 0, 0, 1.5),
                   ystar = c(0.86, 0, 0, 0.86),
                   lab = c("*", "", "", "*"),
                   transcript_name = c("ADAM10-001",
                                       "ADAM10-002",
                                       "ADAM10-008",
                                       "ADAM10-015"))

proportion_plot <- plot_proportions(
  dtelist = dtelist,
  target_gene = target_gene,
  variable = if_else(is.null(dtu_variable), variable, dtu_variable),
  tpm_method = tpm_method,
  path = abundance_path,
  keep_legend = FALSE
) +
  labs(x = if_else(is.null(dtu_x_label), variable, dtu_x_label)) +
  ggtitle(NULL) +
  theme_bw()

proportion_plot +
  scale_y_continuous(limits = c(0, 0.9)) +
  geom_text(data = anno,
            aes(x = xstar,
                y = ystar,
                label = lab),
            size = 7) +
  geom_segment(data = anno,
               aes(x = x1, 
                   xend = x1,
                   y = y1,
                   yend = y2),
               colour = "black") +
  geom_segment(data = anno,
               aes(x = x2,
                   xend = x2,
                   y = y1,
                   yend = y2)) +
  geom_segment(data = anno,
               aes(x = x1,
                   xend = x2, 
                   y = y2,
                   yend = y2),
               colour = "black")

## Now just need to add this into my code
## I'm going to remake the ADAM10 script so I can add this in!
