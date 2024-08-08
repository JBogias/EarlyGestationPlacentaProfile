library(shiny)
library(here)

ui <- shinyUI(
  fluidPage(
    includeHTML(here("www/EarlyGestation_RNAseq_mod.html"))
  )
)

server <- function(input, output) {}

shinyApp(ui, server)
