library(shiny)
library(here)

ui <- shinyUI(
  fluidPage(
    includeHTML(here("www/EarlyGestation_SupplementaryTables.html"))
  )
)

server <- function(input, output) {}

shinyApp(ui, server)