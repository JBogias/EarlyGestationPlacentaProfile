library(shiny)
library(here)

ui <- shinyUI(
  fluidPage(
    includeHTML(here("www/EarlyGestation_transcriptomeProfile.html"))
  )
)

server <- function(input, output) {}

shinyApp(ui, server)
