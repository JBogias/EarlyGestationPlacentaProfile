library(shiny)
library(here)

ui <- shinyUI(
  fluidPage(
    includeHTML(here("www/EarlyGestation_DTU.html"))
  )
)

server <- function(input, output) {}

shinyApp(ui, server)
