library(shiny)
library(here)

ui <- shinyUI(
  fluidPage(
    includeHTML(here("www/EarlyGestation_ADAM10profile.html"))
  )
)

server <- function(input, output) {}

shinyApp(ui, server)
