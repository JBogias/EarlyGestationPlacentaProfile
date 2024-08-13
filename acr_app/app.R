# App for ACR and App Service deployment
# I would just use the docker-compose for ACR but then I'd need to store 
# all seven of the images in the ACR, which would rake up costs
# App service is free though
library(shiny)
library(here)

addResourcePath(prefix = "html", directoryPath = "www/")

ui <- fluidPage(
  class = "fluid-page",
  tags$head(
    tags$link(rel = "stylesheet", type = "text/css", href = "html/styling.css")
  ),
  titlePanel(paste0("Code Companion for 'Placental Transcription Profiling ",
                      "in 6–23 Weeks’ Gestation Reveals Differential",
                      " Transcript Usage in Early Development'",
                      " by Justin Bogias")),
  tabsetPanel(
    tabPanel("Differential Expression",
             class = "tab-panel",
             tags$iframe(src = "html/EarlyGestation_RNAseq.html",
                         style = "width:70vw;height:80vh;",
                         class = "markdown_frame")),
    tabPanel("Differential Transcript Usage",
             class = "tab-panel",
             tags$iframe(src = "html/EarlyGestation_DTU.html",
                         style = "width:70vw;height:80vh;",
                         class = "markdown_frame")),
    tabPanel("Results Differences",
             class = "tab-panel",
             tags$iframe(
               src = "html/EarlyGestation_transcriptomeProfile.html",
               style = "width:70vw;height:80vh;",
               class = "markdown_frame")),
    tabPanel("ADAM10 Profile",
             class = "tab-panel",
             tags$iframe(src = "html/EarlyGestation_ADAM10profile.html",
                         style = "width:70vw;height:80vh;",
                         class = "markdown_frame")),
    tabPanel("Gene Ontology Analysis",
             class = "tab-panel",
             tags$iframe(src = "html/EarlyGestation_GOanalysis.html",
                         style = "width:70vw;height:80vh;",
                         class = "markdown_frame")),
    tabPanel("Supplementary Figures",
             class = "tab-panel",
             tags$iframe(
               src = "html/EarlyGestation_SupplementaryFigures.html",
               style = "width:70vw;height:80vh;",
               class = "markdown_frame")
             ),
    tabPanel("Supplementary Tables",
             class = "tab-panel",
             tags$iframe(
               src = "html/EarlyGestation_SupplementaryTables.html",
               style = "width:70vw;height:80vh;",
               class = "markdown_frame")
    ),
    tabPanel("Publication",
             class = "tab-panel",
             tags$iframe(
               src = "html/ijms-23-04506-v4.pdf",
               style = "width:70vw;height:80vh;",
               class = "markdown_frame"
             ))
  )
)

server <- function(input, output, session) {
}

shinyApp(ui, server)