library(shiny)

ui <- fluidPage(
  titlePanel("Early Gestation"),
  tabsetPanel(
    tabPanel("Differential Expression",
             tags$iframe(src = "http://localhost:4001",
                         style = "width:70vw;height:100vh;")),
    tabPanel("Differential Transcript Usage", 
             tags$iframe(src = "http://localhost:4002",
                         style = "width:70vw;height:100vh;")),
    tabPanel("Results Differences",
             tags$iframe(src = "http://localhost:4003",
                         style = "width:70vw;height:100vh;")),
    tabPanel("ADAM10 Profile",
             tags$iframe(src = "http://localhost:4004",
                         style = "width:70vw;height:100vh;")),
    tabPanel("Gene Ontology Analysis",
             tags$iframe(src = "http://localhost:4005",
                         style = "width:70vw;height:100vh;")),
    tabPanel("Supplementary Figures",
             tags$iframe(src = "http://localhost:4006",
                         style = "width:70vw;height:100vh;")),
    tabPanel("Supplementary Tables",
             tags$iframe(src = "http://localhost:4007",
                         style = "width:70vw;height:100vh;")),
  )
)

server <- function(input, output, session) {
}

shinyApp(ui, server)