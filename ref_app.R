

library(shiny)
library(dplyr)
library(stringr)
source("generate_output.R")

options(shiny.maxRequestSize = 9*1024^2)

# Define UI for application that draws a histogram
ui <- fluidPage(
  
  
  # Application title
  fileInput('image_file', "Choose Image .csv file",
            accept = c(
              'text/csv',
              'text/comma-separated-values',
              'text/tab-separated-values',
              'text/plain',
              '.csv',
              '.tsv'
            )
  ),
  
  fileInput('object_file', "Choose EditedObjects .csv file",
            accept = c(
              'text/csv',
              'text/comma-separated-values',
              'text/tab-separated-values',
              'text/plain',
              '.csv',
              '.tsv'
            )
  ),
  fluidRow(
    column(3,actionButton("run", "Run Script")),
    column(3,uiOutput("download"))),
  
  tableOutput("contents")
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  
  values <- reactiveValues()
  observeEvent(input$run,{
    
    output$contents <- renderTable({
      # input$file1 will be NULL initially. After the user selects
      # and uploads a file, it will be a data frame with 'name',
      # 'size', 'type', and 'datapath' columns. The 'datapath'
      # column will contain the local filenames where the data can
      # be found.
      
      image_file <- input$image_file
      object_file <- input$object_file
      
      if (is.null(image_file) || is.null(object_file))
        return(NULL)
      
      values$res <- generate_csv(object_file$datapath, image_file$datapath)
      
    })
  })
  
  output$download <- renderUI(
    
    if(!is.null(values$res))
    {
      downloadButton("downloadButton", "Download Results")
    }
  )
  
  output$downloadButton <- downloadHandler(
    filename = "Results.csv",
    content = function(file) {
      write.csv(values$res, file, row.names = FALSE)
    }
  )
  
  
}

# Run the application 
shinyApp(ui = ui, server = server)
