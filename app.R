
library(shiny)
library(dplyr)
library(stringr)
source('src/functions.R')
source('src/Main.R')

options(shiny.maxRequestSize = 9*1024^2)

ui <- fluidPage(
    
    
    fileInput('input_images', "Choose tif files",
              accept = c(
                '.tif'
              ),
              multiple=TRUE

    ),
    fluidRow(
        column(3,actionButton("run", "Run Pipeline")),
        column(5,uiOutput("contents"))

   
))

# Define server logic required to draw a histogram
server <- function(input, output) {
    
    values <- reactiveValues()
    
   

        
    
    observeEvent(input$run,
    {
        image_files <<- input$input_images
        
        if (is.null(image_files))
        {
            print('null')
            return(NULL)
        }
        
       
            values$res <- pipeline(image_files,testing=FALSE, gui=FALSE, progress=NULL, interactive = FALSE, factor = 4, chan = gfp_channel, cutoff = 100) 
       
    })
    
    output$contents <- renderUI(
        
        if(is.null(values$res))
            return(NULL)
        else
            return(uiOutput(textOutput("DONE")))
    )
    
 

    
}

# Run the application 
shinyApp(ui = ui, server = server)
