
library(shiny)
library(dplyr)
library(stringr)
source('src/functions.R')
source('src/Main.R')

options(shiny.maxRequestSize = 9*1024^2)

ui <- fluidPage(
    
    titlePanel("Automated Cell Quant - O'Donnell Lab", windowTitle = "Automated Cell Quant - O'Donnell Lab"),
    sidebarLayout(
        sidebarPanel(
            h4("Pipeline Options", align = "center"),
            fluidRow(
                column(1, strong("1.")),
                column(9, fileInput('input_images', NULL, accept = c('.tif'),multiple=FALSE)),
                column(2, actionButton("inputimage_help", "?"))),
            fluidRow(
                column(1, strong("2.")),
                column(3, textInput("cmac_chan", "CMAC", value = "", placeholder = "1")),
                column(3, textInput("gfp_chan", "GFP", value = "", placeholder = "2")),
                column(3, textInput("dic_chan", "DIC", value = "", placeholder = "3")),
                column(2, br(), actionButton("inputchannels_help", "?"))),
            fluidRow(
                column(1, strong("3.")),
                column(9, textInput("cutoff_value", "Cell size cutoff", placeholder="100")),
                column(2, br(), actionButton("cutoffvalue_help", "?"))),
            fluidRow(actionButton("run", "Run Pipeline"), align = "center"),
        width=4),
        mainPanel(h4("main panel", align = "center"), width = 8)
        
        
        
        
    ),
    

    fluidRow(
        
        column(5,uiOutput("contents"))

   
))

# Define server logic required to draw a histogram
server <- function(input, output) {
    
    values <- reactiveValues()
    
    # help buttons
    observeEvent(input$inputimage_help, {
        showModal(modalDialog(
            title = "modal title input",
            "this is the help text for input",
            easyClose = TRUE,
            footer = NULL
        ))
    })
    observeEvent(input$inputchannels_help, {
        showModal(modalDialog(
            title = "modal title channel",
            "this is the help text for channel",
            easyClose = TRUE,
            footer = NULL
        ))
    })
    observeEvent(input$cutoffvalue_help, {
        showModal(modalDialog(
            title = "modal title cutoff",
            "this is the help text for cutoff",
            easyClose = TRUE,
            footer = NULL
        ))
    })
   

        
    
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
    
   # output$contents <- renderUI(
        
      #  if(is.null(values$res))
     #       return(NULL)
        #else
          #  return(uiOutput(textOutput("DONE")))
  #  )
    
 

    
}

# Run the application 
shinyApp(ui = ui, server = server)
