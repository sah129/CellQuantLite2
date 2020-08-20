
library(shiny)
library(dplyr)
library(stringr)
library(gridExtra)

source('src/functions.R')
source('src/Main.R')
source('src/app_functions.R')

options(shiny.maxRequestSize = 100*1024^2)
ui <- fluidPage(
    
    titlePanel("Automated Cell Quant - O'Donnell Lab", windowTitle = "Automated Cell Quant - O'Donnell Lab"),
    sidebarLayout(
        sidebarPanel(
            h4("Pipeline Options", align = "center"),
            fluidRow(
                column(1, strong("1.")),
                column(9, fileInput('input_image', NULL, accept = c('.tif'),multiple=FALSE)),
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
            
            fluidRow(actionButton("show_options", "Demonstate Options for 4 and 5"), align = "center"),
            br(),
            
            fluidRow(
              column(1, strong("4.")),
              column(9, radioButtons("algchoose", "Membrane Detection Algorithm", choices = c("GFP", "DIC"), selected = "GFP", inline = TRUE)),
              column(2, br(), actionButton("algchoose_help", "?"))),
            fluidRow(
              column(1, strong("5.")),
              
              column(9, radioButtons("factorchoose", "Factor:", choices = c("1", "2", "4", "8", "16"), selected = "5", inline = TRUE)),
              
              
              column(2, br(), actionButton("factorchoose_help", "?"))),
            fluidRow(actionButton("run", "Run Pipeline"), align = "center"),
        width=4),
        mainPanel(
          h4("main panel", align = "center"),
          #plotOutput("options_demo"),
          h4("GFP Algorithm"),
          fluidRow(
            column(2,plotOutput("gfp1")),
            column(2,plotOutput("gfp2")),
            column(2,plotOutput("gfp4")),
            column(1,plotOutput("gfp8")),
            column(1,plotOutput("gfp16"))),
          h4("DIC Algorithm"),
          fluidRow(
            column(2,plotOutput("dic1")),
            column(2,plotOutput("dic2")),
            column(2,plotOutput("dic4")),
            column(1,plotOutput("dic8")),
            column(1,plotOutput("dic16"))),
          
          width = 8)
        
        
        
        
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
    observeEvent(input$algchoose_help, {
      showModal(modalDialog(
        title = "modal title algchoose",
        "this is the help text for algchoose",
        easyClose = TRUE,
        footer = NULL
      ))
    })
    observeEvent(input$factorchoose_help, {
      showModal(modalDialog(
        title = "modal title factorchoose",
        "this is the help text for factorchoose",
        easyClose = TRUE,
        footer = NULL
      ))
    })

       
    observeEvent(input$show_options,
    {
       if(is.null(input$input_image) || 
          is.null(input$gfp_chan) ||
          is.null(input$cmac_chan) ||
          is.null(input$dic_chan) ||
          is.null(input$cutoff_value))
        { 
         print("all not selected")
          return(NULL)       
       }
      else
      {
        values$options <- pipeline_options(input$input_image,
                                           gui = FALSE,
                                           progress = NULL,
                                           gfp_chan = input$gfp_chan,
                                           dic_chan = input$dic_chan,
                                           cmac_chan = input$cmac_chan,
                                           cutoff = input$cutoff_value
                                           )
        if(!is.null(values$options))
        {
          output$gfp1 <- get_test(values$options$gfp, 1)
          output$gfp2 <- get_test(values$options$gfp, 2)
          output$gfp4 <- get_test(values$options$gfp, 4)
          output$gfp8 <- get_test(values$options$gfp, 8)
          output$gfp16 <- get_test(values$options$gfp, 16)
          
          output$dic1 <- get_test(values$options$dic, 1)
          output$dic2 <- get_test(values$options$dic, 2)
          output$dic4 <- get_test(values$options$dic, 4)
          output$dic8 <- get_test(values$options$dic, 8)
          output$dic16 <- get_test(values$options$dic, 16)
        }
      }
    }) 
    

    
    observeEvent(input$run,
    {
        image_files <<- input$input_image
        
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
