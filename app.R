

source('src/functions.R')


options(shiny.maxRequestSize = 100*1024^2)
ui <- fluidPage(
    
    titlePanel("PM Detection Algorithm Demo - O'Donnell Lab", windowTitle = "PM Detection Algorithm Demo - O'Donnell Lab"),
   wellPanel(
            #h4("Pipeline Options", align = "center"),
            fluidRow(
                #column(1, strong("1.")),
                column(2, br(), fileInput('input_image', NULL, accept = c('.tif'),multiple=FALSE)),
                column(1, br(), actionButton("inputimage_help", "?")),
            
               # column(1, strong("2.")),
                column(1, textInput("cmac_chan", "CMAC", value = "", placeholder = "1")),
                column(1, textInput("gfp_chan", "GFP", value = "", placeholder = "2")),
                column(1, textInput("dic_chan", "DIC", value = "", placeholder = "3")),
                column(1, br(), actionButton("inputchannels_help", "?")),
          
               # column(1, strong("3.")),
                column(1, textInput("cutoff_value", "Cutoff", placeholder="100")),
                column(1, br(), actionButton("cutoffvalue_help", "?")),
            
            column(3, br(), fluidRow(actionButton("show_options", "Demonstate Options")))),
   ),
       
        
          #plotOutput("options_demo"),
          h6("GFP Algorithm"),
          fluidRow(
            column(2,plotOutput("gfp1")),
            column(2,plotOutput("gfp2")),
            column(2,plotOutput("gfp4")),
            column(2,plotOutput("gfp8")),
            column(2,plotOutput("gfp16"))),
          h6("DIC Algorithm"),
          fluidRow(
            column(2,plotOutput("dic1")),
            column(2,plotOutput("dic2")),
            column(2,plotOutput("dic4"))),
          
          
    
        
        
        
        
 
    


   
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {
    
    values <- reactiveValues()
    
    # help buttons
    observeEvent(input$inputimage_help, {
        showModal(modalDialog(
            title = "Input Image",
            "Select a .tif file.",
            easyClose = TRUE,
            footer = NULL
        ))
    })
    observeEvent(input$inputchannels_help, {
        showModal(modalDialog(
            title = "Input Channels",
            "Frame numbers for the different channels of the .tif file.  If there is no DIC channel, leave the field blank.",
            easyClose = TRUE,
            footer = NULL
        ))
    })
    observeEvent(input$cutoffvalue_help, {
        showModal(modalDialog(
            title = "Cell Area Cutoff",
            HTML(paste0("This should be obtained in imageJ using the \"Measure\" tool.  See the Automated Quant Tutorial for more details.")),
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
      else if( !(str_trim(input$cmac_chan) %in% c("1", "2", "3")))
      {
        showModal(modalDialog(title = "Invalid Input", "Invalid value for CMAC channel.", easyClose = TRUE, footer = NULL))
      }
      else if( !(str_trim(input$gfp_chan) %in% c("1", "2", "3")))
      {
        showModal(modalDialog(title = "Invalid Input", "Invalid value for GFP channel.", easyClose = TRUE, footer = NULL))
      }
      else if( !(str_trim(input$dic_chan) %in% c("1", "2", "3", "")))
      {
        showModal(modalDialog(title = "Invalid Input", "Invalid value for DIC channel.", easyClose = TRUE, footer = NULL))
      }
      else if(any(duplicated(c(input$cmac_chan, input$dic_chan, input$gfp_chan))))
      {
        showModal(modalDialog(title="Invalid Input", "Different channels cannot have the same value", easyClose = TRUE, footer = NULL))
      }
      else
      {
        progress <- shiny::Progress$new()
        
        on.exit(progress$close())
        
        progress$set(message = "Running Test Cases...", value = 0)
        
        values$options <- pipeline_options(input$input_image,
                                           progress = progress,
                                           gfp_chan = input$gfp_chan,
                                           dic_chan = input$dic_chan,
                                           cmac_chan = input$cmac_chan,
                                           cutoff = input$cutoff_value
                                           )
        if(!is.null(values$options))
        {
          if(!is.null(values$options$gfp))
          {
            output$gfp1 <- get_test(values$options$gfp, values$options$channels, 1)
            output$gfp2 <- get_test(values$options$gfp, values$options$channels,2)
            output$gfp4 <- get_test(values$options$gfp, values$options$channels,4)
            output$gfp8 <- get_test(values$options$gfp, values$options$channels,8)
            output$gfp16 <- get_test(values$options$gfp, values$options$channels,16)
          }
          if(!is.null(values$options$dic))
          {
            output$dic1 <- get_test(values$options$dic, values$options$channels,1)
            output$dic2 <- get_test(values$options$dic, values$options$channels,2)
            output$dic4 <- get_test(values$options$dic, values$options$channels,4)
          }
        }
      }
    }) 
    

    

    
    # Close the app when the session completes
    if(!interactive()) {
      session$onSessionEnded(function() {
        stopApp()
        q("no")
      })
    }
 

    
}

# Run the application 
shinyApp(ui = ui, server = server)
