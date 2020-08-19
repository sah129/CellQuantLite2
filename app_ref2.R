library(shiny)
library(shinyjs)
library(shinythemes)
library(shinyFiles)
library(EBImage)
library(ggplot2)
library(DT)

source("src/Main.R")
source("src/shiny_functions.R")

ui <- fluidPage(
  shinyjs::useShinyjs(),
  
  titlePanel( h1( "CellQuant: O'Donnell Lab", align = "center") , windowTitle = "2020 O'Donnell Lab"),
  fluidRow(
    column(11,  shinyDirButton('datasetpath', 'Select a directory containing the images for the pipeline', 'Please select a folder', FALSE, style = "width: 95%"), align = "center"),
    
    column(1, actionButton("process_dataset", "Run"))),
  tags$hr(),
  fluidRow(column(10, verbatimTextOutput("stream_main"))),
  uiOutput("results"),
  tags$hr(),
  fluidRow(column(12, h6("Created by Sarah Hawbaker, O'Donnell Lab, University of Pittsburgh, 2020.", align = "center")))
)


server <- function(input, output,session) 
{
  curr_img = -1
  v <- reactiveValues(res = NULL, log = "", i = -1)
  volumes <- c(def_root = "C:/Users/Sarah's Computer/Documents/GitHub/Cell-Quant-and-Loc/Datasets/", 
               def_dataset = "C:/Users/Sarah's Computer/Documents/GitHub/Cell-Quant-and-Loc/Datasets/Dummy 1B")
  
  shinyDirChoose(input, 'datasetpath',
                 roots = volumes,
                 filetypes=c(',', 'txt'),
                 session = session)
  #return(parseDirPath(getVolumes(),input$datasetpath))
  
  dpath <- reactive(return( parseDirPath(volumes["def_root"], input$datasetpath)))  
  
  observeEvent(input$remove_cells,{
    v$res <- remove_cells_interactive(v$res, v$i, input$cell_selections)
    output$edit_module <- get_edit_module(v$res[[v$i]])
    output$img_result_em <- get_edit_plot(v$res[[v$i]])
  })
  
  observeEvent(input$finish,{
    print("All Finished!")
    finish_up(v$res)
  })
  observeEvent(input$process_dataset, 
               {
                 print(input$datasetpath)
                 progress <- shiny::Progress$new()
                 on.exit(progress$close())
                 withCallingHandlers(
                   {
                     shinyjs::html("stream_main", "")
                     progress$set(message = "Running pipeline...", value = 0)
                     v$res <- pipeline_git1_interactive(dpath(),testing=FALSE, gui=TRUE, progress=progress) #dpath()!!!
                     output$stream_main <- renderText("")
                     output$results <- add_result_ui()
                   },
                   message = function(m) {
                     shinyjs::html(id = "stream_main", html = m$message, add = TRUE) #save the logfile somewhere
                     v$log <- paste0(v$log, m$message)
                   })
                 
               }, once = TRUE)
  output$results_log <- renderText(v$log)
  obsList <- list()
  output$image_select <- renderUI(
    {
      if(is.null(v$res)) 
        return()
      img_buttons <- as.list(1:length(v$res))
      img_buttons <- lapply(img_buttons, function(i)
      {
        img_name <- v$res[[i]]$filename
        
        if( is.null(obsList[[img_name]]))
        {
          obsList[[img_name]] <<- observeEvent(input[[img_name]], 
                                               {
                                                 v$i <- i
                                                 output$img_result_plot <- get_image_plot(v$res[[i]], input$toggle_selections, input$channel_sel)
                                                 output$mpi_table <-get_mpi_table(v$res[[i]])
                                                 output$sum_hist <- get_hist(v$res[[i]])
                                                 output$sum_text <- get_sum_text(v$res[[i]])
                                                 output$labeled_img <- get_final_labeled(v$res[[i]])
                                                 output$labeled_img_ex <- get_final_labeled(v$res[[i]])
                                                 output$img_result_em <- get_edit_plot(v$res[[i]])
                                                 output$edit_module <- get_edit_module(v$res[[i]])
                                                 observeEvent(input$toggle_selections, {
                                                   output$img_result_plot <- get_image_plot(v$res[[i]], input$toggle_selections, input$channel_sel)
                                                   
                                                 },ignoreNULL = FALSE)
                                                 observeEvent(input$channel_sel, {
                                                   output$img_result_plot <- get_image_plot(v$res[[i]], input$toggle_selections, input$channel_sel)
                                                 })
                                               }, ignoreInit = TRUE)
        }
        actionButton(img_name, break_file_name(v$res[[i]]$filename), width="100%")
      }) 
    })
}

shinyApp(ui, server)