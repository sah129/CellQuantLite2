get_mpi_table <- function(res)
{
  renderDT({res$df %>% select(1:7)
    
  }, rownames=FALSE)
}

break_file_name <- function(s)
{
  sp <- strsplit(s, "_")
  sp <- unlist(sp)
  

  s1 <- paste(sp[1], sp[2], sep = "_")
  
  s2 <- sp[-1]
  s2 <- s2[-1]
  s2 <- paste(s2, collapse = "_")
  
  HTML(paste(s1, s2, sep="<br/>"))
}


add_result_ui <-function()
{
  renderUI({
  fluidRow(
    column(2, h4("Results:"),uiOutput("image_select")),
    column(10,
           tabsetPanel(
             tabPanel("Image Result", 
                      br(),
                      sidebarLayout(
                    sidebarPanel(   
                      checkboxGroupInput("toggle_selections",
                                             label = "",
                                             choices = c("Membranes" = "mem_select",
                                                         "Vacuoles" = "vac_select",
                                                         "PM Labels" = "pm_label_select",
                                                         "Vac Labels" = "vac_label_select",
                                                         "Removed" = "rem_select")),
                      tags$hr(),
                      radioButtons("channel_sel", 
                                   label = "", 
                                   choices = c("CMAC" = "cmac",
                                               "GFP" = "gfp",
                                               "DIC" = "dic",
                                               "Grayscale GFP" = 'gray_gfp')), width = 2),
                   
                      #displayOutput("img_result")),
                    mainPanel(plotOutput("img_result_plot"),align="left", width = 4))),
            # tabPanel("Summary", fluidRow(br(),
             #                             column(3, uiOutput("sum_text")), 
              #                            column(7, plotOutput("sum_hist")))),
            tabPanel("Edit Module", fluidRow(br(),
                                             column(4, plotOutput('img_result_em')),
                                             column(3, uiOutput('edit_module')))),
             tabPanel("FIJI Comparison", fluidRow(br(),
                                          column(4, displayOutput('labeled_img_ex')),
                                          column(3, displayOutput('fiji_results')))),
             tabPanel("Table", fluidRow(
                                        br(),
                                        column(5, br(), displayOutput("labeled_img")), 
                                        column(5,DTOutput("mpi_table")))),
             tabPanel("Log", verbatimTextOutput("results_log"))),
    ))
  })
}

get_edit_module <- function(res)
{
  renderUI(
    { 
      fluidRow(
      column(5,checkboxGroupInput("cell_selections",
                         label = "", inline=TRUE,
                         choices = res$df[["CellID"]])),
      column(2, actionButton("remove_cells", "Remove"), actionButton("finish", "Finished Editing (all cells!)")))
    })
}

get_fiji_result <- function(res)
{
  renderDisplay({
    filename = paste0("Demo/Fiji Results/",res$filename, "_fiji.png")
    img_fiji <- readImage(file.path(filename))
    display(img_fiji, method = 'browser')
  })
}

get_edit_plot <- function(res)
{
  renderPlot(
    {
      plot( res$channels$gfp)
      sapply(res$mem_pts, function(x){points(x, type = 'l', col="white")})
      sapply(res$vac_pts, function(x){points(x, type = 'l', col=c("yellow","yellow"))})
      text(x = res$df[,'pm_center_x'],
           y = res$df[, 'pm_center_y'],
           labels = res$df[,'CellID'], 
           col = "red", 
           pos = c(2,3), #(2,3) = to the left of and above
           vfont = c("sans serif", "bold"))
      
    })
  
}
get_image_plot <-function(res, sel, chan)
{
  renderPlot(
    {
     
      if(is.null(res)) 
      {
        return()
      }
      par(mar=rep(0,4))
      if(chan == "cmac")
        plot( res$channels$cmac)
      else if(chan == "gfp")
        plot( res$channels$gfp)
      else if(chan == "dic")
        plot( res$channels$dic)
      else if(chan == "gray_gfp")
        plot(res$img_gray[,,gfp_channel])
      else
        plot(res$channels$gfp)

      par(new=TRUE)
      if(!is.null(sel))
      {
       
        if("mem_select" %in% sel)
          sapply(res$mem_pts, function(x){points(x, type = 'l', col="white")})
        if("vac_select" %in% sel)
          sapply(res$vac_pts, function(x){points(x, type = 'l', col="yellow")})
        if("pm_label_select" %in% sel)
          text(x = res$df[,'pm_center_x'],
               y = res$df[, 'pm_center_y'],
               labels = res$df[,'CellID'], 
               col = "red", 
               pos = c(2,3), #(2,3) = to the left of and above
               vfont = c("sans serif", "bold"))
        if("vac_label_select" %in% sel)
          text(x = res$df[,'pm_center_x'],
               y = res$df[, 'pm_center_y'],
               labels = res$df[,'vacuoles'], 
               col = "orange", 
               pos = c(3,4), # (3,4) = to the right of and above
               vfont = c("sans serif", "bold"))
        
        if("rem_select" %in% sel)
          sapply(res$removed_pts, function(x){points(x, type = 'l', col="red")})
      }
      par(new = TRUE)
    })
}

get_fiji_comparison <- function(res)
{
  renderDisplay({
    t <- paintObjects(res$membranes, tgt = res$channels$gfp, col = c('white','white'))
    tt <- paintObjects(res$vacuoles, tgt = t, col = 'yellow')
    display(tt, method = "browser")
  })
}

get_final_labeled <- function(res)
{
  renderDisplay({
   filename = paste0("Demo/Output/",res$filename, "_final_results.tiff")
   img_final <- readImage(file.path(filename))
   display(img_final, method = 'browser')
  })
}

get_sum_text <- function(res)
{
  renderUI({
    HTML(paste0("Number of cell membranes: ", 
           length(res$FCM[,1]), 
           br(),"Average MPI over all cells: ", 
           format(mean(res$FCM[,"membrane.a.b.mean"]), nsmall = 3),
           br(),"Number budding cells removed: "), length(res$removed_puts))
    
    
  })
}
get_hist <- function(res)
{
  renderPlot({
  df <- as.data.frame(res$FCM[,"membrane.a.b.mean"])
  names(df)[1] <- "mpi"
  ggplot(data = df, mapping = aes(x=mpi))+ 
    geom_histogram(binwidth = .05, 
                   fill="aquamarine4",
                   alpha = 0.5,
                   color = "azure4")+
    labs(title = paste("Mean Pixel Intensity: ",res$filename),
         y = "Frequency",
         x = "Mean Pixel Intensity")+
    theme(plot.title = element_text(hjust = 0.5))
  })
}

