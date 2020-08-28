get_options_output <- function(res)
{
  if(is.null(res))
    return(NULL)
  result<<-res
  renderPlot({
    
    get_display_helper(res$gfp[[1]], res$channels)
  })
}

get_test <- function(res, channels, i) 
{
  title <- switch(as.character(i),
                 "1" = "1",
                 "2" = "2",
                 "4" = "3",
                 "8" = "4",
                 "16" = "5")
  renderPlot({
        if(length(res) > 0 )
          get_display_helper(res[[i]], channels, title = title)
        else
        {
          img <- readImage(file.path('src/res/empty_dic.png'))
          plot(img)
          title(main = title, line = 0)
        }
  })
  
  }