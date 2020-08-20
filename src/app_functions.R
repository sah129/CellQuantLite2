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
  renderPlot({
        get_display_helper(res[[i]], channels, title = i)
  })
  
  }