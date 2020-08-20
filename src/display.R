
# Wrapper function to perform final pass on membranes, vacuoles, and 
# update the data frame.
tidy_up <- function(membranes,vacuoles,res)
{
  final_pm_img <- get_final_pm_img(membranes, res)
  final_vac_img <- get_final_vac_img(vacuoles, res)
  final_df <- renumerate_df(res$df)
  
  list(membranes = final_pm_img, vacuoles = final_vac_img, df = final_df)
}

# Remove fragments and empty cells while taking care to preserve numbering.
get_final_pm_img <- function(mems, res)
{
  r = mems$membranes
  if(!is.null(res$empty_cells))
    r <- rmObjects(mems$membranes, res$empty_cells, reenumerate = FALSE)
  if(!is.null(res$fragments))
    r <- rmObjects(r, res$fragments, reenumerate = FALSE)
  r<-reenumerate(r)
  return(r)
}

# Format the dataframe to conform to O'Donnell lab standards.  Combined vacuoles
# will be assigned label 2x+1, where x=ID of the parent cell.
renumerate_df <- function(df)
{
  if(any(is.na(df['CellID'])))
  {
    df = drop_na(df)
    message('renumerating cell ids')
    if(nrow(df) > 0)
    {
      df['CellID'] = seq(1:nrow(df))
    }
  }
  message('reenumerating vacuoles')
  if(nrow(df) > 0)
  {
    df['vacuoles'] = seq((nrow(df)+1),2*nrow(df))
  }
  return(df)
  
}

# Remove vacuoles belonging to excluded cells and reformat the dataframe.
get_final_vac_img <- function(vacs, res)
{
  vac_df <- drop_na(res$df["vacuoles"])
  if(nrow(vac_df) == 0)
  {
    return(NULL)
  }
  l <- length(table(vacs$vacuoles)[-1])
  vac_list <- vector('list', l)
  i=1
  for(v in vac_df)
  {
    cv <- strsplit(v, ',')
    for(vv in cv)
    {
      for(vvv in vv) 
      {
        vac_list[i]=as.numeric(str_trim(vvv))
        i = i + 1
      }
    }
    vac_list <- unlist(vac_list)
    removedvacs <- which(!(seq(1:length(table(vacs$vacuoles)[-1])) %in% vac_list))
  }
  res_vacs <- rmObjects(vacs$vacuoles, removedvacs, reenumerate = TRUE)
  return(res_vacs)

}




# Wrapper function to generate display image.
# df: dataframe to use
# membranes: membrane objects
# col_membranes: color to paint membranes on image
# vacuoles:  vacuole objects
# col_vacuoles:  color to paint vacuoles on image
# removed:  removed membrane objects
# closed_vacuoles: set to TRUE to paint border+fill of vacuoles
# img: background image to draw on
# showRemoved:  show removed membranes
# showMemLabels:  show membrane labels
# showVacLabels:  show vacuole labels

get_display_img <- function(df,membranes, col_membranes, vacuoles, col_vacuoles, removed,closed_vacuoles, img, showRemoved, showMemLabels, showVacLabels)
{
  if(nrow(df)==0) # nothing detected
  {
    plot(img)
  }
  else
  {
    res_imgA <- paintObjects(membranes, tgt = img, col = c(col_membranes, col_membranes))
    vac_col <- col_vacuoles
    if(closed_vacuoles)
      vac_col <- c(col_vacuoles,col_vacuoles)
    res_img <- paintObjects(vacuoles, tgt = res_imgA, col = vac_col)
    if(showRemoved)
      res_img <- paintObjects(removed, tgt = res_img, col = c('red','red'))
    plot(res_img)
    if(showMemLabels)
    {
      text(x = df[,'pm_center_x'],
           y = df[, 'pm_center_y'],
           labels = df[,'CellID'], 
           col = "red", 
           pos = c(2,3), #(2,3) = to the left of and above
           vfont = c("sans serif", "bold"))
    }
    if(showVacLabels)
    {
      text(x = df[,'pm_center_x'],
           y = df[, 'pm_center_y'],
           labels = df[,'vacuoles'], 
           col = "orange", 
           pos = c(3,4), # (3,4) = to the right of and above
           vfont = c("sans serif", "bold"))
    }
  }
}

get_display_helper <- function(final, channels)
{
  get_display_img(df = final$df,
                  membranes = final$membranes, 
                  col_membranes = 'white', 
                  vacuoles = final$vacuoles, 
                  col_vacuoles ='yellow', 
                  removed = membranes$removed,
                  closed_vacuoles = TRUE, 
                  img = channels$gfp, 
                  showRemoved = FALSE, 
                  showMemLabels = TRUE, 
                  showVacLabels = FALSE)
}


