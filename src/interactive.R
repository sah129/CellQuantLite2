
# Remove cell fragments, empty cells and associated vacuoles without 
# populating the data frame. 
# mems: membranes object
# vacs: vacuole object
# res: object storing all results to be modified here and returned to client
# renum: option to renumber objects,  if desired
get_first_pass <- function(mems,vacs,res, renum)
{
  final_pm_img <- get_final_pm_img(mems, res, renum)
  final_vac_img <- get_final_vac_img(vacs, res, renum)
  list(membranes = final_pm_img, vacuoles = final_vac_img, df = drop_na(res$df))
}

# Remove cells selected by user in interactive interface.  Return updated
# res object.
remove_cells_interactive <- function(res, i, to_remove)
{
  vacs_tr <- res[[i]]$df[ (res[[i]]$df$CellID %in% to_remove), "vacuoles" ]
  res[[i]]$vacuoles <- rmObjects(res[[i]]$vacuoles, vacs_tr, reenumerate = FALSE)
  res[[i]]$membranes <- rmObjects(res[[i]]$membranes, to_remove, reenumerate = FALSE)
  res[[i]]$df <- res[[i]]$df[ !(res[[i]]$df$CellID %in% to_remove), ]
  res[[i]]$mem_pts <- ocontour(res[[i]]$membranes)
  res[[i]]$vac_pts <- ocontour(res[[i]]$vacuoles)
  return(res)
}

# Interactive analog of tidy_up.  When user is done manually pruning results, 
# format the dataframe and write the output. 
finish_up <- function(res)
{

  for(r in res)
  {
    r$df <- renumerate_df(r$df)
    tiff(filename = paste0("FinalOutput/",r$filename, "_final_results.tiff"))
    get_display_img(df = r$df,
                    membranes = r$membranes, 
                    col_membranes = 'white', 
                    vacuoles = r$vacuoles, 
                    col_vacuoles ='yellow', 
                    removed = r$removed,
                    closed_vacuoles = TRUE, 
                    img = channel(r$channels$gfp, "asgreen"), 
                    showRemoved = TRUE, 
                    showMemLabels = TRUE, 
                    showVacLabels = FALSE)
    dev.off()
  }
  sort_data(res)
}

