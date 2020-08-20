source("src/functions.R")

# Main Driver function.  
# datasetpath: path to dataset, see io.R for file/directory structure
# testing: If set set to true, a precomputed .rds file is loaded and
# returned.
# gui:  set to true when using the GUI.  See app.R.
# progress:  holds information for progress updates when using the GUI.
# interactive:  set to manually prune results in GUI
pipeline <- function(image_files, testing, gui, progress, interactive, factor, chan, cutoff)
{
  # defaults for image frames
  cmac_channel = 1
  gfp_channel = 2
  dic_channel = 3
  
  
  
  if(testing)
    return(readRDS("Demo/Saved Results/presentation_results.rds"))
  
  imageset <- read_in_imageset_interactive(image_files)
  results = list()
  for( row in 1:nrow(imageset))
  {
    if(gui)
      progress$inc(1/nrow(imageset), detail = paste0(imageset[row,"filename"], "(", row, "/",nrow(imageset),")" ))
    
    channels <- read_in_channels_interactive(imageset[row,])
    img_gray <- convert_to_grayscale(channels)
   
    membranes <- detect_membranes_new(img_gray, channels, factor, img_gray[,,chan], cutoff)
    vacuoles <- find_vacuoles(membranes, img_gray, channels)
    res <- exclude_and_bind(membranes, vacuoles)
    
    if(interactive)
      final <- get_first_pass(membranes,vacuoles,res)
    else
      final<-tidy_up(membranes,vacuoles,res)
    
    if(nrow(final$df)==0)
    {
      mem_pts = NULL 
      vac_pts = NULL
    }
    else
    {
      mem_pts = ocontour(final$membranes)
      vac_pts = ocontour(final$vacuoles)
    }
    if(length(membranes$removed) > 0)
      removed_pts = ocontour(membranes$removed)
    else
      removed_pts = NULL
    
    tiff(filename = paste0("FinalOutput/Images/",imageset[row, "filename"], "_final_results.tiff"))
    
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
    dev.off()
    
    write.csv(final$df, paste0("FinalOutput/Individual Spreadsheets/",imageset[row, "filename"], '_results.csv'), row.names=FALSE)
    results[[row]] <- list(df = final$df,
                           filename = imageset[row, "filename"])
                        
    
    
    
    
  }
  message("End of main")
  return(results)
}


