
source("src/functions.R")


pipeline <- function(datasetpath, testing, gui, progress)
{
  if(testing)
  {
    return(readRDS("Saved Results/demo_five_results.rds"))
  }
  
  imageset <- read_in_imageset(datasetpath)
  results = list()
  for( row in 1:nrow(imageset))
  {
    if(gui)
      progress$inc(1/nrow(imageset), detail = paste0(imageset[row,"file"], "(", row, "/",nrow(imageset),")" ))
    
    
    
    channels <- read_in_channels(imageset[row,], datasetpath)
    imggray <- convert_to_grayscale(channels)
    # img <- scale_intensities(imggray)  #problem -- refs/grayscale
    cells <- detect_cells(imggray)
    
    
    res <- find_vacuoles(cells, imggray)
    
    
    final <- display_output(res$membranes, channels$gfp, res$vacuoles, labeled = FALSE)
    
    labels <- roi_labels_demo(res$membranes, res$FCM)
    
    writeImage(x = final,
               type = "tiff",
               bits.per.sample = 16,
               compression = "none",
               quality = 100,
               files = paste0("Output/",imageset[row, "file"], "_img.tiff"))
    
    
    write_to_csv(res$FCM, imageset[row, "file"])
    
    
    results[[row]] <- list(channels = channels, 
                           #final = final,
                           labels = labels,
                           filename = imageset[row, "file"],
                           membranes = res$membranes,
                           vacuoles = res$vacuoles,
                           mem_pts = ocontour(res$membranes),
                           vac_pts = ocontour(res$vacuoles),
                           removed_puts = ocontour(cells$removed),
                           FCCI = res$FCCI,
                           #  FCCB = res$FCCB,
                           FV = res$FV,
                           FCM = res$FCM)
    
    gc()
  }
  message("End of main")
  # if(gui)
  return(results)
}

#hres <- pipeline("Datasets/Diploid/Demo_Diploid_Single", testing=FALSE, gui=FALSE, progress=NULL)

#saveRDS(hres, "Saved Results/demo_diploid_single.rds")


#res <- readRDS(file = 'results-04-03-2020.rds')


=======
  
  source("src/functions.R")


pipeline <- function(datasetpath, testing, gui, progress)
{
  if(testing)
  {
    return(readRDS("Saved Results/demo_five_results.rds"))
  }
  
  imageset <- read_in_imageset(datasetpath)
  results = list()
  for( row in 1:nrow(imageset))
  {
    if(gui)
      progress$inc(1/nrow(imageset), detail = paste0(imageset[row,"file"], "(", row, "/",nrow(imageset),")" ))
    
    
    
    channels <- read_in_channels(imageset[row,], datasetpath)
    imggray <- convert_to_grayscale(channels)
    # img <- scale_intensities(imggray)  #problem -- refs/grayscale
    cells <- detect_cells(imggray)
    
    
    res <- find_vacuoles(cells, imggray)
    
    
    final <- display_output(res$membranes, channels$gfp, res$vacuoles, labeled = FALSE)
    
    labels <- roi_labels_demo(res$membranes, res$FCM)
    
    writeImage(x = final,
               type = "tiff",
               bits.per.sample = 16,
               compression = "none",
               quality = 100,
               files = paste0("Output/",imageset[row, "file"], "_img.tiff"))
    
    
    write_to_csv(res$FCM, imageset[row, "file"])
    
    
    results[[row]] <- list(channels = channels, 
                           #final = final,
                           labels = labels,
                           filename = imageset[row, "file"],
                           membranes = res$membranes,
                           vacuoles = res$vacuoles,
                           mem_pts = ocontour(res$membranes),
                           vac_pts = ocontour(res$vacuoles),
                           removed_puts = ocontour(cells$removed),
                           FCCI = res$FCCI,
                           #  FCCB = res$FCCB,
                           FV = res$FV,
                           FCM = res$FCM)
    
    gc()
  }
  message("End of main")
  # if(gui)
  return(results)
}

#hres <- pipeline("Datasets/Diploid/Demo_Diploid_Single", testing=FALSE, gui=FALSE, progress=NULL)

#saveRDS(hres, "Saved Results/demo_diploid_single.rds")


#res <- readRDS(file = 'results-04-03-2020.rds')


