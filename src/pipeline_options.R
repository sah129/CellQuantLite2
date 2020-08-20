
pipeline_options <- function(testpath, gui, progress, cutoff)
{

  #img <- readImage(image_file[1, "filepath"])

  results_gfp = list()
  results_dic = list()
  
  factors = list(1,2,4,8,16)

  #if(gui)
   # progress$inc(1/nrow(imageset), detail = paste0(imageset[row,"filename"], "(", row, "/",nrow(imageset),")" ))
  
  channels <- read_in_channels_single(testpath)
  img_gray <- convert_to_grayscale(channels)
  
  print("##############GFP MEMBRANE DETECTION ALGORITHM###############")
  for(factor in factors)
  {
    membranes <- detect_membranes_new(img_gray, channels, factor, img_gray[,,gfp_channel], cutoff)
    vacuoles <- find_vacuoles(membranes, img_gray, channels)
    res <- exclude_and_bind(membranes, vacuoles)
    final<-tidy_up(membranes,vacuoles,res)

      
    results_gfp[[factor]] <- final
    
  
  }

  print("#################DIC MEMBRANE DETECTION ALGORITHM####################")
  for(factor in factors)
  {
    membranes <- detect_membranes_new(img_gray, channels, factor, img_gray[,,dic_channel], cutoff)
    vacuoles <- find_vacuoles(membranes, img_gray, channels)
    res <- exclude_and_bind(membranes, vacuoles)
    final<-tidy_up(membranes,vacuoles,res)
    

    
    
    results_dic[[factor]] <- final
    
    
  }

    
  
  message("End of pipeline_options")
  list(gfp = results_gfp, dic = results_dic, channels = channels)
}


