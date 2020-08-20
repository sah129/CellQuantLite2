
pipeline_options <- function(image_file, gui, progress, gfp_chan, cmac_chan, dic_chan, cutoff)
{
  
  
  cnum = list(
    cmac_channel = as.numeric(cmac_chan),
    gfp_channel = as.numeric(gfp_chan),
    dic_channel = as.numeric(dic_chan)
    
  )
  


  

  results_gfp = list()
  results_dic = list()

  factors = list(1,2,4,8,16)

  #if(gui)
   # progress$inc(1/nrow(imageset), detail = paste0(imageset[row,"filename"], "(", row, "/",nrow(imageset),")" ))
  
  channels <- read_in_channels_single(image_file)
  img_gray <- convert_to_grayscale(channels)


  
  print("##############GFP MEMBRANE DETECTION ALGORITHM###############")
  for(factor in factors)
  {
    membranes <- detect_membranes_new(img_gray, channels, factor, img_gray[,,cnum$gfp_channel], as.numeric(cutoff), cnum)
    vacuoles <- find_vacuoles(membranes, img_gray, channels, cnum)
    res <- exclude_and_bind(membranes, vacuoles)
    final<-tidy_up(membranes,vacuoles,res)

      
    results_gfp[[factor]] <- final
    
  
  }

  print("#################DIC MEMBRANE DETECTION ALGORITHM####################")
  for(factor in factors)
  {
    membranes <- detect_membranes_new(img_gray, channels, factor, img_gray[,,cnum$dic_channel], as.numeric(cutoff), cnum)
    vacuoles <- find_vacuoles(membranes, img_gray, channels, cnum)
    res <- exclude_and_bind(membranes, vacuoles)
    final<-tidy_up(membranes,vacuoles,res)
    

    
    
    results_dic[[factor]] <- final
    
    
  }

    
  
  message("End of pipeline_options")
  
  result <<- list(gfp = results_gfp, dic = results_dic, channels = channels)
  list(gfp = results_gfp, dic = results_dic, channels = channels)
}


