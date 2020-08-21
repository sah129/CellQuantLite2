
pipeline_options <- function(image_file, gui, progress, gfp_chan, cmac_chan, dic_chan, cutoff)
{
  
  
  cnum = list(
    cmac_channel = as.numeric(cmac_chan),
    gfp_channel = as.numeric(gfp_chan),
    dic_channel = as.numeric(dic_chan)
    
  )
  
  i = 1

  

  results_gfp = list()
  results_dic = list()

  factors = list(1,2,4,8,16)


  
  channels <- read_in_channels_single(image_file)
  img_gray <- convert_to_grayscale(channels)

  
  
  print("##############GFP MEMBRANE DETECTION ALGORITHM###############")
  for(factor in factors)
  {
    out<- tryCatch({
    
      progress$inc(1/10, detail = paste0("GFP DETECTION (iteration ", i, ")" ))
      
      membranes <- detect_membranes_new(img_gray, channels, factor, img_gray[,,cnum$gfp_channel], as.numeric(cutoff), cnum)
      vacuoles <- find_vacuoles(membranes, img_gray, channels, cnum)
      res <- exclude_and_bind(membranes, vacuoles)
      final<-tidy_up(membranes,vacuoles,res)
  
      i = i + 1
      results_gfp[[factor]] <- final

    },
    error = function(cond)
    {
      print(paste0("Error analyzing ", image_file[1,"name"], " in GFP detection for factor=", factor))
      results_gfp[[factor]] <- NULL
      progress$inc(1/10, detail = paste0("GFP DETECTION (iteration ", i, ")" ))
      i = i + 1
    })
    
  
  }

  i = 1
  print("#################DIC MEMBRANE DETECTION ALGORITHM####################")
  for(factor in factors)
  {
    out<- tryCatch({
      
      progress$inc(1/10, detail = paste0("DIC DETECTION (iteration ", i, ")" ))
      
      membranes <- detect_membranes_new(img_gray, channels, factor, img_gray[,,cnum$dic_channel], as.numeric(cutoff), cnum)
      vacuoles <- find_vacuoles(membranes, img_gray, channels, cnum)
      res <- exclude_and_bind(membranes, vacuoles)
      final<-tidy_up(membranes,vacuoles,res)
      
  
      i = i + 1
      
      results_dic[[factor]] <- final
    },
    error = function(cond)
    {
      print(paste0("Error analyzing ", image_file[1,"name"], " in DIC detection for factor=", factor))
      results_dic[[factor]] <- NULL
      progress$inc(1/10, detail = paste0("DIC DETECTION (iteration ", i, ")" ))
      i = i + 1
    })
    
    
  }

 
    
  
  message("End of pipeline_options")
  
  result <<- list(gfp = results_gfp, dic = results_dic, channels = channels)
  list(gfp = results_gfp, dic = results_dic, channels = channels)
}


