setwd("C:/Users/sarah/Desktop/Capstone/Main pipeline")

source("functions.R")

cmac_channel = 1
gfp_channel = 2
dic_channel = 3



pipeline <- function(channels)
{

  img <- convert_to_grayscale(channels)
  
  img <- scale_intensities(img)
  masks <- create_masks(img)
  segments <- create_cell_segments(img)
  edited <- exclude_small_rois(masks$FV, masks, segments)
  
  final <- img_vac_seg(channels$cmac, edited$segments, edited$masks)
 
  list(img = img, 
       masks_Orig = masks, 
       masks_edited = edited$masks,
       segments_orig = segments,
       segments_edited = edited$segments, 
       final = final, 
       FV = masks$FV)
}


pipeline_supp <- function()
{
  print_roi_labels(final,masks)
  
  p <- examine_image_segment(14, img, masks)
  display(p)
  
  irregular_cells = c(45, 58)
  irregular_cells_arr <- lapply(irregular_cells, function(x) examine_image_segment(x, img, masks))
  c <- combine(irregular_cells_arr)
  display(c)
  
  
  FV = computeFeatures(masks$vacuoles,img[,,cmac_channel], xname = "vacuoles")
  
  
  irreg_stats <- lapply(irregular_cells,function(x) vac_stats(x, FV))
  
  mean(FV[,"vacuoles.0.s.area"])
  
  n <- exclude_small_roi(FV, .5, masks)
  display(n)
  
  display(img_vac_seg(channels$cmac, segments, masks$vacuoles))
}
