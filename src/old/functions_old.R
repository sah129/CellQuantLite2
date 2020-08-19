<<<<<<< HEAD

library("EBImage")
library("genefilter")

cmac_channel = 1
gfp_channel = 2
dic_channel = 3

read_in_channels <- function(dic,gfp,cmac)
{

  dic_path=dic
  gfp_path=gfp
  cmac_path=cmac
  
  


  dic = readImage(file.path(dic_path))
  gfp = readImage(file.path(gfp_path))
  cmac = readImage(file.path(cmac_path))
  
  list(cmac = cmac, gfp = gfp, dic = dic)

}
# Function to remove cells that are smaller than 1 standard deviation from the mean size. +-2 standard deviations allows some noise in,
# so far 1 sd is not excluding too many.  Returns cell objects as well as computed feateres.
remove_small_cells <- function(segs, cells)
{
  
  areas <- segs$FC[,"s.area"]
  cutoff <- mean(areas) - sd(areas)
  small <- which(areas < cutoff)
  print(paste0("Cutoff point for cell area (in pixels): ", cutoff))
  print(paste0("Number of cells removed for being below cutoff: ", length(small)))
  segs$segs <- rmObjects(segs$segs, small)
  
  segs$segs <- bwlabel(segs$segs)
  segs$FC = computeFeatures(segs$segs, ref = cells[,,gfp_channel], xname = "seg")
  return(segs)
}


# Wrapper function to find cell objects
find_cell_bodies <- function(cells)
{
  print("################CELLS###############")
  
  segs <- first_pass_new(cells)
  #  segs <- remove_small_cells(segs, cells)
  #segs <- remove_odd_shapes(segs, cells)
  print(paste0("Number of cells after all removals: ", length(table(segs$segs)) ))
  return(segs)
  
}
first_pass <- function(cells)
{
  cb = gblur(cells[,,gfp_channel], sigma = 2)
  ct = thresh(cb, offset = 0.005)
  
  cm = bwlabel(ct)
  FS = computeFeatures.shape(cm)
  sel <- which(FS[,"s.area"] < 500)
  ce = rmObjects(cm, sel)
  
  
  
  cell_segments_full =  fillHull(cell_segments) - cell_segments
  segs <- fillHull(cell_segments_full)
  segs <- bwlabel(segs)
  FC = computeFeatures.shape(segs)
  print(paste0("Number of cells detected on first pass: ", length(table(segs))))
  list(segs = segs, FC = FC)
}
convert_to_grayscale<- function(img)
{
  gfp_gray = channel(img$gfp, 'gray')
  cmac_gray = channel(img$cmac, 'gray')
  
  # RBioformats has .nd2 capability
  img = combine(cmac_gray, gfp_gray, img$dic)
  return(img)
}

scale_intensities <- function(img)
{
  img[,,cmac_channel] = normalize(img[,,cmac_channel])
  img[,,gfp_channel] = normalize(img[,,gfp_channel])
 
  return(img)
}

create_masks <- function(img)
{

  vmask = thresh(img[,,cmac_channel],w=30,h=30,offset = sd(img[,,cmac_channel]))
  vmask = opening(vmask, makeBrush(5,shape="disc"))
  vmask = fillHull(vmask)
  vmask = bwlabel(vmask)
  
  FV = computeFeatures(vmask,img[,,cmac_channel], xname = "vacuoles")
  
  
  ctmask = opening(img[,,gfp_channel] > 0.1, makeBrush(15,shape="disc"))
  cmask = propagate(img[,,gfp_channel], seeds = vmask, mask = ctmask )
  
  list(vacuoles = vmask, cytoplasm = ctmask, cells = cmask, FV = FV)

}


# First pass cell detection.  Objects are detected via blurring, simple threshholding, and noise removal.
# The last 4 lines allow us to return only cells with a complete membrane detectable

detect_cells_old <-function(img)
{
  message("########################CELLS########################")
  
  cells_blurred = gblur(img[,,gfp_channel], sigma = 2)
  ct = thresh(cells_blurred, offset = 0.005)
  
  cm = bwlabel(ct)
  FS = computeFeatures.shape(cm)
  sel <- which(FS[,"s.area"] < 500)
  membranes <-rmObjects(cm, sel)
  
  # membranes <- bwlabel(membranes)
  
  
  cell_bodies <- fillHull(membranes)
  cell_bodies <- bwlabel(cell_bodies)
  
  
  FCCB <- computeFeatures.shape(cell_bodies)
  
  
  
  mask <- list(cell_bodies = cell_bodies, FCCB = FCCB)
  
  mask <- remove_odd_shapes(mask, img)
  
  cell_bodies <- mask$cell_bodies
  FCCB <- mask$FCCB
  
  res <- remove_small_cells(cell_bodies, FCCB, img)
  
  cell_bodies <- res$cell_bodies
  FCCB <- res$FCCB
  
  
  cell_inner <- fillHull(membranes) - membranes
  cell_inner <- bwlabel(cell_inner)
  cell_inner <- fillHull(cell_inner)
  
  cell_inner <- cell_inner*cell_bodies
  cell_inner <- bwlabel(cell_inner)
  
  
  FCI <- computeFeatures.shape(cell_inner)
  sel <- which(FCI[,"s.area"] < 500)
  cell_inner <-rmObjects(cell_inner, sel)
  
  #  cell_inner <- bwlabel(cell_inner)
  
  
  membranes <- cell_bodies - cell_inner
  # membranes <- bwlabel(membranes)
  
  #FCM <- computeFeatures(membranes, ref = img[,,gfp_channel], xname = "membrane")
  
  
  
  FCCI <- computeFeatures.shape(cell_inner)
  message(paste0("Number of cells detected on first pass: ", length(table(cell_bodies))))
  message(paste0("Number of inner detected on first pass: ", length(table(cell_inner))))
  
  list(
    cell_inner = cell_inner,
    removed = mask$removed,
    FCCI = FCCI, 
    membranes = membranes, 
    #FCM = FCM, 
    cell_bodies = cell_bodies, 
    FCCB = FCCB)
  
}


# Function to compute circularity of a cell object
circularity <- function(seg_num, FC) 
{
  4 * pi * ((FC[seg_num,"s.area"]) / (FC[seg_num,"s.perimeter"]^2))
} 

# Removes odd shapes (those with circularity less than 1.)  The circularity cuttoff can be played with,
# 1 generally will be enough to remove budding cells.  
remove_odd_shapes <- function(mask, img)
{
  circ <- lapply(1:length(table(mask$cell_bodies)) -1, circularity, FC = mask$FCCB)
  
  circ<-unlist(circ)
  odd_shapes = which(circ < 1) #budding cells
  message(paste0("Number of cells removed for being oddly shaped: ", length(odd_shapes)))
  
  removed_ind = which(!(seq(1, length(table(mask$cell_bodies))) %in% odd_shapes))
  removed = rmObjects(mask$cell_bodies, removed_ind)
  
  cell_bodies <- rmObjects(mask$cell_bodies, odd_shapes)
  cell_bodies <- bwlabel(cell_bodies)
  FCCB <- computeFeatures.shape(cell_bodies)#(cell_bodies, ref = img[,,gfp_channel], xname = "seg")
  
  list(cell_bodies = cell_bodies, FCCB = FCCB, removed = removed)
  
  
  
  
  
}


remove_small_cells <- function(cell_bodies, FCCB, img)
{
  
  areas <- FCCB[,"s.area"]
  cutoff <- mean(areas) - sd(areas)
  small <- which(areas < cutoff)
  message(paste0("Cutoff point for cell area (in pixels): ", format(cutoff, nsmall = 3)))
  message(paste0("Number of cells removed for being below cutoff: ", length(small)))
  cell_bodies <- rmObjects(cell_bodies, small)
  
  cell_bodies <- bwlabel(cell_bodies)
  FCCB = computeFeatures.shape(cell_bodies)
  list(cell_bodies = cell_bodies, FCCB = FCCB)
  
}



create_cell_segments<- function(cells)
{
  cb = gblur(cells[,,gfp_channel], sigma = 2)
  ct = thresh(cb, offset = 0.005)
  
  cm = bwlabel(ct)
  FS = computeFeatures.shape(cm)
  sel <- which(FS[,"s.area"] < 500)
  ce = rmObjects(cm, sel)
  
  cell_segments <- thresh(ce)
  #cell_segments <- stackObjects(cell_segments,gfp_gray)\
  cell_segments <- bwlabel(cell_segments)
  
  cell_segments_filled <- fillHull(cell_segments)

  return(cell_segments)

}


# If a cell has more than 1 vacuole detected, keeps the biggest one and removes all others.  Also removes cells with
# no detected vacuoles.  After all removals features are recomputed.  This is by far the most computationally expensive function, and can be improved by 
# refining the parameters used to detect the vacuoles.  Still a work in progress.  
remove_multiples <- function(cell_info, vmask, img)
{
  
  l = length(table(cell_info$cell_inner))
  l = l-1
  
  
  maxes <- vector("list", l)
  empty_cells <- vector("list", l)
  
  for(i in seq(1:l))
  {
    
    cseg = cell_info$cell_inner == i
    vcount <- table(cseg*vmask)
    vcount <- vcount[-1]
    #message(vcount)
    if( length(vcount) > 0)
    {
      
      v <- as.numeric(names(as.list(vcount)))
      # message(v)
      vcount <- unname(vcount)
      maxes[i] <- v[which(vcount == max(vcount))]
    }
    else
    {
      # message(paste0(i, " ZERO VCOUNT"))
      empty_cells[i] = i
    }
  }
  
  empty_cells <- unlist(empty_cells)
  
  
  message(paste0("Number of empty cells: ", length(empty_cells)))
  maxes <- unlist(maxes)
  extras <- which(!(seq(1:length(table(vmask)[-1])) %in% maxes))
  
  
  
  vmask <- rmObjects(vmask, extras)
  vmask <- bwlabel(vmask)
  
  
  # Not going to recompute features here for cell bodies/cell inner.  Can add that if needed but as of now it's not necessary and it will slow things down
  # recomputation of membrane features is needed, however
  
  
  
  #remove empty cells
  cell_inner <- rmObjects(cell_info$cell_inner, empty_cells)
  cell_inner <- bwlabel(cell_inner)
  #remove empty cells
  cell_bodies <- rmObjects(cell_info$cell_bodies, empty_cells)
  cell_bodies <- bwlabel(cell_bodies)
  FCCI <- computeFeatures(cell_inner, ref = img[,,gfp_channel])
  FCCB <- computeFeatures(cell_bodies, ref = img[,,gfp_channel])
  
  #remove empty membranes
  #membranes <- rmObjects(cell_info$membranes, empty_cells)
  # membranes <- bwlabel(membranes)
  membranes = cell_bodies - cell_inner
  
  FCM <- computeFeatures(membranes, ref = img[,,gfp_channel], xname = "membrane")
  FV <- computeFeatures(vmask, ref= img[,,gfp_channel], xname = "vac")
  message(paste0("Number of cells after empty cells removed: ", length(table(cell_bodies))))
  
  list(vacuoles = vmask, 
       membranes = membranes,
       FCM = FCM,
       #  FV = FV,
       FCCB = FCCB,
       cell_inner = cell_inner,
       cell_bodies = cell_bodies,
       FCCI = FCCI)
}




img_vac_seg <- function(cmac, cell_segments, masks)
{
  c = paintObjects(cell_segments, tgt = cmac, col = 'brown', closed = TRUE)
  output = paintObjects(masks$vacuoles, c, col='yellow', thick = TRUE)
  return(output)
}



examine_image_segment <- function(segment_number,img,masks)
{

    
  segment = masks$cells == segment_number
  
  gfp_mask = 1- img[,,gfp_channel]
  gfp_mask = gfp_mask*segment
  
  
  vseg = masks$vacuoles == segment_number
  cell_segment = gfp_mask*gfp_gray*20
  
  p <- paintObjects(vseg, tgt = cell_segment, col="red")
  
  return(p)
}

vac_stats <- function(segment_number, FV)
{
  list(majoraxis = FV[segment_number,"vacuoles.0.m.majoraxis"], radiussd = 
  FV[segment_number,"vacuoles.0.s.radius.sd"])
  
}
print_roi_labels <- function(channels,cell_segments, masks)
{
  
  img <- img_vac_seg(channels$cmac,cell_segments, masks$vacuoles)
  labels <- c(1:length(table(masks$vacuoles) - 1))
  center_pts = get_center_points(masks$vacuoles)
  display(img, method = "raster")
  text(center_pts, label = labels, col = "white")
}


get_center_points <- function(vmask)
  {
  vac_centers = computeFeatures.moment(vmask)
  center_pts = vac_centers[,c("m.cx","m.cy")] 
  return(center_pts)
}

exclude_small_rois <- function(FV, masks, cell_segments)
{
  too_small <- which(FV[,"vacuoles.0.s.area"] < 200)
  masks$vacuoles <- rmObjects(masks$vacuoles, too_small)
  masks$cells <- rmObjects(masks$cells, too_small)
  cell_segments = cell_segments*masks$cell
  list(masks = masks, segments = cell_segments)
  
}



roi_labels <- function(vmask)
{
  FV = computeFeatures.moment(vmask)
  vac_centers =FV[,c("m.cx", "m.cy")]
  
  
  
  labels <- c(1:length(table(vmask)) )
  display(vmask, method = "raster")
  text(vac_centers, label = labels, col = "red")
  
  
}


write_to_csv <-function(FCM, filename)
{
  
  mean_pixel_intensity <- FCM[,"membrane.a.b.mean"]
  ids <- as.numeric(names(mean_pixel_intensity))
  
  
  
  mean_pixel_intensity <- as.data.frame( mean_pixel_intensity) 
  mean_pixel_intensity$cell_id = ids
  mean_pixel_intensity = rev(mean_pixel_intensity)
  
  
  write.table(mean_pixel_intensity, 
              file = paste0("Output/",filename, "_quant.csv"),
              sep=",", row.names = FALSE, quote=FALSE)
  
}



detect_membranes_git1_old <-function(img, channels, offset, sigma, cutoff)
{
  message("########################CELLS########################")
  
  cells_blurred = gblur(img[,,gfp_channel], sigma = sigma)
  ct = thresh(cells_blurred, offset = offset)
  
  
  cm = bwlabel(ct)
  FS = computeFeatures.shape(cm)
  sel <- which(FS[,"s.perimeter"] < cutoff)
  membranes <-rmObjects(cm, sel)
  
  
  res <- remove_edge_membranes(membranes, img, channels)
  
  
  
  
  
  message(paste0("Number of cells detected on first pass: ", length(table(membranes))))
  
  list(removed = res$removed, membranes = res$membranes, FM = res$FM)
  
}

read_in_channels <- function(dic,gfp,cmac)
{

  dic_path=dic
  gfp_path=gfp
  cmac_path=cmac
  
  


  dic = readImage(file.path(dic_path))
  gfp = readImage(file.path(gfp_path))
  cmac = readImage(file.path(cmac_path))
  
  list(cmac = cmac, gfp = gfp, dic = dic)

}
# Function to remove cells that are smaller than 1 standard deviation from the mean size. +-2 standard deviations allows some noise in,
# so far 1 sd is not excluding too many.  Returns cell objects as well as computed feateres.
remove_small_cells <- function(segs, cells)
{
  
  areas <- segs$FC[,"s.area"]
  cutoff <- mean(areas) - sd(areas)
  small <- which(areas < cutoff)
  print(paste0("Cutoff point for cell area (in pixels): ", cutoff))
  print(paste0("Number of cells removed for being below cutoff: ", length(small)))
  segs$segs <- rmObjects(segs$segs, small)
  
  segs$segs <- bwlabel(segs$segs)
  segs$FC = computeFeatures(segs$segs, ref = cells[,,gfp_channel], xname = "seg")
  return(segs)
}


# Wrapper function to find cell objects
find_cell_bodies <- function(cells)
{
  print("################CELLS###############")
  
  segs <- first_pass_new(cells)
  #  segs <- remove_small_cells(segs, cells)
  #segs <- remove_odd_shapes(segs, cells)
  print(paste0("Number of cells after all removals: ", length(table(segs$segs)) ))
  return(segs)
  
}
first_pass <- function(cells)
{
  cb = gblur(cells[,,gfp_channel], sigma = 2)
  ct = thresh(cb, offset = 0.005)
  
  cm = bwlabel(ct)
  FS = computeFeatures.shape(cm)
  sel <- which(FS[,"s.area"] < 500)
  ce = rmObjects(cm, sel)
  
  
  
  cell_segments_full =  fillHull(cell_segments) - cell_segments
  segs <- fillHull(cell_segments_full)
  segs <- bwlabel(segs)
  FC = computeFeatures.shape(segs)
  print(paste0("Number of cells detected on first pass: ", length(table(segs))))
  list(segs = segs, FC = FC)
}
convert_to_grayscale<- function(img)
{
  gfp_gray = channel(img$gfp, 'gray')
  cmac_gray = channel(img$cmac, 'gray')
  
  # RBioformats has .nd2 capability
  img = combine(cmac_gray, gfp_gray, img$dic)
  return(img)
}

scale_intensities <- function(img)
{
  img[,,cmac_channel] = normalize(img[,,cmac_channel])
  img[,,gfp_channel] = normalize(img[,,gfp_channel])
 
  return(img)
}

create_masks <- function(img)
{

  vmask = thresh(img[,,cmac_channel],w=30,h=30,offset = sd(img[,,cmac_channel]))
  vmask = opening(vmask, makeBrush(5,shape="disc"))
  vmask = fillHull(vmask)
  vmask = bwlabel(vmask)
  
  FV = computeFeatures(vmask,img[,,cmac_channel], xname = "vacuoles")
  
  
  ctmask = opening(img[,,gfp_channel] > 0.1, makeBrush(15,shape="disc"))
  cmask = propagate(img[,,gfp_channel], seeds = vmask, mask = ctmask )
  
  list(vacuoles = vmask, cytoplasm = ctmask, cells = cmask, FV = FV)

}


# First pass cell detection.  Objects are detected via blurring, simple threshholding, and noise removal.
# The last 4 lines allow us to return only cells with a complete membrane detectable

detect_cells_old <-function(img)
{
  message("########################CELLS########################")
  
  cells_blurred = gblur(img[,,gfp_channel], sigma = 2)
  ct = thresh(cells_blurred, offset = 0.005)
  
  cm = bwlabel(ct)
  FS = computeFeatures.shape(cm)
  sel <- which(FS[,"s.area"] < 500)
  membranes <-rmObjects(cm, sel)
  
  # membranes <- bwlabel(membranes)
  
  
  cell_bodies <- fillHull(membranes)
  cell_bodies <- bwlabel(cell_bodies)
  
  
  FCCB <- computeFeatures.shape(cell_bodies)
  
  
  
  mask <- list(cell_bodies = cell_bodies, FCCB = FCCB)
  
  mask <- remove_odd_shapes(mask, img)
  
  cell_bodies <- mask$cell_bodies
  FCCB <- mask$FCCB
  
  res <- remove_small_cells(cell_bodies, FCCB, img)
  
  cell_bodies <- res$cell_bodies
  FCCB <- res$FCCB
  
  
  cell_inner <- fillHull(membranes) - membranes
  cell_inner <- bwlabel(cell_inner)
  cell_inner <- fillHull(cell_inner)
  
  cell_inner <- cell_inner*cell_bodies
  cell_inner <- bwlabel(cell_inner)
  
  
  FCI <- computeFeatures.shape(cell_inner)
  sel <- which(FCI[,"s.area"] < 500)
  cell_inner <-rmObjects(cell_inner, sel)
  
  #  cell_inner <- bwlabel(cell_inner)
  
  
  membranes <- cell_bodies - cell_inner
  # membranes <- bwlabel(membranes)
  
  #FCM <- computeFeatures(membranes, ref = img[,,gfp_channel], xname = "membrane")
  
  
  
  FCCI <- computeFeatures.shape(cell_inner)
  message(paste0("Number of cells detected on first pass: ", length(table(cell_bodies))))
  message(paste0("Number of inner detected on first pass: ", length(table(cell_inner))))
  
  list(
    cell_inner = cell_inner,
    removed = mask$removed,
    FCCI = FCCI, 
    membranes = membranes, 
    #FCM = FCM, 
    cell_bodies = cell_bodies, 
    FCCB = FCCB)
  
}


# Function to compute circularity of a cell object
circularity <- function(seg_num, FC) 
{
  4 * pi * ((FC[seg_num,"s.area"]) / (FC[seg_num,"s.perimeter"]^2))
} 

# Removes odd shapes (those with circularity less than 1.)  The circularity cuttoff can be played with,
# 1 generally will be enough to remove budding cells.  
remove_odd_shapes <- function(mask, img)
{
  circ <- lapply(1:length(table(mask$cell_bodies)) -1, circularity, FC = mask$FCCB)
  
  circ<-unlist(circ)
  odd_shapes = which(circ < 1) #budding cells
  message(paste0("Number of cells removed for being oddly shaped: ", length(odd_shapes)))
  
  removed_ind = which(!(seq(1, length(table(mask$cell_bodies))) %in% odd_shapes))
  removed = rmObjects(mask$cell_bodies, removed_ind)
  
  cell_bodies <- rmObjects(mask$cell_bodies, odd_shapes)
  cell_bodies <- bwlabel(cell_bodies)
  FCCB <- computeFeatures.shape(cell_bodies)#(cell_bodies, ref = img[,,gfp_channel], xname = "seg")
  
  list(cell_bodies = cell_bodies, FCCB = FCCB, removed = removed)
  
  
  
  
  
}


remove_small_cells <- function(cell_bodies, FCCB, img)
{
  
  areas <- FCCB[,"s.area"]
  cutoff <- mean(areas) - sd(areas)
  small <- which(areas < cutoff)
  message(paste0("Cutoff point for cell area (in pixels): ", format(cutoff, nsmall = 3)))
  message(paste0("Number of cells removed for being below cutoff: ", length(small)))
  cell_bodies <- rmObjects(cell_bodies, small)
  
  cell_bodies <- bwlabel(cell_bodies)
  FCCB = computeFeatures.shape(cell_bodies)
  list(cell_bodies = cell_bodies, FCCB = FCCB)
  
}



create_cell_segments<- function(cells)
{
  cb = gblur(cells[,,gfp_channel], sigma = 2)
  ct = thresh(cb, offset = 0.005)
  
  cm = bwlabel(ct)
  FS = computeFeatures.shape(cm)
  sel <- which(FS[,"s.area"] < 500)
  ce = rmObjects(cm, sel)
  
  cell_segments <- thresh(ce)
  #cell_segments <- stackObjects(cell_segments,gfp_gray)\
  cell_segments <- bwlabel(cell_segments)
  
  cell_segments_filled <- fillHull(cell_segments)

  return(cell_segments)

}


# If a cell has more than 1 vacuole detected, keeps the biggest one and removes all others.  Also removes cells with
# no detected vacuoles.  After all removals features are recomputed.  This is by far the most computationally expensive function, and can be improved by 
# refining the parameters used to detect the vacuoles.  Still a work in progress.  
remove_multiples <- function(cell_info, vmask, img)
{
  
  l = length(table(cell_info$cell_inner))
  l = l-1
  
  
  maxes <- vector("list", l)
  empty_cells <- vector("list", l)
  
  for(i in seq(1:l))
  {
    
    cseg = cell_info$cell_inner == i
    vcount <- table(cseg*vmask)
    vcount <- vcount[-1]
    #message(vcount)
    if( length(vcount) > 0)
    {
      
      v <- as.numeric(names(as.list(vcount)))
      # message(v)
      vcount <- unname(vcount)
      maxes[i] <- v[which(vcount == max(vcount))]
    }
    else
    {
      # message(paste0(i, " ZERO VCOUNT"))
      empty_cells[i] = i
    }
  }
  
  empty_cells <- unlist(empty_cells)
  
  
  message(paste0("Number of empty cells: ", length(empty_cells)))
  maxes <- unlist(maxes)
  extras <- which(!(seq(1:length(table(vmask)[-1])) %in% maxes))
  
  
  
  vmask <- rmObjects(vmask, extras)
  vmask <- bwlabel(vmask)
  
  
  # Not going to recompute features here for cell bodies/cell inner.  Can add that if needed but as of now it's not necessary and it will slow things down
  # recomputation of membrane features is needed, however
  
  
  
  #remove empty cells
  cell_inner <- rmObjects(cell_info$cell_inner, empty_cells)
  cell_inner <- bwlabel(cell_inner)
  #remove empty cells
  cell_bodies <- rmObjects(cell_info$cell_bodies, empty_cells)
  cell_bodies <- bwlabel(cell_bodies)
  FCCI <- computeFeatures(cell_inner, ref = img[,,gfp_channel])
  FCCB <- computeFeatures(cell_bodies, ref = img[,,gfp_channel])
  
  #remove empty membranes
  #membranes <- rmObjects(cell_info$membranes, empty_cells)
  # membranes <- bwlabel(membranes)
  membranes = cell_bodies - cell_inner
  
  FCM <- computeFeatures(membranes, ref = img[,,gfp_channel], xname = "membrane")
  FV <- computeFeatures(vmask, ref= img[,,gfp_channel], xname = "vac")
  message(paste0("Number of cells after empty cells removed: ", length(table(cell_bodies))))
  
  list(vacuoles = vmask, 
       membranes = membranes,
       FCM = FCM,
       #  FV = FV,
       FCCB = FCCB,
       cell_inner = cell_inner,
       cell_bodies = cell_bodies,
       FCCI = FCCI)
}




img_vac_seg <- function(cmac, cell_segments, masks)
{
  c = paintObjects(cell_segments, tgt = cmac, col = 'brown', closed = TRUE)
  output = paintObjects(masks$vacuoles, c, col='yellow', thick = TRUE)
  return(output)
}



examine_image_segment <- function(segment_number,img,masks)
{

    
  segment = masks$cells == segment_number
  
  gfp_mask = 1- img[,,gfp_channel]
  gfp_mask = gfp_mask*segment
  
  
  vseg = masks$vacuoles == segment_number
  cell_segment = gfp_mask*gfp_gray*20
  
  p <- paintObjects(vseg, tgt = cell_segment, col="red")
  
  return(p)
}

vac_stats <- function(segment_number, FV)
{
  list(majoraxis = FV[segment_number,"vacuoles.0.m.majoraxis"], radiussd = 
  FV[segment_number,"vacuoles.0.s.radius.sd"])
  
}
print_roi_labels <- function(channels,cell_segments, masks)
{
  
  img <- img_vac_seg(channels$cmac,cell_segments, masks$vacuoles)
  labels <- c(1:length(table(masks$vacuoles) - 1))
  center_pts = get_center_points(masks$vacuoles)
  display(img, method = "raster")
  text(center_pts, label = labels, col = "white")
}


get_center_points <- function(vmask)
  {
  vac_centers = computeFeatures.moment(vmask)
  center_pts = vac_centers[,c("m.cx","m.cy")] 
  return(center_pts)
}

exclude_small_rois <- function(FV, masks, cell_segments)
{
  too_small <- which(FV[,"vacuoles.0.s.area"] < 200)
  masks$vacuoles <- rmObjects(masks$vacuoles, too_small)
  masks$cells <- rmObjects(masks$cells, too_small)
  cell_segments = cell_segments*masks$cell
  list(masks = masks, segments = cell_segments)
  
}



roi_labels <- function(vmask)
{
  FV = computeFeatures.moment(vmask)
  vac_centers =FV[,c("m.cx", "m.cy")]
  
  
  
  labels <- c(1:length(table(vmask)) )
  display(vmask, method = "raster")
  text(vac_centers, label = labels, col = "red")
  
  
}


>>>>>>> a787f86eca1ccbda411758d1e2b020fbdf89d150
