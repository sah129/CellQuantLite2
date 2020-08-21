# First pass PM detection.  Takes in grayscale image and the reference
# GFP image.  Removes detected membranes that intersect with image edges. 
# Returns a list of removed membrane objects, detected membrane objects,
# and computed features.

detect_membranes <-function(img, channels)
{
  message("########################CELLS########################")
  ct = thresh(img[,,gfp_channel])
  cm = bwlabel(ct)
  fm <- computeFeatures.shape(cm)
  noise <- which(fm[,"s.area"]< 100) # noise removal
  
  membranes <- rmObjects(cm, noise)
  res <- remove_edge_membranes(membranes, img, channels)
  
  message(paste0("Number of cells detected on first pass: ", length(table(membranes))))
  list(removed = res$removed, membranes = res$membranes, FM = res$FM)
}

detect_membranes_new <-function(img, channels, factor, chan, cutoff, cnum)
{
  message("########################CELLS########################")
  
 # chan <- normalize(chan)

  

  g <- gblur(chan*factor, sigma = 2)
  

  
  ct = thresh(g)
  cm = bwlabel(ct)
  fm <- computeFeatures.shape(cm)
  noise <- which(fm[,"s.area"]< cutoff) # noise removal
  
  message(paste0("Number of cells detected on first pass: ", length(table(cm))))
  
  membranes <- rmObjects(cm, noise)

  
 
  res <- remove_edge_membranes(membranes, img, channels, cnum)
  

  message(paste0("Number of cells after noise removal: ", length(table(membranes))))
  list(removed = res$removed, membranes = res$membranes, FM = res$FM)
}
# Removes detected membranes that intersect with image edges.  Saves removed
# membranes and computes complete features on all membranes in reference to
# unaltered GFP channel. 
remove_edge_membranes <-function(membranes,img, channels, cnum)
{
  contours <- ocontour(membranes)
  bound <- list(l = 3, # 3 pixel buffer to account for haze
                r = dim(img[,,cnum$gfp_channel])[1]-3,
                t = 3,
                b = dim(img[,,cnum$gfp_channel])[2]-3)
  
  left <- lapply(contours, function(x){min(x[,1])})
  right <- lapply(contours, function(x){max(x[,1])})
  top <- lapply(contours, function(x){min(x[,2])})
  bottom <- lapply(contours, function(x){max(x[,2])})
  edge_cells <- c(which(left < bound$l),which(right > bound$r),  which(top < bound$t),which(bottom > bound$b))
  edge_cells <- as.numeric(names(edge_cells))
  edge_cells <- unique(edge_cells)
  
  removed_ind = which(!(seq(1, length(table(membranes))) %in% edge_cells))
  removed = rmObjects(membranes, removed_ind)
  
  membranes <- rmObjects(membranes, edge_cells)
  membranes <- bwlabel(membranes)
  
  FM <- computeFeatures(membranes, ref = channels$ref_gfp, xname = "membrane")
  FM<- FM[,c("membrane.a.b.mean", "membrane.0.m.cx", "membrane.0.m.cy", "membrane.0.s.area","membrane.0.s.radius.min")]
  list(removed = removed, membranes = membranes, FM = FM)
}



