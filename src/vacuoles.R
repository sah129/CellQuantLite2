
# First pass vacuole detection.  Creates a moving window that passes over CMAC
# channel.  Returns a list of vacuole objects and associated computed features.
find_vacuoles <- function(cell_info, img, channels)
{
  message("######################VACUOLES#######################")
  
  # Create the moving adaptive window by using the maximum min radius of all 
  # detected cell objects.
  b <- max(cell_info$FM[,"membrane.0.s.radius.min"]) 
  
  # Threshold with an adaptive window and an offset 2 std. dev. from the norm.
  # This is sufficient to pick up only the brightest spots while excluding 
  # haze. 
  vmask = thresh(img[,,cmac_channel],w=b,h=b,offset = 2*sd(img[,,cmac_channel]))
  
  # Removes minor protrosions that are often a result of overflow PM 
  # fluorescence.
  vmask = opening(vmask, makeBrush(5,shape="disc"))
  vmask = bwlabel(vmask)
  
  message(paste0("Number of vacuoles detected on first pass: ", format(length(table(vmask)), nsmall = 4)))
  FV <- computeFeatures(vmask, ref= channels$ref_gfp, xname = "vac")
  FV<-FV[,c("vac.a.b.mean", "vac.0.m.cx", "vac.0.m.cy", "vac.0.s.area")]
  
  res<-list(vacuoles = vmask, FV = FV)
  return(res)
}

# Build the resultant data frame while examining and excluding candidate 
# membranes and their associated vacuoles. Returns data frame, list of 
# fragment PM objects, list of PMs where no vacuoles are found.
exclude_and_bind <- function(mems, vacs)
{

  df <- data.frame(matrix(NA, nrow = length(table(mems$membranes)), ncol = 9))
  names(df) <- c('CellID', 'vacuoles', 'cell_area', 'vac_area', 'PM_vac_ratio', 'cell_mpi', 'vac_mpi', 'pm_center_x', 'pm_center_y')
  
  
  l = length(table(mems$membranes))-1
  empty_cells <- vector("list", l) # cells containing no detected vacuoles
  fragments <- vector("list", l) # PMs that not fully formed
  
  # Loop through list of PMs
  for(i in seq(1:l))
  {
    pm_seg = mems$membranes == i  # isolate PM
    filled_seg = fillHull(pm_seg) # flood fill
    comp <- pm_seg == filled_seg
    
    # If # of pixels in PM is greater than intersection complement, add to 
    # fragment list.  This technique allows whole PMS with some 
    # minor gaps in fluorescence to make the cut. 
    if(length(which(pm_seg == 1)) > length(which(comp == FALSE)))
      fragments[i] = i
    else  
    {
      # If not a fragment, search detected vacuoles within the cell. 
      v_seg <-vacs$vacuoles*(filled_seg-pm_seg)
      vcount <- table(v_seg)
      vcount <- vcount[-1] 
      # If no vacuoles detected, discard
      if( length(vcount) == 0 )
        empty_cells[i] = i
      else
      {
        c_area <-  mems$FM[i, 'membrane.0.s.area']
        v_area <- calc_vac_areas(as.numeric(names(vcount)), vacs$FV)
        
        # Exclude cells containing very small combined vacuole area.
        # The area of the filled membrane would be a much better metric, 
        # however in the interest of space/memory we will use precomputed 
        # areas in lieu of another fillHull operation.
        if(v_area/c_area < .25)
          empty_cells[i] = i
        else  # populate dataframe
        {
          cmpi <- mems$FM[i, 'membrane.a.b.mean']
          vmpi <- calc_vac_mpi(as.numeric(names(vcount)), vacs$FV)
          df[i,] <- list(CellID = i,
                         vacuoles = toString(names(vcount)),
                         cell_area =  c_area,
                         vac_area = v_area,
                         PM_vac_ratio = cmpi/vmpi,
                         cell_mpi = cmpi,
                         vac_mpi = vmpi,
                         pm_center_x = mems$FM[i,'membrane.0.m.cx'],
                         pm_center_y = mems$FM[i,'membrane.0.m.cy'])
        }
      }
    }
  }
  list(df=df, fragments = unlist(fragments), empty_cells = unlist(empty_cells))
}

# Calculate combined areas of all vacuoles in cell
calc_vac_areas <- function(areas, FV)
{
  area = 0
  for(i in areas)
    area = area + FV[i, 'vac.0.s.area']
  return(area)
}

# Calculate mean pixel intensity for all vacuoles in cell.
calc_vac_mpi <- function(mpis, FV)
{
  mpi = 0
  for(i in mpis)
    mpi = mpi + FV[i, 'vac.a.b.mean']
  mpi = mpi/length(mpis)
  return(mpi)
}
