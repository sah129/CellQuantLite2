
# Function to sort batch of imagesets into an N x 4 array where N = # of files.
# Array has structure [ [filename] [cmac image] [gfp image] [dic image] ]

read_in_imageset_dir <- function(datapath)
{
  dirs <- list.dirs(path = datapath)
  dirs <- dirs[-1]
  
  imagesets <- matrix(NA, length(dirs), 4)
  colnames(imagesets) <- c("file", "cmac", "gfp", "dic")
  
  for(dir in seq_along(dirs))
  {
    imagesets[dir,"file"] = str_remove(dirs[dir], paste0(datapath, "/"))
    
    imagesets[dir, "cmac"] <- list.files(
      path = file.path(dirs[dir]),
      pattern='cmac.*?\\.tif',
      ignore.case = TRUE)
    
    imagesets[dir, "gfp"] <- list.files(
      path = file.path(dirs[dir]),
      pattern='gfp.*?\\.tif',
      ignore.case = TRUE)
    
    imagesets[dir, "dic"] <- list.files(
      path = file.path(dirs[dir]),
      pattern='dic.*?\\.tif',
      ignore.case = TRUE)
    
  }
  return(imagesets)
}

 
# Function to read in tiff files from a single directory.  Tifs must be 
# formatted as CMAC-Frame1, GFP-Frame2.  This is the default output when
# batch converting .nd2s to tiffs in FIJI.
# Returns N x 2 array of structure [[filename] [filepath]] where N=number of 
# files.
read_in_imageset_files <- function(datapath)
{
  
  files <- list.files(path = datapath)
  imagesets <- matrix(NA, length(files), 2)
  colnames(imagesets) <- c("filepath", "filename")
  
  for(file in seq_along(files))
  {
    imagesets[file,"filename"] = str_remove(files[file], ".tif")
    imagesets[file,"filepath"]= files[file]
  }
  return(imagesets)
}


read_in_imageset_interactive <- function(image_files)
{
  
  imagesets <- matrix(NA, nrow(image_files), 2)
  colnames(imagesets) <- c("filepath", "filename")
  
  for(file in seq_along(image_files$name))
  {
    imagesets[file,"filename"] = str_remove(image_files[file,"name"], ".tif")
    imagesets[file,"filepath"]= image_files[file, "datapath"]
  }
  return(imagesets)
}

# Function to read and store channels from image.   Also stores unaltered
# pixel intensity matrices to reference when computing features later.
read_in_channels_old <- function(imageset, datasetpath)
{
  message("#####################################################")
  message(paste0("Examining image: ", imageset["filename"]))
 

  gfp = readImage(file.path(paste0(datasetpath, "/", imageset["filepath"])))[,,2]
  cmac = readImage(file.path(paste0(datasetpath, "/", imageset["filepath"])))[,,1]
  
  ref_gfp = readImage(file.path(paste0(datasetpath, "/", imageset["filepath"])),  as.is = TRUE)[,,2]
  ref_cmac = readImage(file.path(paste0(datasetpath, "/", imageset["filepath"])),  as.is = TRUE)[,,1]
  
  list(cmac = cmac, gfp = gfp, ref_cmac = ref_cmac, ref_gfp = ref_gfp)
}

read_in_channels <- function(imageset, datasetpath)
{
  message("#####################################################")
  message(paste0("Examining image: ", imageset["filename"]))
  
  img <- readImage(file.path(paste0(datasetpath, "/", imageset["filepath"])))
  
  gfp <- img[,,2]
  cmac <- img[,,1]
  dic <- NULL
  
  ref_img <- readImage(file.path(paste0(datasetpath, "/", imageset["filepath"])),  as.is = TRUE)
  
  ref_gfp <- ref_img[,,2]
  ref_cmac <- ref_img[,,1]
  ref_dic <- NULL
  
  if(numberOfFrames(img) == 3)
  {
      dic <- img[,,3]
      ref_dic <- ref_img[,,3]
  }
  
  list(cmac = cmac, 
       gfp = gfp, 
       dic = dic,
       ref_cmac = ref_cmac, 
       ref_gfp = ref_gfp,
       ref_dic = ref_dic)
}

read_in_channels_interactive <- function(imageset)
{
  message("#####################################################")
  message(paste0("Examining image: ", imageset["filename"]))
  
  img <- readImage(file.path(imageset["filepath"]))
  
  gfp <- img[,,2]
  cmac <- img[,,1]
  dic <- NULL
  
  ref_img <- readImage(file.path(imageset["filepath"]),  as.is = TRUE)
  
  ref_gfp <- ref_img[,,2]
  ref_cmac <- ref_img[,,1]
  ref_dic <- NULL
  
  if(numberOfFrames(img) == 3)
  {
    dic <- img[,,3]
    ref_dic <- ref_img[,,3]
  }
  
  list(cmac = cmac, 
       gfp = gfp, 
       dic = dic,
       ref_cmac = ref_cmac, 
       ref_gfp = ref_gfp,
       ref_dic = ref_dic)
}


read_in_channels_single <- function(imgpath)
{
  message("#####################################################")

  
  img <- readImage(file.path(imgpath$datapath))
  
  gfp <- img[,,2]
  cmac <- img[,,1]
  dic <- NULL
  
  ref_img <- readImage(file.path(imgpath$datapath),  as.is = TRUE)
  
  ref_gfp <- ref_img[,,2]
  ref_cmac <- ref_img[,,1]
  ref_dic <- NULL
  
  if(numberOfFrames(img) == 3)
  {
    dic <- img[,,3]
    ref_dic <- ref_img[,,3]
  }
  
  list(cmac = cmac, 
       gfp = gfp, 
       dic = dic,
       ref_cmac = ref_cmac, 
       ref_gfp = ref_gfp,
       ref_dic = ref_dic)
}












