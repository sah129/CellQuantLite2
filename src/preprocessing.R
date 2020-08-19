# Converts all channels to grayscale
convert_to_grayscale<- function(img)
{
  colorMode(img$gfp) <- 'Grayscale'
  colorMode(img$cmac) <- 'Grayscale'
  
  gfp_gray <- getFrame(img$gfp, 1, type = 'render')
  cmac_gray <- getFrame(img$cmac, 1, type = 'render')
  dic_gray = NULL
  
  if(!is.null(img$dic))
  {
      colorMode(img$dic) <- 'Grayscale'
      dic_gray <- getFrame(img$dic,1,type='render')  
  }

  imgn = EBImage::combine(cmac_gray, gfp_gray, dic_gray)
  return(imgn)
}