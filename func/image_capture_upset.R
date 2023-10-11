image_capture_upset <- function(image_name){
  #' Upset plot formating for saving as tiff
  pdf(
    image_name,
    width = 2000,
    height = 200,
    bg = 'white',
    paper = 'a4r'
  )
}
