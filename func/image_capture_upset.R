image_capture_upset <- function(image_name){
  #' Upset plot formating for saving as tiff
  tiff(
    image_name,
    width = 1000,
    height = 500,
    res = 200,
    compression = "none"
  )
}
