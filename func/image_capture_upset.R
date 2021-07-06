image_capture_upset <- function(image_name){
  tiff(
    image_name,
    width = 1000,
    height = 500,
    res = 200,
    compression = "none"
  )
}
