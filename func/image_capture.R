image_capture <- function(image_name){
  tiff(
    image_name,
    width = 2400,
    height = 1200,
    res = 200,
    compression = "none"
  )
}
