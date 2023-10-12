image_capture_ht <- function(image_name){
  cairo_pdf(
    image_name,
    width = 13,
    height = 7,
    bg = 'white',
  )
}
