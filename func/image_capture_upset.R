# image_capture_upset <- function(image_name){
#   #' Upset plot formating for saving as tiff
#   library(svglite)
#   svglite(
#     image_name,
#     # width = 1000,
#     # height = 500#,
#     # res = 200,
#     # compression = "none"
#   )
# }

image_capture_upset <- function(image_name){
  #' Upset plot formating for saving as tiff
  pdf(
    image_name,
    width = 2000,
    height = 200,
    bg = 'white',
    paper = 'a4r'
    # res = 200,
    # compression = "none"
  )
}
