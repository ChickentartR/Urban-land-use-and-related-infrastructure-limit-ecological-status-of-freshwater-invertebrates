################
### st_shift ###
################
# Description #### 
# This function shifts single or multiple points x towards another point y 
# by a specified factor of their distance.

# Code ####
st_shift <- function(x, y, factor){

  # disable traceback
  options(error = NULL)
  
  # check for valid inputs
  if(!any(class(x) %in% c("sf", "sfg", "sfc"))) {
    stop("x must be of class sf, sfc or sfg!")
  }
  
  if(!any(class(st_geometry(y)) %in% "sfc_POINT")) {
    stop("y must be of class sf, sfc or sfg, with point geometry!")
  }
  
  if(nrow(st_coordinates(y)) != 1) {
    stop("y must be a single point!")
  }
  
  if(class(factor) != "numeric" | length(factor) != 1) {
    stop("factor must be a single numeric value!")
  }
  
  # extract coordinates and shift
  sh_pts <- st_coordinates(x) + factor * (matrix(rep(st_coordinates(y), nrow(x)), ncol = 2, byrow = TRUE)-st_coordinates(x)) %>% 
  as_tibble()

# convert to sf object
shifted <- st_as_sf(sh_pts, coords = c("X", "Y"), crs = st_crs(x))

return(shifted)
}