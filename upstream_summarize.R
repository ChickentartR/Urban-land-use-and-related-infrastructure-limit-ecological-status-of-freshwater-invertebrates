###################################
### upstream_summarize function ###
###################################
# Description #### 
# This function filters a sfnetwork from a given node to all of its 
# upstream reaches and summarizes information of other specified nodes present 
# in these sections.

# Dependent packages:
# The function requires the dplyr and sfnetworks packages

# Arguments:
# net - the network from which the function calculates upstream reaches. Must be sfnetwork object.

# start - row number of node, from which upstream reache should be calculated.

# IDs - IDs of nodes, which should be counted. Not the row number of nodes!

# node_cols - columns from which the node data should be summarized.

# dist - calculates min, max, or average distance to points specified by IDs.

# area - sf object with polygon geometry for which upstream reach coverage should be calculated. 

# area_cols - columns from which the area data should be summarized.

# threshold - shrinkage value for the area polygon to not intersect with adjacent polygons
#             not covered by upstream reach. Pay attention to CRS units!

# Shreve - should stream magnitude be calculated? Number of sources upstream from start.

# Code ####
upstream_summarize <- function(net, start, IDs, node_cols = NULL, dist = NULL, area = NULL, area_cols = NULL, threshold = NULL, Shreve = NULL) {
  
  # disable traceback
  options(error = NULL)
  
  # Check for valid class of net input
  if (!is.sfnetwork(net)) {
    stop("net must be an sfnetwork object!", call. = F)
  }
  
  # Check for presence of start node in net
  if (!(start %in% rownames(st_as_sf(net, "nodes")))) {
    stop("start node not present in net!", call. = F)
  }
  
  # Check presence of ID-columns in data
  if (any(!(IDs %in% colnames(st_as_sf(net,"nodes"))))) {
    stop("ID-column not found in net!")
  }
  
  # Check presence of columns in data
  if(!is.null(node_cols)){
    if (any(!(node_cols %in% colnames(st_as_sf(net,"nodes"))))) {
      stop("Column name not found in net!")
    }
  }
  
  # Check for valid dist input
  if(!is.null(dist)) {
    if(!(dist %in% c("min", "max", "mean", "all"))) {
      stop("dist must be 'min', 'max', 'mean', or 'all'!", call. = F)
    }
  }
  
  # Check for valid area input
  if (!is.null(area)) {
    if (st_crs(net) != st_crs(area)) {
      stop("net and area are not in the same CRS!", call. = F)}# Check for same net and area CRS
    if (class(area)[1] != "sf" | any(st_geometry_type(area) != "POLYGON")) {
      stop("area must be an sf object which contains only polygon geometries!", call. = F)}  # check for valid network input
  }
  
  if (!is.null(area) & is.null(area_cols)) {
    stop("must provide area_cols if area is given!")
  }
  
  if (!is.null(area_cols)) {
    if (any(!(area_cols %in% colnames(area)))) {
      stop("Column name not found in area!")}
  }
  
  # Check for valid threshold input
  if(!is.null(threshold)) {
    if(!is.numeric(threshold)) {stop("threshold must be a numeric value!", call. = F)}
  }
  
  # Check for valid Shreve input
  if(!is.null(Shreve)) {
    if(!is.logical(Shreve)) {stop("Shreve must be logical!", call. = F)}
  }
  
  # Perform the analysis and filtering steps
  nodes_us <- suppressWarnings(igraph::shortest_paths(net, from = start, to = igraph::V(net), mode = "in")) %>%
    unlist(.$vpath) %>%
    unique()
    
  tab_nodes <- st_as_sf(net, "nodes") %>% mutate(igraph_ID = seq.int(nrow(.)))
  sub_nodes <- tab_nodes[nodes_us,]

  # Count IDs
  tab_sum <- sub_nodes %>% lazy_dt() %>% 
    summarise(across(.cols = all_of(IDs), ~sum(!is.na(.x)), .names = "num_{.col}")) %>% 
    as_tibble()
  
  # summarize node_cols
  if(!is.null(node_cols)){
    col_sum <- sub_nodes %>% lazy_dt() %>% 
      summarise(across(.cols = all_of(node_cols), ~sum(.x, na.rm = T))) %>% 
      as_tibble()
    tab_sum <- bind_cols(tab_sum, col_sum)
  }
  
  # calculate shreve
  if(isTRUE(Shreve)){
    tab_sum <- tab_sum %>% mutate(
      Shreve = st_as_sf(net, "edges") %>% filter(from %in% sub_nodes$igraph_ID & !(from %in% to)) %>% nrow()
    )
  }
  
  # Get OD cost matrix
  if(!is.null(dist)) {
    suppressWarnings(
      tab_dist <- sub_nodes %>% as_tibble() %>% 
        mutate(across(.cols = all_of(IDs), ~ st_network_cost(
        net,
        from = start,
        to = filter(sub_nodes, !is.na(.x)),
        direction = "in"),
        .names = "dist_{.col}"
      )) %>% select(starts_with("dist_")) %>% distinct()
    ) 
    
    # Calculate distances
    if (dist == "min") {
      suppressWarnings(
        tab_dist <- tab_dist %>%  
          summarize(across(.cols = everything(), list(min = ~ min(.x[!is.infinite(.x)]))))
      )
    } else if (dist == "max") {
      suppressWarnings(
      tab_dist <- tab_dist %>%  
        summarize(across(.cols = everything(), list(max = ~ max(.x[!is.infinite(.x)]))))
      )
    } else if (dist == "mean") {
      suppressWarnings(
      tab_dist <- tab_dist %>%  
        summarize(across(.cols = everything(), list(mean = ~ mean(.x[!is.infinite(.x)]))))
      )
    } else if (dist == "all") {
      suppressWarnings(
      tab_dist <- tab_dist %>%  
        summarize(across(.cols = everything(), 
                         list(min = ~ min(.x[!is.infinite(.x)]), 
                              max = ~ max(.x[!is.infinite(.x)]), 
                              mean = ~ mean(.x[!is.infinite(.x)]))
      )) 
      )
    }
    tab_sum <- bind_cols(tab_sum, as_tibble(tab_dist)) 
  }
  
  # select watersheds
  if (!is.null(area)) {
    area_poly <- st_filter(area, st_shift(sub_nodes, st_centroid(st_union(sub_nodes)), 0.001), .predicate = st_intersects) 
    
    if(any(!st_is(area_poly, "POLYGON") | length(st_is(area_poly, "POLYGON")) == 0)) {stop("there were non-polygon geometries created!")}
    else {area_poly <- area_poly %>% st_union() %>% st_remove_holes() %>% st_buffer(dist = -threshold)}
    
    # sample area by filled area polygon and summarize
    area_sum <- st_filter(area, area_poly, .predicate = st_intersects) %>% 
      lazy_dt() %>% summarise(across(.cols = all_of(area_cols), ~sum(.x, na.rm = T))) %>% 
      as_tibble()
    
    tab_sum <- bind_cols(tab_sum, area_sum)
  }
  
  # Return the summarized table with sums
  return(tab_sum)
}