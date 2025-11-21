# README
Daniel Enns
2025-11-21

- [This document…](#this-document)
- [Methods](#methods)
  - [<span class="toc-section-number">0.1</span> Library](#library)
  - [<span class="toc-section-number">0.2</span> Combining Spatial
    Data](#combining-spatial-data)
    - [<span class="toc-section-number">0.2.1</span> Point
      features](#point-features)
    - [<span class="toc-section-number">0.2.2</span> Watershed
      delineation](#watershed-delineation)
    - [<span class="toc-section-number">0.2.3</span> Network
      building](#network-building)
  - [<span class="toc-section-number">0.3</span> Spatial Data
    quantification](#spatial-data-quantification)
  - [<span class="toc-section-number">0.4</span> Data analysis and
    handling](#data-analysis-and-handling)
    - [<span class="toc-section-number">0.4.1</span> EDA](#eda)
    - [<span class="toc-section-number">0.4.2</span> Spatial
      autocorrelation](#spatial-autocorrelation)
    - [<span class="toc-section-number">0.4.3</span> Data cleaning and
      feature engineering](#data-cleaning-and-feature-engineering)
  - [<span class="toc-section-number">0.5</span> Modeling](#modeling)
    - [<span class="toc-section-number">0.5.1</span> Set-up tasks and
      learners](#set-up-tasks-and-learners)
    - [<span class="toc-section-number">0.5.2</span> Tuning and tuning
      evaluation](#tuning-and-tuning-evaluation)
- [Data availability](#data-availability)
- [References](#references)

# This document…

is a detailed explanation of the data and model generating process for
the paper “Urban land use and related infrastructure limit ecological
status of freshwater invertebrates”. It describes the procedures of how
to delineate watersheds, how to build a stream network for routing and
how the variables were quantified. Futher, details about an initial
explorative data analysis, analysis of spatial autocorrelation, as well
as data cleaning and feature engineering are described. Finally, model
construction, which includes task and learner setup and model tuning and
tuning evaluation are explained.

# Methods

## Library

The following packages are needed for the analysis:

``` r
# general
library(readxl)
library(tidyverse)
library(dtplyr)
library(data.table)
library(cowplot)

# Saptial analysis
library(geodata)
library(stars)
library(sf)
library(spdep)
library(nngeo)
library(sfnetworks)
library(tidygraph)
library(whitebox)
library(osmextract)
library(httr2)
library(tmap)
library(units)

# modeling
library(mlr3verse)
library(mlr3spatiotempcv)
library(ranger)
library(iml)
library(pdp)
library(doParallel)
library(parallelMap)
```

## Combining Spatial Data

### Point features

For the spatial analysis all necessary shape and raster files are
projected into EPSG 25832.

The freshwater invertebrate sampling sites, as well as the data on
WWTPs, stormwater overflows (SOs) and barriers can readily be loaded
into R. Sampling sites located on large Rivers (e.g. Rhine, Main, Weser)
were excluded, due to incomplete catchment data.

``` r
# load and filter invertebrate sampling sites
mzb <- read_delim("./MZB.csv", delim = ";") %>% 
  bind_cols(mzb, stream_net[st_nearest_feature(mzb, stream_net),]) %>% 
  select(-geom) %>% 
  filter(!GEWKZ  %in% c("2", "238", "24", "41", "4", "44", "239152", "23932", "2396", "23988", "24779892"))
```

Highway and railroad networks can be read in from OpenStreetMap pbf
files via the `osmextract` package. These can be intersected with the
stream network to create point features. In the osm highway and railroad
networks, each lane / track is represented as an individual line
feature, so that multiple crossings are present in locations with high
lane / track density.

``` r
# Extract highway and railroad networks
transport_net <- oe_read("./hessen-latest.osm.pbf", extra_tags = "railway") %>%
  filter(highway == "motorway" | railway == "rail") %>% 
  st_transform(st_crs(stream_net))

# Intersect with stream network
crossings <- st_intersection(stream_net, transport_net) %>% 
  st_cast("POINT")
```

### Watershed delineation

In order to delineate upstream watersheds for each bio sampling site the
SRTM GL1 30m (OpenTopography 2013) digital elevation model will be used
to extract flow accumulation and pointer grids, from which the
watersheds are calculated. However, the delineation algorithms struggles
in regions with low profiles (namely the Rhine basin), creating an
inaccurate stream raster and watersheds. To improve accuracy, the
existing stream network is burned into the DEM following these steps in
QGIS:

1.  Reproject raster network to stream CRS (EPSG: 25832)

2.  dissolve the stream network

3.  Run the Grass r.carve algorithm with a stream width of 60 m, depth
    of 3 m and with no flat areas in direction of stream flow allowed
    (Note: the algorithm does not work on latitude-longitude CRS, see
    [r.carve GRASS GIS
    manual](https://grass.osgeo.org/grass-stable/manuals/r.carve.html))

The algorithm will take some time to run (approximately 1 hour and 18
minutes with an Intel core i7, 16 GB RAM). Finally, watersheds can be
extracted using the `whitebox` package (Wu and Brown 2025), following
the instructions from [Gannon,
2024](https://vt-hydroinformatics.github.io/Quarto_Book/15-Watershed-Delineation-and-Analysis.html).

``` r
# Breach and fill pits in raster
wbt_breach_depressions_least_cost(
  dem = "./SRTMGL1_30m_2px_burned.tif", 
  output = "./srtm_breach.tif", 
  dist = 7,
  fill = T
)

# Create flow accumulation & pointer grids
wbt_d8_flow_accumulation(
  input = "./srtm_breach.tif",
  output = "./D8FA.tif"
)

wbt_d8_pointer(
  dem = "./srtm_breach.tif",
  output = "./D8pointer.tif"
)

# Extract streams
wbt_extract_streams(
  flow_accum = "./D8FA.tif",
  output = "./stream_raster.tif",
  threshold = 700
)

# Snap points to stream raster
wbt_jenson_snap_pour_points(
  pour_pts = "./MZB.shp",
  streams = "./stream_raster.tif",
  output = "./MZB_snap.shp",
  snap_dist = 300
)

# Delineate watersheds
wbt_watershed(
  d8_pntr = "./D8pointer.tif",
  pour_pts = "./MZB_snap.shp",
  output = "./watersheds.tif"
)
```

Figure 1 demonstrates the structure of the watersheds, which were build
consecutively between multiple pour points.

![Figure 1: Example watershed. Each watershed is build consecutively
between specified
pourpoints.](README_files/figure-commonmark/watershed_example_figure-1.png)

The created watershed raster can be read in and converted to polygon
features, which subsequently be used to intersect the [CORINE Land
Cover](https://land.copernicus.eu/en/products/corine-land-cover)
polygons. The area of each intersection can be calculated using the
`st_area()` function.

``` r
# Read in Watersheds and convert to polygon 
ws_poly <- read_stars("./watersheds.tif") %>% 
  st_as_sf(as_points = F, merge = T)

# Create dataframe with land cover areas
landcover_ws <- st_intersection(landcover, ws_poly) %>% 
  mutate(area = st_area(.)) %>% st_drop_geometry() %>% 
  as.data.frame() %>% filter(!is.na(type))

# Summarize land cover types in each watershed
landcover_sum <- landcover_ws %>% group_by(watersheds.tif, type) %>% 
  summarize(area = sum(area)) %>% 
  pivot_wider(names_from = type, values_from = area)

# Combine watersheds with summarized land cover
watersheds_lc <- left_join(ws_poly, landcover_sum, by = "watersheds.tif")
```

### Network building

The Hessian stream network is a collection of 100 m long line segments,
each containing a stream ID and a segment ID. The stream ID is designed
in such a way, that it consecutively builds upon the ID of a previous
stream segment, adding numbers to it whenever the stream splits in
upstream direction. The segment number increases towards the source of a
stream (Fig. 2). Each 100 m segment already contains information on the
in-stream and surrounding habitat structural quality (1 - pristine to
7 - completly altered).

![Figure 2: Example of stream- and segment ID system. The segment number
increases in the direction of the source. The stream ID of smaller
confluence resembels always that of the bigger stream, in which it flows
into, plus additional
numbers.](README_files/figure-commonmark/ID_example_figure-1.png)

The `sfnetwork` package (van der Meer et al. 2024) provides many useful
functions for network analysis and routing operations, where the main
function `as_sfnetwork()` builds a network by connecting aligning nodes.
An initial build reveals that many stream segments visually seem
connected, but are separated from their nearest node or edge by a few
meters. Thus, to build a useful network for routing the precision of the
stream network `sf` object needs to be lowered (for more information see
Details on `st_set_precision()` and `st_as_binary()` help page). To what
precision to round depends, in general, on the maximum distance of
dangling nodes and if the network can be build by connecting nodes with
edges, or just by connecting nodes. The precision will be set to 70
meters (be aware of the object CRS units!).

In the resulting network, the right streams are connected (Fig. 3) and
it can be used for routing operations. Performance of routing can be
greatly increased by simplifying the network and removing pseudo-nodes
(more info
[here](https://luukvdmeer.github.io/sfnetworks/articles/sfn02_preprocess_clean.html)).

``` r
# Lower precision to 70 meters
str_net_70m <- st_set_precision(stream_net, set_units(70, "m"))

# Create sf_network object, simplify and smooth pseudo-nodes
network <- as_sfnetwork(str_net_70m) %>% convert(to_spatial_simple) %>% 
  convert(to_spatial_smooth)
```

![Figure 3: Example for how precision changes the connectedness of the
network. The numbers indicate from which node a stream
originates.](README_files/figure-commonmark/Precision_example_figure-1.png)

Sampling sites and point stressors can be blend in as nodes into the
created network.

``` r
# Blend in sampling sites, wastewater treatment plants and dams
network_blend <- st_network_blend(network, mzb) %>% st_network_blend(.,wwtp) %>% 
  st_network_blend(.,dams)
```

## Spatial Data quantification

For the quantification of point- and polygon features and their
attributes I wrote two custom functions, `st_shift()` and
`upstream_summarize()`. The first function is required within the latter
and shifts a set of point geometries a proportional distance towards
another specified point. The `upstream_summarize()` function builds a
sub-network from a given point to all its vertices (stream sources) and
extracts the node data, from which it counts the number of specified
nodes and sums up specified attribute values. Further, it can calculate
the minimal, maximal, or average network distance between the start
point and a set of nodes specified by their IDs with the `IDs` argument.
If no such nodes are present in the sub-network, `inf` or `NaN` values
are returned. Additionally the Shreve stream magnitude, defined as the
sum of sources upstream from a specified node, can be calculated. The
function can sum up attribute values of provided polygons the following
way: First, the set of nodes present in the sub-network are shifted
towards their centroid by 0.1% of their length, to avoid including
adjacent polygons. Then, it creates a filtering mask by selecting all
polygons which are touched by the shifted nodes and fills in ‘holes’ in
the set of polygons. By setting a threshold, the mask can be shrunken to
avoid selecting adjacent polygons. This mask is finally used to select
all polygons present in the sub-network, from which their specified
attributes are summarized. The function contains the following set of
arguments:

- `net` : The complete network with blended in points of interest
- `start` : Row number of node from which to rout upstream
- `node_cols` : Names of attribute columns to summarize, present in
  nodes
- `IDs` : ID columns for sets of nodes, which should be unique in each
  set
- `area` : `sf` object with only polygon geometries
- `area_cols` : names of attribute columns to summarize, present in
  `area`
- `dist` : Either `'min'`, `'max'`, `'mean'`, or `'all'`. Calculates
  distances between start and nodes, specified by `IDs`.
- `threshold` : Value (in polygon CRS units) by which the polygon mask
  should be shrunken.
- `Shreve` : logical, should Shreve stream magnitude be calculated?

Before this function can be used, the points from which the function
should rout must be extracted as nodes from the network and their row
name must be saved as a new column.

``` r
# Extract nodes of invertebrate sampling sites and add row names
mzb_nodes <- st_as_sf(network_blend, "nodes") %>% filter(!is.na(ID_SITE)) %>% 
  mutate(ID_NODE = row.names(nodes)[with(nodes, !is.na(ID_SITE))])
```

`upstream_summarize()` can be used in combination with `rowwise()` and
`mutate()` to perform the action over multiple points. Depending on the
number of nodes within the entire network, this task can take quite a
lot of time. Luckily this task can be parallelized to save a lot of
time.

``` r
# Set up cluster
num_cores <- (detectCores()-1)
num_cores %>% makeCluster() %>% registerDoParallel()

# Split data into chunks
data_chunks <- mzb_nodes %>% 
  split(., ceiling(seq_along(row_number(.)) / (length(row_number(.)) %/% num_cores)))

# Apply upstream_summarize parallel and row-wise over all sampling site nodes 
mzb_data_complete <- foreach(chunk = data_chunks, .combine = rbind, .packages = c("dplyr", "dtplyr","sf","sfnetworks","tidygraph","nngeo")) %dopar% {
  chunk %>% rowwise() %>% 
  mutate(upstream_summarize(
    net = network_blend,
    start = ID_NODE,
    node_cols = c("population_equivalents", "dam_discharge"),
    IDs = c("ID_WWTP", "ID_DAMS"),
    area = watersheds_lc,
    area_cols = c("Agriculture", "Urban", "semi-Natural"),
    dist = T,
    threshold = 30)
    )
}
```

## Data analysis and handling

### EDA

Before modeling it is good practice to get a sens of the target and
explanatory variables. Here, multiple models should perform two
different tasks, regression and multi-class classification. The
respective target variables are the multi metric index (MMI) and the
resulting ecological quality class (EQC). The MMI is a summarization of
multiple core metrics, compared to reference conditions (for more
details see (Hering et al. 2006)). It ranges from zero to one,
representing bad to good condition respectively.

![Figure 4: Histogram of the site averaged multi metric
index.](README_files/figure-commonmark/MMI_histogram-1.png)

![Figure 5: Histogram of the site averaged ecological quality
class.](README_files/figure-commonmark/EQC_histogram-1.png)

From the EQC histogram (Fig. 5) a strong class imbalance is observable.
Since the classes are ordered categorical, the minority could be merged
with an adjacent class. Although that would mean a loss of information,
most practitioners are trying to achieve a “good or better” status for
which a combination of the top two classes is legitimate.

Some covariates follow the same method of quantification and are
therefore auto correlated by default. Sites at bigger reaches will have
more point stressors in the upstream watershed compared to sites further
upstream, therefore the different stressor covariates are automatically
correlated. The correlation can be displayed in a correlation plot (Fig.
6).

![Figure 6: Correlation plot of
covariates.](README_files/figure-commonmark/autocorrelation_plot-1.png)

Tree ensemble learning algorithms are robust to autocorrelated
covariates when it comes to overall model performance, but it is good to
keep in mind that these relations between variables exist.

### Spatial autocorrelation

Toblers first law states that *“everything is related to everything
else, but near things are more related than distant things”*. This
degree can be quantified by the Moranes Index, which can be calculated
with distance-based weight matrix for points.

``` r
# K-nearest neighbours 
suppressPackageStartupMessages(require(deldir)) 
data_nb_knn <- knn2nb(knearneigh(mzb_data_complete, k = 1))  

# Distance neighbours 
dsts <- unlist(nbdists(data_nb_knn, mzb_data_complete)) 
summary(dsts) 
max_1nn <- max(dsts)  
data_nb_dst <- dnearneigh(mzb_data_complete, d1 = 0, d2 = 0.75*max_1nn)  
data_lw <- nb2listw(data_nb_dst, style = "W", zero.policy = T)  

# Moran statistics 
moran.test(mzb_data_complete$MMI_mean, data_lw, alternative = "greater") 
moran.plot(mzb_data_complete$MMI_mean, data_lw) 
```

       Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
        1.0   332.7   840.2  1197.1  1663.7  8718.0 


        Moran I test under randomisation

    data:  data_all$MMI_mean  
    weights: data_lw  
    n reduced by no-neighbour observations  

    Moran I statistic standard deviate = 39.571, p-value < 2.2e-16
    alternative hypothesis: greater
    sample estimates:
    Moran I statistic       Expectation          Variance 
         0.4357911968     -0.0005844535      0.0001216123 

The index shows a significant deviation from the expectation, signifying
strong spatial autocorrelation. This can also be visualized with a
Moranes plot (Fig. 7).

![Figure 7: Moranes plot for the site averaged multi metric
index.](README_files/figure-commonmark/moranes_plot-1.png)

Consequently, the spatial heterogenity should be kept in mind for model
tuning and cross validation (Schratz et al. 2019).

### Data cleaning and feature engineering

Tree ensemble learning algorithms can handle covariates of different
types and are robust to auto-correlated covariates. Nevertheless, it is
good practice to clean and modify desired features before modeling. For
each sampling site, land use types in the upstream watershed were
calculated as proportions to the total watershed area.

Distance variables where transformed into relative inverse distances:

$$
dist_{inverse} = 1 - x / (x_{max} + 5 m)
$$

Where $x$ is the distance to the closest point stressor in meters and
$x_{max}$ is the maximal distance to the closest point stressor in
meters.

The Wastewater impact was calculated as the cumulative number of
population equivalents from all WWTPs upstream of a sampling site,
divided by the Shreve stream magnitude of the sampling site.

The set of variables, used in the consecutive section were:

| Variable | Description | Average | Range |
|----|----|----|----|
| count_Dam | count of Dams upstream | 0.48 | 0-11 |
| count_WWTP | count of wastewater treatment plants upstream | 3.1 | 0-40 |
| count_SO | count of stormwater overflows upstream | 30.86 | 0-401 |
| count\_ HWC | count of higway crossings upstream | 3.97 | 0-58 |
| count_RRC | count of railroad crossings upstream | 9.36 | 0-153 |
| dist_barrier | distance to nearest hydrological barrier upstream | 0.91 | 0-1 |
| dist_WWTP | distance to nearest wastewater treatment plant upstream | 0.48 | 0-1 |
| dist_SO | distance to nearest stormwater overflow upstream | 0.84 | 0-1 |
| dist_HWC | distance to nearest highway crossing upstream | 0.29 | 0-1 |
| dist_RRC | distance to nearest railroad crossing upstream | 0.4 | 0-1 |
| prop_agri | proportion of agricultural area in upstream catchment | 0.4 | 0-0.99 |
| prop_urb | proportion of urban area in upstream catchment | 0.1 | 0-0.68 |
| prop_nat | proportion of semi-natural area in upstream catchment | 0.45 | 0-1 |
| WWTP_imp | wastewater impact | 0.02 | 0-1 |
| HQ | habitat quality |  |  |
| ST | stream type |  |  |
| Shreve | Shreve stream magnitude | 32.25 | 1-295 |

Table 1: Variables used for modelling

## Modeling

To understand how the site averaged MMI and the averaged EQC are
influenced by the variables in table 1, the relationships are modelled
using ML algorithms. Since the MMI is a numeric and the EQC is a
categorical variable, the tasks at hand are regression and
classification respectively. The `mlr3` package (Lang, Schratz, and
Becker 2025) offers an extensive framework for model building, tuning,
interpretation, etc.

### Set-up tasks and learners

The regression and classification tasks are set-up by dropping the
respective unnecessary target from the data and specify which features
are coordinates for the spatial resampling technique mentioned later on.
Then, the tasks are partitioned into training and test sets, with a
default ratio of 67% to 33% respectively, using the same splits.

``` r
# seed
set.seed(1234)

# create task
tsk_reg <- data4model_fin %>% select(-EQC_mean) %>% as_task_regr_st(target = "MMI_mean", coordinate_names = c("X", "Y"))
tsk_classif <- data4model_fin %>% select(-MMI_mean) %>% as_task_classif_st(target = "EQC_mean", coordinate_names = c("X", "Y"))

# create partition
splits_reg <- partition(tsk_reg)
tsk_reg_train <- tsk_reg$clone()$filter(splits_reg$train)
tsk_reg_test <- tsk_reg$clone()$filter(splits_reg$test)

tsk_classif_train <- tsk_classif$clone()$filter(splits_reg$train)
tsk_classif_test <- tsk_classif$clone()$filter(splits_reg$test)
```

Data preprocessing, feature selection and modelling will be done using a
graph learner (Fig. 8), a construct where data flows from one operator
to another, potentially branching into different paths. This allows to
examine how preprocessing methods, feature selection and choice of model
algorithm influences the final model performance. The first branch in
the graph is a division between performing a pca on scaled and
standardized data, only scaling the data or performing no operation. The
output is then one-hot encoded, turning categorical features into
numerical ones, so that they can be used in the subsequent models. The
second branch divides the paths into performing a feature selection
based on feature importance, or to perform no operation. The final
branch is the division between the model algorithms used, namely xgboost
and random forest (ranger).

``` r
# regression
# create and plot graph learner
graph <- as_learner(
  po("branch", c("nop", "scale", "pca"), id = "branch_1") %>>% 
    gunion(list(
      po("nop"), 
      po("scale"), 
      po("pca", scale. = T, rank. = to_tune(p_int(1, 10, depends = branch_1.selection == "pca")))
    )) %>>%
    po("unbranch", id = "unbranch_1") %>>%
    po("encode", method = "one-hot") %>>%
    po("branch", c("nop", "filter"), id = "branch_2") %>>%
    gunion(list(
      po("nop", id = "nop_2"),
      po("filter", filter = flt("importance", learner = lrn("regr.xgboost")), filter.nfeat = to_tune(p_int(1, tsk_reg$ncol, depends = branch_2.selection == "filter")))
    )) %>>% 
    po("unbranch", id = "unbranch_2") %>>%
    po("branch", c("xgboost", "ranger"), id = "branch_3") %>>%
    gunion(list(
      lrn("regr.xgboost"),
      lrn("regr.ranger")
    )) %>>%
    po("unbranch", id = "unbranch_3")
)

# set additional params
graph$param_set$set_values(
  branch_1.selection = to_tune(c("nop", "scale", "pca")),
  branch_2.selection = to_tune(c("nop", "filter")),
  branch_3.selection = to_tune(c("xgboost", "ranger")),
  regr.xgboost.eta = to_tune(p_dbl(1e-04, 1, logscale = T, depends = branch_3.selection == "xgboost")),
  regr.xgboost.nrounds = to_tune(p_int(1, 5000, depends = branch_3.selection == "xgboost")),
  regr.xgboost.max_depth = to_tune(p_int(1, 20, depends = branch_3.selection == "xgboost")),
  regr.xgboost.colsample_bytree = to_tune(p_dbl(0.1, 1, depends = branch_3.selection == "xgboost")),
  regr.xgboost.colsample_bylevel = to_tune(p_dbl(0.1, 1, depends = branch_3.selection == "xgboost")),
  regr.xgboost.lambda = to_tune(p_dbl(0.001, 1000, logscale = T, depends = branch_3.selection == "xgboost")),
  regr.xgboost.alpha = to_tune(p_dbl(0.001, 1000, logscale = T, depends = branch_3.selection == "xgboost")),
  regr.xgboost.subsample = to_tune(p_dbl(0.1, 1, depends = branch_3.selection == "xgboost")),
  regr.ranger.mtry.ratio = to_tune(p_dbl(0, 1, depends = branch_3.selection == "ranger")),
  regr.ranger.replace = to_tune(p_lgl(depends = branch_3.selection == "ranger")),
  regr.ranger.sample.fraction = to_tune(p_dbl(0.1, 1, depends = branch_3.selection == "ranger")),
  regr.ranger.num.trees = to_tune(p_int(1, 2000, depends = branch_3.selection == "ranger"))
)

# classification
graph_classif <- as_learner(
  po("branch", c("nop", "scale", "pca"), id = "branch_1") %>>% 
    gunion(list(
      po("nop"), 
      po("scale"), 
      po("pca", scale. = T, rank. = to_tune(p_int(1, 10, depends = branch_1.selection == "pca")))
    )) %>>%
    po("unbranch", id = "unbranch_1") %>>%
    po("encode", method = "one-hot") %>>%
    po("branch", c("nop", "filter"), id = "branch_2") %>>%
    gunion(list(
      po("nop", id = "nop_2"),
      po("filter", filter = flt("importance", learner = lrn("classif.xgboost")), filter.nfeat = to_tune(p_int(1, tsk_classif$ncol, depends = branch_2.selection == "filter")))
    )) %>>% 
    po("unbranch", id = "unbranch_2") %>>%
    po("branch", c("xgboost", "ranger"), id = "branch_3") %>>%
    gunion(list(
      lrn("classif.xgboost"),
      lrn("classif.ranger")
    )) %>>%
    po("unbranch", id = "unbranch_3")
)

graph_classif$param_set$set_values(
  branch_1.selection = to_tune(c("nop", "scale", "pca")),
  branch_2.selection = to_tune(c("nop", "filter")),
  branch_3.selection = to_tune(c("xgboost", "ranger")),
  classif.xgboost.eta = to_tune(p_dbl(1e-04, 1, logscale = T, depends = branch_3.selection == "xgboost")),
  classif.xgboost.nrounds = to_tune(p_int(1, 5000, depends = branch_3.selection == "xgboost")),
  classif.xgboost.max_depth = to_tune(p_int(1, 20, depends = branch_3.selection == "xgboost")),
  classif.xgboost.colsample_bytree = to_tune(p_dbl(0.1, 1, depends = branch_3.selection == "xgboost")),
  classif.xgboost.colsample_bylevel = to_tune(p_dbl(0.1, 1, depends = branch_3.selection == "xgboost")),
  classif.xgboost.lambda = to_tune(p_dbl(0.001, 1000, logscale = T, depends = branch_3.selection == "xgboost")),
  classif.xgboost.alpha = to_tune(p_dbl(0.001, 1000, logscale = T, depends = branch_3.selection == "xgboost")),
  classif.xgboost.subsample = to_tune(p_dbl(0.1, 1, depends = branch_3.selection == "xgboost")),
  classif.ranger.mtry.ratio = to_tune(p_dbl(0, 1, depends = branch_3.selection == "ranger")),
  classif.ranger.replace = to_tune(p_lgl(depends = branch_3.selection == "ranger")),
  classif.ranger.sample.fraction = to_tune(p_dbl(0.1, 1, depends = branch_3.selection == "ranger")),
  classif.ranger.num.trees = to_tune(p_int(1, 2000, depends = branch_3.selection == "ranger"))
)
```

![Figure 8: Graph learner used for
modeling.](README_files/figure-commonmark/graph_learner_plot-1.png)

### Tuning and tuning evaluation

The graph learner can be tuned to find an optimal hyperparameter
constellation. Hyperparameters can be e.g. the number of principal
components or features to keep, internal parameters which control model
building and performance, but also the choice of which path to follow in
the graph. The performance of each model, resulting out of one
hyperparameter constellation, is measured by the root mean squared error
(rmse) and the explained variance for the regression model and by
classification accuracy for the classification model. To get a robust
performance estimation, a 10 fold k-nearest neighbour cross validation
based on coordinates is used. The algorithm used to find the ‘best’
hyperparameter constellation is a random search with 500 evaluation
iterations. It consecutively builds models with random hyperparameter
constellations, evaluates them via the selected resampling method on the
train data and stores the performance and constellation in an archive.

``` r
# setup tuning
# regression
instance <- ti(
  task = tsk_reg_train,
  learner = graph,
  resampling =  rsmp("spcv_coords", folds = 10),
  measure = msrs(c("regr.rsq", "regr.rmse")),
  terminator = trm("evals", n_evals = 500)
)

# classification
instance_classif <- ti(
  task = tsk_classif_train,
  learner = graph_classif,
  resampling =  rsmp("spcv_coords", folds = 10),
  measure = msr("classif.acc"),
  terminator = trm("evals", n_evals = 500)
)

tuner <- tnr("random_search")

# tuning process
# regression
future::plan("multisession", workers = detectCores()-1)  
tuner$optimize(instance)

# classification
future::plan("multisession", workers = detectCores()-1)  
tuner$optimize(instance_classif)
```

The tuning exhaustion can be visualized by plotting the sorted measures
(root mean squared error and classification accuracy) against the
evaluation steps (Fig. 9). Both curves inflection point is situated at
around 200 iterations and is surpassed after 500 evaluation rounds,
where no great improvements in RMSE and Accuracy is observed.

![Figure 9: Tuning exhaustion plot shows the performance measuere
against the number of tuning
iterations](README_files/figure-commonmark/tuning_exhaustion_plot-1.png)

Since the optimization of the regression model uses two performance
measures simultaneously, the relation of both rmse and explained
variance to each other can be visualized via a pareto front plot for
each hyperparameter constellation (fig. 10). The most optimal
constellation maximizes the explained variance and minimizes the rmse.

![Figure 10: Pareto front plot shows the relation of two performance
measures which were tuned
simultaneously](README_files/figure-commonmark/pareto_plot-1.png)

The performance of the path selection for the graph learner is
visualized in figure 11. The lowest rmse results from selecting no data
preprocessing, selecting no filtering as well as selecting a XGBoost
algorithm, whereas the highest classification accuracy results from
selecting feature scaling, selecting no filtering and selecting the
random forest algorithm. The optimal set of hyperparameters is presented
in table 2.

![Figure 11: Path selection plot shows how the selected path influences
the model performance for the top 10 performing tuning
iterations](README_files/figure-commonmark/path_selection_plot-1.png)

| Model          | Hyperparameter            | Value     |
|----------------|---------------------------|-----------|
| Regression     | xgboost.alpha             | -2.990768 |
| Regression     | xgboost.eta               | -3.074318 |
| Regression     | xgboost.lambda            | 4.878687  |
| Regression     | xgboost.max_depth         | 5         |
| Regression     | xgboost.nrounds           | 272       |
| Regression     | xgboost.subsample         | 0.7369354 |
| Regression     | xgboost.colsample_bylevel | 0.3088645 |
| Regression     | xgboost.colsample_bytree  | 0.6463654 |
| Classification | ranger.mtry.ratio         | 0.7950264 |
| Classification | ranger.num.trees          | 279       |
| Classification | ranger.replace            | TRUE      |
| Classification | ranger.sample.fraction    | 0.1824466 |

Table 2: Optimal set of hyper parameters

Finally the model can be trained on the training data, using the optimal
hyperparameter configuration, and then be evaluated on the test data to
get an unbiased performance measure. We can compare the model
performances against a baseline model, which in this case would just be
the average for the regression and the majority class for the
classification model. Model results are easily interpretable by
inspecting feature importance ranking and partial dependence plots (for
more info see
[here](https://mlr3book.mlr-org.com/chapters/chapter12/model_interpretation.html)).

``` r
# evaluate RMSE and explained variance
xgb_tuned <- graph
xgb_tuned$param_set$values = instance$archive$learner_param_vals(uhash = "dd69e4d8-0237-4e3e-8ccd-dbec49be2fa0")
xgb_tuned$train(tsk_reg_train)

prediction <- xgb_tuned$predict(tsk_reg_test)
prediction$score(msrs(c("regr.rsq", "regr.rmse")))

# compare to baseline
reg_baseline <- lrn("regr.featureless")$train(tsk_reg_test)
prediction_baseline <- reg_baseline$predict(tsk_reg_test)
prediction_baseline$score(msrs(c("regr.rsq", "regr.rmse")))

# evaluate Accuracy
ranger_tuned <- graph_classif
ranger_tuned$param_set$values = instance_classif$archive$learner_param_vals(uhash = "0fbd7e51-5420-42e6-9543-f7231006db55")
ranger_tuned$train(tsk_classif_train)

prediction_classif <- ranger_tuned$predict(tsk_classif_test)
prediction_classif$score(msr("classif.acc"))
prediction_classif$confusion

# compare to baseline
classif_baseline <- lrn("classif.featureless")$train(tsk_classif_train)
prediction_cl_baseline <- classif_baseline$predict(tsk_classif_test)
prediction_cl_baseline$score(msr("classif.acc"))
```

The results of model performances, as well as feature importance and
partial dependence are available at the respective paper.

# Data availability

Data on WFD invertebrate sampling and the Hessian stream network were
kindly provided by the Hessian state office for nature, environment and
geology (HLNUG). Data on wastewater treatment plants, stormwater
overflows are available for download from the
[WRRL-viewer](https://wrrl.hessen.de/mapapps/resources/apps/wrrl/index.html?lang=de).
Further, data on hydrological barriers are available from the [Amber
Barrier Atlas](https://amber.international/european-barrier-atlas/).
OpenStreetMap data can be downloaded as pbf files at
[geofabrik](https://download.geofabrik.de/). CORINE land cover
shapefiles and rasters are available at [Copernicus Land Monitoring
Service](https://land.copernicus.eu/en/products/corine-land-cover). The
SRTM GL1 30m digital elevation model can be downloaded from
[OpenTopography](https://portal.opentopography.org/raster?opentopoID=OTSRTM.082015.4326.1).

# References

<div id="refs" class="references csl-bib-body hanging-indent"
entry-spacing="0">

<div id="ref-hering2006" class="csl-entry">

Hering, Daniel, Christian K. Feld, Otto Moog, and Thomas Ofenböck. 2006.
“Cook Book for the Development of a Multimetric Index for Biological
Condition of Aquatic Ecosystems: Experiences from the European AQEM and
STAR Projects and Related Initiatives.” In, 311–24. Springer
Netherlands. <https://doi.org/10.1007/978-1-4020-5493-8_22>.

</div>

<div id="ref-Lang2025" class="csl-entry">

Lang, Michel, Patrick Schratz, and Marc Becker. 2025. *Mlr3verse: Easily
Install and Load the ’Mlr3’ Package Family*.
<https://mlr3verse.mlr-org.com>.

</div>

<div id="ref-opentopography2013" class="csl-entry">

OpenTopography. 2013. “Shuttle Radar Topography Mission (SRTM) Global.”
<https://doi.org/10.5069/G9445JDF>.

</div>

<div id="ref-schratz2019" class="csl-entry">

Schratz, Patrick, Jannes Muenchow, Eugenia Iturritxa, Jakob Richter, and
Alexander Brenning. 2019. “Hyperparameter Tuning and Performance
Assessment of Statistical and Machine-Learning Algorithms Using Spatial
Data.” *Ecological Modelling* 406 (August): 109–20.
<https://doi.org/10.1016/j.ecolmodel.2019.06.002>.

</div>

<div id="ref-vanderMeer2024" class="csl-entry">

van der Meer, Lucas, Lorena Abad, Andrea Gilardi, and Robin Lovelace.
2024. *Sfnetworks: Tidy Geospatial Networks*.
<https://luukvdmeer.github.io/sfnetworks/>.

</div>

<div id="ref-Wu&Brown2025" class="csl-entry">

Wu, Qiusheng, and Andrew Brown. 2025. *’Whitebox’: ’WhiteboxTools’ r
Frontend*. <https://CRAN.R-project.org/package=whitebox>.

</div>

</div>
