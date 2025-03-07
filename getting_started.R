library(openeo)
connect(host = "https://openeo.dataspace.copernicus.eu")

login()

list_collections()
describe_collection("SENTINEL2_L2A")
list_processes()
describe_process("aggregate_temporal")
collection_viewer(x="SENTINEL2_L2A")
process_viewer("load_collection")

get_sample(
  SENTINEL2_L2A,
  replace_aoi = TRUE,
  spatial_extent = NULL,
  execution = "sync",
  immediate = TRUE,
  con = NULL
)

p <- processes()
datacube <- p$load_collection(
  id = "SENTINEL2_L2A",
  spatial_extent=list(west = 5.14, south = 51.17, east = 5.17, north = 51.19),
  temporal_extent=list("2021-02-01", "2021-04-30"),
  bands=list("B08","B04","B02")
)

evi_cube <- p$reduce_dimension(data = datacube, dimension = "bands",
                               reducer = function(data,context) {
  B08 = data[1]
  B04 = data[2]
  B02 = data[3]
  (2.5 * (B08 - B04)) / sum(B08, 6.0 * B04, -7.5 * B02, 1.0)
})

###################################################################
p = processes()

# get the collection list to get easier access to the collection ids, via auto completion
collections = list_collections()

# get the formats
formats = list_file_formats()
# load the initial data collection and limit the amount of data loaded
# note: for the collection id and later the format you can also use the its character value
data = p$load_collection(id = collections$`SENTINEL1_GRD`,
                         spatial_extent = list(west=16.06, 
                                               south=48.06,
                                               east=16.65,
                                               north=48.35),
                         temporal_extent = c("2017-03-01", "2017-06-01"),
                         bands = c("VV"))

# apply the SAR backscatter process to the data
datacube_sar = p$sar_backscatter(data = data, coefficient = "sigma0-ellipsoid")


# create three monthly sub-datasets, which will later be combined back into a single data cube
march = p$filter_temporal(data = datacube_sar,
                          extent = c("2017-03-01", "2017-04-01"))

april = p$filter_temporal(data = datacube_sar,
                          extent = c("2017-04-01", "2017-05-01"))

may = p$filter_temporal(data = datacube_sar,
                        extent = c("2017-05-01", "2017-06-01"))

# The aggregation function for the following temporal reducer
agg_fun_mean = function(data, context) {
  mean(data)
}

march_reduced = p$reduce_dimension(data = march,
                                   reducer = agg_fun_mean,
                                   dimension = "t")

april_reduced = p$reduce_dimension(data = april,
                                   reducer = agg_fun_mean,
                                   dimension = "t")

may_reduced = p$reduce_dimension(data = may,
                                 reducer = agg_fun_mean,
                                 dimension = "t")

# Each band is currently called VV. We need to rename at least the label of one dimension, 
# because otherwise identity of the data cubes is assumed. The bands dimension consists 
# only of one label, so we can rename this to be able to merge those data cubes.
march_renamed = p$rename_labels(data = march_reduced,
                                dimension = "bands",
                                target = c("R"),
                                source = c("VV"))

april_renamed = p$rename_labels(data = april_reduced,
                                dimension = "bands",
                                target = c("G"),
                                source = c("VV"))

may_renamed = p$rename_labels(data = may_reduced,
                              dimension = "bands",
                              target = c("B"),
                              source = c("VV"))

# combine the individual data cubes into one
# this is done one by one, since the dimensionalities have to match between each of the data cubes
merge_1 = p$merge_cubes(cube1 = march_renamed,cube2 = april_renamed)
merge_2 = p$merge_cubes(cube1 = merge_1, cube2 = may_renamed)

# rescale the the back scatter measurements into 8Bit integer to view the results as PNG
rescaled = p$apply(data = merge_2,
                   process = function(data,context) {
                     p$linear_scale_range(x=data, inputMin = -20,inputMax = -5, outputMin = 0, outputMax = 255)
                   })


# store the results using the png format and set the create a job options
result = p$save_result(data = rescaled,format = formats$output$PNG, options = list(red="R",green="G",blue="B"))

# create a job
job = create_job(graph = result, title = "S1 Example R", description = "Getting Started example on openeo.org for R-client")

# then start the processing of the job and turn on logging (messages that are captured on the back-end during the process execution)
start_job(job = job, log = TRUE)

list_results(job = job)


##########################

# If the .glm file is an R object
load("SD.SpeciesDistribution.invasiveAlien.Heracleum-mantegazzianum.glm")

# If it's saved as an RDS file
model <- readRDS("SD.SpeciesDistribution.invasiveAlien.Heracleum-mantegazzianum.glm")
library(multiplex)
data <- read.gml("SD.SpeciesDistribution.invasiveAlien.Heracleum-mantegazzianum.glm",
                 as = c("srt", "array"), directed = TRUE, coords = FALSE)


library(sf)
gml_data <- st_read("Data.glm")
shapefile_data <- st_read("DGURBA_RG_01M_2021.shp")

library(ggplot2)
ggplot(data = shapefile_data) +
  geom_sf() +
  theme_minimal() +
  ggtitle("Map of Shapefile")

