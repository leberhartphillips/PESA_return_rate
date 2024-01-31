#### set up R enviroment ----
# a vector of all the packages needed in the project
packages_required_in_project <- c("sf", "raster", "sp", "rgdal", "tidyverse", "spatstat", "geosphere", "adehabitatLT", "patchwork", "viridis")

# of the required packages, check if some need to be installed
new.packages <- 
  packages_required_in_project[!(packages_required_in_project %in% 
                                   installed.packages()[,"Package"])]

# install all packages that are not locally available
if(length(new.packages)) install.packages(new.packages)

# load all the packages into the current R session
lapply(packages_required_in_project, require, character.only = TRUE)

# set the home directory to where the project is locally based (i.e., to find 
# the relevant datasets to import, etc.
here::set_here()


# Find fonts from computer that you want. Use regular expressions to do this
# For example, load all fonts that are 'verdana' or 'Verdana'
extrafont::font_import(pattern = "[V/v]erdana", prompt = FALSE) 

# check which fonts were loaded
extrafont::fonts()
extrafont::fonttable()
extrafont::loadfonts() # load these into R

# define the plotting theme to be used in subsequent ggplots
luke_theme <- 
  theme_bw() +
  theme(
    text = element_text(family = "Verdana"),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    axis.title.x = element_text(size = 10),
    axis.text.x  = element_text(size = 8), 
    axis.title.y = element_text(size = 10),
    axis.text.y = element_text(size = 8),
    strip.text = element_text(size = 10),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(size = 0.5, colour = "grey40"),
    axis.ticks.length = unit(0.2, "cm"),
    panel.border = element_rect(linetype = "solid", colour = "grey"),
    legend.position = c(0.1, 0.9)
  )

#### data import ----
# import data and convert to a spatial feature with WGS84 coordinate system, 
# then transform to the UTM zone for Barrow (4N) for a planar coordinate system 
# to use spatstat tools
df_sf_utm <- 
  read.csv("data/data_pecs.csv", colClasses = "character") %>% 
  mutate(longit = as.numeric(longit), latit = as.numeric(latit)) %>% 
  st_as_sf(., coords = c("longit", "latit"), 
           crs = 4326) %>% 
  st_transform(., crs = 32604)

#### generate year-specific kernel density surfaces ----
# create an empty list to store the results
kd_list <- list()

# a for loop which cycles through the years of the dataset and generates
# an occupancy probability surface across the study site from which to draw the
# random locations that males could return to.
for(i in 1:length(df_sf_utm %>% pull(year_) %>% unique())){
  year_i <- unique(df_sf_utm$year_)[i]
  
  # extract all data for a given year
  pts_ <- 
    df_sf_utm %>% 
    dplyr::filter(year_ == year_i)
  
  # determine the minimum convex polygon of the data
  mcp_ <- st_convex_hull(st_union(pts_))
  
  # convert mcp and point to a spatstat planar point pattern
  mcp_w <- as.owin(mcp_)
  mcp_w_ppp <- as.ppp(st_coordinates(pts_), mcp_w)
  
  # create a kernal density pixel image
  mcp_kd <- density.ppp(mcp_w_ppp)
  
  # convert the kernal density estimate to a probability
  mcp_kd_df <- 
    as.data.frame(mcp_kd) %>% 
    mutate(probability = value/sum(value),
           year_ = year_i)
  
  # store the results
  kd_list[[i]] <- mcp_kd_df
  
}

# bind all elements of the output list
kd_bind_rows <- do.call(bind_rows, kd_list)

# plot the year-specific kernel densities of the male locations
plot_kd <-
  ggplot(kd_bind_rows, aes(x = x, y = y, color = probability)) +
  geom_point() +
  scale_color_viridis() + 
  facet_wrap(~ year_) +
  labs(x = "easting",
       y = "northing",
       color = "occupancy\nprobability") +
  luke_theme +
  theme(legend.position = c(.8, .2),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

plot_kd

#### return location simulation ----
# specify the number of random points to generate
nrand_pts = 1000

# create an emply list to store results
result_list <- list()

pesa_df <- 
  read.csv("data/data_pecs.csv", colClasses = "character") %>% 
  mutate(date = as.POSIXct(paste(year_, "01", "01", sep = "-")),
         longit = as.numeric(longit),
         latit = as.numeric(latit))

returning_birds <- 
  read.csv("data/data_pecs.csv", colClasses = "character") %>% 
  group_by(ID) %>% 
  summarise(n_obs = n()) %>% 
  filter(n_obs > 1) %>% 
  pull(ID) %>% 
  unique()

# specify coordinate columns
coordinates(pesa_df) <- c("longit","latit")

# define spatial projection as WGS84
proj4string(pesa_df) <- CRS("+init=epsg:4326")

# make a spatial trajectory object of each animal's encounter history
tr_pesa_wp <- 
  as.ltraj(xy = sp::coordinates(pesa_df),
           date = pesa_df$date,
           id = pesa_df$ID)

# convert the trajectory object back to a dataframe, rename columns, and
# calculate the cumulative distance traveled over the observation
dist_y1_y2 <-
  ld(tr_pesa_wp) %>% 
  dplyr::rename(distance = dist,
                ID = id) %>% 
  filter(!is.na(distance)) %>% 
  group_by(ID) %>%
  arrange(date) %>%
  slice(1) %>% 
  mutate(distance = distance * 100000,
         year_ = year(date),
         years_ago = dt/60/60/24/365) %>% 
  dplyr::select(ID, distance, year_, years_ago)

# a for loop that cycles through each male that returned and 1) samples a number 
# random points (specified by the "nrand_pts" object above) according to the 
# year-specific kernel density surfaces shown above, 2) calculates the distance
# between each of these random points and the location of the previous year
for(i in 1:length(unique(returning_birds))){
  
  focal_year <- 
    pesa_df %>% 
    as.data.frame() %>% 
    filter(ID == returning_birds[i]) %>% 
    arrange(date) %>% 
    slice(2) %>% 
    mutate(year_ = year(date)) %>% 
    pull(year_)
  
  focal_bird <- 
    pesa_df %>%
    as.data.frame() %>% 
    filter(ID == returning_birds[i]) %>% 
    arrange(date) %>% 
    slice(2) %>% 
    mutate(year_ = year(date)) %>% 
    dplyr::select(-date) %>% 
    st_as_sf(., coords = c("longit", "latit"), 
             crs = 4326) %>% 
    st_transform(., crs = 32604)
  
  pesa_df_ <- 
    pesa_df %>% 
    as.data.frame() %>% 
    filter(year(date) == focal_bird$year_) %>% 
    st_as_sf(., coords = c("longit", "latit"), 
             crs = 4326) %>% 
    st_transform(., crs = 32604)
  
  # determine the minimum convex polygon of the data
  mcp_ <- st_convex_hull(st_union(pesa_df_))
  
  # convert mcp and point to a spatstat planar point pattern
  mcp_w <- as.owin(mcp_)
  mcp_w_ppp <- as.ppp(st_coordinates(pesa_df_), mcp_w)
  
  # create a kernal density pixel image
  mcp_kd <- density.ppp(mcp_w_ppp)

  # convert the kernal density estimate to a probability
  mcp_kd_df <- 
    as.data.frame(mcp_kd) %>% 
    mutate(probability = value/sum(value))
  
  random_row_index <- sample(seq_len(nrow(mcp_kd_df)), 
                             size = nrand_pts, 
                             prob = mcp_kd_df$probability)
  
  rand_points <- 
    mcp_kd_df %>% 
    slice(random_row_index) %>% 
    st_as_sf(coords = c("x", "y"), crs = 32604)
  
  distances <- 
    data.frame(rand_dist = as.numeric(st_distance(rand_points, focal_bird,  by_element = FALSE)), 
               true_dist = filter(dist_y1_y2, ID == returning_birds[i]) %>% pull(distance), 
               ID = returning_birds[i]) %>% 
    mutate(dist_diff = rand_dist - true_dist)
  
  # Append the result to the list
  result_list[[i]] <- distances
  
}

result_bind_rows <- do.call(bind_rows, result_list)

median(result_bind_rows$rand_dist)
quantile(result_bind_rows$rand_dist, c(0.025, 0.975))

# plot the year-specific kernel densities of the male locations
plot_distance_density <-
  ggplot() +
  geom_density(data = result_bind_rows, aes(rand_dist), fill = "grey70", color = "white") +
  geom_vline(xintercept = median(dist_y1_y2$distance), color = "white") +
  xlab("dispersal distance between consecutive years") +
  annotate(geom = "text", 
           x = median(dist_y1_y2$distance) + 40, 
           y = 0.0002, label = "median distance of males that did return", 
           angle = 90, color = "white", hjust = 0) +
  annotate(geom = "text",
           x = 1000, y = 0.00005, 
           label = "simulated dispersal distances", color = "grey40") +
  luke_theme +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

plot_distance_density
