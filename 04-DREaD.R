---
title: "Speciation mode 2025.Isopogon example"
output: html_document
---

```{r setup, include=FALSE, warnings=F}
knitr::opts_chunk$set(echo = TRUE, eval=T, warning = FALSE, message = FALSE)
```

## Data and libraries
==
Load packages.
```{r}
# Libraries
library(terra)# do everything with terra and sf instead of raster/sp/rgeos/rgdal
library(sf)
library(ape)
library(geiger)
library(phytools)
library(moments)
library(apTreeshape)
library(ggplot2)
library(abc)
library(ENMTools)
library(stringr)
```

Load data.
```{r}
# load phy
#phy <- read.tree("input/clade_data/2025.Isopogon/2025.Isopogon.28.nu.65t.2022.tree")
phy <- read.tree("Isopogon.24t.Secapr.tree.tree")

# load occurrence data
#spp <- read.csv("input/clade_data/2025.Isopogon/2025.Isopogon-All-Oct-clean.COMBINE.69t.csv")
spp <- read.csv("01-Isopogon.cleaned.2023.04.combine.Dread.csv")

# load in the country outline
oz_map <- st_read("input/polygons/Australia.shp")

# load the ecoregions to define sea and swa
ibra7 <- st_read("input/polygons/ibra7_regions.shp")

# Following Cook and Crisp paper
sea <- ibra7[ibra7$REG_CODE_7 %in% c("TSR", "TWE", "TSE", "TNM", "TCH", "BEL", "TNS", "KIN", "FUR",
                                     "SCP", "SVP", "VIM", "NCP", "AUA", "SEC", "AUA", "MDD", "KAN", 
                                     "EYB", "FLB", "RIV", "NSS", "SEH", "SYB", "COP", "DRP", "NAN", 
                                     "NET", "NNC", "SEQ", "BBS"), ]
swa <- ibra7[ibra7$REG_CODE_7 %in% c("MAL", "ESP", "WAR", "JAF", "SWA", "GES", "AVW"), ]

# combine into single polygons
swa <- st_union(st_make_valid(swa))
sea <- st_union(st_make_valid(sea))

swa <- st_transform(swa, "+proj=longlat +datum=WGS84")
sea <- st_transform(sea, "+proj=longlat +datum=WGS84")

plot(st_geometry(oz_map))
plot(st_geometry(sea), add=T, col="dark green")
plot(st_geometry(swa), add=T, col="gold")

```

Clean and match the spatial and phylo data. Specifically the names.

```{r}
# which are mismatched?
phy_not_spp <- phy$tip.label[which(!phy$tip.label %in% spp$Scientific.Name)]
spp_not_phy <- spp$Scientific.Name[which(!spp$Scientific.Name %in% phy$tip.label)]

# drop mismatched
spp <- spp[which(!spp$Scientific.Name %in% spp_not_phy ),]
phy <- drop.tip(phy, phy_not_spp)

# double check they match
all(spp$Scientific.Name %in% phy$tip.label)
all(phy$tip.label %in% spp$Scientific.Name)
```

Get species names.
```{r}
species_names <- unique(spp$Scientific.Name)
head(species_names)
```

## Rasterizing range maps

In order to calculate summary statistics to use in the model selection, we need to convert point occurrence data into raster format. Ideally we want to approximate the species geographic range as closely as possible. There are different methods to do this and they each have pros and cons. We will also mask out the ocean cells using the outline of the Australian coastline.

* Direct raster - here we count only grid cells with an occurrence present. (likely underestimates the total range)

* Buffer raster - here we draw a 0.5 degree buffer around each occurrence record and count all grid cells with their centroid overlapping a buffed region.

* Minimum convex polygon - here we draw a convex hull around X% the occurrence points. (likely overestimates the total range)

```{r}
# create a template raster which will set the resolution and extent of the rasters
raster_template <- rast(resolution=0.1, xmin=100, xmax=190, ymin=-50, ymax= 0)
values(raster_template) <- 0
crs(raster_template) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"# create some empty lists to fill up with rasters or spatial points
spp_stack <- list() # for spatial points

for(species_i in 1:length(species_names)){
  
  spp_i <- spp[which(spp$Scientific.Name==species_names[species_i]),]
  spp_i <-  st_set_crs(st_as_sf(spp_i, coords=c("Longitude", "Latitude")),4326)
  spp_stack[[species_i]] <- spp_i
  if(species_i == 1){
  buffer_stack_01 <- rasterizeWin(x=data.frame(st_coordinates(spp_i), 1), y=raster_template, pars=0.1, fun=any)
  buffer_stack_03 <- rasterizeWin(x=data.frame(st_coordinates(spp_i), 1), y=raster_template, pars=0.3, fun=any)
  buffer_stack_05 <- rasterizeWin(x=data.frame(st_coordinates(spp_i), 1), y=raster_template, pars=0.5, fun=any)
  } else {
  buffer_stack_01 <- c(buffer_stack_01,rasterizeWin(x=data.frame(st_coordinates(spp_i), 1), y=raster_template, pars=0.1, fun=any))
  buffer_stack_03 <- c(buffer_stack_03,rasterizeWin(x=data.frame(st_coordinates(spp_i), 1), y=raster_template, pars=0.3, fun=any))
  buffer_stack_05 <- c(buffer_stack_05,rasterizeWin(x=data.frame(st_coordinates(spp_i), 1), y=raster_template, pars=0.5, fun=any))
  }
}


names(buffer_stack_01) <- species_names
names(buffer_stack_03) <- species_names
names(buffer_stack_05) <- species_names
names(spp_stack)       <- species_names


writeRaster(buffer_stack_01, file="output/2025.Isopogon/2025.Isopogon_buffered_01.tif") 
writeRaster(buffer_stack_03, file="output/2025.Isopogon/2025.Isopogon_buffered_03.tif")
writeRaster(buffer_stack_05, file="output/2025.Isopogon/2025.Isopogon_buffered_04.tif") 
saveRDS(spp_stack, file="output/2025.Isopogon/2025.Isopogon_points.rds")      
```

## Plotting rasters

Plot out the spatial records for a single species (2025.Isopogon achaeta)

```{r}
par(mfrow=c(1,1))
plot(st_geometry(oz_map), main="spatial points")
plot(st_geometry(spp_stack[[1]]),add=T, pch=16, cex=0.2)
```

Plot out the raster's of 2025.Isopogon achaeta to get an impression of how the different methods perform.

```{r, fig.width = 4, fig.height = 10}
# for fun plot the richness of the clade
par(mfrow=c(3,1))
plot(buffer_stack_01[[1]], main="buffer 01")
plot(st_geometry(oz_map), add=T)
plot(buffer_stack_03[[1]], main="buffer 03")
plot(st_geometry(oz_map), add=T)
plot(buffer_stack_05[[1]], main="buffer 05")
plot(st_geometry(oz_map), add=T)
```

We can see they different methods give greatly different estimates of geographic ranges. The buffer is probably the most conservative option going forward. We can also stack the rasters to look at the species richness. Can look at the overall maps of species richness to see an impression of where many species are co-occurring. Can see a very clear hotspot for 2025.Isopogon in SWA.

```{r, fig.width = 6, fig.height = 10}
# for fun plot the richness of the clade
buffer_01_richness <- sum(buffer_stack_01,na.rm=T)
buffer_03_richness <- sum(buffer_stack_03,na.rm=T)
buffer_05_richness <- sum(buffer_stack_05,na.rm=T)

par(mfrow=c(3,1))
plot(buffer_01_richness,  main="buffer 01")
plot(st_geometry(oz_map), add=T)
plot(buffer_03_richness, main="buffer 03")
plot(st_geometry(oz_map), add=T)
plot(buffer_05_richness,  main="buffer 05", col = rev(terrain.colors(100)))
plot(st_geometry(oz_map), add=T)
```

## Identify South West species

The next step once we have the species distributions mapped out, is to identify which species fall in SWA vs SEA. We can use the definitions of Cook and CRisp (2013) to define these regions based on the IBRA ecoregion classifications. Here we ask whether each raster has presence values which overlap the SEA or SWA polygons. Then we create a data frame to store the information.

```{r, echo = F}
# do the species occurence records fall within the SEA or SWA?
swa_presence <- lapply(spp_stack, st_intersection, swa)
sea_presence <- lapply(spp_stack, st_intersection, sea)

# if any occurence falls over return a matching character string
swa_presence_l <- lapply(swa_presence , FUN=function(x){if(any(!is.na(x$Scientific.Name))){return("SWA")}else{return("NSWA")}})
sea_presence_l <- lapply(sea_presence , FUN=function(x){if(any(!is.na(x$Scientific.Name))){return("SEA")}else{return("NSEA")}})

# create a data frame summarising this info
region_df <- data.frame(species =names(swa_presence_l), swa=unlist(swa_presence_l), sea=unlist(sea_presence_l))
region_df$region <- NA

# some If statements for each species - is present in SEA, SWA, both, neither?
for(i in 1:nrow(region_df)){
  if(region_df[i, "swa"] == "NSWA" & region_df[i, "sea"] == "NSEA"){region_df[i, "region"] <- "neither" }
  if(region_df[i, "swa"] == "SWA" & region_df[i, "sea"] == "NSEA"){region_df[i, "region"] <- "SWA" }
  if(region_df[i, "swa"] == "SWA" & region_df[i, "sea"] == "SEA"){region_df[i, "region"] <- "both" }
  if(region_df[i, "swa"] == "NSWA" & region_df[i, "sea"] == "SEA"){region_df[i, "region"] <- "SEA" }
}
```

Plot it out on the phylogeny!
Alysium
```{r, fig.width = 6, fig.height = 10}
# plot the phy
pdf(file = "Phylo.geo.2025.Isopogon.pdf", width = 20, height = 80) 
par(mfrow=c(1,1))
plot(phy, adj=0, label.offset=2, lwd=2, cex=0.6)

# make a vector of colours that matches the order of the tips
colours <- region_df$region[match(phy$tip.label, region_df$species)]
colours <- str_replace_all(colours, "SWA", "gold")
colours <- str_replace_all(colours, "SEA", "dark green")
colours <- str_replace_all(colours, "neither", "grey")
colours <- str_replace_all(colours, "both", "red")

# add time scale
axisPhylo()

# add the tiplabels
tiplabels(pch=21, col="black", adj=1, bg=colours, cex=1)

# add the node labels
nodelabels(cex=0.5)
dev.off()
```


We now know where each species lives, we can split the group into a handful of subclades to look at. For example, node 65 would give us a small clade from neither SWA or SEA (2025.Isopogon from northern Australia), node 70 will give us a primarily SEA clade. node 64 would give us the clade with predominantly non-SWA species. Node 77 will give us the large SWA clade. 

```{r}
# nodes to choose (can edit this for a different subset of clades)
nodes <- c(1025,859,945,999,783,644,605)

# create an empty list
clade_trees <- vector("list",length(nodes))

# fill with phylos
for(i in 1:length(nodes)){
  sp <- tips(phy, nodes[i])
  clade_trees[[i]] <- drop.tip(phy, phy$tip.label[which(!phy$tip.label %in% sp)])
}

names(clade_trees) <- paste0("node", nodes)
```

## Empirical summary statistics for ABC

Now we loop over each subclade and calculate the metrics we need to estimate speciation mode from the simulations using model selection.

```{r, message = FALSE}
# load in required functions
scripts <-paste0("DREaD/", c("getSummaryStats_terrasf_fix.R", "helperFunctions_terrasf.R","findSisters.R","summary_statsitics_functions_terrasf.R", "aoc.R")) 
sapply(scripts, source, echo=FALSE, verbose=FALSE)

#Get summary statistics of each clade for DREaD
emp_df <- data.frame(clade         = names(clade_trees),
                 ntips         = NA,
                 ROmean        = NA,
                 ROsd          = NA,
                 ROslope       = NA,
                 ROintercept   = NA,
                 ROskew        = NA,
                 ROkurtosis    = NA,
                 RO0           = NA,
                 RO50          = NA,
                 RO75          = NA,
                 RO90          = NA,
                 RO100         = NA,
                 TOmean        = NA,
                 TOsd          = NA,
                 asymmean      = NA,
                 asymsd        = NA,
                 asymslope     = NA,
                 asymintercept = NA,
                 RDmean        = NA,
                 RDsd          = NA,
                 RDintercept   = NA,
                 RDslope       = NA,
                 BIMOD50       = NA,
                 BIMODE75      = NA,
                 BIMOD90       = NA,
                 BIMOD100      = NA,
                 RSskew        = NA,
                 RSmean        = NA,
                 RSsd          = NA,
                 CI            = NA,
                 Beta          = NA,
                 Gamma         = NA,
                 SI            = NA,
                 ARCslope      = NA,
                 ARCint        = NA)

sp_rasters <- buffer_stack_01
#alysium

# loop over subclades to get summary stats for model selection
for(sub_clade_i in 1:length(clade_trees)){
  
  print(sub_clade_i)
  
  phy_tmp = clade_trees[[sub_clade_i]]
  
  clade_tmp = names(clade_trees)[sub_clade_i]
  
  sp_rasters_tmp <- sp_rasters[[which(names(sp_rasters) %in% phy_tmp$tip.label)]]
  emp_df[sub_clade_i,] <- getSummaryStats(phy_tmp, clade_tmp, sp_rasters_tmp)
  
}

write.csv(emp_df,file=paste("output/2025.Isopogon/2025.Isopogon_empirical_data_v2.csv",sep=""), row.names=F)

```

Now lets load in the simulated data and get the generating mode of speciation for each run.

## Simulated summary statistics for ABC

```{r}
#Load and prepare DREaD simulation results
sim_df <- read.csv("DREaD/simulation_results_subset.csv",stringsAsFactors = F)

# get mode from which sims were generated 
speciation_modes <- as.character(sim_df$speciation_mode)
sim_df <- sim_df[,which(!colnames(sim_df) %in% "speciation_mode")]

```

## ABC

Now lets perform approximate Bayesian computation to predict the most likely mode of speciation based on the similarity between empirical and simulated summary statistics.

```{r}
#Model selection to infer the predominant geographic mode of speciation: ABC
emp_df_s <- emp_df[, which(colnames(emp_df) %in% colnames(sim_df))]
emp_df_s <- emp_df_s[,match(colnames(sim_df), colnames(emp_df_s))]

# get clade names
clades <- emp_df$clade

# create an empty list to put model results
abc_list <- list()

for(sub_clade_i in 1:length(clades)){
  
  # run ABC with logistic regression correction
  abc_log <- postpr(emp_df_s[sub_clade_i,], speciation_modes, sim_df, tol=.05, method="mnlogistic")
  
    # run ABC with lneural network correction
  abc_neu <- postpr(emp_df_s[sub_clade_i,], speciation_modes, sim_df, tol=.05, method="neuralnet", trace=F)
  
  # combine results into a list and store
  res <- cbind(abc_log$pred,abc_neu$pred)
  colnames(res) <- c("ABC_log", "ABC_nn")
  abc_list[[sub_clade_i]] <- round(res, 3)
  
}

# name list
names(abc_list) <- clades

# print results -> higher posterior probability means more likely predominant mode of speciation accoring to that algorithm (neural network or logistic regression)
print(abc_list)


saveRDS(abc_list, file="output/2025.Isopogon/2025.Isopogon_abc_v2.rds")
```

We can see higher support for a symptric model on the predominantly SWA clade, and greater support for an allopatric model on the predominantly non-SWA clade!
