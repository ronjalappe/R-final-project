##########################################################################################################################################################################################################################
### Burn severity mapping with Sentinel-2
### Author: Ronja Lappe
### Submission: April 2019
##########################################################################################################################################################################################################################

##########################################################################################################################################################################################################################
### Background information 

### Study area:
### - Peatland fire in Meppen, Germany 
### - loacted in the nature reserve "Tinner Dose- Sprakeler Heide" 

### Datasets: 
### - prefire image: S2A 08/06/2018 10:40:21 UTC
### - postfire image: S2B 10/10/2018  10:40:19 UTC
### - preprocessed to 2LA processing level using ESA software "Sen2Cor"
### - AOI shapefile of the nature protected area "Tinner Dose- Sprakeler Heide" (license: https://www.govdata.de/dl-de/by-2-0)

### Script content:
### 1) Data import and conversion 
### 2) Cropping of rasterstack to the extent of AOI
### 3) Calculation of burn indices (NBR, dNBR, MIRBI, dMIRBI)
### 4) Calcuatlion of burned area using predifend thresholds (comparing dNBR and dMIRBI)
### 5) Calculation burn severity levels (according to the USGS burn severity levels) http://un-spider.org/advisory-support/recommended-practices/recommended-practice-burn-severity/in-detail/normalized-burn-ratio
### 6) Quantification of totally burned area and estimation of carbon release 
### 7) Visualization

### Created Output: 
### a) False color composites of pre- and postfire images (B11, B12, B8A)
### b) Burn indices plots for AOI 
### c) Burn area masks
### d) Burn severity rasters
### e) Dataframe including the totally burned area in square meters and the estimated carbon release in tones 
##########################################################################################################################################################################################################################

# load and install packages 
loadandinstall <- function(mypkg) {if (!is.element(mypkg, installed.packages()[,1])){install.packages(mypkg)}; library(mypkg, character.only=TRUE)  }
loadandinstall("sp")
loadandinstall("raster")
loadandinstall("RStoolbox")
loadandinstall("rgdal")
loadandinstall("gdalUtils")
loadandinstall("ggplot2")
loadandinstall("lattice")

##################################### 1) Data import and conversion ##################################### 


# set working directory to the device, where the datasets are saved
# output will be saved here as well 
setwd("/Volumes/FRITZILEIN/R-project")

# relative file-paths for the pre- and post-fire scenes and AOI shapefile 
postfire_img_path <- "./Data/S2B_MSIL2A_20181010T104019_N0206_R008_T31UGU_20181010T161145.SAFE/GRANULE/L2A_T31UGU_A008328_20181010T104018/IMG_DATA/R20m"
prefire_img_path <- "./Data/S2A_MSIL2A_20180806T104021_N0206_R008_T31UGU_20180806T130931.SAFE/GRANULE/L2A_T31UGU_A016307_20180806T104340/IMG_DATA/R20m"
AOI_path <- paste0("./Data/AOI/Tinner_heide.shp")
postfire_tif1 <- "./Data/S2B_MSIL2A_20181010T104019_N0206_R008_T31UGU_20181010T161145.SAFE/GRANULE/L2A_T31UGU_A008328_20181010T104018/IMG_DATA/R20m/T31UGU_20181010T104019_B02_20m.tif"
prefire_tif1 <- "./Data/S2A_MSIL2A_20180806T104021_N0206_R008_T31UGU_20180806T130931.SAFE/GRANULE/L2A_T31UGU_A016307_20180806T104340/IMG_DATA/R20m/T31UGU_20180806T104021_B02_20m.tif"

# convert jp2 to tif files
  # for the postfire image
if(file.exists(postfire_tif1)){                                                   # check if files have already been converted 
  postfire_tif <- list.files(postfire_img_path, pattern = ".tif", full.names = T) # create character vector of files
  postfire_list <- lapply(postfire_tif,function(x) raster(x))                     # list tif-files
  postfire_stack <- stack(postfire_list)                                          # create raster stack from list
} else {
  postfire_jp2 <- list.files(postfire_img_path, pattern = "B.*\\.jp2$", full.names = T) # otherwise convert jp2 to tif files first
  for (file in postfire_jp2) {                                                          
    out_file <- extension(file, 'tif') 
    gdal_translate(src_dataset = file, dst_dataset = out_file, of = "GTiff")           
    postfire_tif <- list.files(postfire_img_path, pattern = ".tif", full.names = T)
    postfire_list <- lapply(postfire_tif,function(x) raster(x))
    postfire_stack <- stack(postfire_list)
  } 
}
  # for the prefire image 
if(file.exists(prefire_tif1)){
  prefire_tif <- list.files(prefire_img_path, pattern = ".tif", full.names = T)
  prefire_list <- lapply(prefire_tif,function(x) raster(x))
  prefire_stack <- stack(prefire_list)
} else {
  prefire_jp2 <- list.files(prefire_img_path, pattern = "B.*\\.jp2$", full.names = T)
  for (file in prefire_jp2) { 
    out_file <- extension(file, 'tif') 
    gdal_translate(src_dataset = file, dst_dataset = out_file, of = "GTiff")
    prefire_tif <- list.files(prefire_img_path, pattern = ".tif", full.names = T)
    prefire_list <- lapply(prefire_tif,function(x) raster(x))
    prefire_stack <- stack(prefire_list)
  } 
}
# rename spectral bands 
names(postfire_stack) <- c("B02","B03","B04","B05","B06","B07","B11","B12","B8A")
names(prefire_stack) <- c("B02","B03","B04","B05","B06","B07","B11","B12","B8A")



##################################### 2) Crop RasterStack to AOI ##################################### 


# load shapefile with AOI 
AOI <- readOGR(AOI_path)  

# transfer projection of RasterStack to AOI 
AOI_proj <- spTransform(AOI,crs(postfire_stack))

# crop and mask RasterStacks to the shape of the AOI polygon
postfire_cropmask <- mask(crop(postfire_stack,AOI_proj),AOI_proj)
prefire_cropmask <- mask(crop(prefire_stack,AOI_proj),AOI_proj)

# rescale pixel values (by factor 10.000 (according to Wegmann et al. 2016)) 
postfire_crop_resc <- postfire_cropmask/10000
prefire_crop_resc <- prefire_cropmask/10000

# plot false color images (B12,B11,B8A)
ggRGB(postfire_crop_resc, r="B12",g="B11",b="B8A", stretch = "lin")+
  ggtitle("FCC post-fire")+
  labs(x="", y="")+
  theme(axis.text.y = element_text(angle = 90))
ggsave("01_FCC_postfire.png")

ggRGB(prefire_crop_resc, r="B12",g="B11",b="B8A", stretch = "lin")+
  ggtitle("FCC pre-fire")+
  labs(x="", y="")+
  theme(axis.text.y = element_text(angle = 90))
ggsave("02_FCC_prefire.png")


##################################### 3) Calculate burn indices ##################################### 

# Function to calculate the Normalized Burn Ratio (NBR)
nbr <- function(img,k,i){
  bk <- img[[k]]
  bi <-img[[i]]
  nbr <- (bk-bi)/(bk+bi)
  return(nbr)
}

# calculate NBR 
NBR_postfire <- nbr(postfire_crop_resc,"B8A", "B12")
NBR_prefire <- nbr(prefire_crop_resc,"B8A", "B12")
dNBR <- NBR_prefire - NBR_postfire

# plot NBR for postfire, prefire and dNBR 
NBR1 <- ggR(NBR_postfire, geom_raster = T)+
  scale_fill_gradient2(low="#67000d", mid="#ef3b2c",high = "#fff5f0", name="NBR")+
  ggtitle("NBR post-fire")+
  labs(x="", y="")+
  theme(axis.text.x = element_text(angle = 45, size = 8))+
  theme(axis.text.y = element_text(angle = 90, size = 8))
ggsave("03_NBR_postfire.png",NBR1)

NBR2 <- ggR(NBR_prefire, geom_raster = T)+
  scale_fill_gradient2(low="#67000d", mid="#ef3b2c",high = "#fff5f0", name="NBR")+
  ggtitle("NBR pre-fire")+
  labs(x="", y="")+
  theme(axis.text.x = element_text(angle = 45, size = 8))+
  theme(axis.text.y = element_text(angle = 90, size = 8))
ggsave("04_NBR_prefire.png",NBR2)

NBR3 <- ggR(dNBR, geom_raster = T)+
  scale_fill_gradient2(low="#fff5f0", mid="#ef3b2c",high ="#67000d", name="NBR")+
  ggtitle("dNBR")+
  labs(x="", y="")+
  theme(axis.text.x = element_text(angle = 45, size = 8))+
  theme(axis.text.y = element_text(angle = 90, size = 8))
ggsave("05_dNBR.png",NBR3)

# Function to calculate the Mid-Infrared Burn Index
mirbi <- function(img,k,i){
  bk <- img[[k]]
  bi <- img[[i]]
  mirbi <- (10*bi-9.8*bk+2)
  return(mirbi)
}

# calculate MIRBI
MIRBI_postfire <- mirbi(postfire_crop_resc,"B11", "B12")
MIRBI_prefire <- mirbi(prefire_crop_resc,"B11", "B12")
dMIRBI <- MIRBI_postfire - MIRBI_prefire

# plot MIRBI postfire, prefire and dMIRBI
MIRBI1 <- ggR(MIRBI_postfire, geom_raster = T)+
  scale_fill_gradient2(low="#fff5f0", mid="#ef3b2c",high = "#67000d", midpoint=1.5,name="MIRBI")+
  ggtitle("MIRBI post-fire")+
  labs(x="", y="")+
  theme(axis.text.x = element_text(angle = 45, size = 8))+
  theme(axis.text.y = element_text(angle = 90, size = 8))
MIRBI1
ggsave("06_MIRBI_postfire.png",MIRBI1)

MIRBI2 <- ggR(MIRBI_prefire, geom_raster = T)+
  scale_fill_gradient2(low="#fff5f0", mid="#ef3b2c",high = "#67000d", midpoint=1.5,name="MIRBI")+  ggtitle("MIRBI pre-fire")+
  labs(x="", y="")+
  theme(axis.text.x = element_text(angle = 45, size = 8))+
  theme(axis.text.y = element_text(angle = 90, size = 8))
MIRBI2
ggsave("07_MIRBI_prefire.png",MIRBI2)

MIRBI3 <- ggR(dMIRBI, geom_raster = T)+
  scale_fill_gradient2(low="#fff5f0", mid="#ef3b2c",high ="#67000d", midpoint=0.5,name="MIRBI")+
  ggtitle("dMIRBI")+
  labs(x="", y="")+
  theme(axis.text.x = element_text(angle = 45, size = 8))+
  theme(axis.text.y = element_text(angle = 90, size = 8))
MIRBI3
ggsave("08_dMIRBI.png",MIRBI3)


################################# 4) Create binary burned area mask ################################# 


# from dNBR, threshold: 0.44 (moderate-high burn severity (according to USGS))
burned_dNBR <- reclassify(dNBR,cbind(-Inf,0.44,NA))
plot(burned_dNBR,col= "red",main="Burned area (dNBR)")
writeRaster(burned_dNBR, filename="09_Burned area mask (dNBR)",format = "GTiff", overwrite = T)

# from dMIRBI, threshold: 1 (according to Lu et al. 2016)
burned_dMIRBI <- reclassify(dMIRBI,cbind(-Inf,1,NA))
burned_dMIRBI_resc <- burned_dMIRBI - 0.8 # to make both scales comparable
plot(burned_dMIRBI_resc,col= "red",main="Burned area (dMIRBI)")
writeRaster(burned_dMIRBI_resc, filename="10_Burned area mask (dMIRBI)",format = "GTiff", overwrite = T)

################################# 5) Create burn severity classes ################################# 


# burn severity classes (based on dNBR)
burned_dNBR_cl <- reclassify(dNBR, c(-Inf,0.27,NA, 0.27,0.44,1, 0.44,0.66,2, 0.66,Inf,3))
writeRaster(burned_dNBR_cl, filename = "11_Burn Severity (dNBR)", format = "GTiff", overwrite = T)

# burn severity classes (based on dMIRBI)
burned_dMIRBI_cl <- reclassify(dMIRBI, c(-Inf,1,NA, 1,1.17,1, 1.17,1.4,2, 1.4,Inf,3))
writeRaster(burned_dMIRBI_cl, filename = "12_Burn Severity (dMIRBI)", format = "GTiff", overwrite = T)

# plots (not as output, but to get an overview of the results)
  # define breaks and colorkey
miat <- c(1,2,3,4)
myColorkey <- list(at = miat,labels = list(labels = c("Low-Moderate", "Moderate-High","High"), at = miat+0.5))
  # plot burn severity classes (dNBR)
spplot(burned_dNBR_cl, maxpixels=1000000, main='Burn severity classes (dNBR)',colorkey= myColorkey)

  # plot burn severity classes (dMIRBI)B
spplot(burned_dMIRBI_cl, maxpixels=1000000, main='Burn severity classes (dMIRBI)',colorkey= myColorkey)


################################# 6) Quantification of burned area ################################# 

# create empty data frame  
area_df <- data.frame(matrix(ncol = 2, nrow = 2))
names(area_df) <- c("total_burned_area(km2)", "total_carbon_release(mio_t)")
row.names(area_df) <- c("dNBR", "dMIRBI")

# area calculations
  # number of pixels of burned area mask 
  calc_pixels1 <- cellStats(burned_dNBR, stat='sum',na.rm=T)
  calc_pixels2 <- cellStats(burned_dMIRBI, stat='sum',na.rm=T)

  # totally burned area in km2 -> multiply with 20^2 (pixel size)
  burned_area1 <- calc_pixels1*20^2/1000000
  burned_area2 <- calc_pixels2*20^2/1000000

# estimation of total CO2 release (according to NABU 50kg/m2, with 1 m depth)
CO2_1 <- burned_area1*180/1000
CO2_2 <- burned_area2*180/1000

# fill data frame with values
area_df$`total_burned_area(km2)` <- c(burned_area1,burned_area2)
area_df$`total_carbon_release(mio_t)` <- c(CO2_1,CO2_2)

# export data frame as csv
write.csv(area_df, file = "burned_area.csv")


####################################### END OF ANALYSIS ####################################### 

