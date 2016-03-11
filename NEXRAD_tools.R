#DOWNLOAD AND OPEN A NEXRAD VOLUME SCAN
#REQUIRES Radar, Year, Month, Day, Time, ragengates, mindist, maxdist

library(ncdf4)
library(ncdf4.helpers)
library(stringr)
library(data.table)
library(plyr)
library(rgdal)
library(raster)
library(sp)
library(geosphere)
library(dismo)
library(rgeos)
library(snow)
library(aspace)
library(abind) 
library(FNN)  
#
# Radar = "KFWS"
# Year = 2015
# Month = 7
# Day = 4
# Time = 112200
# 
# rangegates <- 500 #how many range gates are to be used (assuming radar data start at distance of 0)
# mindist <- 2000   #Minimum distance to process radar data
# maxdist <- 150000 #Maximum distance to process radar data. (we will use the minimum of the ranges defined by this number and the number of range gates)

getRadFileName <- function(fname) {
  Radar <- substr(fname, 1, 4)
  Year <- substr(fname, 5, 8)
  Month <- substr(fname, 9, 10)
  Day <- substr(fname, 11, 12)
  Time <- substr(fname, 14, 19)
  fpath <- paste(Year, Month, Day, Radar, sep = "/")
  
  flist <- system(paste("aws s3 ls --recursive s3://noaa-nexrad-level2", fpath, sep = "/"),
                  intern=TRUE)                #sends command to the shell which creates a list of all the files for the day specified above.
  
  #### Find closest timestamp within the files on the day specified ####
  flist <- as.data.frame(flist, stringsAsFactors = F)            #changes above list into data frame without converting the names to factors
  fext <- str_sub(flist$flist, start = -2, end = -1)
  flist <- data.frame(full = flist[fext == "gz",])               #subsets only those files with the extention '.gz'
  flist$name <- sub(pattern = paste(".*", fpath, "/", sep = ""), replacement = "", flist$full)
  flist$time <-  as.numeric(str_sub(flist$name, start = 14, end = 19))
  time.no <- as.numeric(Time)
  diffs <- (abs(flist$time - time.no))
  mindif <- min(diffs)
  if(mindif > 600){
    print("WARNING: GOOD TIME MATCH NOT FOUND AMONG AVAILABLE NEXRAD FILES")
  }
  closest <- which.min(abs(flist$time - time.no)) 
  Time <- str_pad(flist$time[closest], 6, pad = "0")      # ID closest timestamp
  cat("closest file time is:", Time, "\n")                # print available timestamp
  
  fname <- flist$name[closest]
  return(fname)
}
  
getRadNCDF <- function(fname){
  Radar <- substr(fname, 1, 4)
  Year <- substr(fname, 5, 8)
  Month <- substr(fname, 9, 10)
  Day <- substr(fname, 11, 12)
  Time <- substr(fname, 14, 19)
  fname1 <- paste(Radar, Year, Month, Day, "_", Time, ".gz", sep = "") 
  fname2 <- paste(Radar, Year, Month, Day, "_", Time, sep = "") 
  fname3 <- paste(fname2, "ncdf", sep = '.')
  fpath <- paste("https://noaa-nexrad-level2.s3.amazonaws.com", Year, Month, Day, Radar, fname, sep ="/")
  
  #### Get file ####
  if(!file.exists(fname2)) {
    download.file(fpath, destfile = fname1)                 # download file
    system(paste("gunzip", fname1))                         # unzip radar file
  }
  #You have to maker sure java is installed and that you have this toolsUI thing in an appropriate directory
  ncdftxt <- system(paste("java -classpath /home/rstudio/toolsUI-4.6.jar ucar.nc2.FileWriter -in", fname2, "-out", fname3), intern=TRUE)
  ncdf1 <- nc_open(fname3, write = T) # name ncdf file
  return(ncdf1)
}

extractScan <- function(ncdf1, rangegates, mindist, maxdist, scan, Radar) {
  #### Get variables ####
  refvars <- c("Reflectivity", "distanceR", "elevationR", "azimuthR")
  varlist <- nc.get.variable.list(ncdf1, min.dims = 1)
  if(scan < 1.75 & any(varlist=="Reflectivity_HI")) refvars <- paste(refvars, "HI", sep="_")
  
  reflectivity <- ncvar_get(ncdf1, refvars[1])              # assign reflectivity to a vector
  rangegates <- min(dim(reflectivity)[1], rangegates)
  reflectivity <- reflectivity[1:rangegates,,]
  reflectivity[is.na(reflectivity)] <- -33                # replace 'NA's w/ '-33'
  
  elevation <- ncvar_get(ncdf1,refvars[3])                    # grab elevation
  escan <- min(which(round_any(elevation[10,], accuracy = 0.5)==scan)) # which elevation set best matches the scan
  elevation <- as.vector(elevation[,escan])
  
  distance <- as.vector(ncvar_get(ncdf1,refvars[2])) #grab distance (i.e. range gates)
  distance <- distance[1:rangegates]
  gatesize <- distance[3]-distance[2]
  
  azimuth <-   ncvar_get(ncdf1,refvars[4])       #grab azimuth
  azimuth <-  as.vector(azimuth[,escan])
  
  if(length(azimuth) %% 360 > 1) { #some scans have extra azimuths to them, much reduce them down to a multiple of 360
    azimuth <- azimuth[!is.nan(azimuth)] #usually some Nan's at the end
    extra <- length(azimuth) %% 360
    minind<- extra+1
    maxind <- extra+360
    azimuth <- azimuth[minind:maxind]
    elevation <- elevation[minind:maxind]
    reflectivity <- reflectivity[,minind:maxind,]
  }
  
  CC_length <- length(azimuth)
  
  elevation = rep(elevation, rangegates)          #repeat the elevation vector to generate a value for each range gate
  azimuth = rep(azimuth, each=rangegates)         #same for azimuth
  
  reflectivity <- as.data.frame(reflectivity[1:rangegates,1:CC_length,escan])      #just the low sweep and only the first rangegates range gates
  reflectivity$distance  = distance                     #distance appended to each range gate
  reflectivity = as.data.table(reflectivity)                             
  reflectivity = melt(reflectivity, id=c("distance"))             #melt to three columns 
  reflectivity = as.data.frame(reflectivity)                            
  reflectivity = reflectivity[,c(1,3)]                                  #remove useless column
  reflectivity = cbind(reflectivity, azimuth, elevation)        #add on columns for azimuth and elevation
  reflectivity = subset(reflectivity,elevation!="NaN")              #remove missing data
  remove(distance, azimuth, elevation)                            #scrap bulky variables
  reflectivity = reflectivity[order(reflectivity$azimuth),]      #sort data by azimuth
  reflectivity$azi = round_any(reflectivity$azimuth, accuracy = 0.25)
  reflectivity <- reflectivity[order(reflectivity[,1], reflectivity[,5]),]
  reflectivity <- subset(reflectivity, reflectivity$distance >= mindist & reflectivity$distance <= maxdist)

  Ref <- reflectivity
  
  #set working directory to landcover
  #setwd("/data/landcover")
  site <- read.csv("nexrad_site_list_with_utm.csv") #NEXRAD data file
  site <- subset(site, SITE==substring(Radar, 2)) #SELECT A RADAR HERE
  
  #Define some radar parameters
  BASE=site$BASE_HT_M
  TOWERHT=site$TOWER_HT_M
  LAT=site$LATITUDE_N
  LON=site$LONGITUDE_W
  coords <- matrix(data= c(LON,LAT), nrow=1) #two column matrix of locations - one row only
  #A few constants
  ADJUSTED  = 4/3                   ## 4/3 correction for standard refraction
  EARTHE    = 6378137   	    	      ## Earth equatorial radius (m)
  ANTHT = (TOWERHT)                 ## radar antenna elevation 
  EARTH = (EARTHE*(1 - (0.0033493*(sin_d(LAT))*2)))  #calculate Earth radius at radar latitude
  EARTHRAD= (EARTH+ANTHT)
  #define the minimum and maximum ranges and the size of each range bin 
  minrange <- min(Ref$distance)      #minimum range (hole in the middle)
  maxrange <- max(Ref$distance)      #outermost measurement range
  #sector length - approximate asmuthal length of each segment
  
  #scans <- seq(0.5, 4.5, by=1)  #vector of scan elevations
  
  #Find coordinates for each range bin

  #adjust the ranges to account for beam angle (matters verly little for low scans)
  Ref$centerbeamHts  <- sqrt(( ADJUSTED*EARTHRAD )^2 + (Ref$distance)^2 + 2*( ADJUSTED*EARTHRAD )*
                               (Ref$distance)*sin_d(scan)) - ADJUSTED*EARTHRAD + ANTHT
  Ref$elevAdj <- atan_d((Ref$centerbeamHts-ANTHT)/(Ref$distance))
  Ref$adjDist <- cos_d(Ref$elevAdj)*Ref$distance
  Ref <- subset(Ref, select=-c(centerbeamHts,elevAdj))
  pcoords <- destPoint(matrix(c(LON,LAT), nrow = nrow(Ref), ncol = 2, byrow=TRUE), b=Ref$azi, d=Ref$distance)
  Ref <- cbind(Ref, pcoords)
}

NEXRast <- function(Ref, rastext, Rastres, Radar, hole = T) {
  site <- read.csv("nexrad_site_list_with_utm.csv") #NEXRAD data file
  site <- subset(site, SITE==substring(Radar, 2)) #SELECT A RADAR HERE
  LAT=site$LATITUDE_N
  LON=site$LONGITUDE_W
  coords <- matrix(data= c(LON,LAT), nrow=1) #two column matrix of locations - one row only
  pcoords <- subset(Ref, select= c(lon,lat))
  refRast <- raster(x=rastext, nrows=Rastres, ncols=Rastres, crs = "+proj=longlat +datum=WGS84" )
  refRast <- rasterize(x=Ref[,c("lon", "lat")], y=refRast, field=Ref$value, fun = max)
  refRast@data@values[( ((1:Rastres^2 %/% Rastres)- Rastres/2)^2 + ((1:Rastres^2 %% Rastres) - Rastres/2)^2) >= (Rastres/2)^2] <- -99 #The logic is that the criteria are of the form (x-dx)^2+(y-dy)^2 > r^2 where dx and dy are the center coordinates of the circle and r is the radius 
  if(hole) {
    inCircRatio <- as.numeric(distm(coords, pcoords[1,], fun=distCosine) / distm(coords, pcoords[nrow(pcoords),], fun=distCosine)) #Ratio between the size of the inner circle and the outer circle
    refRast@data@values[( ((1:Rastres^2 %/% Rastres)- Rastres/2)^2 + ((1:Rastres^2 %% Rastres) - Rastres/2)^2) < (Rastres*inCircRatio*0.5)^2] <- -99 #The logic is that the criteria are of the form (x-dx)^2+(y-dy)^2 > r^2 where dx and dy are the center coordinates of the circle and r is the radius 
  }
  missing <- which(is.na(refRast[,]))
  xy<-xyFromCell(refRast,missing)
  rindex <- knnx.index(data=pcoords, query=xy, k=1)
  repl <- Ref$value[rindex]
  refRast[missing] <- repl
  refRast[refRast[,] == -99] <- NA
  return(refRast)
}





