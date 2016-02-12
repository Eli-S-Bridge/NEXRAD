# ###Install packages if necessary.
# install.packages("plyr")
# install.packages("sp")
# install.packages("maptools")
# install.packages("raster")
# install.packages("grDevices")
# install.packages("rgdal")
# install.packages("rgeos")
# install.packages("gdata")
# install.packages("reshape")
# install.packages("cleangeo")
# install.packages("zoom") #not needed for model, but usefull for examining plots
# install.packages("ncdf4")
# install.packages("tidyr")
# install.packages("data.table")
# install.packages("stringr")     # packages should all be installed

#Configure AWS. 
#Open a shell (Tools -> shell)
#type in "aws configure"
#provide the Access Key ID and the Secret Access Key. (For the other entries just hit return)

library(ncdf4)
library(stringr)
library(data.table)
## date loops ##
## June has 30 days
## July and August have 31  



## specify where data will be stored
sweeps = data.frame()


Radar = "KFWS"
Year = 2014
Month = 6
Day = 12
Time = 113000
      
Time <- str_pad(Time, 6, pad = "0")
Month <- str_pad(Month, 2, pad = "0")
Day <- str_pad(Day, 2, pad = "0")
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
fname1 <- paste(Radar, Year, Month, Day, "_", Time, ".gz", sep = "") 
fname2 <- paste(Radar, Year, Month, Day, "_", Time, sep = "") 
fname3 <- paste(fname2, "ncdf", sep = '.')
fpath <- paste("https://noaa-nexrad-level2.s3.amazonaws.com", Year, Month, Day, Radar, fname, sep ="/")

#### Get file ####
download.file(fpath, destfile = fname1)                 # download file
system(paste("gunzip", fname1))                         # unzip radar file
ncdftxt <- system(paste("java -classpath /home/rstudio/toolsUI-4.6.jar ucar.nc2.FileWriter -in", fname2, "-out", fname3), intern=TRUE)
ncdf1 <- nc_open(fname3)                                # name ncdf file

#### Get variables ####
Reflectivity_HI <- ncvar_get(ncdf1, "Reflectivity_HI")              # assign reflectivity to a vector
Reflectivity_HI[is.na(Reflectivity_HI)] <- -33                # replace 'NA's w/ '-33'

distanceR_HI <- as.data.frame(ncvar_get(ncdf1,"distanceR_HI")) #grab distance (i.e. range gates)

elevationR_HI <- ncvar_get(ncdf1,"elevationR_HI")    # grab elevation
elevationR_HI <- as.vector(elevationR_HI[,1])     # get elevation angle of interest

azimuthR_HI <-   ncvar_get(ncdf1,"azimuthR_HI")       #grab azimuth
azimuthR_HI <-  as.vector(azimuthR_HI[,1])

CC_length <- length(azimuthR_HI)
      
      
      
      ######################################################################################################################################
      ####################################### the grab #####################################################################################
      ######################################################################################################################################
      
      
     # if ("elevationR_HI" %in% names(ncdf1$var)==T){ #does this ever fail???

          
          elevationR_HI = rep(elevationR_HI, 480)  #repeat the elevation vector to generate a value for each range gate
          azimuthR_HI = rep(azimuthR_HI, each=480) #same for azimuth
          distanceR_HI = distanceR_HI[[1]][1:480]  #single vector for distance
          
          Reflectivity_HI <- as.data.frame(Reflectivity_HI[1:480,1:CC_length,1])      #just the low sweep and only the first 480 range gates
          Reflectivity_HI = cbind(distanceR_HI, Reflectivity_HI)                      #distance appended to each range gate
          Reflectivity_HI = as.data.table(Reflectivity_HI)                             
          Reflectivity_HI = melt(Reflectivity_HI, id=c("distanceR_HI"))               #melt to three columns 
          Reflectivity_HI = as.data.frame(Reflectivity_HI)                            
          Reflectivity_HI = Reflectivity_HI[,c(1,3)]                                  #remove useless column
          Reflectivity_HI = cbind(Reflectivity_HI, azimuthR_HI, elevationR_HI)        #add on columns for azimuth and elevation
          Reflectivity_HI = subset(Reflectivity_HI,elevationR_HI!="NaN")              #remove missing data
          remove(distanceR_HI, azimuthR_HI, elevationR_HI)                            #scrap bulky variables
          Reflectivity_HI = Reflectivity_HI[order(Reflectivity_HI$azimuthR_HI),]      #sort data by azimuth
          Reflectivity_HI$azimuth = sort(rep(seq(0.5,360,360/CC_length),length(unique(Reflectivity_HI$distanceR_HI))),decreasing = FALSE) #replace actual azimuths with rounded values.


#Draw polygons for x and y of voxels:
          
#Get landcover data within a radial grid that corresponds to radar voxels.

library(rgdal)
library(raster)
library(sp)
library(geosphere)
library(dismo)
library(rgeos)
library(plyr)
library(snow)
library(aspace)
library(abind)          
          
      

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
minrange <- min(Reflectivity_HI$distanceR_HI-125)       #minimum range (hole in the middle) - take desired number and subtract 125
maxrange <- max(Reflectivity_HI$distanceR_HI+125)          #outermost measurement range
secr <- 10000              #sector range - radial distance for each sector
secl <- 10000               #sector length - approximate asmuthal length of each segment
num_azi <- CC_length  #how many azimuths

rangebins <- 250            #range bin length

#scans <- seq(0.5, 4.5, by=1)  #vector of scan elevations
scans <-0.5 #lowest scan only

nbins <- ceiling((maxrange-minrange)/rangebins) #number of range bins (number of concentric rings).
nrings <- nbins+1

plist <- list(NULL) #empty list to hold data

# #list of factors for the number 360
# div <- seq_len(abs(360))
# factors <- div[360 %% div == 0L]

for (e in scans){
  #adjust the ranges to account for beam angle (matters verly little for low scans)
  print(paste("elevation ", e))
  ranges <- seq(from = minrange, to = maxrange, by = 250)
  centerbeamHts  <- sqrt(( ADJUSTED*EARTHRAD )^2 + (ranges+125)^2 + 2*( ADJUSTED*EARTHRAD )*
                           (ranges+125)*sin_d(e)) - ADJUSTED*EARTHRAD + ANTHT
  eadj <- atan_d((centerbeamHts-ANTHT)/(ranges+125))
  adj_ranges <- cos_d(eadj)*ranges
  rm(centerbeamHts,eadj)
  
  circsx <- circsy <- matrix(NA, ncol = nrings, nrow = CC_length+1)
  
  #loop through each bin.
  for(c in 1:nrings) { 
    circ <- circles(p=coords, d=adj_ranges[c], lonlat=T, n=CC_length, dissolve=T)
    cxy <- circ@polygons@polygons[[1]]@Polygons[[1]]@coords
    circsx[,c] <- cxy[,1]
    circsy[,c] <- cxy[,2]
  }
  
  circsx2 <- cbind(circsx[,-1], circsx[,1])    #shift the matrixes to position the x coordinates for each shape above each other in an array.
  circsx3 <- rbind(circsx2[-1,], circsx2[1,])
  circsx4 <- rbind(circsx[-1,], circsx[1,])
  circsx5 <- circsx
 
  
  circsy2 <- cbind(circsy[,-1], circsy[,1])    #shift the matrixes to position the x coordinates for each shape above each other in an array.
  circsy3 <- rbind(circsy2[-1,], circsy2[1,])
  circsy4 <- rbind(circsy[-1,], circsy[1,])
  circsy5 <- circsy
  
  circsxa <- abind(circsx,circsx2,circsx3,circsx4,circsx5, along = 3)
  circsya <- abind(circsy,circsy2,circsy3,circsy4,circsy5, along = 3)
  
  circsa <- abind(circsx,circsx2,circsx3,circsx4,circsx5,circsy,circsy2,circsy3,circsy4,circsy5, along = 3)
  
  tst <- (alply(circsxa,c(1,2)), 
  
  
 tm <- 
  alply(tm,3)
  
  
 A1 <- matrix(runif(12),4,3)
 A2 <- matrix(runif(12),4,3)
 A3 <- matrix(runif(12),4,3)
 
 MyList  <- list(A1,A2,A3)
 
 MyArray <- array(1:36,c(4,3,3))
  
 alply(MyArray,c(1,2))
  
  
  
  
  
  for (r in 1:(nbins)) {
    #for (r in 0:3) {
    #establish inner and outer ring distances
    rinner <- minrange + (r-1)*rangebins #inner range
    router <- minrange + (r)*rangebins #outer range
    rmid <- mean(c(rinner,router)) #middle of range
    # get beam height for given range
    centerbeamHt  <- sqrt(( ADJUSTED*EARTHRAD )^2 + rmid^2 + 2*( ADJUSTED*EARTHRAD )*
                            rmid*sin_d(e)) - ADJUSTED*EARTHRAD + ANTHT
    #calculate the angle based on the center beam height
    eadj <- atan_d((centerbeamHt-ANTHT)/rmid)
    
    #calculate more precise estimates for ground coordinates (range limits) based on angle from beam height.
    
    routeradj <- cos_d(eadj)*router
    rinneradj <- cos_d(eadj)*rinner
    
    
    
    
    
    
    rmidadj <- mean(c(rinneradj,routeradj))
    #Calculate number of wedges for the current ring
    wdeg <- 360/num_azi #wedge angle width 
    print(paste("ring #", r))
    for (i in 1:(num_azi)) {
      #for (i in 0:5) {  
      w <- (i-1)*wdeg
      w2 <- w+wdeg #set the two angles that define the edge
      midw <- mean(c(w,w2))
      
      #determine four corners polygon that defines a vpr wedge based on a radar location
      pt1 <- destPoint(coords,w,rinneradj)
      pt2 <- destPoint(coords,w,routeradj)
      pt3 <- destPoint(coords,w2,routeradj) 
      pt4 <- destPoint(coords,w2,rinneradj)
      
      #make a matrix
      cds1 <- rbind(pt1,pt2,pt3,pt4,pt1) #repeat last point to close the polygon
      #step through creation of spatial polygons object
      nextp <- ifelse(is.null(plist[[1]]), 1, length(plist)+1)
      #name includes the middle azimuth, the elevation, and the midrange
      pname <- paste(nextp, w, w2, rinneradj, routeradj, e, sep = ",")
      plist[nextp] <- Polygons(list(Polygon(cds1)), pname) #compile a long list of polygon data
    }
  }
}  

sps1 = SpatialPolygons(plist) #make the list of polygon data into a spatial polygon
proj4string(sps1) <- "+init=epsg:4202 +proj=longlat +ellps=aust_SA +no_defs" #establish the coordinate system and projection
sps1 = spTransform(sps1, CRS(proj4string(lcmap))) ## Transform projecton to WGS84
comment(sps1) <- as.character(site$SITE) #attach the site name to the spatial object

plot(sps1)

#save the big polygon list.
savename <- paste(site$SITE, secr, secl, "sps1.RData", sep ="_")
save(sps1,file=savename)
load(savename)

