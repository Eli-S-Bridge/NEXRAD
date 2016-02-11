###Install packages if necessary.
install.packages("plyr")
install.packages("sp")
install.packages("maptools")
install.packages("raster")
install.packages("grDevices")
install.packages("rgdal")
install.packages("rgeos")
install.packages("gdata")
install.packages("reshape")
install.packages("cleangeo")
install.packages("zoom") #not needed for model, but usefull for examining plots
install.packages("ncdf4")
install.packages("tidyr")
install.packages("stringr")     # packages should all be installed

#Configure AWS. 
#Open a shell (Tools -> shell)
#type in "aws configure"
#provide the Access Key ID and the Secret Access Key. (For the other entries just hit return)

library(ncdf4)
library(stringr)

## date loops ##
## June has 30 days
## July and August have 31  



## specify file of interest
sweeps = data.frame()


Radar = "KFWS"

#for (y in 2013:2015){}
Year = 2014

for (m in 6:8){
  Month = m  
  
  for (d in 1:3){
    Day = d
    
    times = c(seq(103000,106000,700),seq(110000,116000,700),seq(120000,126000,700))
    
    for (t in times){
      Time = t
      
      ## format date-time accordingly
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
      ncdftxt <- system(paste("java -classpath toolsUI-4.6.jar ucar.nc2.FileWriter -in", fname2, "-out", fname3), intern=TRUE)
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
      
      
      if ("elevationR_HI" %in% names(ncdf1$var)==T){
        
        if(length(elevationR_HI)==361){
          
          elevationR_HI = rep(elevationR_HI, 1192)
          
          azimuthR_HI = rep(azimuthR_HI, each=1192)
          
          distanceR_HI = distanceR_HI[1:1192]
          
          Reflectivity_HI = as.data.frame(Reflectivity_HI[1:1192,1:CC_length,1])      # '1' after CC_length specifies sweep #
          Reflectivity_HI = cbind(distanceR_HI, Reflectivity_HI)
          Reflectivity_HI<-as.data.table(Reflectivity_HI)
          Reflectivity_HI<-melt(Reflectivity_HI, id=c("distanceR_HI"))
          Reflectivity_HI=as.data.frame(Reflectivity_HI)
          Reflectivity_HI<-Reflectivity_HI[,c(1,3)] 
          Reflectivity_HI<-cbind(Reflectivity_HI, azimuthR_HI, elevationR_HI)
          Reflectivity_HI=subset(Reflectivity_HI,elevationR_HI!="NaN")
          remove(distanceR_HI, azimuthR_HI, elevationR_HI) 
        } else if (length(elevationR_HI)==720){
          
          elevationR_HI = rep(elevationR_HI, 1192)
          
          azimuthR_HI = rep(azimuthR_HI, each=1192)
          
          distanceR_HI = distanceR_HI[[1]][1:1192]
          
          Reflectivity_HI = as.data.frame(Reflectivity_HI[1:1192,1:CC_length,1])      #sweep ref
          Reflectivity_HI = cbind(distanceR_HI, Reflectivity_HI)
          Reflectivity_HI = as.data.table(Reflectivity_HI)
          Reflectivity_HI = melt(Reflectivity_HI, id=c("distanceR_HI"))
          Reflectivity_HI = as.data.frame(Reflectivity_HI)
          Reflectivity_HI = Reflectivity_HI[,c(1,3)]
          
          Reflectivity_HI = cbind(Reflectivity_HI, azimuthR_HI, elevationR_HI)
          Reflectivity_HI = subset(Reflectivity_HI,elevationR_HI!="NaN")
          remove(distanceR_HI, azimuthR_HI, elevationR_HI)
          
          Reflectivity_HI = Reflectivity_HI[order(Reflectivity_HI$azimuthR_HI),] 
          Reflectivity_HI$azimuthR_HI = sort(rep(seq(0.5,360,1),length(unique(Reflectivity_HI$distanceR_HI))),decreasing = FALSE)
        }  else if (length(elevationR_HI) > 720){
          
          elevationR_HI = rep(elevationR_HI, 1192)
          
          azimuthR_HI = rep(azimuthR_HI, each=1192)
          
          distanceR_HI = distanceR_HI[[1]][1:1192]
          
          Reflectivity_HI<-as.data.frame(Reflectivity_HI[1:1192,1:CC_length,1])     # nb: third bracket term is sweep elevation
          Reflectivity_HI<-cbind(distanceR_HI, Reflectivity_HI)
          Reflectivity_HI<-as.data.table(Reflectivity_HI)
          Reflectivity_HI<-melt(Reflectivity_HI, id=c("distanceR_HI"))
          Reflectivity_HI=as.data.frame(Reflectivity_HI)
          Reflectivity_HI<-Reflectivity_HI[,c(1,3)]
          
          Reflectivity_HI<-cbind(Reflectivity_HI, azimuthR_HI, elevationR_HI)
          Reflectivity_HI=subset(Reflectivity_HI,elevationR_HI!="NaN")
          remove(distanceR_HI, azimuthR_HI, elevationR_HI)
          
          Reflectivity_HI=Reflectivity_HI[order(Reflectivity_HI$azimuthR_HI),] 
          Reflectivity_HI$azimuthR_HI=sort(rep(seq(0.5,360,1),length(unique(Reflectivity_HI$distanceR_HI))),decreasing = FALSE)
        }  else {
          next
        }
        
      } else {
        
        elevationR_HI = rep(elevationR_HI, 1192)
        
        azimuthR_HI = rep(azimuthR_HI, each=1192)
        
        distanceR_HI = distanceR_HI[[1]][1:1192]
        
        CC_length = length(azimuthR_HI)
        
        Reflectivity_HI<-as.data.frame(Reflectivity[1:1192,1:CC_length,1])
        Reflectivity_HI<-cbind(distanceR_HI, Reflectivity_HI)
        Reflectivity_HI<-as.data.table(Reflectivity_HI)
        Reflectivity_HI<-melt(Reflectivity_HI, id=c("distanceR_HI"))
        Reflectivity_HI=as.data.frame(Reflectivity_HI)
        Reflectivity_HI<-Reflectivity_HI[,c(1,3)]
        
        Reflectivity_HI<-cbind(Reflectivity_HI, azimuthR_HI, elevationR_HI)
        remove(distanceR_HI, azimuthR_HI, elevationR_HI)
        
      }
      
      #########################################################################################################
      ###################################### the math #########################################################
      #########################################################################################################
      
      ### roost domain range boundaries here ###
      Reflectivity_HI=subset(Reflectivity_HI, distanceR_HI>40000 & distanceR_HI<145000)
      ### Azimuthal boundaries ###
      Reflectivity_HI=subset(Reflectivity_HI, azimuthR_HI>30 & azimuthR_HI<85) 
      
      Reflectivity_HI=subset(Reflectivity_HI, value < 35)
      Reflectivity_HI=subset(Reflectivity_HI, value >= -33)
      
      beta=10*log10((10^3*(pi^5)*(0.93))/ ((10.7)^4))
      
      ###cross-section is 13 cm^2 Eastwood 1967 ###
      birds=round(sum((10^((Reflectivity_HI$value+beta)/10)/13) * (((0.35)*(sqrt(2*pi))/(2*log(2)) *(pi*(Reflectivity_HI$distanceR^2)*250 *(rad(0.96)^2))/4)*10^-9)),2)
      ## nb not all radar sweeps have ranges gates of 250 m
      #Change to linear
      Reflectivity_HI$LinearZ= 10^(Reflectivity_HI$value/10)
      
      #aggregate azimuths to generate uniform 1-degree resolution for all sweeps
      #low elevations are 0.5 degree resolution and upper tilts are 1 degree
      
      Reflectivity_HI=as.data.table(Reflectivity_HI)
      
      Reflectivity_HI=Reflectivity_HI[,list(reflectivity = mean(LinearZ)), by = list(Reflectivity_HI$distanceR_HI,Reflectivity_HI$elevationR_HI,Reflectivity_HI$azimuthR_HI)]
      setnames(Reflectivity_HI, c("distanceR","elevationR", "azimuthR", "LinearZ"))
      
      Reflectivity_HI$dBZ= 10*log10(Reflectivity_HI$LinearZ)
      
      elevation=mean(Reflectivity_HI$elevationR)
      MeanZ=mean(Reflectivity_HI$LinearZ)
      SumZ=sum(Reflectivity_HI$LinearZ)
      SDZ=sd(Reflectivity_HI$LinearZ)
      Volumes=length(Reflectivity_HI$LinearZ)
      
      Reflectivity_HI=subset(Reflectivity_HI,dBZ!=-33)
      
      MeanZ_m33=mean(Reflectivity_HI$LinearZ)
      SumZ_m33=sum(Reflectivity_HI$LinearZ)
      SDZ_m33=sd(Reflectivity_HI$LinearZ)
      Volumes_m33=length(Reflectivity_HI$LinearZ)
      
      Reflectivity_HI=cbind(fname2, elevation, birds,MeanZ,SumZ,SDZ, Volumes, MeanZ_m33,SumZ_m33,SDZ_m33, Volumes_m33)
      sweeps=rbind(Reflectivity_HI,sweeps)
      print(fname2)
      
      
      file.remove(fname)
      file.remove(fname1)
      file.remove(fname2)
      file.remove(fname3)
    }}}

#filename="KFWS_reflect.csv"
#setwd("C:/stats/Radar.mat/KFWS/2014")
#write.csv(sweeps,filename, row.names=F )

#############################################################################

