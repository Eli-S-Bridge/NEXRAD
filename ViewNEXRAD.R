#View a composite of the 0.5 and 1.5 degree radar scans

library(raster)

source("NEXRAD_tools.R")

Radar = "KMOB"
Year = 2015
Month = 8
Day = 31
Time = 112650

Time <- str_pad(Time, 6, pad = "0")
Month <- str_pad(Month, 2, pad = "0")
Day <- str_pad(Day, 2, pad = "0")

rangegates <- 500 #how many range gates are to be used (assuming radar data start at distance of 0)
mindist <- 2000   #Minimum distance to process radar data
maxdist <- 150000 #Maximum distance to process radar data. (we will use the minimum of the ranges defined by this number and the number of range gates)

fname <- paste(Radar, Year, Month, Day, "_", Time, sep ="")
fname <- getRadFileName(fname)

ncdf1 <- getRadNCDF(fname)

#Get the 0.5 degree scan
Ref1 <- extractScan(ncdf1, rangegates, mindist, maxdist, 0.5, Radar)

#define an extent for a raster image
rastext <- extent(min(Ref1$lon), max(Ref1$lon), min(Ref1$lat), max(Ref1$lat))
Rast1 <- NEXRast(Ref1, rastext, Rastres=1000, Radar) 
# cuts <- c(-100, seq(-30, 75, by=5), 200) #set breaks for color scheme
# pal <- c("black", "azure", "mistyrose2", "plum", "plum3", "seashell2","seashell3","seashell4","darkturquoise","deepskyblue","blue3","chartreuse1","chartreuse3","chartreuse4","gold","gold3", "darkorange","red","firebrick3","firebrick4", "deeppink1", "mediumorchid3","gray95")
# image(Rast1, axes = T, breaks = cuts, col = pal)


#Get the 1.5 degree scan
Ref2 <-  extractScan(ncdf1, rangegates, mindist, maxdist, 1.5, Radar)
Rast2 <- NEXRast(Ref2, rastext, Rastres=1000, Radar) #use same extent and resolution as before
#image(Rast2, axes = T, breaks = cuts, col = pal)

Rast3 <- overlay(Rast1, Rast2, fun = "max")

try(dev.off)
cuts <- c(-100, seq(-30, 75, by=5), 200) #set breaks for color scheme
pal <- c("black", "azure", "mistyrose2", "plum", "plum3", "seashell2","seashell3","seashell4","darkturquoise","deepskyblue","blue3","chartreuse1","chartreuse3","chartreuse4","gold","gold3", "darkorange","red","firebrick3","firebrick4", "deeppink1", "mediumorchid3","gray95")
image(Rast3, axes = T, breaks = cuts, col = pal)


