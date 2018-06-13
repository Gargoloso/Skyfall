#############################################################################
#############################################################################

##   This program plots maps of luminosities for the H-alpha, [O III] 5007,##
##   and [N II] 6583 emission lines.                                       ##

##   February 25, 2018 A. Robleto-Orús                                     ##

#############################################################################
#############################################################################


##Clean the workspace.
rm(list=ls(all=TRUE))

##Libraries.
library("fields")
require("stringr")
library("png")

########################################################################

##				DATA INPUT			      ##

########################################################################

setwd("~/Rings/ringed_work/") #Directrory with our data.

data <- read.table("lis.dat", header=TRUE)
attach(data)

galaxy <- as.character(name)

####################################################################
####################################################################

##                   ITERATION FOR ALL GALAXIES                   ##

####################################################################
####################################################################

##Loop for all galaxies
for(i in 1:length(galaxy)){
  
  print('****************************')
  print("NEW GALAXY: ")
  print(galaxy[i])
  print('****************************')
  
  ####################################################################
  ####################################################################
  
  ##                EXTRACTION OF EMISSION-LINE DATA                ##
  
  ####################################################################
  ####################################################################
  
  ##Load data
  
  print('Extracting Calzetti A(Ha) extinction.')  
  path_ha <- str_c(galaxy[i],"/Calzetti_Base_Fluxes.dat") #creates path to Ha and [N II]data file for each galaxy
  data0 <- read.table(path_ha, header=TRUE)
  attach(data0)
  

  y <-  (trunc(ID[which(AHa >= 0 & AHa <= 20 & !is.na(AHa))]/100)  ) #Obtaining Y from id.
  x <-  (ID[which(AHa >= 0 & AHa <= 20 & !is.na(AHa))]-(y*100) ) #Obtaining X from id.
  
  ##id of central spaxel
  print('Determining central spaxel postion.')
  xc <- trunc(naxis1[i]/2)
  yc <- trunc(naxis2[i]/2)
  
  x <- x - xc
  y <- y - yc
  
  #x <- match(Calzetti_id,Calzetti_id)
  
  Calzetti_AHa <- AHa[which(AHa >= 0 & AHa <= 20 & !is.na(AHa))]
 
  ##########################################################################
  ##H-alpha Map
  
  print('Plotting H-alpha luminosity map.')
  
  ## Specify color map for z.
  
  rbPal <-  colorRampPalette((rainbow(100,start=0,end=0.9))) 
  Col <- rbPal(100)[as.numeric(cut(Calzetti_AHa,breaks = 100))]
  
  plot(x,y,xlim=c(-40,40),ylim=c(-40,40),xaxs= "i",yaxs="i", col="white", main = galaxy[i],xlab = expression(paste(Delta, alpha," (arcsec)")), ylab = expression(paste(Delta, delta," (arcsec)")), cex.lab=1.3, cex.axis=1.3)
  #rasterImage(ima,-40,-40,40,40,interpolate=F)
  points(x,y, col = (Col), pch=15,cex = 0.7) #Plots spaxels with color map reversed (to have greater values darker and lower lighter).
  grid() #Coordinates'grid.
  abline(h = 0, v = 0, col = "gray60",lwd=2)
  #mtext(side=3,expression(paste(" log"[10]," H",alpha," Luminosity [erg s"^-1,"]")))
  image.plot( -40, -40, legend.only=TRUE, zlim= range(Calzetti_AHa), col = ((rainbow(100,start=0,end=0.9)))) #Plots the color bar. Reverse colors here too.
  
  
  
  #Saving plot.
  
  map <- str_c(galaxy[i],"/",galaxy[i],"_Calzetti_AHa_Map.eps")
  dev.copy2eps(file=map)
  map <- str_c("convert -density 300 ",galaxy[i],"/",galaxy[i],"_Calzetti_AHa_Map.eps ",galaxy[i],"/",galaxy[i],"_Calzetti_AHa_Map.png")
  system(map)
}

