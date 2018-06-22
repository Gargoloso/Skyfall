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
  
  print('Extracting Cardelli A(Ha) extinction.')  
  path_ha <- str_c(galaxy[i],"/Cardelli_Base_Fluxes.dat") #creates path to Ha and [N II]data file for each galaxy
  data0 <- read.table(path_ha, header=TRUE)
  
  Cardelli_id <- data0$id[which(data0$AHa >= 0)]
  
  
  print('Extracting Calzetti A(Ha) extinction.')
  path_hb <- str_c(galaxy[i],"/Calzetti_Base_Fluxes.dat") #creates path to H-beta and [OIII] data file for each galaxy
  data3 <- read.table(path_hb, header=TRUE)
  
  Calzetti_id <- data0$id[which(data3$AHa >= 0)]
  
  
  #x <- match(Cardelli_id,Calzetti_id)
  
  Cardelli_AHa <- data0$AHa[which(data0$AHa >= 0 & data3$AHa >= 0 & !is.na(data0$AHa) & !is.na(data3$AHa))]
  Calzetti_AHa <- data3$AHa[which(data0$AHa >= 0 & data3$AHa >= 0 & !is.na(data0$AHa) & !is.na(data3$AHa))]
  
  plot(Cardelli_AHa, Calzetti_AHa, xaxs= "i",yaxs="i", pch=18, xlab = expression(paste("Cardelli et al. (1989)   A(H",alpha,")")), ylab = expression(paste("Calzetti et al. (2000) A(H",alpha,")")), xlim=c(0,10), ylim=c(0,10))
  #DATA <- merge(data0,data3,by.x = 1, by.y =1)
  #attach(DATA)
  
  #Saving plot.
  
  map <- str_c(galaxy[i],"/",galaxy[i],"_Comp_AHa.eps")
  dev.copy2eps(file=map)
  map <- str_c("convert -density 300 ",galaxy[i],"/",galaxy[i],"_Comp_AHa.eps ",galaxy[i],"/",galaxy[i],"_Comp_AHa.png")
  system(map)
}
  
  