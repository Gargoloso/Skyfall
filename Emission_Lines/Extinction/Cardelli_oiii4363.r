#############################################################################
#############################################################################

##   This calculates extinction corrected fluxes for the [O III] 4363      ##
##   emission line using the CCM 89 extinction law.                        ##    

##   May 15, 2018 A. Robleto-Orús                                          ##
##   June 10, 2018 Optimized.                                              ##

#############################################################################
#############################################################################

#############################################################################

## WARNING!!!!!                                                            ##
## The Cardelli_Base_Fluxes.r script must be runned before this one!!!     ##

#############################################################################

##Clean the workspace.
rm(list=ls(all=TRUE))

##Libraries.
require("stringr")

##################################

##				DATA INPUT			      ##

##################################

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
  
##Load data.
  
print('Extracting [O III] 4363  data.')  
path_oi <- str_c(galaxy[i],"/lis_hgoiii4363.res") #creates path to the [O III] 4363 data file for each galaxy.
data0 <- read.table(path_oi, header=TRUE)
  
print('Extracting H-beta  data.')
path_hb <- str_c(galaxy[i],"/lis_hb.res") #creates path to H-beta and [OIII] 5007 data file for each galaxy.
data1 <- read.table(path_hb, header=TRUE)
  
print('Extracting dereddened H-beta and Cbeta data') #Creates path to a file with dereddened H-beta fluxes and Cbeta normalization constant, from a previous run of Cardelli_Base_Fluxes.dat
pathcb <- str_c(galaxy[i],"/Cardelli_Base_Fluxes.dat")
data2 <- read.table(pathcb, header=TRUE) 

##Merge data.
DATA0 <- merge(data0, data1, by.x = 1, by.y = 1)
DATA  <- merge(DATA0, data2, by.x = 1, by.y = 1)
attach(DATA)
  
  
##Extracting coordinates from ID.
  
print('Extracting spaxel coordinates.')
id <- id[which(fluxhb!=500 & fluxoiii4363!=500)]

##Extracting line surface specific intensities [1e-16 erg cm^-2 s^-1 A^-1 arcsec^-2].

print('Extracting emission lines.')

#H-alpha 6563
print('Extracting H-beta.')
fb <- fluxhb[which(fluxhb!=500 & fluxoiii4363!=500)] #Obtaining the H-beta surface specific intensity.

#[N II] 6583
print('Extracting [O III] 4363.')
foiii4363 <- fluxoiii4363[which(fluxhb!=500 & fluxoiii4363!=500)] #Obtaining the [O III] 4363 surface specific intensity.

#H-beta dereddened flux.
print('Extracting H-beta dereddened flux.')
Fb <- Fb[which(fluxhb!=500 & fluxoiii4363!=500)] #Obtaining the dereddened H-beta surface specific flux.

#Cbeta
print('Extracting H-beta normalization constant.')
Cbeta <- Cbeta[which(fluxhb!=500 & fluxoiii4363!=500)] #Obtaining the H-beta normalization constant for the extinction law.

##Convert from surface specific intensity to standard specific intensity [erg cm^-2 s^-1 A^-1 sr^-1]
#print('Converting surface specific units to specific intensities.')
#foiii4363   <-(1e-16) * foiii4363 / 2.3504e-11 
#fb        <-(1e-16) * fb / 2.3504e-11 

##Multiply all values by 1e-16 (CALIFA reference level) Comment this section if using the previous conversion to specific intensity.

fb   <-(1e-16)*fb
foiii4363  <-(1e-16)*foiii4363


#####################################################

##                    DEREDDENING                 ##  

#####################################################

print('Dereddening...')

Rv <- 3.1
l <-  c(4861,4363) #Lines we are going to unredden, in Angstrom.
l <- l*1e-4 #Convert Angstrom to micrometers.
xl <- 1/l
yl <- xl - 1.82

a_x <- 1 + (0.17699*yl) - (0.50447*yl^2) - (0.02427*yl^3) + (0.72085*yl^4) + (0.01979*yl^5) - (0.77530*yl^6) + (0.32999*yl^7)
b_x <- (1.41338*yl) + (2.28305*yl^2) + (1.07233*yl^3) -(5.38434*yl^4) - (0.62251*yl^5) + (5.30260*yl^6) - (2.09002*yl^7)


rIoiii4363 <- seq(1:length(id))
Foiii4363  <- seq(1:length(id))

##Claculate the extinction law.

for(j in 1:length(id)){
  
  f_l <- seq(1:length(l))  #CCM89 extinction law, normalized to V-band.
  fb_l <- seq(1:length(l)) #Extinction law normalized to H-beta.
  
  for(k in 1:length(l)){
    f_l[k] <- a_x[k] + (b_x[k] / Rv)
  }
  fb_l <- f_l/f_l[1]
  
  ##Deredden each line normalized to I_H-beta.
  
  rIoiii4363[j] <- ((foiii4363[j])/fb[j])*10^(Cbeta[j]*(fb_l[2]-1)) #Dereddened H-beta normalized [O III] 4363 surface specific flux.
  Foiii4363[j]  <- rIoiii4363[j]*Fb[j] #Dereddened [O I] 6300 surface specific flux.
}

##############################################################################
##Save data to files.

print('Saving fluxes to data file.')
resume <- data.frame(id, Foiii4363)
tabla <- str_c(galaxy[i],"/Cardelli_oiii4363_Fluxes.dat")
write.table(resume, tabla, sep="\t",quote=FALSE)

}