#############################################################################
#############################################################################

##   This calculates extinction corrected fluxes for the [SII] 6716,6731   ##
##   emission line doublet using the Calzetti et al. (2000) extinction law.##    

##  June 10, 2018 A. Robleto-Orús                                          ##

#############################################################################
#############################################################################

#############################################################################

## WARNING!!!!!                                                            ##
## The Cardelli_Base_Fluxes.r script must be executed before this one!!!   ##

#############################################################################

##Clean the workspace
rm(list=ls(all=TRUE))

##Libraries
require("stringr")

########################################################################

##				DATA INPUT			      ##

########################################################################

setwd("~/Rings/ringed_work/") #Directrory with our data

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
  
print('Extracting [S II] 6716,6731  data.')  
path_sii <- str_c(galaxy[i],"/lis_sii.res") #creates path to [Ha and [S II] 6716,6731 data file for each galaxy. 
data0 <- read.table(path_sii , header=TRUE)
  
print('Extracting Base data') #Creates path to a file with Base results.
pathcb <- str_c(galaxy[i],"/Calzetti_Base_Fluxes.dat")
data1 <- read.table(pathcb, header=TRUE) 

##Merge data
DATA <- merge(data0, data1, by.x = 1, by.y = 1)
attach(DATA)
  
  
##Extracting coordinates from ID

  
print('Extracting spaxel ID')
id <-    id[which(fluxsii1!=500 & fluxsii2!=500)]

##Extracting line surface specific fluxes[1e-16 erg cm^-2 s^-1 A^-1 arcsec^-2]

print('Extracting emission lines.')

#[S II] 6716
print('Extracting [S II] 6716.')
fsii1 <- fluxsii1[which(fluxsii1!=500 & fluxsii2!=500)] #Obtaining the [O III] 4363 surface specific flux.

#[S II] 6731
print('Extracting [S II] 6731.')
fsii2 <- fluxsii2[which(fluxsii1!=500 & fluxsii2!=500)] #Obtaining the [O III] 4363 flux density (10^-16 erg/(s cm^2 A pix))

#H-beta dereddened flux.
print('Extracting E(B-V).')
EBV <- EBV[which(fluxsii1!=500 & fluxsii2!=500)] #Obtaining 


##Convert from surface specific intensity to standard specific intensity [erg cm^-2 s^-1 A^-1 sr^-1]
#print('Converting surface specific units to specific intensities.')
#foiii4363   <-(1e-16) * foiii4363 / 2.3504e-11 
#fb        <-(1e-16) * fb / 2.3504e-11 

##Multiply all values by 1e-16 (CALIFA reference level) Comment this section if using the previous conversion to specific intensity.

fsii1 <- 1e-16*fsii1
fsii2 <- 1e-16*fsii2


#####################################################

##                    DEREDDENING                 ##  

#####################################################

print('Dereddening...')

Rv <- 4.05 #Value for star-forming or high redshift galaxies, in Calzetti et al=. (2000)
l <- c(6716,6731) #Lines we are going to deredden, in Angstrom.
l <- l*1e-4 #Convert Angstrom to micrometers

#Calzetti et al. (2000) extinction law.
#k_l1 <- (2.659*(-2.156+(1.509/l)-(0.198/l^2)+(0.011/l^3)))+Rv #For  0.12 to 0.63 micrometres
k_l2 <- (2.659*(-1.857+(1.040/l)))+Rv #For 0.63 to 2.20 micrometres.

Fsii1  <-  fsii1*10^(0.4*EBV*k_l2[1]) #[S II] 6716 dereddening
Fsii2  <-  fsii2*10^(0.4*EBV*k_l2[2]) #[S II] 6731 dereddening

##############################################################################
##Save data to files

print('Saving fluxes to data file.')
resume <- data.frame(id, Fsii1, Fsii2)
tabla <- str_c(galaxy[i],"/Calzetti_sii_Fluxes.dat")
write.table(resume, tabla, sep="\t",quote=FALSE)

}