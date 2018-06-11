#############################################################################
#############################################################################

##   This calculates extinction corrected fluxes for the [O III] 4363      ##
##   emission line using the Calzetti et al. (2000) extinction law.        ##    

##  June 10, 2018 A. Robleto-Orús.                                         ## 

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

##################################

##				DATA INPUT			      ##

##################################

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
  
print('Extracting [O III] 4363  data.')  
path_oi <- str_c(galaxy[i],"/lis_hgoiii4363.res") #creates path [O III] 4363 data file for each galaxy.
data0 <- read.table(path_oi, header=TRUE)
  
print('Extracting Base data') #Creates path to a file with Base results.
pathcb <- str_c(galaxy[i],"/Calzetti_Base_Fluxes.dat")
data1 <- read.table(pathcb, header=TRUE) 

##Merge data
DATA <- merge(data0, data1, by.x = 1, by.y = 1)
attach(DATA)
  
  
##Extracting coordinates from ID

print('Extracting spaxels IDs.')
id <- id[which(fluxoiii4363!=500)]

##Extracting line surface specific intensities [1e-16 erg cm^-2 s^-1 A^-1 arcsec^-2].

print('Extracting emission lines.')

#[O III] 4363
print('Extracting [O III] 4363.')
foiii4363 <- fluxoiii4363[which(fluxoiii4363!=500)] #Obtaining the [O III] 4363 surface specific flux.

#H-beta dereddened flux.
print('Extracting E(B-V).')
EBV <- EBV[which(fluxoiii4363!=500)] #Obtaining E(B-V) colour excess from Base.


##Convert from surface specific intensity to standard specific intensity [erg cm^-2 s^-1 A^-1 sr^-1]
#print('Converting surface specific units to specific intensities.')
#foiii4363   <-(1e-16) * foiii4363 / 2.3504e-11 
#fb        <-(1e-16) * fb / 2.3504e-11 

##Multiply all values by 1e-16 (CALIFA reference level) Comment this section if using the previous conversion to specific intensity.

foiii4363  <-(1e-16)*foiii4363


#####################################################

##                    DEREDDENING                 ##  

#####################################################

print('Dereddening...')

Rv <- 4.05  #Value for star-forming or high redshift galaxies, in Calzetti et al=. (2000)
l <-  4363  #Lines we are going to unredden, in Angstrom.
l <- l*1e-4 #Convert Angstrom to micrometers

k_l1 <- (2.659*(-2.156+(1.509/l)-(0.198/l^2)+(0.011/l^3)))+Rv #For  0.12 to 0.63 micrometres
#k_l2 <- (2.659*(-1.857+(1.040/l)))+Rv #For 0.63 to 2.20 micrometres.

Foiii4363  <-  foiii4363*10^(0.4*EBV*k_l1) #[O III] 6583 dereddening

##############################################################################
##Save data to files

print('Saving fluxes to data file.')
resume <- data.frame(id, Foiii4363)
tabla <- str_c(galaxy[i],"/Calzetti_oiii4363_Fluxes.dat")
write.table(resume, tabla, sep="\t",quote=FALSE)

}