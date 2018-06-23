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

print('Extracting spaxel coordinates.')

id2 <-  id[which(!is.na(fluxoiii4363) & snoiii4363 >= 5 & !is.na(eqwoiii4363) & !is.na(EBV) & !is.na(rms_EBV))]

# E(B-V) colour excess.
print('Extracting E(B-V) colour excess.')
EBV2 <- EBV[which(!is.na(fluxoiii4363) & snoiii4363 >= 5 & !is.na(eqwoiii4363) & !is.na(EBV) & !is.na(rms_EBV))]
rms_EBV2 <- rms_EBV[which(!is.na(fluxoiii4363) & snoiii4363 >= 5 & !is.na(eqwoiii4363) & !is.na(EBV) & !is.na(rms_EBV))]

##Extracting line surface specific flux [1e-16 erg cm^-2 s^-1 A^-1 arcsec^-2].
print('Extracting emission lines.')

#[O III] 4363.
print('Extracting [O III] 4363.')
foiii4363    <-  1e-16*fluxoiii4363[which(!is.na(fluxoiii4363) & snoiii4363 >= 5 & !is.na(eqwoiii4363) & !is.na(EBV) & !is.na(rms_EBV))] # Obtaining the [O III] 4959 surface specific flux.
snr_oiii4363 <-    snoiii4363[which(!is.na(fluxoiii4363) & snoiii4363 >= 5 & !is.na(eqwoiii4363) & !is.na(EBV) & !is.na(rms_EBV))] # Signal to noise ratio of the continuum adjacent to the line.
sig_oiii4363 <- sigmaoiii4363[which(!is.na(fluxoiii4363) & snoiii4363 >= 5 & !is.na(eqwoiii4363) & !is.na(EBV) & !is.na(rms_EBV))] # Sigma of the Gaussian curve (width).
ew_oiii4363  <-   eqwoiii4363[which(!is.na(fluxoiii4363) & snoiii4363 >= 5 & !is.na(eqwoiii4363) & !is.na(EBV) & !is.na(rms_EBV))] # Equivalent width of the emission line.
Npix_oiii4363 <- 2*fwhmoiii4363[which(!is.na(fluxoiii4363) & snoiii4363 >= 5 & !is.na(eqwoiii4363) & !is.na(EBV) & !is.na(rms_EBV))] # Number of pixels covered by the line (approx. = 2*FWHM)


##Convert from surface specific intensity to standard specific intensity [erg cm^-2 s^-1 A^-1 sr^-1]
#print('Converting surface specific units to specific intensities.')
#foiii4363   <- foiii4363 / 2.3504e-11 


###################################################

##         FLUX ERROR DETERMINATION              ##

###################################################

# Continuum rms = 1/(Height of the Gaussian*snr of the continuum).
rmsc_oiii4363  <- 1e-16/(sig_oiii4363*sqrt(2*pi)*snr_oiii4363)

# Spectral dispersion (Angstrom/pixel).
D <- 1

# Flux error in the reddened emission line (As in Tresse et al. 1999).
rms_oiii4363  <-  rmsc_oiii4363*D*sqrt((2*Npix_oiii4363)+(ew_oiii4363/D))

# Filtering out NAs.
ID <- id2[which(!is.na(rms_oiii4363) & !is.na(EBV2) & !is.na(rms_EBV2))]

foiii43632    <- foiii4363[which(!is.na(rms_oiii4363) & !is.na(EBV2) & !is.na(rms_EBV2))]
rms_oiii43632 <- rms_oiii4363[which(!is.na(rms_oiii4363) & !is.na(EBV2) & !is.na(rms_EBV2))]

EBV3    <- EBV2[which(!is.na(rms_oiii4363) & !is.na(EBV2) & !is.na(rms_EBV2))]
rms_EBV3 <- rms_EBV2[which(!is.na(rms_oiii4363) & !is.na(EBV2) & !is.na(rms_EBV2))]


#####################################################

##                    DEREDDENING                 ##  

#####################################################

print('Dereddening...')

Rv <- 4.05  #Value for star-forming or high redshift galaxies, in Calzetti et al=. (2000)
l <-  4363  #Lines we are going to unredden, in Angstrom.
l <- l*1e-4 #Convert Angstrom to micrometers

k_l1 <- (2.659*(-2.156+(1.509/l)-(0.198/l^2)+(0.011/l^3)))+Rv #For  0.12 to 0.63 micrometres
#k_l2 <- (2.659*(-1.857+(1.040/l)))+Rv #For 0.63 to 2.20 micrometres.

Foiii4363  <-  foiii43632*10^(0.4*EBV3*k_l1) #[O III] 6583 dereddening
rms_Foiii4363 <-  sqrt((rms_oiii43632*10^(0.4*EBV3*k_l1))^2 + (rms_EBV3*log(10)*foiii43632*exp(log(10)*EBV3*k_l1))^2)

##############################################################################
##Save data to files

print('Saving fluxes to data file.')
resume <- data.frame(ID, Foiii4363, rms_Foiii4363)
tabla <- str_c(galaxy[i],"/Calzetti_oiii4363_Fluxes.dat")
write.table(resume, tabla, sep="\t",quote=FALSE)

}