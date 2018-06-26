#############################################################################
#############################################################################

##   This calculates extinction corrected fluxes for the [O III] 4363        ##
##   emission line using the CCM 89 extinction law.                        ##    

##   May 15, 2018 A. Robleto-Orús                                          ##
##   June 10, 2018 Optimized.                                              ##
##   June 19, 2018 Added error propagation. Not using H-beta normalization ##
##                 anymore. Direct flux dereddening as in Calzetti et al.  ##
##                 (2001) but using the Cardelli et al. (1989) extinction  ##
##                 law.                                                    ##

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

setwd("---") # Directrory with our data.

data <- read.table("---", header=TRUE) # File with general data for all galaxies.
attach(data)

galaxy <- as.character(name)

####################################################################
####################################################################

##                   ITERATION FOR ALL GALAXIES                   ##

####################################################################
####################################################################

##Loop for all galaxies.
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
path_oi <- str_c(---,"/lis_hgoiii4363.res") #creates path to the [O III] 4363 data file for each galaxy.
data0 <- read.table(path_oi, header=TRUE)


print('Extracting complementary data.') #Creates path to a file with complementary data, from a previous run of Cardelli_Base_Fluxes.dat
pathcb <- str_c(---,"/Cardelli_Base_Fluxes.dat")
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

print('Dereddening.')

Rv <- 3.1
l <-  6300 #Lines we are going to unredden, in Angstrom.
l <- l*1e-4 #Convert Angstrom to micrometers.
xl <- 1/l
yl <- xl - 1.82

a_x <- 1 + (0.17699*yl) - (0.50447*yl^2) - (0.02427*yl^3) + (0.72085*yl^4) + (0.01979*yl^5) - (0.77530*yl^6) + (0.32999*yl^7)
b_x <- (1.41338*yl) + (2.28305*yl^2) + (1.07233*yl^3) -(5.38434*yl^4) - (0.62251*yl^5) + (5.30260*yl^6) - (2.09002*yl^7)

f_l <- a_x + (b_x/ Rv)

Foiii4363     <-  foiii43632*10^(0.4*EBV3*f_l[1]) #H-alpha dereddening
rms_Foiii4363 <-  sqrt((rms_oiii43632*10^(0.4*EBV3*f_l[1]))^2 + (rms_EBV3*log(10)*foiii43632*exp(log(10)*EBV3*f_l[1]))^2)


##############################################################################
##Save data to files.

print('Saving fluxes to data file.')
resume <- data.frame(ID, Foiii4363, rms_oiii4363)
tabla <- str_c(---,"/Cardelli_oiii4363_Fluxes.dat") # Output data file.
write.table(resume, tabla, sep="\t",quote=FALSE)

}
