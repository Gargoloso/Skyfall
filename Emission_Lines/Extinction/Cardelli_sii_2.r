#############################################################################
#############################################################################

##   This calculates extinction corrected fluxes for the [S II] 6716,6731      ##
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

setwd("~/Rings/ringed_work/") #Directrory with our data.

data <- read.table("lis.dat", header=TRUE)
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
  
print('Extracting [S II] 6716,6731  data.')  
path_oi <- str_c(galaxy[i],"/lis_sii.res") #creates path to the [S II] 6716,6731 data file for each galaxy.
data0 <- read.table(path_oi, header=TRUE)


print('Extracting complementary data.') #Creates path to a file with complementary data, from a previous run of Cardelli_Base_Fluxes.dat
pathcb <- str_c(galaxy[i],"/Cardelli_Base_Fluxes.dat")
data1 <- read.table(pathcb, header=TRUE) 


##Merge data
DATA <- merge(data0, data1, by.x = 1, by.y = 1)
attach(DATA)
  
  
##Extracting coordinates from ID

print('Extracting spaxel coordinates.')

id2 <-  id[which(!is.na(fluxsii1) & snsii1 >= 5 & !is.na(eqwsii1) & !is.na(fluxsii2) & snsii2 >= 5 & !is.na(eqwsii2) & !is.na(EBV) & !is.na(rms_EBV))]

# E(B-V) colour excess.
print('Extracting E(B-V) colour excess.')
EBV2 <- EBV[which(!is.na(fluxsii1) & snsii1 >= 5 & !is.na(eqwsii1) & !is.na(fluxsii2) & snsii2 >= 5 & !is.na(eqwsii2) & !is.na(EBV) & !is.na(rms_EBV))]
rms_EBV2 <- rms_EBV[which(!is.na(fluxsii1) & snsii1 >= 5 & !is.na(eqwsii1) & !is.na(fluxsii2) & snsii2 >= 5 & !is.na(eqwsii2) & !is.na(EBV) & !is.na(rms_EBV))]

##Extracting line surface specific flux [1e-16 erg cm^-2 s^-1 A^-1 arcsec^-2].
print('Extracting emission lines.')

#[S II] 6716,6731.
print('Extracting [S II] 6716,6731.')
#[S II] 6716
fsii1    <-  1e-16*fluxsii1[which(!is.na(fluxsii1) & snsii1 >= 5 & !is.na(eqwsii1) & !is.na(fluxsii2) & snsii2 >= 5 & !is.na(eqwsii2) & !is.na(EBV) & !is.na(rms_EBV))] # Obtaining the [O III] 4959 surface specific flux.
snr_sii1 <-    snsii1[which(!is.na(fluxsii1) & snsii1 >= 5 & !is.na(eqwsii1) & !is.na(fluxsii2) & snsii2 >= 5 & !is.na(eqwsii2) & !is.na(EBV) & !is.na(rms_EBV))] # Signal to noise ratio of the continuum adjacent to the line.
sig_sii1 <- sigmasii1[which(!is.na(fluxsii1) & snsii1 >= 5 & !is.na(eqwsii1) & !is.na(fluxsii2) & snsii2 >= 5 & !is.na(eqwsii2) & !is.na(EBV) & !is.na(rms_EBV))] # Sigma of the Gaussian curve (width).
ew_sii1  <-   eqwsii1[which(!is.na(fluxsii1) & snsii1 >= 5 & !is.na(eqwsii1) & !is.na(fluxsii2) & snsii2 >= 5 & !is.na(eqwsii2) & !is.na(EBV) & !is.na(rms_EBV))] # Equivalent width of the emission line.
Npix_sii1 <- 2*fwhmsii1[which(!is.na(fluxsii1) & snsii1 >= 5 & !is.na(eqwsii1) & !is.na(fluxsii2) & snsii2 >= 5 & !is.na(eqwsii2) & !is.na(EBV) & !is.na(rms_EBV))] # Number of pixels covered by the line (approx. = 2*FWHM)

#[S II] 6731
fsii2    <-  1e-16*fluxsii2[which(!is.na(fluxsii1) & snsii1 >= 5 & !is.na(eqwsii1) & !is.na(fluxsii2) & snsii2 >= 5 & !is.na(eqwsii2) & !is.na(EBV) & !is.na(rms_EBV))] # Obtaining the [O III] 4959 surface specific flux.
snr_sii2 <-    snsii2[which(!is.na(fluxsii1) & snsii1 >= 5 & !is.na(eqwsii1) & !is.na(fluxsii2) & snsii2 >= 5 & !is.na(eqwsii2) & !is.na(EBV) & !is.na(rms_EBV))] # Signal to noise ratio of the continuum adjacent to the line.
sig_sii2 <- sigmasii2[which(!is.na(fluxsii1) & snsii1 >= 5 & !is.na(eqwsii1) & !is.na(fluxsii2) & snsii2 >= 5 & !is.na(eqwsii2) & !is.na(EBV) & !is.na(rms_EBV))] # Sigma of the Gaussian curve (width).
ew_sii2  <-   eqwsii2[which(!is.na(fluxsii1) & snsii1 >= 5 & !is.na(eqwsii1) & !is.na(fluxsii2) & snsii2 >= 5 & !is.na(eqwsii2) & !is.na(EBV) & !is.na(rms_EBV))] # Equivalent width of the emission line.
Npix_sii2 <- 2*fwhmsii2[which(!is.na(fluxsii1) & snsii1 >= 5 & !is.na(eqwsii1) & !is.na(fluxsii2) & snsii2 >= 5 & !is.na(eqwsii2) & !is.na(EBV) & !is.na(rms_EBV))] # Number of pixels covered by the line (approx. = 2*FWHM)

##Convert from surface specific intensity to standard specific intensity [erg cm^-2 s^-1 A^-1 sr^-1]
#print('Converting surface specific units to specific intensities.')
#fsii1   <- fsii1 / 2.3504e-11 


###################################################

##         FLUX ERROR DETERMINATION              ##

###################################################

# Continuum rms = 1/(Height of the Gaussian*snr of the continuum).
rmsc_sii1  <- 1e-16/(sig_sii1*sqrt(2*pi)*snr_sii1)
rmsc_sii2  <- 1e-16/(sig_sii2*sqrt(2*pi)*snr_sii2)

# Spectral dispersion (Angstrom/pixel).
D <- 1

# Flux error in the reddened emission line (As in Tresse et al. 1999).
rms_sii1  <-  rmsc_sii1*D*sqrt((2*Npix_sii1)+(ew_sii1/D))
rms_sii2  <-  rmsc_sii2*D*sqrt((2*Npix_sii2)+(ew_sii2/D))

# Filtering out NAs.
ID <- id2[which(!is.na(rms_sii1) & !is.na(rms_sii2) & !is.na(EBV2) & !is.na(rms_EBV2))]

fsii12    <- fsii1[which(!is.na(rms_sii1) & !is.na(rms_sii2) & !is.na(EBV2) & !is.na(rms_EBV2))]
rms_sii12 <- rms_sii1[which(!is.na(rms_sii1) & !is.na(rms_sii2) & !is.na(EBV2) & !is.na(rms_EBV2))]

fsii22    <- fsii2[which(!is.na(rms_sii1) & !is.na(rms_sii2) & !is.na(EBV2) & !is.na(rms_EBV2))]
rms_sii22 <- rms_sii2[which(!is.na(rms_sii1) & !is.na(rms_sii2) & !is.na(EBV2) & !is.na(rms_EBV2))]

EBV3    <- EBV2[which(!is.na(rms_sii1) & !is.na(rms_sii2) & !is.na(EBV2) & !is.na(rms_EBV2))]
rms_EBV3 <- rms_EBV2[which(!is.na(rms_sii1) & !is.na(rms_sii2) & !is.na(EBV2) & !is.na(rms_EBV2))]

#####################################################

##                    DEREDDENING                 ##  

#####################################################

print('Dereddening.')

Rv <- 3.1
l <-  c(6716,6731) #Lines we are going to unredden, in Angstrom.
l <- l*1e-4 #Convert Angstrom to micrometers.
xl <- 1/l
yl <- xl - 1.82

a_x <- 1 + (0.17699*yl) - (0.50447*yl^2) - (0.02427*yl^3) + (0.72085*yl^4) + (0.01979*yl^5) - (0.77530*yl^6) + (0.32999*yl^7)
b_x <- (1.41338*yl) + (2.28305*yl^2) + (1.07233*yl^3) -(5.38434*yl^4) - (0.62251*yl^5) + (5.30260*yl^6) - (2.09002*yl^7)

f_l <- a_x + (b_x/ Rv)

Fsii1     <-  fsii12*10^(0.4*EBV3*f_l[1]) #H-alpha dereddening
rms_Fsii1 <-  sqrt((rms_sii12*10^(0.4*EBV3*f_l[1]))^2 + (rms_EBV3*log(10)*fsii12*exp(log(10)*EBV3*f_l[1]))^2)

Fsii2     <-  fsii22*10^(0.4*EBV3*f_l[2]) #H-alpha dereddening
rms_Fsii2 <-  sqrt((rms_sii22*10^(0.4*EBV3*f_l[2]))^2 + (rms_EBV3*log(10)*fsii22*exp(log(10)*EBV3*f_l[2]))^2)

##############################################################################
##Save data to files.

print('Saving fluxes to data file.')
resume <- data.frame(ID, Fsii1, rms_sii1, Fsii2, rms_Fsii2)
tabla <- str_c(galaxy[i],"/Cardelli_sii_Fluxes.dat")
write.table(resume, tabla, sep="\t",quote=FALSE)

}