#############################################################################
#############################################################################

##   This calculates extinction corrected fluxes using the Cardelli        ##
##   extinction law.                                                       ##    

##   February 25, 2018 A. Robleto-Orús                                     ##
##   May 15, 2018 Adapted to work with the last version of fitlines.       ##
##                We now only save the final surface specific fluxes but   ##
##                not line ratios. Also Cbeta and AHa are saved.           ##
##   June 10, 2018 Optimized.                                              ##
##   June 15, 2018 Added error estimation as in Tresse et al. (1999)       ##
##   June 15, 2018 Not using H-beta normalization anymore. Direct flux     ##
##                 dereddening as in Calzetti et al. (2001) but using the  ##
##                 Cardelli et al. (1989) extinction law.                  ##

#############################################################################
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

print('Extracting H-alpha and [N II] 6583  data.')  
path_ha <- str_c(galaxy[i],"/lis_ha.res") #creates path to Ha and [N II] data file for each galaxy.
data0 <- read.table(path_ha, header=TRUE)

print('Extracting H-beta and [O III] 5007 data.')
path_hb <- str_c(galaxy[i],"/lis_hb.res") #creates path to H-beta and [OIII] data file for each galaxy.
data3 <- read.table(path_hb, header=TRUE)
         
##Merge data.

DATA <- merge(data0,data3,by.x = 1, by.y =1)
attach(DATA)


##Extracting spaxels ID

print('Extracting spaxels ID.')
id2 <-  id[which(fluxha!=500 & !is.na(fluxha) & !is.na(fluxnii1) & !is.na(fluxnii2) & !is.na(fluxoiii1) & !is.na(fluxoiii2)  & snhb >= 5 & snha >= 5  & !is.na(eqwha) & !is.na(eqwnii1) & !is.na(eqwnii2) & !is.na(eqwoiii1) & !is.na(eqwoiii2))]

##Extracting line surface specific fluxes [1e-16 erg cm^-2 s^-1 A^-1 arcsec^-2].

print('Extracting emission lines.')

#H-alpha 6563
print('Extracting H-alpha.')
fa    <-  1e-16*fluxha[which(fluxha!=500 & !is.na(fluxha) & !is.na(fluxnii1) & !is.na(fluxnii2) & !is.na(fluxoiii1) & !is.na(fluxoiii2) & snhb >= 5 & snha >= 5  & !is.na(eqwha) & !is.na(eqwnii1) & !is.na(eqwnii2) & !is.na(eqwoiii1) & !is.na(eqwoiii2))] # Obtaining the H-alpha surface specific flux.
snr_a <-    snha[which(fluxha!=500 & !is.na(fluxha) & !is.na(fluxnii1) & !is.na(fluxnii2) & !is.na(fluxoiii1) & !is.na(fluxoiii2) & snhb >= 5 & snha >= 5  & !is.na(eqwha) & !is.na(eqwnii1) & !is.na(eqwnii2) & !is.na(eqwoiii1) & !is.na(eqwoiii2))] # Signal to noise ratio of the continuum adjacent to the line.
sig_a <- sigmaha[which(fluxha!=500 & !is.na(fluxha) & !is.na(fluxnii1) & !is.na(fluxnii2) & !is.na(fluxoiii1) & !is.na(fluxoiii2) & snhb >= 5 & snha >= 5  & !is.na(eqwha) & !is.na(eqwnii1) & !is.na(eqwnii2) & !is.na(eqwoiii1) & !is.na(eqwoiii2))] # Sigma of the Gaussian curve (width).
ew_a  <-   eqwha[which(fluxha!=500 & !is.na(fluxha) & !is.na(fluxnii1) & !is.na(fluxnii2) & !is.na(fluxoiii1) & !is.na(fluxoiii2) & snhb >= 5 & snha >= 5  & !is.na(eqwha) & !is.na(eqwnii1) & !is.na(eqwnii2) & !is.na(eqwoiii1) & !is.na(eqwoiii2))] # Equivalent width of the emission line.
Npix_a  <- 2*fwhmha[which(fluxha!=500 & !is.na(fluxha) & !is.na(fluxnii1) & !is.na(fluxnii2) & !is.na(fluxoiii1) & !is.na(fluxoiii2) & snhb >= 5 & snha >= 5  & !is.na(eqwha) & !is.na(eqwnii1) & !is.na(eqwnii2) & !is.na(eqwoiii1) & !is.na(eqwoiii2))] # Number of pixels covered by the line (approx. = 2*FWHM)

#[N II] 6548
print('Extracting [N II] 6548')
fn1    <-  1e-16*fluxnii1[which(fluxha!=500 & !is.na(fluxha) & !is.na(fluxnii1) & !is.na(fluxnii2) & !is.na(fluxoiii1) & !is.na(fluxoiii2) & snhb >= 5 & snha >= 5  & !is.na(eqwha) & !is.na(eqwnii1) & !is.na(eqwnii2) & !is.na(eqwoiii1) & !is.na(eqwoiii2))] # Obtaining the  [N II] 6548 surface specific flux.
snr_n1 <-    snnii1[which(fluxha!=500 & !is.na(fluxha) & !is.na(fluxnii1) & !is.na(fluxnii2) & !is.na(fluxoiii1) & !is.na(fluxoiii2) & snhb >= 5 & snha >= 5  & !is.na(eqwha) & !is.na(eqwnii1) & !is.na(eqwnii2) & !is.na(eqwoiii1) & !is.na(eqwoiii2))] # Signal to noise ratio of the continuum adjacent to the line.
sig_n1 <- sigmanii1[which(fluxha!=500 & !is.na(fluxha) & !is.na(fluxnii1) & !is.na(fluxnii2) & !is.na(fluxoiii1) & !is.na(fluxoiii2) & snhb >= 5 & snha >= 5  & !is.na(eqwha) & !is.na(eqwnii1) & !is.na(eqwnii2) & !is.na(eqwoiii1) & !is.na(eqwoiii2))] # Sigma of the Gaussian curve (width).
ew_n1  <-   eqwnii1[which(fluxha!=500 & !is.na(fluxha) & !is.na(fluxnii1) & !is.na(fluxnii2) & !is.na(fluxoiii1) & !is.na(fluxoiii2) & snhb >= 5 & snha >= 5  & !is.na(eqwha) & !is.na(eqwnii1) & !is.na(eqwnii2) & !is.na(eqwoiii1) & !is.na(eqwoiii2))] # Equivalent width of the emission line.
Npix_n1 <- 2*fwhmnii1[which(fluxha!=500 & !is.na(fluxha) & !is.na(fluxnii1) & !is.na(fluxnii2) & !is.na(fluxoiii1) & !is.na(fluxoiii2) & snhb >= 5 & snha >= 5  & !is.na(eqwha) & !is.na(eqwnii1) & !is.na(eqwnii2) & !is.na(eqwoiii1) & !is.na(eqwoiii2))] # Number of pixels covered by the line (approx. = 2*FWHM)

#[N II] 6583
print('Extracting [N II] 6583.')
fn    <-  1e-16*fluxnii2[which(fluxha!=500 & !is.na(fluxha) & !is.na(fluxnii1) & !is.na(fluxnii2) & !is.na(fluxoiii1) & !is.na(fluxoiii2) & snhb >= 5 & snha >= 5  & !is.na(eqwha) & !is.na(eqwnii1) & !is.na(eqwnii2) & !is.na(eqwoiii1) & !is.na(eqwoiii2))] # Obtaining the [N II] 6583 surface specific flux.
snr_n <-    snnii2[which(fluxha!=500 & !is.na(fluxha) & !is.na(fluxnii1) & !is.na(fluxnii2) & !is.na(fluxoiii1) & !is.na(fluxoiii2) & snhb >= 5 & snha >= 5  & !is.na(eqwha) & !is.na(eqwnii1) & !is.na(eqwnii2) & !is.na(eqwoiii1) & !is.na(eqwoiii2))] # Signal to noise ratio of the continuum adjacent to the line.
sig_n <- sigmanii2[which(fluxha!=500 & !is.na(fluxha) & !is.na(fluxnii1) & !is.na(fluxnii2) & !is.na(fluxoiii1) & !is.na(fluxoiii2) & snhb >= 5 & snha >= 5  & !is.na(eqwha) & !is.na(eqwnii1) & !is.na(eqwnii2) & !is.na(eqwoiii1) & !is.na(eqwoiii2))] # Sigma of the Gaussian curve (width).
ew_n  <-   eqwnii2[which(fluxha!=500 & !is.na(fluxha) & !is.na(fluxnii1) & !is.na(fluxnii2) & !is.na(fluxoiii1) & !is.na(fluxoiii2) & snhb >= 5 & snha >= 5  & !is.na(eqwha) & !is.na(eqwnii1) & !is.na(eqwnii2) & !is.na(eqwoiii1) & !is.na(eqwoiii2))] # Equivalent width of the emission line.
Npix_n  <- 2*fwhmnii2[which(fluxha!=500 & !is.na(fluxha) & !is.na(fluxnii1) & !is.na(fluxnii2) & !is.na(fluxoiii1) & !is.na(fluxoiii2) & snhb >= 5 & snha >= 5  & !is.na(eqwha) & !is.na(eqwnii1) & !is.na(eqwnii2) & !is.na(eqwoiii1) & !is.na(eqwoiii2))] # Number of pixels covered by the line (approx. = 2*FWHM)

#H-beta 4861
print('Extracting H-beta.')
fb    <-  1e-16*fluxhb[which(fluxha!=500 & !is.na(fluxha) & !is.na(fluxnii1) & !is.na(fluxnii2) & !is.na(fluxoiii1) & !is.na(fluxoiii2) & snhb >= 5 & snha >= 5  & !is.na(eqwha) & !is.na(eqwnii1) & !is.na(eqwnii2) & !is.na(eqwoiii1) & !is.na(eqwoiii2))] # Obtaining the [H-beta] surface specific flux.
snr_b <-    snhb[which(fluxha!=500 & !is.na(fluxha) & !is.na(fluxnii1) & !is.na(fluxnii2) & !is.na(fluxoiii1) & !is.na(fluxoiii2) & snhb >= 5 & snha >= 5  & !is.na(eqwha) & !is.na(eqwnii1) & !is.na(eqwnii2) & !is.na(eqwoiii1) & !is.na(eqwoiii2))] # Signal to noise ratio of the continuum adjacent to the line.
sig_b <- sigmahb[which(fluxha!=500 & !is.na(fluxha) & !is.na(fluxnii1) & !is.na(fluxnii2) & !is.na(fluxoiii1) & !is.na(fluxoiii2) & snhb >= 5 & snha >= 5  & !is.na(eqwha) & !is.na(eqwnii1) & !is.na(eqwnii2) & !is.na(eqwoiii1) & !is.na(eqwoiii2))] # Sigma of the Gaussian curve (width).
ew_b  <-   eqwhb[which(fluxha!=500 & !is.na(fluxha) & !is.na(fluxnii1) & !is.na(fluxnii2) & !is.na(fluxoiii1) & !is.na(fluxoiii2) & snhb >= 5 & snha >= 5  & !is.na(eqwha) & !is.na(eqwnii1) & !is.na(eqwnii2) & !is.na(eqwoiii1) & !is.na(eqwoiii2))] # Equivalent width of the emission line.
Npix_b  <- 2*fwhmhb[which(fluxha!=500 & !is.na(fluxha) & !is.na(fluxnii1) & !is.na(fluxnii2) & !is.na(fluxoiii1) & !is.na(fluxoiii2) & snhb >= 5 & snha >= 5  & !is.na(eqwha) & !is.na(eqwnii1) & !is.na(eqwnii2) & !is.na(eqwoiii1) & !is.na(eqwoiii2))] # Number of pixels covered by the line (approx. = 2*FWHM)

#[O III] 4959
print('Extracting [O III] 4959')
fo1    <-  1e-16*fluxoiii1[which(fluxha!=500 & !is.na(fluxha) & !is.na(fluxnii1) & !is.na(fluxnii2) & !is.na(fluxoiii1) & !is.na(fluxoiii2) & snhb >= 5 & snha >= 5  & !is.na(eqwha) & !is.na(eqwnii1) & !is.na(eqwnii2) & !is.na(eqwoiii1) & !is.na(eqwoiii2))] # Obtaining the [O III] 4959 surface specific flux.
snr_o1 <-    snoiii1[which(fluxha!=500 & !is.na(fluxha) & !is.na(fluxnii1) & !is.na(fluxnii2) & !is.na(fluxoiii1) & !is.na(fluxoiii2) & snhb >= 5 & snha >= 5  & !is.na(eqwha) & !is.na(eqwnii1) & !is.na(eqwnii2) & !is.na(eqwoiii1) & !is.na(eqwoiii2))] # Signal to noise ratio of the continuum adjacent to the line.
sig_o1 <- sigmaoiii1[which(fluxha!=500 & !is.na(fluxha) & !is.na(fluxnii1) & !is.na(fluxnii2) & !is.na(fluxoiii1) & !is.na(fluxoiii2) & snhb >= 5 & snha >= 5  & !is.na(eqwha) & !is.na(eqwnii1) & !is.na(eqwnii2) & !is.na(eqwoiii1) & !is.na(eqwoiii2))] # Sigma of the Gaussian curve (width).
ew_o1  <-   eqwoiii1[which(fluxha!=500 & !is.na(fluxha) & !is.na(fluxnii1) & !is.na(fluxnii2) & !is.na(fluxoiii1) & !is.na(fluxoiii2) & snhb >= 5 & snha >= 5  & !is.na(eqwha) & !is.na(eqwnii1) & !is.na(eqwnii2) & !is.na(eqwoiii1) & !is.na(eqwoiii2))] # Equivalent width of the emission line.
Npix_o1 <- 2*fwhmoiii1[which(fluxha!=500 & !is.na(fluxha) & !is.na(fluxnii1) & !is.na(fluxnii2) & !is.na(fluxoiii1) & !is.na(fluxoiii2) & snhb >= 5 & snha >= 5  & !is.na(eqwha) & !is.na(eqwnii1) & !is.na(eqwnii2) & !is.na(eqwoiii1) & !is.na(eqwoiii2))] # Number of pixels covered by the line (approx. = 2*FWHM)

#[O III] 5007
print('Extracting [O III] 5007.')
fo    <-  1e-16*fluxoiii2[which(fluxha!=500 & !is.na(fluxha) & !is.na(fluxnii1) & !is.na(fluxnii2) & !is.na(fluxoiii1) & !is.na(fluxoiii2) & snhb >= 5 & snha >= 5  & !is.na(eqwha) & !is.na(eqwnii1) & !is.na(eqwnii2) & !is.na(eqwoiii1) & !is.na(eqwoiii2))] # Obtaining the [O III] 5007 surface specific flux.
snr_o <-    snoiii2[which(fluxha!=500 & !is.na(fluxha) & !is.na(fluxnii1) & !is.na(fluxnii2) & !is.na(fluxoiii1) & !is.na(fluxoiii2) & snhb >= 5 & snha >= 5  & !is.na(eqwha) & !is.na(eqwnii1) & !is.na(eqwnii2) & !is.na(eqwoiii1) & !is.na(eqwoiii2))] # Signal to noise ratio of the continuum adjacent to the line.
sig_o <- sigmaoiii2[which(fluxha!=500 & !is.na(fluxha) & !is.na(fluxnii1) & !is.na(fluxnii2) & !is.na(fluxoiii1) & !is.na(fluxoiii2) & snhb >= 5 & snha >= 5  & !is.na(eqwha) & !is.na(eqwnii1) & !is.na(eqwnii2) & !is.na(eqwoiii1) & !is.na(eqwoiii2))] # Sigma of the Gaussian curve (width).
ew_o  <-   eqwoiii2[which(fluxha!=500 & !is.na(fluxha) & !is.na(fluxnii1) & !is.na(fluxnii2) & !is.na(fluxoiii1) & !is.na(fluxoiii2) & snhb >= 5 & snha >= 5  & !is.na(eqwha) & !is.na(eqwnii1) & !is.na(eqwnii2) & !is.na(eqwoiii1) & !is.na(eqwoiii2))] # Equivalent width of the emission line.
Npix_o  <- 2*fwhmoiii2[which(fluxha!=500 & !is.na(fluxha) & !is.na(fluxnii1) & !is.na(fluxnii2) & !is.na(fluxoiii1) & !is.na(fluxoiii2) & snhb >= 5 & snha >= 5  & !is.na(eqwha) & !is.na(eqwnii1) & !is.na(eqwnii2) & !is.na(eqwoiii1) & !is.na(eqwoiii2))] # Number of pixels covered by the line (approx. = 2*FWHM)

###################################################

##         FLUX ERROR DETERMINATION              ##

###################################################

# Continuum rms = 1/(Height of the Gaussian*snr of the continuum).
rmsc_a  <- 1e-16/(sig_a*sqrt(2*pi)*snr_a)
rmsc_n1 <- 1e-16/(sig_n1*sqrt(2*pi)*snr_n1)
rmsc_n  <- 1e-16/(sig_n*sqrt(2*pi)*snr_n)
rmsc_b  <- 1e-16/(sig_b*sqrt(2*pi)*snr_b)
rmsc_o1 <- 1e-16/(sig_o1*sqrt(2*pi)*snr_o1)
rmsc_o  <- 1e-16/(sig_o*sqrt(2*pi)*snr_o)

# Spectral dispersion (Angstrom/pixel)
D <- 1

# Flux error in the reddened emission line (As in Tresse et al. 1999).
rms_a  <-  rmsc_a*D*sqrt((2*Npix_a)+(ew_a/D))
rms_n1 <- rmsc_n1*D*sqrt((2*Npix_n1)+(ew_n1/D))
rms_n  <-  rmsc_n*D*sqrt((2*Npix_n)+(ew_n/D))
rms_b  <-  rmsc_b*D*sqrt((2*Npix_b)+(ew_b/D))
rms_o1 <- rmsc_o1*D*sqrt((2*Npix_o1)+(ew_o1/D))
rms_o  <-  rmsc_o*D*sqrt((2*Npix_o)+(ew_o/D))


##Convert from surface specific intensity to standard specific intensity [erg cm^-2 s^-1 A^-1 sr^-1].

#print('Converting surface specific flux  units to specific intensities.')
#fa   <-fa  /2.3504e-11 
#fn   <-fn  /2.3504e-11 
#fb   <-fb  /2.3504e-11 
#fo   <-fo  /2.3504e-11 

## Comment this section if using the previous conversion to specific intensity.

## Filtering out NAs.
fa2   <- fa[which(!is.na(rms_a) & !is.na(rms_n1) & !is.na(rms_n) & !is.na(rms_b) & !is.na(rms_o1) & !is.na(rms_o))]
fn12  <- fn1[which(!is.na(rms_a) & !is.na(rms_n1) & !is.na(rms_n) & !is.na(rms_b) & !is.na(rms_o1) & !is.na(rms_o))]
fn2   <- fn[which(!is.na(rms_a) & !is.na(rms_n1) & !is.na(rms_n) & !is.na(rms_b) & !is.na(rms_o1) & !is.na(rms_o))] 
fb2   <- fb[which(!is.na(rms_a) & !is.na(rms_n1) & !is.na(rms_n) & !is.na(rms_b) & !is.na(rms_o1) & !is.na(rms_o))] 
fo12  <- fo1[which(!is.na(rms_a) & !is.na(rms_n1) & !is.na(rms_n) & !is.na(rms_b) & !is.na(rms_o1) & !is.na(rms_o))]
fo2   <- fo[which(!is.na(rms_a) & !is.na(rms_n1) & !is.na(rms_n) & !is.na(rms_b) & !is.na(rms_o1) & !is.na(rms_o))] 

rms_a2  <-  rms_a[which(!is.na(rms_a) & !is.na(rms_n1) & !is.na(rms_n) & !is.na(rms_b) & !is.na(rms_o1) & !is.na(rms_o))]
rms_n12 <- rms_n1[which(!is.na(rms_a) & !is.na(rms_n1) & !is.na(rms_n) & !is.na(rms_b) & !is.na(rms_o1) & !is.na(rms_o))]
rms_n2  <-  rms_n[which(!is.na(rms_a) & !is.na(rms_n1) & !is.na(rms_n) & !is.na(rms_b) & !is.na(rms_o1) & !is.na(rms_o))]
rms_b2  <-  rms_b[which(!is.na(rms_a) & !is.na(rms_n1) & !is.na(rms_n) & !is.na(rms_b) & !is.na(rms_o1) & !is.na(rms_o))]
rms_o12 <- rms_o1[which(!is.na(rms_a) & !is.na(rms_n1) & !is.na(rms_n) & !is.na(rms_b) & !is.na(rms_o1) & !is.na(rms_o))]
rms_o2  <-  rms_o[which(!is.na(rms_a) & !is.na(rms_n1) & !is.na(rms_n) & !is.na(rms_b) & !is.na(rms_o1) & !is.na(rms_o))]

ID <- id2[which(!is.na(rms_a) & !is.na(rms_n1) & !is.na(rms_n) & !is.na(rms_b) & !is.na(rms_o1) & !is.na(rms_o))]

###################################################

##                 Dereddening                   ##

###################################################

print('Dereddening.')

Rv <- 3.1 
l <-  c(6563,6548,6583,4861,4959,5007) #Lines we are going to deredden, in Angstroms.
l <- l*1e-4 #Convert Angstroms to micrometres.
xl <- 1/l
yl <- xl - 1.82

## Cardelli et al. (1989) extinction law:
a_x <- 1 + (0.17699*yl) - (0.50447*yl^2) - (0.02427*yl^3) + (0.72085*yl^4) + (0.01979*yl^5) - (0.77530*yl^6) + (0.32999*yl^7)
b_x <- (1.41338*yl) + (2.28305*yl^2) + (1.07233*yl^3) - (5.38434*yl^4) - (0.62251*yl^5) + (5.30260*yl^6) - (2.09002*yl^7)

f_l <- a_x + (b_x/ Rv)

Rab <- (fa2/fb2)/2.86 #Ratio of attenuated Halpha/Hbeta over non-attenuated one.
rms_Rab <-sqrt((rms_a2/(2.86*fb2))^2 + ((-fa2*rms_b2)/(2.86*fb2^2))^2)
  
EBV <- log10(Rab)/(0.4*(f_l[4]-f_l[1])) #Colour excess E(B-V) in magnitudes.
rms_EBV <- (rms_Rab/(0.4*(f_l[4]-f_l[1])*Rab*log(10)))
  
Fa     <-  fa2*10^(0.4*EBV*f_l[1]) #H-alpha dereddening
rms_Fa <-  sqrt((rms_a2*10^(0.4*EBV*f_l[1]))^2 + (rms_EBV*log(10)*fa2*exp(log(10)*EBV*f_l[1]))^2)

Fb  <-  fb2*10^(0.4*EBV*f_l[4]) #H-beta dereddening
rms_Fb <-  sqrt((rms_b2*10^(0.4*EBV*f_l[4]))^2 + (rms_EBV*log(10)*fb2*exp(log(10)*EBV*f_l[4]))^2)

Fn1 <- fn12*10^(0.4*EBV*f_l[2]) #[N II] 6548 dereddening 
rms_Fn1 <-  sqrt((rms_n12*10^(0.4*EBV*f_l[2]))^2 + (rms_EBV*log(10)*fn12*exp(log(10)*EBV*f_l[2]))^2)

Fn  <-  fn2*10^(0.4*EBV*f_l[3]) #[N II] 6583 dereddening
rms_Fn <-  sqrt((rms_n2*10^(0.4*EBV*f_l[3]))^2 + (rms_EBV*log(10)*fn2*exp(log(10)*EBV*f_l[3]))^2)

Fo1 <- fo12*10^(0.4*EBV*f_l[5]) #[O III] 4959 dereddening
rms_Fo1 <-  sqrt((rms_o12*10^(0.4*EBV*f_l[5]))^2 + (rms_EBV*log(10)*fo12*exp(log(10)*EBV*f_l[5]))^2)

Fo  <-  fo2*10^(0.4*EBV*f_l[6]) #[O III] 6583 dereddening
rms_Fo <-  sqrt((rms_o2*10^(0.4*EBV*f_l[6]))^2 + (rms_EBV*log(10)*fo2*exp(log(10)*EBV*f_l[6]))^2)


AHa <- f_l[1]*EBV #H-alpha extinction (in magnitudes)
rms_AHa <- f_l[1]*rms_EBV



##############################################################################
##Save data to files.

print('Saving fluxes to data file.')
resume <- data.frame(ID,EBV,rms_EBV,AHa,rms_AHa,Fa,rms_Fa,Fb,rms_Fb,Fn1,rms_Fn1,Fn,rms_Fn,Fo1,rms_Fo1,Fo,rms_Fo)
#resume <- data.frame(ID,AHa,Fa,Fb,Fn1,Fn,Fo1,Fo)
tabla <- str_c(galaxy[i],"/Cardelli_Base_Fluxes.dat")
write.table(resume, tabla, sep="\t",quote=FALSE)

}
