#############################################################################
#############################################################################

##   This calculates extinction corrected fluxes using the Cardelli        ##
##   extinction law. (Cardelli et al. 1989)                                ##    

##   February 25, 2018 A. Robleto-Orús                                     ##
##   May 15, 2018 Adapted to work with the last version of fitlines.       ##
##                We now only save the final surface specific fluxes but   ##
##                not line ratios. Also Cbeta and AHa are saved.           ##
##   June 10, 2018 Optimized.                                              ##
##   June 15, 2018 Added error estimation as in Arellano-Cordova (2015)    ##
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
id2 <-  id[which(fluxha!=500 & !is.na(fluxha) & !is.na(fluxnii1) & !is.na(fluxnii2) & !is.na(fluxoiii1) & !is.na(fluxoiii2)  & snhb >= 3 & snha >= 3  & !is.na(eqwha) & !is.na(eqwnii1) & !is.na(eqwnii2) & !is.na(eqwoiii1) & !is.na(eqwoiii2))]

##Extracting line surface specific fluxes [1e-16 erg cm^-2 s^-1 A^-1 arcsec^-2].

print('Extracting emission lines.')

#H-alpha 6563
print('Extracting H-alpha.')
fa    <-  1e-16*fluxha[which(fluxha!=500 & !is.na(fluxha) & !is.na(fluxnii1) & !is.na(fluxnii2) & !is.na(fluxoiii1) & !is.na(fluxoiii2) & snhb >= 3 & snha >= 3  & !is.na(eqwha) & !is.na(eqwnii1) & !is.na(eqwnii2) & !is.na(eqwoiii1) & !is.na(eqwoiii2))] #Obtaining the H-alpha surface specific flux.
snr_a <-    snha[which(fluxha!=500 & !is.na(fluxha) & !is.na(fluxnii1) & !is.na(fluxnii2) & !is.na(fluxoiii1) & !is.na(fluxoiii2) & snhb >= 3 & snha >= 3  & !is.na(eqwha) & !is.na(eqwnii1) & !is.na(eqwnii2) & !is.na(eqwoiii1) & !is.na(eqwoiii2))]
sig_a <- sigmaha[which(fluxha!=500 & !is.na(fluxha) & !is.na(fluxnii1) & !is.na(fluxnii2) & !is.na(fluxoiii1) & !is.na(fluxoiii2) & snhb >= 3 & snha >= 3  & !is.na(eqwha) & !is.na(eqwnii1) & !is.na(eqwnii2) & !is.na(eqwoiii1) & !is.na(eqwoiii2))]
ew_a  <-   eqwha[which(fluxha!=500 & !is.na(fluxha) & !is.na(fluxnii1) & !is.na(fluxnii2) & !is.na(fluxoiii1) & !is.na(fluxoiii2) & snhb >= 3 & snha >= 3  & !is.na(eqwha) & !is.na(eqwnii1) & !is.na(eqwnii2) & !is.na(eqwoiii1) & !is.na(eqwoiii2))]
Npix_a  <- 2*fwhmha[which(fluxha!=500 & !is.na(fluxha) & !is.na(fluxnii1) & !is.na(fluxnii2) & !is.na(fluxoiii1) & !is.na(fluxoiii2) & snhb >= 3 & snha >= 3  & !is.na(eqwha) & !is.na(eqwnii1) & !is.na(eqwnii2) & !is.na(eqwoiii1) & !is.na(eqwoiii2))]

#[N II] 6548
print('Extracting [N II] 6548')
fn1    <-  1e-16*fluxnii1[which(fluxha!=500 & !is.na(fluxha) & !is.na(fluxnii1) & !is.na(fluxnii2) & !is.na(fluxoiii1) & !is.na(fluxoiii2) & snhb >= 3 & snha >= 3  & !is.na(eqwha) & !is.na(eqwnii1) & !is.na(eqwnii2) & !is.na(eqwoiii1) & !is.na(eqwoiii2))] #Obtaining the  [N II] 6548 surface specific flux.
snr_n1 <-    snnii1[which(fluxha!=500 & !is.na(fluxha) & !is.na(fluxnii1) & !is.na(fluxnii2) & !is.na(fluxoiii1) & !is.na(fluxoiii2) & snhb >= 3 & snha >= 3  & !is.na(eqwha) & !is.na(eqwnii1) & !is.na(eqwnii2) & !is.na(eqwoiii1) & !is.na(eqwoiii2))]
sig_n1 <- sigmanii1[which(fluxha!=500 & !is.na(fluxha) & !is.na(fluxnii1) & !is.na(fluxnii2) & !is.na(fluxoiii1) & !is.na(fluxoiii2) & snhb >= 3 & snha >= 3  & !is.na(eqwha) & !is.na(eqwnii1) & !is.na(eqwnii2) & !is.na(eqwoiii1) & !is.na(eqwoiii2))]
ew_n1  <-   eqwnii1[which(fluxha!=500 & !is.na(fluxha) & !is.na(fluxnii1) & !is.na(fluxnii2) & !is.na(fluxoiii1) & !is.na(fluxoiii2) & snhb >= 3 & snha >= 3  & !is.na(eqwha) & !is.na(eqwnii1) & !is.na(eqwnii2) & !is.na(eqwoiii1) & !is.na(eqwoiii2))]
Npix_n1 <- 2*fwhmnii1[which(fluxha!=500 & !is.na(fluxha) & !is.na(fluxnii1) & !is.na(fluxnii2) & !is.na(fluxoiii1) & !is.na(fluxoiii2) & snhb >= 3 & snha >= 3  & !is.na(eqwha) & !is.na(eqwnii1) & !is.na(eqwnii2) & !is.na(eqwoiii1) & !is.na(eqwoiii2))]

#[N II] 6583
print('Extracting [N II] 6583.')
fn    <-  1e-16*fluxnii2[which(fluxha!=500 & !is.na(fluxha) & !is.na(fluxnii1) & !is.na(fluxnii2) & !is.na(fluxoiii1) & !is.na(fluxoiii2) & snhb >= 3 & snha >= 3  & !is.na(eqwha) & !is.na(eqwnii1) & !is.na(eqwnii2) & !is.na(eqwoiii1) & !is.na(eqwoiii2))] #Obtaining the [N II] 6583 surface specific flux.
snr_n <-    snnii2[which(fluxha!=500 & !is.na(fluxha) & !is.na(fluxnii1) & !is.na(fluxnii2) & !is.na(fluxoiii1) & !is.na(fluxoiii2) & snhb >= 3 & snha >= 3  & !is.na(eqwha) & !is.na(eqwnii1) & !is.na(eqwnii2) & !is.na(eqwoiii1) & !is.na(eqwoiii2))]
sig_n <- sigmanii2[which(fluxha!=500 & !is.na(fluxha) & !is.na(fluxnii1) & !is.na(fluxnii2) & !is.na(fluxoiii1) & !is.na(fluxoiii2) & snhb >= 3 & snha >= 3  & !is.na(eqwha) & !is.na(eqwnii1) & !is.na(eqwnii2) & !is.na(eqwoiii1) & !is.na(eqwoiii2))]
ew_n  <-   eqwnii2[which(fluxha!=500 & !is.na(fluxha) & !is.na(fluxnii1) & !is.na(fluxnii2) & !is.na(fluxoiii1) & !is.na(fluxoiii2) & snhb >= 3 & snha >= 3  & !is.na(eqwha) & !is.na(eqwnii1) & !is.na(eqwnii2) & !is.na(eqwoiii1) & !is.na(eqwoiii2))]
Npix_n  <- 2*fwhmnii2[which(fluxha!=500 & !is.na(fluxha) & !is.na(fluxnii1) & !is.na(fluxnii2) & !is.na(fluxoiii1) & !is.na(fluxoiii2) & snhb >= 3 & snha >= 3  & !is.na(eqwha) & !is.na(eqwnii1) & !is.na(eqwnii2) & !is.na(eqwoiii1) & !is.na(eqwoiii2))]

#H-beta 4861
print('Extracting H-beta.')
fb    <-  1e-16*fluxhb[which(fluxha!=500 & !is.na(fluxha) & !is.na(fluxnii1) & !is.na(fluxnii2) & !is.na(fluxoiii1) & !is.na(fluxoiii2) & snhb >= 3 & snha >= 3  & !is.na(eqwha) & !is.na(eqwnii1) & !is.na(eqwnii2) & !is.na(eqwoiii1) & !is.na(eqwoiii2))] #Obtaining the [H-beta] surface specific flux.
snr_b <-    snhb[which(fluxha!=500 & !is.na(fluxha) & !is.na(fluxnii1) & !is.na(fluxnii2) & !is.na(fluxoiii1) & !is.na(fluxoiii2) & snhb >= 3 & snha >= 3  & !is.na(eqwha) & !is.na(eqwnii1) & !is.na(eqwnii2) & !is.na(eqwoiii1) & !is.na(eqwoiii2))]
sig_b <- sigmahb[which(fluxha!=500 & !is.na(fluxha) & !is.na(fluxnii1) & !is.na(fluxnii2) & !is.na(fluxoiii1) & !is.na(fluxoiii2) & snhb >= 3 & snha >= 3  & !is.na(eqwha) & !is.na(eqwnii1) & !is.na(eqwnii2) & !is.na(eqwoiii1) & !is.na(eqwoiii2))]
ew_b  <-   eqwhb[which(fluxha!=500 & !is.na(fluxha) & !is.na(fluxnii1) & !is.na(fluxnii2) & !is.na(fluxoiii1) & !is.na(fluxoiii2) & snhb >= 3 & snha >= 3  & !is.na(eqwha) & !is.na(eqwnii1) & !is.na(eqwnii2) & !is.na(eqwoiii1) & !is.na(eqwoiii2))]
Npix_b  <- 2*fwhmhb[which(fluxha!=500 & !is.na(fluxha) & !is.na(fluxnii1) & !is.na(fluxnii2) & !is.na(fluxoiii1) & !is.na(fluxoiii2) & snhb >= 3 & snha >= 3  & !is.na(eqwha) & !is.na(eqwnii1) & !is.na(eqwnii2) & !is.na(eqwoiii1) & !is.na(eqwoiii2))]

#[O III] 4959
print('Extracting [O III] 4959')
fo1    <-  1e-16*fluxoiii1[which(fluxha!=500 & !is.na(fluxha) & !is.na(fluxnii1) & !is.na(fluxnii2) & !is.na(fluxoiii1) & !is.na(fluxoiii2) & snhb >= 3 & snha >= 3  & !is.na(eqwha) & !is.na(eqwnii1) & !is.na(eqwnii2) & !is.na(eqwoiii1) & !is.na(eqwoiii2))] #Obtaining the [O III] 4959 surface specific flux.
snr_o1 <-    snoiii1[which(fluxha!=500 & !is.na(fluxha) & !is.na(fluxnii1) & !is.na(fluxnii2) & !is.na(fluxoiii1) & !is.na(fluxoiii2) & snhb >= 3 & snha >= 3  & !is.na(eqwha) & !is.na(eqwnii1) & !is.na(eqwnii2) & !is.na(eqwoiii1) & !is.na(eqwoiii2))]
sig_o1 <- sigmaoiii1[which(fluxha!=500 & !is.na(fluxha) & !is.na(fluxnii1) & !is.na(fluxnii2) & !is.na(fluxoiii1) & !is.na(fluxoiii2) & snhb >= 3 & snha >= 3  & !is.na(eqwha) & !is.na(eqwnii1) & !is.na(eqwnii2) & !is.na(eqwoiii1) & !is.na(eqwoiii2))]
ew_o1  <-   eqwoiii1[which(fluxha!=500 & !is.na(fluxha) & !is.na(fluxnii1) & !is.na(fluxnii2) & !is.na(fluxoiii1) & !is.na(fluxoiii2) & snhb >= 3 & snha >= 3  & !is.na(eqwha) & !is.na(eqwnii1) & !is.na(eqwnii2) & !is.na(eqwoiii1) & !is.na(eqwoiii2))]
Npix_o1 <- 2*fwhmoiii1[which(fluxha!=500 & !is.na(fluxha) & !is.na(fluxnii1) & !is.na(fluxnii2) & !is.na(fluxoiii1) & !is.na(fluxoiii2) & snhb >= 3 & snha >= 3  & !is.na(eqwha) & !is.na(eqwnii1) & !is.na(eqwnii2) & !is.na(eqwoiii1) & !is.na(eqwoiii2))]

#[O III] 5007
print('Extracting [O III] 5007.')
fo    <-  1e-16*fluxoiii2[which(fluxha!=500 & !is.na(fluxha) & !is.na(fluxnii1) & !is.na(fluxnii2) & !is.na(fluxoiii1) & !is.na(fluxoiii2) & snhb >= 3 & snha >= 3  & !is.na(eqwha) & !is.na(eqwnii1) & !is.na(eqwnii2) & !is.na(eqwoiii1) & !is.na(eqwoiii2))] #Obtaining the [O III] 5007 surface specific flux.
snr_o <-    snoiii2[which(fluxha!=500 & !is.na(fluxha) & !is.na(fluxnii1) & !is.na(fluxnii2) & !is.na(fluxoiii1) & !is.na(fluxoiii2) & snhb >= 3 & snha >= 3  & !is.na(eqwha) & !is.na(eqwnii1) & !is.na(eqwnii2) & !is.na(eqwoiii1) & !is.na(eqwoiii2))]
sig_o <- sigmaoiii2[which(fluxha!=500 & !is.na(fluxha) & !is.na(fluxnii1) & !is.na(fluxnii2) & !is.na(fluxoiii1) & !is.na(fluxoiii2) & snhb >= 3 & snha >= 3  & !is.na(eqwha) & !is.na(eqwnii1) & !is.na(eqwnii2) & !is.na(eqwoiii1) & !is.na(eqwoiii2))]
ew_o  <-   eqwoiii2[which(fluxha!=500 & !is.na(fluxha) & !is.na(fluxnii1) & !is.na(fluxnii2) & !is.na(fluxoiii1) & !is.na(fluxoiii2) & snhb >= 3 & snha >= 3  & !is.na(eqwha) & !is.na(eqwnii1) & !is.na(eqwnii2) & !is.na(eqwoiii1) & !is.na(eqwoiii2))]  
Npix_o  <- 2*fwhmoiii2[which(fluxha!=500 & !is.na(fluxha) & !is.na(fluxnii1) & !is.na(fluxnii2) & !is.na(fluxoiii1) & !is.na(fluxoiii2) & snhb >= 3 & snha >= 3  & !is.na(eqwha) & !is.na(eqwnii1) & !is.na(eqwnii2) & !is.na(eqwoiii1) & !is.na(eqwoiii2))]

###################################################

##         FLUX ERROR DETERMINATION              ##

###################################################

# Continuum rms.
rmsc_a  <- 1e-16/(sig_a*sqrt(2*pi)*snr_a)
rmsc_n1 <- 1e-16/(sig_n1*sqrt(2*pi)*snr_n1)
rmsc_n  <- 1e-16/(sig_n*sqrt(2*pi)*snr_n)
rmsc_b  <- 1e-16/(sig_b*sqrt(2*pi)*snr_b)
rmsc_o1 <- 1e-16/(sig_o1*sqrt(2*pi)*snr_o1)
rmsc_o  <- 1e-16/(sig_o*sqrt(2*pi)*snr_o)

# Spectral dispersion (Angstrom/pixel)
D <- 1

# Flux error in the reddened emission line
rms_a  <-  rmsc_a*D*sqrt((2*Npix_a)+(ew_a/D))
rms_n1 <- rmsc_n1*D*sqrt((2*Npix_n1)+(ew_n1/D))
rms_n  <-  rmsc_n*D*sqrt((2*Npix_n)+(ew_n/D))
rms_b  <-  rmsc_b*D*sqrt((2*Npix_b)+(ew_b/D))
rms_o1 <- rmsc_o1*D*sqrt((2*Npix_o1)+(ew_o1/D))
rms_o  <-  rmsc_o*D*sqrt((2*Npix_o)+(ew_o/D))


##Convert from surface specific intensity to standard specific intensity [erg cm^-2 s^-1 A^-1 sr^-1].

#print('Converting surface specific flux  units to specific intensities.')
#fa   <-(1e-16)*fa  /2.3504e-11 
#fn   <-(1e-16)*fn  /2.3504e-11 
#fb   <-(1e-16)*fb  /2.3504e-11 
#fo   <-(1e-16)*fo  /2.3504e-11 

##Multiply all values by 1e-16 (CALIFA reference level) Comment this section if using the previous conversion to specific intensity.

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

id3 <- id2[which(!is.na(rms_a) & !is.na(rms_n1) & !is.na(rms_n) & !is.na(rms_b) & !is.na(rms_o1) & !is.na(rms_o))]

#####################################################

##                    DEREDDENING                 ##  

#####################################################

print('First dereddening.')

Rv <- 3.1
l <-  c(6563,6548,6583,4861,4959,5007,3729) #Lines we are going to unredden, in Angstrom.
l <- l*1e-4 #Convert Angstrom to micrometres.
xl <- 1/l
yl <- xl - 1.82

a_x <- 1 + (0.17699*yl) - (0.50447*yl^2) - (0.02427*yl^3) + (0.72085*yl^4) + (0.01979*yl^5) - (0.77530*yl^6) + (0.32999*yl^7)
b_x <- (1.41338*yl) + (2.28305*yl^2) + (1.07233*yl^3) - (5.38434*yl^4) - (0.62251*yl^5) + (5.30260*yl^6) - (2.09002*yl^7)


Cbeta1 <- seq(1:length(id3))
rIb   <- seq(1:length(id3))
rIa   <- seq(1:length(id3))
rIo1  <- seq(1:length(id3))
rIo   <- seq(1:length(id3))
rIn1  <- seq(1:length(id3))
rIn   <- seq(1:length(id3))

rms_Cbeta1 <- seq(1:length(id3))
rms_rIb   <- seq(1:length(id3))
rms_rIa   <- seq(1:length(id3))
rms_rIo1  <- seq(1:length(id3))
rms_rIo   <- seq(1:length(id3))
rms_rIn1  <- seq(1:length(id3))
rms_rIn   <- seq(1:length(id3))

##Claculate the extinction law.

for(j in 1:length(id3)){
  
  f_l <- seq(1:length(l))  #CCM89 extinction law, normalized to V-band.
  fb_l <- seq(1:length(l)) #CCM89 extinction law normalized to H-beta.
  
  for(k in 1:length(l)){
    f_l[k] <- a_x[k] + (b_x[k] / Rv)
  }
  fb_l <- f_l/f_l[4]
  
  ##Find constant c(H-beta). 
  Cbeta1[j] <- (log10(2.86) - log10(fa2[j]/fb2[j]))/(fb_l[1] - 1)
  #rms_Cbeta1[j] <- sqrt(((-1/(fa2[j]*(fb_l[1]-1)*log(10)))^2 * rms_a2[j]^2) + ((-1/(fb2[j]*(fb_l[1]-1)*log(10)))^2 * rms_b2[j]^2))
  rms_Cbeta1[j] <-  sqrt(((-1/((fb_l[1]-1)*fa2[j]*log(10)))*rms_a2[j]))^2 + ((1/((fb_l[1]-1)*fb2[j]*log(10)))*rms_b2[j])^2)
z
  ##Deredden each line normalized to I_H-beta.
  rIb[j]     <- ((fb2[j])/fb2[j])*10^(Cbeta1[j]*(fb_l[4]-1))
  rms_rIb[j] <- sqrt(((10^(Cbeta1[j]*(fb_l[4]-1))/fb2[j])^2 * rms_b2[j]^2) + (((-fb2[j]*10^(Cbeta1[j]*(fb_l[4]-1)))/fb2[j]^2)^2 * rms_b2[j]^2) + (((fb2[j]*(fb_l[4]-1)*log(10)*10^(Cbeta1[j]*(fb_l[4]-1)))/fb2[j])^2 * rms_Cbeta1[j]^2))
  
  rIa[j]  <- ((fa2[j])/fb2[j])*10^(Cbeta1[j]*(fb_l[1]-1))
  rms_rIa[j] <- sqrt(((10^(Cbeta1[j]*(fb_l[1]-1))/fb2[j])^2 * rms_a2[j]^2) + (((-fa2[j]*10^(Cbeta1[j]*(fb_l[1]-1)))/fb2[j]^2)^2 * rms_b2[j]^2) + (((fa2[j]*(fb_l[1]-1)*log(10)*10^(Cbeta1[j]*(fb_l[1]-1)))/fb2[j])^2 * rms_Cbeta1[j]^2))
  
  rIn1[j] <- ((fn12[j])/fb2[j])*10^(Cbeta1[j]*(fb_l[2]-1))
  rms_rIn1[j] <- sqrt(((10^(Cbeta1[j]*(fb_l[2]-1))/fb2[j])^2 * rms_n12[j]^2) + (((-fn12[j]*10^(Cbeta1[j]*(fb_l[2]-1)))/fb2[j]^2)^2 * rms_b2[j]^2) + (((fn12[j]*(fb_l[2]-1)*log(10)*10^(Cbeta1[j]*(fb_l[2]-1)))/fb2[j])^2 * rms_Cbeta1[j]^2))
  
  rIn[j]  <- ((fn2[j])/fb2[j])*10^(Cbeta1[j]*(fb_l[3]-1))
  rms_rIn[j] <- sqrt(((10^(Cbeta1[j]*(fb_l[3]-1))/fb2[j])^2 * rms_n2[j]^2) + (((-fn2[j]*10^(Cbeta1[j]*(fb_l[3]-1)))/fb2[j]^2)^2 * rms_b2[j]^2) + (((fn2[j]*(fb_l[3]-1)*log(10)*10^(Cbeta1[j]*(fb_l[3]-1)))/fb2[j])^2 * rms_Cbeta1[j]^2))
  
  rIo1[j] <- ((fo12[j])/fb2[j])*10^(Cbeta1[j]*(fb_l[5]-1))
  rms_rIo1[j] <- sqrt(((10^(Cbeta1[j]*(fb_l[5]-1))/fb2[j])^2 * rms_o12[j]^2) + (((-fo12[j]*10^(Cbeta1[j]*(fb_l[5]-1)))/fb2[j]^2)^2 * rms_b2[j]^2) + (((fo12[j]*(fb_l[5]-1)*log(10)*10^(Cbeta1[j]*(fb_l[5]-1)))/fb2[j])^2 * rms_Cbeta1[j]^2))
  
  rIo[j]  <- ((fo2[j])/fb2[j])*10^(Cbeta1[j]*(fb_l[6]-1))
  rms_rIo[j] <- sqrt(((10^(Cbeta1[j]*(fb_l[6]-1))/fb2[j])^2 * rms_o2[j]^2) + (((-fo2[j]*10^(Cbeta1[j]*(fb_l[6]-1)))/fb2[j]^2)^2 * rms_b2[j]^2) + (((fo2[j]*(fb_l[6]-1)*log(10)*10^(Cbeta1[j]*(fb_l[6]-1)))/fb2[j])^2 * rms_Cbeta1[j]^2))
  
}

########################################################################

##			LINE RATIOS			              ##

########################################################################

print('Calculating BPT-NII line ratios.')

#[N II]6584 / H-alpha 
na1 <- log10(rIn/rIa)
#[O III] 5007 / H-beta
ob1 <- log10(rIo/rIb)

na  <- na1[which(!is.na(na1) & !is.na(ob1))]
ob  <- ob1[which(!is.na(na1) & !is.na(ob1))]



########################################

##			BPT CLASSIFICATION			      ##

########################################

print('BPT-NII classification.')

##Separating by object type  
kewley<-function(x) (0.61/(x-0.47))+1.19 #Kewley 2001 maximum starburst line.
kauffmann<-function(x) (0.61/(x-0.05))+1.3 #Kauffman 2003 mixing line.
cid<-function(x) 1.01*x+0.48 #Cid Fernandes 2010 Sy/LINER division line.
#new <- function(x) (3.45367*x^2) + (3.82356*x) + 0.18961
new <- function(x) (1.80282*x^2) + (2.58101*x) - 0.25599
#modif <- function(x) (m[i]*x)+b[i] #Proposed modified division for Sy/LINER.

type1 <- seq(1,length(na),1)
for(jj in 1:length(na)){
  if(ob[jj] <= kauffmann(na[jj]) & na[jj] < 0.05) {type1[jj] <- "blue"} #Star Forming spaxels.
  #else if(ob[jj] <= modif(na[jj])) {type[jj] <- "orange"} #LINER-like proposed division, uncomment to use it and comment the cid line.
  else if(ob[jj] <= kewley(na[jj]) & na[jj]< 0.47){type1[jj] <- "green"} #Transition Object.
  else if(ob[jj] <= cid(na[jj])) {type1[jj] <- "orange"} #LINER-like, uncomment to use Cid Fernandes division and comment the modif line.
  else {type1[jj] <- "red"} # Seyfert.
} 

###############################################################

##                        DEREDDENING                        ##
##		Dereddening again for AGN spaxels with                 ## 
##		Ha/Hb = 3.1                                            ##
##              Using Cardelli et al. 1989                   ##

###############################################################

print('Second dereddening. Correction for AGN spaxels.')

ID <-  id3[which(!is.na(na1) & !is.na(ob1))]

type <- type1[which(!is.na(na1) & !is.na(ob1))]

fa3  <-  fa2[which(!is.na(na1) & !is.na(ob1))]
fb3  <-  fb2[which(!is.na(na1) & !is.na(ob1))]
fn13 <- fn12[which(!is.na(na1) & !is.na(ob1))]
fn3  <-  fn2[which(!is.na(na1) & !is.na(ob1))]
fo13 <- fo12[which(!is.na(na1) & !is.na(ob1))]
fo3  <-  fo2[which(!is.na(na1) & !is.na(ob1))]

rms_fa3  <-  rms_a2[which(!is.na(na1) & !is.na(ob1))]
rms_fb3  <-  rms_b2[which(!is.na(na1) & !is.na(ob1))]
rms_fn13 <- rms_n12[which(!is.na(na1) & !is.na(ob1))]
rms_fn3  <-  rms_n2[which(!is.na(na1) & !is.na(ob1))]
rms_fo13 <- rms_o12[which(!is.na(na1) & !is.na(ob1))]
rms_fo3  <-  rms_o2[which(!is.na(na1) & !is.na(ob1))]

rIa3  <-  rIa[which(!is.na(na1) & !is.na(ob1))]
rIb3  <-  rIb[which(!is.na(na1) & !is.na(ob1))]
rIn13 <- rIn1[which(!is.na(na1) & !is.na(ob1))]
rIn3  <-  rIn[which(!is.na(na1) & !is.na(ob1))]
rIo13 <- rIo1[which(!is.na(na1) & !is.na(ob1))]
rIo3  <-  rIo[which(!is.na(na1) & !is.na(ob1))]

rms_rIa3  <-  rms_rIa[which(!is.na(na1) & !is.na(ob1))]
rms_rIb3  <-  rms_rIb[which(!is.na(na1) & !is.na(ob1))]
rms_rIn13 <- rms_rIn1[which(!is.na(na1) & !is.na(ob1))]
rms_rIn3  <-  rms_rIn[which(!is.na(na1) & !is.na(ob1))]
rms_rIo13 <- rms_rIo1[which(!is.na(na1) & !is.na(ob1))]
rms_rIo3  <-  rms_rIo[which(!is.na(na1) & !is.na(ob1))]

Cbeta <- Cbeta1[which(!is.na(na1) & !is.na(ob1))]
rms_Cbeta <- rms_Cbeta1[which(!is.na(na1) & !is.na(ob1))]

for(kk in 1:length(ID)){
 if(type[kk] == "red"){

   ## Find constant c(H-beta).
   Cbeta[kk] <- (log10(3.1) - log10(fa3[kk]/fb3[kk]))/(fb_l[1] - 1)
   rms_Cbeta[kk] <- sqrt(((-1/(fa3[kk]*(fb_l[1]-1)*log(10)))^2 * rms_fa3[kk]^2) + ((-1/(fb3[kk]*(fb_l[1]-1)*log(10)))^2 * rms_fb3[kk]^2))
   
  ## Deredden each line normalized to I_H-beta.
  ## Multiply each by 100 if they are used directly.

   rIb3[kk] <- ((fb3[kk])/fb3[kk])*10^(Cbeta[kk]*(fb_l[4]-1))
   rms_rIb3[kk] <- sqrt(((10^(Cbeta[kk]*(fb_l[4]-1))/fb3[kk])^2 * rms_fb3[kk]^2) + (((-fb3[kk]*10^(Cbeta[kk]*(fb_l[4]-1)))/fb3[kk]^2)^2 * rms_fb3[kk]^2) + (((fb3[kk]*(fb_l[4]-1)*log(10)*10^(Cbeta[kk]*(fb_l[4]-1)))/fb3[kk])^2 * rms_Cbeta[kk]^2))
   
   rIa3[kk] <- ((fa3[kk])/fb3[kk])*10^(Cbeta[kk]*(fb_l[1]-1))
   rms_rIa3[kk] <- sqrt(((10^(Cbeta[kk]*(fb_l[1]-1))/fb3[kk])^2 * rms_fa3[kk]^2) + (((-fa3[kk]*10^(Cbeta[kk]*(fb_l[1]-1)))/fb3[kk]^2)^2 * rms_fb3[kk]^2) + (((fa3[kk]*(fb_l[1]-1)*log(10)*10^(Cbeta[kk]*(fb_l[1]-1)))/fb3[kk])^2 * rms_Cbeta[kk]^2))
   
   rIn13[kk] <- ((fn13[kk])/fb3[kk])*10^(Cbeta[kk]*(fb_l[2]-1))
   rms_rIn13[kk] <- sqrt(((10^(Cbeta[kk]*(fb_l[2]-1))/fb3[kk])^2 * rms_fn13[kk]^2) + (((-fn13[kk]*10^(Cbeta[kk]*(fb_l[2]-1)))/fb3[kk]^2)^2 * rms_fb3[kk]^2) + (((fn13[kk]*(fb_l[2]-1)*log(10)*10^(Cbeta[kk]*(fb_l[2]-1)))/fb3[kk])^2 * rms_Cbeta[kk]^2))
   
   rIn3[kk] <- ((fn3[kk])/fb3[kk])*10^(Cbeta[kk]*(fb_l[3]-1))
   rms_rIn3[kk] <- sqrt(((10^(Cbeta[kk]*(fb_l[3]-1))/fb3[kk])^2 * rms_fn3[kk]^2) + (((-fn3[kk]*10^(Cbeta[kk]*(fb_l[3]-1)))/fb3[kk]^2)^2 * rms_fb3[kk]^2) + (((fn3[kk]*(fb_l[3]-1)*log(10)*10^(Cbeta[kk]*(fb_l[kk]-1)))/fb3[kk])^2 * rms_Cbeta[kk]^2))
   
   rIo13[kk] <- ((fo13[kk])/fb3[kk])*10^(Cbeta[kk]*(fb_l[5]-1))
   rms_rIo13[kk] <- sqrt(((10^(Cbeta[kk]*(fb_l[5]-1))/fb3[kk])^2 * rms_fo13[kk]^2) + (((-fo13[kk]*10^(Cbeta[kk]*(fb_l[5]-1)))/fb3[kk]^2)^2 * rms_fb3[kk]^2) + (((fo13[kk]*(fb_l[5]-1)*log(10)*10^(Cbeta[kk]*(fb_l[5]-1)))/fb3[kk])^2 * rms_Cbeta[kk]^2))
   
   rIo3[kk] <- ((fo3[kk])/fb3[kk])*10^(Cbeta[kk]*(fb_l[6]-1))
   rms_rIo3[kk] <- sqrt(((10^(Cbeta[kk]*(fb_l[6]-1))/fb3[kk])^2 * rms_fo3[kk]^2) + (((-fo3[kk]*10^(Cbeta[kk]*(fb_l[6]-1)))/fb3[kk]^2)^2 * rms_fb3[kk]^2) + (((fo3[kk]*(fb_l[6]-1)*log(10)*10^(Cbeta[kk]*(fb_l[6]-1)))/fb3[kk])^2 * rms_Cbeta[kk]^2))
   
 }
}

# Determination of the dereddened surface specific fluxes.

AHa <- seq(1,length(ID),1)
Fa  <- seq(1,length(ID),1)
Fb  <- seq(1,length(ID),1)
Fn1 <- seq(1,length(ID),1)
Fn  <- seq(1,length(ID),1)
Fo1 <- seq(1,length(ID),1)
Fo  <- seq(1,length(ID),1)

rms_AHa <- seq(1,length(ID),1)
rms_Fa  <- seq(1,length(ID),1)
rms_Fb  <- seq(1,length(ID),1)
rms_Fn1 <- seq(1,length(ID),1)
rms_Fn  <- seq(1,length(ID),1)
rms_Fo1 <- seq(1,length(ID),1)
rms_Fo  <- seq(1,length(ID),1)


for (ww in 1:length(ID)){
 if(type[ww] == "red"){
  AHa[ww] <- (f_l[1]/(-0.4*(f_l[1]-f_l[4])))*log10((fa3[ww]/fb3[ww])/3.1) # H-alpha extinction (magnitudes, V-band normalized). Original Cardelli et al. values.
  # AHa[ww] <- (2.53/(-0.4*(2.53-3.61)))*log10((fa[ww]/fb[ww])/3.1) #H-alpha extinction (magnitudes, V-band normalized). Catalan-Torrecilla values. Almost the same result.
  rms_AHa[ww] <- sqrt(((1/(fa3[ww]*log(10)))*rms_fa3[ww])^2 + ((-1/(fb3[ww]*log(10)))*rms_fb3[ww])^2)
    
  Fa[ww]      <- fa3[ww]*10^(0.4*AHa[ww]) #H-alpha dereddened flux.
  rms_Fa[ww]  <- sqrt(((10^(0.4*AHa[ww]))*rms_fa3[ww])^2 + ((0.4*log(10)*10^(0.4*AHa[ww]))*rms_AHa[ww])^2)
  
  Fb[ww]      <- Fa[ww]/3.1      #H-beta dereddened flux.
  rms_Fb[ww]  <- rms_Fa[ww]/3.1
  
  Fn1[ww]     <- rIn13[ww]*Fb[ww] #[N II] 6548 dereddened flux.
  rms_Fn1[ww] <- sqrt((Fb[ww]*rms_rIn13[ww])^2 + (rIn13[ww]*rms_Fb[ww])^2)
  
  Fn[ww]  <- rIn3[ww]*Fb[ww]  #[N II] 6583 dereddenend flux.
  rms_Fn[ww] <- sqrt((Fb[ww]*rms_rIn3[ww])^2 + (rIn3[ww]*rms_Fb[ww])^2)
  
  Fo1[ww] <- rIo13[ww]*Fb[ww] #[O III] 4959 dereddened flux.
  rms_Fo1[ww] <- sqrt((Fb[ww]*rms_rIo13[ww])^2 + (rIo13[ww]*rms_Fb[ww])^2)
  
  Fo[ww]  <- rIo3[ww]*Fb[ww]  #[O III] 5007 dereddened flux.
  rms_Fo[ww] <- sqrt((Fb[ww]*rms_rIo3[ww])^2 + (rIo3[ww]*rms_Fb[ww])^2)
  
 }
 else {
   AHa[ww] <- (f_l[1]/(-0.4*(f_l[1]-f_l[4])))*log10((fa3[ww]/fb3[ww])/2.86) # H-alpha extinction (magnitudes, V-band normalized). Original Cardelli et al. values.
   # AHa[ww] <- (2.53/(-0.4*(2.53-3.61)))*log10((fa[ww]/fb[ww])/3.1) #H-alpha extinction (magnitudes, V-band normalized). Catalan-Torrecilla values. Almost the same result.
   rms_AHa[ww] <- sqrt(((f_l[1]/(-0.4*(f_l[1]-f_l[4])))*(1/(fa3[ww]*log(10)))*rms_fa3[ww])^2 + ((f_l[1]/(-0.4*(f_l[1]-f_l[4])))*(-1/(fb3[ww]*log(10)))*rms_fb3[ww])^2)
   #rms_AHa[ww] <- sqrt(((f_l[1]/(-0.4*(f_l[1]-f_l[4])))*(1/(fa3[ww]*log(10)))*rms_fa3[ww])^2 + ((f_l[1]/(-0.4*(f_l[1]-f_l[4])))*(log(10^0.4)*exp(log(10^0.4)*AHa[ww]))*rms_AHa[ww])^2)
   
   Fa[ww]      <- fa3[ww]*10^(0.4*AHa[ww]) #H-alpha dereddened flux.
   #rms_Fa[ww]  <- sqrt(((10^(0.4*AHa[ww]))*rms_fa3[ww])^2 + ((0.4*log(10)*10^(0.4*AHa[ww]))*rms_AHa[ww])^2)
   #rms_Fa[ww]  <- sqrt(((10^(0.4*AHa[ww]))*rms_fa3[ww])^2 + (fa3[ww]*(0.4*log(10)*10^(0.4*AHa[ww]))*rms_AHa[ww])^2)
   rms_Fa[ww]  <- sqrt(((10^(0.4*AHa[ww]))*rms_fa3[ww])^2 + (fa3[ww]*(log(10^0.4)*exp(AHa[ww]*10^0.4))*rms_AHa[ww])^2)
   
   Fb[ww]      <- Fa[ww]/2.86     #H-beta dereddened flux.
   rms_Fb[ww]  <- sqrt((rms_Fa[ww]/2.86)^2)
   
   Fn1[ww]     <- rIn13[ww]*Fb[ww] #[N II] 6548 dereddened flux.
   rms_Fn1[ww] <- sqrt((Fb[ww]*rms_rIn13[ww])^2 + (rIn13[ww]*rms_Fb[ww])^2)
   
   Fn[ww]  <- rIn3[ww]*Fb[ww]  #[N II] 6583 dereddenend flux.
   rms_Fn[ww] <- sqrt((Fb[ww]*rms_rIn3[ww])^2 + (rIn3[ww]*rms_Fb[ww])^2)
   
   Fo1[ww] <- rIo13[ww]*Fb[ww] #[O III] 4959 dereddened flux.
   rms_Fo1[ww] <- sqrt((Fb[ww]*rms_rIo13[ww])^2 + (rIo13[ww]*rms_Fb[ww])^2)
   
   Fo[ww]  <- rIo3[ww]*Fb[ww]  #[O III] 5007 dereddened flux.
   rms_Fo[ww] <- sqrt((Fb[ww]*rms_rIo3[ww])^2 + (rIo3[ww]*rms_Fb[ww])^2)
   
 }

}


##############################################################################
##Save data to files.

print('Saving fluxes to data file.')
resume <- data.frame(ID,AHa,rms_AHa,Cbeta,rms_Cbeta,Fa,rms_Fa,Fb,rms_Fb,Fn1,rms_Fn1,Fn,rms_Fn,Fo1,rms_Fo1,Fo,rms_Fo)
tabla <- str_c(galaxy[i],"/Cardelli_Base_Fluxes.dat")
write.table(resume, tabla, sep="\t",quote=FALSE)

}
