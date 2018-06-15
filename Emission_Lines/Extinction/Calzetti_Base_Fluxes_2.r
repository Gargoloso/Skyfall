#############################################################################
#############################################################################

##   This calculates extinction corrected fluxes using the Calzetti et al. ##
##   (200) extinction law.                                                 ##    

##   June 10, 2018 by A. Robleto-Orús                                       ##
#############################################################################
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

print('Extracting H-alpha and [N II] 6583  data.')  
path_ha <- str_c(galaxy[i],"/lis_ha.res") #creates path to Ha and [N II] data file for each galaxy
data0 <- read.table(path_ha, header=TRUE)

print('Extracting H-beta and [O III] 5007 data.')
path_hb <- str_c(galaxy[i],"/lis_hb.res") #creates path to H-beta and [OIII] data file for each galaxy
data3 <- read.table(path_hb, header=TRUE)
         
##Merge data

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

ID <- id2[which(!is.na(rms_a) & !is.na(rms_n1) & !is.na(rms_n) & !is.na(rms_b) & !is.na(rms_o1) & !is.na(rms_o))]

#####################################################

##                    DEREDDENING                 ##  

#####################################################

print('Dereddening.')

Rv <- 4.05 #Value for star-forming or high redshift galaxies, in Calzetti et al=. (2000)
l <-  c(6563,6548,6583,4861,4959,5007) #Lines we are going to deredden, in Angstrom.
l <- l*1e-4 #Convert Angstrom to micrometres.

#Calzetti et al. (2000) extinction law.
k_l1 <- (2.659*(-2.156+(1.509/l)-(0.198/l^2)+(0.011/l^3)))+Rv #For 0.12 to 0.63 micrometres
k_l2 <- (2.659*(-1.857+(1.040/l)))+Rv #For 0.63 to 2.20 micrometres.

Rab <- (fa2/fb2)/2.86 #Ratio of attenuated Halpha/Hbeta over non-attenuated one.
rms_Rab <-sqrt((rms_a2/(2.86*fb2))^2 + ((-fa2*rms_b2)/(2.86*fb2^2))^2)

EBV <- log10(Rab)/(0.4*1.163) #Colour excess E(B-V) in magnitudes.
rms_EBV <- (rms_Rab/(0.4*1.163*Rab*log(10)))

Fa  <-  fa2*10^(0.4*EBV*k_l2[1]) #H-alpha dereddening
rms_Fa <-  sqrt((rms_a2*10^(0.4*EBV*k_l2[1]))^2 + (rms_EBV*log(10)*fa2*exp(log(10)*EBV*k_l2[1]))^2)

Fb  <-  fb2*10^(0.4*EBV*k_l1[4]) #H-beta dereddening
rms_Fb <-  sqrt((rms_b2*10^(0.4*EBV*k_l1[4]))^2 + (rms_EBV*log(10)*fb2*exp(log(10)*EBV*k_l1[4]))^2)

Fn1 <- fn12*10^(0.4*EBV*k_l2[2]) #[N II] 6548 dereddening 
rms_Fn1 <-  sqrt((rms_n12*10^(0.4*EBV*k_l2[2]))^2 + (rms_EBV*log(10)*fn12*exp(log(10)*EBV*k_l2[2]))^2)

Fn  <-  fn2*10^(0.4*EBV*k_l2[3]) #[N II] 6583 dereddening
rms_Fn <-  sqrt((rms_n2*10^(0.4*EBV*k_l2[3]))^2 + (rms_EBV*log(10)*fn2*exp(log(10)*EBV*k_l2[3]))^2)

Fo1 <- fo12*10^(0.4*EBV*k_l1[5]) #[O III] 4959 dereddening
rms_Fo1 <-  sqrt((rms_o12*10^(0.4*EBV*k_l1[5]))^2 + (rms_EBV*log(10)*fo12*exp(log(10)*EBV*k_l1[5]))^2)

Fo  <-  fo2*10^(0.4*EBV*k_l1[6]) #[O III] 6583 dereddening
rms_Fo <-  sqrt((rms_o2*10^(0.4*EBV*k_l1[6]))^2 + (rms_EBV*log(10)*fo2*exp(log(10)*EBV*k_l1[6]))^2)

AHa <- k_l2[1]*EBV #H-alpha extinction (in magnitudes)
rms_AHa <- k_l2[1]*rms_EBV

##############################################################################
##Save data to files

print('Saving fluxes to data file.')
#resume <- data.frame(ID,AHa,EBV,Fa,Fb,Fn1,Fn,Fo1,Fo)
resume <- data.frame(ID,AHa,rms_AHa,Fa,rms_Fa,Fb,rms_Fb,Fn1,rms_Fn1,Fn,rms_Fn,Fo1,rms_Fo1,Fo,rms_Fo)
tabla <- str_c(galaxy[i],"/Calzetti_Base_Fluxes.dat")
write.table(resume, tabla, sep="\t",quote=FALSE)

}
