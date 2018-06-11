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


##Extracting coordinates from ID

print('Extracting spaxels IDs.')
id <-    id[which(fluxha!=500)]

##Extracting line surface specific fluxes [1e-16 erg cm^-2 s^-1 A^-1 arcsec^-2]

print('Extracting emission lines.')

#H-alpha 6563
print('Extracting H-alpha.')
fa <-    fluxha[which(fluxha!=500)] #Obtaining the H-alpha surface specific flux.

#[N II] 6548
print('Extracting [N II] 6548')
fn1 <-   fluxnii1[which(fluxha!=500)] #Obtaining the  [N II] 6548 surface specific flux.
 
#[N II] 6583
print('Extracting [N II] 6583.')
fn <-    fluxnii2[which(fluxha!=500)] #Obtaining the[N II] surface specific flux.

#H-beta 4861
print('Extracting H-beta.')
fb <-    fluxhb[which(fluxha!=500)] #Obtaining the[H-beta] surface specific flux.

#[O III] 4959
print('Extracting [O III] 4959')
fo1 <-   fluxoiii1[which(fluxha!=500)] #Obtaining the [O III] 4959 surface specific flux.

#[O III] 5007
print('Extracting [O III] 5007.')
fo <-    fluxoiii2[which(fluxha!=500)] #Obtaining the [O III] 5007 surface specific flux.

##Convert from surface specific flux. to standard specific intensity [erg cm^-2 s^-1 A^-1 sr^-1]

#print('Converting surface specific flux  units to specific intensities.')
#fa   <-(1e-16)*fa  /2.3504e-11 
#fn   <-(1e-16)*fn  /2.3504e-11 
#fb   <-(1e-16)*fb  /2.3504e-11 
#fo   <-(1e-16)*fo  /2.3504e-11 

##Multiply all values by 1e-16 (CALIFA reference level) Comment this section if using the previous conversion to specific intensity.

fa   <-(1e-16)*fa
fn1  <-(1e-16)*fn1
fn   <-(1e-16)*fn 
fb   <-(1e-16)*fb 
fo1  <-(1e-16)*fo1
fo   <-(1e-16)*fo 

#####################################################

##                    DEREDDENING                 ##  

#####################################################

print('Dereddening.')

Rv <- 4.05 #Value for star-forming or high redshift galaxies, in Calzetti et al=. (2000)
l <-  c(6563,6548,6583,4861,4959,5007,3729) #Lines we are going to deredden, in Angstrom.
l <- l*1e-4 #Convert Angstrom to micrometres.

#Calzetti et al. (2000) extinction law.
k_l1 <- (2.659*(-2.156+(1.509/l)-(0.198/l^2)+(0.011/l^3)))+Rv #For 0.12 to 0.63 micrometres
k_l2 <- (2.659*(-1.857+(1.040/l)))+Rv #For 0.63 to 2.20 micrometres.

Rab <- (fa/fb)/2.86 #Ratio of attenuated Halpha/Hbeta over non-attenuated one.

EBV <- (0.4*1.163)*(log10(10)/log10(Rab)) #Colour excess E(B-V) in magnitudes.

Fa  <-  fa*10^(0.4*EBV*k_l2[1]) #H-alpha dereddening
Fb  <-  fb*10^(0.4*EBV*k_l1[4]) #H-beta dereddening
Fn1 <- fn1*10^(0.4*EBV*k_l2[2]) #[N II] 6548 dereddening 
Fn  <-  fn*10^(0.4*EBV*k_l2[3]) #[N II] 6583 dereddening
Fo1 <- fo1*10^(0.4*EBV*k_l1[5]) #[O III] 4959 dereddening
Fo  <-  fo*10^(0.4*EBV*k_l1[6]) #[O III] 6583 dereddening

AHa <- k_l2[1]*EBV #H-alpha extinction (in magnitudes)

##############################################################################
##Save data to files

print('Saving fluxes to data file.')
resume <- data.frame(id,AHa,EBV,Fa,Fb,Fn1,Fn,Fo1,Fo)
tabla <- str_c(galaxy[i],"/Calzetti_Base_Fluxes.dat")
write.table(resume, tabla, sep="\t",quote=FALSE)

}
