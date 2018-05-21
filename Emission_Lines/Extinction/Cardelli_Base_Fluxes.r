#############################################################################
#############################################################################

##   This calculates extinction corrected fluxes using the Cardelli        ##
##   extinction law.                                                       ##    

##   February 25, 2018 A. Robleto-Orús                                     ##
##   May 15, 2018 Adapted to work with the last version of fitlines.       ##
##                We now only save the final surface specific fluxes but   ##
##                not line ratios. Also Cbeta and AHa are saved.           ##
##                
#############################################################################
#############################################################################


##Clean the workspace
rm(list=ls(all=TRUE))

##Libraries
library("fields")
require("stringr")
library("png")

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
path_ha <- str_c(galaxy[i],"/lis_ha.res") #creates path to Ha and [N II]data file for each galaxy
data0 <- read.table(path_ha, header=TRUE)

print('Extracting H-beta and [O III] 5007 data.')
path_hb <- str_c(galaxy[i],"/lis_hb.res") #creates path to H-beta and [OIII] data file for each galaxy
data3 <- read.table(path_hb, header=TRUE)
         
##Merge data

DATA <- merge(data0,data3,by.x = 1, by.y =1)
attach(DATA)


##Extracting coordinates from ID
## Spaxel's coordinates are in the id in format YYXX, with Y = DEC and X = RA.

print('Extracting spaxel coordinates.')
y <-     (trunc(id[which(fluxha!=500)]/100)  ) #Obtaining Y from id.
x <-     (id[which(fluxha!=500)]-(y*100) ) #Obtaining X from id.
id2 <-    id[which(fluxha!=500)]

##Extracting line surface specific intensities [1e-16 erg cm^-2 s^-1 A^-1 arcsec^-2]

print('Extracting emission lines.')

#H-alpha 6563
print('Extracting H-alpha.')
fa <-    fluxha[which(fluxha!=500)] #Obtaining the H-alpha flux density (10^-16 erg/(s cm^2 A pix))

#[N II] 6548
print('Extracting [N II] 6548')
fn1 <-   fluxnii1[which(fluxha!=500)] #Obtaining the  [N II] 6548 flux density (10^-16/erg(s cm^2 A pix))
 
#[N II] 6583
print('Extracting [N II] 6583.')
fn <-    fluxnii2[which(fluxha!=500)] #Obtaining the[N II] flux density (10^-16 erg/(s cm^2 A pix))

#H-beta 4861
print('Extracting H-beta.')
fb <-    fluxhb[which(fluxha!=500)] #Obtaining the[H-beta] flux density (10^-16 erg/(s cm^2 A pix))

#[O III] 4959
print('Extracting [O III] 4959')
fo1 <-   fluxoiii1[which(fluxha!=500)] #Obtaining the [O III] 4959 flux density (10^-16 erg/(s cm^2 A pix))

#[O III] 5007
print('Extracting [O III] 5007.')
fo <-    fluxoiii2[which(fluxha!=500)] #Obtaining the [O III] 5007 flux density (10^-16 erg/(s cm^2 A pix))

#

##id of central spaxel
#print('Determining central spaxel postion.')
#xc <- trunc(naxis1[i]/2)
#yc <- trunc(naxis2[i]/2)

#X3 <- x - xc
#Y3 <- y - yc


##Convert from surface specific intensity to standard specific intensity [erg cm^-2 s^-1 A^-1 sr^-1]

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

print('First dereddening.')

Rv <- 3.1
l <-  c(6563,6548,6583,4861,4959,5007,3729) #Lines we are going to unredden, in Angstrom.
l <- l*1e-4 #Convert Angstrom to micrometers
xl <- 1/l
yl <- xl - 1.82

a_x <- 1 + (0.17699*yl) - (0.50447*yl^2) - (0.02427*yl^3) + (0.72085*yl^4) + (0.01979*yl^5) - (0.77530*yl^6) + (0.32999*yl^7)
b_x <- (1.41338*yl) + (2.28305*yl^2) + (1.07233*yl^3) -(5.38434*yl^4) - (0.62251*yl^5) + (5.30260*yl^6) - (2.09002*yl^7)


Cbeta <- seq(1:length(y))
rIb   <- seq(1:length(y))
rIa   <- seq(1:length(y))
rIo1  <- seq(1:length(y))
rIo   <- seq(1:length(y))
rIn1  <- seq(1:length(y))
rIn   <- seq(1:length(y))

##Claculate the extinction law

for(j in 1:length(y)){
  
  f_l <- seq(1:length(l))  #CCM89 extinction law, normalized to V-band
  fb_l <- seq(1:length(l)) #Extinction law normalized to H-beta
  
  for(k in 1:length(l)){
    f_l[k] <- a_x[k] + (b_x[k] / Rv)
  }
  fb_l <- f_l/f_l[4]
  
  ##Find constant c(H-beta)  
  Cbeta[j] <- (log10(2.86) - log10(fa[j]/fb[j]))/(fb_l[1] - 1)
  
  ##Deredden each line normalized to I_H-beta
  
  rIb[j]  <- ((fb[j])/fb[j])*10^(Cbeta[j]*(fb_l[4]-1))
  
  rIa[j]  <- ((fa[j])/fb[j])*10^(Cbeta[j]*(fb_l[1]-1))
  
  rIn1[j] <- ((fn1[j])/fb[j])*10^(Cbeta[j]*(fb_l[2]-1))
  
  rIn[j]  <- ((fn[j])/fb[j])*10^(Cbeta[j]*(fb_l[3]-1))
  
  rIo1[j] <- ((fo1[j])/fb[j])*10^(Cbeta[j]*(fb_l[5]-1))
  
  rIo[j]  <- ((fo[j])/fb[j])*10^(Cbeta[j]*(fb_l[6]-1))
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



########################################################################

##			BPT CLASSIFICATION			      ##

########################################################################

print('BPT-NII classification.')

##Separating by object type  
kewley<-function(x) (0.61/(x-0.47))+1.19 #Kewley 2001 maximum starburst line
kauffmann<-function(x) (0.61/(x-0.05))+1.3 #Kauffman 2003 mixing line
cid<-function(x) 1.01*x+0.48 #Cid Fernandes 2010 Sy/LINER division line
#new <- function(x) (3.45367*x^2) + (3.82356*x) + 0.18961
new <- function(x) (1.80282*x^2) + (2.58101*x) - 0.25599
#modif <- function(x) (m[i]*x)+b[i] #Proposed modified division for Sy/LINER

type <- seq(1,length(na),1)
for(jj in 1:length(na)){
  if(ob[jj] <= kauffmann(na[jj]) & na[jj] < 0.05) {type[jj] <- "blue"} #Star Forming spaxels
  #else if(ob[jj] <= modif(na[jj])) {type[jj] <- "orange"} #LINER-like proposed division, uncomment to use it and comment the cid line.
  else if(ob[jj] <= kewley(na[jj]) & na[jj]< 0.47){type[jj] <- "green"} #Transition Object
  else if(ob[jj] <= cid(na[jj])) {type[jj] <- "orange"} #LINER-like, uncomment to use Cid Fernandes division and comment the modif line.
  else {type[jj] <- "red"} # Seyfert
} 

###############################################################

##                        DEREDDENING                        ##
##		Dereddening again for AGN spaxels with 
##		Ha/Hb = 3.1
##              Using Cardelli et al. 1989                   ##

###############################################################

print('Second dereddening. Correction for AGN spaxels.')

id <- id2[which(!is.na(na1) & !is.na(ob1))]
X2 <-   x[which(!is.na(na1) & !is.na(ob1))]
Y2 <-   y[which(!is.na(na1) & !is.na(ob1))]

fa  <-  fa[which(!is.na(na1) & !is.na(ob1))]
fb  <-  fb[which(!is.na(na1) & !is.na(ob1))]
fn1 <- fn1[which(!is.na(na1) & !is.na(ob1))]
fn  <-  fn[which(!is.na(na1) & !is.na(ob1))]
fo1 <- fo1[which(!is.na(na1) & !is.na(ob1))]
fo  <-  fo[which(!is.na(na1) & !is.na(ob1))]


rIa  <-  rIa[which(!is.na(na1) & !is.na(ob1))]
rIb  <-  rIb[which(!is.na(na1) & !is.na(ob1))]
rIn1 <- rIn1[which(!is.na(na1) & !is.na(ob1))]
rIn  <-  rIn[which(!is.na(na1) & !is.na(ob1))]
rIo1 <- rIo1[which(!is.na(na1) & !is.na(ob1))]
rIo  <-  rIo[which(!is.na(na1) & !is.na(ob1))]

Cbeta <- Cbeta[which(!is.na(na1) & !is.na(ob1))]

for(kk in 1:length(X2)){
 if(type[kk] == "red"){

   ##Find constant c(H-beta)
   Cbeta[kk] <- (log10(3.1) - log10(fa[kk]/fb[kk]))/(fb_l[1] - 1)

  ##Deredden each line normalized to I_H-beta
  #Multiply each by 100 if they are used directly.

   rIb[kk] <- ((fb[kk])/fb[kk])*10^(Cbeta[j]*(fb_l[4]-1))

   rIa[kk] <- ((fa[kk])/fb[kk])*10^(Cbeta[j]*(fb_l[1]-1))
   
   rIn1[j] <- ((fn1[j])/fb[j])*10^(Cbeta[j]*(fb_l[2]-1))

   rIn[kk] <- ((fn[kk])/fb[kk])*10^(Cbeta[j]*(fb_l[3]-1))
   
   rIo1[j] <- ((fo1[j])/fb[j])*10^(Cbeta[j]*(fb_l[5]-1))

   rIo[kk] <- ((fo[kk])/fb[kk])*10^(Cbeta[j]*(fb_l[6]-1))
 }
}

#Determination of the dereddened surface specific fluxes.

AHa <- seq(1,length(X2),1)
Fa  <- seq(1,length(X2),1)
Fb  <- seq(1,length(X2),1)
Fn1 <- seq(1,length(X2),1)
Fn  <- seq(1,length(X2),1)
Fo1 <- seq(1,length(X2),1)
Fo  <- seq(1,length(X2),1)

for (ww in 1:length(X2)){
 if(type[ww] == 'red'){
  AHa[ww] <- (f_l[1]/(-0.4*(f_l[1]-f_l[4])))*log10((fa[ww]/fb[ww])/3.1) #H-alpha extinction (magnitudes, V-band normalized). Original Cardelli et al. values.
  #AHa[ww] <- (2.53/(-0.4*(2.53-3.61)))*log10((fa[ww]/fb[ww])/3.1) #H-alpha extinction (magnitudes, V-band normalized). Catalan-Torrecilla values. Almost the same result.
  Fa[ww]  <- fa[ww]*10^(0.4*AHa[ww]) #H-alpha dereddened flux.
  Fb[ww]  <- Fa[ww]/3.1      #H-beta dereddened flux.
  Fn1[ww] <- rIn1[ww]*Fb[ww] #[N II] 6548 dereddened flux.
  Fn[ww]  <- rIn[ww]*Fb[ww]  #[N II] 6583 dereddenend flux.
  Fo1[ww] <- rIo1[ww]*Fb[ww] #[O III] 4959 dereddened flux.
  Fo[ww]  <- rIo[ww]*Fb[ww]  #[O III] 5007 dereddened flux.
 }
 else {
  AHa[ww] <- (f_l[1]/(-0.4*(f_l[1]-f_l[4])))*log10((fa[ww]/fb[ww])/2.86) #H-alpha extinction (magnitudes, V-band normalized). Original Cardelli et al. values.
  #AHa[ww] <- (2.53/(-0.4*(2.53-3.61)))*log10((fa[ww]/fb[ww])/2.86) #H-alpha extinction (magnitudes, V-band normalized). Catalan-Torrecilla values. Almost the same result.
  Fa[ww]  <- fa[ww]*10^(0.4*AHa[ww]) #H-alpha dereddened flux.
  Fb[ww]  <- Fa[ww]/2.86      #H-beta dereddened flux.
  Fn1[ww] <- rIn1[ww]*Fb[ww] #[N II] 6548 dereddened flux.
  Fn[ww]  <- rIn[ww]*Fb[ww]  #[N II] 6583 dereddenend flux.
  Fo1[ww] <- rIo1[ww]*Fb[ww] #[O III] 4959 dereddened flux.
  Fo[ww]  <- rIo[ww]*Fb[ww]  #[O III] 5007 dereddened flux.
 }

}



##############################################################################
##Save data to files
print('Saving fluxes to data file.')
resume <- data.frame(id,AHa,Cbeta,Fa,Fb,Fn1,Fn,Fo1,Fo)
tabla <- str_c(galaxy[i],"/Cardelli_Base_Fluxes.dat")
write.table(resume, tabla, sep="\t",quote=FALSE)

}
