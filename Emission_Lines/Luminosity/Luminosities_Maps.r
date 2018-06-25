#############################################################################
#############################################################################

##   This program plots maps of luminosities for the H-alpha, [O III] 5007,##
##   and [N II] 6583 emission lines.                                       ##

##   February 25, 2018 A. Robleto-Orús                                     ##

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

setwd("/data2/CALIFA/") #Directrory with our data

data <- read.table("Seyfert_Centro1.dat", header=TRUE)
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
path_ha <- str_c(galaxy[i],".COMB/lis_ha.res") #creates path to Ha and [N II]data file for each galaxy
data0 <- read.table(path_ha, header=TRUE)

print('Extracting H-beta and [O III] 5007 data.')
path_hb <- str_c(galaxy[i],".COMB/lis_hb.res") #creates path to H-beta and [OIII] data file for each galaxy
data3 <- read.table(path_hb, header=TRUE)
         
##Merge data

DATA <- merge(data0,data3,by.x = 1, by.y =1)
attach(DATA)


##Extracting coordinates from ID
## Spaxel's coordinates are in the id in format YYXX, with Y = DEC and X = RA.

print('Extracting spaxel coordinates.')
y <-     (trunc(id[SN_starlight.x > 5 & snha >= 3 & sigmaha < 400 & fluxha>0.1 & fluxha<1000 & !is.na(fluxha) & snnii2 > 2 & sigmanii2 < 400 & fluxnii2>0.1 & fluxnii2<1000 & !is.na(fluxnii2)& snhb > 2 & sigmahb < 400 & fluxhb>0.1 & fluxhb<1000 & !is.na(fluxhb) & snoiii2 > 2 & sigmaoiii2 < 400 & fluxoiii2>0.1 & fluxoiii2<1000 & !is.na(fluxoiii2)]/100)  ) #Obtaining Y from id.
x <-     (id[SN_starlight.x > 5 & snha >= 3 & sigmaha < 400 & fluxha>0.1 & fluxha<1000 & !is.na(fluxha) & snnii2 > 2 & sigmanii2 < 400 & fluxnii2>0.1 & fluxnii2<1000 & !is.na(fluxnii2)& snhb > 2 & sigmahb < 400 & fluxhb>0.1 & fluxhb<1000 & !is.na(fluxhb) & snoiii2 > 2 & sigmaoiii2 < 400 & fluxoiii2>0.1 & fluxoiii2<1000 & !is.na(fluxoiii2)]-(y*100) ) #Obtaining X from id.

##Extracting line surface specific intensities [1e-16 erg cm^-2 s^-1 A^-1 arcsec^-2]

print('Extracting emission lines.')

#H-alpha 6563
print('Extracting H-alpha.')
fa <-    fluxha[SN_starlight.x > 5 & snha >= 3 & sigmaha < 400 & fluxha>0.1 & fluxha<1000 & !is.na(fluxha) & snnii2 > 2 & sigmanii2 < 400 & fluxnii2>0.1 & fluxnii2<1000 & !is.na(fluxnii2)& snhb > 2 & sigmahb < 400 & fluxhb>0.1 & fluxhb<1000 & !is.na(fluxhb) & snoiii2 > 2 & sigmaoiii2 < 400 & fluxoiii2>0.1 & fluxoiii2<1000 & !is.na(fluxoiii2) ] #Obtaining the H-alpha flux density (10^-16 erg/(s cm^2 A pix))

#[N II] 6583
print('Extracting [N II] 6583.')
fn <-    fluxnii2[SN_starlight.x > 5 & snha >= 3 & sigmaha < 400 & fluxha>0.1 & fluxha<1000 & !is.na(fluxha) & snnii2 > 2 & sigmanii2 < 400 & fluxnii2>0.1 & fluxnii2<1000 & !is.na(fluxnii2)& snhb > 2 & sigmahb < 400 & fluxhb>0.1 & fluxhb<1000 & !is.na(fluxhb) & snoiii2 > 2 & sigmaoiii2 < 400 & fluxoiii2>0.1 & fluxoiii2<1000 & !is.na(fluxoiii2) ] #Obtaining the[N II] flux density (10^-16 erg/(s cm^2 A pix))

#H-beta 4861
print('Extracting H-beta.')
fb <-    fluxhb[SN_starlight.x > 5 & snha >= 3 & sigmaha < 400 & fluxha>0.1 & fluxha<1000 & !is.na(fluxha) & snnii2 > 2 & sigmanii2 < 400 & fluxnii2>0.1 & fluxnii2<1000 & !is.na(fluxnii2) & snhb > 2 & sigmahb < 400 & fluxhb>0.1 & fluxhb<1000 & !is.na(fluxhb) & snoiii2 > 2 & sigmaoiii2 < 400 & fluxoiii2>0.1 & fluxoiii2<1000 & !is.na(fluxoiii2)] #Obtaining the[H-beta] flux density (10^-16 erg/(s cm^2 A pix))

#[O III] 5007
print('Extracting [O III]5007.')
fo <-    fluxoiii2[SN_starlight.x > 5 & snha >= 3 & sigmaha < 400 & fluxha>0.1 & fluxha<1000 & !is.na(fluxha) & snnii2 > 2 & sigmanii2 < 400 & fluxnii2>0.1 & fluxnii2<1000 & !is.na(fluxnii2) & snhb > 2 & sigmahb < 400 & fluxhb>0.1 & fluxhb<1000 & !is.na(fluxhb) & snoiii2 > 2 & sigmaoiii2 < 400 & fluxoiii2>0.1 & fluxoiii2<1000 & !is.na(fluxoiii2)] #Obtaining the[H-beta] flux density (10^-16 erg/(s cm^2 A pix))

#Equivalent width of H-alpha (Angstrom)
print('Extracting H-alpha equivalent width.')
EWa1 <-   eqwha[SN_starlight.x > 5 & snha >= 3 & sigmaha < 400 & fluxha>0.1 & fluxha<1000 & !is.na(fluxha) & snnii2 > 2 & sigmanii2 < 400 & fluxnii2>0.1 & fluxnii2<1000 & !is.na(fluxnii2) & snhb > 2 & sigmahb < 400 & fluxhb>0.1 & fluxhb<1000 & !is.na(fluxhb) & snoiii2 > 2 & sigmaoiii2 < 400 & fluxoiii2>0.1 & fluxoiii2<1000 & !is.na(fluxoiii2)] #Obtaining the[H-beta] flux density (10^-16 erg/(s cm^2 A pix))


##id of central spaxel
print('Determining central spaxel postion.')
xc <- trunc(naxis1[i]/2)
yc <- trunc(naxis2[i]/2)

X3 <- x - xc
Y3 <- y - yc


##Convert from surface specific intensity to standard specific intensity [erg cm^-2 s^-1 A^-1 sr^-1]

print('Converting surface specific units to specific intensities.')
fa   <-(1e-16)*fa  /2.3504e-11 
fn   <-(1e-16)*fn  /2.3504e-11 
fb   <-(1e-16)*fb  /2.3504e-11 
fo   <-(1e-16)*fo  /2.3504e-11 

#####################################################

##                    DEREDDENING                 ##  

#####################################################

print('First dereddening.')

Rv <- 3.1
l <-  c(6563,6548,6584,4861,4959,5007,3729) #Lines we are going to unredden, in Angstrom.
l <- l*1e-4 #Convert Angstrom to micrometers
xl <- 1/l
yl <- xl - 1.82

a_x <- 1 + (0.17699*yl) - (0.50447*yl^2) - (0.02427*yl^3) + (0.72085*yl^4) + (0.01979*yl^5) - (0.77530*yl^6) + (0.32999*yl^7)
b_x <- (1.41338*yl) + (2.28305*yl^2) + (1.07233*yl^3) -(5.38434*yl^4) - (0.62251*yl^5) + (5.30260*yl^6) - (2.09002*yl^7)


Cbeta <- seq(1:length(y))
rIb <- seq(1:length(y))
rIa <- seq(1:length(y))
rIo <- seq(1:length(y))
rIn <- seq(1:length(y))

##Claculate the extinction law

for(j in 1:length(y)){
  
  f_l <- seq(1:length(l))  #CCM89 extinction law, normalized to V-band
  fb_l <- seq(1:length(l)) #Extinction law normalized to H-beta
  
  for(k in 1:length(l)){
    f_l[k] <- (a_x[k] + (b_x[k] / Rv))
  }
  fb_l <- f_l/f_l[4]
  
  ##Find constant c(H-beta)  
  Cbeta[j] <- (log10(2.86) - log10(fa[j]/fb[j]))/(fb_l[1] - 1)
  
  ##Deredden each line normalized to I_H-beta
  
  rIb[j] <- 100*((fb[j])/fb[j])*10^(Cbeta[j]*(fb_l[4]-1))
  
  rIa[j] <- 100*((fa[j])/fb[j])*10^(Cbeta[j]*(fb_l[1]-1))
  
  rIn[j]  <- 100*((fn[j])/fb[j])*10^(Cbeta[j]*(fb_l[3]-1))
  
  rIo[j] <- 100*((fo[j])/fb[j])*10^(Cbeta[j]*(fb_l[6]-1))
}

########################################################################

##			LINE RATIOS			              ##

########################################################################

print('Calculating BPT-NII line ratios.')

#[N II]6584 / H-alpha 
na1 <- log10(rIn/rIa)
#[O III] 5007 / H-beta
ob1 <- log10(rIo/rIb)

na  <- na1[which(!is.na(na1) & !is.na(ob1) & !is.na(EWa1))]
ob  <- ob1[which(!is.na(na1) & !is.na(ob1) & !is.na(EWa1))]
EWa <- EWa1[which(!is.na(na1) & !is.na(ob1) & !is.na(EWa1))]

X2 <- x[which(!is.na(na1) & !is.na(ob1) & !is.na(EWa1))]
Y2 <- y[which(!is.na(na1) & !is.na(ob1) & !is.na(EWa1))]


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


for(kk in 1:length(X2)){
  if(type[kk] == " red"){
    
    ##Find constant c(H-beta)  
    Cbeta[j] <- (log10(3.1) - log10(fa[kk]/fb[kk]))/(fb_l[1] - 1)
    
    ##Deredden each line normalized to I_H-beta
    #Multiply each by 100 if they are used directly.
    
    rIb[kk] <- ((fb[kk])/fb[kk])*10^(Cbeta[j]*(fb_l[4]-1))
    
    rIa[kk] <- ((fa[kk])/fb[kk])*10^(Cbeta[j]*(fb_l[1]-1))
    
    rIn[kk]  <- ((fn[kk])/fb[kk])*10^(Cbeta[j]*(fb_l[3]-1))
    
    rIo[kk] <- ((fo[kk])/fb[kk])*10^(Cbeta[j]*(fb_l[6]-1))
  }
}  

#Determination of the absolute fluxes

AHa <- seq(1,length(X2),1)
Fa  <- seq(1,length(X2),1)
Fb  <- seq(1,length(X2),1)
Fn  <- seq(1,length(X2),1)
Fo  <- seq(1,length(X2),1)

for (ww in 1:length(X2)){
  if(type[ww] == 'red'){
   AHa[ww] <- (f_l[1]/(-0.4*(f_l[1]-f_l[4])))*log10((f_l[1]/f_l[4])/3.1) #H-alpha extinction (magnitudes, V-band normalized).
   Fa[ww] <- fa[ww]*10^(0.4*AHa[ww]) #H-alpha dereddened absolute flux
   Fb[ww] <- Fa[ww]/3.1 #H-beta dereddened flux
   Fn[ww] <- rIn[ww]*Fb[ww] #[N II] 6583 dereddenend flux
   Fo[ww] <- rIo[ww]*Fb[ww] #[O III] 5007 dereddened flux
  
  }
  else {
   AHa[ww] <- (f_l[1]/(-0.4*(f_l[1]-f_l[4])))*log10((f_l[1]/f_l[4])/2.86) #H-alpha extinction (magnitudes, V-band normalized).
   Fa[ww] <- fa[ww]*10^(0.4*AHa[ww]) #H-alpha dereddened absolute flux
   Fb[ww] <- Fa[ww]/2.86 #H-beta dereddened flux
   Fn[ww] <- rIn[ww]*Fb[ww] #[N II] 6583 dereddenend flux
   Fo[ww] <- rIo[ww]*Fb[ww] #[O III] 5007 dereddened flux  
  }
  
}

########################################################################

##				Determination of Luminosities 		      ##

########################################################################

##Calculate distance

H0 <- 70 #Hubble constant (km/s). Our group = 70, Riess(2016) = 73.0, Planck(2015) = 67.8
print(paste('Calculating cosmolgical distance. H0 = ',H0,' [km/s/Mpc].'))
c <- 299792.458 #Speed of light (km/s)
D <- ((c*z[i])/H0)*3.0857e24 #Distance (cm)

## H-alpha luminosity

print('Determining H-alpha luminosity.')
La <- 4*pi*D^2*Fa*2.3504e-11 #Luminosity [erg s^-1 ].


## [O III] 5007 luminosity

print('Determining [O III] 5007 luminosity.')
Lo <- 4*pi*D^2*Fo*2.3504e-11 #Luminosity [erg s^-1 ].


## [N II] 6584

print('Determining [N II] 6583 luminosity.')
Ln <- 4*pi*D^2*Fn*2.3504e-11 #Luminosity [erg s^-1 ].

##########################################################################

##				PLOTS					##

##########################################################################

print("Plotting maps")

par(mar=c(5,5,4,4))
par(bg='white')

####################################
#Set background image. Only once for all maps.

ima <- readPNG(str_c(galaxy[i],".COMB/",galaxy[i],".png")) #Background image file in png.
lim <- par()




##########################################################################
##H-alpha Map

print('Plotting H-alpha luminosity map.')

## Specify color map for z.

rbPal <-  colorRampPalette((rainbow(100,start=0,end=0.9))) 
Col <- rbPal(100)[as.numeric(cut(log10(La),breaks = 100))]

plot(X3,Y3,xlim=c(-40,40),ylim=c(-40,40),xaxs= "i",yaxs="i", col="white", main = galaxy[i],xlab = expression(paste(Delta, alpha," (arcsec)")), ylab = expression(paste(Delta, delta," (arcsec)")), cex.lab=1.3, cex.axis=1.3)
rasterImage(ima,-40,-40,40,40,interpolate=F)
points(X3,Y3, col = (Col), pch=15,cex = 0.7) #Plots spaxels with color map reversed (to have greater values darker and lower lighter).
grid() #Coordinates'grid.
abline(h = 0, v = 0, col = "gray60",lwd=2)
#mtext(side=3,expression(paste(" log"[10]," H",alpha," Luminosity [erg s"^-1,"]")))
image.plot( -40, -40, legend.only=TRUE, zlim= range(log10(La)), col = ((rainbow(100,start=0,end=0.9)))) #Plots the color bar. Reverse colors here too.

#Saving plot.

map <- str_c("/data2/CALIFA/",galaxy[i],".COMB/",galaxy[i],"_Luminosity_Ha.eps")
dev.copy2eps(file=map)
map <- str_c("convert -density 300 /data2/CALIFA/",galaxy[i],".COMB/",galaxy[i],"_Luminosity_Ha.eps /data2/CALIFA/",galaxy[i],".COMB/",galaxy[i],"_Luminosity_Ha.png")
system(map)


##########################################################################
##[O III 5007] Map

print('Plotting [O III] 5007 luminosity map.')

## Specify color map for z.

rbPal <-  colorRampPalette((rainbow(100,start=0,end=0.9))) 
Col <- rbPal(100)[as.numeric(cut(log10(Lo),breaks = 100))]


plot(X3,Y3,xlim=c(-40,40),ylim=c(-40,40),xaxs= "i",yaxs="i", col="white", main = galaxy[i],xlab = expression(paste(Delta, alpha," (arcsec)")), ylab = expression(paste(Delta, delta," (arcsec)")), cex.lab=1.3, cex.axis=1.3)
rasterImage(ima,-40,-40,40,40,interpolate=F)
points(X3,Y3, col = (Col), pch=15,cex = 0.7) #Plots spaxels with color map reversed (to have greater values darker and lower lighter).
grid() #Coordinates'grid.
abline(h = 0, v = 0, col = "gray60",lwd=2)
#mtext(side=3, expression(paste(" log"[10]," [O III] ",lambda,"5007 Luminosity [erg s"^-1,"]")))
image.plot( -40, -40, legend.only=TRUE, zlim= range(log10(Lo)), col = ((rainbow(100,start=0,end=0.9)))) #Plots the color bar. Reverse colors here too.

#Saving plot.

map <- str_c("/data2/CALIFA/",galaxy[i],".COMB/",galaxy[i],"_Luminosity_OIII.eps")
dev.copy2eps(file=map)
map <- str_c("convert -density 300 /data2/CALIFA/",galaxy[i],".COMB/",galaxy[i],"_Luminosity_OIII.eps /data2/CALIFA/",galaxy[i],".COMB/",galaxy[i],"_Luminosity_OIII.png")
system(map)

##########################################################################
##[N II] 6584  Map

print('Plotting [N II] 6583 luminosity map.')

## Specify color map for z.

rbPal <-  colorRampPalette((rainbow(100,start=0,end=0.9))) 
Col <- rbPal(100)[as.numeric(cut(log10(Ln),breaks = 100))]

plot(X3,Y3,xlim=c(-40,40),ylim=c(-40,40),xaxs= "i",yaxs="i", col="white", main=galaxy[i], xlab = expression(paste(Delta, alpha," (arcsec)")), ylab = expression(paste(Delta, delta," (arcsec)")), cex.lab=1.3, cex.axis=1.3)
rasterImage(ima,-40,-40,40,40,interpolate=F)
points(X3,Y3, col = (Col), pch=15,cex = 0.7) #Plots spaxels with color map reversed (to have greater values darker and lower lighter).
grid() #Coordinates'grid.
abline(h = 0, v = 0, col = "gray60",lwd=2)
#mtext(side=3,expression(paste(" log"[10]," [N II] ",lambda," 6584 Luminosity [erg s"^-1,"]")))
image.plot( -40, -40, legend.only=TRUE, zlim= range(log10(Ln)), col = ((rainbow(100,start=0,end=0.9)))) #Plots the color bar. Reverse colors here too.

#Saving plot.

map <- str_c("/data2/CALIFA/",galaxy[i],".COMB/",galaxy[i],"_Luminosity_NII.eps")
dev.copy2eps(file=map)
map <- str_c("convert -density 300 /data2/CALIFA/",galaxy[i],".COMB/",galaxy[i],"_Luminosity_NII.eps /data2/CALIFA/",galaxy[i],".COMB/",galaxy[i],"_Luminosity_NII.png")
system(map)

##########################################################################
##H-alpha equivalent width

print('Plotting H-alpha equivalent width map.')

## Specify color map for z.

rbPal <-  colorRampPalette((rainbow(100,start=0,end=0.9))) 
Col <- rbPal(100)[as.numeric(cut(log10(EWa),breaks = 100))]

plot(X3,Y3,xlim=c(-40,40),ylim=c(-40,40),xaxs= "i",yaxs="i", col="white", main=galaxy[i], xlab = expression(paste(Delta, alpha," (arcsec)")), ylab = expression(paste(Delta, delta," (arcsec)")), cex.lab=1.3, cex.axis=1.3)
rasterImage(ima,-40,-40,40,40,interpolate=F)
points(X3,Y3, col = (Col), pch=15,cex = 0.7) #Plots spaxels with color map reversed (to have greater values darker and lower lighter).
grid() #Coordinates'grid.
abline(h = 0, v = 0, col = "gray60",lwd=2)
#mtext(side=3,expression(paste(" log"[10]," [N II] ",lambda," 6584 Luminosity [erg s"^-1,"]")))
image.plot( -40, -40, legend.only=TRUE, zlim= range(log10(EWa)), col = ((rainbow(100,start=0,end=0.9)))) #Plots the color bar. Reverse colors here too.

#Saving plot.

map <- str_c("/data2/CALIFA/",galaxy[i],".COMB/",galaxy[i],"_WHa.eps")
dev.copy2eps(file=map)
map <- str_c("convert -density 300 /data2/CALIFA/",galaxy[i],".COMB/",galaxy[i],"_WHa.eps /data2/CALIFA/",galaxy[i],".COMB/",galaxy[i],"_WHa.png")
system(map)

copy1 <- str_c("cp ",galaxy[i],".COMB/",galaxy[i],"_Luminosity*.png /home/aitor/Doctorado/Plots/Luminosities/")
system(copy1)

copy2 <- str_c("cp ",galaxy[i],".COMB/",galaxy[i],"_WHa.png /home/aitor/Doctorado/Plots/Equivalent_Width/")
system(copy2)

}