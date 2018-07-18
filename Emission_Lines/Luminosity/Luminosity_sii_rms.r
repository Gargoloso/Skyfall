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
library("astro")

########################################################################

##				DATA INPUT			      ##

########################################################################

setwd("~/Rings/ringed_work") #Directrory with our data

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

print('Extracting base data.')  
path_base <- str_c(galaxy[i],"/Calzetti_sii_Fluxes.dat") #creates path to Ha and [N II]data file for each galaxy
data0 <- read.table(path_base, header=TRUE)
attach(data0)

##id of central spaxel
xc <- trunc(naxis1[i]/2)
yc <- trunc(naxis2[i]/2)


####################################################################
####################################################################

##		  EXTRACTION OF EMISSION-LINE DATA                ##

####################################################################
####################################################################

print('EXTRACTING DEREDDENED LINE SURFACE SPECIFIC FLUXES')

## Extracting coordinates from ID
## Spaxel's coordinates are in the id in format YYXX, with Y = DEC and X = RA.

y <- trunc(ID[which(!is.na(Fsii1) & !is.na(rms_Fsii1) & !is.na(Fsii2) & !is.na(rms_Fsii2) & (rms_Fsii1/Fsii1) <= 1 & (rms_Fsii2/Fsii2) <=1)]/100) # y coordinate.
x <- ID[which(!is.na(Fsii1) & !is.na(rms_Fsii1) & !is.na(Fsii2) & !is.na(rms_Fsii2) & (rms_Fsii1/Fsii1) <= 1 & (rms_Fsii2/Fsii2) <=1)]-(y*100) # x coordinate.

ID2 <- ID[which(!is.na(Fsii1) & !is.na(rms_Fsii1) & !is.na(Fsii2) & !is.na(rms_Fsii2) & (rms_Fsii1/Fsii1) <= 1 & (rms_Fsii2/Fsii2) <=1)] # ID of each spectra.

## Extracting dereddenend surface specific fluxes and rms [erg s^-1 cm^-2 arcsec^-2].
print('Extracting and converting surface specific units to specific intensities.')
# [S II] 6716
FSii1 <- Fsii1[which(!is.na(Fsii1) & !is.na(rms_Fsii1) & !is.na(Fsii2) & !is.na(rms_Fsii2) & (rms_Fsii1/Fsii1) <= 1 & (rms_Fsii2/Fsii2) <=1)]/2.3504e-11
rms_Fii1 <- rms_Fsii1[which(!is.na(Fsii1) & !is.na(rms_Fsii1) & !is.na(Fsii2) & !is.na(rms_Fsii2) & (rms_Fsii1/Fsii1) <= 1 & (rms_Fsii2/Fsii2) <=1)] 

# [S II] 6731
FSii2 <- Fsii2[which(!is.na(Fsii1) & !is.na(rms_Fsii1) & !is.na(Fsii2) & !is.na(rms_Fsii2) & (rms_Fsii1/Fsii1) <= 1 & (rms_Fsii2/Fsii2) <=1)]/2.3504e-11
rms_Fii2 <- rms_Fsii2[which(!is.na(Fsii1) & !is.na(rms_Fsii1) & !is.na(Fsii2) & !is.na(rms_Fsii2) & (rms_Fsii1/Fsii1) <= 1 & (rms_Fsii2/Fsii2) <=1)] 


# Translation of coordinates setting origin to the field centre.
Y <- y - yc
X <- x - xc



########################################################################

##				Determination of Luminosities 		      ##

########################################################################
##Calculate distance


H0 <- 73 #Hubble constant (km/s). Our group = 70, Riess(2016) = 73.0, Planck(2015) = 67.8
#rms_H0 <-
print(paste('Calculating cosmolgical distance. H0 = ',H0,' [km/s/Mpc].'))
c <- 299792458 #Speed of light (km/s)
#D <- ((c*z[i])/H0)*3.0857e24 #Distance (cm)
#rms_D <- (sqrt(((c/H0)*rms_z)^2+((-c*z/H0^2)*rms_H0)^2))**3.0857e24 
d <- lumdist(z[i],c,H0) #Distance (Mpc)
D <- d*3.0857e24 #Distance (cm)
sp <- 4.8414e-6 * 1e6 * d #size of 1 spaxel's side (pc)
asp <- sp^2  #area of 1 spaxel (pc^2)

## H-alpha luminosity

print('Determining [S II] 6716 luminosity.')
LFSii1 <- 4*pi*D^2*FSii1*2.3504e-11 #Luminosity [erg s^-1 ].

## H-beta Luminosity

print('Determining [S II] 6731 luminosity.')
LFSii2 <- 4*pi*D^2*FSii2*2.3504e-11 #Luminosity [erg s^-1 ].


#########################################       MAPS ONLY!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#removing inferior quartile.
remove_outliers <- function(x, na.rm=TRUE, ...){
  qnt <- quantile(x, probs=c(.10,.99, na.rm = na.rm, ...))
  H <- 1.5*IQR(x, na.rm = na.rm)
  clean <- x
  clean[x < (qnt[1] - H)] <- NA
  clean[x > (qnt[2] + H)] <- NA
  clean
}

LFSII1 <- 10^(remove_outliers(log10(LFSii1)))
Lsii1  <- LFSII1[which(!is.na(LFSII1))]/(asp * 3.826e33) #Convert to Solar luminosity / pc^2.
Xsii1 <-  X[which(!is.na(LFSII1))]
Ysii1  <- Y[which(!is.na(LFSII1))]
IDsii1 <-  sprintf("%04d", as.numeric(ID2[which(!is.na(LFSII1))]))


resume <- data.frame(IDsii1,Lsii1)
tabla <- str_c(galaxy[i],"/",galaxy[i],"_Lum_sii6716.dat")
write.table(resume, tabla, sep="\t",quote=FALSE, row.names = FALSE)

LFSII2 <- 10^(remove_outliers(log10(LFSii2)))
Lsii2  <- LFSII2[which(!is.na(LFSII2))]/(asp * 3.826e33) #Convert to Solar luminosity / pc^2.
Xsii2  <- X[which(!is.na(LFSII2))]
Ysii2  <- Y[which(!is.na(LFSII2))]
IDsii2 <-  sprintf("%04d", as.numeric(ID2[which(!is.na(LFSII2))]))


resume <- data.frame(IDsii2,Lsii2)
tabla <- str_c(galaxy[i],"/",galaxy[i],"_Lum_sii6731.dat")
write.table(resume, tabla, sep="\t",quote=FALSE, row.names = FALSE)

##########################################################################

##				PLOTS					##

##########################################################################

print("Plotting maps")

par(mar=c(5,5,4.75,4.75))
par(bg='white')

####################################
#Set background image. Only once for all maps.

#ima <- readPNG(str_c(galaxy[i],".COMB/",galaxy[i],".png")) #Background image file in png.
#lim <- par()

##########################################################################
##H-alpha Map

print('Plotting [H-alpha [SII] 6716 luminosity map.')

## Specify color map for z.

rbPal <-  colorRampPalette(rev(rainbow(100,start=0,end=0.7))) 
Col <- rbPal(100)[as.numeric(cut(log10(Lsii1),breaks = 100))]

title <- str_c(galaxy[i]," [S II] 6716 luminosity")
plot(Xsii1,Ysii1,xlim=c(-40,40),ylim=c(-40,40),xaxs= "i",yaxs="i", col="white", main = title,xlab = expression(paste(Delta, alpha," (arcsec)")), ylab = expression(paste(Delta, delta," (arcsec)")), cex.lab=1.3, cex.axis=1.3)
#rasterImage(ima,-40,-40,40,40,interpolate=F)
points(Xsii1,Ysii1, col = (Col), pch=15,cex = 0.7) #Plots spaxels with color map reversed (to have greater values darker and lower lighter).
grid() #Coordinates'grid.
abline(h = 0, v = 0, col = "gray60",lwd=2)
#mtext(side=3,expression(paste(" log"[10]," H",alpha," Luminosity [erg s"^-1,"]")))
image.plot( -40, -40, legend.only=TRUE, zlim= range(log10(Lsii1)), col = (rev(rainbow(100,start=0,end=0.7)))) #Plots the color bar. Reverse colors here too.

#Saving plot.

map <- str_c(galaxy[i],"/",galaxy[i],"_Luminosity_sii6716.eps")
dev.copy2eps(file=map)
map <- str_c("convert -density 300 ",galaxy[i],"/",galaxy[i],"_Luminosity_sii6716.eps ",galaxy[i],"/",galaxy[i],"_Luminosity_sii6716.png")
system(map)


##########################################################################
##[O III 5007] Map

print('Plotting [S II] 6731 luminosity map.')

## Specify color map for z.

rbPal <-  colorRampPalette(rev(rainbow(100,start=0,end=0.7))) 
Col <- rbPal(100)[as.numeric(cut(log10(Lsii2),breaks = 100))]

title <- str_c(galaxy[i]," [S II] 6731 luminosity")
plot(Xsii2,Ysii2,xlim=c(-40,40),ylim=c(-40,40),xaxs= "i",yaxs="i", col="white", main = title,xlab = expression(paste(Delta, alpha," (arcsec)")), ylab = expression(paste(Delta, delta," (arcsec)")), cex.lab=1.3, cex.axis=1.3)
#rasterImage(ima,-40,-40,40,40,interpolate=F)
points(Xsii2,Ysii2, col = (Col), pch=15,cex = 0.7) #Plots spaxels with color map reversed (to have greater values darker and lower lighter).
grid() #Coordinates'grid.
abline(h = 0, v = 0, col = "gray60",lwd=2)
#mtext(side=3, expression(paste(" log"[10]," [O III] ",lambda,"5007 Luminosity [erg s"^-1,"]")))
image.plot( -40, -40, legend.only=TRUE, zlim= range(log10(Lsii2)), col = (rev(rainbow(100,start=0,end=0.7)))) #Plots the color bar. Reverse colors here too.

#Saving plot.

map <- str_c(galaxy[i],"/",galaxy[i],"_Luminosity_sii6731.eps")
dev.copy2eps(file=map)
map <- str_c("convert -density 300 ",galaxy[i],"/",galaxy[i],"_Luminosity_sii6731.eps ",galaxy[i],"/",galaxy[i],"_Luminosity_sii6731.png")
system(map)



}