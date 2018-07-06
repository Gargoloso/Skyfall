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
path_base <- str_c(galaxy[i],"/Calzetti_Base_Fluxes.dat") #creates path to Ha and [N II]data file for each galaxy
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

y <- trunc(ID[which(!is.na(AHa) & !is.na(rms_AHa) & AHa <= 4)]/100) # y coordinate.
x <- ID[which(!is.na(AHa) & !is.na(rms_AHa) & AHa <= 4)]-(y*100) # x coordinate.

ID2 <- ID[which(!is.na(AHa) & !is.na(rms_AHa) & AHa <= 4)] # ID of each spectra.

## Extracting dereddenend surface specific fluxes and rms [erg s^-1 cm^-2 arcsec^-2].

# H-alpha.
FA <- Fa[which(!is.na(AHa) & !is.na(rms_AHa) & AHa <= 4)] 
rms_FA <- rms_Fa[which(!is.na(AHa) & !is.na(rms_AHa) & AHa <= 4)]

# H-beta.
FB <- Fb[which(!is.na(AHa) & !is.na(rms_AHa) & AHa <= 4)] #H-beta flux.
rms_FB <- rms_Fb[which(!is.na(AHa) & !is.na(rms_AHa) & AHa <= 4)]

# [N II] 6548
FN1 <- Fn1[which(!is.na(AHa) & !is.na(rms_AHa) & AHa <= 4)]
rms_FN1 <- rms_Fn1[which(!is.na(AHa) & !is.na(rms_AHa) & AHa <= 4)]

# [N II] 6583
FN <- Fn[which(!is.na(AHa) & !is.na(rms_AHa) & AHa <= 4)]
rms_FN <- rms_Fn[which(!is.na(AHa) & !is.na(rms_AHa) & AHa <= 4)]

# [O III] 4959
FO1 <- Fo1[which(!is.na(AHa) & !is.na(rms_AHa) & AHa <= 4)]
rms_FO1 <- rms_Fo1[which(!is.na(AHa) & !is.na(rms_AHa) & AHa <= 4)]   

# [O III] 5007
FO <- Fo[which(!is.na(AHa) & !is.na(rms_AHa) & AHa <= 4)]
rms_FO <- rms_Fo[which(!is.na(AHa) & !is.na(rms_AHa) & AHa <= 4)]               

# Translation of coordinates setting origin to the field centre.
Y <- y - yc
X <- x - xc

##Convert from surface specific flux to standard specific intensity [erg cm^-2 s^-1 A^-1 sr^-1]

print('Converting surface specific units to specific intensities.')
fa   <- FA[which(rms_FA/FA <= 1.0)]/2.3504e-11 
fn1  <- FN1[which(rms_FN1/FN1 <= 1.0)]/2.3504e-11 
fn   <- FN[which(rms_FN/FN <= 1.0)]/2.3504e-11 
fb   <- FB[which(rms_FB/FB <= 1.0)]/2.3504e-11
fo1  <- FO1[which(rms_FO1/FO1 <= 1.0)]/2.3504e-11
fo   <- FO[which(rms_FO/FO <= 1.0)]/2.3504e-11 

XA  <- X[which(rms_FA/FA <= 1.0)]
YA  <- Y[which(rms_FA/FA <= 1.0)]
IDA <- ID2[which(rms_FA/FA <= 1.0)]

XN1  <- X[which(rms_FN1/FN1 <= 1.0)]
YN1  <- Y[which(rms_FN1/FN1 <= 1.0)]
IDN1 <- ID2[which(rms_FN1/FN1 <= 1.0)]

XN  <- X[which(rms_FN/FN <= 1.0)]
YN  <- Y[which(rms_FN/FN <= 1.0)]
IDN <- ID2[which(rms_FN/FN <= 1.0)]

XB  <- X[which(rms_FB/FB <= 1.0)]
YB  <- Y[which(rms_FB/FB <= 1.0)]
IDB <- ID2[which(rms_FB/FB <= 1.0)]

XO1  <- X[which(rms_FO1/FO1 <= 1.0)]
YO1  <- Y[which(rms_FO1/FO1 <= 1.0)]
IDO1 <- ID2[which(rms_FO1/FO1 <= 1.0)]

XO  <- X[which(rms_FO/FO <= 1.0)]
YO  <- Y[which(rms_FO/FO <= 1.0)]
IDO <- ID2[which(rms_FO/FO <= 1.0)]

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

print('Determining H-alpha luminosity.')
LA1 <- 4*pi*D^2*fa*2.3504e-11 #Luminosity [erg s^-1 ].

## H-beta Luminosity

print('Determining H-beta luminosity.')
LB1 <- 4*pi*D^2*fb*2.3504e-11 #Luminosity [erg s^-1 ].

## [O III] 4959 luminosity

print('Determining [O III] 4959 luminosity.')
LO11 <- 4*pi*D^2*fo1*2.3504e-11 #Luminosity [erg s^-1 ].

## [O III] 5007 luminosity

print('Determining [O III] 5007 luminosity.')
LO1 <- 4*pi*D^2*fo*2.3504e-11 #Luminosity [erg s^-1 ].

## [N II] 6548

print('Determining [N II] 6548 luminosity.')
LN11 <- 4*pi*D^2*fn1*2.3504e-11 #Luminosity [erg s^-1 ].

## [N II] 6583

print('Determining [N II] 6583 luminosity.')
LN1 <- 4*pi*D^2*fn*2.3504e-11 #Luminosity [erg s^-1 ].

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

LA2 <- 10^(remove_outliers(log10(LA1)))
La  <- LA2[which(!is.na(LA2))]/(asp * 3.826e33) #Convert to Solar luminosity / pc^2.
Xa  <- XA[which(!is.na(LA2))]
Ya  <- YA[which(!is.na(LA2))]
IDa <-  sprintf("%04d", as.numeric(IDA[which(!is.na(LA2))]))

resume <- data.frame(IDa,La)
tabla <- str_c(galaxy[i],"/",galaxy[i],"_Lum_Ha.dat")
write.table(resume, tabla, sep="\t",quote=FALSE, row.names = FALSE)



LN112 <- 10^(remove_outliers(log10(LN11)))
Ln1  <- LN112[which(!is.na(LN112))]/(asp * 3.826e33) #Convert to Solar luminosity.
Xn1  <- XN1[which(!is.na(LN112))]
Yn1  <- YN1[which(!is.na(LN112))]
IDn1 <- IDN1[which(!is.na(LN112))]

LN12 <- 10^(remove_outliers(log10(LN1)))
Ln  <- LN12[which(!is.na(LN12))]/(asp * 3.826e33) #Convert to Solar luminosity.
Xn  <- XN[which(!is.na(LN12))]
Yn  <- YN[which(!is.na(LN12))]
IDn <- IDN[which(!is.na(LN12))]

LB2 <- 10^(remove_outliers(log10(LB1)))
Lb  <- LB2[which(!is.na(LB2))]/(asp * 3.826e33) #Convert to Solar luminosity.
Xb  <- XB[which(!is.na(LB2))]
Yb  <- YB[which(!is.na(LB2))]
IDb <- IDB[which(!is.na(LB2))]

LO112 <- 10^(remove_outliers(log10(LO11)))
Lo1  <- LO112[which(!is.na(LO112))]/(asp * 3.826e33) #Convert to Solar luminosity.
Xo1  <- XO1[which(!is.na(LO112))]
Yo1  <- YO1[which(!is.na(LO112))]
IDo1 <- IDO1[which(!is.na(LO112))]

LO2 <- 10^(remove_outliers(log10(LO1)))
Lo  <- LO2[which(!is.na(LO2))]/(asp * 3.826e33) #Convert to Solar luminosity.
Xo  <- XO[which(!is.na(LO2))]
Yo  <- YO[which(!is.na(LO2))]
IDo <- IDO[which(!is.na(LO2))]


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

print('Plotting H-alpha luminosity map.')

## Specify color map for z.

rbPal <-  colorRampPalette(rev(rainbow(100,start=0,end=0.7))) 
Col <- rbPal(100)[as.numeric(cut(log10(La),breaks = 100))]

title <- str_c(galaxy[i]," H-alpha luminosity")
plot(Xa,Ya,xlim=c(-40,40),ylim=c(-40,40),xaxs= "i",yaxs="i", col="white", main = title,xlab = expression(paste(Delta, alpha," (arcsec)")), ylab = expression(paste(Delta, delta," (arcsec)")), cex.lab=1.3, cex.axis=1.3)
#rasterImage(ima,-40,-40,40,40,interpolate=F)
points(Xa,Ya, col = (Col), pch=15,cex = 0.7) #Plots spaxels with color map reversed (to have greater values darker and lower lighter).
grid() #Coordinates'grid.
abline(h = 0, v = 0, col = "gray60",lwd=2)
#mtext(side=3,expression(paste(" log"[10]," H",alpha," Luminosity [erg s"^-1,"]")))
image.plot( -40, -40, legend.only=TRUE, zlim= range(log10(La)), col = (rev(rainbow(100,start=0,end=0.7)))) #Plots the color bar. Reverse colors here too.

#Saving plot.

map <- str_c(galaxy[i],"/",galaxy[i],"_Luminosity_Ha.eps")
dev.copy2eps(file=map)
map <- str_c("convert -density 300 ",galaxy[i],"/",galaxy[i],"_Luminosity_Ha.eps ",galaxy[i],"/",galaxy[i],"_Luminosity_Harms2.png")
system(map)


##########################################################################
##[O III 5007] Map

print('Plotting [O III] 5007 luminosity map.')

## Specify color map for z.

rbPal <-  colorRampPalette(rev(rainbow(100,start=0,end=0.7))) 
Col <- rbPal(100)[as.numeric(cut(log10(Lo),breaks = 100))]

title <- str_c(galaxy[i]," [O III] 5007 luminosity")
plot(Xo,Yo,xlim=c(-40,40),ylim=c(-40,40),xaxs= "i",yaxs="i", col="white", main = title,xlab = expression(paste(Delta, alpha," (arcsec)")), ylab = expression(paste(Delta, delta," (arcsec)")), cex.lab=1.3, cex.axis=1.3)
#rasterImage(ima,-40,-40,40,40,interpolate=F)
points(Xo,Yo, col = (Col), pch=15,cex = 0.7) #Plots spaxels with color map reversed (to have greater values darker and lower lighter).
grid() #Coordinates'grid.
abline(h = 0, v = 0, col = "gray60",lwd=2)
#mtext(side=3, expression(paste(" log"[10]," [O III] ",lambda,"5007 Luminosity [erg s"^-1,"]")))
image.plot( -40, -40, legend.only=TRUE, zlim= range(log10(Lo)), col = (rev(rainbow(100,start=0,end=0.7)))) #Plots the color bar. Reverse colors here too.

#Saving plot.

map <- str_c(galaxy[i],"/",galaxy[i],"_Luminosity_OIII5007.eps")
dev.copy2eps(file=map)
map <- str_c("convert -density 300 ",galaxy[i],"/",galaxy[i],"_Luminosity_OIII5007.eps ",galaxy[i],"/",galaxy[i],"_Luminosity_OIII5007rms2.png")
system(map)


##########################################################################
##[N II] 6584  Map

print('Plotting [N II] 6583 luminosity map.')

## Specify color map for z.

rbPal <-  colorRampPalette(rev(rainbow(100,start=0,end=0.7))) 
Col <- rbPal(100)[as.numeric(cut(log10(Ln),breaks = 100))]

title <- str_c(galaxy[i]," [N II] 6583 luminosity")
plot(Xn,Yn,xlim=c(-40,40),ylim=c(-40,40),xaxs= "i",yaxs="i", col="white", main= title, xlab = expression(paste(Delta, alpha," (arcsec)")), ylab = expression(paste(Delta, delta," (arcsec)")), cex.lab=1.3, cex.axis=1.3)
#rasterImage(ima,-40,-40,40,40,interpolate=F)
points(Xn,Yn, col = (Col), pch=15,cex = 0.7) #Plots spaxels with color map reversed (to have greater values darker and lower lighter).
grid() #Coordinates'grid.
abline(h = 0, v = 0, col = "gray60",lwd=2)
#mtext(side=3,expression(paste(" log"[10]," [N II] ",lambda," 6584 Luminosity [erg s"^-1,"]")))
image.plot( -40, -40, legend.only=TRUE, zlim= range(log10(Ln)), col = (rev(rainbow(100,start=0,end=0.7)))) #Plots the color bar. Reverse colors here too.

#Saving plot.

map <- str_c(galaxy[i],"/",galaxy[i],"_Luminosity_NII6583.eps")
dev.copy2eps(file=map)
map <- str_c("convert -density 300 ",galaxy[i],"/",galaxy[i],"_Luminosity_NII6583.eps ",galaxy[i],"/",galaxy[i],"_Luminosity_NII6583rms2.png")
system(map)


##########################################################################
##H-alpha equivalent width

#print('Plotting H-alpha equivalent width map.')

## Specify color map for z.

#rbPal <-  colorRampPalette(rev(rainbow(100,start=0,end=0.7))) 
#Col <- rbPal(100)[as.numeric(cut(log10(EWa),breaks = 100))]

#plot(X3,Y3,xlim=c(-40,40),ylim=c(-40,40),xaxs= "i",yaxs="i", col="white", main=galaxy[i], xlab = expression(paste(Delta, alpha," (arcsec)")), ylab = expression(paste(Delta, delta," (arcsec)")), cex.lab=1.3, cex.axis=1.3)
#rasterImage(ima,-40,-40,40,40,interpolate=F)
#points(X3,Y3, col = (Col), pch=15,cex = 0.7) #Plots spaxels with color map reversed (to have greater values darker and lower lighter).
#grid() #Coordinates'grid.
#abline(h = 0, v = 0, col = "gray60",lwd=2)
#mtext(side=3,expression(paste(" log"[10]," [N II] ",lambda," 6584 Luminosity [erg s"^-1,"]")))
#image.plot( -40, -40, legend.only=TRUE, zlim= range(log10(EWa)), col = (rev(rainbow(100,start=0,end=0.7)))) #Plots the color bar. Reverse colors here too.

#Saving plot.

#map <- str_c("/data2/CALIFA/",galaxy[i],".COMB/",galaxy[i],"_WHa.eps")
#dev.copy2eps(file=map)
#map <- str_c("convert -density 300 /data2/CALIFA/",galaxy[i],".COMB/",galaxy[i],"_WHa.eps /data2/CALIFA/",galaxy[i],".COMB/",galaxy[i],"_WHa.png")
#system(map)

#copy1 <- str_c("cp ",galaxy[i],".COMB/",galaxy[i],"_Luminosity*.png /home/aitor/Doctorado/Plots/Luminosities/")
#system(copy1)

#copy2 <- str_c("cp ",galaxy[i],".COMB/",galaxy[i],"_WHa.png /home/aitor/Doctorado/Plots/Equivalent_Width/")
#system(copy2)

}