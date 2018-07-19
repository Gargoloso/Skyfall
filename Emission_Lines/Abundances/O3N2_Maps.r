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

gametxy <- as.character(name)

####################################################################
####################################################################

##                   ITERATION FOR ALL GAmetXIES                   ##

####################################################################
####################################################################

##Loop for all gametxies
for(i in 1:length(gametxy)){
  
print('****************************')
print("NEW GAmetXY: ")
print(gametxy[i])
print('****************************')
  
####################################################################
####################################################################
  
##                EXTRACTION OF EMISSION-LINE DATA                ##
  
####################################################################
####################################################################

##Load data

print('Extracting base data.')  
path_base <- str_c(gametxy[i],"/Calzetti_Base_Fluxes.dat") #creates path to Ha and [N II]data file for each gametxy
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

# [N II] 6583
FN <- Fn[which(!is.na(AHa) & !is.na(rms_AHa) & AHa <= 4)]
rms_FN <- rms_Fn[which(!is.na(AHa) & !is.na(rms_AHa) & AHa <= 4)]

# [O III] 5007
FO <- Fo[which(!is.na(AHa) & !is.na(rms_AHa) & AHa <= 4)]
rms_FO <- rms_Fo[which(!is.na(AHa) & !is.na(rms_AHa) & AHa <= 4)]               

# Transmettion of coordinates setting origin to the field centre.
Y <- y - yc
X <- x - xc

##Convert from surface specific flux to standard specific intensity [erg cm^-2 s^-1 A^-1 sr^-1]

print('Converting surface specific units to specific intensities.')

crit1 <- rms_FA/FA <= 1.0 & rms_FN/FN <= 1.0 & rms_FB/FB <= 1.0 & rms_FO/FO <= 1.0
fa   <- FA[which(crit1)]/2.3504e-11 
fn   <- FN[which(crit1)]/2.3504e-11 
fb   <- FB[which(crit1)]/2.3504e-11
fo   <- FO[which(crit1)]/2.3504e-11 

XA  <- X[which(crit1)]
YA  <- Y[which(crit1)]
IDA <- ID2[which(crit1)]



################################################

##				Determination Metalicities 		      ##

################################################

#O3N2
print('Calculating O3N2 parametre')
O3N2 <- log10((fo/fb)*((fa/fn)))

print('Calculating 12 + log(O/H)')
met <- 8.533-(0.214*O3N2[which(O3N2 >= -1.1 & O3N2 <= 1.7)])
O3N2m <- O3N2[which(O3N2 >= -1.1 & O3N2 <= 1.7)]
Xm <- XA[which(O3N2 >= -1.1 & O3N2 <= 1.7)]
Ym <- YA[which(O3N2 >= -1.1 & O3N2 <= 1.7)]
IDm <- IDA[which(O3N2 >= -1.1 & O3N2 <= 1.7)]


#########################################       
#MAPS ONLY!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#removing inferior quartile.
# remove_outliers <- function(x, na.rm=TRUE, ...){
#   qnt <- quantile(x, probs=c(.10,.99, na.rm = na.rm, ...))
#   H <- 1.5*IQR(x, na.rm = na.rm)
#   clean <- x
#   clean[x < (qnt[1] - H)] <- NA
#   clean[x > (qnt[2] + H)] <- NA
#   clean
# }

# met2 <- 10^(remove_outliers(log10(met1)))
# met  <- met2[which(!is.na(met2))]/(asp * 3.826e33) #Convert to Sometr luminosity / pc^2.
# Xa  <- XA[which(!is.na(met2))]
# Ya  <- YA[which(!is.na(met2))]
# IDa <-  sprintf("%04d", as.numeric(IDA[which(!is.na(met2))]))

print('Saving data to file')

resume <- data.frame(IDm,O3N2m,met)
tabmet <- str_c(gametxy[i],"/",gametxy[i],"_O3N2.dat")
write.table(resume, tabmet, sep="\t",quote=FALSE, row.names = FALSE)


##########################################################################

##				PLOTS					##

##########################################################################

print("Plotting maps")

par(mar=c(5,5,4.75,4.75))
par(bg='white')

####################################
#Set background image. Only once for all maps.

#ima <- readPNG(str_c(gametxy[i],".COMB/",gametxy[i],".png")) #Background image file in png.
#lim <- par()

##########################################################################
##H-alpha Map

print('Plotting 12 + O/H map.')

## Specify color map for z.

rbPal <-  colorRampPalette(rev(rainbow(100,start=0,end=0.7))) 
Col <- rbPal(100)[as.numeric(cut((met),breaks = 100))]

title <- str_c(gametxy[i]," 12 + log(O/H)")
plot(Xm,Ym,xlim=c(-40,40),ylim=c(-40,40),xaxs= "i",yaxs="i", col="white", main = title,xmetb = expression(paste(Delta, alpha," (arcsec)")), ymetb = expression(paste(Delta, delta," (arcsec)")), cex.metb=1.3, cex.axis=1.3)
#rasterImage(ima,-40,-40,40,40,interpomette=F)
points(Xm,Ym, col = (Col), pch=15,cex = 0.7) #Plots spaxels with color map reversed (to have greater values darker and lower lighter).
grid() #Coordinates'grid.
abline(h = 0, v = 0, col = "gray60",lwd=2)
#mtext(side=3,expression(paste(" log"[10]," H",alpha," Luminosity [erg s"^-1,"]")))
image.plot( -40, -40, legend.only=TRUE, zlim= range((met)), col = (rev(rainbow(100,start=0,end=0.7)))) #Plots the color bar. Reverse colors here too.

#Saving plot.

map <- str_c(gametxy[i],"/",gametxy[i],"_O3N2.eps")
dev.copy2eps(file=map)
map <- str_c("convert -density 300 ",gametxy[i],"/",gametxy[i],"_O3N2.eps ",gametxy[i],"/",gametxy[i],"_O3N2.png")
system(map)

}