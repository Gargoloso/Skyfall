#########################################################################
#########################################################################


##	This programs finds draws BPT diagrams  and maps for all galaxies. ##
                                        
##  07/02/2018 by A. Robleto-Orús                                      ##

#########################################################################
#########################################################################

##Clean the workspace
rm(list=ls(all=TRUE))

##Libraries
require("stringr")
library("png")


########################################################################

##				DATA INPUT			      ##

########################################################################

setwd("~/Rings/ringed_work/") # Directrory with our data.

data <- read.table("lis.dat", header=TRUE) #File with general information for all galaxies.
attach(data)
galaxy<-as.character(name)


####################################################################
####################################################################

##		     ITERATION FOR ALL GALAXIES                   ##

####################################################################
####################################################################



nSy <- seq(1,length(galaxy),1)
nSF <- seq(1,length(galaxy),1)
nTO <- seq(1,length(galaxy),1)
nAn <- seq(1,length(galaxy),1)
tot <- seq(1,length(galaxy),1)


for (i in 1:length(galaxy)){
print("NEW GALAXY: ")
print(galaxy[i])

##id of central spaxel
xc <- trunc(naxis1[i]/2)
yc <- trunc(naxis2[i]/2)


####################################################################
####################################################################

##		  EXTRACTION OF EMISSION-LINE DATA                ##

####################################################################
####################################################################

print('Extracting Base data') #Creates path to a file with Base results.
pathcb <- str_c(galaxy[i],"/Calzetti_Base_Fluxes.dat")
data1 <- read.table(pathcb, header=TRUE) 
attach(data1)

print('EXTRACTING DESyDENED LINE SURFACE SPECIFIC FLUXES')

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

# Translation of coordinates setting origin to the field centre.
Y <- y - yc
X <- x - xc



########################################################################

##			LINE RATIOS			              ##

########################################################################

print('CALCULATING LINE RATIOS')
## Note: rms are calculated from the non-logarythmic line ratios and 
# applying logarythm after error propagation. Direct logarythmic propagation
# should be on the form: #rms_na1 <- sqrt(((1/(FN*log(10)))*rms_FN)^2 + ((-1/(FA*log(10)))*rms_FA)^2)

#[N II]6584 / H-alpha
na1 <- log10(FN/FA)
rms_na2 <- sqrt(((1/FA)*rms_FN)^2 + ((-FN/(FA^2))*rms_FA)^2)
rms_na3 <- log10(rms_na2)

#[O III] 5007 / H-beta
ob1 <- log10(FO/FB)
rms_ob2 <- sqrt(((1/FB)*rms_FO)^2 + ((-FO/(FB^2))*rms_FB)^2)
rms_ob3 <- log10(rms_ob2)

## Filtering out NAs.

na     <- na1[which(!is.na(na1) & !is.na(ob1))]
rms_na <- rms_na3[which(!is.na(na1) & !is.na(ob1))]

ob     <- ob1[which(!is.na(na1) & !is.na(ob1))]
rms_ob <- rms_ob3[which(!is.na(na1) & !is.na(ob1))]

X2 <- X[which(!is.na(na1) & !is.na(ob1))]
Y2 <- Y[which(!is.na(na1) & !is.na(ob1))]


########################################################################

##			BPT CLASSIFICATION			      ##

########################################################################

print('BPT CLASSIFICATION')

##Separating by object type  
kewley<-function(x) (0.61/(x-0.47))+1.19 #Kewley 2001 maximum starburst line
kauffmann<-function(x) (0.61/(x-0.05))+1.3 #Kauffman 2003 mixing line
cid<-function(x) 1.01*x+0.48 #Cid Fernandes 2010 Sy/LINER division line
#new <- function(x) (3.45367*x^2) + (3.82356*x) + 0.18961
new <- function(x) (1.80282*x^2) + (2.58101*x) - 0.25599
#modif <- function(x) (m[i]*x)+b[i] #Proposed modified division for Sy/LINER
     
type <- seq(1,length(na),1)
for(jj in 1:length(na)){
  if(ob[jj] <= kauffmann(na[jj]) & na[jj] < 0.05) {type[jj] <- "SF"} #Star Forming spaxels
  #else if(ob[jj] <= modif(na[jj])) {type[jj] <- "LINER"} #LINER-like proposed division, uncomment to use it and comment the cid line.
  else if(ob[jj] <= kewley(na[jj]) & na[jj]< 0.47){type[jj] <- "TO"} #Transition Object
  else if(ob[jj] <= cid(na[jj])) {type[jj] <- "LINER"} #LINER-like, uncomment to use Cid Fernandes division and comment the modif line.
  else {type[jj] <- "Sy"} # Seyfert
   } 
   
##########################

##				PLOTS					##

###########################
print("PLOTTING")

#############################################################################
##BPT diagram

##Plot BPT diagram
par(mar=c(5,5,4.75,4.75))
plot(na,ob,pch=19, ylab= expression(paste("log"[10]," [O III] ",lambda,"5007 / H",beta)), main = galaxy[i], xlab= expression(paste("log"[10]," [N II] ",lambda,"6584 / H",alpha)),xlim=c(-2,1),ylim=c(-1.5,1.5), col="white", cex.lab=1.4, cex.axis=1.4, asp=1)
points(na[type=="SF"],ob[type== "SF"], col = "blue", pch=19, cex=0.7)
points(na[type=="TO"],ob[type=="TO"], col = "green", pch=19, cex=0.7)
points(na[type=="LINER"],ob[type=="LINER"], col = "orange", pch=19, cex=0.7)
points(na[type=="Sy"],ob[type=="Sy"], col = "red", pch=19, cex=0.7)
#contour(kde_BPT, levels= levels,labels = c(4,3,2.5,2,1.5,1,0.5),add=TRUE, lwd=1.5, labcex=1)
curve(kewley,add=TRUE,lty=1,col="black",from=-2.3,to=0.3,lwd=2.5)
curve(kauffmann,add=TRUE,lty=2,col="black",from=-1.28,to=0,lwd=2.5)
curve(cid,add=TRUE,lty=3,col="black",from=-0.43,to=1.3,lwd=2.5)
curve(new,add=TRUE,lty=5,col="purple",from=-0.25,to=1.3,lwd=2.5)
#curve(Model4,add=TRUE,lty=5,col="Sy",from=-0.43,to=1.3,lwd=1.5)
#lines(fitted(model),X_na,lty=5,col="black",lwd=1.5)
text(-0.5,1.4, label="Seyfert",cex=1.5)
text(-2,-1.5, label="SF",cex=1.5)
text(0,-1.5, label="TO",cex=1.5)
text(0.7,-1.5, label="LINER",cex=1.5) #Anomalous


##############################################################################
### Save plots as eps and png files

print('SAVING PLOT')

map <- str_c(galaxy[i],"/",galaxy[i],"_BPTdiag_NII.eps")
dev.copy2eps(file=map)
map <- str_c("convert -density 300 ",galaxy[i],"/",galaxy[i],"_BPTdiag_NII.eps ",galaxy[i],"/",galaxy[i],"_BPTdiag_NII.png")
system(map)

##############################################################################
### Plot BPT-NII map.
print("PLOTTING MAP")

colour <- seq(1:length(X2))
for (cc in 1:length(X2)) {
  if (type[cc] == "Sy"){colour[cc] <- "red"}
  if (type[cc] == "SF"){colour[cc] <- "blue"}
  if (type[cc] == "LINER"){colour[cc] <- "orange"}
  if (type[cc] == "TO"){colour[cc] <- "green"}
}

#Set background image
#conv <-str_c("convert -density -300 ",galaxy[i],"/",galaxy[i],".jpg ",galaxy[i],"/",galaxy[i],".png")
#system(conv)
#ima <- readPNG(str_c(galaxy[i],"/",galaxy[i],".png")) #Background image file in png.
#lim <- par()

par(mar=c(5,5,4,4))
par(bg="white")
plot(X2,Y2,col="white", main = galaxy[i], xlim=c(-40,40),ylim=c(-40,40),xaxs= "i",yaxs="i",xlab = expression(paste(Delta, alpha, " (arcsec)")), ylab = expression(paste(Delta, delta, " (arcsec)")), cex.lab=1.3, cex.axis=1.3,pty='s')
#rasterImage(ima,-40,-40,40,40,interpolate=F) #Inserts background image between specified coordinates.
points(X2,Y2,pch=15, col = colour, cex = 0.8)

grid()
abline(h = 0, v = 0, col = "gray60",lwd=2)


##############################################################################
### Save plots as eps and png files

print('SAVING PLOT')

map <- str_c(galaxy[i],"/",galaxy[i],"_BPT_NII_Map.eps")
dev.copy2eps(file=map)
map <- str_c("convert -density 300  ",galaxy[i],"/",galaxy[i],"_BPT_NII_Map.eps ",galaxy[i],"/",galaxy[i],"_BPT_NII_Map.png")
system(map)

##############################################################################
##Save data to files.

print('Saving line ratios and BPT types to data file.')
resume <- data.frame(ID2, na, rms_na, ob, rms_ob, type)
tabla <- str_c(galaxy[i],"/BPT-NII_data.dat")
write.table(resume, tabla, sep="\t",quote=FALSE)

##############################################################################
##Count number of spaxels of each type  

nSy[i] <- length(na[which(type=="Sy")])
nSF[i] <- length(na[which(type=="SF")])
nTO[i] <- length(na[which(type=="TO")])
nAn[i] <- length(na[which(type=="LINER")])
tot[i] <- nSy[i] + nSF[i] + nTO[i] + nAn[i]
}

##############################################################################
##Save summary data to a file.

Cuentas <- data.frame(galaxy,nSy,nSF,nTO,nAn, tot)
tabla <- str_c('BPT-NII_SpaxCount.dat')
write.table(Cuentas, tabla, sep="\t",quote=FALSE)

