#########################################################################
#########################################################################


##	This programs finds draws BPT-SII maps for all galaxies.

##	02/23/2017 by A. Robleto-Orus
##      02/13/2018 Modified for more general purpose by A. Robleto-Orus

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

data <- read.table("lis.dat", header=TRUE) # File with general information for all galaxies.
attach(data)
galaxy<-as.character(name)


####################################################################
####################################################################

##		     ITERATION FOR ALL GALAXIES                   ##

####################################################################
####################################################################

nSy <- seq(1,length(galaxy),1)
nSF <- seq(1,length(galaxy),1)
nAn <- seq(1,length(galaxy),1)
tot <- seq(1,length(galaxy),1)

#Iteration for all galaxies.
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
  
  ##Load data
  
  path_base <- str_c(galaxy[i],'/Calzetti_Base_Fluxes.dat') #creates path to Ha and [N II]data file for each galaxy
  data0 <- read.table(path_base, header=TRUE)
  
  path_oi6300 <- str_c(galaxy[i],"/Calzetti_oi6300_Fluxes.dat") #creates path to Ha and [S II]data file for each galaxy
  data1 <- read.table(path_oi6300, header=TRUE)
  
  ##Merge data
  DATA <- merge(data0, data1,by.x = 1, by.y =1)
  attach(DATA)
  
  ID2 <- ID[which(!is.na(AHa) & !is.na(rms_AHa) & AHa <= 4& !is.na(Foi6300) & !is.na(rms_Foi6300))] # ID of each spectra.
  
  ##Extracting coordisates from ID
  ## Spaxel's coordisates are in the id in format YYXX, with Y = DEC and X = RA.
  y <- (trunc(ID[which(!is.na(AHa) & !is.na(rms_AHa) & AHa <= 4 & !is.na(Foi6300) & !is.na(rms_Foi6300))]/100)  ) #Obtaining Y from id.
  
  x <- (ID[which(!is.na(AHa) & !is.na(rms_AHa) & AHa <= 4 & !is.na(Foi6300) & !is.na(rms_Foi6300))]-(y*100) ) #Obtaining X from id.
  
  ##Extracting line surface specific intensities [1e-16 erg cm^-2 s^-1 A^-1 arcsec^-2] Conditions must be met by all lines.
  
  #H-alpha 6563
  FA <-    Fa[which(!is.na(AHa) & !is.na(rms_AHa) & AHa <= 4 & !is.na(Foi6300) & !is.na(rms_Foi6300))] #Obtaining the[N II] flux density (10^-16 erg/(s cm^2 A pix))
  rms_FA <- rms_Fa[which(!is.na(AHa) & !is.na(rms_AHa) & AHa <= 4 & !is.na(Foi6300) & !is.na(rms_Foi6300))]
  
  #[S II] 6716
  FOI <-    Foi6300[which(!is.na(AHa) & !is.na(rms_AHa) & AHa <= 4 & !is.na(Foi6300) & !is.na(rms_Foi6300))] #Obtaining the[N II] flux density (10^-16 erg/(s cm^2 A pix))
  rms_FOI <- rms_Foi6300[which(!is.na(AHa) & !is.na(rms_AHa) & AHa <= 4 & !is.na(Foi6300) & !is.na(rms_Foi6300))]
  
   #H-beta 4861
  FB <-    Fb[which(!is.na(AHa) & !is.na(rms_AHa) & AHa <= 4 & !is.na(Foi6300) & !is.na(rms_Foi6300))] #Obtaining the[N II] flux density (10^-16 erg/(s cm^2 A pix))
  rms_FB <- rms_Fb[which(!is.na(AHa) & !is.na(rms_AHa) & AHa <= 4 & !is.na(Foi6300) & !is.na(rms_Foi6300))]
  
  #[O III] 5007
  FO <-    Fo[which(!is.na(AHa) & !is.na(rms_AHa) & AHa <= 4 & !is.na(Foi6300) & !is.na(rms_Foi6300))] #Obtaining the[N II] flux density (10^-16 erg/(s cm^2 A pix))
  rms_FO <- rms_Fo[which(!is.na(AHa) & !is.na(rms_AHa) & AHa <= 4 & !is.na(Foi6300) & !is.na(rms_Foi6300))]
  
  # Translation of coordinates setting origin to the field centre.
  Y <- y - yc
  X <- x - xc
  
  
  ########################################################################
  
  ##			LINE RATIOS			              ##
  
  ########################################################################
  ## Note: rms are calculated from the non-logarythmic line ratios and 
  # applying logarythm after error propagation. Direct logarythmic propagation
  # should be on the form: #rms_na1 <- sqrt(((1/(FN*log(10)))*rms_FN)^2 + ((-1/(FA*log(10)))*rms_FA)^2)
  
  #[N II]6584 / H-alpha
  o1a1 <- log10(FOI/FA)
  rms_o1a1 <- sqrt(((1/FA)*rms_FOI)^2 + ((-FOI/(FA^2))*rms_FA)^2)
  rms_o1a2 <- log10(rms_o1a1)
  
  #[O III] 5007 / H-beta
  ob1 <- log10(FO/FB)
  rms_ob2 <- sqrt(((1/FB)*rms_FO)^2 + ((-FO/(FB^2))*rms_FB)^2)
  rms_ob3 <- log10(rms_ob2)
  
  ## Filtering out NAs and relative errors >= 1.
  
  o1a <- o1a1[which(!is.na(o1a1) & !is.na(ob1) & (rms_o1a1/(FOI/FA)) <= 1 & (rms_ob2/(FO/FB)) <= 1)]
  rms_o1a <- rms_o1a2[which(!is.na(o1a1) & !is.na(ob1) & (rms_o1a1/(FOI/FA)) <= 1 & (rms_ob2/(FO/FB)) <= 1)]
  
  ob <- ob1[which(!is.na(o1a1) & !is.na(ob1) & (rms_o1a1/(FOI/FA)) <= 1 & (rms_ob2/(FO/FB)) <= 1)]
  rms_ob <- rms_ob3[which(!is.na(o1a1) & !is.na(ob1) & (rms_o1a1/(FOI/FA)) <= 1 & (rms_ob2/(FO/FB)) <= 1)]
  
  X2 <- X[which(!is.na(o1a1) & !is.na(ob1) & (rms_o1a1/(FOI/FA)) <= 1 & (rms_ob2/(FO/FB)) <= 1)]
  Y2 <- Y[which(!is.na(o1a1) & !is.na(ob1) & (rms_o1a1/(FOI/FA)) <= 1 & (rms_ob2/(FO/FB)) <= 1)]
  
  ID3 <- ID2[which(!is.na(o1a1) & !is.na(ob1) & (rms_o1a1/(FOI/FA)) <= 1 & (rms_ob2/(FO/FB)) <= 1)]
  
  ########################################################################
  
  ##			BPT CLASSIFICATION			      ##}
  
  ########################################################################
  
  ##Separating by object type  
  mainAGN  <- function(x) (0.72/(x -0.32)) + 1.30 
  Kewley06 <- function(x) (1.89*x) + 0.76
  new <-   function (x) (0.56115520*x^2) + (2.02573341*x) - 0.02954284
  
  type <- seq(1,length(o1a),1)
  for(jj in 1:length(o1a)){
    if(ob[jj] <= mainAGN(o1a[jj]) & o1a[jj] < 0.05) {type[jj] <- "SF"} #Star Forming spaxels
    else if(ob[jj] <= Kewley06(o1a[jj])) {type[jj] <- "LINER"} #LINER-like, uncomment to use Cid Fersandes division and comment the modif line.
    else {type[jj] <- "Sy"} # Seyfert
  } 
  
  ##########################################################################
  
  ##				PLOTS					##
  
  ##########################################################################
  #map <- str_c("Maps/BPT/",galaxy[i],"_BPT_DEP_1.0_radius.png")
  #png(map,width=2200, height=2200,bg="white")
  print("Plotting maps")
  
  #############################################################################
  ##BPT Map
  
  ##Plot BPT diagram
  par(mar=c(5,5,4.75,4.75))
  plot(o1a,ob,pch=19, main = galaxy[i], ylab= expression(paste("log"[10]," ([O III] ",lambda,"5007 / H",beta,")")), xlab= expression(paste("log"[10]," ([O I] ",lambda,"6300 / H",alpha,")")),xlim=c(-1.5,1.5),ylim=c(-1.5,1.5), col="white", cex.lab=1.4, cex.axis=1.4,asp=1)
  points(o1a[type=="SF"],ob[type== "SF"], col = "blue", pch=19, cex=0.7)
  points(o1a[type=="LINER"],ob[type=="LINER"], col = "orange", pch=19, cex=0.7)
  points(o1a[type=="Sy"],ob[type=="Sy"], col = "red", pch=19, cex=0.7)
  #contour(kde_BPT, levels= levels,labels = c(4,3,2.5,2,1.5,1,0.5),add=TRUE, lwd=1.5, labcex=1)
  curve(mainAGN,add=TRUE,lty=1,col="black",from=-2.3,to=0.3,lwd=2.5)
  curve(Kewley06,add=TRUE,lty=3,col="black",from=-0.28,to=1.5,lwd=2.5)
  curve(new,add=TRUE,lty=5,col="purple",from=-0.1,to=1.3,lwd=2.5)
  text(-1,1.4, label="Seyfert",cex=1.5)
  text(-1.3,-1.5, label="SF",cex=1.5)
  text(0.8,-1.5, label="Anomalous",cex=1.5)
  #legend('topright',legend = c('order 2, 2-sigma fit','order 4, 2-sigma fit'), lty = c(5,5),lwd=c(1.5,1.5), col=c("purple","Sy"),bg = 'white')
  
  ##############################################################################
  ### Save plots as eps and png files
  
  print('SAVING PLOT')
  
  map <- str_c(galaxy[i],"/",galaxy[i],"_BPTdiag_OI.eps")
  dev.copy2eps(file=map)
  map <- str_c("convert -density 300 ",galaxy[i],"/",galaxy[i],"_BPTdiag_OI.eps ",galaxy[i],"/",galaxy[i],"_BPTdiag_OI.png")
  system(map)
  
  ##############################################################################
  ### Plot BPT-SII map.
  print("PLOTTING MAP")
  
  colour <- seq(1:length(X2))
  for (cc in 1:length(X2)) {
    if (type[cc] == "Sy"){colour[cc] <- "red"}
    if (type[cc] == "SF"){colour[cc] <- "blue"}
    if (type[cc] == "LINER"){colour[cc] <- "orange"}
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
  
  map <- str_c(galaxy[i],"/",galaxy[i],"_BPT_OI_Map.eps")
  dev.copy2eps(file=map)
  map <- str_c("convert -density 300  ",galaxy[i],"/",galaxy[i],"_BPT_OI_Map.eps ",galaxy[i],"/",galaxy[i],"_BPT_OI_Map.png")
  system(map)
  
  ##############################################################################
  ##Save data to files.
  
  print('Saving line ratios and BPT types to data file.')
  resume <- data.frame(ID3, o1a, rms_o1a, ob, rms_ob, type)
  tabla <- str_c(galaxy[i],"/BPT-SII_data.dat")
  write.table(resume, tabla, sep="\t",quote=FALSE)
  
  ##############################################################################
  ##Count number of spaxels of each type  
  
  nSy[i] <- length(o1a[which(type=="Sy")])
  nSF[i] <- length(o1a[which(type=="SF")])
  nAn[i] <- length(o1a[which(type=="LINER")])
  tot[i] <- nSy[i] + nSF[i] + nAn[i]
  
}

##############################################################################
##Save summary data to a file.

Cuentas <- data.frame(galaxy,nSy,nSF,nAn, tot)
tabla <- str_c('BPT-OI_SpaxCount.dat')
write.table(Cuentas, tabla, sep="\t",quote=FALSE)

