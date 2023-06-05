############################
# spatial preference heterogeneity analysis of the data of the BonaRes soil ecosystem services CE
# n=1500, conducted online between 21 June and 6 July 2021
# authors: Bartosz Bartkowski, Nele Lienhoop, Lukas Mahlich, Julian Massenberg 
# contact: bartosz.bartkowski@ufz.de
############################

############## Outline ###############
# 1. Packages and data
# 2. Calculate WTP means for PLZ areas
# 3. Spatial weighting matrices
# 4. Calculate Moran's I
# 5. Robustness check
#  5.1. Jitter coordinates
#  5.2. Recalculate spatial weighting matrix
#  5.3. Recalculate Moran's I
######################################

######################
# 1. Packages and data
######################

#load packages
require(spdep)
require(ape) #in the most cases, the more efficient option compared to spdep
require(here)

#import dataset from CE including individual WTP scores (code for calculation of individual WTP scores from mixed logit in WTP space available in soil_CE_analysis_main_study_apollo.R), coordinates of the postal codes (PLZ) centroids as well as additional data on indicators of drought impacts, flood impacts, groundwater quality and soil quality
ce_data <- read.csv(here("Spatial data/ce_results_spatial_data_clean.csv"), header = T)

######################
# 2. Calculate WTP means for PLZ areas
######################

#calculate mean WTP scores for each PLZ (postcode)
wtp_plz <- as.data.frame(matrix(nrow=0,ncol=7)) #create empty dataframe
unique_plz <- unique(ce_data$plz) #create vector of unique PLZ
unique_x <- unique(ce_data$xcoor) #same for x coordinates
unique_y <- unique(ce_data$ycoor) #...and for y coordinates

#generate a 1-row matrix for each PLZ with mean WTP scores and add it to dataframe
for (i in 1:length(unique_plz)){
  this.plz <- matrix(nrow = 1, ncol = 7)
  this.plz[,1] <- unique_plz[i]
  this.plz[,2] <- mean(subset(ce_data,plz==unique_plz[i])$drought)
  this.plz[,3] <- mean(subset(ce_data,plz==unique_plz[i])$flood)
  this.plz[,4] <- mean(subset(ce_data,plz==unique_plz[i])$climate)
  this.plz[,5] <- mean(subset(ce_data,plz==unique_plz[i])$water)
  this.plz[,6] <- unique_x[i]
  this.plz[,7] <- unique_y[i]
  wtp_plz <- rbind(wtp_plz, this.plz)
}
rm(this.plz,unique_plz,unique_x,unique_y,i)
colnames(wtp_plz) <- c("plz","drought","flood","climate","water","xcoor","ycoor") #rename columns

######################
# 3. Spatial weighting matrices
######################

#extract coordinates
coords <- wtp_plz[,c(6:7)]

#full inverse-distance weighting
dists <- as.matrix(dist(coords)) #make matrix of distances
dists_i <- 1/dists #invert
diag(dists_i) <- 0 #add diagonal
row_sums <- apply(dists_i,1,sum) #calculate row sums
dists_i <- dists_i/row_sums #normalize values using row sums
#weights_all <- mat2listw(dists_i,style="W") #reformat to have format readable by spdep

#inverse-distance weighting for 41494 m range (in line with default settings in ArcGIS)
dists_40 <- 1/dists
diag(dists_40) <- 0
dists_40[dists_40<1/41494] <- 0
row_sums <- apply(dists_40,1,sum)
dists_40 <- dists_40/row_sums
#weights_40 <- mat2listw(dists_40,style="W")

#inverse-distance weighting within (arbitrary) 200000 m range
dists_200 <- 1/dists
diag(dists_200) <- 0
dists_200[dists_200<1/200000] <- 0
row_sums <- apply(dists_200,1,sum)
dists_200 <- dists_200/row_sums
#weights_200 <- mat2listw(dists_200,style="W")

#k-nearest neighbours with k=8
col.knn <- knearneigh(as.matrix(coords), k = 8)
weights_knn <- nb2listw(knn2nb(col.knn), zero.policy = T)

######################
# 4. Calculate Moran's I
######################

#spatial autocorrelation with full weighting matrix
Moran.I(wtp_plz$drought, dists_i)
Moran.I(wtp_plz$flood, dists_i)
Moran.I(wtp_plz$climate, dists_i)
Moran.I(wtp_plz$water, dists_i)

#spatial autocorrelation with 40-km matrix
Moran.I(wtp_plz$drought, dists_40)
Moran.I(wtp_plz$flood, dists_40)
Moran.I(wtp_plz$climate, dists_40)
Moran.I(wtp_plz$water, dists_40)

#spatial autocorrelation with 200-km matrix
Moran.I(wtp_plz$drought, dists_200)
Moran.I(wtp_plz$flood, dists_200)
Moran.I(wtp_plz$climate, dists_200)
Moran.I(wtp_plz$water, dists_200)

#spatial autocorrelation with k-nearest neighbours (k=8)
moran.test(wtp_plz$drought, weights_knn, alternative = "two.sided")
moran.test(wtp_plz$flood, weights_knn, alternative = "two.sided")
moran.test(wtp_plz$climate, weights_knn, alternative = "two.sided")
moran.test(wtp_plz$water, weights_knn, alternative = "two.sided")

##spdep example
#spatial autocorrelation with full weighting matrix
#moran.test(wtp_plz$drought, weights_all, alternative = "two.sided")
#moran.test(wtp_plz$flood, weights_all, alternative = "two.sided")
#moran.test(wtp_plz$climate, weights_all, alternative = "two.sided")
#moran.test(wtp_plz$water, weights_all, alternative = "two.sided")

######################
# 5. Robustness check
######################

######################
# 5.1. Jitter coordinates
######################

#read in data again
ce_data <- read.csv(here("Spatial data/ce_results_spatial_data_clean.csv"),header = T)

#identify duplicates (i.e. observations in same PLZ area)
duplicates <- which(duplicated(ce_data$plz))

#check the range of coordinates
summary(ce_data$xcoor)
summary(ce_data$ycoor)

#manually jitter the coordinates of the duplicates (manually because I can't find a good alternative to the spJitter function from the discontinued spsann package)
set.seed(2104)
for (i in 1:length(duplicates)) {
  ce_data$xcoor[duplicates[i]] <- ce_data$xcoor[duplicates[i]] + runif(1,-20,20)
  ce_data$ycoor[duplicates[i]] <- ce_data$ycoor[duplicates[i]] + runif(1,-20,20)
}
#check results (if everything worked well, there shouldn't be any observations left that have the same coordinates)
which(duplicated(ce_data$xcoor) & duplicated(ce_data$ycoor))

#save modified data
write.csv(ce_data, here("Spatial data/ce_results_spatial_data_jittered.csv"), fileEncoding = "UTF-8")

######################
# 5.2. Recalculate weighting matrix
######################

#extract coordinates
coords <- ce_data[,c("xcoor","ycoor")]

#full inverse-distance weighting
dists <- as.matrix(dist(coords)) #make matrix of distances
dists_i <- 1/dists #invert
diag(dists_i) <- 0 #add diagonal
row_sums <- apply(dists_i,1,sum) #calculate row sums
dists_i <- dists_i/row_sums #normalize values using row sums

#inverse-distance weighting for 41494 m range (in line with default settings in ArcGIS)
dists_40 <- 1/dists
diag(dists_40) <- 0
dists_40[dists_40<1/41494] <- 0
row_sums <- apply(dists_40,1,sum)
dists_40 <- dists_40/row_sums

######################
# 5.3. Recalculate Moran's I
######################

#spatial autocorrelation with jittered coordinates and full inverse-distance weighting
Moran.I(ce_data$drought, dists_i)
Moran.I(ce_data$flood, dists_i)
Moran.I(ce_data$climate, dists_i)
Moran.I(ce_data$water, dists_i)

#spatial autocorrelation with jittered coordinates and "ArcGIS weighting"
Moran.I(ce_data$drought, dists_40)
Moran.I(ce_data$flood, dists_40)
Moran.I(ce_data$climate, dists_40)
Moran.I(ce_data$water, dists_40)
