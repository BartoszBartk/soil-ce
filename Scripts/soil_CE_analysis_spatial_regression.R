############################
# spatial preference heterogeneity analysis of the data of the BonaRes soil ecosystem services CE
# n=1500, conducted online between 21 June and 6 July 2021
# authors: Bartosz Bartkowski, Nele Lienhoop, Lukas Mahlich, Julian Massenberg
# contact: bartosz.bartkowski@ufz.de
############################

############## Outline ###############
# 1. Packages and data
# 2. Weighting matrix
# 3. Spatial lag models
# 4. Spatial error models
# 5. Spatial Durbin models
######################################

######################
# 1. Packages and data
######################

#optional: clear environment
rm(list=ls())

#load packages
require(here)
require(spatialreg)
require(spdep)
require(rgdal)
require(sp)
require(dplyr)
require(modelsummary)

#import dataset from CE including individual WTP scores (code for calculation of individual WTP scores from mixed logit in WTP space available in soil_CE_analysis_main_study_apollo.R), coordinates of the postal codes (PLZ) centroids as well as additional data on indicators of drought impacts, flood impacts, groundwater quality and soil quality
ce_data <- read.csv(here("Spatial data/ce_results_spatial_data_clean.csv"),header = T)

######################
# 2. Weighting matrix
######################

##create weighting matrix for spatial lags
##note: spatial weighting matrices from soil_CE_analysis_spatial_autocorrelation.R cannot be used here, as there, there is only one aggregated observation per PLZ area
#full weighting
coords <- as.matrix(ce_data[,c(27:28)]) #extract coordinates
dists <- as.matrix(dist(coords)) #make matrix of distances
dists[dists == 0] <- 1 #correct to avoid dividing by 0
dists_i <- 1/dists #invert
diag(dists_i) <- 0 #add diagonal
row_sums <- apply(dists_i,1,sum) #calculate row sums
dists_i <- dists_i/row_sums #normalize values using row sums
weights_all <- mat2listw(dists_i,style="W") #reformat to have format readable by spatialreg
#neighbs <- listw2sn(weights_list)

#only include neighbours within 41494 m range (in line with default settings in ArcGIS)
dists_40 <- 1/dists
diag(dists_40) <- 0
dists_40[dists_40<1/41494] <- 0
row_sums <- apply(dists_40,1,sum)
dists_40 <- dists_40/row_sums
weights_40 <- mat2listw(dists_40,style="W")

#only include neighbours within (arbitrary) 200000 m range
dists_200 <- 1/dists
diag(dists_200) <- 0
dists_200[dists_200<1/200000] <- 0
row_sums <- apply(dists_200,1,sum)
dists_200 <- dists_200/row_sums
weights_200 <- mat2listw(dists_200,style="W")

#k-nearest neighbours with k=8
col.knn <- knearneigh(as.matrix(coords), k = 8, use_kd_tree = F)
weights_knn <- nb2listw(knn2nb(col.knn), zero.policy = T)

######################
# 3. Spatial lag models
######################

######################
# 3.1 Hypothesis 1 Drought
######################

#first, use SMI_Min as indicator to see differences between different spatial weighting schemes
#inverse distances, full matrix
splag_drought_smi_dist_all <- lagsarlm(drought ~ SMI_Min + age + gender + edu + income_ + member + donation + no_ag,
                                   data = ce_data,
                                   listw = weights_all,
                                   zero.policy = T, 
                                   na.action = na.omit)
#save model output to docx file
modelsummary(models = splag_drought_smi_dist_all,
             output = "Output/splag_drought_smi_dist_all.docx",
             fmt = 3,
             estimate = "{estimate} [{conf.low}, {conf.high}]",
             statistic = "{p.value}")

#inverse distances, only points within 40 km
splag_drought_smi_dist_40 <- lagsarlm(drought ~ SMI_Min + age + gender + edu + income_ + member + donation + no_ag,
                                       data = ce_data,
                                       listw = weights_40,
                                       zero.policy = T, 
                                       na.action = na.omit)
#save model output to docx file
modelsummary(models=splag_drought_smi_dist_40,
             output="Output/splag_drought_smi_dist_40.docx",
             fmt = 3,
             estimate = "{estimate} [{conf.low}, {conf.high}]",
             statistic = "{p.value}")

#inverse distances, only points within 200 km
splag_drought_smi_dist_200 <- lagsarlm(drought ~ SMI_Min + age + gender + edu + income_ + member + donation + no_ag,
                                       data = ce_data,
                                       listw = weights_200,
                                       zero.policy = T, 
                                       na.action = na.omit)
#save model output to docx file
modelsummary(models=splag_drought_smi_dist_200,
             output="Output/splag_drought_smi_dist_200.docx",
             fmt = 3,
             estimate = "{estimate} [{conf.low}, {conf.high}]",
             statistic = "{p.value}")

#k-nearest neighbours (k=8)
splag_drought_smi_dist_knn <- lagsarlm(drought ~ SMI_Min + age + gender + edu + income_ + member + donation + no_ag,
                                       data = ce_data,
                                       listw = weights_knn,
                                       zero.policy = T, 
                                       na.action = na.omit)
#save model output to docx file
modelsummary(models=splag_drought_smi_dist_knn,
             output="Output/splag_drought_smi_dist_knn.docx",
             fmt = 3,
             estimate = "{estimate} [{conf.low}, {conf.high}]",
             statistic = "{p.value}")

#alternative indicators only with full inverse distances
#DIM since 2014
splag_drought_dim14_dist_all <- lagsarlm(drought ~ DIM2014 + age + gender + edu + income_ + member + donation + no_ag,
                                       data = ce_data,
                                       listw = weights_all,
                                       zero.policy = T, 
                                       na.action = na.omit)
#save model output to docx file
modelsummary(models=splag_drought_dim14_dist_all,
             output="Output/splag_drought_dim14_dist_all.docx",
             fmt = 3,
             estimate = "{estimate} [{conf.low}, {conf.high}]",
             statistic = "{p.value}")

#DIM since 2007
splag_drought_dim07_dist_all <- lagsarlm(drought ~ DIM2007 + age + gender + edu + income_ + member + donation + no_ag,
                                         data = ce_data,
                                         listw = weights_all,
                                         zero.policy = T, 
                                         na.action = na.omit)
#save model output to docx file
modelsummary(models=splag_drought_dim07_dist_all,
             output="Output/splag_drought_dim07_dist_all.docx",
             fmt = 3,
             estimate = "{estimate} [{conf.low}, {conf.high}]",
             statistic = "{p.value}")

#DIM since 2000
splag_drought_dim00_dist_all <- lagsarlm(drought ~ DIM2000 + age + gender + edu + income_ + member + donation + no_ag,
                                         data = ce_data,
                                         listw = weights_all,
                                         zero.policy = T, 
                                         na.action = na.omit)
#save model output to docx file
modelsummary(models=splag_drought_dim00_dist_all,
             output="Output/splag_drought_dim00_dist_all.docx",
             fmt = 3,
             estimate = "{estimate} [{conf.low}, {conf.high}]",
             statistic = "{p.value}")

######################
# 3.2 Hypothesis 2 Flood
######################

#continuous indicator based on HOWAS database
splag_flood_dist_all <- lagsarlm(flood ~ howas + age + gender + edu + income_ + member + donation + no_ag,
                                 data = ce_data,
                                 listw = weights_all,
                                 zero.policy = T, 
                                 na.action = na.omit)
#save model output to docx file
modelsummary(models=splag_flood_dist_all,
             output="Output/splag_flood_dist_all.docx",
             fmt = 3,
             estimate = "{estimate} [{conf.low}, {conf.high}]",
             statistic = "{p.value}")

#binary indicator based on HOWAS database
splag_flood_binary_dist_all <- lagsarlm(flood ~ flooded + age + gender + edu + income_ + member + donation + no_ag,
                                        data = ce_data,
                                        listw = weights_all,
                                        zero.policy = T, 
                                        na.action = na.omit)
#save model output to docx file
modelsummary(models=splag_flood_binary_dist_all,
             output="Output/splag_flood_binary_dist_all.docx",
             fmt = 3,
             estimate = "{estimate} [{conf.low}, {conf.high}]",
             statistic = "{p.value}")

######################
# 3.3 Hypothesis 3 Climate
######################

#combine SMI and HOWAS as indicators of climate change impacts
splag_climate_dist_all <- lagsarlm(climate ~ SMI_Min + howas + age + gender + edu + income_ + member + donation + no_ag,
                                   data = ce_data,
                                   listw = weights_all,
                                   zero.policy = T,
                                   na.action = na.omit)
#save model output to docx file
modelsummary(models=splag_climate_dist_all,
             output="Output/splag_climate_dist_all.docx",
             fmt = 3,
             estimate = "{estimate} [{conf.low}, {conf.high}]",
             statistic = "{p.value}")

######################
# 3.4 Hypothesis 4 Water
######################

#binary indicator based on Rote Gebiete
splag_water_dist_all <- lagsarlm(water ~ Nitrate_RG + age + gender + edu + income_ + member + donation + no_ag,
                                 data = ce_data,
                                 listw = weights_all,
                                 zero.policy = T, 
                                 na.action = na.omit)
#save model output to docx file
modelsummary(models=splag_water_dist_all,
             output="Output/splag_water_dist_all.docx",
             fmt = 3,
             estimate = "{estimate} [{conf.low}, {conf.high}]",
             statistic = "{p.value}")

######################
# 3.5 Hypothesis 5 Opportunity costs
######################

#SQR with 33 km buffer
splag_opp_33km_dist_all <- lagsarlm(wtp_all ~ SQR_33km + age + gender + edu + income_ + member + donation + no_ag,
                                    data = ce_data,
                                    listw = weights_all,
                                    zero.policy = T, 
                                    na.action = na.omit)
#save model output to docx file
modelsummary(models=splag_opp_33km_dist_all,
             output="Output/splag_opp_33km_dist_all.docx",
             fmt = 3,
             estimate = "{estimate} [{conf.low}, {conf.high}]",
             statistic = "{p.value}")

######### Those would only work for subsets of the dataset (due to NAs in the SQR indicators), which doesn't seem to make sense
#SQR only with overlap
splag_opp_dist_all <- lagsarlm(wtp_all ~ SQR + age + gender + edu + income_ + member + donation + no_ag,
                               data = ce_data,
                               listw = weights_sqr,
                               zero.policy = T, 
                               na.action = na.omit)

#SQR with 1-km buffer
splag_opp_1km_dist_all <- lagsarlm(wtp_all ~ SQR_1km + age + gender + edu + income_ + member + donation + no_ag,
                                   data = ce_data,
                                   listw = weights_all,
                                   zero.policy = T, 
                                   na.action = na.omit)

#SQR with 10-km buffer
splag_opp_10km_dist_all <- lagsarlm(wtp_all ~ SQR_10km + age + gender + edu + income_ + member + donation + no_ag,
                                    data = ce_data,
                                    listw = weights_all,
                                    zero.policy = T, 
                                    na.action = na.omit)

#SQR with 20 km buffer
splag_opp_20km_dist_all <- lagsarlm(wtp_all ~ SQR_20km + age + gender + edu + income_ + member + donation + no_ag,
                                    data = ce_data,
                                    listw = weights_all,
                                    zero.policy = T, 
                                    na.action = na.omit)

######################
# 4. Spatial error models
######################

######################
# 4.1 Hypothesis 1 Drought
######################

#SMI_Min indicator, inverse distances, full matrix
sperr_drought_smi_dist_all <- errorsarlm(drought ~ SMI_Min + age + gender + edu + income_ + member + donation + no_ag,
                                        data = ce_data,
                                        listw = weights_all,
                                        zero.policy = T, 
                                        na.action = na.omit)

######################
# 4.2 Hypothesis 2 Flood
######################

#continuous indicator based on HOWAS database
sperr_flood_dist_all <- errorsarlm(flood ~ howas + age + gender + edu + income_ + member + donation + no_ag,
                                  data = ce_data,
                                  listw = weights_all,
                                  zero.policy = T, 
                                  na.action = na.omit)

######################
# 4.3 Hypothesis 3 Climate
######################

#combine SMI and HOWAS as indicators of climate change impacts
sperr_climate_dist_all <- errorsarlm(climate ~ SMI_Min + howas + age + gender + edu + income_ + member + donation + no_ag,
                                    data = ce_data,
                                    listw = weights_all,
                                    zero.policy = T,
                                    na.action = na.omit)

######################
# 4.4 Hypothesis 4 Water
######################

#binary indicator based on Rote Gebiete
sperr_water_dist_all <- errorsarlm(water ~ Nitrate_RG + age + gender + edu + income_ + member + donation + no_ag,
                                  data = ce_data,
                                  listw = weights_all,
                                  zero.policy = T, 
                                  na.action = na.omit)

######################
# 4.5 Hypothesis 5 Opportunity costs
######################

#SQR with 33 km buffer
sperr_opp_33km_dist_all <- errorsarlm(wtp_all ~ SQR_33km + age + gender + edu + income_ + member + donation + no_ag,
                                      data = ce_data,
                                      listw = weights_all,
                                      zero.policy = T, 
                                      na.action = na.omit)

######################
# 5. Spatial Durbin models
######################

######################
# 5.1 Hypothesis 1 Drought
######################

#SMI_Min indicator, inverse distances, full matrix
spdurb_drought_smi_dist_all <- lagsarlm(drought ~ SMI_Min + age + gender + edu + income_ + member + donation + no_ag,
                                        data = ce_data,
                                        listw = weights_all,
                                        zero.policy = T, 
                                        na.action = na.omit,
                                        Durbin = T)

######################
# 5.2 Hypothesis 2 Flood
######################

#continuous indicator based on HOWAS database
spdurb_flood_dist_all <- lagsarlm(flood ~ howas + age + gender + edu + income_ + member + donation + no_ag,
                                  data = ce_data,
                                  listw = weights_all,
                                  zero.policy = T, 
                                  na.action = na.omit,
                                  Durbin = T)

######################
# 5.3 Hypothesis 3 Climate
######################

#combine SMI and HOWAS as indicators of climate change impacts
spdurb_climate_dist_all <- lagsarlm(climate ~ SMI_Min + howas + age + gender + edu + income_ + member + donation + no_ag,
                                    data = ce_data,
                                    listw = weights_all,
                                    zero.policy = T,
                                    na.action = na.omit,
                                    Durbin = T)

######################
# 5.4 Hypothesis 4 Water
######################

#binary indicator based on Rote Gebiete
spdurb_water_dist_all <- lagsarlm(water ~ Nitrate_RG + age + gender + edu + income_ + member + donation + no_ag,
                                  data = ce_data,
                                  listw = weights_all,
                                  zero.policy = T, 
                                  na.action = na.omit,
                                  Durbin = T)

######################
# 5.5 Hypothesis 5 Opportunity costs
######################

#SQR with 33 km buffer
spdurb_opp_33km_dist_all <- lagsarlm(wtp_all ~ SQR_33km + age + gender + edu + income_ + member + donation + no_ag,
                                     data = ce_data,
                                     listw = weights_all,
                                     zero.policy = T, 
                                     na.action = na.omit,
                                     Durbin = T)