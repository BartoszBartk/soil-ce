##### Clean up and prepare dataset for spatial analyses

require(here)

#load dataset prepared by Lukas
ce_data <- read.csv(here("Spatial data/04 Hypothesis Evaluation/04_06 Results/ce-results.csv"),header = T)
#add modified variables
#make educational achievement binary
ce_data$edu <- ifelse(ce_data$edu == "abitur" | ce_data$edu == "university",1,0)
#income: take middle points of brackets as proxies for continuous income (for highest bracket, ca. double the mean net income in Germany [3580 according to Destatis])
ce_data$income_ <- c(500,1250,1750,2250,3000,4250,7000,ce_data$income_)[match(ce_data$income,c("below_1000","1000_1500","1500_2000","2000_2500","2500_3500","3500_5000","5000_above",ce_data$income))]
#sum of WTP
ce_data$wtp_all <- ce_data$drought + ce_data$flood + ce_data$climate + ce_data$water
#binary flood variable
ce_data$flooded <- ifelse(ce_data$howas > 0, TRUE, FALSE)
#change some variable names
ce_data <- rename(ce_data, howas = Insg_Howas,
                  DIM2014 = X2014.20__1,
                  DIM2007 = X2007.20__1,
                  DIM2000 = X2000.20__1,
                  xcoor = x.coordina,
                  ycoor = y.coordina)
#remove irrelevant variables
ce_data <- subset(ce_data, select = c("ID","gender","age","plz","urban","edu","member","donation","no_ag","income_","drought","flood","climate","water","wtp_all","SQR","SQR_1km","SQR_10km","SQR_20km","SQR_33km","howas","DIM2014","DIM2007","DIM2000","SMI_Min","Nitrate_RG","xcoor","ycoor"))
#force DIM variables to numeric
ce_data$DIM2014 <- as.numeric(gsub(",", "\\.", ce_data$DIM2014))
ce_data$DIM2007 <- as.numeric(gsub(",", "\\.", ce_data$DIM2007))
ce_data$DIM2000 <- as.numeric(gsub(",", "\\.", ce_data$DIM2000))
#recode Nitrate_RG
ce_data$Nitrate_RG <- ifelse(ce_data$Nitrate_RG == "schlecht",TRUE,FALSE)
#write out the final file
write.csv(ce_data, here("Spatial data/ce_results_spatial_data_clean.csv"), row.names = F, fileEncoding = "UTF-8")
