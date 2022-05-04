############################
# Analysis of the data from the BonaRes soil ecosystem services CE
# n=1500, conducted online between 21 June and 6 July 2021
# authors: Bartosz Bartkowski, Julian Massenberg, Nele Lienhoop
# contact: bartosz.bartkowski@ufz.de
# This script can be found at https://github.com/BartoszBartk/soil-ce 
############################

#optional: clear environment
rm(list=ls())

#load packages
require(apollo)
require(here)

#read in dataset already formatted for analysis
#the original dataset is available at https://maps.bonares.de/mapapps/resources/apps/bonares/index.html?lang=en&mid=cc787645-47d2-41a8-af10-b699a64d7b3b, citable as Bartkowski et al. (2021): German public's preferences for soil-based ecosystem services (discrete choice experiment). DOI:10.13039/501100002347
#the dataset used here has already been cleaned by removing protest votes (see script cleaning_protest_votes.R in the GitHub repository)
ce_main_dataset <- read.csv(paste0(here(),"/Data/BonaRes_CE_results_clean.csv"),header=T)

#initialize apollo package
apollo_initialise()

#create database by reshaping dataset to wide format as required by apollo
database <- reshape(ce_main_dataset,
                    idvar="STR",
                    timevar="ALT",
                    v.names=c("drought","flood","climate","water","price","RES"),
                    drop=c("X","ASC"),
                    direction="wide")

#add choice variable indicating which alternative was chosen
database$choice <- ifelse(
  database$RES.1==TRUE,1,ifelse(
    database$RES.2==TRUE,2,3
  ))

#add modified variables
#abi (binary variable of educational achievement)
database$abi <- ifelse(database$edu=="abitur" | database$edu=="university",1,0)
#income: take middle points of brackets as proxies for continuous income (for highest bracket, ca. double the mean net income in Germany [3580 according to Destatis])
database$income_ <- c(500,1250,1750,2250,3000,4250,7000,database$income_)[match(database$income,c("below_1000","1000_1500","1500_2000","2000_2500","2500_3500","3500_5000","5000_above",database$income))]

#remove original dataset
rm(ce_main_dataset)

##################################### multinomial logit
##implement a simple MNL
#set required parameters of project
apollo_control <- list(
  modelName="soil_CE_mnl",
  modelDescr="Discrete choice experiment into preferences for soil-based ecosystem services",
  indivID="ID"
)

#define coefficients to be estimated and starting values for each
apollo_beta <- c(asc_1=0,
                 asc_2=0,
                 asc_3=0,
                 b_drought=0,
                 b_flood=0,
                 b_climate=0,
                 b_water=0,
                 b_price=0)

#set fixed 0 coefficients (here: only one of the alternative-specific constants)
apollo_fixed <- "asc_1"

#validate inputs
apollo_inputs <- apollo_validateInputs()

#define model (for details of the required function elements, see apollo documentation)
apollo_probabilities <- function(apollo_beta,apollo_inputs,functionality="estimate"){
  apollo_attach(apollo_beta,apollo_inputs)
  on.exit(apollo_detach(apollo_beta,apollo_inputs))
  P = list()
  V = list()
  #define utilities for each alternative
  V[['alt1']] = asc_1 + b_drought * drought.1 + b_flood * flood.1 + b_climate * climate.1 + b_water * water.1 + b_price * price.1
  V[['alt2']] = asc_2 + b_drought * drought.2 + b_flood * flood.2 + b_climate * climate.2 + b_water * water.2 + b_price * price.2
  V[['alt3']] = asc_3 + b_drought * drought.3 + b_flood * flood.3 + b_climate * climate.3 + b_water * water.3 + b_price * price.3
  mnl_settings = list(
    alternatives = c(alt1=1,alt2=2,alt3=3),
    choiceVar = choice,
    V = V
  ) #model components
  P[["model"]] = apollo_mnl(mnl_settings,functionality)
  P = apollo_panelProd(P,apollo_inputs,functionality)
  P = apollo_prepareProb(P,apollo_inputs,functionality)
  return(P)
}

#estimate model
mnl <- apollo_estimate(apollo_beta,apollo_fixed,apollo_probabilities,apollo_inputs)

#show outputs
apollo_modelOutput(mnl,
                   modelOutput_settings=list(printPVal=1))

#write outputs
apollo_saveOutput(mnl,
                  saveOutput_settings=list(printPVal=1))

############################## mixed logit
##implement a simple mixed logit model to allow for inter-individual preference heterogeneity
apollo_control = list(
  modelName = "soil_CE_mxl",
  modelDescr = "Discrete choice experiment into preferences for soil-based ecosystem services",
  indivID = "ID",
  mixing = TRUE,
  nCores = 3,
  seed=2104
)

#define starting values (using previous results from analyses with gmnl package)
apollo_beta <- c(
  asc_1 = 0,
  asc_2 = 0,
  asc_3 = 0,
  mu_b_drought = 0.01,
  mu_b_flood = 0.01,
  mu_b_climate = 0.01,
  mu_b_water = 0.01,
  mu_log_b_price = -1,
  sigma_b_drought = 0.01,
  sigma_b_flood = 0.01,
  sigma_b_climate = 0.01,
  sigma_b_water = 0.01,
  sigma_log_b_price = 1.5
)

#set fixed 0 coefficients (here: only one of the alternative-specific constants)
apollo_fixed <- "asc_1"

#define parameters for the simulation
apollo_draws = list(
  interDrawsType = "sobol",
  interNDraws = 1000,
  interNormDraws = c("draws_price_inter","draws_drought_inter","draws_flood_inter","draws_climate_inter","draws_water_inter")
)

#define random coefficients
apollo_randCoeff = function(apollo_beta,apollo_inputs){
  randcoeff = list()
  randcoeff[["b_drought"]] = mu_b_drought + sigma_b_drought * draws_drought_inter
  randcoeff[["b_flood"]] = mu_b_flood + sigma_b_flood * draws_flood_inter
  randcoeff[["b_climate"]] = mu_b_climate + sigma_b_climate * draws_climate_inter
  randcoeff[["b_water"]] = mu_b_water + sigma_b_water * draws_water_inter
  randcoeff[["b_price"]] = -exp(mu_log_b_price + sigma_log_b_price * draws_price_inter)
  return(randcoeff)
}

#validate inputs
apollo_inputs <- apollo_validateInputs()

#define model (for details, see apollo documentation)
apollo_probabilities = function(apollo_beta,apollo_inputs,functionality="estimate"){
  apollo_attach(apollo_beta,apollo_inputs)
  on.exit(apollo_detach(apollo_beta,apollo_inputs))
  P = list()
  V = list()
  V[['alt1']] = asc_1 + b_drought * drought.1 + b_flood * flood.1 + b_climate * climate.1 + b_water * water.1 + b_price * price.1
  V[['alt2']] = asc_2 + b_drought * drought.2 + b_flood * flood.2 + b_climate * climate.2 + b_water * water.2 + b_price * price.2
  V[['alt3']] = asc_3 + b_drought * drought.3 + b_flood * flood.3 + b_climate * climate.3 + b_water * water.3 + b_price * price.3
  mnl_settings = list(
    alternatives = c(alt1=1,alt2=2,alt3=3),
    choiceVar = choice,
    V = V
  ) #model components
  P[["model"]] = apollo_mnl(mnl_settings,functionality)
  P = apollo_panelProd(P,apollo_inputs,functionality)
  P = apollo_avgInterDraws(P,apollo_inputs,functionality)
  P = apollo_prepareProb(P,apollo_inputs,functionality)
  return(P)
}

#estimate model
mxl <- apollo_estimate(apollo_beta,apollo_fixed,apollo_probabilities,apollo_inputs)

#show output
apollo_modelOutput(mxl,
                   modelOutput_settings=list(printPVal=1))
#write output
apollo_saveOutput(mxl,
                  saveOutput_settings=list(printPVal=1))

############################## mixed logit in WTP space
##implement the same simple mixed logit, but in WTP space
apollo_control = list(
  modelName = "soil_CE_mxl_wtp",
  modelDescr = "Discrete choice experiment into preferences for soil-based ecosystem services",
  indivID = "ID",
  mixing = TRUE,
  nCores = 3,
  seed=2104
)

#define starting values (using previous results from analyses with gmnl package)
apollo_beta <- c(
  asc_1 = 0,
  asc_2 = 0,
  asc_3 = 0,
  mu_b_drought = 0.01,
  mu_b_flood = 0.01,
  mu_b_climate = 0.01,
  mu_b_water = 0.01,
  mu_log_b_price = -1,
  sigma_b_drought = 0.01,
  sigma_b_flood = 0.01,
  sigma_b_climate = 0.01,
  sigma_b_water = 0.01,
  sigma_log_b_price = 1.5
)

#set fixed 0 coefficients (here: only one of the alternative-specific constants)
apollo_fixed = c("asc_1")

#define parameters for the simulation
apollo_draws = list(
  interDrawsType = "sobol",
  interNDraws = 1000,
  interNormDraws = c("draws_price_inter","draws_drought_inter","draws_flood_inter","draws_climate_inter","draws_water_inter")
)

#define random coefficients
apollo_randCoeff = function(apollo_beta,apollo_inputs){
  randcoeff = list()
  randcoeff[["b_drought"]] = mu_b_drought + sigma_b_drought * draws_drought_inter
  randcoeff[["b_flood"]] = mu_b_flood + sigma_b_flood * draws_flood_inter
  randcoeff[["b_climate"]] = mu_b_climate + sigma_b_climate * draws_climate_inter
  randcoeff[["b_water"]] = mu_b_water + sigma_b_water * draws_water_inter
  randcoeff[["b_price"]] = -exp(mu_log_b_price + sigma_log_b_price * draws_price_inter)
  return(randcoeff)
}

#validate inputs
apollo_inputs <- apollo_validateInputs()

#define model
apollo_probabilities = function(apollo_beta,apollo_inputs,functionality="estimate"){
  apollo_attach(apollo_beta,apollo_inputs)
  on.exit(apollo_detach(apollo_beta,apollo_inputs))
  P = list()
  V = list()
  #multiply all coefficients by the price coefficient
  V[['alt1']] = asc_1 + b_price * (b_drought * drought.1 + b_flood * flood.1 + b_climate * climate.1 + b_water * water.1 + price.1)
  V[['alt2']] = asc_2 + b_price * (b_drought * drought.2 + b_flood * flood.2 + b_climate * climate.2 + b_water * water.2 + price.2)
  V[['alt3']] = asc_3 + b_price * (b_drought * drought.3 + b_flood * flood.3 + b_climate * climate.3 + b_water * water.3 + price.3)
  mnl_settings = list(
    alternatives = c(alt1=1,alt2=2,alt3=3),
    choiceVar = choice,
    V = V
  ) #model components
  P[["model"]] = apollo_mnl(mnl_settings,functionality)
  P = apollo_panelProd(P,apollo_inputs,functionality)
  P = apollo_avgInterDraws(P,apollo_inputs,functionality)
  P = apollo_prepareProb(P,apollo_inputs,functionality)
  return(P)
}

#estimate model
mxl_wtp <- apollo_estimate(apollo_beta,apollo_fixed,apollo_probabilities,apollo_inputs)

#show output
apollo_modelOutput(mxl_wtp,
                   modelOutput_settings=list(printPVal=1))
#write output
apollo_saveOutput(mxl_wtp,
                  saveOutput_settings=list(printPVal=1))

#compute conditional individual estimates (to be used later for spatial analyses)
conditionals <- apollo_conditionals(mxl_wtp,apollo_probabilities,apollo_inputs)
#combine in a dataframe
condits <- cbind(conditionals$b_drought,conditionals$b_flood,conditionals$b_climate,conditionals$b_water)
#remove unnecessary ID columns
condits <- condits[,c(1,2,3,5,6,8,9,11,12)]
#rename columns
colnames(condits) <- c("ID","drought.m","drought.sd","flood.m","flood.sd","climate.m","climate.sd","water.m","water.sd")
#write
write.csv(condits,paste0(here(),"/conditional_wtps.csv"),row.names=F)

############################## mixed logit with individual-variable interactions to explain heterogeneity in ES preference estimates
apollo_control = list(
  modelName = "soil_CE_mxl_het",
  modelDescr = "Discrete choice experiment into preferences for soil-based ecosystem services",
  indivID = "ID",
  mixing = TRUE,
  nCores = 3,
  seed=2104
)

#define starting values (using previous results from analyses with gmnl package)
apollo_beta <- c(
  asc_1 = 0,
  asc_2 = 0,
  asc_3 = 0,
  mu_b_drought = 0.01,
  mu_b_flood = 0.01,
  mu_b_climate = 0.01,
  mu_b_water = 0.01,
  mu_log_b_price = -1,
  sigma_b_drought = 0.01,
  sigma_b_flood = 0.01,
  sigma_b_climate = 0.01,
  sigma_b_water = 0.01,
  sigma_log_b_price = 1.5,
  b_exp_drought = 0,
  b_aware_drought = 0,
  b_know_drought = 0,
  b_urban_drought = 0,
  b_ag_drought = 0,
  b_exp_flood = 0,
  b_aware_flood = 0,
  b_know_flood = 0,
  b_urban_flood = 0,
  b_ag_flood = 0,
  b_aware_clim = 0,
  b_know_clim = 0,
  b_urban_clim = 0,
  b_ag_clim = 0,
  b_aware_water = 0,
  b_know_water = 0,
  b_urban_water = 0,
  b_ag_water = 0,
  b_aware_price = 0,
  b_know_price = 0,
  b_urban_price = 0,
  b_ag_price = 0
)

#set fixed 0 coefficients (here: only one of the alternative-specific constants)
apollo_fixed <- "asc_1"

#define parameters for the simulation
apollo_draws = list(
  interDrawsType = "sobol",
  interNDraws = 1000,
  interNormDraws = c("draws_price_inter","draws_drought_inter","draws_flood_inter","draws_climate_inter","draws_water_inter")
)

#define random coefficients
apollo_randCoeff = function(apollo_beta,apollo_inputs){
  randcoeff = list()
  randcoeff[["b_drought"]] = mu_b_drought + sigma_b_drought * draws_drought_inter
  randcoeff[["b_flood"]] = mu_b_flood + sigma_b_flood * draws_flood_inter
  randcoeff[["b_climate"]] = mu_b_climate + sigma_b_climate * draws_climate_inter
  randcoeff[["b_water"]] = mu_b_water + sigma_b_water * draws_water_inter
  randcoeff[["b_price"]] = -exp(mu_log_b_price + sigma_log_b_price * draws_price_inter)
  return(randcoeff)
}

#validate inputs
apollo_inputs <- apollo_validateInputs()

#define model (for details, see apollo documentation)
apollo_probabilities = function(apollo_beta,apollo_inputs,functionality="estimate"){
  apollo_attach(apollo_beta,apollo_inputs)
  on.exit(apollo_detach(apollo_beta,apollo_inputs))
  P = list()
  
  #define interaction terms
  b_drought_value = b_drought + b_exp_drought * exp_drought + b_aware_drought * awareness + b_know_drought * knowledge + b_urban_drought * urban + b_ag_drought * no_ag
  b_flood_value = b_flood + b_exp_flood * exp_flood + b_aware_flood * awareness + b_know_flood * knowledge + b_urban_flood * urban + b_ag_flood * no_ag
  b_climate_value = b_climate + b_aware_clim * awareness + b_know_clim * knowledge + b_urban_clim * urban + b_ag_clim * no_ag
  b_water_value = b_water + b_aware_water * awareness + b_know_water * knowledge + b_urban_water * urban + b_ag_water * no_ag
  b_price_value = b_price + b_aware_price * awareness + b_know_price * knowledge + b_urban_price * urban + b_ag_price * no_ag
  
  V = list()
  V[['alt1']] = asc_1 + b_drought_value * drought.1 + b_flood_value * flood.1 + b_climate_value * climate.1 + b_water_value * water.1 + b_price_value * price.1
  V[['alt2']] = asc_2 + b_drought_value * drought.2 + b_flood_value * flood.2 + b_climate_value * climate.2 + b_water_value * water.2 + b_price_value * price.2
  V[['alt3']] = asc_3 + b_drought_value * drought.3 + b_flood_value * flood.3 + b_climate_value * climate.3 + b_water_value * water.3 + b_price_value * price.3
  
  mnl_settings = list(
    alternatives = c(alt1=1,alt2=2,alt3=3),
    choiceVar = choice,
    V = V
  ) #model components
  P[["model"]] = apollo_mnl(mnl_settings,functionality)
  P = apollo_panelProd(P,apollo_inputs,functionality)
  P = apollo_avgInterDraws(P,apollo_inputs,functionality)
  P = apollo_prepareProb(P,apollo_inputs,functionality)
  return(P)
}

#estimate model
mxl_het <- apollo_estimate(apollo_beta,apollo_fixed,apollo_probabilities,apollo_inputs)

#show output
apollo_modelOutput(mxl_het,
                   modelOutput_settings=list(printPVal=1))
#write output
apollo_saveOutput(mxl_het,
                  saveOutput_settings=list(printPVal=1))

############################## mixed logit with SQ interactions
apollo_control = list(
  modelName = "soil_CE_mxl_sq",
  modelDescr = "Discrete choice experiment into preferences for soil-based ecosystem services",
  indivID = "ID",
  mixing = TRUE,
  nCores = 3,
  seed=2104
)

#define starting values (using previous results from analyses with gmnl package)
apollo_beta <- c(
  asc_1 = 0,
  asc_2 = 0,
  asc_3 = 0,
  asc_3_gender = 0,
  asc_3_age = 0,
  asc_3_urban = 0,
  asc_3_abi = 0,
  asc_3_income = 0,
  asc_3_awareness = 0,
  asc_3_knowledge = 0,
  asc_3_no_ag = 0,
  asc_3_donation = 0,
  asc_3_member = 0,
  mu_b_drought = 0.01,
  mu_b_flood = 0.01,
  mu_b_climate = 0.01,
  mu_b_water = 0.01,
  mu_log_b_price = -1,
  sigma_b_drought = 0.01,
  sigma_b_flood = 0.01,
  sigma_b_climate = 0.01,
  sigma_b_water = 0.01,
  sigma_log_b_price = 1.5
)

#set fixed 0 coefficients (here: only one of the alternative-specific constants; asc_1 rather than the SQ alternative in order to be able to later interact the SQ alternative with individual variables)
apollo_fixed <- "asc_1"

#define parameters for the simulation
apollo_draws = list(
  interDrawsType = "sobol",
  interNDraws = 1000,
  interNormDraws = c("draws_price_inter","draws_drought_inter","draws_flood_inter","draws_climate_inter","draws_water_inter")
)

#define random coefficients
apollo_randCoeff = function(apollo_beta,apollo_inputs){
  randcoeff = list()
  randcoeff[["b_drought"]] = mu_b_drought + sigma_b_drought * draws_drought_inter
  randcoeff[["b_flood"]] = mu_b_flood + sigma_b_flood * draws_flood_inter
  randcoeff[["b_climate"]] = mu_b_climate + sigma_b_climate * draws_climate_inter
  randcoeff[["b_water"]] = mu_b_water + sigma_b_water * draws_water_inter
  randcoeff[["b_price"]] = -exp(mu_log_b_price + sigma_log_b_price * draws_price_inter)
  return(randcoeff)
}

#validate inputs
apollo_inputs <- apollo_validateInputs()

#define model (for details, see apollo documentation)
apollo_probabilities = function(apollo_beta,apollo_inputs,functionality="estimate"){
  apollo_attach(apollo_beta,apollo_inputs)
  on.exit(apollo_detach(apollo_beta,apollo_inputs))
  P = list()
  
  #define SQ ASC interactions
  asc_3_value = asc_3 + asc_3_abi * abi + asc_3_age * age + asc_3_awareness * awareness + asc_3_donation * donation + asc_3_gender * gender + asc_3_income * income_ + asc_3_knowledge * knowledge + asc_3_member * member + asc_3_no_ag * no_ag + asc_3_urban * urban
  
  V = list()
  V[['alt1']] = asc_1 + b_drought * drought.1 + b_flood * flood.1 + b_climate * climate.1 + b_water * water.1 + b_price * price.1
  V[['alt2']] = asc_2 + b_drought * drought.2 + b_flood * flood.2 + b_climate * climate.2 + b_water * water.2 + b_price * price.2
  V[['alt3']] = asc_3_value + b_drought * drought.3 + b_flood * flood.3 + b_climate * climate.3 + b_water * water.3 + b_price * price.3
  
  mnl_settings = list(
    alternatives = c(alt1=1,alt2=2,alt3=3),
    choiceVar = choice,
    V = V
  ) #model components
  P[["model"]] = apollo_mnl(mnl_settings,functionality)
  P = apollo_panelProd(P,apollo_inputs,functionality)
  P = apollo_avgInterDraws(P,apollo_inputs,functionality)
  P = apollo_prepareProb(P,apollo_inputs,functionality)
  return(P)
}

#estimate model
mxl_sq <- apollo_estimate(apollo_beta,apollo_fixed,apollo_probabilities,apollo_inputs)

#show output
apollo_modelOutput(mxl_sq,
                   modelOutput_settings=list(printPVal=1))
#write output
apollo_saveOutput(mxl_sq,
                  saveOutput_settings=list(printPVal=1))