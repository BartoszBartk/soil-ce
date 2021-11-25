############################
# Analysis of the data from the BonaRes soil ecosystem services CE
# n=1500, conducted online between 21 June and 6 July 2021
# authors: Bartosz Bartkowski, Julian Massenberg, Nele Lienhoop
# contact: bartosz.bartkowski@ufz.de
# This script can be found at https://github.com/BartoszBartk/soil-ce 
############################

#manually install older version of mlogit due to issues with newer version ('gmnl' requires mlogit.data(), which has been deprecated by 'mlogit')
install.packages("https://cran.r-project.org/src/contrib/Archive/mlogit/mlogit_1.0-2.tar.gz", repos=NULL,type="source")
require(mlogit)
require(dplyr)
require(gmnl)
require(support.CEs)
require(modelsummary)

#tidy and glance for modelsummary(), thanks to Maksym Polyakov
source("tidy.gmnl.R") #can be downloaded from Maksym's GitHub rep https://github.com/polyama/tidy.gmnl/tree/main/R
gm <- tibble::tribble(
  ~raw,        ~clean,          ~fmt,
  "nind",     "Num.Ind",       0,
  "nobs",      "Num.Obs",        0,
  "AIC",       "AIC",            1,
  "BIC",       "BIC",            1,
  "logLik",    "Log.Lik",        3
)

#read in dataset already formatted for clogit/mlogit analysis
#the original dataset is available at https://maps.bonares.de/mapapps/resources/apps/bonares/index.html?lang=en&mid=cc787645-47d2-41a8-af10-b699a64d7b3b, citable as Bartkowski et al. (2021): German public's preferences for soil-based ecosystem services (discrete choice experiment). DOI:10.13039/501100002347
#the dataset used here has already been cleaned by removing protest votes (see script cleaning_protest_votes.R in the GitHub repository)
ce_main_dataset <- read.csv("Data/BonaRes_CE_results_clean.csv",header=T)

#create dataset for mixed logit
data_gmnl <- mlogit.data(ce_main_dataset,
                         choice="RES",
                         shape="long",
                         id.var="ID",
                         group.var="BLOCK",
                         alt.var="ALT",
                         varying = 7:11) #starting with 'drought'

#####################################
#run simple conditional logit
clogit <- gmnl(RES ~ drought + flood + climate + water + price | 1,
               data=data_gmnl,
               reflevel=3)
summary(clogit)

#save results
modelsummary(models=clogit,
             output="Output/clogit.docx",
             estimate="{estimate} ({std.error}){stars}",
             gof_map = gm,
             statistic=NULL)
#with exact p-values & confidence intervals
modelsummary(models=clogit,
             output="Output/clogit2.docx",
             estimate="{estimate}",
             gof_map = gm,
             statistic="{p.value} [{conf.low}, {conf.high}]")

#######################################
#create dataset for mixed logit with opposite of price
data_gmnl2 <- mlogit.data(ce_main_dataset,
                         choice="RES",
                         shape="long",
                         id.var="ID",
                         group.var="BLOCK",
                         alt.var="ALT",
                         varying = 7:11,
                         opposite=c("price"))

##run simple mixed logit
#set starting values for mixed logit using coefficient estimates from clogit except for 'price' because of the assumed lognormal distribution
start <- c(coef(clogit)[1:6],-3,rep(1,4),2)

#run model
mxl1 <- gmnl(RES ~ drought + flood + climate + water + price | 1,
             data=data_gmnl2,
             model="mixl",
             panel=T,
             ranp=c(drought="n",flood="n",climate="n",water="n",price="ln"),
             R=1000,
             haltons=NA,
             start=start,
             reflevel=3,
             method="bhhh")
summary(mxl1)

#save results
modelsummary(models=mxl1,
             output="Output/mxl1.docx",
             estimate="{estimate} ({std.error}){stars}",
             gof_map = gm,
             statistic=NULL)
#with exact p-values & confidence intervals
modelsummary(models=mxl1,
             output="Output/mxl1_2.docx",
             estimate="{estimate}",
             gof_map = gm,
             statistic="{p.value} [{conf.low}, {conf.high}]")

##calculate marginal WTP estimates
#transform the price coefficient to be on the same scale as the ES coefficients
mxl1_transformed <- mxl1
mu <- mxl1_transformed$coefficients[7]
sigma <- mxl_transformed$coefficients[12]
mxl1_transformed$coefficients[7] <- -exp(mu + (sigma^2)/2)

#calculate mean mWTP (option 1, using 'support.CEs' [WTPs + CIs])
wtp_mxl1_mean <- mwtp(mxl1_transformed,
                      monetary.variables = "price",
                      method = "delta")
wtp_mxl1_mean
#option 2, using 'gmnl' [WTPs + SEs]
wtp_mxl1_mean2 <- wtp.gmnl(mxl1_transformed,
                           wrt="price")

#calculate median mWTP
mxl1_transformed$coefficients[7] <- -exp(mu)
wtp_mxl1_median <- mwtp(mxl1_transformed,
                        monetary.variables = "price",
                        method = "delta")
wtp_mxl1_median

############################
#reestimate model in WTP space to check WTP estimates
#set starting values for S-MNL model (procedure largely based on section 3.5 by Sarrias and Daziano (2017))
stval <- c(0, 0, 0.1, 0, 0, 0, 0, 0, 0) #using 1 instead of 0.1 leads to a singular matrix

#estimate S-MNL model in WTP space
smnl_wtps <- gmnl(RES ~ price + drought + flood + climate + water | 1 | 0 | 0 | 1,
                  data=data_gmnl2,
                  model="smnl",
                  R=1,
                  fixed=c(F,F,T,F,F,F,F,T,F),
                  panel=T,
                  start=stval,
                  method="bhhh",
                  reflevel=3,
                  interlim=500)
summary(smnl_wtps)

#now the same in a G-MNL model, using S-MNL parameters for starting values
stval2 <- c(coef(smnl_wtps)[1:2],-3,coef(smnl_wtps)[3:7],2,rep(0.1,5),0)

#estimate generalized multinomial logit model in WTP space
mxl_wtps <- gmnl(RES ~ price + drought + flood + climate + water | 1 | 0 | 0 | 1,
                 data=data_gmnl2,
                 model="gmnl",
                 R=1000,
                 fixed=c(F,F,T,rep(F,11),T), 
                 panel=T,
                 start=stval2,
                 method="bhhh",
                 ranp=c(price="ln",drought="n",flood="n",climate="n",water="n"),
                 reflevel=3)
summary(mxl_wtps)

#save results
modelsummary(models=mxl_wtps,
             output="Output/mxl_wtps.docx",
             estimate="{estimate} ({std.error}){stars}",
             gof_map = gm,
             statistic=NULL)
#with exact p-values & confidence intervals
modelsummary(models=mxl_wtps,
             output="Output/mxl_wtps_2.docx",
             estimate="{estimate}",
             gof_map = gm,
             statistic="{p.value} [{conf.low}, {conf.high}]")

############################
#run mixed logit with explanatory variables for mean of random parameters
#set starting values for mixed logit, including for 22 interactions
start2 <- c(coef(clogit)[1:6],-3,rep(0,22),rep(1,4),2)

#run model
mxl2 <- gmnl(RES ~ drought + flood + climate + water + price | 1 | 0 | exp_drought + exp_flood + awareness + knowledge + urban + no_ag,
                  data=data_gmnl2,
                  model="mixl",
                  panel=T,
                  ranp=c(drought="n",flood="n",climate="n",water="n",price="ln"),
                  mvar=list(drought=c("exp_drought","awareness","knowledge","urban","no_ag"),flood=c("exp_flood","awareness","knowledge","urban","no_ag"),climate=c("awareness","knowledge","urban","no_ag"),water=c("awareness","knowledge","urban","no_ag"),price=c("awareness","knowledge","urban","no_ag")),
                  R=1000,
                  haltons=NA,
                  start=start2,
                  reflevel=3,
                  method="bhhh")
summary(mxl2)

#save results
modelsummary(models=mxl2,
             output="Output/mxl2.docx",
             estimate="{estimate} ({std.error}){stars}",
             gof_map = gm,
             statistic=NULL)
#with exact p-values & confidence intervals
modelsummary(models=mxl2,
             output="Output/mxl2_2.docx",
             estimate="{estimate}",
             gof_map = gm,
             statistic="{p.value} [{conf.low}, {conf.high}]")

#################################
###run mixed logit with exaplanatory variables for the choice of status-quo alternatives
##simplify the variables edu & income
#edu: create binary variable with Abitur and above as 1, below Abitur as 0
data_gmnl2$abi <- ifelse(data_gmnl2$edu=="abitur" | data_gmnl2$edu=="university",1,0)
#income: take middle points of brackets as proxies for continuous income (for highest bracket, ca. double the mean net income in Germany [3580 according to Destatis])
data_gmnl2$income_ <- c(500,1250,1750,2250,3000,4250,7000,data_gmnl2$income_)[match(data_gmnl2$income,c("below_1000","1000_1500","1500_2000","2000_2500","2500_3500","3500_5000","5000_above",data_gmnl2$income))]

#set starting values for mixed logit, including for 10 interactions
start3 <- c(coef(clogit)[1:6],-3,rep(0,20),rep(1,4),2)

#run model
mxl3 <- gmnl(RES ~ drought + flood + climate + water + price | gender + age + urban + abi + income_ + awareness + knowledge + no_ag + donation + member,
                  data=data_gmnl2,
                  model="mixl",
                  panel=T,
                  ranp=c(drought="n",flood="n",climate="n",water="n",price="ln"),
                  R=1000,
                  haltons=NA,
                  start=start3,
                  method="bhhh")
summary(mxl3)

#save results
modelsummary(models=mxl3,
             output="Output/mxl3.docx",
             estimate="{estimate} ({std.error}){stars}",
             gof_map = gm,
             statistic=NULL)
#with exact p-values & confidence intervals
modelsummary(models=mxl3_gmnl,
             output="Output/mxl3_2.docx",
             estimate="{estimate}",
             gof_map = gm,
             statistic="{p.value} [{conf.low}, {conf.high}]")

################################
##rescale ES attributes and rerun the basic mixed logit
#rescale ES attributes
data_gmnl3 <- data_gmnl2
data_gmnl3$drought <- as.numeric(data_gmnl2$drought/100)
data_gmnl3$flood <- as.numeric(data_gmnl2$flood/100)
data_gmnl3$climate <- as.numeric(data_gmnl2$climate/100)
data_gmnl3$water <- as.numeric(data_gmnl2$water/100)

start4 <- c(coef(clogit)[1:6]*100,coef(clogit)[7],rep(1,5))

#run simple mixed logit
mxl4 <- gmnl(RES ~ drought + flood + climate + water + price | 1,
                  data=data_gmnl3,
                  model="mixl",
                  panel=T,
                  ranp=c(drought="n",flood="n",climate="n",water="n",price="ln"),
                  R=1000,
                  haltons=NA,
                  start=start4,
                  reflevel=3,
                  method="bhhh")
summary(mxl4)

#save results
modelsummary(models=mxl4,
             output="Output/mxl4.docx",
             estimate="{estimate} ({std.error})",
             gof_map = gm,
             statistic="{p.value} [{conf.low}, {conf.high}]")

###############################
##run mixed logit with fixed price parameter
#run model
mxl1_alt <- gmnl(RES ~ price + drought + flood + climate + water | 1,
                      data=data_gmnl2,
                      model="mixl",
                      panel=T,
                      ranp=c(drought="n",flood="n",climate="n",water="n"),
                      R=1000,
                      haltons=NA,
                      reflevel=3,
                      method="bhhh")
summary(mxl1_alt)

#save results
modelsummary(models=mxl1_alt,
             output="Output/mxl1_alt.docx",
             estimate="{estimate} ({std.error})",
             gof_map = gm,
             statistic="{p.value} [{conf.low}, {conf.high}]")
