#analysis of pretest data of a CE (n=50, conducted between 10 and 14 June 2021)
require(survival)
require(support.CEs)
pretest <- read.csv("Data/UFZ_DiscreteChoiceExperiment_Pretest_N50.csv",header=T,sep=";")
summary(pretest)

##create a dataset readable for clogit() function
#extract choice responses
pretest_ce <- pretest[c(1,76:84)]
colnames(pretest_ce) <- c("ID","q1","q2","q3","q4","q5","q6","q7","q8","BLOCK")

#create a design matrix
design_pretest <- read.csv("Data/design_pretest.csv",skip=1,header=T,sep=";")
#add status quo alternative
design_pretest$sq.drought <- rep(70,8)
design_pretest$sq.flood <- rep(70,8)
design_pretest$sq.climate <- rep(50,8)
design_pretest$sq.water <- rep(30,8)
design_pretest$sq.price <- rep(0,8)
#extract alternatives 1 to 3, add variable ALT = 1, rename columns
design_pretest_alt1 <- design_pretest[,1:6]
design_pretest_alt1$ALT <- rep(1,8)
colnames(design_pretest_alt1) <- c("QES","drought","flood","climate","water","price","ALT")
design_pretest_alt2 <- design_pretest[,c(1,7:11)]
design_pretest_alt2$ALT <- rep(2,8)
colnames(design_pretest_alt2) <- c("QES","drought","flood","climate","water","price","ALT")
design_pretest_alt3 <- design_pretest[,c(1,12:16)]
design_pretest_alt3$ALT <- rep(3,8)
colnames(design_pretest_alt3) <- c("QES","drought","flood","climate","water","price","ALT")
#combine by rows
design_matrix <- rbind(design_pretest_alt1,design_pretest_alt2,design_pretest_alt3)
#add BLOCK variable
design_matrix$BLOCK <- rep(1,24)
#add ASC variable (0 for status quo alternative)
design_matrix$ASC <- ifelse(design_matrix$ALT == 3,0,1)

#create dataset
pretest_dataset <- make.dataset(respondent.dataset=pretest_ce,design.matrix=design_matrix,choice.indicators=c("q1","q2","q3","q4","q5","q6","q7","q8"))

##run simple conditional logit
fm <- RES ~ ASC + drought + flood + climate + water + price + strata(STR)
cl <- clogit(fm,pretest_dataset)
cl
gofm(cl)
#add confidence intervals
ci <- confint(cl,level=0.95)
