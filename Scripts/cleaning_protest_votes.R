############################
# Removing protest votes and cleaning data from the BonaRes soil ecosystem services CE
# n=1500, conducted online between 21 June and 6 July 2021
# authors: Bartosz Bartkowski, Julian Massenberg, Nele Lienhoop
# contact: bartosz.bartkowski@ufz.de
# This script can be found at https://github.com/BartoszBartk/soil-ce 
############################

require(support.CEs)
require(stringr)

#read in survey results
#note: this script uses the original dataset with one row per respondent, rather than the final, publicly available dataset from BonaRes repository; the procedure to remove protest votes should work nonetheless
#setwd()
ce_main_results <- read.csv("Data/BonaRes_CE_main_results.csv",header=T)

##remove protest votes
#add variable identifying respondents who always choose the status quo option
ce_main_results$sq_bias <- ifelse(ce_main_results$q1 + 
                                    ce_main_results$q2 + 
                                    ce_main_results$q3 + 
                                    ce_main_results$q4 + 
                                    ce_main_results$q5 + 
                                    ce_main_results$q6 +
                                    ce_main_results$q7 +
                                    ce_main_results$q8 +
                                    ce_main_results$q9==27,
                                    1,
                                    0)
#identify lowest median among response times for CE questions (excluding 1st question)
threshold <- min(median(ce_main_results$sys_pagetime_Random2),
                 median(ce_main_results$sys_pagetime_Random3),
                 median(ce_main_results$sys_pagetime_Random4),
                 median(ce_main_results$sys_pagetime_Random5),
                 median(ce_main_results$sys_pagetime_Random6),
                 median(ce_main_results$sys_pagetime_Random7),
                 median(ce_main_results$sys_pagetime_Random8),
                 median(ce_main_results$sys_pagetime_Fixed1))
#add variable identifying respondents with consistent below threshold response times
ce_main_results$resp_time <- ifelse(ce_main_results$sys_pagetime_Random2<threshold &
                                      ce_main_results$sys_pagetime_Random3<threshold &
                                      ce_main_results$sys_pagetime_Random4<threshold &
                                      ce_main_results$sys_pagetime_Random5<threshold &
                                      ce_main_results$sys_pagetime_Random6<threshold &
                                      ce_main_results$sys_pagetime_Random7<threshold &
                                      ce_main_results$sys_pagetime_Random8<threshold &
                                      ce_main_results$sys_pagetime_Fixed1<threshold,
                                      1,
                                      0)
#add variable identifying respondents with at least one "strongly agree" response to motivation questions 1 to 5
ce_main_results$protest <- ifelse(ce_main_results$motive_1==5 |
                                    ce_main_results$motive_2==5 |
                                    ce_main_results$motive_3==5 |
                                    ce_main_results$motive_4==5 |
                                    ce_main_results$motive_5==5,
                                    1,
                                    0)
#apply all three criteria to dataset
ce_main_results <- subset(ce_main_results,
                          sq_bias==0 |
                          resp_time==0 |
                          protest==0)
#remove unnecessary criteria variables and the threshold variable
ce_main_results <- ce_main_results[1:95]
rm(threshold)

####### further cleaning procedures for the original raw data (not relevant when using the publicly available dataset from BonaRes repository)
#correct variable 'plz'
ce_main_results$plz <- str_pad(ce_main_results$plz,5,pad="0")
ce_main_results$plz <- as.factor(ce_main_results$plz)

#save cleaned dataset for analysis of auxiliary questions
write.csv(ce_main_results,"Data/BonaRes_CE_main_results_clean.csv")

#extract choice responses
results_ce <- ce_main_results[c(2,77:86)]

#read in design file and create a design matrix
design_ce_main <- read.csv("Data/BonaRes_CE_main_design.csv",skip=23,header=T,sep=";")
#reorder data.frame according to blocks
design_ce_main <- design_ce_main %>% arrange(Block)
#give each choice task a new, block-specific number corresponding with final implemented design
design_ce_main$Choice.situation <- rep(c(1,2,3,4,6,7,8,9),30)
#add fixed task (first task from pretest design) to each block (as choice task 5)
fixed_task_original <- read.csv("Data/design_pretest.csv",skip=1,header=T,sep=";")
fixed_task <- fixed_task_original[1,] #extract first choice task
fixed_task <- rbind(fixed_task, fixed_task[rep(1,29),]) #create data.frame with 30 identical tasks
fixed_task$Block <- seq(1:30) #add block numbers
fixed_task$Choice.situation <- rep(5,30) #add choice task number
design_ce_main <- rbind(design_ce_main,fixed_task)
design_ce_main <- design_ce_main %>% arrange(Block)
#add status quo alternative to all tasks
design_ce_main$sq.drought <- rep(70,270)
design_ce_main$sq.flood <- rep(70,270)
design_ce_main$sq.climate <- rep(50,270)
design_ce_main$sq.water <- rep(30,270)
design_ce_main$sq.price <- rep(0,270)
#save the full design, including fixed task and status quo alternative
write.csv(design_ce_main,"Data/BonaRes_CE_main_design_full.csv")
#extract alternatives 1 to 3, add variable ALT = 1, 2, 3, rename columns
design_alt1 <- design_ce_main[,c(1:6,12)]
design_alt1$ALT <- rep(1,270)
colnames(design_alt1) <- c("QES","drought","flood","climate","water","price","BLOCK","ALT")
design_alt2 <- design_ce_main[,c(1,7:12)]
design_alt2$ALT <- rep(2,270)
colnames(design_alt2) <- c("QES","drought","flood","climate","water","price","BLOCK","ALT")
design_alt3 <- design_ce_main[,c(1,12:17)]
design_alt3$ALT <- rep(3,270)
colnames(design_alt3) <- c("QES","BLOCK","drought","flood","climate","water","price","ALT")
#combine by rows
design_matrix <- rbind(design_alt1,design_alt2,design_alt3)
#add ASC variable (0 for status quo alternative)
design_matrix$ASC <- ifelse(design_matrix$ALT == 3,0,1)
rownames(design_matrix) <- seq(1:810)
#remove unnecessary data frames used in the process
rm(design_alt1,design_alt2,design_alt3,fixed_task,fixed_task_original)

#create full dataset
ce_main_dataset <- make.dataset(respondent.dataset=results_ce,
                                design.matrix=design_matrix,
                                choice.indicators=c("q1","q2","q3","q4","q5","q6","q7","q8","q9"))
#add all other variables
rest_data <- ce_main_results[c(1,3:73,75,76,87:95)]
rest_data <- rest_data %>% slice(rep(1:n(),each=27))
rest_data <- rest_data[2:83]
ce_main_dataset <- cbind(ce_main_dataset,rest_data)
write.csv(ce_main_dataset,"Data/BonaRes_CE_results_clean.csv")
