############################
# Economic valuation of BODIUM-modelled soil management scenarios using data from a DCE
# DCE citation: @article{bartkowski_investigating_2022,
#     title = {Investigating preferences for soil-based ecosystem services},
#     volume = {2},
#     doi = {10.1093/qopen/qoac035},
#     number = {2},
#     journal = {Q Open},
#     author = {Bartkowski, Bartosz and Massenberg, Julian R and Lienhoop, Nele},
#     month = jan,
#     year = {2022},
#     pages = {qoac035}
#     }
# Autumn 2024
# Contributors: Bartosz Bartkowski, Samuel Fischer
# Contact: bartosz.bartkowski@ufz.de
# This script can be found at https://github.com/BartoszBartk/soil-ce/tree/main/BonaRes_test_case
############################

require(here)
require(dplyr)
require(reshape2)
require(ggplot2)

## Read first chunk of BODIUM data
bodium_data <- read.csv(here("Bodium/notill_diverse_mineral_26_sum.csv"), header = T, sep = "")

# add scenario column 
bodium_data$Scenario <- 'notill_diverse_mineral_26'

# add datasets with other management scenarios
d2 <- read.csv(here("Bodium/notill_diverse_mineral_85_sum.csv"), header = T, sep = "")
d2$Scenario <- 'notill_diverse_mineral_85'
d3 <- read.csv(here("Bodium/notill_simple_mineral_26_sum.csv"), header = T, sep = "")
d3$Scenario <- 'notill_simple_mineral_26'
d4 <- read.csv(here("Bodium/notill_simple_mineral_85_sum.csv"), header = T, sep = "")
d4$Scenario <- 'notill_simple_mineral_85'
d5 <- read.csv(here("Bodium/till_diverse_mineral_26_sum.csv"), header = T, sep = "")
d5$Scenario <- 'till_diverse_mineral_26'
d6 <- read.csv(here("Bodium/till_diverse_mineral_85_sum.csv"), header = T, sep = "")
d6$Scenario <- 'till_diverse_mineral_85'
d7 <- read.csv(here("Bodium/till_diverse_mixed_26_sum.csv"), header = T, sep = "")
d7$Scenario <- 'till_diverse_mixed_26'
d8 <- read.csv(here("Bodium/till_diverse_mixed_85_sum.csv"), header = T, sep = "")
d8$Scenario <- 'till_diverse_mixed_85'
d9 <- read.csv(here("Bodium/till_simple_mineral_26_sum.csv"), header = T, sep = "")
d9$Scenario <- 'till_simple_mineral_26'
d10 <- read.csv(here("Bodium/till_simple_mineral_85_sum.csv"), header = T, sep = "")
d10$Scenario <- 'till_simple_mineral_85'
d11 <- read.csv(here("Bodium/till_simple_mixed_26_sum.csv"), header = T, sep = "")
d11$Scenario <- 'till_simple_mixed_26'
d12 <- read.csv(here("Bodium/till_simple_mixed_85_sum.csv"), header = T, sep = "")
d12$Scenario <- 'till_simple_mixed_85'
bodium_data <- rbind(bodium_data,d2,d3,d4,d5,d6,d7,d8,d9,d10,d11,d12)
rm(d2,d3,d4,d5,d6,d7,d8,d9,d10,d11,d12)

# remove irrelevant columns
bodium_data <- subset(bodium_data, select = c("Scenario", "Year", "Plant", "Plant_biomass_SO", "Nitrate_leach", "Corg", "AWC"))

# restrict dataset to projection part
bodium_dataset <- subset(bodium_data, Year > 2020)

# recalculate plant biomass values to have them per ha and in tons
bodium_dataset$Plant_biomass_SO <- bodium_dataset$Plant_biomass_SO * 10000 / 1000
############# TEMPORARY SOLUTION BASED ON ROUGH ESTIMATE OF 1000 kg dry soil / m3
# recalculate Corg from kg/kg to kg C / ha
bodium_dataset$Corg <- bodium_dataset$Corg * 20000000

# add columns soil management (conventional vs. conservation), crop rotation (simple vs. diverse) and fertilization (mineral vs. mixed)
bodium_dataset <- bodium_dataset %>%
  mutate(Soil_mngmt = ifelse(substr(Scenario, 1, 4) == "till", "conventional", "conservation"), #if scenario starts with 'till', 'conventional' is printed
         Rotation = ifelse(grepl("diverse", Scenario), "diverse", "simple"),
         Fertlzr = ifelse(grepl("mineral", Scenario), "mineral", "mixed")
  )

## Read WTP data (recalculated for relevant units)
wtp_data <- read.csv("WTPs_recalculated.csv", header = T, sep = ";")

## Read price/cost data for food provision service
price_data <- read.csv("costs_prices_food.csv", header = T, sep = ";", dec = ",")

## Add gross margins based on price_data 
full_dataset <- merge(bodium_dataset, price_data, by.x = c("Plant","Soil_mngmt"), by.y = c("Crop","Soil_mngmt"))
# add costs of cover crops to winter wheat years in diverse scenarios
full_dataset <- full_dataset %>% mutate(Costs_ha = ifelse(Rotation == "diverse" & Plant == "winterwheat", Costs_ha + 202.22, Costs_ha))
full_dataset$GM <- full_dataset$Plant_biomass_SO * full_dataset$Price_ton - full_dataset$Costs_ha

## Add WTP estimates
full_dataset$WTP_Corg <- full_dataset$Corg * as.numeric(wtp_data[wtp_data$ES == "climate",]$WTP_per_unit)
full_dataset$WTP_AWC <- full_dataset$AWC * as.numeric(wtp_data[wtp_data$ES == "drought",]$WTP_per_unit)
full_dataset$WTP_nitrate <- full_dataset$Nitrate_leach * as.numeric(wtp_data[wtp_data$ES == "water",]$WTP_per_unit)

## Calculate total value
full_dataset$tot_wtp <- full_dataset$WTP_AWC + full_dataset$WTP_Corg + full_dataset$WTP_nitrate

## Add total values per year to new scenario dataframe
# simplify full_dataset
temp_data <- subset(full_dataset, select = c("Scenario", "Year", "GM", "tot_wtp"))

##################################
# optional: discount the values
r <- 0.04 # set discount rate for GM (private good)
p <- 0.01 # set discount rate to tot_wtp (public goods)
temp_data$GM <- temp_data$GM/(1 + r)^(temp_data$Year - min(temp_data$Year))
temp_data$tot_wtp <- temp_data$tot_wtp/(1 + p)^(temp_data$Year - min(temp_data$Year))
##################################

# reshape to long format to have all monetary values in one column
temp_data_ <- melt(temp_data,
                  id.vars = c("Scenario", "Year"))
# reorder by year
temp_data_ <- temp_data_[order(temp_data_$Year),]
# merge Scenario and variable (GM vs. tot_wtp)
temp_data_$Scenario <- paste0(temp_data_$Scenario, "_", temp_data_$variable)
# remove redundant "variable" column
temp_data_ <- temp_data_ %>% select(-variable)

# reshape to create target dataset
val_data_ <- reshape(temp_data_, direction = "wide",
                    idvar = "Scenario",
                    timevar = "Year",
                    v.names = "value")
colnames(val_data_) <- c("Scenario", c(2021:2049))
val_data_ <- val_data_[order(val_data_$Scenario),]

# add (discounted) sum for each scenario (for comparability, reduce to 24 data points, as crop rotation length is 6 in "diverse")
val_data_$NPV <- rowSums(val_data_[, c(2:25)])
val_data_$mean <- rowMeans(val_data_[, c(2:25)])
val_data_$sd <- apply(val_data_[, c(2:25)], 1, sd)

# ALTERNATIVELY: do the same with temp_data to work with sums of GM & tot_wtp
temp_data$tot_val <- temp_data$GM + temp_data$tot_wtp
temp_data <- temp_data %>% select(-c(GM, tot_wtp))
# reorder by year
temp_data <- temp_data[order(temp_data$Year),]
# reshape to create target dataset
val_data <- reshape(temp_data, direction = "wide",
                    idvar = "Scenario",
                    timevar = "Year",
                    v.names = "tot_val")
colnames(val_data) <- c("Scenario", c(2021:2049))
val_data <- val_data[order(val_data$Scenario),]
# add (discounted) sum for each scenario (24 years only, see explanation above)
val_data$NPV <- rowSums(val_data[, c(2:25)])
val_data$mean <- rowMeans(val_data[, c(2:25)])
val_data$sd <- apply(val_data[, c(2:25)], 1, sd)

## Create plots
# create plots for each scenario
# sums
boxplots <- ggplot(temp_data, aes(x = Scenario, y = tot_val)) + 
  geom_boxplot()
boxplots
# separate GM / tot_wtp
boxplots_ <- ggplot(temp_data_, aes(x = Scenario, y = value)) + 
  geom_boxplot()
boxplots_
# sums
timeline <- ggplot(temp_data, aes(x = Year, y = tot_val, color = Scenario)) +
  geom_line()
timeline
# separate
timeline_ <- ggplot(temp_data_, aes(x = Year, y = value, color = Scenario)) +
  geom_line()
timeline_

# create NPV bar chart across scenarios [CURRENTLY BASED ON val_data; WHICH DOESN'T SUM UP GM AND val_wtp]
scenario_bars <- ggplot(val_data, aes(x = Scenario, y = NPV)) +
  geom_bar(stat = "identity")
scenario_bars
