############################
# Economic valuation of BODIUM-modelled soil management scenarios using data from a DCE
# Main DCE paper: Bartkowski et al., 2022, DOI: 10.1093/qopen/qoac035
# December 2025
# Contributors: Bartosz Bartkowski, Samuel Fischer
# Contact: bartosz.bartkowski@ufz.de
# This script can be found at https://github.com/BartoszBartk/soil-ce/tree/main/BonaRes_test_case
############################

require(here)
require(dplyr)
require(reshape2)
require(ggplot2)
require(tidyverse)
require(ggnewscale)
require(ggh4x)
require(patchwork)

##########
# 1. BODIUM data
##########

## Read summary Data
bodium_data <- read.csv(here("sum_data.csv"), header = T, sep = ",")

# remove irrelevant columns
bodium_data <- subset(bodium_data, select = c("Year", "Plant", "Plant_biomass_SO", "Nitrate_leach", "Corg_200", 
                                              "AWC", "Weight_200","site","prec","tillage","cr","rcp"))

# merge values from cover crop years to keep cash crop yields and all other values from another crop
# first, add unique ID
bodium_data$ID <- paste0(bodium_data$Year, bodium_data$site, bodium_data$rcp, bodium_data$cr, bodium_data$tillage, bodium_data$prec)
# code to take non-zero values for selected variables
nz_remove <- function(x, default = 0) {
  x <- x[!is.na(x) & x != 0]
  if (length(x)) x[1] else default
}
# define variables to coalesce by non-zero
vars_rel <- c("Nitrate_leach", "Corg_200", "AWC")
# apply to dataset
bodium_dataset <- bodium_data %>%
  group_by(ID) %>%
  summarise(
    # keep wheat and its biomass
    Plant = {
      yp <- nz_remove(Plant[Plant != "phacelia"], default = NA_real_)
      ifelse(is.na(yp) || yp == 0, nz_remove(Plant), yp)
    },
    Plant_biomass_SO = {
      yw <- nz_remove(Plant_biomass_SO[Plant != "phacelia"], default = NA_real_)
      ifelse(is.na(yw) || yw == 0, nz_remove(Plant_biomass_SO), yw)
    },
    # keep variables that don't vary between cash crop and cover crop
    across(where(is.character) & !any_of("Plant"), ~ first(.x[.x != ""], default = "")),
    # coalesce relevant variables (and all numeric variables that don't vary)
    across(where(is.numeric) & !any_of("Plant_biomass_SO"), ~ nz_remove(.x)),
    .groups = "drop"
    )
# rearrange dataframe
bodium_dataset <- arrange(bodium_dataset, site, rcp, cr, tillage, prec, Year)

# recalculate plant biomass values to have them per ha and in tons
bodium_dataset$Plant_biomass_SO <- bodium_dataset$Plant_biomass_SO * 10000 / 1000

# recalculate Corg from kg/kg to kg C / ha based on weight variable
bodium_dataset$Corg_200 <- bodium_dataset$Corg_200 * 10000 * bodium_dataset$Weight_200

##########
# 2. Economic value data
##########

## Read WTP data (recalculated for relevant units)
wtp_data <- read.csv(here("WTPs_recalculated.csv"), header = T, sep = ",")

## Read price/cost data for food provision service
price_data <- read.csv(here("costs_prices_food_meanmedian.csv"), header = T)  

##########
# 3. Merge datasets
##########

## Add gross margins based on price_data 
full_dataset <- merge(bodium_dataset, price_data, by.x = c("Plant","tillage"), by.y = c("Crop","tillage"))
# add costs of cover crops to winter wheat years in diverse scenarios
full_dataset <- full_dataset %>% mutate(Costs_ha = ifelse(cr == "diverse" & Plant == "winterwheat", Costs_ha + 177, Costs_ha))
full_dataset$GM <- full_dataset$Plant_biomass_SO * full_dataset$Price_ton - full_dataset$Costs_ha

## Add WTP estimates
full_dataset$WTP_Corg <- full_dataset$Corg_200 * as.numeric(wtp_data[wtp_data$ES == "climate",]$WTP_per_unit)
full_dataset$WTP_AWC <- full_dataset$AWC * as.numeric(wtp_data[wtp_data$ES == "drought",]$WTP_per_unit)
full_dataset$WTP_nitrate <- full_dataset$Nitrate_leach * as.numeric(wtp_data[wtp_data$ES == "water",]$WTP_per_unit) * -1

## Calculate total value
full_dataset$tot_wtp <- full_dataset$WTP_AWC + full_dataset$WTP_Corg + full_dataset$WTP_nitrate

##########
# 4. Yearly discounted economic values per soil function/ES
##########

# remove irrelevant variables from dataset
temp_data <- subset(full_dataset, select = c("tillage","cr","rcp","site","prec","Year", "GM", "WTP_Corg", "WTP_AWC", "WTP_nitrate",
                                              "tot_wtp"))

# create a Scenario variable
temp_data$rcp <- as.character(temp_data$rcp)
temp_data$Scenario <- interaction(temp_data$tillage,temp_data$cr,temp_data$rcp,temp_data$site,temp_data$prec)
temp_data <- subset(temp_data, select = c("Scenario","Year", "GM", "WTP_Corg", "WTP_AWC", "WTP_nitrate", "tot_wtp"))

# discount the values
r <- 0.04 # set discount rate for GM (private good)
p <- 0.01 # set discount rate for WTP (public goods)
temp_data$GM <- temp_data$GM/(1 + r)^(temp_data$Year - min(temp_data$Year))
temp_data$tot_wtp <- temp_data$tot_wtp/(1 + p)^(temp_data$Year - min(temp_data$Year))
temp_data$WTP_Corg <- temp_data$WTP_Corg/(1 + p)^(temp_data$Year - min(temp_data$Year))
temp_data$WTP_AWC <- temp_data$WTP_AWC/(1 + p)^(temp_data$Year - min(temp_data$Year))
temp_data$WTP_nitrate <- temp_data$WTP_nitrate/(1 + p)^(temp_data$Year - min(temp_data$Year))

# calculate total economic value per scenario
temp_data$tot_val <- temp_data$GM + temp_data$tot_wtp

# order by year
temp_data <- temp_data[order(temp_data$Year),]

# create a dataset with an economic value for each scenario/year/soil function combination
temp_data_ <- melt(temp_data,
                    id.vars = c("Scenario", "Year"))
# order by year again
temp_data_ <- temp_data_[order(temp_data_$Year),]

# reshape back to wide for target dataset with yearly values for each scenario/soil function combo
val_data <- reshape(temp_data_, direction = "wide",
                     idvar = c("Scenario", "variable"),
                     timevar = "Year",
                     v.names = "value")

# rename columns for later visualization purposes
colnames(val_data) <- c("Scenario", "Group", c(2021:2050))

# order by scenario
val_data <- val_data[order(val_data$Scenario),]

# add net present value (sum of discounted values)
val_data$NPV <- rowSums(val_data[, c(3:32)])
val_data2 <- subset(val_data, Group != "tot_wtp" & Group != "tot_val") #leave out tot_wtp and tot_val (only for visualization purposes)

## create dataset for final overview of economic analysis
# recreate the individual scenario components
val_data3 <- val_data2 %>%
  separate(Scenario, into = c("Management", "Rotation", "RCP", "Site", "Precipitation"), sep = "\\.", remove = F)
# select only relevant variables, first restricting them to tot_val
val_data4 <- val_data %>%
  separate(Scenario, into = c("Management", "Rotation", "RCP", "Site", "Precipitation"), sep = "\\.", remove = F)
val_data4 <- subset(val_data3, Group == "tot_val")
val_data4 <- subset(val_data4, select = c("Scenario", "Site", "Management", "Rotation", "RCP", "Precipitation", "NPV"))
# create and order individual subsets for each site
bad_lauchst <- subset(val_data4, Site == "bl")
bad_lauchst <- bad_lauchst[order(-bad_lauchst$NPV),]
muencheb <- subset(val_data4, Site == "mb")
muencheb <- muencheb[order(-muencheb$NPV),]

##########
# 5. Plotting
##########

# stacked bar chart using NPV, across scenarios
bar <- ggplot(val_data2, aes(x = Scenario, y = NPV, fill = Group)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("GM" = "darkred","WTP_Corg" = "darkgreen","WTP_AWC" = "orange","WTP_nitrate" = "darkblue")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 
bar


# Stacked bar chart with table below graph
sites <- unique(val_data3$Site)

# BL
df1 <- val_data3 %>% filter(Site == sites[1])

# 1) bar chart as before
df1_tot <- df1 %>% group_by(Scenario) %>% summarise(NPV_total = sum(NPV), .groups="drop")
p_bar <- ggplot(df1, aes(Scenario, NPV, fill = Group)) + geom_bar(stat="identity", position="stack") +
  geom_point(data = df1_tot,aes(x = Scenario, y = NPV_total, shape = "Total NPV"),inherit.aes = FALSE,color = "black",size = 2.5) +
  scale_fill_manual(values=c("GM"="darkred","WTP_Corg"="darkgreen","WTP_AWC"="orange","WTP_nitrate"="darkblue")) +
  scale_shape_manual(values = c("Total NPV" = 16), name = NULL) + theme_classic() +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1)) + labs(title = paste("Site:", sites[1]))

# 2) overview: first distinct, then pivot_longer
df_axis_long <- df1 %>%
  distinct(Scenario, Management, Rotation, RCP, Precipitation) %>%
  pivot_longer(cols = c(Precipitation, RCP, Rotation, Management),names_to = "row",values_to = "label") %>%
  mutate(row = factor(row, levels = c("Precipitation", "RCP", "Rotation", "Management")),
         row_label = interaction(row, label, drop = TRUE))

pal <- c("Precipitation.dry" = "#c7e9b4","Precipitation.normal" = "#7fcdbb","Precipitation.wet" = "#2c7fb8",
         "RCP.26" = "#fee08b","RCP.85" = "#fdae61","Rotation.diverse" = "#EFE6FA","Rotation.simple" = "#D8C7F3",
         "Management.notill" = "#FAD4D4","Management.till"    = "#F4A6A6")

p_overview <- ggplot(df_axis_long, aes(x = Scenario, y = row)) + geom_tile(aes(fill = row_label), height = 0.9) +
  geom_text(aes(label = label), size = 3) + scale_fill_manual(values = pal, guide = "none") +  theme_minimal(base_size = 11) +
  theme(axis.title = element_blank(),axis.text.x = element_blank(),panel.grid = element_blank(),
        plot.margin = margin(t = 0, r = 10, b = 0, l = 10)) + coord_cartesian(clip="off")

# 3) join
(p_bar / p_overview) + plot_layout(heights = c(4, 1.3))


# same for MB
df2 <- val_data3 %>% filter(Site == sites[2])

# 1) bar chart
df2_tot <- df2 %>% group_by(Scenario) %>% summarise(NPV_total = sum(NPV), .groups="drop")
p_bar_mb <- ggplot(df2, aes(Scenario, NPV, fill = Group)) + geom_bar(stat="identity", position="stack") +
  geom_point(data = df2_tot,aes(x = Scenario, y = NPV_total, shape = "Total NPV"),inherit.aes = FALSE,color = "black",size = 2.5) +
  scale_fill_manual(values=c("GM"="darkred","WTP_Corg"="darkgreen","WTP_AWC"="orange","WTP_nitrate"="darkblue")) +
  scale_shape_manual(values = c("Total NPV" = 16), name = NULL) + theme_classic() +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1)) + labs(title = paste("Site:", sites[2]))

# 2) overview
df_axis_long_mb <- df2 %>%
  distinct(Scenario, Management, Rotation, RCP, Precipitation) %>%
  pivot_longer(cols = c(Precipitation, RCP, Rotation, Management),names_to = "row",values_to = "label") %>%
  mutate(row = factor(row, levels = c("Precipitation", "RCP", "Rotation", "Management")),
         row_label = interaction(row, label, drop = TRUE))

p_overview_mb <- ggplot(df_axis_long_mb, aes(x = Scenario, y = row)) + geom_tile(aes(fill = row_label), height = 0.9) +
  geom_text(aes(label = label), size = 3) + scale_fill_manual(values = pal, guide = "none") +  theme_minimal(base_size = 11) +
  theme(axis.title = element_blank(),axis.text.x = element_blank(),panel.grid = element_blank(),
        plot.margin = margin(t = 0, r = 10, b = 0, l = 10)) + coord_cartesian(clip="off")

# 3) join
(p_bar_mb / p_overview_mb) + plot_layout(heights = c(4, 1.3))



