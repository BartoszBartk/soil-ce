############################
# Economic valuation of BODIUM-modelled soil management scenarios using data from a DCE
# KTBL Standarddeckungsbeitragsrechner webscraping code + inflation correction + mean prices
# December 2025
# Contributors: Bartosz Bartkowski, Samuel Fischer
# Contact: bartosz.bartkowski@ufz.de
# This script can be found at https://github.com/BartoszBartk/soil-ce/tree/main/BonaRes_test_case
############################

require(stringr)
require(rvest)
require(here)
require(dplyr)

ktbl_prices <- as.data.frame(matrix(ncol=4, nrow=0))
colnames(ktbl_prices) <- c("region","year","crop","price")

page <- read_html(here("KTBL-SDB-Druckansicht.html"))

# extract information from saved print page table
# 6 crop types, 2 regions, 20 years --> 240 entries
for (i in 3:243) {
  ktbl_prices[i,3] <- html_text2(html_element(page,paste0("table.printTable > tbody:nth-child(3) > tr:nth-child(1) > td:nth-child(1) > table:nth-child(1) > tbody:nth-child(1) > tr:nth-child(2) > td:nth-child(1) > table:nth-child(1) > tbody:nth-child(",i,") > tr:nth-child(2) > td:nth-child(1) > b:nth-child(1)")))
  ktbl_prices[i,2] <- html_text2(html_element(page,paste0("table.printTable > tbody:nth-child(3) > tr:nth-child(1) > td:nth-child(1) > table:nth-child(1) > tbody:nth-child(1) > tr:nth-child(2) > td:nth-child(1) > table:nth-child(1) > tbody:nth-child(",i,") > tr:nth-child(3) > td:nth-child(1)")))
  ktbl_prices[i,1] <- html_text2(html_element(page,paste0("table.printTable > tbody:nth-child(3) > tr:nth-child(1) > td:nth-child(1) > table:nth-child(1) > tbody:nth-child(1) > tr:nth-child(2) > td:nth-child(1) > table:nth-child(1) > tbody:nth-child(",i,") > tr:nth-child(4) > td:nth-child(1)")))
  ktbl_prices[i,4] <- html_text2(html_element(page,paste0("table.printTable > tbody:nth-child(3) > tr:nth-child(1) > td:nth-child(1) > table:nth-child(1) > tbody:nth-child(1) > tr:nth-child(2) > td:nth-child(1) > table:nth-child(1) > tbody:nth-child(",i,") > tr:nth-child(3) > td:nth-child(4)")))
}
ktbl_prices <- ktbl_prices[3:242,] # remove 2 empty rows at the top and 1 at the bottom
rownames(ktbl_prices) <- c(1:240) # correct row numbering

# force German coded numbers to numeric
ktbl_prices$price <- as.numeric(ktbl_prices$price)

# extract harvest year from year variable
ktbl_prices[,2] <- gsub(pattern = " ", replacement = "", ktbl_prices[,2])
ktbl_prices[,2] <- paste0(20,str_sub(ktbl_prices[,2], start = -2))
ktbl_prices[,2] <- as.numeric(ktbl_prices[,2])

# remove post-2021 years (since WTP data from 2021)
ktbl_prices <- subset(ktbl_prices, year < 2022)

# inflation-correct prices
infl <- 0.015 # rough estimate of mean inflation rate in Germany from 2004 till 2020, based on https://www.statista.com/statistics/262859/inflation-rate-in-germany-changes-of-the-cpi-compared-to-the-previous-year/
ktbl_prices$price <- ktbl_prices$price * (1 + infl)^(2021 - ktbl_prices$year)

# save output just in case
write.csv(ktbl_prices, here("ktbl_prices_halle_brandenburg_2004_2021.csv"), row.names = F, fileEncoding = "UTF-8")

# load again if needed
# ktbl_prices <- read.csv(here("ktbl_prices_halle_brandenburg_2004_2021.csv"))

# create a new, simple dataframe with median values for prices per crop for each of the two regions
median_prices <- as.data.frame(matrix(ncol = 4, nrow = 6))
colnames(median_prices) <- c("crop", "Brandenburg", "Halle", "mean_medians")

# either take ktbl_prices as is or, optionally, reduce number of included years
ktbl_p <- ktbl_prices
# ktbl_p <- subset(ktbl_prices, year > 2012)

# add crop names
median_prices$crop <- unique(ktbl_p$crop)

# add medians
x <- subset(ktbl_p, region == "Brandenburg")
y <- subset(ktbl_p, region == "Halle")
for (i in 1:nrow(median_prices)) {
  a <- subset(x, crop == median_prices$crop[i])
  median_prices$Brandenburg[i] <- median(a$price)
  b <- subset(y, crop == median_prices$crop[i])
  median_prices$Halle[i] <- median(b$price)
}

# calculate mean across regions
median_prices$mean_medians <- rowMeans(median_prices[, c(2:3)])

# save output
write.csv(median_prices, here("median_prices_ktbl.csv"), row.names= F, fileEncoding = "UTF-8")
