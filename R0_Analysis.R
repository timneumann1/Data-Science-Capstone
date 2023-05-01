####################################
#
#  Estimating the basic reproduction number (R0) of COVID-19 using a Monte Carlo method on the SIR model of differential equations  
#  
#  Author: Tim Neumann
#  
#  Reference and Copyright Information can be found at the end of the script
#
##################################################################################

##### R_0 Estimation with Poisson Error Function ######

####################
### Dependencies ###
####################

install.packages("tidyverse")
install.packages("readr")
install.packages("dplyr")
install.packages("estimatr")
install.packages("modelsummary")
install.packages("ggplot2")
install.packages("gridExtra") 
install.packages("chron")
install.packages("sfsmisc")
install.packages("deSolve")
install.packages("plotly")
install.packages("pacman")
install.packages("R0")
pacman::p_load(deSolve,chron,sfsmisc,gridExtra,ggplot2,tidyverse,readr,estimatr,modelsummary,dplyr,plotly,R0)

setwd("/Users/tim/Projects/DataScience_Capstone/")

data = data.frame(countries = c("Germany","India","South Africa", "Peru", "Russia"), median_percentile = c(2.044,1.105,1.123,1.135,1.240),subspace = c(1.809,1.094,1.131,1.136,1.244), package = c(2.490,1.149,1.267,1.081,1.522))

R0s = data.frame(x = data$countries, y = c(data$median_percentile, data$subspace,data$package), `Estimation Method` = c(rep("Median", nrow(data)),
                                                                                                          rep("Subspaces", nrow(data)),
                                                                                                          rep("R0 Package", nrow(data))))
 
plot = R0s %>% ggplot(aes(x, y, col= Estimation.Method)) + geom_point() + theme_bw() + theme(axis.text.x = element_text( size=10, angle=15),
                                                                                 axis.text.y = element_text( size=10, angle=0)) +    
                ggtitle(paste0("R0 estimates per country")) + xlab("Country") + 
                                                              ylab("R0 Estimate") 
plot

### Plotting the infection data including beginning of lockdowns

countries = c("Germany","India", "South Africa", "Peru", "Russia")
peaks = c(65,230,135,85,100)
names(peaks) = countries
lockdowns = c(as.Date("2020/03/22"), as.Date("2020/03/25"), as.Date("2020/03/27"), as.Date("2020/03/15"), as.Date("2020/03/28") ) 
names(lockdowns) = countries

# in some cases, there might not be a specific day at which a lockdown was implemented, but there is no ambiguity about the timing relative to the occurrence of the first infection peak
for (country in countries) {
  plot_data = data.frame(date = time_series_global[[paste0("time_series_",country)]]$date, cases = time_series_global[[paste0("time_series_",country)]]$cases )
  #print(min(plot_data$date) + peaks[country])
  plot = plot_data %>% ggplot(aes(x = date, y = cases)) + geom_point()  + 
    geom_vline(xintercept = min(plot_data$date) + peaks[country], color = "red", linewidth=0.5) +
    geom_vline(xintercept = lockdowns[country], linewidth = 0.5) + 
    ggtitle(paste0("Infection curve with peak and lockdown location"))
  print(plot)
  #ggsave(sprintf("./runs/%s/R0_Subspaces.png",country),ratioPlot)
}













