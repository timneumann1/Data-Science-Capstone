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

setwd("/Users/tim/Projects/DS_Capstone/")

data = data.frame(countries = c("Germany","India","South Africa", "Peru", "Russia"), median_percentile = c(2.044,1.105,1.123,1.135,1.240),subspace = c(1.809,1.094,1.131,1.136,1.244), package = c(2.490,1.149,1.267,1.081,1.522))

R0s = data.frame(x = data$countries, y = c(data$median_percentile, data$subspace,data$package), group = c(rep("Median", nrow(data)),
                                                                                                          rep("Subspaces", nrow(data)),
                                                                                                          rep("Package", nrow(data))))
 
plot = R0s %>% ggplot(aes(x, y, col= group)) + geom_point() + theme_bw() + theme(axis.text.x = element_text( size=10, angle=15),
                                                                                 axis.text.y = element_text( size=10, angle=0))
plot











