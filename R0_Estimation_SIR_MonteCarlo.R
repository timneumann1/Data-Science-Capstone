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

###############################
######## FUNCTIONS ############
###############################

#################################
#### Derivatives Calculation ####
#################################

SIR_derivatives=function(t, x, vparameters){
  S = x[1]  
  I = x[2]  
  R = x[3]  
  I[I<0] = 0
  S[S<0] = 0
  R[R<0] = 0
  
  with(as.list(vparameters),{
    npop = S+I+R   
    # SIR model
    dS = -beta*S*I/npop            
    dI = +beta*S*I/npop - gamma*I  
    dR = +gamma*I                  

    vout = c(dS,dI,dR)
    list(vout)
  })
}

#################################
#### Monte Carlo ####
#################################
monte_carlo = function(data,iterations,range_beta,range_gamma,name,population,initialI,initialR,directory_name){
  set.seed(123)
  incidence_observed = data$cases
  removed_observed = data$removed
  times_of_observed = data$days_rel_start
  time_binning = min(diff(times_of_observed))
  
  vBeta = numeric(0)
  vGamma = numeric(0)
  vpois_negloglike = numeric(0) 

  best_pois_negloglike_so_far = 1e10 # normally replaced by first error term
  vbest_pois_negloglike_fit_incidence_prediction = rep(0,length(incidence_observed))
  vbest_pois_negloglike_fit_removed_prediction = rep(0,length(removed_observed))

  b_range = range_beta 
  g_range = range_gamma 
  
  beta_A = 1e10
  gamma_A = 1e10
  
  ###################################
  ##### Monte Carlo iterations ######
  ###################################
  
  cat("Running... \n")  
  
  for (iter in 1:iterations){
    
    npop = population   # size of population
    I_0 = initialI      
    R_0 = initialR      # initial recovered and immune
    S_0 = npop-I_0-R_0
    
    # randomly sample beta and gamma uniformly
    beta = runif(1,min(range_beta),max(range_beta)) 
    gamma = runif(1,min(range_gamma),max(range_gamma))         
    
    t0 = -1   # t0 is the time-of-introduction of the disease to the population, by construction 
              # of our dataframes, this will be the day before the first infection is reported

    vparameters = c(beta=beta, gamma=gamma)
    inits = c(S=S_0,I=I_0,R=R_0)
    
    tmin = t0
    tmax = max(times_of_observed)
    tmin = min(t0,min(times_of_observed))
    
    vt = seq(tmin,tmax)
    
    solved_model = as.data.frame(lsoda(inits, vt, SIR_derivatives, vparameters)) # remember: the COVID data is incidence, not prevalence 
    tmin_data = min(times_of_observed)-time_binning
    tmax_data = max(times_of_observed)
    vtime_data = seq(tmin_data,tmax_data,time_binning)
    
    susceptible_predicted = solved_model$S[solved_model$time%in%vtime_data]  
    removed_predicted = solved_model$R[solved_model$time%in%vtime_data]  
    incidence_predicted = -diff(susceptible_predicted)
    removed_predicted = diff(removed_predicted)
    
    frac_confirmed = sum(incidence_observed)/sum(incidence_predicted)
    frac_confirmed_removed = sum(removed_observed)/sum(removed_predicted)
    incidence_predicted = incidence_predicted*frac_confirmed 
    removed_predicted = removed_predicted*frac_confirmed_removed 
    
    # Poisson neglog likelihood statistic that compares the data to this model calculated
    # under a particular hypothesis of beta and gamma
    
    if ((length(incidence_predicted)==length(incidence_observed))
        & (!is.na(sum(incidence_predicted)))  & ( 
          length(removed_predicted)==length(removed_observed))
        & (!is.na(sum(removed_predicted))) & (min(removed_predicted) > 0) ) #
      # last restriction is needed in order to work with the Poisson negative log likelihood,
      
    {
      
      ###########################################################################
      # calculate the Poisson negative log-likelihood statistic
      ###########################################################################
      
      pois_negloglike = sum(-incidence_observed*log(incidence_predicted)+incidence_predicted) +
        sum(-removed_observed*log(removed_predicted)+removed_predicted)
      vBeta = append(vBeta,beta)
      vGamma = append(vGamma,gamma)

      vpois_negloglike = append(vpois_negloglike,pois_negloglike)

      if (pois_negloglike<best_pois_negloglike_so_far){
        best_pois_negloglike_so_far = pois_negloglike
        vbest_pois_negloglike_fit_incidence_prediction = incidence_predicted
        vbest_pois_negloglike_fit_removed_prediction = removed_predicted
        beta_best = beta
        gamma_best = gamma
      }

    } # end check that the predicted incidence vector is the same length as the observed and doesn't contain NA's
  } # end loop over the Monte Carlo iterations
      
  
  #################################
  #### Plots (Poisson) ############
  #################################  
  
  num_points_to_show = 10000
  ymax = max(vpois_negloglike)  # filter inf values
  if (length(vpois_negloglike)>num_points_to_show) ymax = sort(vpois_negloglike)[num_points_to_show]
  l = which(vpois_negloglike<=ymax)
  lmin = which.min(vpois_negloglike)

  # Save in directory
  
  #### R0 #####
  
  png(filename=paste0(directory_name,"/R0.png"))

  plot(vBeta[l]/vGamma[l]
       ,vpois_negloglike[l]
       ,ylab="Poisson negative log-likelihood"
       ,xlab="R0 hypothesis"
       ,main=paste("Best-fit (beta,gamma,R0):",round(beta_best,3),round(gamma_best,3),round(beta_best/gamma_best,3))
       ,col.main=3
       ,xaxp = c(-1, 9, 1))
  points(vBeta[lmin]/vGamma[lmin],vpois_negloglike[lmin],col=3,cex=2)

  dev.off()

  #### Incidence #####

  png(filename=paste0(directory_name,"/incidence.png"))

  ymax = max(c(incidence_observed,vbest_pois_negloglike_fit_incidence_prediction))
  plot(times_of_observed
       ,incidence_observed
       ,ylim=c(0,1.2*ymax)
       ,xlab="Time in days relative to outbreak"
       ,ylab="Incidence"
       ,cex=2)
  lines(times_of_observed
        ,vbest_pois_negloglike_fit_incidence_prediction
        ,col=2
        ,lwd=5)

  dev.off()

  #### Removed #####

  png(filename=paste0(directory_name,"/removed.png"))

  ymax = max(c(removed_observed,vbest_pois_negloglike_fit_removed_prediction))
  plot(times_of_observed
       ,removed_observed
       ,ylim=c(0,1.2*ymax)
       ,xlab="Time in days relative to outbreak"
       ,ylab="Removed"
       ,cex=2)
  lines(times_of_observed
        ,vbest_pois_negloglike_fit_removed_prediction
        ,col=2
        ,lwd=5)

  dev.off()

  ### error terms 3d plot ###
  
  error_terms = vpois_negloglike
  betas = vBeta
  gammas = vGamma
  .GlobalEnv$plot_data = data.frame(error_terms[l],betas[l],gammas[l])
  threedee = plot_ly(plot_data, x = ~gammas.l., y = ~betas.l., z = ~error_terms.l., size = .75, type="scatter3d",mode="markers")
  htmlwidgets::saveWidget(as_widget(threedee), paste0(directory_name,"/Errors_3d.html"))

  ######## Iteration Result ########

  print(sprintf("The best pair (beta,gamma) is (%.3f,%.3f), resulting in a best-fit R0 of %.3f.",beta_best,gamma_best,beta_best/gamma_best) )
  cat("\n")
}

##########################
######### Data ###########
##########################

# Global Covid Data (Johns Hopkins University): https://github.com/CSSEGISandData

####################################
#### Compartment: Infected data #####
####################################

df = read_csv("data/time_series_covid19_confirmed_global.csv") %>% dplyr::select(-c("Province/State","Lat","Long"))

time_series_global = list()  # list of data frames per country
country_labels = list()
i = 1

for (country in c("Germany","India", "South Africa", "Peru", "Russia")) {  
  print(country)
  print(i)
  # determine start date
  data = df %>% filter(`Country/Region` == country) %>% dplyr::select(-c("Country/Region")) %>%
    pivot_longer(
      everything(),
      names_to = "date",
      values_to = "cases"
    )
  
  data$date = as.Date(data$date, "%m/%d/%y")
  temp = data$cases[1]
  for (x in 2:nrow(data)) {
    temp2 = data$cases[x] 
    data$cases[x] = temp2-temp
    temp = temp2
  }
  index = which(data$cases!=0)[1]
  first_day = data$date[index]
  data$days_rel_start = julian(data$date,as.Date(first_day,"%m-%d-%y")) # start date of pandemic in country
  data_AoI = data %>% filter(days_rel_start>=0 & days_rel_start < 250) # 90-130
  
  country_labels[[i]] = paste0("time_series_",country)
  time_series_global[[i]] = data_AoI
  i = i+1
}

names(time_series_global) = country_labels

####################################
#### Compartment: Removed data #####
####################################

df = read_csv("data/time_series_covid19_recovered_global.csv") %>% dplyr::select(-c("Province/State","Lat","Long"))
df2 = read_csv("data/time_series_covid19_deaths_global.csv") %>% dplyr::select(-c("Province/State","Lat","Long"))

df = df %>%  mutate_all(~ ifelse(is.na(.), 0, .))
df2 = df2 %>%  mutate_all(~ ifelse(is.na(.), 0, .))  

time_series_removed_global = list()  # list of data frames per country

country_labels = list()

i = 1
for (country in c("Germany","India", "South Africa", "Peru", "Russia")) {
  print(country)
  print(i)
  
  #### RECOVERIES ###### 
  
  data = df %>% filter(`Country/Region` == country) %>% dplyr::select(-c("Country/Region")) %>%
    pivot_longer(
      everything(),
      names_to = "date",
      values_to = "recoveries"
    )
  
  data$date = as.Date(data$date, "%m/%d/%y")
  temp = data$recoveries[1]
  for (x in 2:nrow(data)) {
    temp2 = data$recoveries[x] 
    data$recoveries[x] = temp2-temp
    temp = temp2
  }
  
  #### DEATHS ######  
  
  data2 = df2 %>% filter(`Country/Region` == country) %>% dplyr::select(-c("Country/Region")) %>%
    pivot_longer(
      everything(),
      names_to = "date",
      values_to = "deaths"
    )
  
  data2$date = as.Date(data2$date, "%m/%d/%y")
  temp = data2$deaths[1]
  for (x in 2:nrow(data2)) {
    temp2 = data2$deaths[x] 
    data2$deaths[x] = temp2-temp
    temp = temp2
  }
  
  df3 = merge(data,data2, by = "date")
  df3$removed = df3$recoveries + df3$deaths
  df3 = df3 %>% dplyr::select(date, removed)
  
  index = which(df3$removed!=0)[1]
  first_day = df3$date[index]
  df3$days_rel_start = julian(df3$date,as.Date(first_day,"%m-%d-%y")) # start date of pandemic in country
  data_AoI = df3 %>% filter(days_rel_start>=0 & days_rel_start < 250) # must contain the first peak
  data_AoI = data_AoI %>%  mutate(removed =ifelse(removed < 0, 0, removed)) 
  country_labels[[i]] = paste0("time_series_removed_",country)
  time_series_removed_global[[i]] = data_AoI
  
  i = i+1
}

names(time_series_removed_global) = country_labels

#####################################
############ ANALYSIS ###############
#####################################

population_numbers = c(81000000,1396000000,58800000,33300000,144100000)
names(population_numbers) = c("Germany","India", "South Africa", "Peru", "Russia") 

for (country in c("Germany","India","South Africa", "Peru", "Russia")) {  
  vR0 = list()
  print(country)
  dir_name = paste0("./runs/",country)
  dir.create(dir_name, showWarnings = FALSE)
  sdir_name = paste0(dir_name,"/plots/")
  dir.create(sdir_name,showWarnings = FALSE)
  
  population_in = population_numbers[country]  #  population of country in January 2020
  
  df_cases = time_series_global[paste0("time_series_",country)][[1]]
  df_removed = time_series_removed_global[paste0("time_series_removed_",country)][[1]] %>% dplyr::select(-"days_rel_start")
  
  df = merge(df_cases,df_removed, by = "date", all.x = TRUE) %>% dplyr::select(date,days_rel_start,cases,removed) 
  df = df %>%  mutate (cases = ifelse(is.na(cases)|cases<0, 0, cases), removed = ifelse(is.na(removed)|removed<0, 0, removed))
  
  print(df %>% ggplot(aes(x = days_rel_start, y = cases)) + geom_point() + theme_bw() + scale_x_continuous("Days from start of the pandemic", breaks=c(0:50)*5) + theme(axis.text.x = element_text(angle = 90)))
  
  # after the plot is displayed, the user inputs the peak location in days from start of the pandemic
  print(sprintf("Based on the plot, please estimate the time of the first peak, relative to the outbreak on %s", as.character(df$date[1])))
  peak = as.integer(readline())
  
  # Germany: 65, India: 230, South Africa: 135, Peru: 85, Russia: 100
  
  nrows_range = c((peak-5):(peak+5)) # sample 11 intervals centered at peak estimate
  print(paste0("Running Monte Carlo simulations for ",country))
  
  .GlobalEnv$errorBetaGamma_data = data.frame()
  
  for (timeframe in nrows_range) {
    print(paste0("Starting Monte carlo simulation for subset of ", timeframe, " rows"))
    subdirectory_name = paste0(dir_name,"/plots/",timeframe)
    dir.create(subdirectory_name, showWarnings = FALSE)
    
    .GlobalEnv$plot_data = data.frame()
    monte_carlo(data = df[1:timeframe,],iterations = 1500, range_beta = c(1/16,3), range_gamma = c(1/8,1/2) ,name = country, population = population_in[[1]], initialI = df$cases[1], initialR = df$removed[1], directory_name = subdirectory_name) 
    # we inform gamma by reality (infectious period), then let beta such that R_0 could be between 1/2 and 6
    
    R0_estimates = plot_data[order(plot_data$error_terms.l.),][1:15,]
    R0_estimates$R0 = R0_estimates$betas.l. / R0_estimates$gammas.l.
    R0sTemp = R0_estimates$R0
    vR0 = append(vR0,R0sTemp)
    
    .GlobalEnv$errorBetaGamma_data = rbind(.GlobalEnv$errorBetaGamma_data,R0_estimates)
  }
  
  ######################################################################################################
  ## Estimating R0 with percentiles (determining non-parametric median and 90% confidence interval)  ###
  ######################################################################################################
  
  R0s = sort(unlist(vR0))
  a=quantile(R0s,c(0,.025,.05,.25,.5,.75,.95,.975,1))
  hist(R0s,breaks = 10)
  abline(v=a[5], col='green', lwd=3, lty='dashed')
  abline(v=a[3], col='blue', lwd=3)
  abline(v=a[7], col='blue', lwd=3)
  arrows(a[3], 2, a[7], 2)
  arrows(a[7], 2, a[3], 2)
  
  text(x=(a[3]+a[7])/2, y=4,'90 % CI', cex=.75)
  
  print(sprintf("The final estimate on R0 for %s is %.3f with a 90 percent confidence interval (%.3f,%.3f). ", country, a[5], a[3], a[7]))
  
  fileConn = file(sprintf("./runs/%s/R0_Percentiles.txt",country))
  writeLines(sprintf("The final estimate on R0 for %s is %.3f with a 90 percent confidence interval (%.3f,%.3f). ", country, a[5], a[3], a[7]), fileConn)
  close(fileConn)
  
  ######################################################################################################
  ##################  Estimating R0 with best-fit subspace #############################################
  ######################################################################################################

  regressionR0 = lm_robust(errorBetaGamma_data$betas.l. ~ 0 + errorBetaGamma_data$gammas.l.)
  coefs = summary(regressionR0)$coefficients
  write.table(coefs, file=sprintf("./runs/%s/R0_Subspaces-Regression.csv",country), sep = ",", col.names=NA)
  ratioPlotData = data.frame(beta = errorBetaGamma_data$betas.l., gamma = errorBetaGamma_data$gammas.l.)
  ratioPlot = ratioPlotData %>% ggplot(aes(x = gamma, y = beta)) + geom_point() + xlim(0,0.5) +ylim(0,1.2) + 
                                    geom_abline(intercept = 0, slope = coefs[1]) + ggtitle(paste0("R0 = ", coefs[1]))
  ggsave(sprintf("./runs/%s/R0_Subspaces.png",country),ratioPlot)
  
  ######################################################################################################
  ##################  Estimating R0 with R0 Estimation Package #########################################
  ######################################################################################################
  
  data_covid = time_series_global[paste0("time_series_",country)][[1]]
  incid = data_covid$cases
  names(incid) = seq(from=as.Date(min(data_covid$date)), to=as.Date(max(data_covid$date)), by=1)
  
  # create generation time : gamma distribution, with mean 5 time units and standard deviation 1 time unit> 
  GT.cov = generation.time("gamma", c(5,1)) # assuming knowledge about infectious period
  res.R = estimate.R(incid, GT.cov, methods=c("EG"), pop.size = population_in ) # EG = exponential growth
  print(sprintf("R0 based on R0 Estimation Package in %s:", country))
  print(res.R)
  plot(res.R)
}



##################################################################################
# an R script to fit the parameters of a simple SIR model to influenza epidemic
# data by minimizing the Least Squares statistic
#
# Author:  Sherry Towers
#          smtowers@asu.edu
# Created: Dec 6, 2012
#
# Copyright Sherry Towers, 2012
#
# This script is not guaranteed to be free of bugs and/or errors
#
# This script can be freely used and shared as long as the author and
# copyright information in this header remain intact.
##################################################################################

