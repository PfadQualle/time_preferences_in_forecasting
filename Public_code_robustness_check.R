## ---------------------------
##
## Script name: Public_code_empirical_analysis.R
##
## Purpose of script: Runs the empirical analysis for the paper "Time preferences in forecasting"
##
## Author: Niklas Valentin Lehmann
##
## Date Created: 2025-08-18
##
## 
## Email: Niklas-Valentin.Lehmann@vwl.tu-freiberg.de
##
## ---------------------------
##
## Notes: The code will not run because the Metaculus API undergoes continuous change and the data can no longer be pulled. 
##        Furthermore, the code is not written to be run top to bottom; some parts are optional (such as plots throughout).
##
## ---------------------------

### furthermore, this project includes a robustness check involving averaging within-event predictions. We kick out the bootstrap.

#requires: 

library(tidyverse)
library(broom)
library(forestplot)
library(lmtest)# for heteroskedasticity-consistent variance estimation
library(sandwich)# for heteroskedasticity-consistent variance estimation
library(metafor) # for random-effects meta analysis
library(httr2) # needed to make requests from Metaculus API



#-----

# DATA : 


# We just load the finished workflow data from the empirical analysis.


#-----

### Structure of document: 


#   We can cut most steps in the beginning because we already have the "finished data" - we start with step 5, where we collapse the forecasting data to exclude within-event dependence.


# STEP 5: ESTIMATE THE CALIBRATIONS FOR EACH SUBJECT (OLS/WLS APPROXIMATION)

# STEP 6: CALCULATE TREATMENT EFFECT FOR EACH SUBJECT (LEGACY)

# STEP 7: ANALYZE DIFFERENCES IN TREATMENT EFFECT IN THE POPULATION

# STEP 8: (OPTIONAL) PLOT CUMULATIVE CALIBRATION CURVES



###########################################################################
##################### # CHANGE TO SCRIPT : COLLAPSE UPDATED FORECASTS  #################
############################################################################

## The important dataframe is "forecast_binary_resolved_all", so we filter this one 

# Here we collapse the forecasts
forecast_binary_resolved_all <- forecast_binary_resolved_all |> group_by(author_id,question_id) |> # Make groups of forecasters 
  summarise(across(everything(), first),  # take all values from the first entry 
    .groups="drop")    



###########################################################################
##################### # STEP 5: ESTIMATE THE CALIBRATIONS FOR EACH SUBJECT (OLS/WLS APPROXIMATION)  #################
############################################################################


### For each subject, calculate the average effect of tvarying using a regression analysis

## Create a routine that creates a dataframe for every subject - necessary for the regression analysis and saves the regression results

preliminary_results <- data.frame()


# Counter 
for (i in 1:nrow(subjects)) {
  
  
  
  
  # Create dataframe with predictions on questions with fixed resolution
  temp_1 <- forecast_binary_resolved_all |>
    filter(author_id == subjects$author_id[i]) |> # take forecasts from one subject
    filter(tvarying == 0) #|> # Add the question information to each prediction
  #mutate(outcome = if_else(resolution == "yes", 1, 0))  #make a new variable "outcome" that is binary 
  #select(probability_yes,outcome,tvarying) # get rid of unnecessary variables
  
  
  # Create dataframe with predictions on questions with varying resolution
  temp_2 <- forecast_binary_resolved_all |>
    filter(author_id == subjects$author_id[i]) |> # take forecasts from one subject
    filter(tvarying == 1) #|> # Add the question information to each prediction
  #mutate(outcome = if_else(resolution == "yes", 1, 0))  #make a new variable "outcome" that is binary 
  #select(probability_yes,outcome,tvarying) # get rid of unnecessary variables
  
  

## Calculate frequencies for subject (make "plot")

# build the dataset for the calibration
calibration <- data.frame(prediction = double(),
                          frequency= double(),
                          tvarying = integer(),
                          weight = integer()) # initialize dataframe




for (k in 1:99){
  
  ##calculate frequency of events in dataset tfixed 
  
  temp_3 <- temp_1 |> 
    filter(probability_yes >= ((k-0.5)/100) & probability_yes < ((k+0.5)/100))
                      
  
  # Get the number of positive and negative outcomes conditional on prediction, weighted by IPTW
  temp_3 <- temp_3 |>
    mutate(weighted_outcome = outcome*IPTW) 

  frequency_tfixed <- as.numeric(sum(temp_3$weighted_outcome)/sum(temp_3$IPTW)) # divides the number of positive outcomes with the total number of outcomes
  
  #calculate frequency of events in dataset tvarying
  temp_4 <- temp_2 |> 
    filter(probability_yes >= ((k-0.5)/100) & probability_yes < ((k+0.5)/100))
  
  
  # Get the number of positive and negative outcomes conditional on prediction
  temp_4 <- temp_4 |>
    mutate(weighted_outcome = outcome*IPTW) 
  
  frequency_tvarying <- as.numeric(sum(temp_4$weighted_outcome)/sum(temp_4$IPTW)) # divides the number of positive outcomes with the total number of outcomes
  
  
  
  
  # Store the frequency and prediction for tfixed
  calibration <- calibration |> 
    add_row(prediction = k/100 , frequency = frequency_tfixed, tvarying=0,  weight = sum(temp_3$IPTW))
    
  # Store the frequency for tvarying 
  calibration <- calibration |> 
    add_row(prediction = k/100 , frequency = frequency_tvarying, tvarying= 1, weight = sum(temp_4$IPTW))
  
  
} 


#Sys.sleep(1) # reduce thermal load 

# estimate the calibration using WLS

lm <- lm(frequency ~ prediction + tvarying + prediction*tvarying, data = calibration, weights = weight)


# tidy results

lm_tidy <- tidy(lm)

# calculate heteroskedasticity-robust std. errors 

temp <- coeftest(lm, vcov = vcovHC(lm, type = "HC1")) # using McKinnon&White (1985) std. error calculation

#save heteroskedasticity-robust std. errors 

lm_tidy$std.error[1] <- temp[1,2]
lm_tidy$std.error[2] <- temp[2,2]
lm_tidy$std.error[3] <- temp[2,2]
lm_tidy$std.error[4] <- temp[2,2]

#kick t-values and p-values (which are no longer valid) 


lm_tidy <- lm_tidy |> 
  select(term, estimate, std.error)

## save results to data frame 


preliminary_results <- preliminary_results |>
  bind_rows(lm_tidy) 

#End routine
}




  
  ###########################################################################
  ##################### # STEP 7: ANALYZE DIFFERENCES IN TREATMENT EFFECT IN THE POPULATION  #################
  ############################################################################
  
  

# Re-name row 1 (Intercept) to Intercept for use in the analysis

#preliminary_results <- rename(preliminary_results, Intercept = "(Intercept)" )

# We need separate dataframes for each variable (regression coefficient) to plot them 

# Create a dataframe for intercepts
preliminary_results_intercepts <- tribble(
   ~estimate, ~conf.low, ~conf.high, ~std_err
)

  # Save intercepts to dataframe
for (i in 1:nrow(subjects)) {
  preliminary_results_intercepts <- preliminary_results_intercepts |>
    add_row(std_err = preliminary_results$std.error[1+4*(i-1)], estimate = preliminary_results$estimate[1+4*(i-1)], conf.low = preliminary_results$estimate[1+4*(i-1)] - 1.96*preliminary_results$std.error[1+4*(i-1)], conf.high = preliminary_results$estimate[1+4*(i-1)] + 1.96*preliminary_results$std.error[1+4*(i-1)])
}

# Create a dataframe for prediction
preliminary_results_probability_yes <- tribble(
  ~estimate, ~conf.low, ~conf.high, ~std_err
)

for (i in 1:nrow(subjects)) {
  preliminary_results_probability_yes <- preliminary_results_probability_yes |>
    add_row(std_err = preliminary_results$std.error[2+4*(i-1)], estimate = preliminary_results$estimate[2+4*(i-1)], conf.low = preliminary_results$estimate[2+4*(i-1)] - 1.96*preliminary_results$std.error[2+4*(i-1)], conf.high = preliminary_results$estimate[2+4*(i-1)] + 1.96*preliminary_results$std.error[2+4*(i-1)])
}



# Create a dataframe for tvarying
preliminary_results_tvarying <- tribble(
  ~estimate, ~conf.low, ~conf.high, ~std_err
)

for (i in 1:nrow(subjects)) {
  preliminary_results_tvarying <- preliminary_results_tvarying |>
    add_row(std_err = preliminary_results$std.error[3+4*(i-1)], estimate = preliminary_results$estimate[3+4*(i-1)], conf.low = preliminary_results$estimate[3+4*(i-1)] - 1.96*preliminary_results$std.error[3+4*(i-1)], conf.high = preliminary_results$estimate[3+4*(i-1)] + 1.96*preliminary_results$std.error[3+4*(i-1)])
}


# Create a dataframe for probability_yes*tvarying
preliminary_results_probability_yes_tvarying <- tribble(
  ~estimate, ~conf.low, ~conf.high, ~std_err
)

for (i in 1:nrow(subjects)) {
  preliminary_results_probability_yes_tvarying <- preliminary_results_probability_yes_tvarying |>
    add_row(std_err = preliminary_results$std.error[4*i], estimate = preliminary_results$estimate[4*i], conf.low = preliminary_results$estimate[4*i] - 1.96*preliminary_results$std.error[4*i], conf.high = preliminary_results$estimate[4*i] + 1.96*preliminary_results$std.error[4*i])
}



########## Calculate mean effect sizes


#### Calculate the mean and standard error of the average individuals coefficient.

### Calculate mean 

temp_1 <- rma(yi = preliminary_results_intercepts$estimate, sei = preliminary_results_intercepts$std_err, data = preliminary_results_intercepts, method = "DL")
temp_2 <- rma(yi = preliminary_results_probability_yes$estimate, sei = preliminary_results_probability_yes$std_err, data = preliminary_results_probability_yes, method = "DL") 
# the z-test needs to be adjusted to test against 1 - but the estimate is still significantly different from 1

# Run:  2-2*pnorm(1.0409, mean=1, sd=0.0097)

temp_3 <- rma(yi = preliminary_results_tvarying$estimate, sei = preliminary_results_tvarying$std_err, data = preliminary_results_tvarying, method = "DL")
temp_4 <- rma(yi = preliminary_results_probability_yes_tvarying$estimate, sei = preliminary_results_probability_yes_tvarying$std_err, data = preliminary_results_probability_yes_tvarying, method = "DL")


#summary(temp_1)



###### Make forestplots  #---------------------#

# 


# Add number of predictions to dataset


preliminary_results_intercepts <- preliminary_results_intercepts |> 
  bind_cols(subjects) 
  
#Re-order the subjects by most predictions to least predictions (Default would probably be old to new, but this is not too insightful.)
preliminary_results_intercepts <- preliminary_results_intercepts |> 
  arrange(desc(n_tfixed+n_tvarying))

# Kick out unnecessary variables again
preliminary_results_intercepts <- preliminary_results_intercepts |> 
  select(estimate,conf.high,conf.low,std_err,author_id) 



# Add mean effect to data to plot in forestplot
preliminary_results_intercepts <- preliminary_results_intercepts |>
  bind_rows(tibble(estimate = -0.0303, conf.low = -0.0390, conf.high = -0.0216))


# Generate the forest plot for intercepts # Save as A4 8.27 *11.69 inches
forestplot(
  labeltext= rep(NA,203),
  mean = preliminary_results_intercepts$estimate,
  lower = preliminary_results_intercepts$conf.low,
  upper = preliminary_results_intercepts$conf.high,
  zero = 0,  # This denotes the reference line at Risk Ratio = 1
  xlog = FALSE,  # Log-transform x-axis
  col = fpColors(box = "brown", lines = "brown1", zero = "black"),
  xticks = c( -0.5, 0, 0.5),
  is.summary = c(rep(FALSE, 202), TRUE),  # Make mean row more prominent
  #new_page = TRUE,
  grid= TRUE, 
)



# Add number of predictions to dataset
preliminary_results_probability_yes <- preliminary_results_probability_yes |> 
  bind_cols(subjects) 

#Re-order the subjects by most predictions to least predictions (Default would probably be old to new, but this is not too insightful.)
preliminary_results_probability_yes <- preliminary_results_probability_yes |>  
  arrange(desc(n_tfixed+n_tvarying)) 
  #select(estimate,conf.high,conf.low,std_err,author_id) # Kick out unnecessary variables again


# Add mean effect to data to plot in forestplot
preliminary_results_probability_yes <- preliminary_results_probability_yes |>
   bind_rows(tibble(estimate = 1.0409, conf.low = 1.0219, conf.high = 1.0599))


# Generate the forest plot for probability_yes
forestplot(
  labeltext= rep(NA,203),
  mean = preliminary_results_probability_yes$estimate,
  lower = preliminary_results_probability_yes$conf.low,
  upper = preliminary_results_probability_yes$conf.high,
  zero = 1,  # This denotes the reference line at Risk Ratio = 1
  xlog = FALSE,  # Log-transform x-axis
  col = fpColors(box = "brown", lines = "brown1", zero = "black"),
  xticks = c(0, 0.5, 1, 1.5, 2),
  is.summary = c(rep(FALSE, 202), TRUE),  # Make mean row more prominent
  #new_page = TRUE,
  grid= TRUE, 
)




# Add number of predictions to dataset
preliminary_results_tvarying <- preliminary_results_tvarying |> 
  bind_cols(subjects)
  

#Re-order the subjects by most predictions to least predictions (Default would probably be old to new, but this is not too insightful.)
  preliminary_results_tvarying <- preliminary_results_tvarying |> arrange(desc(n_tfixed+n_tvarying)) 

# Kick out unnecessary variables again
preliminary_results_tvarying <- preliminary_results_tvarying |> 
  select(estimate,conf.high,conf.low,std_err,author_id) 


# Add mean effect to data to plot in forestplot
preliminary_results_tvarying <- preliminary_results_tvarying |>
  bind_rows(tibble(estimate = -0.0129, conf.low = -0.0223, conf.high = -0.0035))



# Generate the forest plot for tvarying
forestplot(
  labeltext= rep(NA,203),
  mean = preliminary_results_tvarying$estimate,
  lower = preliminary_results_tvarying$conf.low,
  upper = preliminary_results_tvarying$conf.high,
  zero = 0,  # This denotes the reference line at Risk Ratio = 1
  xlog = FALSE,  # Log-transform x-axis
  col = fpColors(box = "darkblue", lines = "cadetblue3", zero = "black"),
  xticks = c( -0.5, 0, 0.5),
  is.summary = c(rep(FALSE, 202), TRUE),  # Make mean row more prominent
  #new_page = TRUE,
  grid= TRUE, 
)




# Add number of predictions to dataset
preliminary_results_probability_yes_tvarying <- preliminary_results_probability_yes_tvarying|> 
  bind_cols(subjects) 

#Re-order the subjects by most predictions to least predictions (Default would probably be old to new, but this is not too insightful.)
preliminary_results_probability_yes_tvarying <- preliminary_results_probability_yes_tvarying|> 
  arrange(desc(n_tfixed+n_tvarying)) 

# Kick out unnecessary variables again
preliminary_results_probability_yes_tvarying <- preliminary_results_probability_yes_tvarying |> 
  select(estimate,conf.high,conf.low,std_err,author_id) 


# Add mean effect to data to plot in forestplot
preliminary_results_probability_yes_tvarying <- preliminary_results_probability_yes_tvarying |>
  bind_rows(tibble(estimate = -0.1292, conf.low =  -0.1520, conf.high = -0.1064))




# Generate the forest plot for probability_yes:tvarying
forestplot(
  labeltext= rep(NA,203),
  mean = preliminary_results_probability_yes_tvarying$estimate,
  lower = preliminary_results_probability_yes_tvarying$conf.low,
  upper = preliminary_results_probability_yes_tvarying$conf.high,
  zero = 0,  # This denotes the reference line at Risk Ratio = 1
  xlog = FALSE,  # Log-transform x-axis
  col = fpColors(box = "darkblue", lines = "cadetblue3", zero = "black"),
  xticks = c(-1, -0.5, 0, 0.5, 1),
  is.summary = c(rep(FALSE, 202), TRUE),  # Make mean row more prominent
  #new_page = TRUE,
  grid= TRUE, 
  boxsize = 0.05/preliminary_results_probability_yes_tvarying$std_err
)




###########################################################################
##################### # STEP 8: PLOT CUMULATIVE CALIBRATION CURVES########## #################
############################################################################

########## 






cumulative_calibration <- data.frame(prediction = double(),
                                     frequency= double(),
                                     tvarying = integer()) # initialize dataframe

# Get a dataset with the predictions and outcomes - we take the prediction dataset and add outcomes





for (i in 1:99){# for tfixed
  
  #calculate frequency of events in dataset tfixed 
  temp <- forecast_binary_resolved_filtered_tfixed |># temp_1 and temp_2 are taken from STEP 5
    filter(probability_yes == i/100) |>
    inner_join(questions_binary_resolved_tfixed, by = join_by(question_id)) |>
    mutate(outcome = if_else(resolution == "yes", 1, 0)) 
    
  
  frequency <- as.numeric(count(temp, outcome)[2,2]/count(temp)) # divides the number of positive outcomes with the total number of outcomes
  
  
  # add row to the dataset for the cumulative calibration
  cumulative_calibration <- cumulative_calibration |> 
    add_row(prediction = i/100 , frequency = frequency, tvarying=0)
  
}
  

for (i in 1:99){# for tvarying
  

  #calculate frequency of events in dataset tvarying
  temp <- forecast_binary_resolved_filtered_tvarying |> 
    filter(probability_yes == i/100) |>
    inner_join(questions_binary_resolved_tvarying, by = join_by(question_id)) |>
    mutate(outcome = if_else(resolution == "yes", 1, 0)) 
  
  frequency <- as.numeric(count(temp, outcome)[2,2]/count(temp)) # divides the number of positive outcomes with the total number of outcomes
  
  # build the dataset for the cumulative calibration
  cumulative_calibration <- cumulative_calibration |> 
    add_row(prediction = i/100 , frequency = frequency, tvarying = 1)
  
}
  

for (i in 1:99){# for dynamic_minus (tvarying = -1)
  
  
  #calculate frequency of events in dataset dynamic_minus
  temp <- forecast_binary_resolved_filtered_dynamic_minus |> 
    filter(probability_yes == i/100) |>
    inner_join(questions_binary_resolved_dynamic_minus, by = join_by(question_id)) |>
    mutate(outcome = if_else(resolution == "yes", 1, 0)) 
  
  frequency <- as.numeric(count(temp, outcome)[2,2]/count(temp)) # divides the number of positive outcomes with the total number of outcomes
  
  
  # build the dataset for the cumulative calibration
  cumulative_calibration <- cumulative_calibration |> 
    add_row(prediction = i/100 , frequency = frequency, tvarying = -1)
}



### --------------------------------
### make the actual calibration plots 


## Plot all predictions and all states 

ggplot(cumulative_calibration, aes(x=prediction, y=frequency, color= as.factor(tvarying))) +
  geom_point() +
  theme_minimal() +
  scale_color_manual(values = c( "green2","red","cadetblue3"),
                    name = "Time-varying" ) +
  labs(,
       x = "Predictions",
       y = "Frequency")




## Plot all predictions and 2 states 


  
ggplot(filter(cumulative_calibration, tvarying != -1), aes(x=prediction, y=frequency, color= as.factor(tvarying))) +
  geom_point() +
  theme_minimal() +
  scale_color_manual(values = c( "red","cadetblue3"),
                     name = "Event type",
                     labels = c("Time-fixed", "Time-varying")) +
  labs(,
       x = "Predictions",
       y = "Frequency")



## Plot subjects predictions and 2 states  ( No change ) + OLS Analysis


# Calculate calibration anew



cumulative_calibration <- data.frame(prediction = double(),
                                     frequency= double(),
                                     tvarying = integer(), 
                                     count = integer()) # initialize dataframe

# Get a dataset with the predictions and outcomes - we take the prediction dataset and add outcomes





for (i in 1:99){# for tfixed
  
  #calculate frequency of events in dataset tfixed 
  temp <- forecast_binary_resolved_filtered_tfixed |># temp_1 and temp_2 are taken from STEP 5
    filter(probability_yes == i/100) |>
    filter(author_id %in% subjects$author_id) |>
    inner_join(questions_binary_resolved_tfixed, by = join_by(question_id)) |>
    mutate(outcome = if_else(resolution == "yes", 1, 0)) 
  
  
  frequency <- as.numeric(count(temp, outcome)[2,2]/count(temp)) # divides the number of positive outcomes with the total number of outcomes
  
  
  # add row to the dataset for the cumulative calibration
  cumulative_calibration <- cumulative_calibration |> 
    add_row(prediction = i/100 , frequency = frequency, tvarying=0, count=count(temp)[1,1])
  
}


for (i in 1:99){# for tvarying
  
  
  #calculate frequency of events in dataset tvarying
  temp <- forecast_binary_resolved_filtered_tvarying |> 
    filter(probability_yes == i/100) |>
    filter(author_id %in% subjects$author_id) |>
    inner_join(questions_binary_resolved_tvarying, by = join_by(question_id)) |>
    mutate(outcome = if_else(resolution == "yes", 1, 0)) 
  
  frequency <- as.numeric(count(temp, outcome)[2,2]/count(temp)) # divides the number of positive outcomes with the total number of outcomes
  
  # build the dataset for the cumulative calibration
  cumulative_calibration <- cumulative_calibration |> 
    add_row(prediction = i/100 , frequency = frequency, tvarying = 1, count=count(temp)[1,1])
  
}




ggplot(cumulative_calibration, aes(x=prediction, y=frequency, color= as.factor(tvarying))) +
  geom_point() +
  theme_minimal() +
  scale_color_manual(values = c( "red","cadetblue3"),
                     name = "Event type",
                     labels = c("Time-fixed", "Time-varying")) +
  labs(,
       x = "Predictions",
       y = "Frequency")


temp <- lm(frequency ~ prediction + tvarying + prediction*tvarying, data = cumulative_calibration)
# summary(temp)




############################################################################
################# Effects on accuracy #####################################


## the following script assesses the changes in Brier score as resulting from the difference in calibration or resolution of the forecasts on time-varying events and time-fixed events.


###---- Get weights over the state space 0.01, 0.02, ... 0.98, 0.99 for tvarying events (because these would change)

# from forecast_binary_resolved_filtered_tvarying make a plot 

# 
# ggplot(data = forecast_binary_resolved_filtered_tvarying, aes(x = probability_yes)) +
#   geom_histogram(aes(y = ..count..), bins = 999, fill = "lightblue", color = "black") 

# We take the dataframe 'cumulative calibration' as calculated in the last chunk of code prior to this section
# and use the weights/counts from tvarying only (see supplementary materials)


cumulative_calibration$count[1:99] <- cumulative_calibration$count[100:198]
  
  

##--- Calculate divergence between ideal calibration and actual frequency over the state space 


cumulative_calibration <- cumulative_calibration |>
  mutate(d = (prediction - frequency)^2)



## --- Sum up divergence 

# for tfixed 

temp_1 <- 0 #initialize
temp_2 <- 0 #initialize

for (k in c(1:99)) {
  temp_1 <- temp_1 + (cumulative_calibration$count[k]/sum(cumulative_calibration$count[1:99])) * cumulative_calibration$d[k]
}

# for tvarying

for (k in c(100:198)) {
  temp_2 <- temp_2 + (cumulative_calibration$count[k]/sum(cumulative_calibration$count[1:99])) * cumulative_calibration$d[k]
}

# the effect of time preferences on accuracy, as defined by the Brier score is temp_2 - temp_1 


