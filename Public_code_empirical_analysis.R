## ---------------------------
##
## Script name: Public_code_empirical_analysis.R
##
## Purpose of script: Runs the empirical analysis for the paper "Time preferences in forecasting"
##
## Author: Niklas Valentin Lehmann
##
## Date Created: 2025-06-13
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

# This script requires information on which questions induce a potential time preference bias, and which do not. 

# This includes the dataframe "questions_binary_resolved" - which lists questions, IDs and tvarying labels. 

# Furthermore, this script requires information on all Metaculus forecasts on these questions. 


#-----

### Structure of document: 

# STEP 1: IDENTIFY INDIVIDUALS FOR THE WITHIN-SUBJECT EXPERIMENT

# STEP 2: IDENTIFY DATA/QUESTIONS/EVENTS USED AS CONTROL TREATMENT 

# STEP 3: USE IPTW TO CREATE BALANCED SAMPLES 

# STEP 4: CHECK FOR BALANCE

# STEP 5: ESTIMATE THE CALIBRATIONS FOR EACH SUBJECT (OLS/WLS APPROXIMATION)

# STEP 6: CALCULATE TREATMENT EFFECT FOR EACH SUBJECT (LEGACY)

# STEP 7: ANALYZE DIFFERENCES IN TREATMENT EFFECT IN THE POPULATION

# STEP 8: (OPTIONAL) PLOT CUMULATIVE CALIBRATION CURVES


############
#     STEP 1: IDENTIFY INDIVIDUALS FOR THE WITHIN-SUBJECT EXPERIMENT
############


### we need to identify all forecasters who have made at least 30 predictions on each group (tvarying 0,1)

# lets first group the questions into tvarying 0,1

questions_binary_resolved_tvarying <- q_binary_res_100225 |>
  filter(tvarying == 1) 

questions_binary_resolved_tfixed <- q_binary_res_100225 |>
  filter(tvarying == 0)



# This dataset concerns questions with the opposing incentive (i.e. can only resolve "0/No" early)
questions_binary_resolved_dynamic_minus <- q_binary_res_100225 |>
  filter(tvarying == -1)



## The analysis of questions which can only resolve negatively early on is left out here - but can be added later (the group would not be very large)



### Now we subset the predictions such that there are only predictions on these questions

## filter predictions by predictions made on questions of interest


forecast_binary_resolved_filtered_tfixed <- forecast_binary |>
  filter(question_id %in% questions_binary_resolved_tfixed$question_id)


forecast_binary_resolved_filtered_tvarying <- forecast_binary |>
  filter(question_id %in% questions_binary_resolved_tvarying$question_id)


forecast_binary_resolved_filtered_dynamic_minus <- forecast_binary |>
  filter(question_id %in% questions_binary_resolved_dynamic_minus$question_id)



# Optional: get rid of unnecessary variables

forecast_binary_resolved_filtered_tfixed <- forecast_binary_resolved_filtered_tfixed |>
  select(question_id,author_id,probability_yes)

forecast_binary_resolved_filtered_tvarying <- forecast_binary_resolved_filtered_tvarying |>
  select(question_id,author_id,probability_yes)

# identify subjects that made 50+ predictions in EACH group

temp_1 <- forecast_binary_resolved_filtered_tfixed |>
  count(author_id) |>
  filter(n > 80) 
  
temp_2 <- forecast_binary_resolved_filtered_tvarying |>
  count(author_id) |>
  filter(n > 80) 


subjects <- filter(temp_1, author_id %in% temp_2$author_id) 


# It turns out that there are XXX individuals that made XX+ predictions for each group



# Save the number of predictions for the check for balance later

subjects <- left_join(subjects, temp_2, by = join_by(author_id)) 


subjects <- subjects |> 
  rename(n_tfixed = n.x) |>
  rename(n_tvarying = n.y)


### Filter for subjects that made predictions on at least 20 different questions in each group


for (i in c(1:nrow(subjects))) {
  
  temp_1 <- forecast_binary_resolved_filtered_tfixed |> #look at all forecasts tfixed
    filter(author_id == subjects$author_id[i]) |> # look at one subject
    distinct(question_id)  # how many different questions did the subject forecast on?

  subjects$tfixed_distinct_q[i] <- nrow(temp_1) # store the number of distinct questions the forecaster has forecasted on
  
  temp_1 <- forecast_binary_resolved_filtered_tvarying |> #look at all forecasts tvarying
    filter(author_id == subjects$author_id[i]) |> # look at one subject
    distinct(question_id)  # how many different questions did the subject forecast on?
  
  subjects$tvarying_distinct_q[i] <- nrow(temp_1) # store the number of distinct questions the forecaster has forecasted on
  
  i <- i+1
}


## Now filter out subjects who do not have predictions on at least 39 questions in each group

subjects <- subjects |> 
  filter(tvarying_distinct_q > 39) |>
  filter(tfixed_distinct_q > 39)
  
# Amazingly, 202 forecasters are left! 




###########################################################################
############### STEP 2: IDENTIFY DATA/QUESTIONS/EVENTS USED AS CONTROL TREATMENT  #################
############################################################################

# this part of the script is mostly a slightly changed version from the file control_topics_120325.R


## We first condition on fixed-event (tfixed) and dynamic events (tvarying)

save_question_id_tvarying  <- questions_binary_resolved_tvarying|> select(question_id) # save question_id vector (needed later)
save_question_id_tfixed  <- questions_binary_resolved_tfixed |> select(question_id) # save question_id vector (needed later)



## Prepare questions data for compatibility with API 

# Group questions topics cannot be accessed via their question_id, but a separate identifier
# We insert these identifiers manually as question_ids <- it is important to undo this later in order to not corrupt the question id's


questions_binary_resolved_tfixed$question_id[7:11] <- 15657
questions_binary_resolved_tfixed$question_id[20:23] <- 12551
questions_binary_resolved_tfixed$question_id[38:39] <- 21607
questions_binary_resolved_tfixed$question_id[49:56] <- 14484
questions_binary_resolved_tfixed$question_id[57:60] <- 13442
questions_binary_resolved_tfixed$question_id[63:93] <- 19494
questions_binary_resolved_tfixed$question_id[94:98] <- 17326
questions_binary_resolved_tfixed$question_id[99:102] <- 20920
questions_binary_resolved_tfixed$question_id[103:133] <- 22313
questions_binary_resolved_tfixed$question_id[134:137] <- 11590
questions_binary_resolved_tfixed$question_id[138:139] <- 10963
questions_binary_resolved_tfixed$question_id[140:151] <- 11654
questions_binary_resolved_tfixed$question_id[152:160] <- 15767
questions_binary_resolved_tfixed$question_id[161:164] <- 15755
questions_binary_resolved_tfixed$question_id[165:171] <- 19093
questions_binary_resolved_tfixed$question_id[165:171] <- 19093
questions_binary_resolved_tfixed$question_id[172:178] <- 15931
questions_binary_resolved_tfixed$question_id[179:185] <- 15923
questions_binary_resolved_tfixed$question_id[186:217] <- 19594
questions_binary_resolved_tfixed$question_id[218:220] <- 10988
questions_binary_resolved_tfixed$question_id[221:231] <- 19039
questions_binary_resolved_tfixed$question_id[232:248] <- 11370
questions_binary_resolved_tfixed$question_id[249:260] <- 14236
questions_binary_resolved_tfixed$question_id[261:268] <- 14236
questions_binary_resolved_tfixed$question_id[269:275] <- 20000
questions_binary_resolved_tfixed$question_id[276:291] <- 19115
questions_binary_resolved_tfixed$question_id[292:299] <- 18393
questions_binary_resolved_tfixed$question_id[300:307] <- 11087
questions_binary_resolved_tfixed$question_id[308:312] <- 11809
questions_binary_resolved_tfixed$question_id[313:316] <- 11921
questions_binary_resolved_tfixed$question_id[317:331] <- 12436
questions_binary_resolved_tfixed$question_id[332:335] <- 14493
questions_binary_resolved_tfixed$question_id[336:355] <- 15036
questions_binary_resolved_tfixed$question_id[356:359] <- 18596
questions_binary_resolved_tfixed$question_id[361:376] <- 11358
questions_binary_resolved_tfixed$question_id[407:408] <- 20274
questions_binary_resolved_tfixed$question_id[592:596] <- 25148
questions_binary_resolved_tfixed$question_id[701:702] <- 11780
questions_binary_resolved_tfixed$question_id[718:723] <- 19166
questions_binary_resolved_tfixed$question_id[724:743] <- 19178
questions_binary_resolved_tfixed$question_id[744:756] <- 14369
questions_binary_resolved_tfixed$question_id[757:759] <- 15140
questions_binary_resolved_tfixed$question_id[770:775] <- 12694

questions_binary_resolved_tvarying$question_id[137:144] <- 10860
questions_binary_resolved_tvarying$question_id[169]     <- 19393
questions_binary_resolved_tvarying$question_id[321]     <- 20106
questions_binary_resolved_tvarying$question_id[322:324] <- 20044
questions_binary_resolved_tvarying$question_id[327:330] <- 25088
questions_binary_resolved_tvarying$question_id[333:334] <- 18178
questions_binary_resolved_tvarying$question_id[365:366] <- 14006
questions_binary_resolved_tvarying$question_id[365:366] <- 14006
questions_binary_resolved_tvarying$question_id[369:370] <- 13125
questions_binary_resolved_tvarying$question_id[439]     <- 20212
questions_binary_resolved_tvarying$question_id[468:469] <- 13522
questions_binary_resolved_tvarying$question_id[486]     <- 15426
questions_binary_resolved_tvarying$question_id[550:554] <- 20048
questions_binary_resolved_tvarying$question_id[657]     <- 17472
questions_binary_resolved_tvarying$question_id[760:764] <- 18008




########### ADD TOPICS USING THE METACULUS API 



# We collect topics and tags from the metaculus API and add them to the QUESTIONS DATA

##### first for fixed questions, then dynamic questions 


for(q in 1:nrow(questions_binary_resolved_tfixed)){# For each question, we get topics and tags and store them in the questions data
  
  ## Get data from Metaculus API 
  
  temp_1 <- request("https://www.metaculus.com") # build request 
  temp_1 <- temp_1 |> req_url_path(paste("/api/posts/",questions_binary_resolved_tfixed$question_id[q],sep="")) # Adds the necessary path and takes correct question_id value from the data
  
  temp_2 <- temp_1 |> req_perform() # get the data from metaculus API
  
  temp_3 <- temp_2 |> resp_body_json() # Save the data from metaculus API 
  
  ### Process the data 
  
  # get projects
  
  for(i in length(temp_3$projects$category)){ #get every project one by one
    
    temp_topic <- temp_3$projects$category[[i]]$name #save category/topic that 
    
    #CHECK WHETHER A TAG IS ASSIGNED AT ALL
    if(is.null(temp_topic) == TRUE){#do nothing
    } else { # continue
      
      # Is the topic already in the table ?
      if(temp_topic  %in% colnames(questions_binary_resolved_tfixed)){ #if yes, then simply annotate this
        
        questions_binary_resolved_tfixed[q,temp_topic] <- 1 #and annotate the question (has to be done in base R)
        
      } else{ # if not, add the variable - and annotate the question
        
        questions_binary_resolved_tfixed <- questions_binary_resolved_tfixed |>
          mutate(!!temp_topic := 0)#add the variable (unquoted and define is necessary )
        
        questions_binary_resolved_tfixed[q,temp_topic] <- 1 #and annotate the question (has to be done in base R)
      }
      
    }
  }
  
  # get tags
  
  for(i in length(temp_3$projects$tag)){ #get every project one by one
    
    temp_topic <- temp_3$projects$tag[[i]]$name #save category/topic that 
    
    #CHECK WHETHER A TAG IS ASSIGNED AT ALL
    if(is.null(temp_topic) == TRUE){#do nothing
    } else { # continue
      
      
      # Is the topic already in the table ?
      if(temp_topic  %in% colnames(questions_binary_resolved_tfixed)){ #if yes, then simply annotate this
        
        questions_binary_resolved_tfixed[q,temp_topic] <- 1 #and annotate the question (has to be done in base R)
        
      } else{ # if not, add the variable - and annotate the question
        
        questions_binary_resolved_tfixed <- questions_binary_resolved_tfixed |>
          mutate(!!temp_topic := 0)#add the variable (unquoted and define is necessary )
        
        questions_binary_resolved_tfixed[q,temp_topic] <- 1 #and annotate the question (has to be done in base R)
      }
      
    }
  }
  
} # End for-loop


##### Now dynamic questions 


for(q in 1:nrow(questions_binary_resolved_tvarying)){# For each question, we get topics and tags and store them in the questions data
  
  ## Get data from Metaculus API 
  
  temp_1 <- request("https://www.metaculus.com") # build request 
  temp_1 <- temp_1 |> req_url_path(paste("/api/posts/",questions_binary_resolved_tvarying$question_id[q],sep="")) # Adds the necessary path and takes correct question_id value from the data
  
  temp_2 <- temp_1 |> req_perform() # get the data from metaculus API
  
  temp_3 <- temp_2 |> resp_body_json() # Save the data from metaculus API 
  
  ### Process the data 
  
  # get projects
  
  for(i in length(temp_3$projects$category)){ #get every project one by one
    
    temp_topic <- temp_3$projects$category[[i]]$name #save category/topic that 
    
    #CHECK WHETHER A TAG IS ASSIGNED AT ALL
    if(is.null(temp_topic) == TRUE){#do nothing
    } else { # continue
      
      # Is the topic already in the table ?
      if(temp_topic  %in% colnames(questions_binary_resolved_tvarying)){ #if yes, then simply annotate this
        
        questions_binary_resolved_tvarying[q,temp_topic] <- 1 #and annotate the question (has to be done in base R)
        
      } else{ # if not, add the variable - and annotate the question
        
        questions_binary_resolved_tvarying <- questions_binary_resolved_tvarying |>
          mutate(!!temp_topic := 0)#add the variable (unquoted and define is necessary )
        
        questions_binary_resolved_tvarying[q,temp_topic] <- 1 #and annotate the question (has to be done in base R)
      }
      
    }
  }
  
  # get tags
  
  for(i in length(temp_3$projects$tag)){ #get every project one by one
    
    temp_topic <- temp_3$projects$tag[[i]]$name #save category/topic that 
    
    #CHECK WHETHER A TAG IS ASSIGNED AT ALL
    if(is.null(temp_topic) == TRUE){#do nothing
    } else { # continue
      
      
      # Is the topic already in the table ?
      if(temp_topic  %in% colnames(questions_binary_resolved_tvarying)){ #if yes, then simply annotate this
        
        questions_binary_resolved_tvarying[q,temp_topic] <- 1 #and annotate the question (has to be done in base R)
        
      } else{ # if not, add the variable - and annotate the question
        
        questions_binary_resolved_tvarying <- questions_binary_resolved_tvarying |>
          mutate(!!temp_topic := 0)#add the variable (unquoted and define is necessary )
        
        questions_binary_resolved_tvarying[q,temp_topic] <- 1 #and annotate the question (has to be done in base R)
      }
      
    }
  }
  
} # End for-loop


###



###  undo question_id changes

questions_binary_resolved_tfixed <- questions_binary_resolved_tfixed |> 
  select(1,3:302) |> #kick out question ids
  add_column(save_question_id_tfixed, .before = 2) #add question_ids



questions_binary_resolved_tvarying <- questions_binary_resolved_tvarying |> 
  select(1,3:302) |> #kick out question ids
  add_column(save_question_id_tvarying, .before = 2) #add question_ids




###########################################################################
############### STEP 3: USE IPTW TO CREATE BALANCED SAMPLES  #################
############################################################################

### Identify covariates/topics to control for 

## Which topics come up at least 5 times across the entire sample? 

## --- Get a dataset with topics and tags going 

# Get the number of tags and topics in fixed
temp_1 <- questions_binary_resolved_tfixed |> 
  summarise_all(function(x){sum(x==1)}) |> # this function simply counts all rows across all columns
  pivot_longer(everything())

# Get the number of tags and topics in dynamic and add to fixed
temp_2 <- questions_binary_resolved_tvarying |> 
  summarise_all(function(x){sum(x==1)}) |> # this function simply counts all rows across all columns
  pivot_longer(everything()) |>
  full_join(temp_1, by= join_by(name))  # Add fixed question topics

temp_2[is.na(temp_2)] <- 0 # Kick NAs to be able to sum across datasets

temp_2 <- temp_2 |> 
  mutate(sum= value.x + value.y) |> # sum up across datasets
  filter(sum >= 5) |> # filter for topics that come up at least 5 times
  select(name,sum) |>
  pivot_wider(names_from=name, values_from = sum)
  

# 80 topics left, pretty neat!

### --- Get a dataset with all forecasts going


forecast_binary_resolved_all <- full_join(forecast_binary_resolved_filtered_tfixed,forecast_binary_resolved_filtered_tvarying)

questions_binary_resolved_all <- full_join(questions_binary_resolved_tfixed,questions_binary_resolved_tvarying)

## Add the topics and tvarying data


forecast_binary_resolved_all <- forecast_binary_resolved_all |>
  filter(author_id %in% subjects$author_id) |> # take forecasts from subjects only
  left_join(questions_binary_resolved_all, by = "question_id") |> # Add the question information to each prediction
  mutate(outcome = if_else(resolution == "yes", 1, 0)) 


#make tvarying numeric
forecast_binary_resolved_all$tvarying <- as.numeric(forecast_binary_resolved_all$tvarying) 
  
# get rid of NAs
forecast_binary_resolved_all[is.na(forecast_binary_resolved_all)] <- 0


### Estimate the propensity score for each question by running a logit regression


log_model <- glm(tvarying ~ Technology + `Social Sciences` + `Economy & Business` + `S&P 500 Index` + `Natural Sciences`
                 + California + Bioethics + `Computing and Math` + Politics + Russia + Geopolitics + `Health & Pandemics` +
                   Epidemiology + SpaceX + Cryptocurrencies + Medicine + China + Ukraine + `US Senate` + Elections + 
                   `Artificial Intelligence` + `Joe Biden` + Virology + `Elon Musk` + `Donald Trump` + `President of the US` +
                   DeepMind + `Google Scholar` + `Sports & Entertainment` + LIGO + `United States` +  `North Korea` +
                   `Planetary science` + Space + `Tesla Model 3` + `Bureau of Labor Statistics` + Law + `Aerospace engineering` +
                   `Environment & Climate` + Microsoft + `Benjamin Netanyahu` + `Boris Johnson` + Twitter + Chatbot + `Silicon Valley Bank` +
                   SEC + OpenAI + `Hunter Biden` + `New York (state)` + `PredictIt` + Tesla + `Washington, D.C.` + `Republican Party (US)` +
                   `Australia women's national soccer team` + `Climate change denial` + `Democratic Left Alliance` + 
                   `FIDE world rankings` + `Saudi Professional League` + `Yesh Atid` + `Kukiz'15` + `Pieter Omtzigt` + 
                   `European United Left–Nordic Green Left` + `Roy Jenkins` + `S.C. Braga` + `Phil Murphy` + `Donald Trump Jr.` +
                   `Ilya Sutskever` + `Ecuador` + `Richárd Rapport` + `Roman Baber` + `Param Singh` + `Meral Akşener` +
                   `Liz Truss` + `Nate Silver` + `African National Congress` + `Jake Sherman` + `Brian Fitzpatrick` + `Patrick McHenry` +
                   `List of Major League Baseball progressive single-season home run leaders`, data = forecast_binary_resolved_all, family = "binomial") 

#Note: We leave the 2024 Leaderboard Topic, because it is probably useless for the analysis

#summary(log_model)







### -- Begin function that calculates (treatment) effects beta_0 through beta_3 ---


## Initialize dataframe for final results

bootstrap_results <- tribble(
  ~beta_0, ~beta_1, ~beta_2, ~beta_3, ~std_err_beta_0, ~std_err_beta_1, ~std_err_beta_2, ~std_err_beta_3
)


##### FUNCTION to resample IPTW 


resample_log_model <- function(){ # We simply resample the log-model randomly to produce a new log_model that we use in the meta-analysis

  # Sample new coefficients from standard normal distribution
  temp_1 <- as_tibble(log_model$coefficients) |> # make into dataframe
    add_column(sqrt(diag(vcov(log_model)))) |>
    rename(coefficients = value) |>
    rename(std_err = "sqrt(diag(vcov(log_model)))") |>
    rowwise() |>
    mutate(resampled_coefficients = coefficients + rnorm(1,mean=0, sd = std_err)) |> # resample coefficient estimates
    mutate(resampled_coefficients = ifelse(coefficients >= 3, 100, resampled_coefficients)) |> # manually set propensity scores of determinants of fixed questions to 100 
    mutate(resampled_coefficients = ifelse(coefficients <= -3, -100, resampled_coefficients))   # manually set propensity scores of determinants of dynamic questions to -100
    
  # Return a new model object with tweaked coefficients
  log_model$coefficients <- temp_1$resampled_coefficients
  return(log_model)
}
#####



###
#BEGIN BOOTSTRAP 
###

set.seed(142857) # set seed for replicability


for (boot in c(1:30)) {

  ### Resample propensity scores
  
resampled_log_model <- resample_log_model() ## We keep "log_model" as the "original" best estimate of propensity scores

### Calculate inverse probability to be treated weights 

# For each forecast in the dataset all  - calculate the propensity score

forecast_binary_resolved_all <- forecast_binary_resolved_all |>
  mutate(propensity_score = predict(resampled_log_model, select(forecast_binary_resolved_all,9:524), type = "response")) ## calculate probability to be a prediction on a dynamic question

#make propensity score numeric

forecast_binary_resolved_all$propensity_score <- as.numeric(forecast_binary_resolved_all$propensity_score)

forecast_binary_resolved_all <- forecast_binary_resolved_all |> 
  mutate(IPTW = if_else(tvarying == 1,1/propensity_score,1/(1-propensity_score))) |> # calculate weights
  relocate(IPTW,.after=tvarying) |>
  relocate(propensity_score,.after=tvarying) 


###
###
# We cut out step 4 (Check for balance) in the bootstrap version
###
###

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
    filter(probability_yes == k/100)
                      
  
  # Get the number of positive and negative outcomes conditional on prediction, weighted by IPTW
  temp_3 <- temp_3 |>
    mutate(weighted_outcome = outcome*IPTW) 

  frequency_tfixed <- as.numeric(sum(temp_3$weighted_outcome)/sum(temp_3$IPTW)) # divides the number of positive outcomes with the total number of outcomes
  
  #calculate frequency of events in dataset tvarying
  temp_4 <- temp_2 |> 
    filter(probability_yes == k/100)
  
  
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


Sys.sleep(1) # reduce thermal load 

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





### Return 

bootstrap_results <- bootstrap_results |> 
  add_row(beta_0 = temp_1$beta, beta_1 = temp_2$beta, beta_2 = temp_3$beta, beta_3 = temp_4$beta, std_err_beta_0 = temp_1$se, std_err_beta_1 = temp_2$se, std_err_beta_2 = temp_3$se, std_err_beta_3 = temp_4$se)
}
## ---------- Here we end the bootstrapping function


########
####### Collect results 


set.seed(428571)


beta_0 <- bootstrap_results |>
  rowwise() |>
  mutate(beta_0_dist = list(rnorm(1000, mean =beta_0, sd = std_err_beta_0))) |>
  unnest(beta_0_dist) |>
  select(beta_0_dist)

beta_1 <- bootstrap_results |>
  rowwise() |>
  mutate(beta_1_dist = list(rnorm(1000, mean =beta_1, sd = std_err_beta_1))) |>
  unnest(beta_1_dist) |>
  select(beta_1_dist)

beta_2 <- bootstrap_results |>
  rowwise() |>
  mutate(beta_2_dist = list(rnorm(1000, mean =beta_2, sd = std_err_beta_2))) |>
  unnest(beta_2_dist) |>
  select(beta_2_dist)

beta_3 <- bootstrap_results |>
  rowwise() |>
  mutate(beta_3_dist = list(rnorm(1000, mean =beta_3, sd = std_err_beta_3))) |>
  unnest(beta_3_dist) |>
  select(beta_3_dist)

## To get results

# mean(beta_0$beta_0_dist)
# sd(beta_0)



# Plot each result

ggplot(beta_3, aes(x = beta_3_dist)) +
  geom_histogram(aes(y = ..density..), bins = 50, fill = "lightblue", color = "black") +
  geom_density(color = "red3", linewidth = 1) +
  ggtitle("Mixture of Normal Distributions") +
  xlab("Coefficient estimate")



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
  arrange(desc(n_tfixed+n_tvarying)) +
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


