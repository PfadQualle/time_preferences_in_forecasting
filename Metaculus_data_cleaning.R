## ---------------------------
##
## Script name: R script metaculus data_modified.R
##
## Purpose of script: Takes in Metaculus data and cleans it for the analysis "Time preferences in forecasting"
##
## Author: Niklas Valentin Lehmann
##
## Date Created: 2026-06-02
##
## 
## Email: Niklas-Valentin.Lehmann@vwl.tu-freiberg.de
##
## ---------------------------
##
## Notes: This code is used to clean the Metaculus data, and has mostly been supplied by Metaculus staff. 
##          The author takes full responsibility for mistakes and errors.
##
## ---------------------------

### furthermore, this project includes a robustness check involving averaging within-event predictions. We kick out the bootstrap.

#requires: 

library(data.table)
library(magrittr)
library(lubridate)

folder <- "PATH_TO_FILES" #Instead of specifying a path, I just loaded the data in the same working directory.

## Read in the data ------------------------------------------------------------
q <- fread("questions.csv", header = TRUE, sep = ",")
setnames(q, old = "id", new = "question_id")
f <- fread("forecasts.csv", header = TRUE, sep = ",")


## Separating binary and multiple choice data ----------------------------------
binary_question_ids <- q[type == "binary", question_id]
mc_question_ids <- q[type == "multiple_choice", question_id]

questions_binary <- q[question_id %in% binary_question_ids]
forecast_binary <- f[question_id %in% binary_question_ids]

questions_mc <- q[question_id %in% mc_question_ids]
forecast_mc <- f[question_id %in% mc_question_ids]



## Extracting the continuous_cdf column for continuous forecasts ---------------
# convert continous_cdf column from string to a list of numerical values
numeric_question_ids <- q[type == "numeric", question_id]
date_question_ids <- q[type == "date", question_id]
f[question_id %in% c(numeric_question_ids, date_question_ids),
  continuous_cdf := continuous_cdf %>%
    gsub("\\{", "", .) %>%
    gsub("\\}", "", .) |>
    strsplit(",") |>
    lapply(as.numeric)
    ]

## Define the grid of x-values corresponding to the cdf ------------------------
# Helper functions
# These helper functions convert from the internal presentation of x-values
# to the actual values required for data analysis.
zero_point_to_deriv_ratio <- function(zero_point, lower_bound, upper_bound) {
  return ((upper_bound - zero_point) / (lower_bound - zero_point))
}
internal_to_actual <- function(x, zero_point, lower_bound, upper_bound) {
  actual_range <- upper_bound - lower_bound
  # if zero_point is not specified, return the linear transformation
  if (is.na(zero_point)) {
    return (lower_bound + actual_range * x)
  }
  dr <- zero_point_to_deriv_ratio(zero_point, lower_bound, upper_bound)
  if (dr == 1) {
    return (lower_bound + actual_range * x)
  } else {
    return (lower_bound + actual_range * (dr^x - 1) / (dr - 1))
  }
}

# apply helper function to compute actual grid of x values
q[question_id %in% c(numeric_question_ids, date_question_ids),
  x_grid := list(list(internal_to_actual(
    x = seq(0, 1, length.out = 201),
    zero_point = zero_point,
    lower_bound = range_min, upper_bound = range_max
  ))), by = question_id]


## Separating numeric and date forecasts ---------------------------------------
questions_numeric <- q[question_id %in% numeric_question_ids]
forecast_numeric <- f[question_id %in% numeric_question_ids]

questions_date <- q[question_id %in% date_question_ids]
questions_date[, xgrid := lapply(x_grid, lubridate::as_datetime)]
forecast_date <- f[question_id %in% date_question_ids]
