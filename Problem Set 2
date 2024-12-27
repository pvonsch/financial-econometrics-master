
#Preamble 

# Install and load required packages
if (!require("tidyverse")) install.packages("tidyverse") # For data manipulation and visualization               
if (!require("psych")) install.packages("psych")         # For descriptive statistics
if(!requireNamespace("devtools")) install.packages("devtools") # Contains financial data, including DJ_d and DJ_w
devtools::install_github("yukai-yang/FE")

library(tidyverse)
library(FE)
library(psych)


#Problem 1.1.1

acf(DJ_d$r_Dow_Jones, lag.max = 44, main = "ACF for daily returns")
acf(DJ_w$r_close, lag.max = 31, main = "ACF for weekly returns")

pacf(DJ_d$r_Dow_Jones, lag.max = 44, main = "PACF for daily returns")
pacf(DJ_w$r_close, lag.max = 231, main = "PACF for weekly returns")

#1.1.2

# Ljung-Box test for daily returns
for (i in 1:44) {
  print(paste("Lag", i, ":"))
  print(Box.test(DJ_d$r_Dow_Jones, lag = i, type = "Ljung-Box"))
}

# Ljung-Box test for weekly returns
for (i in 1:31) {
  print(paste("Lag", i, ":"))
  print(Box.test(DJ_w$r_close, lag = i, type = "Ljung-Box"))
}
#Problem 1.2
# Load necessary libraries
library(stats)

# Function to aggregate returns over n periods
aggregate_returns <- function(returns, n) {
  aggregated <- sapply(seq(1, length(returns), by = n), function(i) {
    if ((i + n - 1) <= length(returns)) {
      sum(returns[i:(i + n - 1)])
    } else {
      NA
    }
  })
  return(na.omit(aggregated))
}

# Aggregating returns
# Assuming DJ_d$r_Dow_Jones and DJ_w$r_close are the daily and weekly return series
two_day_returns <- aggregate_returns(DJ_d$r_Dow_Jones, 2)
two_week_returns <- aggregate_returns(DJ_w$r_close, 2)

# Aggregating higher periods (e.g., monthly and quarterly)
monthly_returns <- aggregate_returns(DJ_d$r_Dow_Jones, 21)  # Approx. 21 trading days
quarterly_returns <- aggregate_returns(DJ_w$r_close, 13)    # Approx. 13 weeks per quarter

# ACF and PACF for two-day returns
acf(two_day_returns, lag.max = 37, main = "ACF for Two-Day Returns")
pacf(two_day_returns, lag.max = 37, main = "PACF for Two-Day Returns")

# ACF and PACF for two-week returns
acf(two_week_returns, lag.max = 26, main = "ACF for Two-Week Returns")
pacf(two_week_returns, lag.max = 26, main = "PACF for Two-Week Returns")

# ACF and PACF for monthly returns
acf(monthly_returns, lag.max = 20, main = "ACF for Monthly Returns")
pacf(monthly_returns, lag.max = 20, main = "PACF for Monthly Returns")

# ACF and PACF for quarterly returns
acf(quarterly_returns, lag.max =16, main = "ACF for Quarterly Returns")
pacf(quarterly_returns, lag.max = 16, main = "PACF for Quarterly Returns")

# Variance Ratio Test Function
variance_ratio_test <- function(returns) {
  # Compute the variance of the aggregated returns and the variance of a random walk
  aggregated_var <- var(returns)
  random_walk_var <- var(diff(returns))
  ratio <- aggregated_var / random_walk_var
  return(ratio)
}

# Ljung-Box test function for selected lags (1, 11, 22, 44)
perform_ljung_box <- function(returns, lags, period_name) {
  results <- data.frame(Lag = integer(), TestStatistic = numeric(), PValue = numeric())
  for (lag in lags) {
    test <- Box.test(returns, lag = lag, type = "Ljung-Box")
    results <- rbind(results, data.frame(Lag = lag, TestStatistic = test$statistic, PValue = test$p.value))
  }
  cat(paste("\nLjung-Box Test Results for", period_name, "\n"))
  print(results)
  return(results)
}

# Perform Ljung-Box test for aggregated returns at specified lags (1, 11, 22, 44)
ljung_box_two_day <- perform_ljung_box(two_day_returns, c(1,11,22,37), "Two-Day Returns")
ljung_box_two_week <- perform_ljung_box(two_week_returns, c(1,11,22,26), "Two-Week Returns")
ljung_box_monthly <- perform_ljung_box(monthly_returns, c(1,11,20), "Monthly Returns")
ljung_box_quarterly <- perform_ljung_box(quarterly_returns, c(1,11,16), "Quarterly Returns")

# Calculate variance ratio for each aggregated return series
cat("\nVariance Ratio Test for Two-Day Returns:", variance_ratio_test(two_day_returns), "\n")
cat("Variance Ratio Test for Two-Week Returns:", variance_ratio_test(two_week_returns), "\n")
cat("Variance Ratio Test for Monthly Returns:", variance_ratio_test(monthly_returns), "\n")
cat("Variance Ratio Test for Quarterly Returns:", variance_ratio_test(quarterly_returns), "\n")


#1.3
#Dividing into sub samples
#Daily data

nD <- length(DJ_d$r_Dow_Jones)
sizeD <- floor (nD/3)

DailySubsample1 <- DJ_d$r_Dow_Jones[1:sizeD]
DailySubsample2 <- DJ_d$r_Dow_Jones[(sizeD + 1):(2 * sizeD)]
DailySubsample3 <- DJ_d$r_Dow_Jones[(2 * sizeD + 1):nD]

#Weekly data
nW <- length(DJ_w$r_close)
sizeW <- floor (nW/3)

WeeklySubsample1 <- DJ_w$r_close[1:sizeW]
WeeklySubsample2 <- DJ_w$r_close[(sizeW + 1):(2 * sizeW)]
WeeklySubsample3 <- DJ_w$r_close[(2 * sizeW + 1):nW]

#Ljung-Box on subsamples
#Daily subsamples
#1
for ( i in 1:34) {
  print ( paste ( "Lag", i, ":" ) )
  print ( Box.test (DailySubsample1, lag = i , type = "Ljung-Box"))
}
#2
for ( i in 1:34) {
  print ( paste ( "Lag", i, ":" ) )
  print ( Box.test (DailySubsample2, lag = i , type = "Ljung-Box"))
}
#3
for ( i in 1:34) {
  print ( paste ( "Lag", i, ":" ) )
  print ( Box.test (DailySubsample3, lag = i , type = "Ljung-Box"))
}

# Weekly subsamples
#1
for ( i in 1:24) {
  print ( paste ( "Lag", i, ":" ) )
  print ( Box.test (WeeklySubsample1, lag = i , type = "Ljung-Box"))
}
#2
for ( i in 1:24) {
  print ( paste ( "Lag", i, ":" ) )
  print ( Box.test (WeeklySubsample2, lag = i , type = "Ljung-Box"))
}
#3
for ( i in 1:24) {
  print ( paste ( "Lag", i, ":" ) )
  print ( Box.test (WeeklySubsample3, lag = i , type = "Ljung-Box"))
}
#Variance ratio test for the 6 subsamples

VDR(vr=DailySubsample1,iq=5)
VDR(vr=DailySubsample2,iq=5)
VDR(vr=DailySubsample3,iq=5)

VDR(vr=WeeklySubsample1,iq=5)
VDR(vr=WeeklySubsample2,iq=5)
VDR(vr=WeeklySubsample3,iq=5)


#Section B 

#ACF
# Generate ACF plots
for (i in 14:19) {
  portfolio_name <- colnames(portcap_m)[i]
  acf(portcap_m[[i]], main = paste("ACF for", portfolio_name), lag.max = 21)
}


#PACF 
# Generate PACF plots
for (i in 2:21) {
  portfolio_name <- colnames(portcap_m)[i]
  pacf(portcap_m[[i]], main = paste("PACF for", portfolio_name), lag.max = 21)
}

# instead summarzing all the acf pacf into a table




# Define the series to exclude
exclude_series <- c("Date", "Rf", "Market")

# Create a named list of portfolio time series by excluding the unwanted series
ts_list <- portcap_m %>%
  select(-all_of(exclude_series))  # Exclude Date, Rf, Market

# Verify the list
print(names(ts_list))




# Function to extract significant ACF lags
get_significant_acf <- function(ts_data, max_lag = 20, conf_level = 0.95) {
  # Compute ACF without plotting
  acf_res <- Acf(ts_data, plot = FALSE, lag.max = max_lag)
  
  # Calculate the confidence limit
  conf_limit <- qnorm((1 + conf_level)/2) / sqrt(length(ts_data))
  
  # Identify significant lags (exclude lag 0)
  sig_lags <- which(abs(acf_res$acf[-1]) > conf_limit)
  
  return(sig_lags)
}

# Function to extract significant PACF lags
get_significant_pacf <- function(ts_data, max_lag = 20, conf_level = 0.95) {
  # Compute PACF without plotting
  pacf_res <- Pacf(ts_data, plot = FALSE, lag.max = max_lag)
  
  # Calculate the confidence limit
  conf_limit <- qnorm((1 + conf_level)/2) / sqrt(length(ts_data))
  
  # Identify significant lags
  sig_lags <- which(abs(pacf_res$acf) > conf_limit)
  
  return(sig_lags)
}



# Initialize an empty data frame to store results
summary_table <- data.frame(
  Series = character(),
  Significant_ACF_Lags = character(),
  Significant_PACF_Lags = character(),
  stringsAsFactors = FALSE
)

# Define maximum lag based on your analysis needs
# For consistency, we'll use max_lag = 20 for all portfolios
max_lag_acf <- 21
max_lag_pacf <- 21

# Iterate over each portfolio and extract significant lags
for (series_name in names(ts_list)) {
  ts_data <- ts_list[[series_name]]
  
  # Extract significant ACF lags
  sig_acf <- get_significant_acf(ts_data, max_lag = max_lag_acf)
  
  # Extract significant PACF lags
  sig_pacf <- get_significant_pacf(ts_data, max_lag = max_lag_pacf)
  
  # Convert lags to comma-separated strings or 'None' if no significant lags
  sig_acf_str <- if(length(sig_acf) > 0) paste(sig_acf, collapse = ", ") else "None"
  sig_pacf_str <- if(length(sig_pacf) > 0) paste(sig_pacf, collapse = ", ") else "None"
  
  # Append to the summary table
  summary_table <- rbind(summary_table, data.frame(
    Series = series_name,
    Significant_ACF_Lags = sig_acf_str,
    Significant_PACF_Lags = sig_pacf_str,
    stringsAsFactors = FALSE
  ))
}

# View the summary table
print(summary_table)

# Ljung-Box test for each portfolio
for (i in 1:ncol(portcap_m)) {
  portfolio_name <- colnames(portcap_m)[i]
  test <- Box.test(portcap_m[, i], lag = 10, type = "Ljung-Box")
  cat(paste("Ljung-Box Test for", portfolio_name, "p-value:", test$p.value, "\n"))
}

#for more lags
# Load necessary library
library(xtable)

# Define lags to test
lags <- c(1, 5, 10, 21)

# Initialize a data frame to store results
results <- data.frame(
  Portfolio = rep(colnames(portcap_m), each = length(lags)),
  Lag = rep(lags, ncol(portcap_m)),
  P_Value = NA
)

# Fill the data frame with Ljung-Box test results
row_index <- 1
for (i in 1:ncol(portcap_m)) {
  for (lag in lags) {
    test <- Box.test(portcap_m[, i], lag = lag, type = "Ljung-Box")
    results$P_Value[row_index] <- test$p.value
    row_index <- row_index + 1
  }
}

# Round p-values for better readability
results$P_Value <- round(results$P_Value, 4)

# View results
print(results)



VDR <- function(vr,iq){
  iTT = length(vr)
  im = floor(iTT/iq)
  iT = im*iq
  
  rr = vr[1:iT]
  mu = mean(rr)
  sa2 = var(rr)
  
  arr = NULL
  for(iter in 1:(iT-iq+1))
    arr = c(arr,sum(rr[iter:(iter+iq-1)]))
  
  sc2 = sum((arr-mu*iq)**2)/iq/(iT-iq+1)/(1-(1/im))
  
  VD = sc2 - sa2
  VR = sc2/sa2
  tmp = sqrt(2*(2*iq-1)*(iq-1)/3/iq)
  
  VD = VD*sqrt(iT)/sa2/tmp
  VR = (VR-1)*sqrt(iT)/tmp
  
  p_VD = 2*(1-pnorm(abs(VD)))
  p_VR = 2*(1-pnorm(abs(VR)))
  
  return(list(VD=VD, VR=VR))
}

view(VDR)

VDR(portcap_m$Lo_30, iq=5)

# Initialize a list to store results
vdr_results <- list()

# Loop through columns 2 to 19
for (i in 2:19) {
  portfolio_name <- colnames(portcap_m)[i]  # Get the column name
  vdr_result <- VDR(portcap_m[[i]], iq = 5)  # Run the VDR function
  vdr_results[[portfolio_name]] <- vdr_result  # Store results in a list
}

# View results
print(vdr_results)




#B question 2, splitting the data set in two halfs to see satbility over time

# Determine the midpoint
mid_point <- floor(nrow(portcap_m) / 2)

# Split the dataset into two halves
first_half <- portcap_m[1:mid_point, ]
second_half <- portcap_m[(mid_point + 1):nrow(portcap_m), ]

# Verify the results
print(head(first_half))  # First few rows of the first half
print(head(second_half))  # First few rows of the second half


#ACF
# Generate ACF plots first half
for (i in 2:4) {
  portfolio_name <- colnames(portcap_m)[i]
  acf(first_half[[i]], main = paste("ACF for first half", portfolio_name), lag.max = 17)
}
#ACF
# Generate ACF plots second half
for (i in 2:4) {
  portfolio_name <- colnames(portcap_m)[i]
  acf(second_half[[i]], main = paste("ACF for second half", portfolio_name), lag.max = 17)
}

#PACF 
# Generate PACF plots for first half
for (i in 2:4) {
  portfolio_name <- colnames(portcap_m)[i]
  pacf(first_half[[i]], main = paste("PACF for first half", portfolio_name), lag.max = 17)
}

# Generate PACF plots for second half
for (i in 2:4) {
  portfolio_name <- colnames(portcap_m)[i]
  pacf(second_half[[i]], main = paste("PACF for second half", portfolio_name), lag.max = 17)
}


#Ljung-Box on Portfolio subsamples
#First period
#1
for ( i in 1:17) {
  print ( paste ( "Lag", i, ":" ) )
  print ( Box.test (first_half[[2]], lag = i , type = "Ljung-Box"))
}
#2
for ( i in 1:17) {
  print ( paste ( "Lag", i, ":" ) )
  print ( Box.test (first_half[[3]], lag = i , type = "Ljung-Box"))
}
#3
for ( i in 1:17) {
  print ( paste ( "Lag", i, ":" ) )
  print ( Box.test (first_half[[4]], lag = i , type = "Ljung-Box"))
}

# Second period

#1
for ( i in 1:17) {
  print ( paste ( "Lag", i, ":" ) )
  print ( Box.test (second_half[[2]], lag = i , type = "Ljung-Box"))
}
#2
for ( i in 1:17) {
  print ( paste ( "Lag", i, ":" ) )
  print ( Box.test (second_half[[3]], lag = i , type = "Ljung-Box"))
}
#3
for ( i in 1:17) {
  print ( paste ( "Lag", i, ":" ) )
  print ( Box.test (second_half[[4]], lag = i , type = "Ljung-Box"))
}

#VR-test on the 3 portfolios in the 2 periods
VDR(vr=first_half[[2]],iq=5)
VDR(vr=first_half[[3]],iq=5)
VDR(vr=first_half[[4]],iq=5)

VDR(vr=second_half[[2]],iq=5)
VDR(vr=second_half[[3]],iq=5)
VDR(vr=second_half[[4]],iq=5)

#B question 3


portcap_m <- as.data.frame(portcap_m)
# Test indexing
ncol(portcap_m)            # Ensure it returns the correct number of columns
colnames(portcap_m)        # Ensure it returns valid column names
portcap_m[[1]]             # Test accessing the first column

####Compute acf ####
# Plot ACF for each portfolio

lapply(seq_along(portcap_m), function(i) {
  portfolio_name <- colnames(portcap_m)[i]  # Get the column name
  acf(portcap_m[[i]], main = paste("ACF for", portfolio_name), lag.max = 21)
})


# Save plots to a PNG file
png("ACF_Plots3.png", width = 1600, height = 1200)  # Set larger dimensions
par(mfrow = c(3, 2))  # Layout for 20 plots

# Generate ACF plots
for (i in 14:19) {
  portfolio_name <- colnames(portcap_m)[i]
  acf(portcap_m[[i]], main = paste("ACF for", portfolio_name), lag.max = 21)
}

# Close the graphics device
dev.off()


####generate PACF####


# Save PACF plots to a PNG file
png("PACF_Plots3.png", width = 1600, height = 1200)  # Set larger dimensions
par(mfrow = c(3, 2))  # Layout for 20 plots (5 rows, 4 columns)

# Generate PACF plots
for (i in 14:19) {
  portfolio_name <- colnames(portcap_m)[i]
  pacf(portcap_m[[i]], main = paste("PACF for", portfolio_name), lag.max = 21)
}

# Close the graphics device
dev.off()



####Ljung Box ####
# Ljung-Box test for each portfolio
for (i in 1:ncol(portcap_m)) {
  portfolio_name <- colnames(portcap_m)[i]
  test <- Box.test(portcap_m[, i], lag = 10, type = "Ljung-Box")
  cat(paste("Ljung-Box Test for", portfolio_name, "p-value:", test$p.value, "\n"))
}

#for more lags
# Load necessary library
library(xtable)

# Define lags to test
lags <- c(1, 5, 10, 21)

# Initialize a data frame to store results
results <- data.frame(
  Portfolio = rep(colnames(portcap_m), each = length(lags)),
  Lag = rep(lags, ncol(portcap_m)),
  P_Value = NA
)

# Fill the data frame with Ljung-Box test results
row_index <- 1
for (i in 1:ncol(portcap_m)) {
  for (lag in lags) {
    test <- Box.test(portcap_m[, i], lag = lag, type = "Ljung-Box")
    results$P_Value[row_index] <- test$p.value
    row_index <- row_index + 1
  }
}

# Round p-values for better readability
results$P_Value <- round(results$P_Value, 4)

# View results
print(results)

library(tidyr)

# Pivot the data frame to wide format
results_wide <- results %>%
  pivot_wider(names_from = Lag, values_from = P_Value, names_prefix = "Lag_")

# Rename columns for clarity
colnames(results_wide) <- c("Portfolio", paste0("Lag ", lags))

# View the wide-format table
print(results_wide)

library(xtable)
# Create LaTeX table
latex_table <- xtable(results_wide, caption = "Ljung-Box Test Results", label = "tab:ljung_box")

# Export to LaTeX
print(latex_table, include.rownames = FALSE, file = "ljung_box_results.tex")



#New table
library(tidyr)  # For reshaping the table

# Define lags to test
lags <- c(1, 5, 10, 21)

# Initialize data frame to store results
results <- data.frame(
  Portfolio = rep(colnames(portcap_m), each = length(lags)),
  Lag = rep(lags, ncol(portcap_m)),
  Statistic = NA,
  P_Value = NA
)

# Compute Ljung-Box test statistics and p-values
row_index <- 1
for (i in 1:ncol(portcap_m)) {
  for (lag in lags) {
    test <- Box.test(portcap_m[, i], lag = lag, type = "Ljung-Box")
    results$Statistic[row_index] <- test$statistic
    results$P_Value[row_index] <- test$p.value
    row_index <- row_index + 1
  }
}

# Round for better readability
results$Statistic <- round(results$Statistic, 4)
results$P_Value <- round(results$P_Value, 4)

# View the results
print(results)


#now for latex
# Pivot to wide format for better presentation
results_wide <- results %>%
  pivot_wider(
    names_from = Lag,
    values_from = c(Statistic, P_Value),
    names_glue = "{.value}_Lag{Lag}"
  )

# View the reshaped table
print(results_wide)
library(xtable)
# Convert to LaTeX table
latex_table <- xtable(results_wide, 
                      caption = "Ljung-Box Test Statistics and P-Values", 
                      label = "tab:ljung_box")

# Export the table to a .tex file
print(latex_table, include.rownames = FALSE, file = "results_wide_latex.tex")

####variance ratio test ####


VDR <- function(vr,iq){
  iTT = length(vr)
  im = floor(iTT/iq)
  iT = im*iq
  
  rr = vr[1:iT]
  mu = mean(rr)
  sa2 = var(rr)
  
  arr = NULL
  for(iter in 1:(iT-iq+1))
    arr = c(arr,sum(rr[iter:(iter+iq-1)]))
  
  sc2 = sum((arr-mu*iq)**2)/iq/(iT-iq+1)/(1-(1/im))
  
  VD = sc2 - sa2
  VR = sc2/sa2
  tmp = sqrt(2*(2*iq-1)*(iq-1)/3/iq)
  
  VD = VD*sqrt(iT)/sa2/tmp
  VR = (VR-1)*sqrt(iT)/tmp
  
  p_VD = 2*(1-pnorm(abs(VD)))
  p_VR = 2*(1-pnorm(abs(VR)))
  
  return(list(VD=VD, VR=VR))
}

view(VDR)

VDR(portcap_m$Lo_30, iq=5)

# Initialize a list to store results
vdr_results <- list()

# Loop through columns 2 to 19
for (i in 2:19) {
  portfolio_name <- colnames(portcap_m)[i]  # Get the column name
  vdr_result <- VDR(portcap_m[[i]], iq = 5)  # Run the VDR function
  vdr_results[[portfolio_name]] <- vdr_result  # Store results in a list
}

# View results
print(vdr_results)


#Calculate the p-values for VD and VR
p_VD = 2*(1-pnorm(abs(VD)))
p_VR = 2*(1-pnorm(abs(VR)))

# Initialize a data frame to store results
vdr_results <- data.frame(
  Portfolio = character(),  # Portfolio names
  VD = numeric(),           # Variance difference
  VR = numeric(),           # Variance ratio
  p_VD = numeric(),         # p-value for VD
  p_VR = numeric()          # p-value for VR
)

# Loop through columns 2 to 19
for (i in 2:19) {
  portfolio_name <- colnames(portcap_m)[i]  # Get the column name
  vdr_result <- VDR(portcap_m[[i]], iq = 5)  # Run the VDR function
  
  # Append results to the data frame
  vdr_results <- rbind(vdr_results, data.frame(
    Portfolio = portfolio_name,
    VD = vdr_result$VD,
    VR = vdr_result$VR,
    p_VD = 2 * (1 - pnorm(abs(vdr_result$VD))),  # Compute p-value for VD
    p_VR = 2 * (1 - pnorm(abs(vdr_result$VR)))   # Compute p-value for VR
  ))
}

# View the results
print(vdr_results)


return(list(VD=VD, VR=VR, p_VD=p_VD, p_VR=p_VR))

library(xtable)

# Convert to LaTeX table
latex_table <- xtable(vdr_results, caption = "VDR Results with p-values", label = "tab:vdr_results")

# Export the LaTeX table
print(latex_table, include.rownames = FALSE, file = "vdr_results.tex")


#B question 4

#Q4
#Splitting high30 and low30
ccf(first_half[[2]], first_half[[4]], lag.max = 17, main= "Cross-autocorrelation Lo_30, Hi_30")
ccf(second_half[[2]], second_half[[4]], lag.max = 17, main= "Cross-autocorrelation Lo_30, Hi_30")

#LjungBoxing in the first period
#### Ljung-Box test for cross-correlations####
for (lag in c(1, 5, 10, 15, 17)) {
  cat("Lag:", lag, "\n")
  test <- Box.test(ccf(first_half[[2]], first_half[[4]], lag.max = lag)$acf[-1], lag = lag, type = "Ljung-Box")
  print(test)
}

# Define the lags to test
lags <- c(1, 5, 10, 15, 17)

# Initialize a data frame to store the results
ljung_box_results <- data.frame(
  Lag = lags,
  X_squared = numeric(length(lags)),
  df = numeric(length(lags)),
  p_value = numeric(length(lags))
)

# Perform the Ljung-Box test for cross-correlations
for (i in seq_along(lags)) {
  lag <- lags[i]
  # Compute the cross-correlation and Ljung-Box test
  ccf_vals <- ccf(first_half[[2]], first_half[[4]], lag.max = lag, plot = FALSE)$acf[-1]
  test <- Box.test(ccf_vals, lag = lag, type = "Ljung-Box")
  
  # Store results
  ljung_box_results$X_squared[i] <- test$statistic
  ljung_box_results$df[i] <- test$parameter
  ljung_box_results$p_value[i] <- test$p.value
}

# Round for better readability
ljung_box_results <- round(ljung_box_results, 4)

# View the results
print(ljung_box_results)


#export it
library(xtable)

#In period 2

#### Ljung-Box test for cross-correlations####
for (lag in c(1, 5, 10, 15, 17)) {
  cat("Lag:", lag, "\n")
  test <- Box.test(ccf(second_half[[2]], second_half[[4]], lag.max = lag)$acf[-1], lag = lag, type = "Ljung-Box")
  print(test)
}

# Define the lags to test
lags <- c(1, 5, 10, 15, 17)

# Initialize a data frame to store the results
ljung_box_results <- data.frame(
  Lag = lags,
  X_squared = numeric(length(lags)),
  df = numeric(length(lags)),
  p_value = numeric(length(lags))
)

# Perform the Ljung-Box test for cross-correlations
for (i in seq_along(lags)) {
  lag <- lags[i]
  # Compute the cross-correlation and Ljung-Box test
  ccf_vals <- ccf(second_half[[2]], second_half[[4]], lag.max = lag, plot = FALSE)$acf[-1]
  test <- Box.test(ccf_vals, lag = lag, type = "Ljung-Box")
  
  # Store results
  ljung_box_results$X_squared[i] <- test$statistic
  ljung_box_results$df[i] <- test$parameter
  ljung_box_results$p_value[i] <- test$p.value
}

# Round for better readability
ljung_box_results <- round(ljung_box_results, 4)

# View the results
print(ljung_box_results)


#export it
library(xtable)
