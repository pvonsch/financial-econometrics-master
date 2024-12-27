
# Install and load required packages
if (!require("tidyverse")) install.packages("tidyverse") # For data manipulation and visualization               
if (!require("psych")) install.packages("psych")         # For descriptive statistics
if(!requireNamespace("devtools")) install.packages("devtools") # Contains financial data, including DJ_d and DJ_w
devtools::install_github("yukai-yang/FE")

library(tidyverse)
library(FE)
library(psych)

#############################################################
## Part 1: Statistical Properties of Asset Returns
#############################################################

## A. Distributional properties of Dow Jones index returns

# 1. Plot the log returns

# Plot daily log returns
plot(DJ_d$r_Dow_Jones, type='l', main='Log returns of the daily Dow Jones index',
     ylab='Log return', xlab='Time horizon')

# ggplot alternative for daily log returns
DJ_d %>% ggplot() + geom_line(aes(y = r_Dow_Jones, x = 1:nrow(DJ_d))) +
  labs(x='Time horizon', y='Log return')

# Plot weekly log returns
plot(DJ_w$r_close, type='l', main='Log returns of the weekly Dow Jones index',
     ylab='Log return', xlab='Time horizon')

# Descriptive statistics for daily and weekly log returns
describe(DJ_d$r_Dow_Jones)
describe(DJ_w$r_close)

# 2. Evaluate empirical distributions using QQ-plots
fit_daily <- fitdistr(DJ_d$r_Dow_Jones, "t")
fit_daily$estimate  # Check estimated df for daily # value: 5.6824958249 

# Fit for weekly returns
fit_weekly <- fitdistr(DJ_w$r_close, "t")
fit_weekly$estimate  # Check estimated df for weekly # value: 3.515671117 

# Standardize the data
daily_standardized <- scale(DJ_d$r_Dow_Jones)
weekly_standardized <- scale(DJ_w$r_close)

# Compare daily log returns to normal and t-distributions
qqnorm(daily_standardized, main = "QQ Plot - Standardized Daily Log Returns vs Normal", ylab = "Sample Quantiles")
qqline(daily_standardized, col = "blue", lty = 2)  # Line for normal distribution

qqplot(rt(length(daily_standardized), df = 5.6824958249), daily_standardized,
       main = "QQ Plot - Standardized Daily Log Returns vs t(5.6824)", 
       xlab = "Theoretical Quantiles (t)", ylab = "Sample Quantiles")
qqline(daily_standardized, distribution = function(p) qt(p, df = 5.6824958249), col = "green", lty = 2)  # Line for t(5.6824)

# Compare weekly log returns to normal and t-distributions
qqnorm(weekly_standardized, main = "QQ Plot - Standardized Weekly Log Returns vs Normal", ylab = "Sample Quantiles")
qqline(weekly_standardized, col = "blue", lty = 2)  # Line for normal distribution

qqplot(rt(length(weekly_standardized), df = 3.515671117), weekly_standardized,
       main = "QQ Plot - Standardized Weekly Log Returns vs t(3.515)", 
       xlab = "Theoretical Quantiles (t)", ylab = "Sample Quantiles")
qqline(weekly_standardized, distribution = function(p) qt(p, df = 3.515671117), col = "green", lty = 2)  # Line for t(3.515)




# 3. Chi-square goodness-of-fit tests

# Normalize data
vdata = DJ_d$r_Dow_Jones
vdata = (vdata - mean(vdata)) / sd(vdata)

# Chi-square test for normality
ik = 20
grids = 1:ik / ik
vq = pnorm(vdata)
vn = NULL
for (val in grids) vn = c(vn, sum(vq <= val))
vn = c(vn[1], diff(vn))
test = sum((vn - length(vdata) / ik)^2 / (length(vdata) / ik))
cat("Normality test =", test, "df =", ik - 3, "p-value =", 1 - pchisq(test, df = ik - 3), "\n")

# Chi-square test for t-distribution
df = 3
ndata = vdata * sqrt(df / (df - 2))
vq = pt(ndata, df = 5)
vn = NULL
for (val in grids) vn = c(vn, sum(vq <= val))
vn = c(vn[1], diff(vn))
test = sum((vn - length(vdata) / ik)^2 / (length(vdata) / ik))
cat("t-distribution test =", test, "df =", ik - 3, "p-value =", 1 - pchisq(test, df = ik - 3), "\n")

# Chi-square test for mixture of normals
func <- function(vx) {
  alpha = vx[1]
  sigma = vx[2]
  ndata = vdata * sqrt((1 - alpha) + alpha * sigma^2)
  vq = (1 - alpha) * pnorm(ndata) + alpha * pnorm(ndata, sd = sigma)
  vn = NULL
  for (val in grids) vn = c(vn, sum(vq <= val))
  vn = c(vn[1], diff(vn))
  return(sum((vn - length(vdata) / ik)^2 / (length(vdata) / ik)))
}

optim(par = c(0.1, 4), fn = func, method = "BFGS")

#############################################################
## B. Dynamical properties of financial return series
#############################################################

# 1. Generate log returns and plot autocorrelations

lret = apply(log(index_d), 2, diff)
summary(lret)
matplot(lret, type = 'l', ylab = 'Log returns', xlab = 'Time horizon')

# Autocorrelation analysis for a specific index

par(mfrow = c (2,1))

tmp = 1
colnames(lret)[tmp]
acf(lret[!is.na(lret[, tmp]), tmp], main = "ACF for DAXINDX")
pacf(lret[!is.na(lret[, tmp]), tmp], main = "PACF for DAXINDX")
tmp = 2
colnames(lret)[tmp]
acf(lret[!is.na(lret[, tmp]), tmp], main = "ACF for FRCAC40")
pacf(lret[!is.na(lret[, tmp]), tmp], main = "PACF for FRCAC40")
tmp = 3
colnames(lret)[tmp]
acf(lret[!is.na(lret[, tmp]), tmp], main = "ACF for FTSE100")
pacf(lret[!is.na(lret[, tmp]), tmp], main = "PACF for FTSE100")
tmp = 4
colnames(lret)[tmp]
acf(lret[!is.na(lret[, tmp]), tmp], main = "ACF for HNGKNGI")
pacf(lret[!is.na(lret[, tmp]), tmp], main = "PACF for HNGKNI")
tmp = 5
colnames(lret)[tmp]
acf(lret[!is.na(lret[, tmp]), tmp], main = "ACF for NIKKEI")
pacf(lret[!is.na(lret[, tmp]), tmp], main = "PACF for NIKKEI")
tmp = 6
colnames(lret)[tmp]
acf(lret[!is.na(lret[, tmp]), tmp], main = "ACF for SNGALLS")
pacf(lret[!is.na(lret[, tmp]), tmp], main = "PACF for SNGALLS")
tmp = 7
colnames(lret)[tmp]
acf(lret[!is.na(lret[, tmp]), tmp], main = "ACF for SPCOMP")
pacf(lret[!is.na(lret[, tmp]), tmp], main = "PACF for SPCOMP")













# 2. Ljung-Box test

LB <- function(vx, lag, ip) {
  tmp = acf(vx, lag.max = lag, plot = FALSE)$acf
  tmp = tmp[2:(lag + 1)]^2
  test = sum(tmp / (length(vx) - 1:lag)) * length(vx) * (length(vx) + 2)
  return(list(test = test, pval = 1 - pchisq(test, df = lag - ip)))
}

LB(vx = lret[!is.na(lret[, tmp]), tmp], lag = 100, ip = 0)

tmp = 1
LB(vx = lret[!is.na(lret[, tmp]), tmp], lag = 100, ip = 0)
tmp = 2
LB(vx = lret[!is.na(lret[, tmp]), tmp], lag = 100, ip = 0)
tmp = 3
LB(vx = lret[!is.na(lret[, tmp]), tmp], lag = 100, ip = 0)
tmp = 4
LB(vx = lret[!is.na(lret[, tmp]), tmp], lag = 100, ip = 0)
tmp = 5
LB(vx = lret[!is.na(lret[, tmp]), tmp], lag = 100, ip = 0)
tmp = 6
LB(vx = lret[!is.na(lret[, tmp]), tmp], lag = 100, ip = 0)
tmp = 7
LB(vx = lret[!is.na(lret[, tmp]), tmp], lag = 100, ip = 0)
tmp = 1
LB(vx = lret[!is.na(lret[, tmp]), tmp], lag = 50, ip = 0)
tmp = 2
LB(vx = lret[!is.na(lret[, tmp]), tmp], lag = 50, ip = 0)
tmp = 3
LB(vx = lret[!is.na(lret[, tmp]), tmp], lag = 50, ip = 0)
tmp = 4
LB(vx = lret[!is.na(lret[, tmp]), tmp], lag = 50, ip = 0)
tmp = 5
LB(vx = lret[!is.na(lret[, tmp]), tmp], lag = 50, ip = 0)
tmp = 6
LB(vx = lret[!is.na(lret[, tmp]), tmp], lag = 50, ip = 0)
tmp = 7
LB(vx = lret[!is.na(lret[, tmp]), tmp], lag = 50, ip = 0)
tmp = 1
LB(vx = lret[!is.na(lret[, tmp]), tmp], lag = 10, ip = 0)
tmp = 2
LB(vx = lret[!is.na(lret[, tmp]), tmp], lag = 10, ip = 0)
tmp = 3
LB(vx = lret[!is.na(lret[, tmp]), tmp], lag = 10, ip = 0)
tmp = 4
LB(vx = lret[!is.na(lret[, tmp]), tmp], lag = 10, ip = 0)
tmp = 5
LB(vx = lret[!is.na(lret[, tmp]), tmp], lag = 10, ip = 0)
tmp = 6
LB(vx = lret[!is.na(lret[, tmp]), tmp], lag = 10, ip = 0)
tmp = 7
LB(vx = lret[!is.na(lret[, tmp]), tmp], lag = 10, ip = 0)



### 3. Cross-correlations

```{r crosscor}
corr_mat <- function(data, lag=1) {
  # column names
  cols <- colnames(data)
  # correlation matrix
  corr_mat <- matrix(nrow = ncol(data), ncol = ncol(data),
                     dimnames=list(lags=cols, cols))
  for (c1 in cols){
    for (c0 in cols){
      # periods with either are NA
      isna <- (is.na(data[, c1]) | is.na(data[, c0]))
      # non-na returns of pair
      pair_data <- data[!isna, c(c1, c0)]
      T <- nrow(pair_data)
      # add period lag autocorrelation to matrix
      corr_mat[c1, c0] <- cor(pair_data[1:(T-lag),c1], pair_data[(1+lag):T,c0])
      # acf(pair_data)
      # pacf(pair_data)
    }
  }
  # re-label past/current period returns
  # rownames(corr_mat) <- paste(cols, '(t-1)', sep='')
  # colnames(corr_mat) <- paste(cols, '(t)', sep='')
  return(corr_mat)
}

ncols = length(colnames(lret))

kbl(corr_mat(lret, lag=0), caption = "Cross-correlations with Period $t$",
    booktabs = T, digits = 3) %>%
  add_header_above(c("", "Period t"=ncols)) %>%
  kable_styling(latex_options = c("repeat_header"))

kbl(corr_mat(lret, lag=0), caption = "Cross-correlations with Period $t-1$",
    booktabs = T, digits = 3) %>%
  add_header_above(c("", "Period t"=ncols)) %>%
  kable_styling(latex_options = c("repeat_header"))

#creating acf-plots for each pair

# Number of time series
n <- ncol(lret)

# Set maximum lag
max_lags <- 30

# Creating directory to save plots
dir.create("CrossCorrelationPlots", showWarnings = FALSE)

# Loop through all unique pairs
for (i in 1:(n - 1)) {
  for (j in (i + 1):n) {
    # Select pair of time series
    pair <- c(i, j)
    x <- lret[, i]
    y <- lret[, j]
    
    # Remove missing values
    data <- na.omit(data.frame(x, y))
    
    # Skip if there's not enough data after removing NAs
    if (nrow(data) < max_lags) {
      cat(sprintf("Skipping pair (%d, %d) due to insufficient data\n", i, j))
      next
    }
    
    # Extract cleaned series
    x <- data$x
    y <- data$y
    
    # Generate and save cross-correlation plot
    plot_name <- sprintf("CrossCorrelation_Series%d_Series%d.png", i, j)
    png(file.path("CrossCorrelationPlots", plot_name))
    ccf(x, y, lag.max = max_lags, main = sprintf("Cross-correlation: Series %d vs Series %d", i, j))
    dev.off()
  }
}










# 4. Analyze squared returns

#Square log returns
lret2 = lret**2

tmp = 7; colnames(lret2)[tmp]; acf(lret2[!is.na(lret2[,tmp]),tmp])
summary(lret2)

#Generate and export stats table
TableB4 <- as.data.frame(describe(lret2))
print(xtable(TableB4, type = "latex", digits=10), file = "TableB4")

