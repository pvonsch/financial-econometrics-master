
#Preamble 

# Install and load required packages
if (!require("tidyverse")) install.packages("tidyverse") # For data manipulation and visualization               
if (!require("psych")) install.packages("psych")         # For descriptive statistics
if(!requireNamespace("devtools")) install.packages("devtools") # Contains financial data, including DJ_d and DJ_w
devtools::install_github("yukai-yang/FE")
if (!require("xtable")) install.packages("xtable")
library(xtable)

library(tidyverse)
library(FE)
library(psych)

############################
####Question A1.############# 
############################

#table 1, table 2. We extracted the values for the ME and MEBE portfolios from the summary, using the summary function

summary(portfolio_m)

############################
####Question A2.############# 
############################

##############A2(a)##############

## This creates the table in the appendix with all the values)

alpha <- ret$alpha
Sigma <- diag(crossprod(mZ - vZm %*% ret$beta) / nrow(mZ)) # Estimate residual variance
alpha_se <- sqrt(Sigma / nrow(mZ)) # Standard errors of alpha

t_stat <- alpha / alpha_se
p_values <- 2 * (1 - pnorm(abs(t_stat))) # Two-tailed p-values

cat("Alpha Estimates and Tests:\n")
for (i in 1:length(alpha)) {
  cat(sprintf("Portfolio %d: Alpha = %.4f, t-stat = %.4f, p-value = %.4f\n", 
              i, alpha[i], t_stat[i], p_values[i]))
}

# Load the xtable library
library(xtable)

# Create a data frame with the results
results_df <- data.frame(
  Portfolio = 1:length(alpha),
  Alpha = alpha,
  `t-stat` = t_stat,
  `p-value` = p_values
)

# Convert the data frame to a LaTeX table
latex_table <- xtable(results_df, caption = "Alpha Estimates and Tests", label = "tab:alpha_tests")

# Print the LaTeX table to the console (you can copy it from here)
print(latex_table, include.rownames = FALSE)


### This filters out the signicant values - table 3


# Rename columns to make them easier to work with
colnames(results_df) <- c("Portfolio", "Alpha", "t_stat", "p_value")

# Filter significant results (p_value <= 0.05)
significant_results <- results_df %>%
  filter(p_value <= 0.05)

# Print significant results to the console
cat("Significant Alpha Estimates:\n")
print(significant_results)

# Convert the significant results to a LaTeX table
significant_latex_table <- xtable(significant_results, 
                                  caption = "Significant Alpha Estimates and Tests (p-value <= 0.05)", 
                                  label = "tab:significant_alpha_tests")

# Print the LaTeX table for significant results
print(significant_latex_table, include.rownames = FALSE)


##############A2(b)##############
#table 4#

# Define the deciles
deciles <- 1:10

# Define the ME-beta values
ME_beta <- c(1.227, 1.269, 1.256, 1.211, 1.178, 1.155, 1.169, 1.100, 1.012, 0.900)

# Define the MEBE-beta values
MEBE_beta <- c(1.075, 1.046, 1.026, 1.025, 0.923, 0.912, 0.878, 0.910, 0.952, 1.082)

# Create the data frame
beta_df <- data.frame(
  Decile = deciles,
  ME_beta = ME_beta,
  MEBE_beta = MEBE_beta
)

# View the data frame
print(beta_df)


#figure 1#

# Scatter plot for ME-beta
plot_ME <- ggplot(beta_df, aes(x = Decile, y = ME_beta)) +
  geom_point(color = "blue", size = 3) +
  geom_line(color = "blue", linetype = "dashed") +
  labs(
    title = "ME-beta Across Deciles",
    x = "Decile",
    y = "ME-beta"
  ) +
  theme_minimal()

# Scatter plot for MEBE-beta
plot_MEBE <- ggplot(beta_df, aes(x = Decile, y = MEBE_beta)) +
  geom_point(color = "darkgreen", size = 3) +
  geom_line(color = "darkgreen", linetype = "dashed") +
  labs(
    title = "MEBE-beta Across Deciles",
    x = "Decile",
    y = "MEBE-beta"
  ) +
  theme_minimal()

# Arrange plots side by side using gridExtra
grid.arrange(plot_ME, plot_MEBE, ncol = 2)

# Alternatively, arrange plots using patchwork
# combined_plot <- plot_ME + plot_MEBE + plot_layout(ncol = 2)
# print(combined_plot)


#BETA FOR ME AND MEBE AND R2
FamaMacbethR2 <- function(mZ, vZm)
{
  iT = nrow(mZ)
  iN = ncol(mZ)
  
  # First Pass
  mX = cbind(1,vZm)
  pars = chol2inv(chol(crossprod(mX)))%*%crossprod(mX,mZ)
  alpha = c(pars[1,])
  # estimated betas
  beta = c(pars[2,])
  # sample mean of each excess return
  mu = apply(mZ,2,mean)
  # sample mean of the market excess return
  mum = mean(vZm)
  # variance of the market excess return
  sigm2 = c(crossprod(vZm-mum))/iT
  # residules
  me = mZ-mX%*%pars
  # covariance matrix of the residules
  Sigma = crossprod(me)/iT


##############A2(c)##############

  
  #Calculating R2
  SST = apply(mZ, 2, function(x) sum((x - mean(x))^2)) # Total Sum of Squares
  SSR = apply(me, 2, function(x) sum(x^2))             # Residual Sum of Squares
  R2 = 1 - SSR / SST    
  
  # do the test on pp.24 & 25 in slides for lecture 4
  tmp = c(t(alpha)%*%chol2inv(chol(Sigma))%*%alpha)/(1+mum^2/sigm2)
  
  # Gibbons/Ross/Shanken (1989)
  J1 = tmp*(iT-iN-1)/iN
  jpval = 1-pf(J1,df1=iN,df2=iT-iN-1)
  
  # Wald
  W = tmp*iT
  wpval = 1-pchisq(W,df=iN)
  
  # modified LR, Jobson/Korkie (1982)
  mX = matrix(vZm,iT,1)
  pars = chol2inv(chol(crossprod(mX)))%*%crossprod(mX,mZ)
  me = mZ-mX%*%pars
  # covariance matrix of the residules
  SigmaR = crossprod(me)/iT
  LR = (sum(log(eigen(SigmaR)$values))-sum(log(eigen(Sigma)$values)))*(iT-iN/2-2)
  lrpval = 1-pchisq(LR,df=iN)
  
  
  # Second Pass
  mX = cbind(1,beta)
  xxinv = chol2inv(chol(crossprod(mX)))
  xxinvx = tcrossprod(xxinv,mX)
  gamma = NULL; gamsd = NULL
  for(iter in 1:iT){
    vy = mZ[iter,]
    tmp = xxinvx%*%vy
    # gamma
    gamma = rbind(gamma, c(tmp))
    s2 = c(crossprod(vy - mX%*%tmp))/iN
    # gammas' standard error
    gamsd = rbind(gamsd, sqrt(diag(s2 * xxinv)))
  }
  
  # gamma0 and market risk premia
  gamma0 = c(gamma[,1]); g0sd = c(gamsd[,1])
  mrprem = c(gamma[,2]); prsd = c(gamsd[,2])
  
  # w tests on pp.30 in slides for lecture 4
  # should be compared with student t with T-1 degrees of freedom
  # or with standard normal if T is large
  wgamma0 = mean(gamma0)/sqrt(sum((gamma0-mean(gamma0))**2)/iT/(iT-1))
  wgamma1 = mean(mrprem)/sqrt(sum((mrprem-mean(mrprem))**2)/iT/(iT-1))
  
  return(list(
    alpha = alpha, beta = beta, R2 = R2, J1 = c(J1, jpval), W = c(W, wpval), 
    mLR = c(LR, lrpval), gamma0 = gamma0, g0sd = g0sd, mrprem = mrprem, prsd = prsd, 
    wgamma0 = wgamma0, wgamma1 = wgamma1
  ))
}

subsample = 501:700
mR = as.matrix(portfolio_m[subsample,5:24])
Rf = as.matrix(portfolio_m[subsample,'Tbill'])
Rm = as.matrix(portfolio_m[subsample,'Market'])
mZ = sweep(mR,1,Rf)
mR - c(Rf)

vZm = Rm - Rf

tmp = EstCAPM(mZ, vZm)
tmp$beta

Fama_MacBeth(mZ, vZm)
ret = FamaMacbethR2(mZ=mZ, vZm=vZm)
ret





##############A3##############

# Create a data frame with your alpha estimates
alpha_table <- data.frame(
  
  t_stat = c(-1.29,  0.49,  0.41,  1.65,  1.21,  1.14,  2.30,  2.17,  2.82,  2.53,
             -0.55,  0.75,  1.74,  0.29,  1.08,  1.94,  2.74,  1.52,  2.74,  2.77,
             -1.09,  0.75,  1.01,  0.42,  1.23,  2.42,  2.19,  2.91,  3.39,  2.41,
             -0.54,  0.15,  1.87,  3.83,  1.37,  2.81,  2.82,  2.68,  3.39,  2.41,
             -1.04,  0.57,  1.22,  0.83,  1.20,  2.42,  1.64,  1.91,  2.26,  1.07,
             -0.67, -0.22, -0.33, -0.83,  2.75,  1.80,  1.48,  2.24,  2.24,  1.34,
             -2.00,  0.06, -0.44, -0.03,  0.39,  0.52,  0.64,  1.14,  0.81,  1.27, 
             0.98,  1.05,  1.51,  1.73,  1.99,  1.18,  1.33,  1.45,  0.14, -1.06, 
             -1.79, -0.20,  0.48,  0.49,  0.26,  1.89,  1.53,  2.96,  2.79,  2.57)
)

# Verify the structure of the data frame
str(alpha_table)

alpha_table <- data.frame(
  t_stat = t_stat,
  stringsAsFactors = FALSE  # Prevents automatic conversion to factors
)


# Calculate the sum of squared t-statistics
test_statistic <- sum(alpha_table$t_stat^2)

# Determine degrees of freedom
degrees_freedom <- nrow(alpha_table)  # 100

# Calculate the p-value from the chi-squared distribution
p_value_joint <- 1 - pchisq(test_statistic, df = degrees_freedom)

# Display the results
cat("Joint F-test for all intercepts (alphas) being zero:\n")
cat("Test Statistic (Sum of t^2):", test_statistic, "\n")
cat("Degrees of Freedom:", degrees_freedom, "\n")
cat("p-value:", p_value_joint, "\n")


##############A4##############



FamaMacbeth0 <- function(mZ, vZm)
{
  iT = nrow(mZ)
  iN = ncol(mZ)
  
  # First Pass
  mX = cbind(1,vZm)
  pars = chol2inv(chol(crossprod(mX)))%*%crossprod(mX,mZ)
  alpha = c(pars[1,])
  # estimated betas
  beta = c(pars[2,])
  # sample mean of each excess return
  mu = apply(mZ,2,mean)
  # sample mean of the market excess return
  mum = mean(vZm)
  # variance of the market excess return
  sigm2 = c(crossprod(vZm-mum))/iT
  # residules
  me = mZ-mX%*%pars
  # covariance matrix of the residules
  Sigma = crossprod(me)/iT
  
  #Calculating R2
  SST = apply(mZ, 2, function(x) sum((x - mean(x))^2)) # Total Sum of Squares
  SSR = apply(me, 2, function(x) sum(x^2))             # Residual Sum of Squares
  R2 = 1 - SSR / SST    
  
  # do the test on pp.24 & 25 in slides for lecture 4
  tmp = c(t(alpha)%*%chol2inv(chol(Sigma))%*%alpha)/(1+mum^2/sigm2)
  
  # Gibbons/Ross/Shanken (1989)
  J1 = tmp*(iT-iN-1)/iN
  jpval = 1-pf(J1,df1=iN,df2=iT-iN-1)
  
  # Wald
  W = tmp*iT
  wpval = 1-pchisq(W,df=iN)
  
  # modified LR, Jobson/Korkie (1982)
  mX = matrix(vZm,iT,1)
  pars = chol2inv(chol(crossprod(mX)))%*%crossprod(mX,mZ)
  me = mZ-mX%*%pars
  # covariance matrix of the residules
  SigmaR = crossprod(me)/iT
  LR = (sum(log(eigen(SigmaR)$values))-sum(log(eigen(Sigma)$values)))*(iT-iN/2-2)
  lrpval = 1-pchisq(LR,df=iN)
  
  
  # Second Pass
  mX = cbind(1,beta)
  xxinv = chol2inv(chol(crossprod(mX)))
  xxinvx = tcrossprod(xxinv,mX)
  gamma = NULL; gamsd = NULL
  for(iter in 1:iT){
    vy = mZ[iter,]
    tmp = xxinvx%*%vy
    # gamma
    gamma = rbind(gamma, c(tmp))
    s2 = c(crossprod(vy - mX%*%tmp))/iN
    # gammas' standard error
    gamsd = rbind(gamsd, sqrt(diag(s2 * xxinv)))
  }
  
  # gamma0 and market risk premia
  gamma0 = c(gamma[,1]); g0sd = c(gamsd[,1])
  mrprem = c(gamma[,2]); prsd = c(gamsd[,2])
  
  # w tests on pp.30 in slides for lecture 4
  # should be compared with student t with T-1 degrees of freedom
  # or with standard normal if T is large
  wgamma0 = mean(gamma0)/sqrt(sum((gamma0-mean(gamma0))**2)/iT/(iT-1))
  wgamma1 = mean(mrprem)/sqrt(sum((mrprem-mean(mrprem))**2)/iT/(iT-1))
  
  return(list(alpha=alpha,beta=beta,J1=c(J1,jpval),W=c(W,wpval),mLR=c(LR,lrpval),
              gamma0=gamma0,g0sd=g0sd,mrprem=mrprem,prsd=prsd,wgamma0=wgamma0,wgamma1=wgamma1, SSR=SSR))
}

subsample = 1:942
mR = as.matrix(portfolio_m[subsample,25:124])
Rf = as.matrix(portfolio_m[subsample,'Tbill'])
Rm = as.matrix(portfolio_m[subsample,'Market'])
mZ = sweep(mR,1,Rf)
mR - c(Rf)

vZm = Rm - Rf

tmp = EstCAPM(mZ, vZm)
tmp$beta

FamaMacbeth0(mZ, vZm)
ret = FamaMacbeth0(mZ=mZ, vZm=vZm)
alphaFull <- ret$alpha
betaFull <- ret$beta
SSRFull <- ret$SSR
ret


# Determine the midpoint
mid_point <- floor(nrow(portfolio_m) / 2)

# Split the dataset into two halves
first_halfPM <- portfolio_m[1:mid_point, ]
second_halfPM <- portfolio_m[(mid_point + 1):nrow(portcap_m), ]

#FIRST SUBSAMPLE
FamaMacbethSS1 <- function(mZ, vZm)
{
  iT = nrow(mZ)
  iN = ncol(mZ)
  
  # First Pass
  mX = cbind(1,vZm)
  pars = chol2inv(chol(crossprod(mX)))%*%crossprod(mX,mZ)
  alpha = c(pars[1,])
  # estimated betas
  beta = c(pars[2,])
  # sample mean of each excess return
  mu = apply(mZ,2,mean)
  # sample mean of the market excess return
  mum = mean(vZm)
  # variance of the market excess return
  sigm2 = c(crossprod(vZm-mum))/iT
  # residules
  me = mZ-mX%*%pars
  # covariance matrix of the residules
  Sigma = crossprod(me)/iT
  
  #Calculating R2
  SST = apply(mZ, 2, function(x) sum((x - mean(x))^2)) # Total Sum of Squares
  SSR = apply(me, 2, function(x) sum(x^2))             # Residual Sum of Squares
  R2 = 1 - SSR / SST    
  
  # do the test on pp.24 & 25 in slides for lecture 4
  tmp = c(t(alpha)%*%chol2inv(chol(Sigma))%*%alpha)/(1+mum^2/sigm2)
  
  # Gibbons/Ross/Shanken (1989)
  J1 = tmp*(iT-iN-1)/iN
  jpval = 1-pf(J1,df1=iN,df2=iT-iN-1)
  
  # Wald
  W = tmp*iT
  wpval = 1-pchisq(W,df=iN)
  
  # modified LR, Jobson/Korkie (1982)
  mX = matrix(vZm,iT,1)
  pars = chol2inv(chol(crossprod(mX)))%*%crossprod(mX,mZ)
  me = mZ-mX%*%pars
  # covariance matrix of the residules
  SigmaR = crossprod(me)/iT
  LR = (sum(log(eigen(SigmaR)$values))-sum(log(eigen(Sigma)$values)))*(iT-iN/2-2)
  lrpval = 1-pchisq(LR,df=iN)
  
  
  # Second Pass
  mX = cbind(1,beta)
  xxinv = chol2inv(chol(crossprod(mX)))
  xxinvx = tcrossprod(xxinv,mX)
  gamma = NULL; gamsd = NULL
  for(iter in 1:iT){
    vy = mZ[iter,]
    tmp = xxinvx%*%vy
    # gamma
    gamma = rbind(gamma, c(tmp))
    s2 = c(crossprod(vy - mX%*%tmp))/iN
    # gammas' standard error
    gamsd = rbind(gamsd, sqrt(diag(s2 * xxinv)))
  }
  
  # gamma0 and market risk premia
  gamma0 = c(gamma[,1]); g0sd = c(gamsd[,1])
  mrprem = c(gamma[,2]); prsd = c(gamsd[,2])
  
  # w tests on pp.30 in slides for lecture 4
  # should be compared with student t with T-1 degrees of freedom
  # or with standard normal if T is large
  wgamma0 = mean(gamma0)/sqrt(sum((gamma0-mean(gamma0))**2)/iT/(iT-1))
  wgamma1 = mean(mrprem)/sqrt(sum((mrprem-mean(mrprem))**2)/iT/(iT-1))
  
  return(list(
    alpha = alpha, beta = beta, R2 = R2, J1 = c(J1, jpval), W = c(W, wpval), 
    mLR = c(LR, lrpval), gamma0 = gamma0, g0sd = g0sd, mrprem = mrprem, prsd = prsd, 
    wgamma0 = wgamma0, wgamma1 = wgamma1, SSR=SSR
  ))
}

subsample = 1:471
mR = as.matrix(first_halfPM[subsample,25:124])
Rf = as.matrix(first_halfPM[subsample,'Tbill'])
Rm = as.matrix(first_halfPM[subsample,'Market'])
mZ = sweep(mR,1,Rf)
mR - c(Rf)

vZm = Rm - Rf

tmp = EstCAPM(mZ, vZm)
tmp$beta

FamaMacbethSS1(mZ, vZm)
retSS1 = FamaMacbethSS1(mZ=mZ, vZm=vZm)
alphaSS1<-retSS1$alpha
betaSS1 <- retSS1$beta
SSRSS1 <- retSS1$SSR
retSS1


#SECOND SUBSAMPLE
FamaMacbethSS2 <- function(mZ, vZm)
{
  iT = nrow(mZ)
  iN = ncol(mZ)
  
  # First Pass
  mX = cbind(1,vZm)
  pars = chol2inv(chol(crossprod(mX)))%*%crossprod(mX,mZ)
  alpha = c(pars[1,])
  # estimated betas
  beta = c(pars[2,])
  # sample mean of each excess return
  mu = apply(mZ,2,mean)
  # sample mean of the market excess return
  mum = mean(vZm)
  # variance of the market excess return
  sigm2 = c(crossprod(vZm-mum))/iT
  # residules
  me = mZ-mX%*%pars
  # covariance matrix of the residules
  Sigma = crossprod(me)/iT
  
  #Calculating R2
  SST = apply(mZ, 2, function(x) sum((x - mean(x))^2)) # Total Sum of Squares
  SSR = apply(me, 2, function(x) sum(x^2))             # Residual Sum of Squares
  R2 = 1 - SSR / SST    
  
  # do the test on pp.24 & 25 in slides for lecture 4
  tmp = c(t(alpha)%*%chol2inv(chol(Sigma))%*%alpha)/(1+mum^2/sigm2)
  
  # Gibbons/Ross/Shanken (1989)
  J1 = tmp*(iT-iN-1)/iN
  jpval = 1-pf(J1,df1=iN,df2=iT-iN-1)
  
  # Wald
  W = tmp*iT
  wpval = 1-pchisq(W,df=iN)
  
  # modified LR, Jobson/Korkie (1982)
  mX = matrix(vZm,iT,1)
  pars = chol2inv(chol(crossprod(mX)))%*%crossprod(mX,mZ)
  me = mZ-mX%*%pars
  # covariance matrix of the residules
  SigmaR = crossprod(me)/iT
  LR = (sum(log(eigen(SigmaR)$values))-sum(log(eigen(Sigma)$values)))*(iT-iN/2-2)
  lrpval = 1-pchisq(LR,df=iN)
  
  
  # Second Pass
  mX = cbind(1,beta)
  xxinv = chol2inv(chol(crossprod(mX)))
  xxinvx = tcrossprod(xxinv,mX)
  gamma = NULL; gamsd = NULL
  for(iter in 1:iT){
    vy = mZ[iter,]
    tmp = xxinvx%*%vy
    # gamma
    gamma = rbind(gamma, c(tmp))
    s2 = c(crossprod(vy - mX%*%tmp))/iN
    # gammas' standard error
    gamsd = rbind(gamsd, sqrt(diag(s2 * xxinv)))
  }
  
  # gamma0 and market risk premia
  gamma0 = c(gamma[,1]); g0sd = c(gamsd[,1])
  mrprem = c(gamma[,2]); prsd = c(gamsd[,2])
  
  # w tests on pp.30 in slides for lecture 4
  # should be compared with student t with T-1 degrees of freedom
  # or with standard normal if T is large
  wgamma0 = mean(gamma0)/sqrt(sum((gamma0-mean(gamma0))**2)/iT/(iT-1))
  wgamma1 = mean(mrprem)/sqrt(sum((mrprem-mean(mrprem))**2)/iT/(iT-1))
  
  return(list(
    alpha = alpha, beta = beta, R2 = R2, J1 = c(J1, jpval), W = c(W, wpval), 
    mLR = c(LR, lrpval), gamma0 = gamma0, g0sd = g0sd, mrprem = mrprem, prsd = prsd, 
    wgamma0 = wgamma0, wgamma1 = wgamma1, SSR=SSR
  ))
}

subsample = 1:471
mR = as.matrix(second_halfPM[subsample,25:124])
Rf = as.matrix(second_halfPM[subsample,'Tbill'])
Rm = as.matrix(second_halfPM[subsample,'Market'])
mZ = sweep(mR,1,Rf)
mR - c(Rf)

vZm = Rm - Rf

tmp = EstCAPM(mZ, vZm)
tmp$beta

FamaMacbethSS2(mZ, vZm)
retSS2 = FamaMacbethSS2(mZ=mZ, vZm=vZm)
alphaSS2<-retSS2$alpha
betaSS2 <- retSS2$beta
SSRSS2 <- retSS2$SSR
retSS2



###CHOW/F-TEST AND P-VALUES###


#Prepare to store results
f_test_result <- list()
p_value_result <- list()

# Loop through each Value and run test

  
  # Run the regression
  for (i in 1:length(alphaSS1)) {
  f_stat <- ((SSRFull[i] - (SSRSS1[i]+SSRSS2[i]))/2)/
    ((SSRSS1[i] + SSRSS2[i]) / (471 - 2 * 2))
  f_test_result[[paste0("f_stat_", i)]] <- f_stat
  
  # Store the result in the list
p_value <- pf(f_stat,df1 = 2,df2 = 471-2*2,lower.tail = FALSE)

p_value_result[[paste0("p_value_", i)]] <- p_value
}
# Inspect results
for (name in names(f_test_result)) {
  cat("\nSummary for", name, ":\n")
  print(summary(f_test_result[[name]]))
}
for (name in names(p_value_result)) {
  cat("\nSummary for", name, ":\n")
  print(summary(p_value_result[[name]]))
}

#Plotting p-value frequencie
ggplot(df, aes(x = P_value)) +
  geom_histogram(binwidth = 0.05, fill = "lightgreen", color = "black", alpha = 0.7) +
  geom_vline(xintercept = 0.075, linetype = "solid", color = "red", size=1.1) +  # Add vertical line at p = 0.05
  labs(title = "Histogram of P-values with Significance Threshold",
       x = "P-value",
       y = "Frequency") +
  theme_minimal() +
  theme(text = element_text(size = 12))

  #Plotting first diff for alpha
# Calculate the first differences
first_differenceAlphas <- alphaSS1 - alphaSS2

# Plot the first differences
plot(first_differenceAlphas, type = "o", main = "First Differences Between Alpha period 1 & 2",
     xlab = "Index", ylab = "First Difference", col = "blue", lwd = 1.5)

# Add horizontal line at 0 for reference
abline(h = 0, col = "red", lty = 5, lwd=2)



#For beta
# Calculate the first differences
first_differenceBetas <- betaSS1 - betaSS2

# Plot the first differences
plot(first_differenceBetas, type = "o", main = "First Differences Between Beta period 1 & 2",
     xlab = "Index", ylab = "First Difference", col = "blue", lwd = 1.5)

# Add horizontal line at 0 for reference
abline(h = 0, col = "red", lty = 5, lwd=2)


##############A4(a)##############

#Time Subsamples Alpha and Beta

FamaMacbeth0 <- function(mZ, vZm)
{
  iT = nrow(mZ)
  iN = ncol(mZ)
  
  # First Pass
  mX = cbind(1,vZm)
  pars = chol2inv(chol(crossprod(mX)))%*%crossprod(mX,mZ)
  alpha = c(pars[1,])
  # estimated betas
  beta = c(pars[2,])
  # sample mean of each excess return
  mu = apply(mZ,2,mean)
  # sample mean of the market excess return
  mum = mean(vZm)
  # variance of the market excess return
  sigm2 = c(crossprod(vZm-mum))/iT
  # residules
  me = mZ-mX%*%pars
  # covariance matrix of the residules
  Sigma = crossprod(me)/iT
  
  #Calculating R2
  SST = apply(mZ, 2, function(x) sum((x - mean(x))^2)) # Total Sum of Squares
  SSR = apply(me, 2, function(x) sum(x^2))             # Residual Sum of Squares
  R2 = 1 - SSR / SST    
  
  # do the test on pp.24 & 25 in slides for lecture 4
  tmp = c(t(alpha)%*%chol2inv(chol(Sigma))%*%alpha)/(1+mum^2/sigm2)
  
  # Gibbons/Ross/Shanken (1989)
  J1 = tmp*(iT-iN-1)/iN
  jpval = 1-pf(J1,df1=iN,df2=iT-iN-1)
  
  # Wald
  W = tmp*iT
  wpval = 1-pchisq(W,df=iN)
  
  # modified LR, Jobson/Korkie (1982)
  mX = matrix(vZm,iT,1)
  pars = chol2inv(chol(crossprod(mX)))%*%crossprod(mX,mZ)
  me = mZ-mX%*%pars
  # covariance matrix of the residules
  SigmaR = crossprod(me)/iT
  LR = (sum(log(eigen(SigmaR)$values))-sum(log(eigen(Sigma)$values)))*(iT-iN/2-2)
  lrpval = 1-pchisq(LR,df=iN)
  
  
  # Second Pass
  mX = cbind(1,beta)
  xxinv = chol2inv(chol(crossprod(mX)))
  xxinvx = tcrossprod(xxinv,mX)
  gamma = NULL; gamsd = NULL
  for(iter in 1:iT){
    vy = mZ[iter,]
    tmp = xxinvx%*%vy
    # gamma
    gamma = rbind(gamma, c(tmp))
    s2 = c(crossprod(vy - mX%*%tmp))/iN
    # gammas' standard error
    gamsd = rbind(gamsd, sqrt(diag(s2 * xxinv)))
  }
  
  # gamma0 and market risk premia
  gamma0 = c(gamma[,1]); g0sd = c(gamsd[,1])
  mrprem = c(gamma[,2]); prsd = c(gamsd[,2])
  
  # w tests on pp.30 in slides for lecture 4
  # should be compared with student t with T-1 degrees of freedom
  # or with standard normal if T is large
  wgamma0 = mean(gamma0)/sqrt(sum((gamma0-mean(gamma0))**2)/iT/(iT-1))
  wgamma1 = mean(mrprem)/sqrt(sum((mrprem-mean(mrprem))**2)/iT/(iT-1))
  
  return(list(alpha=alpha,beta=beta,J1=c(J1,jpval),W=c(W,wpval),mLR=c(LR,lrpval),
              gamma0=gamma0,g0sd=g0sd,mrprem=mrprem,prsd=prsd,wgamma0=wgamma0,wgamma1=wgamma1, SSR=SSR))
}

subsample = 1:942
mR = as.matrix(portfolio_m[subsample,25:124])
Rf = as.matrix(portfolio_m[subsample,'Tbill'])
Rm = as.matrix(portfolio_m[subsample,'Market'])
mZ = sweep(mR,1,Rf)
mR - c(Rf)

vZm = Rm - Rf

tmp = EstCAPM(mZ, vZm)
tmp$beta

FamaMacbeth0(mZ, vZm)
ret = FamaMacbeth0(mZ=mZ, vZm=vZm)
alphaFull <- ret$alpha
betaFull <- ret$beta
SSRFull <- ret$SSR
ret




##############A4(b)##############
 
FamaMacbeth0 <- function(mZ, vZm)
{
  iT = nrow(mZ)
  iN = ncol(mZ)
  
  
  mX = cbind(1, vZm)
  pars = chol2inv(chol(crossprod(mX))) %*% crossprod(mX, mZ)
  
  # First Pass
  mX = cbind(1,vZm)
  pars = chol2inv(chol(crossprod(mX)))%*%crossprod(mX,mZ)
  alpha = c(pars[1,])
  # estimated betas
  beta = c(pars[2,])
  # sample mean of each excess return
  mu = apply(mZ,2,mean)
  # sample mean of the market excess return
  mum = mean(vZm)
  # variance of the market excess return
  sigm2 = c(crossprod(vZm-mum))/iT
  # residules
  me = mZ-mX%*%pars
  # covariance matrix of the residules
  Sigma = crossprod(me)/iT
  
  #Calculating R2
  SST = apply(mZ, 2, function(x) sum((x - mean(x))^2)) # Total Sum of Squares
  SSR = apply(me, 2, function(x) sum(x^2))             # Residual Sum of Squares
  R2 = 1 - SSR / SST    
  
  # do the test on pp.24 & 25 in slides for lecture 4
  tmp = c(t(alpha)%*%chol2inv(chol(Sigma))%*%alpha)/(1+mum^2/sigm2)
  
  # Gibbons/Ross/Shanken (1989)
  J1 = tmp*(iT-iN-1)/iN
  jpval = 1-pf(J1,df1=iN,df2=iT-iN-1)
  
  # Wald
  W = tmp*iT
  wpval = 1-pchisq(W,df=iN)
  
  # modified LR, Jobson/Korkie (1982)
  mX = matrix(vZm,iT,1)
  pars = chol2inv(chol(crossprod(mX)))%*%crossprod(mX,mZ)
  me = mZ-mX%*%pars
  # covariance matrix of the residules
  SigmaR = crossprod(me)/iT
  LR = (sum(log(eigen(SigmaR)$values))-sum(log(eigen(Sigma)$values)))*(iT-iN/2-2)
  lrpval = 1-pchisq(LR,df=iN)
  
  
  # Second Pass
  mX = cbind(1,beta)
  xxinv = chol2inv(chol(crossprod(mX)))
  xxinvx = tcrossprod(xxinv,mX)
  gamma = NULL; gamsd = NULL
  for(iter in 1:iT){
    vy = mZ[iter,]
    tmp = xxinvx%*%vy
    # gamma
    gamma = rbind(gamma, c(tmp))
    s2 = c(crossprod(vy - mX%*%tmp))/iN
    # gammas' standard error
    gamsd = rbind(gamsd, sqrt(diag(s2 * xxinv)))
  }
  
  # gamma0 and market risk premia
  gamma0 = c(gamma[,1]); g0sd = c(gamsd[,1])
  mrprem = c(gamma[,2]); prsd = c(gamsd[,2])
  
  # w tests on pp.30 in slides for lecture 4
  # should be compared with student t with T-1 degrees of freedom
  # or with standard normal if T is large
  wgamma0 = mean(gamma0)/sqrt(sum((gamma0-mean(gamma0))**2)/iT/(iT-1))
  wgamma1 = mean(mrprem)/sqrt(sum((mrprem-mean(mrprem))**2)/iT/(iT-1))
  
  return(list(alpha=alpha,beta=beta,J1=c(J1,jpval),W=c(W,wpval),mLR=c(LR,lrpval),
              gamma0=gamma0,g0sd=g0sd,mrprem=mrprem,prsd=prsd,wgamma0=wgamma0,wgamma1=wgamma1, SSR=SSR, R2=R2))
}

subsample = 1:942
mR = as.matrix(portfolio_m[subsample,25:124])
Rf = as.matrix(portfolio_m[subsample,'Tbill'])
Rm = as.matrix(portfolio_m[subsample,'Market'])
mZ = sweep(mR,1,Rf)
mR - c(Rf)

vZm = Rm - Rf

tmp = EstCAPM(mZ, vZm)
tmp$beta

FamaMacbeth0(mZ, vZm)
ret = FamaMacbeth0(mZ=mZ, vZm=vZm)
alphaFull <- ret$alpha
betaFull <- ret$beta
SSRFull <- ret$SSR
R2Full <- ret$R2
ret

######## TEST WITH TIME SPECIFIC EFFECTS
# Convert input matrices
mR = as.matrix(portfolio_m[subsample, 25:124])
Rf = as.matrix(portfolio_m[subsample, 'Tbill'])
Rm = as.matrix(portfolio_m[subsample, 'Market'])
mZ = sweep(mR, 1, Rf)
vZm = Rm - Rf

# Create period-specific effects (dummy variables)
year <- as.factor(portfolio_m[subsample, 'year', drop = TRUE])
# Create dummy variables for periods
period_dummies <- model.matrix(~ year - 1)  # Ensure 'year' is a factor
# Include period dummies as additional explanatory variables
vZm_with_period <- cbind(vZm, period_dummies)

# Estimate the CAPM model with period-specific effects
tmp <- EstCAPM(mZ, vZm_with_period)
tmp$beta

# Fama-Macbeth regression with period-specific effects
ret <- FamaMacbeth0(mZ = mZ, vZm = vZm_with_period)
alphaFullPS <- ret$alpha
betaFullPS <- ret$beta
SSRFullPS <- ret$SSR
R2FullPS <- ret$R2

avg_R2_no_period <- mean(R2Full)
avg_R2_with_period <- mean(R2FullPS)

# Output the results
ret
avg_R2_no_period
avg_R2_with_period

#Analyzing betas for ME and MEBE in two periods
mid_point <- floor(nrow(portfolio_m) / 2)

# Split the dataset into two halves
first_halfPM <- portfolio_m[1:mid_point, ]
second_halfPM <- portfolio_m[(mid_point + 1):nrow(portcap_m), ]

#FIRST SUBSAMPLE
FamaMacbethSS1 <- function(mZ, vZm)
{
  iT = nrow(mZ)
  iN = ncol(mZ)
  
  # First Pass
  mX = cbind(1,vZm)
  pars = chol2inv(chol(crossprod(mX)))%*%crossprod(mX,mZ)
  alpha = c(pars[1,])
  # estimated betas
  beta = c(pars[2,])
  # sample mean of each excess return
  mu = apply(mZ,2,mean)
  # sample mean of the market excess return
  mum = mean(vZm)
  # variance of the market excess return
  sigm2 = c(crossprod(vZm-mum))/iT
  # residules
  me = mZ-mX%*%pars
  # covariance matrix of the residules
  Sigma = crossprod(me)/iT
  
  #Calculating R2
  SST = apply(mZ, 2, function(x) sum((x - mean(x))^2)) # Total Sum of Squares
  SSR = apply(me, 2, function(x) sum(x^2))             # Residual Sum of Squares
  R2 = 1 - SSR / SST    
  
  # do the test on pp.24 & 25 in slides for lecture 4
  tmp = c(t(alpha)%*%chol2inv(chol(Sigma))%*%alpha)/(1+mum^2/sigm2)
  
  # Gibbons/Ross/Shanken (1989)
  J1 = tmp*(iT-iN-1)/iN
  jpval = 1-pf(J1,df1=iN,df2=iT-iN-1)
  
  # Wald
  W = tmp*iT
  wpval = 1-pchisq(W,df=iN)
  
  # modified LR, Jobson/Korkie (1982)
  mX = matrix(vZm,iT,1)
  pars = chol2inv(chol(crossprod(mX)))%*%crossprod(mX,mZ)
  me = mZ-mX%*%pars
  # covariance matrix of the residules
  SigmaR = crossprod(me)/iT
  LR = (sum(log(eigen(SigmaR)$values))-sum(log(eigen(Sigma)$values)))*(iT-iN/2-2)
  lrpval = 1-pchisq(LR,df=iN)
  
  
  # Second Pass
  mX = cbind(1,beta)
  xxinv = chol2inv(chol(crossprod(mX)))
  xxinvx = tcrossprod(xxinv,mX)
  gamma = NULL; gamsd = NULL
  for(iter in 1:iT){
    vy = mZ[iter,]
    tmp = xxinvx%*%vy
    # gamma
    gamma = rbind(gamma, c(tmp))
    s2 = c(crossprod(vy - mX%*%tmp))/iN
    # gammas' standard error
    gamsd = rbind(gamsd, sqrt(diag(s2 * xxinv)))
  }
  
  # gamma0 and market risk premia
  gamma0 = c(gamma[,1]); g0sd = c(gamsd[,1])
  mrprem = c(gamma[,2]); prsd = c(gamsd[,2])
  
  # w tests on pp.30 in slides for lecture 4
  # should be compared with student t with T-1 degrees of freedom
  # or with standard normal if T is large
  wgamma0 = mean(gamma0)/sqrt(sum((gamma0-mean(gamma0))**2)/iT/(iT-1))
  wgamma1 = mean(mrprem)/sqrt(sum((mrprem-mean(mrprem))**2)/iT/(iT-1))
  
  return(list(
    alpha = alpha, beta = beta, R2 = R2, J1 = c(J1, jpval), W = c(W, wpval), 
    mLR = c(LR, lrpval), gamma0 = gamma0, g0sd = g0sd, mrprem = mrprem, prsd = prsd, 
    wgamma0 = wgamma0, wgamma1 = wgamma1, SSR=SSR
  ))
}

subsample = 1:471
mR = as.matrix(first_halfPM[subsample,5:24])
Rf = as.matrix(first_halfPM[subsample,'Tbill'])
Rm = as.matrix(first_halfPM[subsample,'Market'])
mZ = sweep(mR,1,Rf)
mR - c(Rf)

vZm = Rm - Rf

tmp = EstCAPM(mZ, vZm)
tmp$beta

FamaMacbethSS1(mZ, vZm)
retSS1 = FamaMacbethSS1(mZ=mZ, vZm=vZm)
alphaSS1<-retSS1$alpha
betaSS1 <- retSS1$beta
SSRSS1 <- retSS1$SSR
retSS1


#SECOND SUBSAMPLE
FamaMacbethSS2 <- function(mZ, vZm)
{
  iT = nrow(mZ)
  iN = ncol(mZ)
  
  # First Pass
  mX = cbind(1,vZm)
  pars = chol2inv(chol(crossprod(mX)))%*%crossprod(mX,mZ)
  alpha = c(pars[1,])
  # estimated betas
  beta = c(pars[2,])
  # sample mean of each excess return
  mu = apply(mZ,2,mean)
  # sample mean of the market excess return
  mum = mean(vZm)
  # variance of the market excess return
  sigm2 = c(crossprod(vZm-mum))/iT
  # residules
  me = mZ-mX%*%pars
  # covariance matrix of the residules
  Sigma = crossprod(me)/iT
  
  #Calculating R2
  SST = apply(mZ, 2, function(x) sum((x - mean(x))^2)) # Total Sum of Squares
  SSR = apply(me, 2, function(x) sum(x^2))             # Residual Sum of Squares
  R2 = 1 - SSR / SST    
  
  # do the test on pp.24 & 25 in slides for lecture 4
  tmp = c(t(alpha)%*%chol2inv(chol(Sigma))%*%alpha)/(1+mum^2/sigm2)
  
  # Gibbons/Ross/Shanken (1989)
  J1 = tmp*(iT-iN-1)/iN
  jpval = 1-pf(J1,df1=iN,df2=iT-iN-1)
  
  # Wald
  W = tmp*iT
  wpval = 1-pchisq(W,df=iN)
  
  # modified LR, Jobson/Korkie (1982)
  mX = matrix(vZm,iT,1)
  pars = chol2inv(chol(crossprod(mX)))%*%crossprod(mX,mZ)
  me = mZ-mX%*%pars
  # covariance matrix of the residules
  SigmaR = crossprod(me)/iT
  LR = (sum(log(eigen(SigmaR)$values))-sum(log(eigen(Sigma)$values)))*(iT-iN/2-2)
  lrpval = 1-pchisq(LR,df=iN)
  
  
  # Second Pass
  mX = cbind(1,beta)
  xxinv = chol2inv(chol(crossprod(mX)))
  xxinvx = tcrossprod(xxinv,mX)
  gamma = NULL; gamsd = NULL
  for(iter in 1:iT){
    vy = mZ[iter,]
    tmp = xxinvx%*%vy
    # gamma
    gamma = rbind(gamma, c(tmp))
    s2 = c(crossprod(vy - mX%*%tmp))/iN
    # gammas' standard error
    gamsd = rbind(gamsd, sqrt(diag(s2 * xxinv)))
  }
  
  # gamma0 and market risk premia
  gamma0 = c(gamma[,1]); g0sd = c(gamsd[,1])
  mrprem = c(gamma[,2]); prsd = c(gamsd[,2])
  
  # w tests on pp.30 in slides for lecture 4
  # should be compared with student t with T-1 degrees of freedom
  # or with standard normal if T is large
  wgamma0 = mean(gamma0)/sqrt(sum((gamma0-mean(gamma0))**2)/iT/(iT-1))
  wgamma1 = mean(mrprem)/sqrt(sum((mrprem-mean(mrprem))**2)/iT/(iT-1))
  
  return(list(
    alpha = alpha, beta = beta, R2 = R2, J1 = c(J1, jpval), W = c(W, wpval), 
    mLR = c(LR, lrpval), gamma0 = gamma0, g0sd = g0sd, mrprem = mrprem, prsd = prsd, 
    wgamma0 = wgamma0, wgamma1 = wgamma1, SSR=SSR
  ))
}

subsample = 1:471
mR = as.matrix(second_halfPM[subsample,5:24])
Rf = as.matrix(second_halfPM[subsample,'Tbill'])
Rm = as.matrix(second_halfPM[subsample,'Market'])
mZ = sweep(mR,1,Rf)
mR - c(Rf)

vZm = Rm - Rf

tmp = EstCAPM(mZ, vZm)
tmp$beta

FamaMacbethSS2(mZ, vZm)
retSS2 = FamaMacbethSS2(mZ=mZ, vZm=vZm)
alphaSS2<-retSS2$alpha
betaSS2 <- retSS2$beta
SSRSS2 <- retSS2$SSR
retSS2

#Plotting the difference

# Extract the ME and MEBE
MEFirstsPerB <- betaSS1[1:10]
MEBEFirstsPerB <- betaSS1[11:20]

MESecondPerB <- betaSS2[1:10]
MEBESecondPerB <- betaSS2[11:20]

#FOR ME1-10
# Create a line plot
plot(MEFirstsPerB, type = "o", main = "Betas for ME 1-10", 
     xlab = "Index", ylab = "Betas", col = "blue", xaxt = "n",size=1, ylim = range(c(MEFirstsPerB, MESecondPerB)))

# Add custom x-axis labels
axis(1, at = 1:10, labels = names(MEFirstsPerB))

# Add the second line
lines(MESecondPerB, type = "o", col = "red")

# Add a legend to distinguish the two lines
legend("topright", legend = c("First Period", "Second Period"), col = c("blue", "red"), lty = 1, pch = 1)



#FOR MEBE1-10
plot(MEBEFirstsPerB, type = "o", main = "Betas for MEBE 1-10", 
     xlab = "Index", ylab = "Betas", col = "blue", xaxt = "n",size=1, ylim = range(c(MEBEFirstsPerB, MEBESecondPerB)))

# Add custom x-axis labels
axis(1, at = 1:10, labels = names(MEBEFirstsPerB))

# Add the second line
lines(MESecondPerB, type = "o", col = "red")

# Add a legend to distinguish the two lines
legend("topright", legend = c("First Period", "Second Period"), col = c("blue", "red"), lty = 1, pch = 1)



##############B1(a)##############

#cross sectional beta 
betas <- ret$beta
hist(betas)
library(psych)
describe(betas)
#test if betas are differnt than one 
t_test_result <- t.test(betas, mu = 1)
print(t_test_result)


#relationship between excess return and market betas 
# Calculate the average of each column in mZ
mean_mZ <- colMeans(mZ, na.rm = TRUE)

# View the result
print(mean_mZ)
length(mean_mZ)
mean(mean_mZ)

##############B1(c)##############


# Scatter plot of average excess returns vs market betas
plot(betas, mean_mZ, 
     xlab = "Market Betas", 
     ylab = "Average Excess Returns", 
     main = "Market Betas and Average Excess Returns",
     pch = 16, col = "black")

# Add a regression line to the plot
reg_model <- lm(mean_mZ ~ betas)  # Linear regression model
abline(reg_model, col = "red", lwd = 2)  # Add regression line

# Display regression summary
summary(reg_model)




##############B2(a)##############

a.\text{Test for } \gamma_{0t} = 0 \text{ for all } t

####gamma0####
# Loop over each month (row of mZ)
gamma0_results <- data.frame(Time = numeric(), Gamma0 = numeric(), SE = numeric(), Tstat = numeric(), Pvalue = numeric())

class(mZ)
dim(mZ)
View(mZ)
for (t in 1:nrow(mZ)) {
# Get Z_t (excess returns for all portfolios in month t)
  
  Z_t <- mZ[t, ]
  
  # Regress Z_t on constant and betas
  model <- lm(Z_t ~ ret$beta)
  
  # Extract regression results
  gamma0 <- coef(model)[1]  # Intercept (gamma0)
  gamma0_se <- summary(model)$coefficients[1, 2]  # Standard error of gamma0
  t_stat <- summary(model)$coefficients[1, 3]  # t-statistic
  p_value <- summary(model)$coefficients[1, 4]  # p-value
  
  # Store the results
  gamma0_results <- rbind(gamma0_results, data.frame(
    Time = t, 
    Gamma0 = gamma0, 
    SE = gamma0_se, 
    Tstat = t_stat, 
    Pvalue = p_value
  ))
}

# View results
print(gamma0_results)

# Identify months where gamma0 is significantly different from 0
significant_months <- subset(gamma0_results, Pvalue < 0.05)
cat("Number of months with significant gamma0:", nrow(significant_months), "\n")
cat("Proportion of months with significant gamma0:", nrow(significant_months) / nrow(mZ), "\n")


# Histogram of gamma0 p-values
ggplot(gamma0_results, aes(x = Pvalue)) +
  geom_histogram(binwidth = 0.05, fill = "blue", color = "black", alpha = 0.7) +
  labs(title = "Histogram of Gamma0 P-values",
       x = "P-value",
       y = "Frequency") +
  theme_minimal() +
  geom_vline(xintercept = 0.05, color = "red", linetype = "dashed") +
  annotate("text", x = 0.05, y = max(table(cut(gamma0_results$Pvalue, breaks=seq(0,1,0.05)))), 
           label = "Significance Threshold (0.05)", vjust = -0.5, hjust = -0.1, color = "red")



##############B2(b)##############



#gamma1 mrprem

# Run the FamaMacbeth0 function to get second-pass regression results
ret <- FamaMacbeth0(mZ = mZ, vZm = vZm)

# Extract the monthly estimates of gamma1 (mrprem) and their standard errors
mrprem <- ret$mrprem    # Vector of gamma1 (market risk premia) for each time period
dim(mrprem)

mrprem_se <- ret$prsd   # Vector of standard errors for gamma1

# Calculate t-statistics for gamma1
t_stat <- mrprem / mrprem_se

# Calculate p-values for one-sided test (H0: gamma1 <= 0, H1: gamma1 > 0)
p_values <- 1 - pnorm(t_stat)

# Combine results into a data frame
gamma1_results <- data.frame(
  Time = seq_along(mrprem),  # Time periods (e.g., months)
  Gamma1 = mrprem,
  SE = mrprem_se,
  Tstat = t_stat,
  Pvalue = p_values
)

# View the results
print(gamma1_results)

# Identify months where gamma1 > 0 and statistically significant
significant_gamma1 <- subset(gamma1_results, Gamma1 > 0 & Pvalue < 0.05)


# Histogram of gamma1 p-values
ggplot(gamma1_results, aes(x = Pvalue)) +
  geom_histogram(binwidth = 0.05, fill = "blue", color = "black", alpha = 0.7) +
  labs(title = "Histogram of Gamma1 P-values",
       x = "P-value",
       y = "Frequency") +
  theme_minimal() +
  geom_vline(xintercept = 0.05, color = "red", linetype = "dashed") +
  annotate("text", x = 0.05, y = max(table(cut(gamma0_results$Pvalue, breaks=seq(0,1,0.05)))), 
           label = "Significance Threshold (0.05)", vjust = -0.5, hjust = -0.1, color = "red")






# Summary of results
cat("Number of months with gamma1 > 0 and significant:", nrow(significant_gamma1), "\n")
cat("Proportion of months with gamma1 > 0 and significant:", nrow(significant_gamma1) / nrow(gamma1_results), "\n")


##############B2(c)##############


c.\text{Test for } \bar{\gamma}_0 = 0 \text{ and } \bar{\gamma}_1 > 0
##B2c 

# Calculate averages of gamma1 (mrprem) and gamma0
gamma1_mean <- mean(ret$mrprem)  # Average of gamma1
gamma1_se <- sqrt(var(ret$mrprem) / length(ret$mrprem))  # Standard error of gamma1
gamma0_mean <- mean(ret$gamma0)  # Average of gamma0
gamma0_se <- sqrt(var(ret$gamma0) / length(ret$gamma0))  # Standard error of gamma0

# Test for gamma1 > 0
t_stat_gamma1 <- gamma1_mean / gamma1_se
p_value_gamma1 <- 1 - pnorm(t_stat_gamma1)  # One-sided test

# Test for gamma0 = 0
t_stat_gamma0 <- gamma0_mean / gamma0_se
p_value_gamma0 <- 2 * (1 - pnorm(abs(t_stat_gamma0)))  # Two-sided test

# Print the results
cat("Results for gamma1:\n")
cat("Mean of gamma1 =", gamma1_mean, "\n")
cat("Standard Error =", gamma1_se, "\n")
cat("t-statistic =", t_stat_gamma1, "\n")
cat("p-value =", p_value_gamma1, "\n\n")

cat("Results for gamma0:\n")
cat("Mean of gamma0 =", gamma0_mean, "\n")
cat("Standard Error =", gamma0_se, "\n")
cat("t-statistic =", t_stat_gamma0, "\n")
cat("p-value =", p_value_gamma0, "\n")


# Create a data frame with results
results_table_gammas <- data.frame(
  Statistic = c("Mean of Gamma0", "Standard Error (Gamma0)", "t-statistic (Gamma0)", "p-value (Gamma0)",
                "Mean of Gamma1", "Standard Error (Gamma1)", "t-statistic (Gamma1)", "p-value (Gamma1)"),
  Value = c(gamma0_mean, gamma0_se, t_stat_gamma0, p_value_gamma0,
            gamma1_mean, gamma1_se, t_stat_gamma1, p_value_gamma1)
)

# Create the LaTeX table using xtable
library(xtable)
latex_table <- xtable(results_table_gammas, caption = "Results for Gamma0 and Gamma1 Testing", label = "tab:gamma_tests")
print(latex_table, include.rownames = FALSE)



#R2
# Calculate R-squared for second-pass regression
r_squared <- 1 - (sum((mZ - (ret$gamma0 + ret$mrprem * ret$beta))^2) / sum((mZ - mean(mZ))^2))
cat("Goodness-of-Fit (R-squared):", r_squared, "\n")

length(ret$beta)
length(ret$mrprem)

#fixing code 
# Ensure compatibility between gamma0, mrprem, and beta
predicted_mZ <- matrix(NA, nrow = nrow(mZ), ncol = ncol(mZ))  # Initialize predicted values

for (t in 1:nrow(mZ)) {
  predicted_mZ[t, ] <- ret$gamma0[t] + ret$mrprem[t] * ret$beta  # Predicted returns for each portfolio in period t
}

# Calculate R-squared
ss_total <- sum((mZ - mean(mZ))^2)  # Total sum of squares
ss_residual <- sum((mZ - predicted_mZ)^2)  # Residual sum of squares
r_squared <- 1 - (ss_residual / ss_total)

cat("R-squared for second-pass regression:", r_squared, "\n")


#To preform the robustness over time, I ran the same code
#but changed the subsample (1:471 or 472:942)



##############B3(a)##############


# Extract column names
portfolio_names <- colnames(portfolio_m)[25:124]
print(portfolio_names)

# Extract size and book-to-market deciles
size_deciles <- as.numeric(gsub("R(\\d{1,2})(\\d{1,2})$", "\\1", portfolio_names))
bm_deciles <- as.numeric(gsub("R(\\d{1,2})(\\d{1,2})$", "\\2", portfolio_names))
print(size_deciles)
size_deciles[size_deciles == 91] <- 9
print(bm_deciles)
bm_deciles[bm_deciles == 0] <- 10

###plot size and Bm to avg excess return 

# Scatter plot of average excess returns vs size decile
plot(size_deciles, mean_mZ, 
     xlab = "Size decile", 
     ylab = "Average Excess Returns", 
     main = "Size Decile and Average Excess Returns",
     pch = 16, col = "black")

# Add a regression line to the plot
reg_model <- lm(mean_mZ ~ size_deciles)  # Linear regression model
abline(reg_model, col = "red", lwd = 2)  # Add regression line

# Display regression summary
summary(reg_model)

# Scatter plot of average excess returns vs BM decile
plot(bm_deciles, mean_mZ, 
     xlab = "BM decile", 
     ylab = "Average Excess Returns", 
     main = "BM decile and Average Excess Returns",
     pch = 16, col = "black")

# Add a regression line to the plot
reg_model <- lm(mean_mZ ~ bm_deciles)  # Linear regression model
abline(reg_model, col = "red", lwd = 2)  # Add regression line

# Display regression summary
summary(reg_model)


##############B3(b)##############
#adding controls and then taking log of them
FamaMacbethWithControls <- function(mZ, vZm, size_control, bm_control) {
  iT <- nrow(mZ)
  iN <- ncol(mZ)
  
  # First Pass
  mX <- cbind(1, vZm)
  pars <- chol2inv(chol(crossprod(mX))) %*% crossprod(mX, mZ)
  alpha <- c(pars[1, ])
  beta <- c(pars[2, ])
  mu <- apply(mZ, 2, mean)
  mum <- mean(vZm)
  sigm2 <- c(crossprod(vZm - mum)) / iT
  me <- mZ - mX %*% pars
  Sigma <- crossprod(me) / iT
  
  # Gibbons/Ross/Shanken (1989)
  tmp <- c(t(alpha) %*% chol2inv(chol(Sigma)) %*% alpha) / (1 + mum^2 / sigm2)
  J1 <- tmp * (iT - iN - 1) / iN
  jpval <- 1 - pf(J1, df1 = iN, df2 = iT - iN - 1)
  
  # Wald
  W <- tmp * iT
  wpval <- 1 - pchisq(W, df = iN)
  
  # Modified LR, Jobson/Korkie (1982)
  mX <- matrix(vZm, iT, 1)
  pars <- chol2inv(chol(crossprod(mX))) %*% crossprod(mX, mZ)
  me <- mZ - mX %*% pars
  SigmaR <- crossprod(me) / iT
  LR <- (sum(log(eigen(SigmaR)$values)) - sum(log(eigen(Sigma)$values))) * (iT - iN / 2 - 2)
  lrpval <- 1 - pchisq(LR, df = iN)
  
  # Second Pass with Controls
  gamma <- NULL
  gamsd <- NULL
  
  for (iter in 1:iT) {
    vy <- mZ[iter, ]
    mX <- cbind(1, beta, size_control[iter, ], bm_control[iter, ])  # Add controls to design matrix
    xxinv <- chol2inv(chol(crossprod(mX)))
    xxinvx <- tcrossprod(xxinv, mX)
    
    tmp <- xxinvx %*% vy
    # Store gamma and standard errors
    gamma <- rbind(gamma, c(tmp))
    s2 <- c(crossprod(vy - mX %*% tmp)) / iN
    gamsd <- rbind(gamsd, sqrt(diag(s2 * xxinv)))
  }
  
  # Extract results
  gamma0 <- c(gamma[, 1])  # Intercept
  mrprem <- c(gamma[, 2])  # Market risk premia
  size_coeff <- c(gamma[, 3])  # Size control coefficients
  bm_coeff <- c(gamma[, 4])    # Book-to-market control coefficients
  
  g0sd <- c(gamsd[, 1])   # SE for intercept
  prsd <- c(gamsd[, 2])   # SE for market risk premia
  size_sd <- c(gamsd[, 3])  # SE for size control
  bm_sd <- c(gamsd[, 4])    # SE for book-to-market control
  
  # w-tests
  wgamma0 <- mean(gamma0) / sqrt(sum((gamma0 - mean(gamma0))^2) / iT / (iT - 1))
  wgamma1 <- mean(mrprem) / sqrt(sum((mrprem - mean(mrprem))^2) / iT / (iT - 1))
  w_size <- mean(size_coeff) / sqrt(sum((size_coeff - mean(size_coeff))^2) / iT / (iT - 1))
  w_bm <- mean(bm_coeff) / sqrt(sum((bm_coeff - mean(bm_coeff))^2) / iT / (iT - 1))
  
  return(list(
    alpha = alpha,
    beta = beta,
    J1 = c(J1, jpval),
    W = c(W, wpval),
    mLR = c(LR, lrpval),
    gamma0 = gamma0, g0sd = g0sd,
    mrprem = mrprem, prsd = prsd,
    size_coeff = size_coeff, size_sd = size_sd,
    bm_coeff = bm_coeff, bm_sd = bm_sd,
    wgamma0 = wgamma0, wgamma1 = wgamma1,
    w_size = w_size, w_bm = w_bm
  ))
}


subsample = 1:942
mR = as.matrix(portfolio_m[subsample,25:124])
Rf = as.matrix(portfolio_m[subsample,'Tbill'])
Rm = as.matrix(portfolio_m[subsample,'Market'])
mZ = sweep(mR,1,Rf)
mR - c(Rf)

vZm = Rm - Rf

tmp = EstCAPM(mZ, vZm)
tmp$beta

Fama_MacBeth(mZ, vZm)
ret = FamaMacbeth0(mZ=mZ, vZm=vZm)
ret
#with controls
ret_with_controls <- FamaMacbethWithControls(mZ = mZ, vZm = vZm, size_control = size_control, bm_control = bm_control)

# View results
print(ret_with_controls$size_coeff)  # Coefficients for size control
print(ret_with_controls$bm_coeff)    # Coefficients for book-to-market control
print(ret_with_controls$w_size)      # w-statistic for size
print(ret_with_controls$w_bm)        # w-statistic for book-to-market

length(ret_with_controls$size_coeff)

# Check p-values for the w-tests
p_value_size <- 2 * (1 - pnorm(abs(ret_with_controls$w_size)))
p_value_bm <- 2 * (1 - pnorm(abs(ret_with_controls$w_bm)))

cat("P-value for size control:", p_value_size, "\n")
cat("P-value for book-to-market control:", p_value_bm, "\n")


#now, with ln controls 
print(dim(size_control))  # Should return 942 x 100
print(dim(bm_control))    # Should return 942 x 100
size_control_ln <- log(size_control + 1e-6)  # Add a small constant to avoid log(0)
bm_control_ln <- log(bm_control + 1e-6)     # Same for book-to-market control

print(dim(size_control_ln))  # Should return 942 x 100
print(dim(bm_control_ln))    # Should return 942 x 100


# Run the Fama-MacBeth regression with transformed controls
ret_with_ln_controls <- FamaMacbethWithControls(
  mZ = mZ,
  vZm = vZm,
  size_control = size_control_ln,
  bm_control = bm_control_ln
)

# View results
print(ret_with_ln_controls$size_coeff)  # Coefficients for ln(size_control)
print(ret_with_ln_controls$bm_coeff)    # Coefficients for ln(bm_control)

# Extract w-statistics and p-values for the log controls
w_size_ln <- ret_with_ln_controls$w_size
p_value_size_ln <- 2 * (1 - pnorm(abs(w_size_ln)))

w_bm_ln <- ret_with_ln_controls$w_bm
p_value_bm_ln <- 2 * (1 - pnorm(abs(w_bm_ln)))

# Display results
cat("Test for log(size) control (gamma_size = 0):\n")
cat("w-statistic =", w_size_ln, "\n")
cat("p-value =", p_value_size_ln, "\n\n")

cat("Test for log(book-to-market) control (gamma_bm = 0):\n")
cat("w-statistic =", w_bm_ln, "\n")
cat("p-value =", p_value_bm_ln, "\n")

mean(ret_with_ln_controls$bm_coeff)
mean(ret_with_ln_controls$size_coeff)


-
