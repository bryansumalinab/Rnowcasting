source('functions/cubicbs.R')
source('functions/simR.R')
source('functions/Rnowcast.R')

delayprob.list <- list(
  c(0, 0.30, 0.25, 0.20, 0.1, 0.1, 0.05, 0.05),
  c(0.25, 0.20, 0.15, 0.15, 0.1, 0.05, 0.05, 0.05),
  c(0.5, 0.1, 0.1, 0.1, 0.05, 0.05, 0.05, 0.05)
)

# set delay probabilities
p <- delayprob.list[[1]]

# set nowcast Date
T.now <- 124

endepi <- 200
N <- 250 # number of simulations
max.delay <- 7 # maximum delay
si <- c(0.05, 0.05, 0.1, 0.1, 0.1, 0.1, 0.1, 0.05, 0.05, 0.1, 0.1, 0.1)

bias <- mape <- ciwdth <- coverage <- bias.low <- c()

for (j in 1:N) {
  tryCatch({
    # Simulate incidence data
    simulepi <- simR(si = si, endepi = endepi,
                     dist = "negbin", overdisp = 15)
    simcases <- matrix(nrow = endepi, ncol = max.delay + 1)
    for (i in 1:endepi) {
      simcases[i, ] <- rmultinom(n = 1, size = simulepi$y[i], prob = p)
    }
    # Create a data frame of Cases for all t, d combinations
    data <- data.frame(
      t = rep(1:endepi, times = max.delay + 1),
      d = rep(0:max.delay, each = endepi),
      Cases = as.vector(simcases)
    )
    
    data <- data[data$t <= T.now, ]
    Reported <- ifelse(data$t + data$d <= T.now, "Reported", "Not yet reported")
    Reported <- factor(Reported, levels = c("Reported", "Not yet reported"))
    data$Reported <- Reported
    
    # Rnowcasting
    Rnow <- Rnowcast(data = data,
                     serial_interval = si)
    Rnowcast_T.now <- Rnow[Rnow$t == T.now, ]
    
    bias[j] <- Rnowcast_T.now$R_estim - simulepi$Rtrue(T.now)
    mape[j] <- abs((Rnowcast_T.now$R_estim - simulepi$Rtrue(T.now)) / simulepi$Rtrue(T.now))
    ciwdth[j] <- Rnowcast_T.now$CIRt_up - Rnowcast_T.now$CIRt_low
    coverage[j] <- ifelse(simulepi$Rtrue(T.now) >= Rnowcast_T.now$CIRt_low & simulepi$Rtrue(T.now) <= Rnowcast_T.now$CIRt_up, 1, 0)
    bias.low[j] <- ifelse(Rnowcast_T.now$R_estim < simulepi$Rtrue(T.now), 1, 0)
    
  }, error = function(e) {
    cat("ERROR :", conditionMessage(e), "\n")
  })
  print(j)
}

median.bias <- median(bias, na.rm = TRUE) # BIAS
mean.biaslow <- mean(bias.low, na.rm = TRUE) * 100 # % (Rhat < Rtrue)
median.mape <- median(mape, na.rm = TRUE) # MAPE
median.ciwdth <- median(ciwdth, na.rm = TRUE) # CI width
mean.cov <- mean(coverage, na.rm = TRUE) * 100 # CI coverage

result <- data.frame(median.bias, 
                     mean.biaslow, 
                     median.mape, 
                     median.ciwdth, 
                     mean.cov)
result