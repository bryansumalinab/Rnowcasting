source('functions/cubicbs.R')
source('functions/simR.R')
source('functions/Nowcasting.R')

library("EpiLPS")

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

set.seed(0519)
bias_rep <- mape_rep <- ciwdth_rep <- coverage_rep <- bias.low_rep <- c()
bias_enow <- mape_enow <- ciwdth_enow <- coverage_enow <- bias.low_enow <- c()

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
    
    #### Rt estimation using nowcasted incidence
    nowcast <- Nowcasting(data = data)
    epifit_nowcast <- estimR(incidence = nowcast$cases.now$y, si = si)
    
    value_enow <- subset(epifit_nowcast$RLPS, Time == T.now)
    bias_enow[j] <- value_enow$R - simulepi$Rtrue(T.now)
    mape_enow[j] <- abs((value_enow$R - simulepi$Rtrue(T.now)) / simulepi$Rtrue(T.now))
    ciwdth_enow[j] <- value_enow$Rq0.975 - value_enow$Rq0.025
    coverage_enow[j] <- ifelse(simulepi$Rtrue(T.now) >= value_enow$Rq0.025 & simulepi$Rtrue(T.now) <= value_enow$Rq0.975, 1, 0)
    bias.low_enow[j] <- ifelse(value_enow$R < simulepi$Rtrue(T.now), 1, 0)
    
    #### Rt estimation using reported cases only
    simy_rep <- subset(nowcast$data, Reported == "Reported")
    simy_rep <- aggregate(Cases ~ t, data = simy_rep, FUN = sum)
    simy_rep <- subset(simy_rep, t <= T.now)
    
    epifit_simy_rep <- estimR(incidence = simy_rep$Cases, si = si)
    
    value_rep <- subset(epifit_simy_rep$RLPS, Time == T.now)
    bias_rep[j] <- value_rep$R - simulepi$Rtrue(T.now)
    mape_rep[j] <- abs((value_rep$R - simulepi$Rtrue(T.now)) / simulepi$Rtrue(T.now))
    ciwdth_rep[j] <- value_rep$Rq0.975 - value_rep$Rq0.025
    coverage_rep[j] <- ifelse(simulepi$Rtrue(T.now) >= value_rep$Rq0.025 & simulepi$Rtrue(T.now) <= value_rep$Rq0.975, 1, 0)
    bias.low_rep[j] <- ifelse(value_rep$R < simulepi$Rtrue(T.now), 1, 0)
    
  }, error = function(e) {
    cat("ERROR :", conditionMessage(e), "\n")
  })
  print(j)
}

######## Calculate results ########

######## Using nowcasted incidence ########
# (a) BIAS
bias_enow <- median(bias_enow, na.rm = TRUE)
# (b) BIAS % (Rhat < Rtrue)
biaslow_enow <- mean(bias.low_enow, na.rm = TRUE) * 100
# (c) MAPE
mape_enow <- median(mape_enow, na.rm = TRUE)* 100
# (d) CI width
ciwdth_enow <- median(ciwdth_enow, na.rm = TRUE)
# (e) CI COVERAGE
cov_enow <- mean(coverage_enow, na.rm = TRUE) * 100

######## Using reported cases ########
# (a) BIAS
bias_rep <- median(bias_rep, na.rm = TRUE)
# (b) % (Rhat < Rtrue)
biaslow_rep <- mean(bias.low_rep, na.rm = TRUE) * 100
# (c) MAPE
mape_rep <- median(mape_rep, na.rm = TRUE)* 100
# (d) CI width
ciwdth_rep <- median(ciwdth_rep, na.rm = TRUE)
# (e) CI COVERAGE
cov_rep <- mean(coverage_rep, na.rm = TRUE) * 100

result_enow <- data.frame(bias_enow,
                          biaslow_enow,
                          mape_enow,
                          ciwdth_enow,
                          cov_enow)

result_rep <- data.frame(bias_rep,
                         biaslow_rep,
                         mape_rep,
                         ciwdth_rep,
                         cov_rep)

######## Using nowcasted incidence
result_enow

######## Using reported cases
result_rep
