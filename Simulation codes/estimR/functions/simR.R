simR <- function(si, endepi = 50, dist = c("poiss", "negbin"), overdisp = 1,
                 verbose = FALSE, plotsim = FALSE) {
  
  Rtrue <- function(t) {

    A <- c(1.2, 0.75, 1.0)  # Height of the peak (Rt)
    mu <- c(50, 100, 150)    # Peak (time)
    sigma <- c(15,8,10)    # width of the curve (sd)
    
    # Parameters for the decreasing component
    B <- 0.75               # Initial value
    t0 <- 170               # Inflection point
    decrease_rate <- 0.0005   # Decrease rate after t0
    
    # Calculate the bell-shaped components
    bell_shapes <- sapply(1:length(A), function(i) A[i] * exp(-((t - mu[i])^2) / (2 * sigma[i]^2)))
    
    # Calculate the decreasing component
    decreasing_component <- B - (t - t0) * decrease_rate
    
    # Combine the components and sum them up
    R <- sum(bell_shapes) + decreasing_component
    
    return(R)
  }
    
    smax <- length(si)
    mu_y <- c()
    y <- c()
    
    sampling_dist <- match.arg(dist)
    
    for (t in 1:endepi) {
      if (t == 1) {
        mu_y[t] <- 10
        y[t] <- 10
      } else if (t >= 2 && t <= smax) {
        mu_y[t] <- Rtrue(t) *
          sum(rev(y[1:(smax - 1)][1:(t - 1)]) *
                si[1:(smax - 1)][1:(t - 1)])
        if (sampling_dist == "poiss") {
          y[t] <- stats::rpois(n = 1, lambda = mu_y[t])
        } else if (sampling_dist == "negbin") {
          y[t] <- stats::rnbinom(n = 1, mu = mu_y[t], size = overdisp)
        }
      } else if (t > smax && t <= endepi) {
        mu_y[t] <-
          Rtrue(t) * sum(rev(y[(t - smax):(t - 1)]) * si)
        if (sampling_dist == "poiss") {
          y[t] <- stats::rpois(n = 1, lambda = mu_y[t])
        } else if (sampling_dist == "negbin") {
          y[t] <- stats::rnbinom(n = 1, mu = mu_y[t], size = overdisp)
        }
      }
    }
    
    
    #-- Plot simulation-based results
    if (plotsim == TRUE) {
      # Plot 1 (incidence data)
      daysvec <- seq_len(endepi)
      incidata <- data.frame(daysvec = daysvec, y = y)
      
      plotinci <- ggplot2::ggplot(data = incidata,
                                  ggplot2::aes(x = daysvec, y = y)) +
        ggplot2::geom_bar(stat = "identity", width = 0.35, color = "steelblue",
                          fill = "steelblue") +
        ggplot2::xlab("Days") +
        ggplot2::ylab("Incidence") +
        ggplot2::theme(
          axis.title.x = ggplot2::element_text(size = 14),
          axis.title.y = ggplot2::element_text(size = 14),
          axis.text.x = ggplot2::element_text(size = 14),
          axis.text.y = ggplot2::element_text(size = 14)
        )
      
      # Plot 2 (serial interval)
      silen <- seq_len(smax)
      sispec <- data.frame(silen = silen, si = si)
      
      plotsint <- ggplot2::ggplot(data = sispec,
                                  ggplot2::aes(x = silen, y = si)) +
        ggplot2::scale_x_discrete(name = "",
                                  limits = as.character(silen)) +
        ggplot2::geom_bar(stat = "identity", width = 0.35,
                          color = "forestgreen", fill = "forestgreen") +
        ggplot2::xlab("Serial interval index") +
        ggplot2::ylab("Serial interval distribution") +
        ggplot2::theme(
          axis.title.x = ggplot2::element_text(size = 14),
          axis.title.y = ggplot2::element_text(size = 14),
          axis.text.x =  ggplot2::element_text(size = 14),
          axis.text.y =  ggplot2::element_text(size = 14)
        )
      
      # Plot 3 (true R function)
      tdom <- seq(1, endepi, length = 500)
      Rtdom <- sapply(tdom, Rtrue)
      Rfunc <- data.frame(tdom = tdom, Rtdom = Rtdom)
      
      plotR <- ggplot2::ggplot(data = Rfunc, ggplot2::aes(x = tdom, y = Rtdom)) +
        ggplot2::geom_line(size = 1.1) +
        ggplot2::xlab("Days") +
        ggplot2::ylab("R") +
        ggplot2::theme(
          axis.title.x = ggplot2::element_text(size = 14),
          axis.title.y = ggplot2::element_text(size = 14),
          axis.text.x = ggplot2::element_text(size = 14),
          axis.text.y = ggplot2::element_text(size = 14)
        )
      plot_simsummary <- gridExtra::grid.arrange(plotinci, plotsint, plotR,
                                                 nrow = 1)
    }
    
    
    outlist <- list(y = y, mu_y = mu_y, Rtrue = Rtrue,
                    si = si)
    attr(outlist, "class") <- "episim"
    outlist
}
