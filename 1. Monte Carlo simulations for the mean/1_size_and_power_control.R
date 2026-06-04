#size vs power

#Mean-type statistic
#X = Data matrix, mu = Null hypothesis, Rn = numbers of resamples
Mean_stat <- function(X,mu,Rn){
  n <- length(X)
  
  sigma_hat <- var(X) * (n-1)/n
  
  indices <- sample(1:n,Rn,replace=TRUE)
  W <- (X[indices] - mu) / sqrt(sigma_hat)
  
  T_M_tilde <- sum(W) / sqrt(Rn)
  T_M <- T_M_tilde^2
  
  return(T_M)
}

#U-type statistic
#X = Data, mu = Null hypothesis, Rn = numbers of resamples  
U_stat <- function(X,mu,Rn){
  n <- length(X)
  
  sigma_hat <- var(X) * (n-1)/n
  
  i <- sample.int(n, Rn, replace = TRUE)
  j <- sample.int(n, Rn, replace = TRUE)
  
  duplicates <- which(i == j)
  while (length(duplicates)) {
    j[duplicates] <- sample.int(n, length(duplicates), replace = TRUE)
    duplicates <- duplicates[i[duplicates] == j[duplicates]]
  }
  
  W <- (X[i] - mu) * (X[j] - mu) / sigma_hat
  T_U <- sum(W) / sqrt(Rn)
  
  return(T_U)
}

#Simulation

#parameters
alpha <- 0.05
iterations <- 2000
rho <- 1
mu_0 <- 0
mu_alt <- 0.5

#clustering scenarios
scenarios <- list(c(100,2),
             c(200,5),
             c(5,200),
             c(30,30))

all_type1 <- matrix(NA, nrow = length(scenarios), ncol = 3)
all_power <- matrix(NA, nrow = length(scenarios), ncol = 3)
scenario_labels <- character(length(scenarios))

for (idx in seq_along(scenarios)){
  num_clusters <- scenarios[[idx]][1]
  cluster_size <- scenarios[[idx]][2]
  n <- num_clusters*cluster_size
  cluster_id   <- rep(1:num_clusters, each = cluster_size)
  Rn_M         <- ceiling(sqrt(n))
  Rn_U         <- ceiling((n / 2)^(4 / 3))
  
  scenario_labels[idx] <- paste0("(", num_clusters, ",", cluster_size, ")")
  cat("Running scenario", scenario_labels[idx],
      "| n =", n, "| Rn_M =", Rn_M, "| Rn_U =", Rn_U, "\n")
  
  results <- data.frame(
    Type1_Mean = logical(iterations), Type1_U = logical(iterations),
    Type1_t    = logical(iterations), Power_Mean = logical(iterations),
    Power_U    = logical(iterations), Power_t    = logical(iterations)
  )
  
  for (i in 1:iterations){
    #Generate clustered data
    alpha_c_null <- rnorm(num_clusters, mu_0,   sqrt(rho))
    eps_null     <- rnorm(n)
    X_null       <- alpha_c_null[cluster_id] + eps_null
    
    alpha_c_alt  <- rnorm(num_clusters, mu_alt, sqrt(rho))
    eps_alt      <- rnorm(n)
    X_alt        <- alpha_c_alt[cluster_id] + eps_alt
    
    # Type-I error
    tm_null <- Mean_stat(X_null, mu = mu_0, Rn = Rn_M)
    tu_null <- U_stat(X_null,    mu = mu_0, Rn = Rn_U)
    t_null  <- t.test(X_null, mu = mu_0)$p.value
    
    results$Type1_Mean[i] <- tm_null > qchisq(1 - alpha, df = 1)
    results$Type1_U[i]    <- tu_null > qnorm(1 - alpha)
    results$Type1_t[i]    <- t_null  < alpha
    
    # Power
    tm_alt <- Mean_stat(X_alt, mu = mu_0, Rn = Rn_M)
    tu_alt <- U_stat(X_alt,    mu = mu_0, Rn = Rn_U)
    t_alt  <- t.test(X_alt, mu = mu_0)$p.value
    
    results$Power_Mean[i] <- tm_alt > qchisq(1 - alpha, df = 1)
    results$Power_U[i]    <- tu_alt > qnorm(1 - alpha)
    results$Power_t[i]    <- t_alt  < alpha
  }
  
  all_type1[idx, ] <- c(mean(results$Type1_Mean),
                          mean(results$Type1_U),
                          mean(results$Type1_t))
  all_power[idx, ] <- c(mean(results$Power_Mean),
                          mean(results$Power_U),
                          mean(results$Power_t))
}

#Plot
scenario_cols <- c("#B03E3E", "#2E7D5C", "#B07B20", "#2A5A8C")

stat_labels <- c("Mean-type", "U-type", "t-test")
stat_pch    <- c(15, 17, 19)   

xlim_range <- c(0, max(all_type1, na.rm = TRUE) * 1.35)

par(mar = c(5.5, 4.5, 2, 1), mgp = c(3, 0.7, 0), tcl = -0.3)

plot(NA,
     xlab = "Type I error (rejection rate under H0)",
     ylab = "Power (rejection rate under H1)",
     xlim = xlim_range,
     ylim = c(0, 1),
     las  = 1,          
     bty  = "o")        

grid(nx = NULL, ny = NULL, col = "grey85", lty = 1, lwd = 0.7)

abline(v = 0.05, col = "#A01010", lty = 2, lwd = 1.2)

for (s_idx in seq_along(scenarios)) {
  col    <- scenario_cols[s_idx]
  type1s <- all_type1[s_idx, ]
  powers <- all_power[s_idx, ]
  
  
  ord <- order(type1s)
  lines(type1s[ord], powers[ord],
        col = col, lty = 3, lwd = 0.9)
  
  points(type1s, powers,
         pch = stat_pch,
         col = col,
         bg  = col,
         cex = 1.5,
         lwd = 1.8)
}


legend("topright",
       legend = scenario_labels,
       col    = scenario_cols,
       pch    = 19,
       pt.cex = 1.3,
       lty    = 3,
       lwd    = 0.9,
       title  = expression(bold("Scenario ("*n[c]*", "*c[s]*")")),
       bty    = "n",
       cex    = 0.82)

legend("bottomright",
       legend = stat_labels,
       pch    = stat_pch,
       col    = "grey30",
       pt.cex = 1.3,
       pt.lwd = 1.8,
       title  = expression(bold("Statistic")),
       bty    = "n",
       cex    = 0.82)



