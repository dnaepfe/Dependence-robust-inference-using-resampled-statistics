#power curve (power vs sample size)

#sanity check
iid = FALSE

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

#Parameters
num_clusters <- c(10,20,30,50,100,200,400)
cluster_size <- 5

alpha <- 0.05
rho <- 1
iterations <- 2000

mu_0 <- 0
mu_alt <- 0.5

results <- data.frame(n=num_clusters * cluster_size,
                      Power_Mean = NA_real_,
                      Power_U    = NA_real_,
                      Power_t    = NA_real_
)

for (k in seq_along(num_clusters)){
  num_cluster <- num_clusters[k]
  n <- num_cluster * cluster_size
  cluster_id <- rep(1:num_cluster, each=cluster_size)
  Rn_M <- ceiling(sqrt(n))
  Rn_U <- ceiling((n / 2)^(4 / 3)) 
  
  rej_Mean <- logical(iterations)
  rej_U    <- logical(iterations)
  rej_t    <- logical(iterations)
  
  for (i in 1:iterations){
    
    #generate clustered data
    alpha_f_alt <- rnorm(num_cluster, mu_alt, sqrt(rho))
    eps_alt <- rnorm(n)
    X_alt <- alpha_f_alt[cluster_id] + eps_alt
    
    if(iid){
      X_alt <- rnorm(n,mu_alt,sqrt(rho))
    }
    
    tm_alt <- Mean_stat(X_alt, mu = mu_0, Rn = Rn_M)
    tu_alt <- U_stat(X_alt,    mu = mu_0, Rn = Rn_U)
    t_p    <- t.test(X_alt, mu = mu_0)$p.value
    
    rej_Mean[i] <- tm_alt > qchisq(1 - alpha, df = 1)
    rej_U[i]    <- tu_alt > qnorm(1 - alpha)
    rej_t[i]    <- t_p    < alpha
  }
  
  results$Power_Mean[k] <- mean(rej_Mean)
  results$Power_U[k] <- mean(rej_U)
  results$Power_t[k] <- mean(rej_t)
  
  cat(sprintf("n = %4d | Mean: %.3f | U: %.3f | t: %.3f\n",
              n,
              results$Power_Mean[k],
              results$Power_U[k],
              results$Power_t[k]))
}

#plot 
plot(results$n, results$Power_Mean,
     type = "b", pch = 19, col = "#2A5A8C", lty = 3,
     xlab = "Sample size n",
     ylab = "Power",
     ylim = c(0, 1))
grid(nx = NULL, ny = NULL, col = "grey85", lty = 1, lwd = 0.7)
lines(results$n, results$Power_U,
      type = "b", pch = 17, col = "#2E7D5C", lty = 3)
lines(results$n, results$Power_t,
      type = "b", pch = 15, col = "#B03E3E", lty = 3)
abline(h = alpha, col = "#A01010", lty = 2)
legend("bottomright",
       legend = c("Mean-type", "U-type", "t-test"),
       col = c("#2A5A8C", "#2E7D5C", "#B03E3E"),
       pch = c(19, 17, 15), lty = 2, bty = "n")
