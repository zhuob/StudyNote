
source("H:/Projects/CA20737/protocol_revision/FM_score.R")




## this is the Chan's exact method. max(p_exact) is the corresponding p-value
tail_prob1 <- function(x1, x2, n1, n2, delta, p2_search, z1){
  
  niter <- length(p2_search)
  z_test <- test_stat_fm(x1, x2, n1, n2, delta)$z_eq
  
  z1 <- z_table(n1, n2, delta)  # calculate z for each combination of i and j
  z_index <- z1 <= z_test   # a matrix of TRUE or FALSE 
  
  
  
  p_exact <- rep(NA, niter)
  for(i in 1:niter){
    p_exact[i] <- exact_prob(n1, n2, p2_search[i], delta = delta, z_index = z_index)
  }
  
  return(data.frame(p2 = p2_search, p_exact = p_exact))
  
}



## the true significance levels, where the critical value is the normal quantile (that's why it's called normal approxmiation) 
true_alpha_normal_app1 <- function( n1, n2, delta, p2, alpha = 0.05, z1){
  
  z_critical <- qnorm(alpha)
  z_index <- z1 <= z_critical  # the critical value if normal approximation is used
  
  # summing up the joint probabilities for all combinations of (i, j), such that z[i, j] < z_critical
  true_alpha <- exact_prob(n1, n2, p2, delta, z_index) 
  
  return(true_alpha)
  
}

## the true significant level for exact test. I removed the calculation of z1 to make the simulation run faster
true_alpha_exact1 <- function(n1, n2, delta, p2, alpha = 0.05, z1){
  library(dplyr)

  # the joint probability for every combination of x=i and y = j;
  d1 <- dbinom(0:n1, size = n1, prob = delta + p2)
  d2 <- dbinom(0:n2, size = n2, prob = p2)
  # probability matrix, prob_mat[i, j] = p(x = i|n1, p2 + delta)*p(y = j|n2, p2)
  prob_mat <- d1 %*% t(d2) 
  
  # find the critical z value such that the tail is smaller than and 
  # as close as possible to the desired significance level
  z1_vec <- as.vector(z1)  # vectorize the table by column
  prob_vec <- as.vector(prob_mat)
  
  
  # find the more "extreme z values" such that their corresponding tail p value is <= alpha 
  ztab <- data.frame(z1_vec, prob_vec) %>% arrange(z1_vec) %>% # sort the z values
    mutate(cumprob = cumsum(prob_vec)) %>%    # calculate cummulative probabilities
    mutate(indicator = cumprob <= alpha)      # decide a cut off z value
  z_critical <- ztab %>% filter(indicator == TRUE)          # the cumprob in last row is the corresponding exact significant level
  
  row_id <- nrow(z_critical)
  true_alpha <- z_critical[row_id, 3]
  return(true_alpha)
}



type1_error_simu <- function(p1, p2, n_cohort1= 46, delta, alpha = 0.05){
  
  # step 1: p1 and p2 specified
  n1 <- n_cohort1/2              # number of subjects in each Test/Reference group
  
  
  # step 2: simulate x1 and y1
  x1 <- rbinom(1, n1, prob = p1)  # number of ADA+ in Test group 
  y1 <- rbinom(1, n1, prob = p2)  # number of ADA+ in Reference group 
  
  # using Jeffreys prior, the posterior mean is (x1 + x2 + 0.5)/(n1 + 1)
  p_jeffreys <- (x1 + y1 + 0.5)/(n_cohort1 + 1)
  p_ini <- p_jeffreys
  
  # Step 3: calculate the sample size needed for cohort 2
  n_calc_size <- sample_size(p1 = p_ini, p2 = p_ini, theta = 1,  delta, alpha = alpha, beta = 0.2)
  n2 <- ceiling(n_calc_size$n/2) - n1  # subjects needed in cohrt 2
  
  # step 4: Simulate the outcome of Cohort 2
  x2 <- rbinom(1, n2, p1)
  y2 <- rbinom(1, n2, p2)
  
  # step 5: unblind the data
  x <- x1 + x2 
  y <- y1 + y2
  n_per_arm <- n1 + n2
  
  p2_search <- seq(0.01, min(1-p1, 1-p2, 1-delta), by = 0.001)
  
  # the observed z_table for every combination of x=i and y = j
  z1 <- z_table(n_per_arm, n_per_arm, delta)
  
  ## using Chan's exact method (1998) to calcualte the p-value
  result <- tail_prob1(x, y, n_per_arm, n_per_arm, delta = delta,
                       p2_search = p2_search, z1 = z1)
  pval_exact <- max(result$p_exact)
  
  ## using normal approximation 
  pval_norm_apr <- pnorm(test_stat_fm(x, y, n_per_arm, n_per_arm, delta = delta)$z_eq)
  
  
  ## get the actual significance level 
  alpha_norm <- true_alpha_normal_app1(n1 = n_per_arm, n2 = n_per_arm,
                                       delta = delta, p2 = p2, alpha = alpha, z1 = z1) 
  alpha_exact <- true_alpha_exact1(n1 = n_per_arm, n2 = n_per_arm, 
                                   delta = delta, p2 = p2, alpha = alpha, z1 = z1) 
  
  result <- data.frame(n_per_arm, n1, n2, x1, x2, y1, y2,  
                       pval_exact, alpha_exact, pval_norm_apr, alpha_norm,
                       pt_true = p1, pr_true = p2, alpha_nominal = alpha)
  return(result)
}





report_result <- function(pt, pr, k, alpha_nominal = 0.05, simu_type = "Type1error"){
 
  fre <- length(pt)  # how many parameters there are
  sim1 <- data.frame(matrix(ncol =14, nrow = k*fre))
  names(sim1) <- c("n_per_arm", "n1", "n2", "x1", "x2", "y1", "y2",
                   "pval_exact",  "alpha_exact", "pval_norm_apr", "alpha_norm", 
                   "pt_true", 	"pr_true", 	"alpha_nominal")
  
  for (s in 1:fre){
    for (i in 1:k){
      row_in <- (s-1)*k + i 
      temp1 <- type1_error_simu(pt[s], pr[s], delta = 0.1, alpha = alpha_nominal)
      sim1[row_in, ] <- temp1 
      cat('\r', row_in)
    }

  }
 
  
  saveobj <- paste("H:/Projects/CA20737/protocol_revision/Simu_results/", simu_type, ".csv", sep = "")
  write.csv(sim1, file = saveobj,  row.names = FALSE)
 
}



### run the simulations 

pt <- seq(0.14, 0.24, by = 0.01); pr <- rep(0.04, length(pt))  # this corresponds to type 1 error
k <- 1000  # number of simulations for each pt and pr

report_result(pt, pr, k, alpha_nominal = 0.05, simu_type = "type1error_setup1")


pt <- seq(0.14, 0.24, by = 0.01); pr <- pt - 0.1  # this corresponds to type 1 error
report_result(pt, pr, k, alpha_nominal = 0.05, simu_type = "type1error_setup2")


## power

pt <- seq(0.01, 0.10, by = 0.01); pr <- rep(0.04, length(pt))  
report_result(pt, pr, k, alpha_nominal = 0.05, simu_type = "power_setup1")


pr <- seq(0.01, 0.10, by = 0.01); pt <- pr  + 0.05
report_result(pt, pr, k, alpha_nominal = 0.05, simu_type = "power_setup2")
















