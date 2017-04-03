## sample size calculation using FM-score method (Farrington and Manning 1990)


# http://www.ncss.com/wp-content/themes/ncss/pdf/Procedures/PASS/Non-Inferiority_Tests_for_Two_Proportions.pdf
## this function is used to calculate sample size 

sample_size <- function(p1, p2, theta = 1,  delta, alpha = 0.05, beta = 0.2){
  # parameters
  
  # Group        Success     Failure      Total
  # Treatment    n11          n12          n1         
  # Control      n21          n22          n2
  # Total        n.1          n.2           n
  
  # p10: group 1 proportion tested by the null 
  # delta0: non-inferiority margin
  # p1: binomial proportions = n11/n1
  # p2: binomial proportions = n21/n2
  # p : overall proportion  = m1/n
  # delta = p1 - p2
  
  #         H0: p10 - p2 >= delta  vs  H1: p10 - p2 < delta
  
  a <- 1 + theta
  b <- -(1 + theta + p1 + theta*p2 + delta*(theta + 2))
  c <- delta^2 + delta*(2*p1 + theta + 1) + p1 + theta*p2
  d <- -p1*delta*(1+delta)
  
  v <- (b/(3 * a))^3 - b * c/6/a^2 + d/2/a
  u <- (sign(v) + (v == 0)) * sqrt((b/3/a)^2 - c/3/a)
  
  # to avoid acos(1+1e-15) == NaN
  s <- v/(u^3)
  s[s>1] <- 1
  
  w <- (pi + acos(s))/3 
  
  p1_tilt <- 2*u*cos(w) - b/(3*a)
  p2_tilt <- p1_tilt - delta
  
  p2_tilt[p2_tilt < 0] <- 0  # 0.1 - 0.1 = -1.94289e-16
  # print(c(p1_tilt, p2_tilt))
  se0 <- sqrt( (p1_tilt*(1-p1_tilt) + p2_tilt*(1-p2_tilt)/theta)*(1 + theta))
  se1 <- sqrt((p1*(1-p1) + p2*(1-p2)/theta)*(1 + theta))
  
  z_alpha <- qnorm(1 - alpha)  
  z_beta <- qnorm(1 - beta)
  
  nsize <- ( (z_alpha*se0 + z_beta*se1)/(p1 - p2 - delta) )^2
  # print(nsize)
  #  print(c(se0, se1))
  result <- data.frame(n = nsize, n1 = nsize/(theta + 1), n2 = nsize*theta/(theta + 1),
                       p1_tilt, p2_tilt)
  
  
  return(result)
}



## estimate the MLE of two proportions using formula given by Farrington and Manning
## Tomorrow --> handle the cases of (0, 0), (0, 1), (1, 0), (1, 1)
prop_mle <- function(p1_hat, p2_hat, delta, theta){
  
  # the boundary p values are handled based on Chan's 1999 paper: 
  # TEST-BASED EXACT CONFIDENCE INTERVALS FOR THE DIFFERENCES OF TWO BINOMIAL PROPORTIONS
    
    a <- 1 + theta
    b <- -(1 + theta + p1_hat + theta*p2_hat + delta*(theta + 2))
    c <- delta^2 + delta*(2*p1_hat + theta + 1) + p1_hat + theta*p2_hat
    d <- -p1_hat*delta*(1+delta)
    
    v <- (b/(3 * a))^3 - b * c/6/a^2 + d/2/a
    u <- (sign(v) + (v == 0)) * sqrt((b/3/a)^2 - c/3/a)
    
    # to avoid acos(1+1e-15) == NaN
    s <- v/(u^3)
    s[s>1] <- 1
  
    w <- (pi + acos(s))/3 
    
    p1_tilt <- 2*u*cos(w) - b/(3*a)
    p2_tilt <- p1_tilt - delta
    
    p2_tilt[p2_tilt < 0] <- 0  # 0.1 - 0.1 = -1.94289e-16
  
  
  return(c(p1_tilt, p2_tilt))
}



## test statistics by formula 1 of Farrington and Manning

# the FM-statistic 
test_stat_fm <- function(x1, x2, n1, n2, delta){
  # x1: the number of events in treatment
  # x2: the nubmer of events in control
  # n1: total number of events in treatment
  # n2: total number of events in control
  
# this statistic is calculated based on the hypothesis that 
#     H0: p1-p2 >= delta    vs    H1: p1 - p2 < delta
  
  p1_hat <- x1/n1
  p2_hat <- x2/n2
  n <- n1 + n2
  theta <- n2/n1
  
  ## Tomorrow --> handle the cases of (0, 0), (0, 1), (1, 0), (1, 1)
  # the boundary z values are handled based on Chan's 1999 paper: 
  # TEST-BASED EXACT CONFIDENCE INTERVALS FOR THE DIFFERENCES OF TWO BINOMIAL PROPORTIONS
  

  
  if( (p1_hat == 0 & p2_hat == 0 ) | (p1_hat == 1 & p2_hat == 1 ) ){
    z_eq <- 0
    p1_tilt <- p2_tilt <- 0
    }
  else {
      if(p1_hat ==0 & p2_hat == 1) { 
          p1_hat <- 1/(2*n1); p2_hat = 1- 1/(2*n2) }
    
      if (p1_hat == 1 & p2_hat == 0){ 
          p1_hat <- 1- 1/(2*n1); p2_hat = 1/(2*n2) }
  
  p_mle <- prop_mle(p1_hat = p1_hat, p2_hat = p2_hat, delta = delta, theta = theta)
  p1_tilt <- p_mle[1]; 
  p2_tilt <- p_mle[2]
  
  s_eq <- sqrt( (p1_tilt*(1-p1_tilt) + p2_tilt*(1-p2_tilt)/theta)*(1 + theta))/sqrt(n)
  
  z_eq <- (p1_hat-p2_hat - delta)/s_eq
  
  }
  
  
  return(list(p1=p1_tilt, p2= p2_tilt, z_eq = z_eq))
  
}

# the observed test statistic

test_stat_obs <- function(x1, x2, n1, n2){
  
  p1_hat <- x1/n1
  p2_hat <- x2/n2
  n <- n1 + n2
  theta <- n2/n1
  
  s_obs <- sqrt( (p1_hat*(1-p1_hat) + p2_hat*(1-p2_hat)/theta )*(1+theta))/sqrt(n)
  
  z_obs <- (p1_hat-p2_hat)/s_obs

  return(list(p1_hat = p1_hat, p2_hat = p2_hat, z_obs = z_obs))
}



# the sum of probabilites whose combination of i and j results a test statistic that is less than or equal to z_eq 
exact_prob <- function(n1, n2, p2, delta, z_index){
  
  d1 <- dbinom(0:n1, size = n1, prob = delta + p2)
  d2 <- dbinom(0:n2, size = n2, prob = p2)
  
  # probability matrix, prob_mat[i, j] = p(x = i|n1, p2 + delta)*p(y = j|n2, p2)
  prob_mat <- d1 %*% t(d2)  
  
  p_tail <- sum(prob_mat[z_index])
  
  return(p_tail)
  
}


# the z table for all possible combinations of x = i and y = j
# the z values are calculated based on FM-Score test.
z_table <- function(n1, n2, delta){
  z_1 <- matrix(NA, nrow= n1 + 1, ncol = n2 + 1)
  
  for(i in 1:nrow(z_1)){  # calculate the z_eq for each pair of (i, j)
    for(j in 1:ncol(z_1)){
      temp <- test_stat_fm(i-1, j-1, n1, n2, delta)
      z_1[i, j] <- temp$z_eq
    }
    
  }
  
  return(z_1)
  
}



## this is the Chan's exact method. max(p_exact) is the corresponding p-value
tail_prob <- function(x1, x2, n1, n2, delta, p2_search){
  
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
true_alpha_normal_app <- function( n1, n2, delta, p2, alpha = 0.05){
  
  z1 <- z_table(n1, n2, delta)
  z_critical <- qnorm(alpha)
  z_index <- z1 <= z_critical  # the critical value if normal approximation is used
  
  # summing up the joint probabilities for all combinations of (i, j), such that z[i, j] < z_critical
  true_alpha <- exact_prob(n1, n2, p2, delta, z_index) 
  
  return(true_alpha)
  
}

## the true significant level for exact test.
true_alpha_exact <- function(n1, n2, delta, p2, alpha = 0.05){
  library(dplyr)
  # the observed z_table for every combination of x=i and y = j
  z1 <- z_table(n1, n2, delta)
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


