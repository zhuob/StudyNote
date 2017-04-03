er1 <- read.csv("H:/Projects/CA20737/protocol_revision/Simu_results/type1error_setup1.csv", header =T)
pwr1 <- read.csv("H:/Projects/CA20737/protocol_revision/Simu_results/power_setup1.csv", header =T)
pwr2 <- read.csv("H:/Projects/CA20737/protocol_revision/Simu_results/power_setup2.csv", header =T)


library(dplyr)
table_summary <- function(data){
  result <- data %>% group_by(pt_true, pr_true) %>%
                summarize(r_exact_exact_level = mean(pval_exact <= alpha_exact), 
                          r_norm_exact_level  = mean(pval_norm_apr <= alpha_norm),
                          r_exact_norminal_level = mean(pval_exact <= alpha_nominal), 
                          r_norm_norminal_level = mean(pval_norm_apr <= alpha_nominal))
  return(result)
  }

ty1_rate_setup1 <- table_summary(er1)
pwr_rate_setup1 <- table_summary(pwr1)
pwr_rate_setup2 <- table_summary(pwr2)


library(ggplot2)

ggplot(data = er1 %>% mutate(pt_true =(as.factor(pt_true)))) + 
              stat_qq(aes(sample = pval_exact, color = pt_true, 
                                 linetype  = pt_true), 
                             distribution = qunif, size = 1, geom= "line") + 
                     geom_abline(slope = 1, intercept = 0) + 
      labs(y = "calculated p values", title = "QQ plot of the p values (exact method)", 
           subtitle = "ADA+ rate for reference = 0.04")

ggplot(data = er1 %>% mutate(pt_true =(as.factor(pt_true)))) + 
  stat_qq(aes(sample = pval_norm_apr, color = pt_true, 
              linetype  = pt_true), 
          distribution = qunif, size = 1, geom= "line") + 
  geom_abline(slope = 1, intercept = 0) + 
  labs(y = "calculated p values", title = "QQ plot of the p values (normal approximation)",
       subtitle = "ADA+ rate for reference = 0.04")

