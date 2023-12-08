source("prowl.R")
library(tidyverse)
stard_wide <- read.csv("Data/stard_wide.csv") # Local

# Keep only level 2 and the second visits of level 1
# Remove the dem variables which are majorly missing
# If someone is in level 2.1, take that treatment rather than the level 2 treatment
stard_sub <-  stard_wide %>%
                select(-contains("L4."), 
                       -contains("L3.")) %>%
                filter(!is.na(txassign.L2.1) | 
                         !is.na(txassign.L2.2) | 
                         !is.na(txassign.L2.1.1) | 
                         !is.na(txassign.L2.1.2)) %>% 
                select(-c(white, black, marital, educat, empl, publica, privins)) %>% 
                mutate(
                  db.0 = days_baseline.L1.2,
                  db.1 = ifelse(is.na(days_baseline.L2.1.1), days_baseline.L2.1, days_baseline.L2.1.1),
                  db.2 = ifelse(is.na(days_baseline.L2.1.2), days_baseline.L2.2, days_baseline.L2.1.2),
                  week.0 = week.L1.2,
                  week.1 = ifelse(is.na(week.L2.1.1), week.L2.1, week.L2.1.1),
                  week.2 = ifelse(is.na(week.L2.1.2), week.L2.2, week.L2.1.2),
                  trtmt.0 = trtmt.L1.2,
                  trtmt.1 = ifelse(is.na(trtmt.L2.1.1), trtmt.L2.1, trtmt.L2.1.1),
                  qccur_r.0 = qccur_r.L1.2,
                  qccur_r.1 = ifelse(is.na(qccur_r.L2.1.1), qccur_r.L2.1, qccur_r.L2.1.1),
                  qccur_r.2 = ifelse(is.na(qccur_r.L2.1.2), qccur_r.L2.2, qccur_r.L2.1.2),
                  qscur_r.0 = qscur_r.L1.2,
                  qscur_r.1 = ifelse(is.na(qscur_r.L2.1.1), qscur_r.L2.1, qscur_r.L2.1.1),
                  qscur_r.2 = ifelse(is.na(qscur_r.L2.1.2), qscur_r.L2.2, qscur_r.L2.1.2),
                  qcimp_r.0 = qcimp_r.L1.2,
                  qcimp_r.1 = ifelse(is.na(qscur_r.L2.1.1),qcimp_r.L2.1,qcimp_r.L2.1.1),
                  qcimp_r.2 = ifelse(is.na(qscur_r.L2.1.2),qcimp_r.L2.2,qcimp_r.L2.1.2),
                  action_call.0 = action_call.L1.2,
                  action_call.1 = ifelse(is.na(qscur_r.L2.1.1),action_call.L2.1,action_call.L2.1.1),
                  action_call.2 = ifelse(is.na(qscur_r.L2.1.2),action_call.L2.2,action_call.L2.1.2),
                  txassign.0 = txassign.L1.1,
                  txassign.1 = ifelse(is.na(qscur_r.L2.1.1),txassign.L2.1,txassign.L2.1.1),
                  medsw.0 = medswl2.L1.2,
                  medsw.1 = ifelse(is.na(qscur_r.L2.1.1),medswl2.L2.1,medswl2.L2.1.1),
                  medsw.2 = ifelse(is.na(qscur_r.L2.1.2),medswl2.L2.2,medswl2.L2.1.2),
                  medaug.0 = medaugl2.L1.2,
                  medaug.1 = ifelse(is.na(qscur_r.L2.1.1),medaugl2.L2.1,medaugl2.L2.1.1),
                  medaug.2 = ifelse(is.na(qscur_r.L2.1.2),medaugl2.L2.2,medaugl2.L2.1.2),
                  cogsw.0 = cogswl2.L1.2,
                  cogsw.1 = ifelse(is.na(qscur_r.L2.1.1),cogswl2.L2.1,cogswl2.L2.1.1),
                  cogsw.2 = ifelse(is.na(qscur_r.L2.1.2),cogswl2.L2.2,cogswl2.L2.1.2),
                  cogaug.0 = cogaugl2.L1.2,
                  cogaug.1 = ifelse(is.na(qscur_r.L2.1.1),cogaugl2.L2.1,cogaugl2.L2.1.1),
                  cogaug.2 = ifelse(is.na(qscur_r.L2.1.2),cogaugl2.L2.2,cogaugl2.L2.1.2),
                  grseb.0 = grseb.L1.2,
                  grseb.1 = ifelse(is.na(qscur_r.L2.1.1),grseb.L2.1,grseb.L2.1.1),
                  grseb.2 = ifelse(is.na(qscur_r.L2.1.2),grseb.L2.2,grseb.L2.1.2),
                  fisin.0 = fisin.L1.2,
                  fisin.1 = ifelse(is.na(qscur_r.L2.1.1),fisin.L2.1,fisin.L2.1.1),
                  fisin.2 = ifelse(is.na(qscur_r.L2.1.2),fisin.L2.2,fisin.L2.1.2)
                ) %>% 
                mutate(Q1 = ifelse(is.na(qccur_r.1), qscur_r.1, qccur_r.1),
                       Q2 = ifelse(is.na(qccur_r.2), qscur_r.2, qccur_r.2),
                       S1 = ifelse(is.na(qcimp_r.1), qcimp_r.0, qcimp_r.1)/100,
                       MS = medsw.1,
                       MA = medaug.1,
                       CS = cogsw.1,
                       CA = cogaug.1) %>%
                select(-contains("L1"), 
                       -contains("L2"), 
                       -contains("week"), 
                       -contains("trtmt"), 
                       -contains("cur_r"),
                       -contains("action_call"),
                       -contains("imp_r"),
                       -contains("medsw"),
                       -contains("medaug"),
                       -contains("cogsw"),
                       -contains("cogaug"),
                       -c("grseb.2", "fisin.2")) %>% 
                filter(complete.cases(.))


## Randomization Probabilities
# We are considering single versus combination therapies
# Combination therapies are available to those who said that they were okay with augmentation
# The only individuals capable of receiving a combined therapy are those with one of [MA, CA] == 1
# The only individuals capable of receiving a mono therapy are those with one of [MS, CS] == 1

# As a result, we consider only the subset upon which there is no positivity violation
# That is, filter to individuals who have both: MA+CA < 4 and MS+CS < 4.
stard_sub %>% filter(MA == 1 | CA == 1, MS == 1 | CS == 1) %>% nrow()
stard_sub %>% filter(MS == 1 | (CS == 1 & (MA == 1 | CA == 1))) %>% nrow()

# There is only n = 109 who satisfy the requirements we have set out.
# This may not be ideal. We can certainly try though? 
# If it is bad then we can instead turn to doing SSRIs vs. Non-SSRIs.
owl_stard <- stard_sub %>% filter(MA == 1 | CA == 1, MS == 1 | CS == 1)

## Begin PrOWL
library(personalized)

# Randomization Probabilities
owl_stard <- stard_sub %>% 
              mutate(A = ifelse(grepl("\\+", txassign.1), 1, -1)) %>% 
              mutate(totalTrtOptions = 3 * (MS == 1) + 1 * (CS == 1) + 2 * (MA == 1) + 1 * (CA == 1),
                     totalMonoOptions = 3 * (MS == 1) + 1 * (CS == 1),
                     totalComboOptions = 2 * (MA == 1) + 1 * (CA == 1)) %>%
              mutate(pA = totalComboOptions / totalTrtOptions,
                     Y = 27-Q2) %>%
              filter(pA != 0, pA != 1) %>% 
              select(-contains("total"), -contains('txassign'), -c("MS", "MA", "CS", "CA"))

## Set a seed for reproducibility of PrOWL
set.seed(31415)

d <- 7
diam <- sqrt(1^2 + (900-216)^2 + 6^2 + 6^2 + 27^2 + 1^2)
xi <- 26
zeta <- 0 
cLower <- 0.2
cUpper <- 0.75
L <- 1
cost <- 0.01
n <- 109
Delta <- 3
epsilon <- 1

privacy_parms <- compute_prowl_parms(kernel_parms = list(F = Inf, 
                                                         kappa = 1, 
                                                         sigma_p = sqrt(d), 
                                                         diam = diam),
                                     xi,
                                     zeta,
                                     cLower,
                                     cUpper,
                                     L,
                                     cost,
                                     d,
                                     n = n,
                                     beta = NULL,
                                     Delta = Delta,
                                     epsilon = epsilon)

# The Formula
stard_form <- A ~ gender + interview_age + grseb.1 + fisin.1 + Q1 + S1

priv_est <- PrOWL(df_formula = stard_form, 
                     data = owl_stard, 
                     cost = cost,
                     kernel_type = "rbfdot", 
                     treatment_var = "A", 
                     reward_var = "Y", 
                     prob_var = "pA",
                     lambda = privacy_parms$lambda,
                     epsilon = epsilon,
                     d = privacy_parms$set_parms$d,
                     d_tilde = privacy_parms$d_tilde) 

# Try to get OWL working well
standard_owl_with_cv <- OWL(df_formula = stard_form, 
                            data = owl_stard, 
                            cost = seq(0.01, 1, length.out = 10),
                            cv_folds = 10,
                            kernel_type = "rbfdot", 
                            treatment_var = "A", 
                            reward_var = "Y",
                            prob_var = "pA")

private_predictions <- priv_est$predict_treatment(model.matrix(stard_form, owl_stard))
non_private_predictions <- standard_owl_with_cv$predict_treatment(model.matrix(stard_form, owl_stard))

sum(non_private_predictions == private_predictions)

## Non parametric bootstrap for value estimation
set.seed(9265)
bs_reps <- 10000
value_results <- matrix(nrow = bs_reps, ncol = 2) 

for(ii in 1:bs_reps) {
  bs_sample <- sample(1:nrow(owl_stard), nrow(owl_stard), replace = TRUE)
  
  bs_df <- owl_stard[bs_sample, ]
  
  bs_preds_priv <- priv_est$predict_treatment(model.matrix(stard_form, bs_df))
  bs_preds_nonpriv <- standard_owl_with_cv$predict_treatment(model.matrix(stard_form, bs_df))
  
  value_results[ii,] <- c(estimated_value(bs_preds_priv, bs_df$A, bs_df$Y, bs_df$pA),
                            estimated_value(bs_preds_nonpriv, bs_df$A, bs_df$Y, bs_df$pA))
  
}

round(quantile(value_results[,1], probs = c(0.025, 0.975)), 2)
round(quantile(value_results[,2], probs = c(0.025, 0.975)), 2)
round(priv_est$value, 2)
round(standard_owl_with_cv$value, 2)

## Because these results look very promising, it is worth considering whether we got a lucky seed
## or similar. To assess this we will consider nonparametric Bootstrap for the Full Estimation procedure
## These results seem to confirm the single analysis was not a fluke with these data
## Not reporting in the analysis in the paper, but contained here for confirmation
# set.seed(358979)
# bs_all_reps <- 500
# all_results <- matrix(nrow = bs_all_reps, ncol = 3) 

# for(ii in 1:bs_all_reps) {
#   tryCatch({
#     bs_sample <- sample(1:nrow(owl_stard), nrow(owl_stard), replace = TRUE)
    
#     bs_df <- owl_stard[bs_sample, ]
    
#     # Do Prediction
#     bs_privacy_parms <- compute_prowl_parms(kernel_parms = list(F = Inf, 
#                                                              kappa = 1, 
#                                                              sigma_p = sqrt(d), 
#                                                              diam = diam),
#                                          xi,
#                                          zeta,
#                                          cLower,
#                                          cUpper,
#                                          L,
#                                          cost,
#                                          d,
#                                          n = n,
#                                          beta = NULL,
#                                          Delta = Delta,
#                                          epsilon = epsilon)
    
    
#     bs_priv_est <- PrOWL(df_formula = stard_form, 
#                           data = bs_df, 
#                           cost = cost,
#                           kernel_type = "rbfdot", 
#                           treatment_var = "A", 
#                           reward_var = "Y", 
#                           prob_var = "pA",
#                           lambda = bs_privacy_parms$lambda,
#                           epsilon = epsilon,
#                           d = bs_privacy_parms$set_parms$d,
#                           d_tilde = bs_privacy_parms$d_tilde) 
    
#     # Try to get OWL working well
#     bs_standard_owl_with_cv <- OWL(df_formula = stard_form, 
#                                 data = bs_df, 
#                                 cost = seq(0.01, 1, length.out = 10),
#                                 cv_folds = 10,
#                                 kernel_type = "rbfdot", 
#                                 treatment_var = "A", 
#                                 reward_var = "Y",
#                                 prob_var = "pA")
    
#     bs_preds_priv <- priv_est$predict_treatment(model.matrix(stard_form, bs_df))
#     bs_preds_nonpriv <- standard_owl_with_cv$predict_treatment(model.matrix(stard_form, bs_df))
    
#     all_results[ii,] <- c(estimated_value(bs_preds_priv, bs_df$A, bs_df$Y, bs_df$pA),
#                             estimated_value(bs_preds_nonpriv, bs_df$A, bs_df$Y, bs_df$pA),
#                             sum(bs_preds_priv == bs_preds_nonpriv))
#   }, error=function(e){cat(ii, " - ERROR :",conditionMessage(e), "\n")})
# }

# all_results <- all_results[complete.cases(all_results),]
# round(quantile(all_results[,1], probs = c(0.025, 0.975)), 2)
# round(quantile(all_results[,2], probs = c(0.025, 0.975)), 2)
# round(quantile(all_results[,3], probs = c(0.025, 0.975)), 2)


