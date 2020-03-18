# The following code snippet has company oriented information taken away
# and serve only as a simple construct to demonstrate the idea

library(tidyverse)
library(ompr)
library(ROI)
library(ROI.plugin.glpk)
library(ompr.roi)

setwd("~/CreditRisk/Test Linear Program")

df_fraud_wrangle <- readRDS("df_fraud_wrangle.rds")

glimpse(df_fraud_wrangle)

df_fraud_wrangle <- df_fraud_wrangle %>%
  mutate(
    INTERNET_IND = as.factor(case_when(
      as.character(SLS_DIST_CHNL_TYPE_CD_ORIG) == "N" ~ "Internet",
      TRUE ~ "Non-Internet"
    )),
    REAL_OUTCOME_IND = case_when(
      as.character(REAL_OUTCOME) == "Bad" ~ 1,
      TRUE ~ 0
    )
    ) %>%
  select(CREDIT_APP_ID, FLAG_STATUS, FIQ_SCORE, INTERNET_IND, REAL_OUTCOME_IND)


### first conceptually setup the maximization problem

# define following variables:
#   f_ij is whether the i-th point in the j-th channel is flagged
#   s_j is the score threshold for channel j 

# define the following coefficients/parameters:
#   t_ij is the true binary outcome for the i-th point in the j-th channel
#   s_ij is the score for the i-th point in the j-th channel

# for i: 1 for "positive", 0 for "negative"
# for j: 1 for internet, 2 for non-internet

# MAXIMIZE (from an overall population perspective)
#   true_positive_rate - false_positive_rate

# SUBJECT TO

# A. if person i's score is greater than or equal to threshold score, then flag is 1, 0 otherwise 
#    
#   s_ij >= s_j - M * (1 - f_ij) for all i, j
#   s_ij <= s_j + M * f_ij for all i, j

# B. each flag rate should be less than or equal to 5%
#   (f_11 + f_21 + ... + f_n1) / n <= 0.05  
#   (f_12 + f_22 + ... + f_m2) / m <= 0.05

# C. variable ranges
#   f_ij in {0, 1} for all i, j
#   s_j between 0 and 1000 for all j

# where true_positive_rate = (1 / (t_11 + t_21 + ... + t_n1)) * (f_11 * t_11 + f_21 * t_21 + ... + f_n1 * t_n1) +
#                            (1 / (t_12 + t_22 + ... + t_m2)) * (f_12 * t_12 + f_22 * t_22 + ... + f_m2 * t_m2)
# and false_positive_rate =  (1 / ((1 - t_11) + (1 - t_21) + ... + (1 - t_n1)) * (f_11 * (1 - t_11) + f_21 * (1 - t_21) + ... + f_n1 * (1 - t_n1)) +
#                            (1 / ((1 - t_12) + (1 - t_22) + ... + (1 - t_m2)) * (f_12 * (1 - t_12) + f_22 * (1 - t_22) + ... + f_n2 * (1 - t_n2))


# in ompr...
df_fraud_wrangle_sample <- df_fraud_wrangle %>% sample_n(1000)
c <- 2  # number of channels
n <- df_fraud_wrangle_sample %>% filter(as.character(INTERNET_IND) == 'Internet') %>% count() %>% pull()  # number of data points in Internet
m <- df_fraud_wrangle_sample %>% filter(as.character(INTERNET_IND) != 'Internet') %>% count() %>% pull()  # number of data points in Non-Internet
B <- 1000000

totals_by_channel <- df_fraud_wrangle_sample %>%
  group_by(INTERNET_IND) %>%
  summarize(TOTAL_POSITIVES = sum(REAL_OUTCOME_IND),
            TOTAL_NEGATIVES = sum(1 - REAL_OUTCOME_IND))
  
df_fraud_wrangle_agg <- df_fraud_wrangle_sample %>%
  group_by(INTERNET_IND) %>%
  mutate(TOTAL_POSITIVES = sum(REAL_OUTCOME_IND),
         TOTAL_NEGATIVES = sum(1 - REAL_OUTCOME_IND),
         TP_WEIGHT = (1 / TOTAL_POSITIVES) * REAL_OUTCOME_IND,
         FP_WEIGHT = (1 / TOTAL_NEGATIVES) * (1 - REAL_OUTCOME_IND)) %>%
  select(INTERNET_IND, FIQ_SCORE, REAL_OUTCOME_IND, TP_WEIGHT, FP_WEIGHT) %>%
  arrange(INTERNET_IND)

tp_weight <- df_fraud_wrangle_agg$TP_WEIGHT  # weights of true positives for objective function
fp_weight <- df_fraud_wrangle_agg$FP_WEIGHT  # weights of false positives for objective function
score <- df_fraud_wrangle_agg$FIQ_SCORE


mip_model <- MIPModel() %>%  # 0. start with an "empty" mixed-integer program
  
  # 1. add the flag variables
  add_variable(f[i], i = 1:(n + m), type = "binary") %>%  
  
  # 2. add the score threshold variables which are integers in [0, 1000]
  add_variable(s[j], j = 1:c, type = "integer", lb = 0, ub = 1000) %>%  
  
  # 3. add the objective function
  set_objective(sum_expr((tp_weight[i] - fp_weight[i]) * f[i], i = 1:(n + m))) %>%  
  
  # 4. finally add all the constraints
  
  #  4.a. conditional constraint on score for Internet (j = 1)
  add_constraint(score[i] >= s[1] - B * (1 - f[i]), i = 1:n) %>%
  add_constraint(score[i] <= s[1] + B * f[i], i = 1:n) %>%
  
  #  4.b. conditional constraint on score for Non-Internet (j = 2)
  add_constraint(score[i] >= s[2] - B * (1 - f[i]), i = (n + 1):(n + m)) %>%
  add_constraint(score[i] <= s[2] + B * f[i], i = (n + 1):(n + m)) %>%
  
  # 4.c. flag constraints in each group
  add_constraint(sum_expr(f[i], i = 1:n) / n <= 0.05) %>%
  add_constraint(sum_expr(f[i], i = (n + 1):(n + m)) / m <= 0.05) %>%
  
  # finally, solve the MIP!
  solve_model(with_ROI("glpk", verbose = TRUE))

get_solution(mip_model, s[j])
#   variable j value
# 1        s 1   908
# 2        s 2   592

# now wrap into a function and run it many times
milp_solver_func <- function(df) {
  df_fraud_wrangle_sample <- df %>% sample_n(1000)
  c <- 2  # number of channels
  n <- df_fraud_wrangle_sample %>% filter(as.character(INTERNET_IND) == 'Internet') %>% count() %>% pull()  # number of data points in Internet
  m <- df_fraud_wrangle_sample %>% filter(as.character(INTERNET_IND) != 'Internet') %>% count() %>% pull()  # number of data points in Non-Internet
  B <- 1000000
  
  totals_by_channel <- df_fraud_wrangle_sample %>%
    group_by(INTERNET_IND) %>%
    summarize(TOTAL_POSITIVES = sum(REAL_OUTCOME_IND),
              TOTAL_NEGATIVES = sum(1 - REAL_OUTCOME_IND))
  
  df_fraud_wrangle_agg <- df_fraud_wrangle_sample %>%
    group_by(INTERNET_IND) %>%
    mutate(TOTAL_POSITIVES = sum(REAL_OUTCOME_IND),
           TOTAL_NEGATIVES = sum(1 - REAL_OUTCOME_IND),
           TP_WEIGHT = (1 / TOTAL_POSITIVES) * REAL_OUTCOME_IND,
           FP_WEIGHT = (1 / TOTAL_NEGATIVES) * (1 - REAL_OUTCOME_IND)) %>%
    select(INTERNET_IND, FIQ_SCORE, REAL_OUTCOME_IND, TP_WEIGHT, FP_WEIGHT) %>%
    arrange(INTERNET_IND)
  
  tp_weight <- df_fraud_wrangle_agg$TP_WEIGHT  # weights of true positives for objective function
  fp_weight <- df_fraud_wrangle_agg$FP_WEIGHT  # weights of false positives for objective function
  score <- df_fraud_wrangle_agg$FIQ_SCORE
  
  
  mip_model <- MIPModel() %>%  # 0. start with an "empty" mixed-integer program
    
    # 1. add the flag variables
    add_variable(f[i], i = 1:(n + m), type = "binary") %>%  
    
    # 2. add the score threshold variables which are integers in [0, 1000]
    add_variable(s[j], j = 1:c, type = "integer", lb = 0, ub = 1000) %>%  
    
    # 3. add the objective function
    set_objective(sum_expr((tp_weight[i] - fp_weight[i]) * f[i], i = 1:(n + m))) %>%  
    
    # 4. finally add all the constraints
    
    #  4.a. conditional constraint on score for Internet (j = 1)
    add_constraint(score[i] >= s[1] - B * (1 - f[i]), i = 1:n) %>%
    add_constraint(score[i] <= s[1] + B * f[i], i = 1:n) %>%
    
    #  4.b. conditional constraint on score for Non-Internet (j = 2)
    add_constraint(score[i] >= s[2] - B * (1 - f[i]), i = (n + 1):(n + m)) %>%
    add_constraint(score[i] <= s[2] + B * f[i], i = (n + 1):(n + m)) %>%
    
    # 4.c. flag constraints in each group
    add_constraint(sum_expr(f[i], i = 1:n) / n <= 0.05) %>%
    add_constraint(sum_expr(f[i], i = (n + 1):(n + m)) / m <= 0.05) %>%
    
    # finally, solve the MIP!
    solve_model(with_ROI("glpk", verbose = TRUE))
  
  get_solution(mip_model, s[j]) %>% as_tibble()
  
  #c("s1" = results$value[1], "s2" = results$value[2])
}

results_1 <- milp_solver_func(df = df_fraud_wrangle)

results_10 <- map_df(1:10, ~milp_solver_func(df = df_fraud_wrangle))
results_10 %>% filter(j == 1) %>% pull(value) %>% summary()
results_10 %>% filter(j == 2) %>% pull(value) %>% summary()

results_100 <- map_df(1:100, ~milp_solver_func(df = df_fraud_wrangle))
results_100 %>% filter(j == 1) %>% pull(value) %>% summary()
results_100 %>% filter(j == 2) %>% pull(value) %>% summary()
