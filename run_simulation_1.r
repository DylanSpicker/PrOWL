
# Test DF 
source("src/prowl.R")
source("src/simulation1.R")
# library(future.apply)
# library("doFuture")
# library("doRNG")

# Define a helper function to simplify the complexity of the simulations
extract_results <- function(results, df_formula, val_df, jj, compare_to = NULL, compare_at = NULL) {
    fitted_value <- results$value

    # Return the estimated value
    val_predicted_vals <- results$predict_treatment(model.matrix(df_formula, val_df))
    val_value <- estimated_value(val_predicted_vals, val_df$A, unlist(val_df[paste0("R",jj)]), val_df$a_p)
    val_acc <- mean(factor(val_predicted_vals) == unlist(val_df[paste0("Astar",jj)]))
    max_value <- estimated_value(unlist(val_df[paste0("Astar",jj)]), val_df$A, unlist(val_df[paste0("R",jj)]), val_df$a_p)

    agree <- NA
    delta.agree <- NA
    if(!is.null(compare_to)){
        compare_vals <- compare_to$predict_classifier(model.matrix(df_formula, val_df))
        compare_preds <- sign(compare_vals)

        d_filter <- which(compare_vals >= compare_at)

        agree <- mean(val_predicted_vals == compare_preds)
        
        if(length(d_filter) == 0) {
            delta.agree <- 1
        } else {
            delta.agree <- mean(val_predicted_vals[d_filter] == compare_preds[d_filter])
        }
    }

    c(fitted_value,         # Value on the training set
      val_acc,              # Accuracy on the validation
      val_value,            # Value on the validation
      max_value,            # Maximum possible value
      sum(val_predicted_vals == 1), # Number of '1' diagnoses
      sum(val_predicted_vals == -1), # Number of '-1' diagnoses
      agree,                # Agreement with 'compare_to'
      delta.agree           # Agreement with 'compare_to', filtered by those >= 'compare_at'
      )
}

# Define Simulation Parameters
nruns <- 1000
which_scenarios <- c(1, 4)
scenario_labels <- c("_A", "_B")

# SVM Parameters
df_formula <- A ~ -1 + X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10
ns <- c(100, 200, 400, 1000)
epsilons <- c(0.5, 1, 5, 10)
cost <- 1/100
failed_eps_setting <- 1/100 # A parameter to help with stability of estimators;
                            # The larger the value of this, the less sensitive the numeric algorithm is.
                            # 1/100 should suffice for all scenarios consdiered
        
fit_kernel <- TRUE      # Should the kernel version of OWL be run?
fit_linear <- TRUE      # Should the lienar version of OWL be run?

# Set up results list
privacy_parms <- compute_sim1_privacy_parms(ns, epsilons)
all_results <- vector(length(ns), mode = "list")

# Run this in parallel ?
# No. Strangely: setting more than 1 worker dramatically slows down the performance,
#            however, multisession with 1 worker is better than sequential.
#            unsure of why this is happening on my machine, but could be adjusted otherwise
# registerDoFuture()
# plan(multisession, workers=1L)

results_per <- 8
num_scenarios <- fit_kernel + fit_linear + length(epsilons) * 2

start_time <- Sys.time()
zz <- 1

while (zz <= length(ns)) {
    
    set.seed(31415 * zz)

    sc_n <- ns[zz]
    cat(paste0("Working on ", sc_n, "\n"))

    # mat_results <- matrix(nrow = nruns, ncol = length(which_scenarios) * 4 * 2)
    # mat_results_LOWL <- matrix(nrow = nruns, ncol = length(which_scenarios) * 4)

    cat(paste0("Starting at ", Sys.time(), ". It has been ", difftime(Sys.time(), start_time, units="mins"), " minutes since starting.\n"))

    # mat_results <- foreach(ii = 1:nruns, .combine=rbind) %dorng% {
    mat_results <- matrix(nrow = nruns, ncol = length(which_scenarios) * results_per * num_scenarios)
    for(ii in 1:nruns){
        cat(paste0("\r", ii, " / ", nruns))
        # Compute privacy parameters
        my_df <- generate_single_run_data(n = sc_n)
        val_df <- generate_single_run_data(n = 10000)

        save(my_df, val_df, file = paste0("src/backups/", sc_n,"/", stringr::str_pad(ii, 4, pad="0"), "_df.Rda"))

        run_results <- rep(0, length(which_scenarios) * results_per * num_scenarios)
        
        for(jj in which_scenarios){
            jj_idx <- which(which_scenarios == jj)
            idx_offset <- (jj_idx-1) * results_per * num_scenarios

            # Run the kernel OWL estimator [non-private]
            if(fit_kernel) {
                kernel_owl_results <- tryCatch(
                    OWL(df_formula = df_formula, 
                            data = my_df, 
                            cost = cost,
                            cv_folds = 1,
                            kernel_type = "rbfdot", 
                            treatment_var = "A", 
                            reward_var = paste0("R",jj), 
                            prob_var = "a_p"),
                    error=function(e) e
                )
                if(!inherits(kernel_owl_results, "error")) {
                    run_results[(idx_offset + 1):(idx_offset + results_per)] <- extract_results(results = kernel_owl_results, df_formula, val_df, jj)
                } else {
                    # If the original algorithm fails; try again overriding epsilon.
                    kernel_owl_results <- tryCatch(
                        OWL(df_formula = df_formula, 
                                data = my_df, 
                                cost = cost,
                                cv_folds = 1,
                                kernel_type = "rbfdot", 
                                treatment_var = "A", 
                                reward_var = paste0("R",jj), 
                                prob_var = "a_p", 
                                eps = failed_eps_setting),
                        error=function(e) e
                    )

                    if(!inherits(kernel_owl_results, "error")) {
                        run_results[(idx_offset + 1):(idx_offset + results_per)] <- extract_results(results = kernel_owl_results, df_formula, val_df, jj)
                    } else {
                        run_results[(idx_offset + 1):(idx_offset + results_per)] <- rep(NA, results_per)
                    } 
                }
            }
            

            # Run the OWL estimator [non-private]
            if(fit_linear){
                linear_owl_results <- tryCatch(
                    OWL(df_formula = df_formula, 
                            data = my_df, 
                            cost = cost,
                            cv_folds = 1,
                            kernel_type = "vanilladot", 
                            treatment_var = "A", 
                            reward_var = paste0("R",jj), 
                            prob_var = "a_p",
                            min_weight = 0),
                    error=function(e) e
                )
                if(!inherits(linear_owl_results, "error")) {
                    run_results[(idx_offset + results_per + 1):(idx_offset + 2*results_per)] <- extract_results(results = linear_owl_results, df_formula, val_df, jj)
                } else {
                    # If the original algorithm fails; try again overriding epsilon.
                    linear_owl_results <- tryCatch(
                        OWL(df_formula = df_formula, 
                                data = my_df, 
                                cost = cost,
                                cv_folds = 1,
                                kernel_type = "vanilladot", 
                                treatment_var = "A", 
                                reward_var = paste0("R",jj), 
                                prob_var = "a_p",
                                min_weight = 0,
                                eps = failed_eps_setting),
                        error=function(e) e
                    )

                    if(!inherits(linear_owl_results, "error")) {
                        run_results[(idx_offset + results_per + 1):(idx_offset + 2*results_per)] <- extract_results(results = linear_owl_results, df_formula, val_df, jj)
                    } else {
                        run_results[(idx_offset + results_per + 1):(idx_offset + 2*results_per)] <- rep(NA, results_per)
                    }
                }
            }
            
            # Run the kernel prowl over values of epsilon
            for(tt in 1:length(epsilons)){
                my_p <- privacy_parms[[as.character(sc_n)]][[as.character(epsilons[tt])]][[paste0("infinite", scenario_labels[jj_idx])]]
                possibleError <- tryCatch(
                    PrOWL(df_formula = df_formula, 
                            data = my_df, 
                            cost = cost,
                            kernel_type = "rbfdot", 
                            treatment_var = "A", 
                            reward_var = paste0("R",jj), 
                            prob_var = "a_p",
                            lambda = my_p$lambda,
                            epsilon = epsilons[tt],
                            d = my_p$set_parms$d,
                            d_tilde = my_p$d_tilde),
                    error=function(e) e
                )

                if(!inherits(possibleError, "error")) {
                    run_results[(idx_offset + (tt-1)*results_per + 2*results_per + 1):(idx_offset + (tt-1)*results_per + 3*results_per)] <- extract_results(results = possibleError, df_formula, val_df, jj, kernel_owl_results, my_p$Delta)
                } else {
                    possibleError <- tryCatch(
                        PrOWL(df_formula = df_formula, 
                                data = my_df, 
                                cost = cost,
                                kernel_type = "rbfdot", 
                                treatment_var = "A", 
                                reward_var = paste0("R",jj), 
                                prob_var = "a_p",
                                lambda = my_p$lambda,
                                epsilon = epsilons[tt],
                                d = my_p$set_parms$d,
                                d_tilde = my_p$d_tilde,
                                eps = failed_eps_setting),
                        error=function(e) e
                    )
                    if(!inherits(possibleError, "error")) {
                        run_results[(idx_offset + (tt-1)*results_per + 2*results_per + 1):(idx_offset + (tt-1)*results_per + 3*results_per)] <- extract_results(results = possibleError, df_formula, val_df, jj, kernel_owl_results, my_p$Delta)
                    } else {
                        run_results[(idx_offset + (tt-1)*results_per + 2*results_per + 1):(idx_offset + (tt-1)*results_per + 3*results_per)] <- rep(NA, results_per)
                    }
                }

            }

            # Run the prowl
            for(tt in 1:length(epsilons)){
                my_p <- privacy_parms[[as.character(sc_n)]][[as.character(epsilons[tt])]][[paste0("finite", scenario_labels[jj_idx])]]
                possibleError <- tryCatch(
                    PrOWL(df_formula = df_formula, 
                            data = my_df, 
                            cost = cost,
                            kernel_type = "vanilladot", 
                            treatment_var = "A", 
                            reward_var = paste0("R",jj), 
                            prob_var = "a_p",
                            lambda = my_p$lambda,
                            epsilon = epsilons[tt]),
                    error=function(e) e
                )

                if(!inherits(possibleError, "error")) {
                    run_results[(idx_offset + (tt-1)*results_per + 6*results_per + 1):(idx_offset + (tt-1)*results_per + 7*results_per)] <- extract_results(results = possibleError, df_formula, val_df, jj, linear_owl_results, my_p$Delta)
                } else {
                    possibleError <- tryCatch(
                        PrOWL(df_formula = df_formula, 
                                data = my_df, 
                                cost = cost,
                                kernel_type = "vanilladot", 
                                treatment_var = "A", 
                                reward_var = paste0("R",jj), 
                                prob_var = "a_p",
                                lambda = my_p$lambda,
                                epsilon = epsilons[tt],
                                eps = failed_eps_setting),
                        error=function(e) e
                    )
                    if(!inherits(possibleError, "error")) {
                        run_results[(idx_offset + (tt-1)*results_per + 6*results_per + 1):(idx_offset + (tt-1)*results_per + 7*results_per)] <- extract_results(results = possibleError, df_formula, val_df, jj, linear_owl_results, my_p$Delta)
                    } else {
                        run_results[(idx_offset + (tt-1)*results_per + 6*results_per + 1):(idx_offset + (tt-1)*results_per + 7*results_per)] <- rep(NA, results_per)
                    }
                }

            }
        }
        mat_results[ii, ] <- run_results 
    }

    cat(paste0("\nCompleted in ", difftime(Sys.time(), start_time, units="mins"), " minutes. Saving.\n"))
    all_results[[zz]] <- mat_results
    zz <- zz + 1
    save(all_results, file = paste0("src/sim_mat_results_FULL_", stringr::str_replace_all(Sys.Date(), "-", "."), ".RDa"))
}


# Work through all of the analysis
library(tidyverse)
library(xtable)

# Write out the analysis for a single result df 

for(ii in 1:length(all_results)) {
    single_result_df <- all_results[[ii]]
    
    # Write the column names
    # Order goes: kernel_owl_A, linear_owl_A, kernel_prowl_A [x4], linear_prowl_A [x4]
    #             kernel_owl_B, linear_owl_B, kernel_prowl_B [x4], linear_prowl_B [x4]
    # Inside each order goes:
    #   fitted_value, accuracy, max_value, treatment_A, treatment_B, agreement, delta_agreement

    estring <- stringr::str_replace(epsilons, "\\.", "")

    names_vector <- apply(expand.grid(
        c("fvalue", "accuracy", "vvalue", "mvalue", "nA", "nB", "agree", "dagree"),
        c("kernel_owl", "linear_owl", paste0("kernel_prowl", estring), paste0("linear_prowl", estring)),
        c("A", "B")), 1, paste, collapse="_")

    colnames(single_result_df) <- names_vector

    single_result_df %>% as.data.frame() %>% select(starts_with("nA")) %>% colMeans()
    single_result_df %>% as.data.frame() %>% select(starts_with("fvalue")) %>% colMeans()
    single_result_df %>% as.data.frame() %>% select(starts_with("accuracy")) %>% colMeans()
    single_result_df %>% as.data.frame() %>% select(starts_with("vvalue")) %>% colMeans()
    single_result_df %>% as.data.frame() %>% select(starts_with("agree")) %>% colMeans()
    single_result_df %>% as.data.frame() %>% select(starts_with("dagree")) %>% colMeans()

    ## Build out the table results
    # Each row represents a different statistic; each column will be a method
    ordered_df <- as.data.frame(single_result_df) %>% 
                    mutate(run = row_number()) %>%
                    pivot_longer(cols = -run, names_sep = "_", names_to = c("stat", "kernel", "technique", "scenario")) %>%
                    pivot_wider(id_cols = c("stat", "scenario", "run"), names_from = c("kernel", "technique"), names_sep = "_") %>% 
                    group_by(stat, scenario) %>% 
                    summarise(across(-run, ~ mean(.x, na.rm = TRUE))) %>% 
                    select(
                        "stat",
                        "scenario",
                        "linear_owl",
                        "kernel_owl",
                        "linear_prowl05",
                        "linear_prowl1",
                        "linear_prowl5",
                        "linear_prowl10",
                        "kernel_prowl05",
                        "kernel_prowl1",
                        "kernel_prowl5",
                        "kernel_prowl10"
                    ) %>% 
                    arrange(scenario, stat)

    # Convert vvalue into regrets
    ordered_df[which(ordered_df$stat == "vvalue"), -c(1,2)] <- ordered_df[which(ordered_df$stat == "mvalue"), -c(1,2)] - ordered_df[which(ordered_df$stat == "vvalue"), -c(1,2)]

    # Convert agreement into nA for OWL
    ordered_df[which(ordered_df$stat == "agree"), c(3, 4)] <- ordered_df[which(ordered_df$stat == "nA"), c(3, 4)]

    # Select only the right statistics
    ordered_df[which(ordered_df$stat %in% c("agree", "accuracy", "dagree", "vvalue")), ] %>% xtable(x = ., type = "latex") %>% print(include.rownames = FALSE, NA.string = "--")

}
