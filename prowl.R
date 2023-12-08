# source("src/weighted_ksvm.R")
library(personalized)


# Split function for CV
CV_chunk <- function(x,n) split(x, cut(seq_along(x), n, labels = FALSE)) 

# Step 1: Sample rho vectors from the Fourier transform;
#         Requires: d, Kernel, d-tilde
#         Returns: List of rho vectors and varphi
approximate_feature_mapping <- function(d, d_tilde, sample_fn) {

    # Rho will be a matrix, each row representing a d-vector
    # with d-tilde columns.
    rho <- sample_fn(d, d_tilde) # Sample from the correct distribution

    varphi <- function(x) {
        inner_prod <- t(rho) %*% x

        cos_trans <- cos(inner_prod)
        sin_trans <- sin(inner_prod)

        (1/sqrt(d_tilde)) * c(rbind(cos_trans, sin_trans))
    }

    list(
        "mapping" = varphi,
        "rho" = rho
    )   
}

# Step 2: Define the kernel from the feature mapping
define_kernel <- function(varphi) {
    function(x, y) { t(varphi(x)) %*% varphi(y) }
}

# Step 3: Compute the individual weights
# Step 3A: Find the propensity scores
ridge_logistic_regression <- function(formula, data) {
    # This function call is an interface to glmnet

    pred_mat <- model.matrix(formula, data)
    outcome_mat <- data[all.vars(formula)[1], ]

    glmnet::glmnet(
        x = pred_mat,           # Matrix of predictors
        y = otucome_mat,        # Outcome variable (treatment)
        family = "binomial",    # Family is binomial for logistic regression
        alpha = 0               # Alpha = 0 forces ridge penalty
    )
}

# CBSR.R > ridge_logistic_cbsr() provides the capacity to calculate weights using
#                                the CBSR setup in place of logistic regression

# Step 3B: Determine the actual weights
individual_weights <- function(formula, data, rewards, treatments) {
    # OWL weights are R / propensity score
    # As a result we simply need to take rewards over propensity scores 
    propensity_score_model <- ridge_logistic_regression(formula, data)
    probabilities <- predict(propensity_score_model)

    rewards / ((treatments == 1)*probabilities + (treatments == -1)*(1 - probabilities))
}

# Step 6: Privatize the Vector
privatize_vector <- function(w, lambda, epsilon) { w + VGAM::rlaplace(length(w), 0, lambda/epsilon) }


# Convert the predicted values, treatments, rewards, and probabilities into the estimated value 
# function. This is for use in cross validation
estimated_value <- function(predicted_vals, treatments, rewards, probs) {
    
    m_probs <- probs * (treatments == 1) + (1 - probs) * (treatments == -1)
  
    numerator <- as.numeric(mean((predicted_vals == treatments)*rewards/m_probs))
    denominator <- as.numeric(mean((predicted_vals == treatments)/m_probs))
    numerator/denominator
}


# OWL: Outcome-Weighted Learning
OWL <- function(df_formula,
                data,
                kernel_type = "rbfdot",
                cv_folds = 5,
                cost = 1/100,
                kpar = "automatic",
                treatment_var,
                reward_var,
                prob_var,
                min_weight = 0,
                ...) {
                    
    # Specify the kernel type
    if(startsWith(kernel_type, "p")) {
        kernel_type <- "polydot"
    } else if(startsWith(kernel_type, "v")) {
        kernel_type <- "vanilladot"
    } else { 
        # Infinite Dimensional Kernels
        kernel_type <- "rbfdot"
        if(startsWith(kernel_type, "l")) {
            kernel_type <- "laplacedot"
        } else if(! startsWith(kernel_type, "r")) { 
            warning("The provided kernel did not match any of `rbfdot`, `polydot`, `vanilladot`, or `laplacedot`. The results use `rbfdot`.") 
        }
    }

    # Extract 
    treatments <- as.numeric(unlist(data[treatment_var]))
    rewards <- as.numeric(unlist(data[reward_var]))
    probs <- as.numeric(unlist(data[prob_var]))
    weight <- rewards / ((treatments == 1)*probs + (treatments == -1)*(1 - probs))
    
    # Model Matrix of Predictors
    model_x <- model.matrix(df_formula, data)
    
    # Try to get the weighted KSVM fit
    w_svm <- weighted.ksvm(y = treatments, 
                           x = model_x, 
                           weights = pmax(min_weight, weight),
                           C = cost,
                           kernel = kernel_type,
                           kpar = kpar,
                           nfolds = cv_folds,
                           ...)

    # Return the estimated value
    predicted_vals <- predict(w_svm, newx = model_x, type = 'class')
    ev <- estimated_value(predicted_vals, factor(treatments), rewards, probs)

    # Build Object to Return
    list(
        method = "OWL",
        decision_function_model = df_formula,
        fitted_svm = w_svm,
        value = ev,
        predict_treatment = function(ind_data) {
            predict(w_svm, newx = ind_data, type = 'class')
        },
        predict_classifier = function(ind_data) {
            predict(w_svm, newx = ind_data, type = 'linear.predictor')
        }
    )
}

# Define Cross Validation Function
self_CV.OWL <- function(df_formula,
                        data,
                        kernel_type = "rbfdot",
                        grid_size = 50,
                        cv_folds = 5,
                        c_max = 10,
                        treatment_var,
                        reward_var,
                        prob_var) {

    cost_values <- seq(1/ncol(data), c_max, length.out = grid_size)
    cv_table <- expand.grid(cost_values)

    colnames(cv_table) <- c("cost")
    
    # Perform Cross Validation
    cv_vals <- CV_chunk(sample(1:nrow(data), nrow(data), FALSE), cv_folds)

    # Go through the different CV scenarios
    cv_results_list <- future_lapply(1:nrow(cv_table), function(ii) {
        
        parms <- cv_table[ii, ]

        cost <- parms[["cost"]]
        # kpar <- as.list(parms[names(parms)[-which(names(parms) == "cost")]])
        # names(kpar) <- names(parms)[-which(names(parms) == "cost")]
        kpar <- "automatic"
        
        results <- rep(0, cv_folds)
        kernels <-  vector(mode = "list", length = cv_folds)
        
        for(foldidx in 1:cv_folds){
            
            # Split data
            test_data <- data[cv_vals[[foldidx]], ]
            train_data <- data[-cv_vals[[foldidx]], ]
            
        
            # Model Matrix of Predictors
            # model_x <- model.matrix(df_formula, train_data)

            # Try to get the weighted KSVM fit
            cur_OWL <- OWL(
                df_formula,
                train_data,
                kernel_type = kernel_type,
                grid_size = grid_Size,
                cv_folds = cv_folds,
                cost = cost,
                kpar = kpar,
                treatment_var,
                reward_var,
                prob_var)

            # Return the estimated value
            predicted_vals <- cur_OWL$predict_treatment(model.matrix(df_formula, test_data))


            results[foldidx] <- estimated_value(predicted_vals, 
                                                factor(as.numeric(unlist(data[treatment_var]))), 
                                                as.numeric(unlist(data[reward_var])),
                                                as.numeric(unlist(data[prob_var])))
            kernels[[foldidx]] <- cur_OWL$fitted_svm$kernel
            
        }

        list("results" = c(results, mean(results)), "kernels" = kernels)

    })

    
    cv_results <- do.call(rbind, lapply(cv_results_list, function(ii){ ii$results }))
    best_idx <- which.max(cv_results[nrow(cv_results),])
    best_parms <- cv_table[best_idx,]

    if(is.list(best_parms)) {
        best_cost <- best_parms[["cost"]]
        best_kpar <- as.list(best_parms[names(best_parms)[-which(names(best_parms) == "cost")]])
        names(best_kpar) <- names(best_parms)[-which(names(best_parms) == "cost")]
    } else {
        best_cost <- best_parms
        best_kpar <- list()
    }

    cv_folds <- 1
    cv_result_table <- list(parms = cv_table, results = cv_results)

}


# Compute privacy noise
compute_prowl_parms <- function(
    kernel_parms,   # List including F; diam; sigma_p; kappa
    xi,             # Bounds on the (modified) rewards
    zeta,           # Bounds on the sensitivity of propensity scores
    cLower,            # Lower bound truncation of propensity score
    cUpper,            # Upper bound truncation of propensity score
    L,              # Lipschitz parameter for the loss function
    cost,           # Cost parameter
    d = NULL,       # Dimension
    n = NULL,       # Sample size
    beta = NULL,    # Agreement Probability
    Delta = NULL,   # Indifference parameter
    epsilon = NULL  # Privacy parameter
) {

    # Compute lower bound probability.
    c_L <- min(cLower, 1-cUpper)

    # Kernel Parameters
    if(! "F" %in% names(kernel_parms)) stop("`F` must be provided in the kernel parameters.")
    F <- kernel_parms$F

    # Define helper quantities for function progress
    nulls <- c(n = is.null(n), beta = is.null(beta), Delta = is.null(Delta), epsilon = is.null(epsilon))
    return_range <- FALSE
    
    # Run check based on the nulls
    if(sum(nulls) > 1) {
        stop("Only one of `n`, `beta`, `Delta`, or `epsilon` can be null.")
    } else if(sum(nulls) == 0) { 
        warning("The full set of parameters was provided. Bounds are provided based on the given parameters, and then lambda is solved for to minimize `epsilon`.")
        solve_for <- "epsilon"
        
        # Compute based on given
        lower.lambda <- 2*L*xi*kappa*sqrt(F)/(epsilon*n*c_L)*((n-1)*zeta/c_L + 2)
        upper.lambda <- ifelse(is.infinite(F), -Delta/(sqrt(2)*log(2*(1-sqrt(beta)))), -Delta/(sqrt(2)*kappa*log(2*(1-beta))))
        return_range <- TRUE 
    } else {
        solve_for <- names(nulls)[which(nulls)]
    }
    

    # Check on infinite kernel parameters
    if(is.infinite(F)) {
        # Check on the infinite kernel parameters
        if(! "diam" %in% names(kernel_parms)) stop("Since `F` is infinite, `diam` must be provided in the kernel parameters.")
        diam <- kernel_parms$diam

        if(! "sigma_p" %in% names(kernel_parms)) stop("Since `F` is infinite, `sigma_p` must be provided in the kernel parameters.")
        sigma_p <- kernel_parms$sigma_p

        if(is.null(d)) stop("The dimension of the database (`d`) is required for infinite kernels.")

        # Set `kappa`
        kappa <- 1

        # Solving for 'n' or 'epsilon' does not require any further dependence 
        # between d_tilde and the parameters. 
        # In the event that `beta` is solved for, beta needs to be used in both
        # the computation of d_tilde and vice-versa
        if(solve_for == "n") {
            theta.Delta <- min(1, Delta^2/cost*(1+sqrt(6)))^2
            d_tilde <- ceiling(4*(d+2)/theta.Delta*log(2^8*(sigma_p*diam)^2/((1-sqrt(beta))*theta.Delta)))
            F <- 2*d_tilde

            LBV <- log(2*(1-sqrt(beta)))
            A <- 2*L*xi*kappa*sqrt(F)/(epsilon*c_L)
            B <- zeta/c_L
            C <- -Delta/(sqrt(2)*kappa*LBV)

            n <- 2*A*(1-B)/(C-A*B)
        } else if (solve_for == "beta") {
            # Function to solve for the joint zeros
            theta.Delta <- min(1, Delta^2/cost*(1+sqrt(6)))^2

            C1 <- (-Delta*epsilon*n*c_L/(2*sqrt(2)*kappa^2*L*((n-1)*xi/c_L + 2)))^2/2
            C2 <- 4*(d+2)/theta.Delta
            C3 <- 2^8*(sigma_p*diam)^2/theta.Delta

            cubic_coefs <- c(-1, log(C3)-2*log(2), (2*log(2)*log(C3) - log(2)^2), log(2)^2*log(C3)-C1/C2)
            possible_beta <- (1-exp(RConics::cubic(cubic_coefs))/2)^2
            
            beta <- tryCatch({
                max(possible_beta[which(possible_beta > 0.5 & possible_beta < 1)])
            }, warning = function(cond){
                warning("There was no `beta` value which could be established based on the provided parameters. As a result, `beta` is set to 0.501, and a new value for `epsilon` will be solved for.")
                0.501
            })

            d_tilde <- ceiling(4*(d+2)/theta.Delta*log(2^8*(sigma_p*diam)^2/((1-sqrt(beta))*theta.Delta)))
            F <- 2*d_tilde

            if(beta == 0.501) {
                
                LBV <- log(2*(1-sqrt(beta)))
                A <- 2*L*xi*kappa*sqrt(F)/(n*c_L)*((n-1)*zeta/c_L + 2)
                B <- -Delta/(sqrt(2)*kappa*LBV)

                epsilon <- A/B
            }

        } else if (solve_for == "Delta") {
            # Two step process
            # First:  Assume that theta.Delta = 1. Solve for d_tilde and Delta, and then check if Delta corresponds to Theta = 1
            # Second: If not, the solution corresponds to Theta = Delta^4/C^2(1+sqrt(6))^2; use lamW to solve.
            LBV <- log(2*(1-sqrt(beta)))

            # Set theta=1 and run through calculation
            theta.Delta <- 1
            d_tilde <- ceiling(4*(d+2)/theta.Delta*log(2^8*(sigma_p*diam)^2/((1-sqrt(beta))*theta.Delta)))
            F <- 2*d_tilde

            A <- 2*L*xi*kappa*sqrt(F)/(epsilon*n*c_L)*((n-1)*zeta/c_L + 2)
            B <- -sqrt(2)*kappa*LBV

            Delta <- A*B

            if (Delta < sqrt(cost*(1+sqrt(6)))) {
                # Theta will be less than 1, which will further inflate d_tilde. 
                Omega1 <- 4*(d+2)*cost^2*(1+sqrt(6))^2
                Omega2 <- 2^8*sigma_p^2*diam^2*cost^2*(1+sqrt(6))^2/(1-sqrt(beta))
                Omega3 <- sqrt(2)*kappa*LBV*2*L*xi*kappa/(epsilon*n*c_L)*((n-1)*zeta/c_L + 2)

                const1 <- Omega1 / (2 * Omega3^4)
                const2 <- Omega1/(4*Omega3 ^ 4) * (log(Omega2) - 4*log(Omega3) - 2*log(2))

                cubed <- const1 / 3 * lamW::lambertW0(3 * exp( 3 * const2 / const1) / const1)

                if(cubed < 0) stop("The provided parameters cannot be resolved. These imply a negative size for `d_tilde`.")

                d_tilde <- ceiling(cubed^(1/3))
                F <- 2*d_tilde

                A <- 2*L*xi*kappa*sqrt(F)/(epsilon*n*c_L)*((n-1)*zeta/c_L + 2)
                B <- -sqrt(2)*kappa*LBV

                Delta <- A*B
            }

        } else {
            theta.Delta <- min(1, Delta^2/cost*(1+sqrt(6)))^2
            d_tilde <- ceiling(4*(d+2)/theta.Delta*log(2^8*(sigma_p*diam)^2/((1-sqrt(beta))*theta.Delta)))
            F <- 2*d_tilde

            LBV <- log(2*(1-sqrt(beta)))
            A <- 2*L*xi*kappa*sqrt(F)/(n*c_L)*((n-1)*zeta/c_L + 2)
            B <- -Delta/(sqrt(2)*kappa*LBV)

            epsilon <- A/B
        }

    } else {
        if(! "kappa" %in% names(kernel_parms)) stop("Since `F` is finite, `kappa` must be provided in the kernel parameters.")
        kappa <- kernel_parms$kappa

        # Solving things is easy
        if(solve_for == "n") {
            LBV <- log(2*(1-beta))
            A <- 2*L*xi*kappa*sqrt(F)/(epsilon*c_L)
            B <- zeta/c_L
            C <- -Delta/(sqrt(2)*kappa*LBV)

            n <- 2*A*(1-B)/(C-A*B)
        } else if (solve_for == "beta") {
            A <- 2*L*xi*kappa*sqrt(F)/(epsilon*n*c_L)*((n-1)*zeta/c_L + 2)
            B <- -Delta/(sqrt(2)*kappa)

            beta <- 1 - exp(B/A)/2
        } else if (solve_for == "Delta") {
            LBV <- log(2*(1-beta))
            A <- 2*L*xi*kappa*sqrt(F)/(epsilon*n*c_L)*((n-1)*zeta/c_L + 2)
            B <- -sqrt(2)*kappa*LBV

            Delta <- A*B
        } else {
            LBV <- log(2*(1-beta))
            A <- 2*L*xi*kappa*sqrt(F)/(n*c_L)*((n-1)*zeta/c_L + 2)
            B <- -Delta/(sqrt(2)*kappa*LBV)

            epsilon <- A/B
        }

        d_tilde <- NULL

    }

    if(return_range) {
        range <- c(lower.lambda, upper.lambda)
    } else {
        range <- c(2*L*xi*kappa*sqrt(F)/(epsilon*n*c_L)*((n-1)*zeta/c_L + 2), -Delta/(sqrt(2)*kappa*log(2*(1-beta))))
    }

    # kernel_parms,   # List including F; diam; sigma_p; kappa
    #     xi,             # Bounds on the (modified) rewards
    #     zeta,           # Bounds on the sensitivity of propensity scores
    #     c_L,            # Lower bound truncation of propensity score
    #     c_H,            # Upper bound truncation of propensity score
    #     L,              # Lipschitz parameter for the loss function
    #     cost,           # Cost parameter

    list(
        lambda = 2*L*xi*kappa*sqrt(F)/(epsilon*n*c_L)*((n-1)*zeta/c_L + 2),
        lambda.range = range,
        solve_for = solve_for, n = n, beta = beta, Delta = Delta, epsilon = epsilon, d_tilde = d_tilde,
        set_parms = list(
            used_F = F, used_kappa = kappa, kernel_parms = kernel_parms, xi = xi, zeta = zeta, c_L = cLower, c_H = cUpper, L = L, cost = cost, d = d
        ))
}

# Predict from a PrOWL Vector Helper
predict_prowl <- function(w, newx) { 
    as.numeric(sign(cbind(1, newx) %*% w)) 
}

make_transform_function <- function(fm) {
    function(model_x) {
        t(apply(model_x, MARGIN=1, function(x){fm$mapping(x)}))
    }
}

# PROWL: 
# Allowable Kernels:
    # rbfdot      # Radial Basis kernel "Gaussian"
    # polydot     # Polynomial kernel
    # vanilladot  # Linear kernel
    # laplacedot  # Laplacian kernel
PrOWL <- function(df_formula,
                  data,
                  kernel_type = "rbfdot",
                  cost = NULL,
                  treatment_var,
                  reward_var,
                  prob_var,
                  lambda,
                  epsilon,
                  d = NULL,
                  d_tilde = NULL,
                  ...) {
    
    # Model Matrix of Predictors
    model_x <- model.matrix(df_formula, data)
    
    # Map kernel_type to correct kernel_type
    kernel_type <- tolower(kernel_type)
    if(startsWith(kernel_type, "p")) {
        kernel_type <- "polydot"
        transform_function <- function(model_x){ model_x }
    } else if(startsWith(kernel_type, "v")) {
        kernel_type <- "vanilladot"
        transform_function <- function(model_x){ model_x }
    } else { 
        # Infinite Dimensional Kernels
        kernel_type <- "rbfdot"
        sample_fn <- function(d, d_tilde){
            do.call(rbind, lapply(1:d, function(ii) { rnorm(d_tilde) }))
        }


        if(startsWith(kernel_type, "l")) {
            kernel_type <- "laplacedot"
            
            sample_fn <- function(d, d_tilde){
                do.call(rbind, lapply(1:d, function(ii) { rcauchy(d_tilde) }))
            }

        } else if(! startsWith(kernel_type, "r")) { 
            warning("The provided kernel did not match any of `rbfdot`, `polydot`, `vanilladot`, or `laplacedot`. The results use `rbfdot`.") 
        }
        # We need to do the kernel approximation here.
        fm <- approximate_feature_mapping(d, d_tilde, sample_fn)

        transform_function <- make_transform_function(fm)

        # kernel_fn <- define_kernel(fm$varphi)
        model_x <- transform_function(model_x)
    }


    # Extract 
    treatments <- as.numeric(unlist(data[treatment_var]))
    rewards <- as.numeric(unlist(data[reward_var]))
    probs <- as.numeric(unlist(data[prob_var]))
    weight <- rewards / ((treatments == 1)*probs + (treatments == -1)*(1 - probs))

    # Fit the SVM
    w_svm <- weighted.ksvm(y = treatments, 
                           x = model_x, 
                           weights = weight,
                           C = cost,
                           kernel = "vanilladot",
                           nfolds = 1,
                           ...)

    # Grab the Vector
    w_tilde <- c(-1*w_svm$dual, colSums(w_svm$y * w_svm$primal * model_x))

    # Privatize it 
    w <- privatize_vector(w_tilde, lambda, epsilon)

    # Return the estimated value
    predicted_vals <- predict_prowl(w, model_x)
    ev <- estimated_value(predicted_vals, treatments, rewards, probs)

    # Build Object to Return
    list(
        method = "PrOWL",
        decision_function_model = df_formula,
        fitted_svm_vector = w,
        value = ev,
        transform_function = transform_function, 
        predict_treatment = function(ind_data) {
            predict_prowl(w, newx = transform_function(ind_data))
        }
    )
}

