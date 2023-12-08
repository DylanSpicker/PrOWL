# Simulations take inspiration from: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3636816/pdf/nihms379905.pdf

generate_single_run_data <- function(n, treatment_p = 0.5, np = 10) {
    # Truncate the data
    cut_off_1 <- 9
    cut_off_4 <- 17.5


    # Generate a single run of the data for a given sample size 'n'
    # Generate np U[-1,1] variates and convert to DF
    gen_data <- as.data.frame(do.call(cbind, lapply(1:np, function(idx){ runif(n, -1, 1) })))
    colnames(gen_data) <- paste0("X", 1:np)

    gen_data$a_p <- treatment_p

    # Generate Treatment Independent of Covariates
    gen_data$A <- rbinom(n, 1, gen_data$a_p)*2 - 1

    # Generate the Set of Reward Outcomes;  Draw from N(Rj, 1)
    offset <- 1 + 2*gen_data$X1 + gen_data$X2 + 0.5*gen_data$X3

    # R1: T0(X,A) = 0.442(1 − X1 − X2)A
    T1 <- 0.442*(1 - gen_data$X1 - gen_data$X2)

    gen_data$R1 <- offset + T1*gen_data$A  + rnorm(n, 0, 1)
    gen_data$R1 <- pmin(pmax(0, gen_data$R1), cut_off_1)

    gen_data$Astar1 <- sign(T1)

    # R4: T0(X,A) = (1 - X1^3 + exp(X3^2 + X5)+0.6X6-(X7+X8)^2)A
    T4 <- (1 - gen_data$X1^3 + exp(gen_data$X3^2 + gen_data$X5) + 0.6*gen_data$X6 - (gen_data$X7 + gen_data$X8)^2)
    gen_data$R4 <- offset + T4*gen_data$A + rnorm(n, 0, 1)
    gen_data$R4 <- pmin(pmax(0, gen_data$R4), cut_off_4)
    gen_data$Astar4 <- sign(T4)

    # Return data frame
    gen_data
}

compute_sim1_privacy_parms <- function(ns, epsilons, beta = 0.85, cost = 0.01) {

    xi_1 <- 9 # [0, 9]
    xi_2 <- 17.5 # [0, 17.5]

    # Other Parameters
    zeta <- 0       # Known probabilities
    c_L <- 0.5
    c_H <- 0.5
    L <- 1          # Hinge is 1-Lipschitz
    d <- 10

    linear_kernel <- list(F = d, kappa = sqrt(d))
    gaussian_kernel <- list(F = Inf, kappa = 1, sigma_p = sqrt(d), diam = 2*sqrt(d))

    privacy_parms <- lapply(ns, function(n){
        each_eps <- lapply(epsilons, function(epsilon){
            # Finite 1
            # Finite 2
            # Infinite 1
            # Infinite 2
            list(
                "finite_A" = compute_prowl_parms(
                    kernel_parms = linear_kernel,
                    xi_1,
                    zeta,
                    cLower = c_L,
                    cUpper = c_H,
                    L,
                    cost,
                    d,
                    n = n,
                    beta = beta,
                    Delta = NULL,
                    epsilon = epsilon),
                "finite_B" = compute_prowl_parms(
                    kernel_parms = linear_kernel,
                    xi_2,
                    zeta,
                    cLower = c_L,
                    cUpper = c_H,
                    L,
                    cost,
                    d,
                    n = n,
                    beta = beta,
                    Delta = NULL,
                    epsilon = epsilon),
                "infinite_A" = compute_prowl_parms(
                    kernel_parms = gaussian_kernel,
                    xi_1,
                    zeta,
                    cLower = c_L,
                    cUpper = c_H,
                    L,
                    cost,
                    d,
                    n = n,
                    beta = beta,
                    Delta = NULL,
                    epsilon = epsilon),
                "infinite_B" = compute_prowl_parms(
                    kernel_parms = gaussian_kernel,
                    xi_2,
                    zeta,
                    cLower = c_L,
                    cUpper = c_H,
                    L,
                    cost,
                    d,
                    n = n,
                    beta = beta,
                    Delta = NULL,
                    epsilon = epsilon))
        })

        names(each_eps) <- epsilons
        each_eps
    })

    names(privacy_parms) <- ns

    privacy_parms
}
