
get_initial_MV_log_chisq <- function(params, y){
  variance <- 1 / ((1 - params$phi^2) / params$sigma^2 + 1 / trigamma(.5))
  mean <- variance * ((log(y^2) - log(2) - digamma(.5) - 2 * log(params$beta)) / trigamma(.5))
  c(mean = mean, variance = variance)
}

get_initial_particles <- function(N, 
                                  params,
                                  y, 
                                  method = "bootstrap"){
  if(method == "bootstrap"){
    return(params$sigma / sqrt(1 - params$phi^2) * rnorm(N))
  }
  if(method == "log_chisq"){
    mean_variance <- get_initial_MV_log_chisq(params, y)
    return(rnorm(N, mean_variance["mean"], sqrt(mean_variance["variance"])))
  }
}


get_initial_weights <- function(y, 
                                particles,
                                params,
                                method = "bootstrap", 
                                log_ = TRUE, normed = TRUE){
  if(method == "bootstrap"){
    log_weights <- dnorm(y, 0, params$beta * exp(.5 * particles), 
                         log = TRUE)
  }
  else if(method == "log_chisq"){
    mean_variance <- get_initial_MV_log_chisq(params, y)
    log_weights <- dnorm(y, 0, params$beta * exp(.5 * particles), 
                         log = TRUE) + 
      dnorm(particles, 0, params$sigma / sqrt(1 - params$phi^2), log = TRUE) -
      dnorm(particles, mean_variance["mean"], sqrt(mean_variance["variance"]),
            log = TRUE)
  }
  if(!log_){
    if(normed){
      return(exp(log_weights) / sum(exp(log_weights)))
    }
    return(exp(log_weights))
  }
  else
    return(log_weights)
}

get_next_MV_log_chisq <- function(params, current_particles, y){
  variance <- params$sigma^2 * trigamma(.5) / (params$sigma^2 + trigamma(.5)) 
  means <- variance * (params$phi * current_particles / params$sigma^2 + 
                         (log(y^2) - log(2) - digamma(.5) - 2 * log(params$beta)) / trigamma(.5))
  list(mean = means, variance = variance)
}

get_next_particles <- function(current_particles,
                               params,
                               y,
                               method = "bootstrap"){
  N <- length(current_particles)
  if(method == "bootstrap"){
    return(rnorm(N, params$phi * current_particles, params$sigma))
  }
  if(method == "log_chisq"){
    mean_variance <- get_next_MV_log_chisq(params, current_particles, y)
    return(rnorm(N, mean_variance[["mean"]], sqrt(mean_variance[["variance"]])))
  }
}

get_next_weights <- function(y,
                             current_particles,
                             next_particles,
                             params,
                             method = "bootstrap", 
                             log_ = TRUE, normed = TRUE){
  if(method == "bootstrap"){
    log_weights <- dnorm(y, 0, 
                         params$beta * exp(.5 * next_particles), log = TRUE)
  }
  else if(method == "log_chisq"){
    mean_variance <- get_next_MV_log_chisq(params, current_particles, y)
    log_weights <- dnorm(y, 0, 
                         params$beta * exp(.5 * next_particles), log = TRUE) +
      dnorm(next_particles, params$phi * current_particles, 
            params$sigma, log = TRUE) -
      dnorm(next_particles,  mean_variance[["mean"]], 
            sqrt(mean_variance[["variance"]]), log = TRUE)
  }
  if(!log_){
    if(normed){
      return(exp(log_weights) / sum(exp(log_weights)))
    }
    return(exp(log_weights))
  }
  else
    return(log_weights)
}

# cll stands for complete log likelihood

get_initial_cll <- function(y, initial_particles, params){
  dnorm(initial_particles, 0, params$sigma / sqrt(1 - params$phi^2), 
        log = TRUE) +
    dnorm(y, 0, params$beta * exp(.5 * initial_particles), 
          log = TRUE)
}

get_transition_cll <- function(y, parent_particles, child_particle, params){
  result <- dnorm(y, 0, params$beta * exp(.5 * child_particle), 
                  log = TRUE) +
    dnorm(child_particle, params$phi * parent_particles, params$sigma, 
          log = TRUE)
}

get_backward_weights <- function(parent_particles, child_particle, params, 
                                 normed = TRUE){
  weights <- dnorm(child_particle, params$phi * parent_particles, params$sigma)
  if(normed)
    weights <- weights / sum(weights)
  return(weights)
}

get_parameter_offspring <- function(current_param, innovation_variance,
                                    known_beta = NULL, known_sigma = NULL){
  ok <- FALSE
  while(!ok){
    candidate <- map(current_param,
                     function(x) x * exp(rnorm(1, - .5 * innovation_variance, 
                                               sqrt(innovation_variance))))
    if(!is.null(known_sigma)){
      candidate$sigma <- known_sigma
    }
    if(!is.null(known_beta)){
      candidate$beta <- known_beta
    }
    ok <- (abs(candidate$phi) < 1) & (candidate$sigma > 0)
  }
  return(candidate)
}

get_initial_guess <- function(ys, known_beta = NULL, known_sigma = NULL){
  zs <- log(ys^2) - (log(2) + digamma(.5))
  log_beta <- ifelse(is.null(known_beta), mean(zs), log(known_beta))
  sigma_phi_ratio <- var(zs) - trigamma(.5)
  phi <- min(cov(zs[-1], zs[-length(zs)]) / sigma_phi_ratio, .99)
  sigma <- ifelse(is.null(known_sigma), sqrt(sigma_phi_ratio * (1 - phi^2)),
                  known_sigma)
  list(beta = exp(log_beta), phi = phi, sigma = sigma)
}

get_particle_filter <- function(Ys, params, method = "bootstrap",
                                N = 50, 
                                smoothing = c("PoorMan", "FFBSi", "FFBSm", "BIS"),
                                IS_method = "IS",
                                get_N_tilde = function(N) floor(N^(.6))){
  # Initialization
  n_obs <- length(Ys)
  particles <- filtering_weights <- unnormed_weights <- genealogy <- matrix(NA, nrow = N, 
                                                                            ncol = n_obs)
  E_statistics <- NULL
  loglik <- 0
  for(t in 1:n_obs){
    if(t == 1){ # Initialization step
      particles[, t] <- get_initial_particles(N = N, params = params, y = Ys[t], method = method)
      log_unnormed_weights <- get_initial_weights(Ys[t], 
                                                  particles[, t], 
                                                  params, 
                                                  log = TRUE, 
                                                  normed = FALSE,
                                                  method = method)
      genealogy[, t] <- 1:N
    }
    else if(t > 1){ # Propagation step
      selected_indexes <- sample(1:N, 
                                 size = N,
                                 prob = filtering_weights[, t - 1],
                                 replace = TRUE)
      selected_particles <- particles[selected_indexes, t - 1]
      particles[, t] <- get_next_particles(current_particles = selected_particles, 
                                           params = params, 
                                           y = Ys[t],
                                           method = method)
      log_unnormed_weights <- get_next_weights(y = Ys[t], 
                                               current_particles = selected_particles, 
                                               next_particles = particles[, t], 
                                               params = params, 
                                               log_ = TRUE, normed = FALSE, 
                                               method = method)
      new_genealogy <- genealogy[selected_indexes,]
      new_genealogy[, n_obs] <- 1:N
      new_genealogy[, t - 1] <- selected_indexes
      genealogy <- new_genealogy
    }
    max_log_uw <- max(log_unnormed_weights)
    unnormed_weights[, t] <- exp(log_unnormed_weights - max_log_uw)
    loglik <- loglik + max_log_uw + log(mean(unnormed_weights[, t]))
    filtering_weights[, t] <- unnormed_weights[, t] / sum(unnormed_weights[, t])
  }
  smoothing_particles <- smoothing_weights <- backward_indices <- backward_weights <- NULL
  if(smoothing == "FFBSi"){
    # Smoothing
    smoothing_particles <- particles
    for(t in n_obs:1){
      if(t == n_obs){
        smoothing_particles[, t] <- sample(particles[, n_obs], size = N, replace = TRUE,
                                           prob = filtering_weights[, t])
      }
      else{
        smoothing_particles[, t] <- sapply(smoothing_particles[, t + 1],
                                           function(x){
                                             log_u_weights <- log(unnormed_weights[, t]) + 
                                               dnorm(x, params$phi * particles[, t], 
                                                     params$sigma, log = TRUE)
                                             max_log_uw <- max(log_unnormed_weights)
                                             unnormed_weights <- exp(log_unnormed_weights - max_log_uw)
                                             sample(particles[, t], size = 1,
                                                    prob = unnormed_weights / sum(unnormed_weights))
                                           })
      }
    }
    smoothing_weights <- matrix(1 / N, nrow = N, ncol = n_obs)
    E_statistics <- get_E_statistics(smoothing_particles,
                                     y = Ys, W = smoothing_weights[, n_obs])
  }
  else if(smoothing == "FFBSm"){
    # Smoothing
    smoothing_weights <- filtering_weights
    for(t in (n_obs - 1):1){
      smoothing_weights[, t] <- sapply(smoothing_particles[, t + 1],
                                       function(x){
                                         log_u_weights <- log(filtering_weights[, t]) + 
                                           dnorm(x, params$phi * particles[, t], params$sigma, log = TRUE)
                                         max_log_uw <- max(log_unnormed_weights)
                                         unnormed_weights <- exp(log_unnormed_weights - max_log_uw)
                                         return(sum(unnormed_weights))
                                       })
    }
    smoothing_weights <- apply(smoothing_weights, 2, function(x) x / sum(x))
    smoothing_particles <- particles
  }
  else if(smoothing == "BIS"){
    N_tilde <- get_N_tilde(N)
    E_statistics <- get_E_statistics_BIS(particles, Ys, filtering_weights,
                                         N_tilde, params, method = IS_method)
  }
  else if(smoothing == "PoorMan"){
    smoothing_weights <- matrix(filtering_weights[, n_obs], ncol = n_obs,
                                nrow = N)
    smoothing_particles <- apply(genealogy, 1, 
                                 function(gen_index) 
                                   mapply(gen_index, 1:n_obs,
                                          FUN = function(i, t) particles[i, t])) %>% 
      t()
    E_statistics <- get_E_statistics(smoothing_particles,
                                     y = Ys, W = smoothing_weights[, n_obs])
  }
  return(list(filtering_particles = particles,
              filtering_weights = filtering_weights,
              smoothing_particles = smoothing_particles,
              smoothing_weights = smoothing_weights,
              backward_indices = backward_indices,
              backward_weights = backward_weights,
              loglik = loglik,
              E_statistics = E_statistics))
}

get_complete_log_likelihood <- function(x, y, params){
  beta <- params[1]
  phi <- params[2]
  sigma <- params[3]
  n_obs <- length(x)
  (
    sum(dnorm(x, phi * dplyr::lag(x), sigma, log = TRUE),
        na.rm = TRUE) +
      sum(dnorm(y, 0, beta * exp(.5 * x), log = TRUE)) +
      dnorm(x[1], 0, sigma / sqrt(1 - phi^2), log = TRUE)) / n_obs
}

get_cll <- function(x, y, params){
  beta <- params[1]
  phi <- params[2]
  sigma <- params[3]
  n_obs <- length(x)
  s0 <- x[1]^2
  s1 <- sum(x[-n_obs]^2)
  s2 <- sum(x[-1]^2)
  s3 <- sum(x[-1] * x[-n_obs])
  s4 <- sum(y^2 * exp(-x))
  (-.5 * -n_obs * (log(beta^2) + log(sigma^2)) +
      .5 * log(1 - phi^2) -
      .5 * (s4 / beta^2 + (1 - phi^2) * s0 / sigma^2 + (s2 - 2 * phi * s3 + phi^2 * s1))) / n_obs
}

get_E_statistics <- function(X, y, W = NULL){
  if(is.null(W))
    W <- rep(1/nrow(X), nrow = nrow(X), ncol = ncol(X))
  n_obs <- length(y)
  s0 <- sum(W * X[, 1]^2)
  s1 <- sum(W * rowSums(X[,-n_obs]^2))
  s2 <- sum(W * rowSums(X[,-1]^2))
  s3 <- sum(W * apply(X, 1, function(x) sum(x[-1] * x[-n_obs])))
  s4 <- sum(W * apply(X, 1, function(x) sum(y^2 * exp(-x))))
  c(s0 = s0, s1 = s1, s2 = s2, s3 = s3, s4 = s4)
}

get_E_statistics_BIS <- function(particles, Ys, filtering_weights, 
                                 N_tilde, current_param,
                                 method = c("IS", "Exact")){
  N <- nrow(particles)
  n_obs <- length(Ys)
  tau_s0 <- tau_s1 <- tau_s2 <- tau_s3 <- tau_s4  <- matrix(0, nrow = nrow(particles), 
                                                            ncol = ncol(particles))
  for(t in 1:(n_obs - 1)){
    for(i in 1:N){
      if(method == "IS"){
        backward_indexes <- sample(1:N, 
                                   size = N_tilde,
                                   prob = filtering_weights[, t],
                                   replace = TRUE)
        backward_parts <- particles[backward_indexes, t]
        backward_weights <- get_backward_weights(backward_parts,  
                                                 particles[i, t + 1], 
                                                 current_param)
      }
      else if(method == "Exact"){
        backward_indexes <- 1:nrow(particles)
        backward_parts <- particles[backward_indexes, t]
        u_bw <- filtering_weights[, t] * dnorm(particles[i, t + 1],
                                               current_param$phi * backward_parts,
                                               current_param$sigma)
        backward_weights <- u_bw / sum(u_bw)
      }
      tau_s0[i, t + 1] <- sum(backward_weights * tau_s0[backward_indexes, t])
      if(t == 1){
        tau_s0[i, t + 1] <- tau_s0[i, t + 1] + 
          sum(backward_weights * backward_parts^2)
      }
      tau_s1[i, t + 1] <- sum(backward_weights * (tau_s1[backward_indexes, t] + 
                                                    backward_parts^2)) 
      tau_s2[i, t + 1] <- sum(backward_weights * (tau_s2[backward_indexes, t] +
                                                    particles[i, t + 1]^2))
      tau_s3[i, t + 1] <- sum(backward_weights * (tau_s3[backward_indexes, t] + 
                                                    particles[i, t + 1] * backward_parts))
      tau_s4[i, t + 1] <- sum(backward_weights * (tau_s4[backward_indexes, t] + 
                                                    Ys[t]^2 * exp(-backward_parts)))
      if(t == n_obs - 1){
        tau_s4[i, t + 1] <- tau_s4[i, t + 1] + 
          Ys[t + 1]^2 * exp(-particles[i, t + 1])
      }
    }
  }
  s0 <- sum(filtering_weights[, n_obs] * tau_s0[, n_obs])
  s1 <- sum(filtering_weights[, n_obs] * tau_s1[, n_obs])
  s2 <- sum(filtering_weights[, n_obs] * tau_s2[, n_obs])
  s3 <- sum(filtering_weights[, n_obs] * tau_s3[, n_obs])
  s4 <- sum(filtering_weights[, n_obs] * tau_s4[, n_obs])
  c(s0 = s0, s1 = s1, s2 = s2, s3 = s3, s4 = s4)
}

get_E_step <- function(vec_params, X, y, method = "FFBSi"){
  if(method == "FFBSi"){
    return(-mean(apply(X, 1, get_complete_log_likelihood, 
                       y = y, params = vec_params)))
  }
}

get_exact_M_step <- function(n_obs, statistics){
  s0 <- statistics["s0"]; s1 <- statistics["s1"]
  s2 <- statistics["s2"]; s3 <- statistics["s3"];
  s4 <- statistics["s4"];
  beta <- sqrt(s4 / n_obs)
  phi_cands <- polyroot(c(n_obs * s3, # cst
             -s2 + (n_obs - 1) * s0 - n_obs * s1, # linear
             -(n_obs - 2) * s3,
             (n_obs - 1) * (s1 - s0)))
  phi <- Re(phi_cands[abs(Re(phi_cands)) < 1 & round(Im(phi_cands), digits = 12) == 0])
  sigma <- sqrt((s0 + s2 - 2 * phi * s3 + phi^2 * (s1 - s0)) / n_obs)
  names(beta) <- names(sigma) <- NULL
  list(beta = beta, phi = phi, sigma = sigma)
}

get_approx_M_step <- function(current_param, X, y, method = "FFBSi"){
  foo <- function(param){
    get_E_step(param, X, y, method)
  }
  optim(unlist(current_param))
}
# 
# foo <- function(params, X, y){
#   -mean(apply(X, 1, get_complete_log_likelihood, y = y, params = params))
# }
# optim(unlist(param_true), foo, 
#       lower = c(1e-5, 1e-5, 1e-5),
#       upper = c(10, 0.99999, 10),
#       y = Ys, X = filter(my_data, signal == "True") %>% pull() %>% 
#         matrix(nrow = 2, byrow = TRUE),
#       method = "L-BFGS-B")


# get_EM  <- function(Ys,
#                     method = "log_chisq",
#                     current_guess = NULL,
#                     n_initial_guess = 200,
#                     GEM_variance = 0.005,
#                     known_beta = NULL,
#                     known_sigma = NULL,
#                     N = 100,
#                     n_EM_iterations = 10,
#                     n_params = 50){
#   N_tilde <- floor(N^(.6))
#   n_obs <- length(Ys)
#   if(is.null(current_guess)){
#     first_initial_guess <- get_initial_guess(Ys, known_beta = known_beta,
#                                              known_sigma = known_sigma)
#     initial_guess_candidates <- c(list(first_initial_guess),
#                                   rerun(n_initial_guess, 
#                                         get_parameter_offspring(first_initial_guess, 
#                                                                 .1, known_beta, 
#                                                                 known_sigma = known_sigma)))
#     initial_lls <- mclapply(initial_guess_candidates,
#                             function(param){
#                               get_particle_filter(Ys, param, method, N = 100,
#                                                   get_log_lik = TRUE)$loglik
#                             }, 
#                             mc.cores = max(min(n_initial_guess, detectCores() - 2), 1)) %>% 
#       unlist()
#     current_guess <- initial_guess_candidates[[which.max(initial_lls)]]
#   }
#   old_guess <- current_guess
#   estimated_params <- as.data.frame(current_guess) %>% 
#     mutate(loglik = NA)
#   print("After random search, first param is:")
#   print(unlist(current_guess))
#   for(step_ in 1:n_EM_iterations){
#     # print(step_)
#     tested_params <- c(list(current_guess),
#                        list(old_guess),
#                        rerun(n_params, 
#                              get_parameter_offspring(current_param = current_guess, 
#                                                      GEM_variance,
#                                                      known_beta,
#                                                      known_sigma = known_sigma)))
#     particles <- filtering_weights <- unnormed_weights <- matrix(NA, nrow = N, 
#                                                                  ncol = n_obs)
#     E_step_taus <- array(0, dim = c(length(tested_params),
#                                     N, n_obs))
#     loglik <- 0
#     tic()
#     for(t in 1:n_obs){
#       if(t > 1){ # Propagation step
#         selected_particles <- sample(particles[, t - 1], 
#                                      size = N,
#                                      prob = filtering_weights[, t - 1],
#                                      replace = TRUE)
#         particles[, t] <- get_next_particles(current_particles = selected_particles, 
#                                              params = current_guess, 
#                                              y = Ys[t],
#                                              method = method)
#         log_unnormed_weights <- get_next_weights(y = Ys[t], 
#                                                  current_particles = selected_particles, 
#                                                  next_particles = particles[, t], 
#                                                  params = current_guess, 
#                                                  log_ = TRUE, normed = FALSE, 
#                                                  method = method)
#       }
#       else if(t == 1){ # Initialization step
#         particles[, t] <- get_initial_particles(N = N, 
#                                                 params = current_guess, 
#                                                 y = Ys[t], 
#                                                 method = method)
#         log_unnormed_weights <- get_initial_weights(y = Ys[t], 
#                                                     particles = particles[, t], 
#                                                     params = current_guess, 
#                                                     log = TRUE, 
#                                                     normed = FALSE,
#                                                     method = method)
#       }
#       max_log_uw <- max(log_unnormed_weights)
#       unnormed_weights[, t] <- exp(log_unnormed_weights - max_log_uw)
#       loglik <- loglik + max_log_uw + log(mean(unnormed_weights[, t]))
#       filtering_weights[, t] <- exp(unnormed_weights[, t]) / sum(exp(unnormed_weights[, t]))
#       if(t > 1){
#         E_step_taus[,, t] <- sapply(1:N, function(i){
#           backward_indices <- sample(1:N,
#                                      size = N_tilde,
#                                      prob = filtering_weights[, t - 1],
#                                      replace = TRUE)
#           importance_weights <- get_backward_weights(parent_particles = particles[backward_indices, t - 1], 
#                                                      child_particle = particles[i, t], 
#                                                      params = current_guess)
#           sapply(1:length(tested_params),
#                  function(l){
#                    param <- tested_params[[l]]
#                    htilde <- get_transition_cll(y = Ys[t - 1],
#                                                 parent_particles = particles[backward_indices, t - 1],
#                                                 child_particle = particles[i, t],
#                                                 params = param)
#                    if(t == 2)
#                      htilde <- htilde + get_initial_cll(y = Ys[1],
#                                                         initial_particles = particles[backward_indices, t - 1],
#                                                         params = param)
#                    old_tau <- E_step_taus[l, backward_indices, t - 1]
#                    sum(importance_weights * (htilde + old_tau))
#                  })
#         })
#       }
#     } # End loop over t
#     toc()
#     print("Likelihood of the old parameter")
#     print(loglik / n_obs)
#     E_quantity <- colSums(filtering_weights[, n_obs] * t(E_step_taus[,, n_obs]))
#     old_guess <- current_guess
#     current_guess <- tested_params[[which.max(E_quantity)]]
#     estimated_params <- bind_rows(estimated_params,
#                                   as.data.frame(current_guess) %>% 
#                                     mutate(loglik = NA))
#     
#     estimated_params[step_, "loglik"] <- loglik / n_obs
#     print("New parameter")
#     print(unlist(current_guess))
#   }
#   estimated_params[step_ + 1, "loglik"] <- get_particle_filter(Ys, current_guess,
#                                                                "log_chisq", N = 100,
#                                                                get_log_lik = TRUE)$loglik
#   estimated_params
# }

