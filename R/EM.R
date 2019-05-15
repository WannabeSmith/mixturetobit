# Density of tobit model
dtobit <- function(y, xTbeta, sigma, left = -1, log = FALSE)
{
  d <- y > left

  if (log)
  {
    return(d * (dnorm((y - xTbeta)/sigma, log = TRUE) - log(sigma)) +
             (1 - d) * pnorm((left - xTbeta)/sigma, log.p = TRUE))
  } else
  {
    return((1/sigma * dnorm((y - xTbeta)/sigma))^d * pnorm((xTbeta - left)/sigma, lower.tail = FALSE)^(1-d))
  }
}

# Log-likelihood of mixture of tobit
llMixtureTobit <- function(y, K, w_ik, lambda, xTbeta, sigma)
{
  # A list of the K components of the log-likelihood
  ll.components <- lapply(1:K, function(k){
    sum(w_ik[[k]] * (log(lambda[k]) +
                       dtobit(y = y, xTbeta = xTbeta[[k]],
                              sigma = sigma[k], log = TRUE)))
  })

  ll <- do.call(sum, ll.components)

  return(ll)
}

devectorize.params <- function(theta, K)
{
  # length(theta) - K had better be divisible by K
  stopifnot((length(theta) - K) %% K == 0)

  # dimension of beta
  d <- (length(theta) - K) / K

  theta.split <- split(theta, ceiling(seq_along(theta)/d))

  # The last list will contain our sigmas
  beta <- theta.split[-(K + 1)]
  names(beta) <- NULL
  sigma <- theta.split[K + 1][[1]]

  return(list("beta" = beta,
              "sigma" = sigma))
}

vectorize.params <- function(beta, sigma)
{
  theta <- c(do.call(c, beta), sigma)
  return(theta)
}

# Function to initialize J. J takes in theta and produces a negative log-likelihood.
# init.J takes in K, delta, lambda.tplus1, y, and produces J.
# J must be updated at each E-step but must be held constant for the M-step.
# J is the negative of Q.
init.J <- function(y, X, K, delta, lambda.tplus1)
{
  J <- function(theta)
  {
    devec <- devectorize.params(theta, K)

    beta <- devec$beta
    sigma <- devec$sigma

    xTbeta <- lapply(beta, function(b){
      return(X %*% b)
    })

    return(-llMixtureTobit(y = y, K = K, w_ik = delta, lambda = lambda.tplus1, xTbeta = xTbeta, sigma = sigma))
  }

  return(J)
}

EM <- function(y, start.beta, start.sigma, start.lambda, K, ll.prev, X,
               theta.lower = NULL, theta.upper = NULL, method = "L-BFGS-B", tol = 1e-5)
{
  xTbeta <- lapply(start.beta, function(b){
    return(X %*% b)
  })

  # Components of f(y | theta). That is, f(y | theta) = \sum_k^K[\lambda_k \times f_k(y | theta)]
  f.y.theta.comps <- lapply(1:K, function(k){
    return(start.lambda[k] * dtobit(y, xTbeta = xTbeta[[k]], sigma = start.sigma[[k]]))
  })

  # add these vectors together element-wise to get a single vector of size nrow(X)
  # note: the `+` function is simply the + operator.
  f.y.theta <- do.call(`+`, f.y.theta.comps)

  # Estimated probabilities of each observation belonging to group k
  delta <- lapply(1:K, function(k){
    start.lambda[k] * dtobit(y = y, xTbeta[[k]], start.sigma[[k]]) / f.y.theta
  })

  deltasum <- as.vector(do.call(`+`, delta))
  # Have to use all.equal because, amazingly, == does not recognize
  # min(deltasum) as 1 even though it prints exactly 1.
  # More here: https://www.reddit.com/r/rstats/comments/34p0bm/logical_operators_not_working_as_expected_r_says/
  stopifnot(all.equal(min(deltasum), 1) && all.equal(max(deltasum), 1))

  lambda.tplus1 <- sapply(delta, mean)

  theta.init <- vectorize.params(start.beta, start.sigma)

  print("Trying...")
  print(theta.init)

  # Create the Q function (the negative of the objective function)
  # Doing this inside the EM function so that we don't need
  # to pass parameters that don't change ()
  J <- init.J(y = y, X = X, K = K, delta = delta, lambda.tplus1 = lambda.tplus1)
  J <- compiler::cmpfun(J) # Compile for speed

  if(is.null(theta.lower) && is.null(theta.upper))
  {
    optimum <- optim(theta.init, J, method = "Nelder-Mead")
  } else
  {
    optimum <- optim(theta.init, J, method = "L-BFGS-B", lower = theta.lower,
                     upper = theta.upper)
  }

  devec <- devectorize.params(optimum$par, K = K)
  beta.tplus1 <- devec$beta
  sigma.tplus1 <- devec$sigma

  xTbeta.tplus1 <- lapply(beta.tplus1, function(b){
    return(X %*% b)
  })

  # Get predicted latent class for each observation
  pred.latent.class <- max.col(do.call(cbind, delta))

  # Create list of K vectors where the k'th vector is a vector of logicals
  # indicating class membership to latent class k. This is required for the
  # log-likelihood to be evaluated.
  w_ik <- lapply(1:K, function(k){
    return(pred.latent.class == k)
  })

  ll <- llMixtureTobit(y = y, K = K, w_ik = w_ik, lambda = lambda.tplus1, xTbeta = xTbeta.tplus1, sigma = sigma.tplus1)
  print(paste("log-likelihood:", ll))

  if(abs(ll/ll.prev - 1) < tol)
  {
    print("converged")
    return(list(beta = beta.tplus1,
                sigma = sigma.tplus1,
                lambda = lambda.tplus1,
                w_ik = w_ik, ll = ll,
                delta = delta))
  } else
  {
    return(EM(y = y, start.beta = beta.tplus1, start.sigma = sigma.tplus1,
              start.lambda = lambda.tplus1, K = K, ll.prev = ll, X = X,
              method = method, theta.lower = theta.lower,
              theta.upper = theta.upper, tol = tol))
  }
}

mixturetobit <- function(formula, data, K = 2, start.beta = NULL,
                         start.sigma = NULL, start.lambda = NULL,
                         left = -1, tol = 1e-5, theta.lower = NULL,
                         theta.upper = NULL, method = "L-BFGS-B")
{
  stopifnot(K >= 2)

  X <- model.matrix(formula, data = data)
  d <- ncol(X) # dimension of beta

  response.varname <- all.vars(formula)[1]
  y <- data[, response.varname]

  if(is.null(start.beta) || is.null(start.sigma) || is.null(start.lambda))
  {
    print("Using naive tobit regression model to initialize optimization")
    # create a list with K elements, each of which is a d-dimensional numeric vector
    start.beta <- rep(list(numeric(d)), K)

    naive.model <- survreg(Surv(y, y > left, type = "left") ~
                                           X[, -1], dist = "gaussian")
    naive.beta <- naive.model$coefficients
    names(naive.beta) <- colnames(X)

    start.beta <- lapply(start.beta, function(x){
      return(rnorm(n = d, mean = naive.beta, sd = mean(abs(naive.beta)))) # Assign some random starting points
    })

    start.sigma <- rep(naive.model$scale, K) # Give same sigma to each group
  }

  MLE <- EM(y = y, start.beta = start.beta, start.sigma = start.sigma,
            start.lambda = start.lambda, K = K, ll.prev = Inf, X = X,
            theta.lower = theta.lower, theta.upper = theta.upper, method = method, tol = tol)
}

# Incorporate the following as an example later

# formula <- tto ~ mo + sc + ua + pd + ad
# theta.lower <- c(rep(-3, 2 * 21), rep(1e-16, 2))
# theta.upper <- c(rep(3, 2 * 21), rep(5, 2))
#
#
# # start.beta <- list(beta.tto.true[[1]] + 0.2, beta.tto.true[[2]] - 0.2)
# # start.sigma <- sigma.tto + c(0.1, -0.1)
#
# start.beta <- NULL
# mixturetobit(formula, data = eqdata.tto, K = 2, start.beta = start.beta,
#              start.sigma = start.sigma, start.lambda = c(0.5, 0.5),
#              theta.lower = theta.lower, theta.upper = theta.upper, method = "L", tol = 0.0001)
