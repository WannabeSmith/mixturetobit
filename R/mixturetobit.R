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
llMixtureTobit <- function(y, K, w_ik, lambda, xTbeta, sigma, n_i)
{
  # A list of the K components of the log-likelihood
  ll.components <- lapply(1:K, function(k){
    sum(w_ik[[k]] * (log(lambda[k])/n_i +
                       dtobit(y = y, xTbeta = xTbeta[[k]],
                              sigma = sigma[k], log = TRUE)))
  })

  ll <- do.call(sum, ll.components)

  return(ll)
}

prep.theta.optim <- function(beta, sigma, K)
{
  theta.init <- lapply(1:K, function(k){
    return(c(beta[[k]], sigma[k]))
  })

  return(theta.init)
}

recover.params <- function(theta, K)
{
  beta <- lapply(1:K, function(k){
    length.theta_k <- length(theta[[k]])
    return(theta[[k]][-length.theta_k])
  })

  sigma <- lapply(1:K, function(k){
    length.theta_k <- length(theta[[k]])
    return(theta[[k]][length.theta_k])
  })

  sigma <- do.call(c, sigma)

  return(list(beta = beta, sigma = sigma))
}

Q_k.tobit <- function(theta_k, y, X, delta_k, lambda.tplus1_k, n_i)
{
  beta <- theta_k[-length(theta_k)]
  sigma <- theta_k[length(theta_k)]

  xTbeta <- X %*% beta

  ll <- sum(delta_k * (log(lambda.tplus1_k)/n_i +
                         dtobit(y = y, xTbeta = xTbeta,
                                sigma = sigma, log = TRUE)))

  return(ll)
}

Q_k.tobit <- cmpfun(Q_k.tobit)

EM.tobit <- function(y, start.beta, start.sigma, start.lambda, K, ll.prev, X, id.vec, n_i,
               theta.lower = NULL, theta.upper = NULL, method = "L-BFGS-B", tol = 1e-5,
               best.beta, best.sigma, best.lambda, best.delta = NULL, best.ll)
{
  xTbeta <- lapply(start.beta, function(b){
    return(X %*% b)
  })

  # Components of f(y | theta). That is, f(y | theta) = \sum_k^K[\lambda_k \times f_k(y | theta)]
  f.y.theta.comps <- lapply(1:K, function(k){
    f.y.theta_k <- data.table(aggregate(dtobit(y, xTbeta = xTbeta[[k]], sigma = start.sigma[[k]]),
                             by = list(id.vec), function(x){prod(x) * start.lambda[k]}))
    names(f.y.theta_k) <- c("id", "f.y.theta_k")

    merged <- merge(data.table("id" = id.vec), f.y.theta_k, by = "id")
    return(merged$f.y.theta_k)
  })

  # add these vectors together element-wise to get a single vector of size nrow(X)
  f.y.theta <- rowSums(do.call(cbind, f.y.theta.comps))

  # Estimated probabilities of each observation belonging to group k
  delta <- lapply(1:K, function(k){
    f.y.theta.comps[[k]] / f.y.theta
  })

  deltasum <- rowSums(do.call(cbind, delta))
  # Have to use all.equal because, amazingly, == does not recognize
  # min(deltasum) as 1 even though it prints exactly 1.
  # More here: https://www.reddit.com/r/rstats/comments/34p0bm/logical_operators_not_working_as_expected_r_says/
  stopifnot(all.equal(min(deltasum), 1) && all.equal(max(deltasum), 1))

  lambda.tplus1 <- sapply(delta, mean)

  theta.init <- prep.theta.optim(beta = start.beta, sigma = start.sigma, K = K)

  print("Trying beta...")
  print(theta.init)

  print("Previous lambda...")
  print(start.lambda)

  print(paste("Previous log-likelihood:", ll.prev))

  # Create the objective function (the negative of the Q function)
  # Doing this inside the EM function so that we don't need
  # to pass parameters that don't change (X, K, delta, ...)
  Js <- lapply(1:K, function(k){
    J <- function(theta_k){-Q_k.tobit(theta_k = theta_k, y = y, X = X,
                                delta_k = delta[[k]],
                                lambda.tplus1_k = lambda.tplus1[k],
                                n_i = n_i)}
    return(cmpfun(J))
  })


  if(method == "Nelder-Mead")
  {
    optims <- lapply(1:K, function(k){
      optim(theta.init[[k]], Js[[k]], method = "Nelder-Mead")
    })
  } else if(method == "L-BFGS-B")
  {
    if(!is.null(theta.lower) && !is.null(theta.upper))
    {
      optims <- lapply(1:K, function(k){
        optim(theta.init[[k]], Js[[k]], method = "L-BFGS-B",
              lower = theta.lower, upper = theta.upper)
      })
    } else
    {
      optims <- lapply(1:K, function(k){
        optim(theta.init[[k]], Js[[k]], method = method)
      })
    }
  }

  optimal.params <- lapply(optims, function(o){
    return(o$par)
  })

  optima <- lapply(optims, function(o){
    return(o$value)
  })

  params <- recover.params(optimal.params, K = K)
  beta.tplus1 <- params$beta
  sigma.tplus1 <- params$sigma

  xTbeta.tplus1 <- lapply(beta.tplus1, function(b){
    return(X %*% b)
  })

  # Get predicted latent class for each observation
  pred.latent.class <- max.col(do.call(cbind, delta))

  # Create list of K vectors where the k'th vector is a vector of logicals
  # indicating class membership to latent class k. This is required for the
  # log-likelihood to be evaluated.
  # w_ik <- lapply(1:K, function(k){
  #   return(pred.latent.class == k)
  # })

  ll <- -do.call(sum, optima)

  if(ll < ll.prev)
  {
    warning("E-step did not increase log-likelihood!")
  }

  if(ll > best.ll)
  {
    best.ll <- ll
    best.beta <- beta.tplus1
    best.lambda <- lambda.tplus1
    best.sigma <- sigma.tplus1
    best.delta <- delta
  }

  print(paste("log-likelihood:", ll))

  if(abs(ll/ll.prev - 1) < tol)
  {
    print("converged")
    return(list(beta = best.beta,
                sigma = best.sigma,
                lambda = best.lambda,
                delta = best.delta,
                ll = best.ll))
  } else
  {
    return(EM.tobit(y = y, start.beta = beta.tplus1, start.sigma = sigma.tplus1,
              start.lambda = lambda.tplus1, K = K, id.vec = id.vec, n_i = n_i, ll.prev = ll, X = X,
              method = method, theta.lower = theta.lower,
              theta.upper = theta.upper, tol = tol, best.beta = best.beta, best.sigma = best.sigma,
              best.lambda = best.lambda, best.delta = best.delta, best.ll = best.ll))
  }
}

#' Perform maximum-likelihood estimation for a mixture of tobit regression models
#'
#' This function performs maximum-likelihood estimation via the E-M algorithm to obtain
#' estimates of regression coefficients in a mixture of tobit regression models.
#'
#' @importFrom stats model.matrix dnorm pnorm rnorm optim aggregate
#' @importFrom survival survreg
#' @importFrom compiler cmpfun
#' @import data.table
#' @param formula a regression formula describing the relationship between the response and the covariates
#' @param data the data.frame containing the responses and covariates
#' @param K the number of mixtures (or latent classes)
#' @param start.beta a list of length K of starting values for each mixture's beta coefficients
#' @param start.sigma a vector of length K of starting values for each mixture's sigma value
#' @param start.lambda a vector of length K of starting values for the mixing proportions
#' @param id a string specifying the name of the column that identifies subjects
#' @param left a number specifying where left-censoring occurred
#' @param tol a number specifying the tolerance used to determine convergence
#' @param theta.lower a numeric vector of lower bounds for the theta parameters
#' @param theta.upper a numeric vector of upper bounds for the theta parameters
#' @param method a string specifying the optimization routine to be used by optim
#' @return a list containing the following elements:\cr
#' \item{beta}{a list containing the estimated regression coefficients}
#' \item{sigma}{a vector containing the estimated values of sigma}
#' \item{lambda}{a vector containing the estimated mixing proportions}
#' \item{delta}{a list of length K containing the estimated class membership probabilities for each observation}
#' \item{ll}{the log-likelihood function evaluated at the MLE}
#' @export
mixturetobit <- function(formula, data, K = 2, start.beta = NULL,
                         start.sigma = NULL, start.lambda = NULL, id = NULL,
                         left = -1, tol = 1e-5, theta.lower = NULL,
                         theta.upper = NULL, method = "L-BFGS-B")
{
  data <- data.table(data)

  stopifnot(K >= 2)

  X <- model.matrix(formula, data = data)
  d <- ncol(X) # dimension of beta

  response.varname <- all.vars(formula)[1]
  y <- data[[response.varname]]

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
      return(rnorm(n = d, mean = naive.beta, sd = 0.5)) # Assign some random starting points
    })

    start.sigma <- rep(naive.model$scale, K)/(2*K) # Give same sigma to each group

    start.lambda <- rep(1/K, K)
  }

  id.vec <- data[[id]]

  n_i <- data[, n_i := .N, by = eval(id)]$n_i

  MLE <- EM.tobit(y = y, start.beta = start.beta, start.sigma = start.sigma,
            start.lambda = start.lambda, K = K, ll.prev = -Inf, X = X, id.vec = id.vec, n_i = n_i,
            theta.lower = theta.lower, theta.upper = theta.upper, method = method, tol = tol,
            best.beta = start.beta, best.sigma = start.sigma, best.lambda = start.lambda,
            best.ll = -Inf)

  return(MLE)
}

# Incorporate the following as an example later
#
# K=3
# formula <- tto ~ mo + sc + ua + pd + ad
# theta.lower <- c(rep(-Inf, 1 * 21), rep(1e-16, 1))
# theta.upper <- c(rep(Inf, 1 * 21), rep(Inf, 1))
# start.lambda <- rep(1/K, K)
#
# start.beta <- list(beta.tto.true[[1]] + 1, beta.tto.true[[2]], beta.tto.true[[3]])
# start.sigma <- sigma.tto + c(0.1, -0.1, 0)
#
# # start.beta <- beta.tto.true
# # start.sigma <- sigma.tto
# set.seed(75)
# system.time(MLE.KSO <- mixturetobit(formula, data = data, K = K, start.beta = start.beta,
#                     start.sigma = start.sigma, start.lambda = start.lambda, id = "SubjectID",
#                     theta.lower = theta.lower, theta.upper = theta.upper, method = "L-BFGS-B", tol = 1e-8))
