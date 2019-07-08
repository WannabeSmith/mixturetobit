# Density of logit model
dlogit <- function(y, xTbeta, log = FALSE)
{
  p <- exp(xTbeta) / (1 + exp(xTbeta))
  return(dbinom(y, size = 1, prob = p, log = log))
}

# Log-likelihood of mixture of logit
llMixtureLogit <- function(y, K, w_ik, lambda, xTbeta)
{
  # A list of the K components of the log-likelihood
  ll.components <- lapply(1:K, function(k){
    sum(w_ik[[k]] * (log(lambda[k]) +
                       dlogit(y = y, xTbeta = xTbeta[[k]], log = TRUE)))
  })

  ll <- do.call(sum, ll.components)

  return(ll)
}

Q_k <- function(beta_k, y, X, delta_k, lambda.tplus1_k)
{
  xTbeta_k <- X %*% beta_k

  ll <- sum(delta_k * (log(lambda.tplus1_k) +
                         dlogit(y = y, xTbeta = xTbeta_k, log = TRUE)))

  return(ll)
}

Q_k <- cmpfun(Q_k)

EM <- function(y, start.beta, start.sigma, start.lambda, K, ll.prev, X, id.vec = NULL,
               theta.lower = NULL, theta.upper = NULL, method = "L-BFGS-B", tol = 1e-5,
               best.beta, best.lambda, best.ll)
{
  xTbeta <- lapply(start.beta, function(b){
    return(X %*% b)
  })

  # f_k(y | theta) computed for all subjects (not all observations. That is next)
  f.y.theta.comps.subj <- lapply(1:K, function(k){
    f.y.theta_k <- aggregate(dlogit(y, xTbeta = xTbeta[[k]], log = FALSE),
                             by = list(id.vec),
                             function(x){prod(x) * start.lambda[k]})
    return(f.y.theta_k$x)
  })

  # Components of f(y | theta). That is, f(y | theta) = \sum_k^K[\lambda_k \times f_k(y | theta)]
  # This returns for all observations, not all subjects
  f.y.theta.comps <- lapply(1:K, function(k){
    f.y.theta_k <- data.table("f.y.theta_k" = f.y.theta.comps.subj[[k]])
    f.y.theta_k$id <- unique(id.vec)

    merged <- merge(data.table("id" = id.vec), f.y.theta_k, by = "id")
    return(merged$f.y.theta_k)
  })

  # add these vectors together element-wise to get a single vector of size nrow(X)
  f.y.theta.subj <- rowSums(do.call(cbind, f.y.theta.comps.subj))
  f.y.theta <- rowSums(do.call(cbind, f.y.theta.comps))

  delta.subj <- lapply(1:K, function(k){
    f.y.theta.comps.subj[[k]] / f.y.theta.subj
  })

  # Estimated probabilities of each observation belonging to group k
  delta <- lapply(1:K, function(k){
    f.y.theta.comps[[k]] / f.y.theta
  })

  deltasum <- rowSums(do.call(cbind, delta))
  # Have to use all.equal because, amazingly, == does not recognize
  # min(deltasum) as 1 even though it prints exactly 1.
  # More here: https://www.reddit.com/r/rstats/comments/34p0bm/logical_operators_not_working_as_expected_r_says/
  stopifnot(all.equal(min(deltasum), 1) && all.equal(max(deltasum), 1))

  start.params <- c(start.lambda, 1)

  lambda.tplus1 <- sapply(delta, mean)

  if(min(lambda.tplus1) < 0.01)
  {
    stop("Smallest lambda is less than 1%, likely indicating convergence issues")
  }

  theta.init <- start.beta

  print("Trying beta...")
  print(theta.init)

  print("Previous lambda...")
  print(start.lambda)

  print(paste("Previous log-likelihood:", ll.prev))

  # Create the objective function (the negative of the Q function)
  # Doing this inside the EM function so that we don't need
  # to pass parameters that don't change (X, K, delta, ...)
  Js <- lapply(1:K, function(k){
    J <- function(theta_k){-Q_k(beta_k = theta_k, y = y, X = X,
                                delta_k = delta[[k]], lambda.tplus1_k = lambda.tplus1[k])}
    return(cmpfun(J))
  })


  optimal.params <- lapply(1:K, function(k){
    optim(theta.init[[k]], Js[[k]], method = "Nelder-Mead")$par
  })

  optima <- lapply(1:K, function(k){
    optim(theta.init[[k]], Js[[k]], method = "Nelder-Mead")$value
  })

  beta.tplus1 <- optimal.params

  xTbeta.tplus1 <- lapply(beta.tplus1, function(b){
    return(X %*% b)
  })

  # Get predicted latent class for each observation
  # pred.latent.class <- max.col(do.call(cbind, delta))

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
  }

  if(abs(ll/ll.prev - 1) < tol)
  {
    print("converged")
    return(list(beta = beta.tplus1,
                lambda = lambda.tplus1,
                delta = delta,
                best.ll = best.ll))
  } else
  {
    return(EM(y = y, start.beta = beta.tplus1, start.lambda = lambda.tplus1,
              K = K, id.vec = id.vec, ll.prev = ll, X = X,
              method = method, theta.lower = theta.lower,
              theta.upper = theta.upper, tol = tol, best.beta = best.beta,
              best.lambda = best.lambda, best.ll = best.ll))
  }
}

#' Perform maximum-likelihood estimation for a latent class logistic regression model
#'
#' This function performs maximum-likelihood estimation via the E-M algorithm to obtain
#' estimates of regression coefficients in a latent class logistic regression model.
#'
#' @importFrom stats model.matrix dnorm pnorm rnorm optim
#' @importFrom survival survreg
#' @importFrom compiler cmpfun
#' @param formula a regression formula describing the relationship between the response and the covariates
#' @param data the data.frame containing the responses and covariates
#' @param K the number of mixtures (or latent classes)
#' @param start.beta a list of length K of starting values for each mixture's beta coefficients
#' @param start.lambda a vector of length K of starting values for the mixing proportions
#' @param id the (character) name of the column containing subject IDs
#' @param tol a numeric tolerance used to determine convergence
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
latentclasslogit <- function(formula, data, K = 2, start.beta = NULL,
                         start.lambda = NULL, id = NULL, tol = 1e-5, theta.lower = NULL,
                         theta.upper = NULL, method = "L-BFGS-B")
{
  data <- data.table(data)

  stopifnot(K >= 2)

  X <- model.matrix(formula, data = data)
  d <- ncol(X) # dimension of beta

  response.varname <- all.vars(formula)[1]
  y <- data[[response.varname]]

  id.vec <- data[[id]]

  MLE <- EM(y = y, start.beta = start.beta, start.lambda = start.lambda, K = K,
            ll.prev = -Inf, X = X, id.vec = id.vec,theta.lower = theta.lower,
            theta.upper = theta.upper, method = method, tol = tol, best.beta = start.beta,
            best.lambda = start.lambda, best.ll = -Inf)

  return(MLE)
}

# Incorporate the following as an example later
# #
# K=3
# formula <- formula.dce
# theta.lower <- c(rep(-Inf, 1 * 21))
# theta.upper <- c(rep(Inf, 1 * 21))
# # theta.lower <- NULL
# # theta.upper <- NULL
# start.lambda <- rep(1/K, K)
#
# start.beta <- beta.dce.true
#
# # start.beta <- beta.dce.true
# # start.sigma <- sigma.dce
# set.seed(75)
# system.time(MLE.KSO <- latentclasslogit(formula, data = data, K = K, start.beta = start.beta,
#                                         start.lambda = start.lambda, id = "SubjectID",
#                                         theta.lower = theta.lower, theta.upper = theta.upper,
#                                         method = "L", tol = 1e-8))
