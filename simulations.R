# Number of simulations
K = 10

glmFamily = poisson(link = "log")
n = 512
yFinal = base::matrix(data = NA_integer_,
                      ncol = K,
                      nrow = n)
seedSample = base::rep(x = NA_integer_,
                       times = K)
poisSimple = list()

p = 8 # Total number of variables
betaNonZero = c(2, -0.5, -1.5) # Nonzero coefficients
p1 = base::length(betaNonZero)
p0 = p - p1

beta = c(betaNonZero, rep(x = 0, times = p0))

x = mvtnorm::rmvnorm(n = n,
                     mean = rep(x = 0,
                                times = p),
                     sigma = base::diag(x = p))

base::colnames(x) <- base::paste0("X", 1:p)

eta = x %*% beta
mu = exp(eta)

for (i in 1:K){
  seed = base::sample(x = 1:100000,
                      size = 1)

  y = stats::rpois(n = n, lambda = mu)
  yFinal[,i] = y

  data = data.frame(x,y)

  simRun = paste0("PoisSimple",i)

  poisSimple[[i]] = run_everything(data = data,
                        glmFamily = glmFamily,
                        saveName = simRun)

}


######################################################################
K = 1

glmFamily = poisson(link = "log")

n = 512 # Number of observations

yFinal = matrix(data = NA, ncol = K, nrow = n)
poisMed = list()

p = 15 # Total number of variables
betaNonZero = c(2, -0.5, -1.5, 1, 2, -0.5) # Nonzero coefficients
p1 = length(betaNonZero)
p0 = p - 2 * p1

beta = c(betaNonZero, rep(x = 0, times = p0), betaNonZero)

x = mvtnorm::rmvnorm(n = n, mean = rep(x = 0,times = p), sigma = stats::toeplitz(0.2^(1:p)))

colnames(x) <- paste0("X", 1:p)

eta = x %*% beta
mu = exp(eta)

for (i in 1:K){
  seed = base::sample(x = 1:100000,
                      size = 1)

  y = stats::rpois(n = n, lambda = mu)
  y[which(n + 1 - rank(x[,5]) <= 20)] = rpois(n = 20, lambda = 10)
  yFinal[,i] = y

  data = data.frame(x,y)

  simRun = paste0("PoisMed",i)

  poisMed[[i]] = run_everything(data = data,
                                   glmFamily = glmFamily,
                                   saveName = simRun)

}


