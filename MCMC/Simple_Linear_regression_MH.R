# Purpose:
#   perform iid draws from posterior of simple linear regression model using
#     Metropolis-Hastings method.
#
# Model:
#   y = beta_0 + X * beta_1 + e  e ~N(0, sigmasq)
#          y is n x 1
#          X is n x 2
#          beta is 2 x 1 vector of coefficients
#
# Priors:  beta ~ N(betabar, sigmasq*A^-1)
#          beta_0 ~ N(0, 1)
#          beta_1 ~ N(1, 2)
#          sigmasq ~ gamma(1, 2)
# 
#
# check arguments
#
#########################################################################
# Generate test dataset 
set.seed(100)
N <- 1000; betas <- c(2, 4); true_sd <- 10
X <<- cbind(1, (-(N - 1) / 2):((N - 1)/2))
y <<- X %*% betas + rnorm(n = N, mean = 0, sd = true_sd)
plot(X[, 2], y, main="Test Data")

#########################################################################
# calculate likelihood value

likelihood <- function(betas, sigma) {
	y_pred <- X %*% betas 
	poster_y <- sum(dnorm(y, mean = y_pred, sd = sigma, log = TRUE))
	poster_beta <- sum(dnorm(betas, mean = 2, sd = 0.5, log = TRUE))
	poster_sigma <- dgamma(sigma, shape = 2, log = TRUE)
	#poster_sigma = dunif(sigma, min=0, max=30, log = T)
	return(poster_y + poster_beta + poster_sigma)
}

#########################################################################
# Proposal function

proposal <- function(betas, sigma) {
	new_betas <- rnorm(2, mean = betas, sd = 0.1)
	new_sigma <- rnorm(1, mean = sigma, sd = 0.5)
	return(c(new_betas, new_sigma))
}

#########################################################################
# Run Metropolis-Hastings

MH_estimate <- function(betas, sigma, maxiter = 1000) {
	chains <- array(dim = c(maxiter + 1, 3))
	chains[1, ] <- c(betas, sigma)
	
	for (ii in 1 : maxiter) {
		new_params <- proposal(betas, sigma)
		ratio <- exp(likelihood(new_params[1:2], new_params[3]) - likelihood(betas, sigma))
		print(ratio)
		if (runif(1) < min(ratio, 1)) {
			betas <- new_params[1 : 2]
			sigma <- new_params[3]
		} else {
			betas <- betas
			sigma <- sigma
		}		
		chains[ii + 1, ] <- new_params
	}
	return(chains)
}

#########################################################################
# Initial params and estimate.

sigma <- 1; betas <- c(0.1, 0.1)

chains <- MH_estimate(betas, sigma, maxiter = 10000)
