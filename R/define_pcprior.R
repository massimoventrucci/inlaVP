
# given 'u' and 'alpha', the function 'pcprior.interaction.lambda'
# compute the 'lambda' for the PC prior for mixing parameter f12 vs (f1+f2)
pcprior.interaction.lambda <- function(u, alpha)
{
  stopifnot(!missing(u) && !missing(alpha))
  alpha.min <- sqrt(u)
  if (!(alpha > alpha.min))
  {
    #stop(paste("inla.pc.cor1.b.lambda: alpha >", alpha.min))
    print(paste("alpha < alpha.min", alpha.min))
    fac <- 0.98
    alpha <- alpha.min * fac + 1 * (1-fac)
    print(paste("alpha set to ", alpha))
  }
  fun <- function(lam, u, alpha) {
    F <- (1 - exp(-lam * sqrt(u)))/(1 - exp(-lam))
    return((F - alpha)^2)
  }
  lambdas <- unique(c(seq(1e-07, 10, len = 100), seq(10, 100, len = 10)))
  idx <- which.min(fun(lambdas, u, alpha))
  stopifnot(idx > 1 && idx < length(lambdas))
  lambda <- optimize(fun, interval = lambdas[c(idx - 1, idx + 1)],
                    maximum = FALSE, u = u, alpha = alpha)$minimum
  stopifnot(fun(lambda, u, alpha) < 1e-08)
  return(lambda)
}

# # given lambda, compute the log density of the pc prior for 'sqrt(gamma)'
# # this was implemented in the augmented approach, not in the joint
# pc.sqrt.gamma <- function(sqrt.gamma, lambda, u, alpha, log = FALSE)
# {
#   if(missing(lambda))  lambda = pcprior.interaction.lambda(u,alpha)
#   log.dens = log(lambda) - lambda * sqrt.gamma - log(1 - exp(-lambda))
#   return(INLA:::inla.ifelse(log, log.dens, exp(log.dens)))
# }

# given lambda, compute the log density of the pc prior for 'gamma'
pc.gamma <- function(gamma, lambda, u, alpha, log = FALSE)
{
  if(missing(lambda))  lambda = pcprior.interaction.lambda(u,alpha)
  log.dens = log(lambda) - lambda * sqrt(gamma) - log(2*sqrt(gamma)) - log(1 - exp(-lambda))
  return(INLA:::inla.ifelse(log, log.dens, exp(log.dens)))
}



######### model (2) in the paper
# functions to map the 5-dim 'theta=log(prec)'
# into the hyperpar tau, gamma, phi, psi1, psi2 of model (2) in the paper (striid model)
# DONE and corrected on 12/07/21
# inputs: theta: 5-dim vector with log(c(tau1,tau2,tau3,tau4,tau5) values
#        tau1 = prec f(time.str)
#        tau2 = prec f(space.str)
#        tau3 = prec f(space.time interaction)
#        tau4 = prec f(time.iid)
#        tau5 = prec f(space.iid)
# returns: the hyperpar 'tau': total (generalized) precision
theta.to.tau.striid <- function(theta){
  theta1 <- theta[1]
  theta2 <- theta[2]
  theta3 <- theta[3]
  theta4 <- theta[4]
  theta5 <- theta[5]
  return( exp(theta1+theta2+theta3+theta4+theta5) / (exp(theta1+theta2+theta4+theta5)+
                                                     exp(theta2+theta3+theta5)*(exp(theta1)+exp(theta4))+
                                                     exp(theta1+theta3+theta4)*(exp(theta2)+exp(theta5))))
}

# inputs: theta: 5-dim vector with log(c(tau1,tau2,tau3,tau4,tau5) values
# returns: the hyperpar 'gamma': mixing int vs main
theta.to.gamma.striid <- function(theta){
  theta1 <- theta[1]
  theta2 <- theta[2]
  theta3 <- theta[3]
  theta4 <- theta[4]
  theta5 <- theta[5]
  return( exp(theta1+theta2+theta4+theta5) / (exp(theta1+theta2+theta4+theta5)+
                                                exp(theta2+theta3+theta5)*(exp(theta1)+exp(theta4))+
                                                exp(theta1+theta3+theta4)*(exp(theta2)+exp(theta5))))
}

# returns: the hyperpar 'phi': mixing time vs space
theta.to.phi.striid <- function(theta){
  theta1 <- theta[1]
  theta2 <- theta[2]
  # theta3 <- theta[3]
  theta4 <- theta[4]
  theta5 <- theta[5]
  return((exp(theta1+theta4)*(exp(theta2)+exp(theta5))) / (exp(theta1+theta4)*(exp(theta2)+exp(theta5)) + exp(theta2+theta5)*(exp(theta1)+exp(theta4))))
}

# returns: the hyperpar 'psi1': mixing str vs iid for time
theta.to.psi1.striid <- function(theta){
  theta1 <- theta[1]
  # theta2 <- theta[2]
  # theta3 <- theta[3]
  theta4 <- theta[4]
  # theta5 <- theta[5]
  return(exp(theta1) / (exp(theta1)+exp(theta4)))
}

# returns: the hyperpar 'psi2': mixing str vs iid for space
theta.to.psi2.striid <- function(theta){
  # theta1 <- theta[1]
  theta2 <- theta[2]
  # theta3 <- theta[3]
  # theta4 <- theta[4]
  theta5 <- theta[5]
  return(exp(theta2) / (exp(theta2)+exp(theta5)))
}


# inputs: taugammaphipsi1psi2 is a 5-dim vector with hyperpars values (tau, gamma, phi, psi1, psi2)
# returns: theta=log(prec)
taugammaphipsi1psi2.to.theta <- function(taugammaphipsi1psi2){
  tau <- taugammaphipsi1psi2[1]
  gama <- taugammaphipsi1psi2[2]
  phi <- taugammaphipsi1psi2[3]
  psi1 <- taugammaphipsi1psi2[4]
  psi2 <- taugammaphipsi1psi2[5]
  theta <- rep(NA, 5)
  theta[1] <- log(tau/((1-gama)*(1-phi)*(1-psi1)))
  theta[2] <- log(tau/((1-gama)*phi*(1-psi2)))
  theta[3] <- log(tau/gama)
  theta[4] <- log(tau/((1-gama)*(1-phi)*psi1))
  theta[5] <- log(tau/((1-gama)*phi*psi2))
  return(theta)
}



######### model (1) in the paper
# functions to map the 3-dim 'theta=log(prec)'
# into the hyperpar tau, gamma, phi of model (1) in the paper
theta.to.tau <- function(theta){
  theta1 <- theta[1]
  theta2 <- theta[2]
  theta3 <- theta[3]
  return(exp(theta1+theta2+theta3)/(exp(theta2+theta3) + exp(theta1+theta3) + exp(theta1+theta2)))
}

theta.to.gamma <- function(theta){
  theta1 <- theta[1]
  theta2 <- theta[2]
  theta3 <- theta[3]
  return(exp(theta1+theta2)/(exp(theta2+theta3) + exp(theta1+theta3) + exp(theta1+theta2)))
}

theta.to.phi <- function(theta){
  theta1 <- theta[1]
  theta2 <- theta[2]
  #theta3 <- theta[3]
  return(exp(theta1)/(exp(theta2) + exp(theta1)))
}

# compute jacobian, this goers inside the joint prior function
ljac <- function(theta){
  theta1 <- theta[1]
  theta2 <- theta[2]
  theta3 <- theta[3]
  detjac <- (exp(theta1)^3 * exp(theta2)^3 * exp(theta3)^2)/((exp(theta1) + exp(theta2))*(exp(theta2+theta3) + exp(theta1+theta3) + exp(theta1+theta2))^3)
  return(log(abs(detjac)))
}

# from tau, gamma, phi to the 3-dim theta
taugammaphi.to.theta <- function(taugammaphi){
  tau <- taugammaphi[1]
  gama <- taugammaphi[2]
  phi <- taugammaphi[3]
  theta <- rep(NA, 3)
  theta[1] <- log(tau/((1-gama)*(1-phi)))
  theta[2] <- log(tau/((1-gama)*phi))
  theta[3] <- log(tau/gama)
  return(theta)
}
