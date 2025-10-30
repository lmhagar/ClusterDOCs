## BEGIN SETUP ##

## load necessary packages
require(qrng)
require(foreach)
require(doParallel)
require(doSNOW)
require(rjags)
require(coda)
require(purrr)

## define these functions for later use
expit <- function(x){1/(1 + exp(-x))}
logit <- function(x){log(x) - log(1-x)}

## set up parallelization
cores=detectCores()
cl <- makeSOCKcluster(min(47, cores[1]-1))

## set up progress bar
m <- 10000
registerDoSNOW(cl)
pb <- txtProgressBar(max = m, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

## set up equivalence margin (delta_U in text)
d_ae <- 0.04

## set up array of cluster counts to consider
Ms <- seq(80, 160, 10)

## start with the low ICC setting
## set up data generating coefficients for the four scenarios
beta0s <- c(-3.9746, -3.9746, -3.9746, -3.9746)
beta1s <- c(0.0000, 0.4176, 0.9529, 1.1477)
## set variance for the random intercepts
sig_u <- sqrt(pi^2/57)

## save MCMC settings for use later
n.chains = 1
n.burnin = 500
n.draws = 2000
n.thin = 1

## loop over the four scenarios
for (ii in 1:4){
  ## loop over each cluster count considered
  for (j in 1:length(Ms)){
    ## extract revelant cluster count
    M <- Ms[j]
    ## get the sampling distribution estimate in parallel
    sim.res <- foreach(k=1:m, .packages=c('rjags', 'coda', 'purrr'), .combine=rbind,
                       .options.snow=opts) %dopar% {
                         
                         ## generate the cluster sizes
                         set.seed((ii-1)*j*m + (j-1)*m + k)
                         nclust <- M
                         mclust <- rdunif(nclust, 4, 6)
                         
                         ## extract the relevant coefficients
                         beta0 <- beta0s[ii]
                         beta1 <- beta1s[ii]
                      
                         dat <- data.frame(ae = 0, x = 0, clust = 0)
                         
                         ## for each cluster
                         for (i in 1:nclust){
                           
                           ## get the random intercept and treatment assignment
                           ae.b0 <- rnorm(1, 0, sig_u)
                           x.temp <- rep(runif(1)< 1/2, mclust[i])
                           
                           ## get the linear predictor and binary data
                           ae.eta <- expit(beta0 + ae.b0 + beta1*x.temp)
                           ae.temp <- rbinom(mclust[i], 1, ae.eta)
                           
                           ## save data
                           dat <- rbind(dat,
                                        data.frame(ae = ae.temp, x = x.temp,
                                                   clust = rep(i, mclust[i])))
                         }
                         dat <- dat[-1,]
                         
                         ## fit model using JAGS
                         log.int.wd <- paste(getwd(), '/jags_logistic_int.txt', sep='')
                         model1.fit <- jags.model(file=log.int.wd,
                                                  data=list(N=nrow(dat), Y=dat$ae, clust = dat$clust,
                                                            nclust = nclust, X = dat$x,
                                                            p0= 0.0001, p1 = 0.001,
                                                            sig.up = 25),
                                                  n.chains = n.chains)
                         
                         update(model1.fit, n.burnin)
                         
                         model1.samples.ae <- coda.samples(model1.fit, c("beta0", "beta1", "sigma_int", "b0"), 
                                                        n.iter=n.draws, thin=n.thin)
                         
                         ## extract posterior draws
                         beta0.post.ae <- unlist(model1.samples.ae[,ncol(model1.samples.ae[[1]]) - 2])
                         beta1.post.ae <- unlist(model1.samples.ae[,ncol(model1.samples.ae[[1]]) - 1])
                         
                         save0.ae <- NULL
                         save1.ae <- NULL
                         mm <- nrow(model1.samples.ae[[1]])
                         for (i in 1:mm){
                           ## get Dirichlet weights for treatment arm
                           mu1.ae <- expit(beta0.post.ae[i] + beta1.post.ae[i] + model1.samples.ae[[1]][i, 1:nclust])
                           w1 <- rexp(nclust, 1)
                           w1 <- w1/sum(w1)
                           
                           ## get Dirichlet weights for control arm
                           mu0.ae <- expit(beta0.post.ae[i] + model1.samples.ae[[1]][i, 1:nclust])
                           w0 <- rexp(nclust, 1)
                           w0 <- w0/sum(w0)
                           
                           ## sum to get the marginal estimands in each group
                           save1.ae[i] <- sum(w1*mu1.ae) 
                           save0.ae[i] <- sum(w0*mu0.ae)
                         }
                         ## get the difference in marginal estimands
                         ldiff.ae <- save1.ae - save0.ae
                         
                         ## get the posterior probability using kernel density estimation
                         kd.ae <- density(na.omit(ldiff.ae))
                         np.prob.ae <- mean(pnorm(d_ae, ifelse(is.finite(ldiff.ae), ldiff.ae, 
                                                               max(na.omit(ldiff.ae)) + 1), kd.ae$bw, lower.tail = FALSE))
                         
                         ## use normal approximation to posterior if probability is really small (for stability)
                         if (np.prob.ae < 0.00004){
                           np.prob.ae <- pnorm(d_ae, mean(na.omit(ldiff.ae)), sd(na.omit(ldiff.ae)), lower.tail = FALSE)
                         }
      
                         ## output posterior probability
                         np.prob.ae
                       }
    
    ## save results for each scenario
    write.csv(sim.res, paste0("scen1",ii,"_log_c_", M, ".csv"), row.names = FALSE)
  }
}

## continue with the moderate ICC setting
## set up data generating coefficients for the four scenarios
beta0s <- c(-4.1669, -4.1669, -4.1669, -4.1669)
beta1s <- c(0.0000, 0.4231, 0.9688, 1.1686)
## set variance for the random intercepts
sig_u <- sqrt(pi^2/17)

## loop over the four scenarios
for (ii in 1:4){
  ## loop over each cluster count considered
  for (j in 1:length(Ms)){
    ## extract revelant cluster count
    M <- Ms[j]
    ## get the sampling distribution estimate in parallel
    sim.res <- foreach(k=1:m, .packages=c('rjags', 'coda', 'purrr'), .combine=rbind,
                       .options.snow=opts) %dopar% {
                         
                         ## generate the cluster sizes
                         set.seed(440000 + (ii-1)*j*m + (j-1)*m + k)
                         nclust <- M
                         mclust <- rdunif(nclust, 4, 6)
                         
                         ## extract the relevant coefficients
                         beta0 <- beta0s[ii]
                         beta1 <- beta1s[ii]
                         
                         dat <- data.frame(ae = 0, x = 0, clust = 0)
                         
                         ## for each cluster
                         for (i in 1:nclust){
                           
                           ## get the random intercept and treatment assignment
                           ae.b0 <- rnorm(1, 0, sig_u)
                           x.temp <- rep(runif(1)< 1/2, mclust[i])
                           
                           ## get the linear predictor and binary data
                           ae.eta <- expit(beta0 + ae.b0 + beta1*x.temp)
                           ae.temp <- rbinom(mclust[i], 1, ae.eta)
                           
                           ## save data
                           dat <- rbind(dat,
                                        data.frame(ae = ae.temp, x = x.temp,
                                                   clust = rep(i, mclust[i])))
                         }
                         dat <- dat[-1,]
                         
                         ## fit model using JAGS
                         log.int.wd <- paste(getwd(), '/jags_logistic_int.txt', sep='')
                         model1.fit <- jags.model(file=log.int.wd,
                                                  data=list(N=nrow(dat), Y=dat$ae, clust = dat$clust,
                                                            nclust = nclust, X = dat$x,
                                                            p0= 0.0001, p1 = 0.001,
                                                            sig.up = 25),
                                                  n.chains = n.chains)
                         
                         update(model1.fit, n.burnin)
                         
                         model1.samples.ae <- coda.samples(model1.fit, c("beta0", "beta1", "sigma_int", "b0"), 
                                                           n.iter=n.draws, thin=n.thin)
                         
                         ## extract posterior draws
                         beta0.post.ae <- unlist(model1.samples.ae[,ncol(model1.samples.ae[[1]]) - 2])
                         beta1.post.ae <- unlist(model1.samples.ae[,ncol(model1.samples.ae[[1]]) - 1])
                         
                         save0.ae <- NULL
                         save1.ae <- NULL
                         mm <- nrow(model1.samples.ae[[1]])
                         for (i in 1:mm){
                           ## get Dirichlet weights for treatment arm
                           mu1.ae <- expit(beta0.post.ae[i] + beta1.post.ae[i] + model1.samples.ae[[1]][i, 1:nclust])
                           w1 <- rexp(nclust, 1)
                           w1 <- w1/sum(w1)
                           
                           ## get Dirichlet weights for control arm
                           mu0.ae <- expit(beta0.post.ae[i] + model1.samples.ae[[1]][i, 1:nclust])
                           w0 <- rexp(nclust, 1)
                           w0 <- w0/sum(w0)
                           
                           ## sum to get the marginal estimands in each group
                           save1.ae[i] <- sum(w1*mu1.ae) 
                           save0.ae[i] <- sum(w0*mu0.ae)
                         }
                         ## get the difference in marginal estimands
                         ldiff.ae <- save1.ae - save0.ae
                         
                         ## get the posterior probability using kernel density estimation
                         kd.ae <- density(na.omit(ldiff.ae))
                         np.prob.ae <- mean(pnorm(d_ae, ifelse(is.finite(ldiff.ae), ldiff.ae, 
                                                               max(na.omit(ldiff.ae)) + 1), kd.ae$bw, lower.tail = FALSE))
                         
                         ## use normal approximation to posterior if probability is really small (for stability)
                         if (np.prob.ae < 0.00004){
                           np.prob.ae <- pnorm(d_ae, mean(na.omit(ldiff.ae)), sd(na.omit(ldiff.ae)), lower.tail = FALSE)
                         }
                         
                         ## output posterior probability
                         np.prob.ae
                       }
    
    ## save results for each scenario
    write.csv(sim.res, paste0("scen2",ii,"_log_c_", M, ".csv"), row.names = FALSE)
  }
}

## continue with the high ICC setting
## set up data generating coefficients for the four scenarios
beta0s <- c(-4.4041, -4.4041, -4.4041, -4.4041)
beta1s <- c(0.0000, 0.4330, 0.9937, 1.2008)
## set variance for the random intercepts
sig_u <- sqrt(pi^2/9)

## loop over the four scenarios
for (ii in 1:4){
  ## loop over each cluster count considered
  for (j in 1:length(Ms)){
    ## extract revelant cluster count
    M <- Ms[j]
    ## get the sampling distribution estimate in parallel
    sim.res <- foreach(k=1:m, .packages=c('rjags', 'coda', 'purrr'), .combine=rbind,
                       .options.snow=opts) %dopar% {
                         
                         ## generate the cluster sizes
                         set.seed(880000 + (ii-1)*j*m + (j-1)*m + k)
                         nclust <- M
                         mclust <- rdunif(nclust, 4, 6)
                         
                         ## extract the relevant coefficients
                         beta0 <- beta0s[ii]
                         beta1 <- beta1s[ii]
                         
                         dat <- data.frame(ae = 0, x = 0, clust = 0)
                         
                         ## for each cluster
                         for (i in 1:nclust){
                           
                           ## get the random intercept and treatment assignment
                           ae.b0 <- rnorm(1, 0, sig_u)
                           x.temp <- rep(runif(1)< 1/2, mclust[i])
                           
                           ## get the linear predictor and binary data
                           ae.eta <- expit(beta0 + ae.b0 + beta1*x.temp)
                           ae.temp <- rbinom(mclust[i], 1, ae.eta)
                           
                           ## save data
                           dat <- rbind(dat,
                                        data.frame(ae = ae.temp, x = x.temp,
                                                   clust = rep(i, mclust[i])))
                         }
                         dat <- dat[-1,]
                         
                         ## fit model using JAGS
                         log.int.wd <- paste(getwd(), '/jags_logistic_int.txt', sep='')
                         model1.fit <- jags.model(file=log.int.wd,
                                                  data=list(N=nrow(dat), Y=dat$ae, clust = dat$clust,
                                                            nclust = nclust, X = dat$x,
                                                            p0= 0.0001, p1 = 0.001,
                                                            sig.up = 25),
                                                  n.chains = n.chains)
                         
                         update(model1.fit, n.burnin)
                         
                         model1.samples.ae <- coda.samples(model1.fit, c("beta0", "beta1", "sigma_int", "b0"), 
                                                           n.iter=n.draws, thin=n.thin)
                         
                         ## extract posterior draws
                         beta0.post.ae <- unlist(model1.samples.ae[,ncol(model1.samples.ae[[1]]) - 2])
                         beta1.post.ae <- unlist(model1.samples.ae[,ncol(model1.samples.ae[[1]]) - 1])
                         
                         save0.ae <- NULL
                         save1.ae <- NULL
                         mm <- nrow(model1.samples.ae[[1]])
                         for (i in 1:mm){
                           ## get Dirichlet weights for treatment arm
                           mu1.ae <- expit(beta0.post.ae[i] + beta1.post.ae[i] + model1.samples.ae[[1]][i, 1:nclust])
                           w1 <- rexp(nclust, 1)
                           w1 <- w1/sum(w1)
                           
                           ## get Dirichlet weights for control arm
                           mu0.ae <- expit(beta0.post.ae[i] + model1.samples.ae[[1]][i, 1:nclust])
                           w0 <- rexp(nclust, 1)
                           w0 <- w0/sum(w0)
                           
                           ## sum to get the marginal estimands in each group
                           save1.ae[i] <- sum(w1*mu1.ae) 
                           save0.ae[i] <- sum(w0*mu0.ae)
                         }
                         ## get the difference in marginal estimands
                         ldiff.ae <- save1.ae - save0.ae
                         
                         ## get the posterior probability using kernel density estimation
                         kd.ae <- density(na.omit(ldiff.ae))
                         np.prob.ae <- mean(pnorm(d_ae, ifelse(is.finite(ldiff.ae), ldiff.ae, 
                                                               max(na.omit(ldiff.ae)) + 1), kd.ae$bw, lower.tail = FALSE))
                         
                         ## use normal approximation to posterior if probability is really small (for stability)
                         if (np.prob.ae < 0.00004){
                           np.prob.ae <- pnorm(d_ae, mean(na.omit(ldiff.ae)), sd(na.omit(ldiff.ae)), lower.tail = FALSE)
                         }
                         
                         ## output posterior probability
                         np.prob.ae
                       }
    
    ## save results for each scenario
    write.csv(sim.res, paste0("scen3",ii,"_log_c_", M, ".csv"), row.names = FALSE)
  }
}