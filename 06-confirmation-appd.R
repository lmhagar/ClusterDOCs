## BEGIN SETUP ##

## load necessary packages
## run 04-simulations-appd.R first to get the .csv files used in this code file
## run 05-simulations-appd.R first to get the necessary functions
require(qrng)
require(foreach)
require(doParallel)
require(doSNOW)
require(rjags)
require(coda)
require(purrr)

## code to get the recommended cluster counts based on the linear approximations
gam <- 0.94
cs <- NULL
for (s in 1:3){
  ## consider clearly acceptable treatment
  j <- 1
  
  ## extract .csv files and construct linear approximations
  assign(paste0("res", s, j, "_25"), read.csv(paste0("scen", s, j, "_log_c_25.csv")))
  assign(paste0("res", s, j, "_45"), read.csv(paste0("scen", s, j, "_log_c_45.csv")))
  
  oc.temp <- lin_app(get(paste0("res", s,j,"_25")), get(paste0("res", s,j,"_45")), 
                         25, 45, 20, 50, 0.97)
  ## get the smallest cluster count where power is attained
  cs[s] <- oc.temp[which(oc.temp[,2] > 0.8)][1]
}

## get the operating characteristics based on linear approximations
res.lin <- NULL
## loop for ICC setting and scenario
for (s in 1:3){
  M <- cs[s]
  for (j in 1:4){
    
    ## extract posterior probabilities
    assign(paste0("res", s, j, "_25"), read.csv(paste0("scen", s, j, "_log_c_25.csv")))
    assign(paste0("res", s, j, "_45"), read.csv(paste0("scen", s, j, "_log_c_45.csv")))
    
    ## construct linear approximations and get operating characteristics for recommended c value
    res.lin <- c(res.lin, lin_app(get(paste0("res", s,j,"_25")), get(paste0("res", s,j,"_45")), 
                                          25, 45, 20, 50, gam)[M-19, 2])
  }
}

## format into table (rows for ICC setting, columnds for scenario)
res.lin <- matrix(res.lin, nrow = 3, byrow = TRUE)

## define these functions for later use
expit <- function(x){1/(1 + exp(-x))}
logit <- function(x){log(x) - log(1-x)}

## set up parallelization
cores=detectCores()
cl <- makeSOCKcluster(cores[1]-1)

## set up progress bar
m <- 10000
registerDoSNOW(cl)
pb <- txtProgressBar(max = m, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

## set up equivalence margin (delta_U in text)
d_ae <- 0.06

## start with the low ICC setting
## set up data generating coefficients for the four scenarios
beta0s <- c(-4.6798, -4.6798, -4.6798, -4.6798)
beta1s <- c(0.0000, 1.1228, 1.6581, 2.0195)
## set variance for the random intercepts
sig_u <- sqrt(pi^2/57)

## save MCMC settings for use later
n.chains = 1
n.burnin = 1000
n.draws = 4000
n.thin = 1

## get cluster count recommendation
Ms <- cs[1]
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
                         mclust <- rdunif(nclust, 10, 12)
                         
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
                                                            p0= 0.001, p1 = 0.01,
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
beta0s <- c(-4.8777, -4.8777, -4.8777, -4.8777)
beta1s <- c(0.0000, 1.1339, 1.6796, 2.0509)
## set variance for the random intercepts
sig_u <- sqrt(pi^2/17)

## get cluster count recommendation
Ms <- cs[2]
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
                         mclust <- rdunif(nclust, 10, 12)
                         
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
                                                            p0= 0.001, p1 = 0.01,
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
beta0s <- c(-5.1247, -5.1247, -5.1247, -5.1247)
beta1s <- c(0.0000, 1.1526, 1.7143, 2.0999)
## set variance for the random intercepts
sig_u <- sqrt(pi^2/9)

## get cluster count recommendation
Ms <- cs[3]
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
                         mclust <- rdunif(nclust, 10, 12)
                         
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
                                                            p0= 0.001, p1 = 0.01,
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

## get the estimated operating characteristics based on the confirmatory simulations
res.dat <- NULL
for (s in 1:3){
  for (j in 1:4){
    M <- cs[s]
    ## read in confirmatory estimates of sampling distribution (csv files)
    res_temp <- read.csv(paste0("scen",s,j,"_log_c_", M, ".csv"))
    ## compute operating characteristics using empirical averages
    res.dat <- c(res.dat, mean(1-res_temp[,1] > gam[1]))
  }
}

## format into table
res.dat <- matrix(res.dat, nrow = 3, byrow = TRUE)

## format into final tables
tab.app <- rbind(res.lin, res.dat)
tab.app <- tab.app[c(1, 4, 2, 5, 3, 6),]

## save results to .csv file
write.csv(tab.app, "table_app.csv", row.names = FALSE)