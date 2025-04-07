## BEGIN SETUP ##

## load necessary packages
## run 01-simulations.R first to get the .csv files used in this code file
require(qrng)
require(foreach)
require(doParallel)
require(doSNOW)
require(rjags)
require(coda)
require(purrr)

## functions to convert results from .csv files to logits
## of posterior probabilities
logit <- function(x){log(x) - log(1-x)}
## try to get "power" curves
l1 <- function(x){
  -1*logit(x)
}

## function to obtain operating characteristics for the interim analysis
lin_app <- function(ml, mu, cl, cu, lb, ub, gam){
  # ml is lower matrix
  # mu is upper matrix
  # cl is lower c
  # cu is upper c
  # lb is lower bound
  # ub is upper bound
  # gam are the thresholds
  
  # get logits for the smaller and larger sample sizes
  ll <- apply(ml, 2, l1)
  lu <- apply(mu, 2, l1)
  
  # adjust an infinite logits
  for (j in 1:ncol(ll)){
    ll[,j] <- ifelse(ll[,j] == -Inf, min(subset(ll, is.finite(ll[,j]))[,j]) - 1, ll[,j])
    ll[,j] <- ifelse(ll[,j] == Inf, max(subset(ll, is.finite(ll[,j]))[,j]) + 1, ll[,j])
  }
  
  for (j in 1:ncol(lu)){
    lu[,j] <- ifelse(lu[,j] == -Inf, min(subset(lu, is.finite(lu[,j]))[,j]) - 1, lu[,j])
    lu[,j] <- ifelse(lu[,j] == Inf, max(subset(lu, is.finite(lu[,j]))[,j]) + 1, lu[,j])
  }
  
  # get indexes to combine marginal logits later
  ll <- cbind(ll, seq(1, nrow(ll), 1))
  lu <- cbind(lu, seq(1, nrow(lu), 1))
  
  # process for AE event
  # sort logits
  ll_s <- ll[order(ll[,1]),1]
  lu_s <- lu[order(lu[,1]),1]
  
  # get slopes and intercepts of linear approximations
  l_slope <- (lu_s - ll_s)/(cu-cl)
  l_int <- ll_s - l_slope*cl
  
  # reorder according to smaller sample size
  l_slope[ll[order(ll[,1]),4]] <- l_slope 
  l_int[ll[order(ll[,1]),4]] <- l_int 
  
  # save results
  slopes <- NULL
  slopes <- cbind(slopes, l_slope)
  
  ints <- NULL
  ints <- cbind(ints, l_int)
  
  # repeat process for completion outcome
  ll_s <- ll[order(ll[,2]),2]
  lu_s <- lu[order(lu[,2]),2]
  
  l_slope <- (lu_s - ll_s)/(cu-cl)
  l_int <- ll_s - l_slope*cl
  
  l_slope[ll[order(ll[,2]),4]] <- l_slope 
  l_int[ll[order(ll[,2]),4]] <- l_int 
  
  slopes <- cbind(slopes, l_slope)
  ints <- cbind(ints, l_int)
  
  # repeat process for tolerability outcome
  ll_s <- ll[order(ll[,3]),3]
  lu_s <- lu[order(lu[,3]),3]
  
  l_slope <- (lu_s - ll_s)/(cu-cl)
  l_int <- ll_s - l_slope*cl
  
  l_slope[ll[order(ll[,3]),4]] <- l_slope 
  l_int[ll[order(ll[,3]),4]] <- l_int 
  
  slopes <- cbind(slopes, l_slope)
  ints <- cbind(ints, l_int)
  
  # get matrices of new posterior probabilities corresponding to
  # linear approximations
  cs <- seq(lb, ub,1)
  lae <- matrix(0, nrow = nrow(slopes), ncol=length(cs))
  lcomp <- matrix(0, nrow = nrow(slopes), ncol=length(cs))
  ltol <- matrix(0, nrow = nrow(slopes), ncol=length(cs))
  for (i in 1:length(cs)){
    lae[,i] <- ints[,1] + slopes[,1]*cs[i]
    lcomp[,i] <- ints[,2] + slopes[,2]*cs[i]
    ltol[,i] <- ints[,3] + slopes[,3]*cs[i]
  }
  
  # create indicator nu to determine if all criteria are satisified
  ylae <- lae > logit(gam[1])
  ylcomp <- lcomp > logit(gam[2])
  yltol <- ltol > logit(gam[3])
  
  all3 <- ylae + ylcomp + yltol
  
  return(cbind(cs, colMeans(all3 >= 3)))
  
}

## function to obtain operating characteristics for the final analysis
lin_app_seq <- function(ml, mu, cl, cu, lb, ub, gam){
  # ml is lower matrix
  # mu is upper matrix
  # cl is lower c
  # cu is upper c
  # lb is lower bound
  # ub is upper bound
  # gam are the thresholds
  
  # get logits for the smaller and larger sample sizes
  ll <- apply(ml, 2, l1)
  lu <- apply(mu, 2, l1)
  
  # adjust an infinite logits
  for (j in 1:ncol(ll)){
    ll[,j] <- ifelse(ll[,j] == -Inf, min(subset(ll, is.finite(ll[,j]))[,j]) - 1, ll[,j])
    ll[,j] <- ifelse(ll[,j] == Inf, max(subset(ll, is.finite(ll[,j]))[,j]) + 1, ll[,j])
  }
  
  for (j in 1:ncol(lu)){
    lu[,j] <- ifelse(lu[,j] == -Inf, min(subset(lu, is.finite(lu[,j]))[,j]) - 1, lu[,j])
    lu[,j] <- ifelse(lu[,j] == Inf, max(subset(lu, is.finite(lu[,j]))[,j]) + 1, lu[,j])
  }
  
  # get indexes to combine marginal logits later
  ll <- cbind(ll, seq(1, nrow(ll), 1))
  lu <- cbind(lu, seq(1, nrow(lu), 1))
  
  # process for AE event
  # sort logits
  ll_s <- ll[order(ll[,1]),1]
  lu_s <- lu[order(lu[,1]),1]
  
  # get slopes and intercepts of linear approximations
  l_slope <- (lu_s - ll_s)/(cu-cl)
  l_int <- ll_s - l_slope*cl
  
  # reorder according to smaller sample size
  l_slope[ll[order(ll[,1]),5]] <- l_slope 
  l_int[ll[order(ll[,1]),5]] <- l_int 
  
  # save results
  slopes <- NULL
  slopes <- cbind(slopes, l_slope)
  
  ints <- NULL
  ints <- cbind(ints, l_int)
  
  # repeat process for completion outcome
  ll_s <- ll[order(ll[,2]),2]
  lu_s <- lu[order(lu[,2]),2]
  
  l_slope <- (lu_s - ll_s)/(cu-cl)
  l_int <- ll_s - l_slope*cl
  
  l_slope[ll[order(ll[,2]),5]] <- l_slope 
  l_int[ll[order(ll[,2]),5]] <- l_int 
  
  slopes <- cbind(slopes, l_slope)
  ints <- cbind(ints, l_int)
  
  # repeat process for tolerability outcome
  ll_s <- ll[order(ll[,3]),3]
  lu_s <- lu[order(lu[,3]),3]
  
  l_slope <- (lu_s - ll_s)/(cu-cl)
  l_int <- ll_s - l_slope*cl
  
  l_slope[ll[order(ll[,3]),5]] <- l_slope 
  l_int[ll[order(ll[,3]),5]] <- l_int 
  
  slopes <- cbind(slopes, l_slope)
  ints <- cbind(ints, l_int)
  
  # repeat process for second AE analysis
  ll_s <- ll[order(ll[,4]),4]
  lu_s <- lu[order(lu[,4]),4]
  
  cuu <- 1.25*cu
  cll <- 1.25*cl
  
  l_slope <- (lu_s - ll_s)/(cuu-cll)
  l_int <- ll_s - l_slope*cll
  
  l_slope[ll[order(ll[,4]),5]] <- l_slope 
  l_int[ll[order(ll[,4]),5]] <- l_int 
  
  slopes <- cbind(slopes, l_slope)
  ints <- cbind(ints, l_int)
  
  # get matrices of new posterior probabilities corresponding to
  # linear approximations
  cs <- seq(lb, ub,1)
  lae <- matrix(0, nrow = nrow(slopes), ncol=length(cs))
  lcomp <- matrix(0, nrow = nrow(slopes), ncol=length(cs))
  ltol <- matrix(0, nrow = nrow(slopes), ncol=length(cs))
  lae2 <- matrix(0, nrow = nrow(slopes), ncol=length(cs))
  for (i in 1:length(cs)){
    lae[,i] <- ints[,1] + slopes[,1]*cs[i]
    lcomp[,i] <- ints[,2] + slopes[,2]*cs[i]
    ltol[,i] <- ints[,3] + slopes[,3]*cs[i]
    lae2[,i] <- ints[,4] + slopes[,4]*1.25*cs[i]
  }
  
  # create indicator nu to determine if all criteria are satisified
  ylae <- lae > logit(gam[1])
  ylcomp <- lcomp > logit(gam[2])
  yltol <- ltol > logit(gam[3])
  ylae2 <- lae2 > logit(gam[4])
  
  all4 <- ylae + ylcomp + yltol + ylae2
  
  return(cbind(cs, colMeans(all4 >= 4)))
  
}

## code to get the recommended sample sizes based on the linear approximations
gam <- c(0.8, 0.55, 0.55)
cs <- NULL
for (s in 1:3){
  ## consider acceptable treatment
  j <- 2
  assign(paste0("res", s, j, "_120"), read.csv(paste0("scen", s, j, "_log_c_120.csv")))
  assign(paste0("res", s, j, "_180"), read.csv(paste0("scen", s, j, "_log_c_180.csv")))
  
  oc.temp <- lin_app_seq(get(paste0("res", s,j,"_120")), get(paste0("res", s,j,"_180")), 
                         120, 180, 100, 200, c(gam, 0.94))
  cs[s] <- oc.temp[which(oc.temp[,2] > 0.8)][1]
}

## get the operating characteristics for the interim analysis based on linear approximations
res.int.lin <- NULL
for (s in 1:3){
  M <- cs[s]
  for (j in 1:4){
    
    assign(paste0("res", s, j, "_120"), read.csv(paste0("scen", s, j, "_log_c_120.csv"))[,1:3])
    assign(paste0("res", s, j, "_180"), read.csv(paste0("scen", s, j, "_log_c_180.csv"))[,1:3])
    
    res.int.lin <- c(res.int.lin, lin_app(get(paste0("res", s,j,"_120")), get(paste0("res", s,j,"_180")), 
                            120, 180, 100, 200, gam)[M-99, 2])
  }
}

## format into table
res.int.lin <- matrix(res.int.lin, nrow = 3, byrow = TRUE)

## get the operating characteristics for the final analysis based on linear approximations
res.fin.lin <- NULL
for (s in 1:3){
  M <- cs[s]
  for (j in 1:4){
    
    assign(paste0("res", s, j, "_120"), read.csv(paste0("scen", s, j, "_log_c_120.csv")))
    assign(paste0("res", s, j, "_180"), read.csv(paste0("scen", s, j, "_log_c_180.csv")))
    
    res.fin.lin <- c(res.fin.lin, lin_app_seq(get(paste0("res", s,j,"_120")), get(paste0("res", s,j,"_180")), 
                                          120, 180, 100, 200, c(gam, 0.94))[M-99, 2])
  }
}

## format into table
res.fin.lin <- matrix(res.fin.lin, nrow = 3, byrow = TRUE)

## define these functions for later use
expit <- function(x){1/(1 + exp(-x))}
logit <- function(x){log(x) - log(1-x)}

## set up parallelization
cores=detectCores()
cl <- makeSOCKcluster(cores[1]-1)

m <- 10000
registerDoSNOW(cl)
pb <- txtProgressBar(max = m, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

## define non-inferiority margins
d_ae <- 0.04
d_comp <- 0.1
d_ntol <- 0.1
Ms <- cs[1]
Mfs <- ceiling(1.25*Ms)

## parameters for low ICC setting
beta0s <- c(-3.9746, -3.9746, -3.9746, -3.9746)
beta1s <- c(0.0000, 0.4176, 0.9529, 1.1477)
phi0s <- c(1.4916, 1.4916, 1.4916, 1.4916)
phi1s <- c(0.0000, -0.1550, -0.1891, -0.4675)
for (ii in 1:4){
  for (j in 1:length(Ms)){
    M <- Ms[j]
    Mf <- Mfs[j]
    
    sim.res <- foreach(k=1:m, .packages=c('rjags', 'coda', 'purrr'), .combine=rbind,
                       .options.snow=opts) %dopar% {
                         
                         set.seed((ii-1)*j*m + (j-1)*m + k)
                         nclust <- M
                         mclust <- rdunif(nclust, 4, 6)
                         
                         beta0 <- beta0s[ii]
                         beta1 <- beta1s[ii]
                         phi0 <- phi0s[ii]
                         phi1 <- phi1s[ii]
                         
                         sig_u <- sqrt(pi^2/57)
                         sig_v <- sqrt(pi^2/7)
                         rho <- 0.9
                         
                         dat <- data.frame(ae = 0, comp = 0, ntol = 0, 
                                           x = 0, clust = 0)
                         for (i in 1:nclust){
                           ae.b0 <- rnorm(1, 0, sig_u)
                           x.temp <- rep(runif(1)< 1/2, mclust[i])
                           ae.eta <- expit(beta0 + ae.b0 + beta1*x.temp)
                           ae.temp <- rbinom(mclust[i], 1, ae.eta)
                           
                           comp.b0 <- rnorm(1, 0, sig_v)
                           comp.eta <- expit(phi0 + comp.b0 + phi1*x.temp)
                           comp.temp <- rbinom(mclust[i], 1, comp.eta)
                           comp.temp <- ifelse(ae.temp == 1, 0, comp.temp)
                           
                           tol.b0 <- rnorm(1, rho*comp.b0, sqrt(1-rho^2)*sig_v)
                           tol.eta <- expit(phi0 + tol.b0 + phi1*x.temp)
                           tol.temp <- rbinom(mclust[i], 1, tol.eta)
                           ntol.temp <- ifelse(ae.temp == 1, 1, 1 - tol.temp)
                           
                           dat <- rbind(dat,
                                        data.frame(ae = ae.temp, comp = comp.temp,
                                                   ntol = ntol.temp, x = x.temp,
                                                   clust = rep(i, mclust[i])))
                         }
                         dat <- dat[-1,]
                         
                         n.chains = 1
                         n.burnin = 500
                         n.draws = 2000
                         n.thin = 1
                         
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
                         
                         beta0.post.ae <- unlist(model1.samples.ae[,ncol(model1.samples.ae[[1]]) - 2])
                         beta1.post.ae <- unlist(model1.samples.ae[,ncol(model1.samples.ae[[1]]) - 1])
                         
                         ## repeat for completion
                         model1.fit <- jags.model(file=log.int.wd,
                                                  data=list(N=nrow(dat), Y=dat$comp, clust = dat$clust,
                                                            nclust = nclust, X = dat$x,
                                                            p0= 0.0001, p1 = 0.001,
                                                            sig.up = 25),
                                                  n.chains = n.chains)
                         
                         update(model1.fit, n.burnin)
                         
                         model1.samples.comp <- coda.samples(model1.fit, c("beta0", "beta1", "sigma_int", "b0"), 
                                                             n.iter=n.draws, thin=n.thin)
                         
                         beta0.post.comp <- unlist(model1.samples.comp[,ncol(model1.samples.comp[[1]]) - 2])
                         beta1.post.comp <- unlist(model1.samples.comp[,ncol(model1.samples.comp[[1]]) - 1])
                         
                         ## repeat for non-tolerability
                         model1.fit <- jags.model(file=log.int.wd,
                                                  data=list(N=nrow(dat), Y=dat$ntol, clust = dat$clust,
                                                            nclust = nclust, X = dat$x,
                                                            p0= 0.0001, p1 = 0.001,
                                                            sig.up = 25),
                                                  n.chains = n.chains)
                         
                         update(model1.fit, n.burnin)
                         
                         model1.samples.tol <- coda.samples(model1.fit, c("beta0", "beta1", "sigma_int", "b0"), 
                                                            n.iter=n.draws, thin=n.thin)
                         
                         beta0.post.tol <- unlist(model1.samples.tol[,ncol(model1.samples.tol[[1]]) - 2])
                         beta1.post.tol <- unlist(model1.samples.tol[,ncol(model1.samples.tol[[1]]) - 1])
                         
                         save0.ae <- NULL
                         save1.ae <- NULL
                         save0.comp <- NULL
                         save1.comp <- NULL
                         save0.tol <- NULL
                         save1.tol <- NULL
                         mm <- nrow(model1.samples.ae[[1]])
                         for (i in 1:mm){
                           mu1.ae <- expit(beta0.post.ae[i] + beta1.post.ae[i] + model1.samples.ae[[1]][i, 1:nclust])
                           mu1.comp <- expit(beta0.post.comp[i] + beta1.post.comp[i] + model1.samples.comp[[1]][i, 1:nclust])
                           mu1.tol <- expit(beta0.post.tol[i] + beta1.post.tol[i] + model1.samples.tol[[1]][i, 1:nclust])
                           w1 <- rexp(nclust, 1)
                           w1 <- w1/sum(w1)
                           
                           mu0.ae <- expit(beta0.post.ae[i] + model1.samples.ae[[1]][i, 1:nclust])
                           mu0.comp <- expit(beta0.post.comp[i] + model1.samples.comp[[1]][i, 1:nclust])
                           mu0.tol <- expit(beta0.post.tol[i] + model1.samples.tol[[1]][i, 1:nclust])
                           w0 <- rexp(nclust, 1)
                           w0 <- w0/sum(w0)
                           
                           save1.ae[i] <- sum(w1*mu1.ae) 
                           save0.ae[i] <- sum(w0*mu0.ae)
                           save1.comp[i] <- sum(w1*mu1.comp) 
                           save0.comp[i] <- sum(w0*mu0.comp)
                           save1.tol[i] <- sum(w1*mu1.tol) 
                           save0.tol[i] <- sum(w0*mu0.tol)
                         }
                         ldiff.ae <- save1.ae - save0.ae
                         
                         kd.ae <- density(na.omit(ldiff.ae))
                         np.prob.ae <- mean(pnorm(d_ae, ifelse(is.finite(ldiff.ae), ldiff.ae, 
                                                               max(na.omit(ldiff.ae)) + 1), kd.ae$bw, lower.tail = FALSE))
                         
                         if (np.prob.ae < 0.00004){
                           np.prob.ae <- pnorm(d_ae, mean(na.omit(ldiff.ae)), sd(na.omit(ldiff.ae)), lower.tail = FALSE)
                         }
                         
                         ldiff.comp <- save0.comp - save1.comp
                         
                         kd.comp <- density(na.omit(ldiff.comp))
                         np.prob.comp <- mean(pnorm(d_comp, ifelse(is.finite(ldiff.comp), ldiff.comp, 
                                                                   max(na.omit(ldiff.comp)) + 1), kd.comp$bw, lower.tail = FALSE))
                         
                         if (np.prob.comp < 0.00004){
                           np.prob.comp <- pnorm(d_comp, mean(na.omit(ldiff.comp)), sd(na.omit(ldiff.comp)), lower.tail = FALSE)
                         }
                         
                         ldiff.tol <- save1.tol - save0.tol
                         
                         kd.tol <- density(na.omit(ldiff.tol))
                         np.prob.tol <- mean(pnorm(d_ntol, ifelse(is.finite(ldiff.tol), ldiff.tol, 
                                                                  max(na.omit(ldiff.tol)) + 1), kd.tol$bw, lower.tail = FALSE))
                         
                         if (np.prob.tol < 0.00004){
                           np.prob.tol <- pnorm(d_ntol, mean(na.omit(ldiff.tol)), sd(na.omit(ldiff.tol)), lower.tail = FALSE)
                         }
                         
                         dat <- dat[, c(-2, -3)]
                         ## now get the probability for the final analysis
                         nclust <- Mf
                         mclust <- c(mclust, rdunif(Mf - M, 4, 6))
                         
                         for (i in (M+1):Mf){
                           ae.b0 <- rnorm(1, 0, sig_u)
                           x.temp <- rep(runif(1)< 1/2, mclust[i])
                           ae.eta <- expit(beta0 + ae.b0 + beta1*x.temp)
                           ae.temp <- rbinom(mclust[i], 1, ae.eta)
                           
                           dat <- rbind(dat,
                                        data.frame(ae = ae.temp, x = x.temp,
                                                   clust = rep(i, mclust[i])))
                         }
                         
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
                         
                         beta0.post.ae <- unlist(model1.samples.ae[,ncol(model1.samples.ae[[1]]) - 2])
                         beta1.post.ae <- unlist(model1.samples.ae[,ncol(model1.samples.ae[[1]]) - 1])
                         
                         save0.ae <- NULL
                         save1.ae <- NULL
                         mm <- nrow(model1.samples.ae[[1]])
                         for (i in 1:mm){
                           mu1.ae <- expit(beta0.post.ae[i] + beta1.post.ae[i] + model1.samples.ae[[1]][i, 1:nclust])
                           w1 <- rexp(nclust, 1)
                           w1 <- w1/sum(w1)
                           
                           mu0.ae <- expit(beta0.post.ae[i] + model1.samples.ae[[1]][i, 1:nclust])
                           w0 <- rexp(nclust, 1)
                           w0 <- w0/sum(w0)
                           
                           save1.ae[i] <- sum(w1*mu1.ae) 
                           save0.ae[i] <- sum(w0*mu0.ae)
                         }
                         ldiff.ae2 <- save1.ae - save0.ae
                         
                         kd.ae2 <- density(na.omit(ldiff.ae2))
                         np.prob.ae2 <- mean(pnorm(d_ae, ifelse(is.finite(ldiff.ae2), ldiff.ae2, 
                                                                max(na.omit(ldiff.ae2)) + 1), kd.ae2$bw, lower.tail = FALSE))
                         
                         if (np.prob.ae2 < 0.00004){
                           np.prob.ae2 <- pnorm(d_ae, mean(na.omit(ldiff.ae2)), sd(na.omit(ldiff.ae2)), lower.tail = FALSE)
                         } 
                         
                         c(np.prob.ae, np.prob.comp, np.prob.tol, np.prob.ae2)
                       }
    write.csv(sim.res, paste0("scen1",ii,"_log_c_", M, ".csv"), row.names = FALSE)
  }
}

## parameters for moderate ICC setting
Ms <- cs[2]
Mfs <- ceiling(1.25*Ms)

beta0s <- c(-4.1669, -4.1669, -4.1669, -4.1669)
beta1s <- c(0.0000, 0.4231, 0.9688, 1.1686)
phi0s <- c(1.6322, 1.6322, 1.6322, 1.6322)
phi1s <- c(0.0000, -0.1658, -0.2048, -0.5105)
for (ii in 1:4){
  for (j in 1:length(Ms)){
    M <- Ms[j]
    Mf <- Mfs[j]
    sim.res <- foreach(k=1:m, .packages=c('rjags', 'coda', 'purrr'), .combine=rbind,
                       .options.snow=opts) %dopar% {
                         
                         set.seed(440000 + (ii-1)*j*m + (j-1)*m + k)
                         nclust <- M
                         mclust <- rdunif(nclust, 4, 6)
                         
                         beta0 <- beta0s[ii]
                         beta1 <- beta1s[ii]
                         phi0 <- phi0s[ii]
                         phi1 <- phi1s[ii]
                         
                         sig_u <- sqrt(pi^2/17)
                         sig_v <- sqrt(2*pi^2/9)
                         rho <- 0.9
                         
                         dat <- data.frame(ae = 0, comp = 0, ntol = 0, 
                                           x = 0, clust = 0)
                         for (i in 1:nclust){
                           ae.b0 <- rnorm(1, 0, sig_u)
                           x.temp <- rep(runif(1)< 1/2, mclust[i])
                           ae.eta <- expit(beta0 + ae.b0 + beta1*x.temp)
                           ae.temp <- rbinom(mclust[i], 1, ae.eta)
                           
                           comp.b0 <- rnorm(1, 0, sig_v)
                           comp.eta <- expit(phi0 + comp.b0 + phi1*x.temp)
                           comp.temp <- rbinom(mclust[i], 1, comp.eta)
                           comp.temp <- ifelse(ae.temp == 1, 0, comp.temp)
                           
                           tol.b0 <- rnorm(1, rho*comp.b0, sqrt(1-rho^2)*sig_v)
                           tol.eta <- expit(phi0 + tol.b0 + phi1*x.temp)
                           tol.temp <- rbinom(mclust[i], 1, tol.eta)
                           ntol.temp <- ifelse(ae.temp == 1, 1, 1 - tol.temp)
                           
                           dat <- rbind(dat,
                                        data.frame(ae = ae.temp, comp = comp.temp,
                                                   ntol = ntol.temp, x = x.temp,
                                                   clust = rep(i, mclust[i])))
                         }
                         dat <- dat[-1,]
                         
                         n.chains = 1
                         n.burnin = 500
                         n.draws = 2000
                         n.thin = 1
                         
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
                         
                         beta0.post.ae <- unlist(model1.samples.ae[,ncol(model1.samples.ae[[1]]) - 2])
                         beta1.post.ae <- unlist(model1.samples.ae[,ncol(model1.samples.ae[[1]]) - 1])
                         
                         ## repeat for completion
                         model1.fit <- jags.model(file=log.int.wd,
                                                  data=list(N=nrow(dat), Y=dat$comp, clust = dat$clust,
                                                            nclust = nclust, X = dat$x,
                                                            p0= 0.0001, p1 = 0.001,
                                                            sig.up = 25),
                                                  n.chains = n.chains)
                         
                         update(model1.fit, n.burnin)
                         
                         model1.samples.comp <- coda.samples(model1.fit, c("beta0", "beta1", "sigma_int", "b0"), 
                                                             n.iter=n.draws, thin=n.thin)
                         
                         beta0.post.comp <- unlist(model1.samples.comp[,ncol(model1.samples.comp[[1]]) - 2])
                         beta1.post.comp <- unlist(model1.samples.comp[,ncol(model1.samples.comp[[1]]) - 1])
                         
                         ## repeat for non-tolerability
                         model1.fit <- jags.model(file=log.int.wd,
                                                  data=list(N=nrow(dat), Y=dat$ntol, clust = dat$clust,
                                                            nclust = nclust, X = dat$x,
                                                            p0= 0.0001, p1 = 0.001,
                                                            sig.up = 25),
                                                  n.chains = n.chains)
                         
                         update(model1.fit, n.burnin)
                         
                         model1.samples.tol <- coda.samples(model1.fit, c("beta0", "beta1", "sigma_int", "b0"), 
                                                            n.iter=n.draws, thin=n.thin)
                         
                         beta0.post.tol <- unlist(model1.samples.tol[,ncol(model1.samples.tol[[1]]) - 2])
                         beta1.post.tol <- unlist(model1.samples.tol[,ncol(model1.samples.tol[[1]]) - 1])
                         
                         save0.ae <- NULL
                         save1.ae <- NULL
                         save0.comp <- NULL
                         save1.comp <- NULL
                         save0.tol <- NULL
                         save1.tol <- NULL
                         mm <- nrow(model1.samples.ae[[1]])
                         for (i in 1:mm){
                           mu1.ae <- expit(beta0.post.ae[i] + beta1.post.ae[i] + model1.samples.ae[[1]][i, 1:nclust])
                           mu1.comp <- expit(beta0.post.comp[i] + beta1.post.comp[i] + model1.samples.comp[[1]][i, 1:nclust])
                           mu1.tol <- expit(beta0.post.tol[i] + beta1.post.tol[i] + model1.samples.tol[[1]][i, 1:nclust])
                           w1 <- rexp(nclust, 1)
                           w1 <- w1/sum(w1)
                           
                           mu0.ae <- expit(beta0.post.ae[i] + model1.samples.ae[[1]][i, 1:nclust])
                           mu0.comp <- expit(beta0.post.comp[i] + model1.samples.comp[[1]][i, 1:nclust])
                           mu0.tol <- expit(beta0.post.tol[i] + model1.samples.tol[[1]][i, 1:nclust])
                           w0 <- rexp(nclust, 1)
                           w0 <- w0/sum(w0)
                           
                           save1.ae[i] <- sum(w1*mu1.ae) 
                           save0.ae[i] <- sum(w0*mu0.ae)
                           save1.comp[i] <- sum(w1*mu1.comp) 
                           save0.comp[i] <- sum(w0*mu0.comp)
                           save1.tol[i] <- sum(w1*mu1.tol) 
                           save0.tol[i] <- sum(w0*mu0.tol)
                         }
                         ldiff.ae <- save1.ae - save0.ae
                         
                         kd.ae <- density(na.omit(ldiff.ae))
                         np.prob.ae <- mean(pnorm(d_ae, ifelse(is.finite(ldiff.ae), ldiff.ae, 
                                                               max(na.omit(ldiff.ae)) + 1), kd.ae$bw, lower.tail = FALSE))
                         
                         if (np.prob.ae < 0.00004){
                           np.prob.ae <- pnorm(d_ae, mean(na.omit(ldiff.ae)), sd(na.omit(ldiff.ae)), lower.tail = FALSE)
                         }
                         
                         ldiff.comp <- save0.comp - save1.comp
                         
                         kd.comp <- density(na.omit(ldiff.comp))
                         np.prob.comp <- mean(pnorm(d_comp, ifelse(is.finite(ldiff.comp), ldiff.comp, 
                                                                   max(na.omit(ldiff.comp)) + 1), kd.comp$bw, lower.tail = FALSE))
                         
                         if (np.prob.comp < 0.00004){
                           np.prob.comp <- pnorm(d_comp, mean(na.omit(ldiff.comp)), sd(na.omit(ldiff.comp)), lower.tail = FALSE)
                         }
                         
                         ldiff.tol <- save1.tol - save0.tol
                         
                         kd.tol <- density(na.omit(ldiff.tol))
                         np.prob.tol <- mean(pnorm(d_ntol, ifelse(is.finite(ldiff.tol), ldiff.tol, 
                                                                  max(na.omit(ldiff.tol)) + 1), kd.tol$bw, lower.tail = FALSE))
                         
                         
                         if (np.prob.tol < 0.00004){
                           np.prob.tol <- pnorm(d_ntol, mean(na.omit(ldiff.tol)), sd(na.omit(ldiff.tol)), lower.tail = FALSE)
                         }
                         
                         
                         dat <- dat[, c(-2, -3)]
                         ## now get the probability for the final analysis
                         nclust <- Mf
                         mclust <- c(mclust, rdunif(Mf - M, 4, 6))
                         
                         for (i in (M+1):Mf){
                           ae.b0 <- rnorm(1, 0, sig_u)
                           x.temp <- rep(runif(1)< 1/2, mclust[i])
                           ae.eta <- expit(beta0 + ae.b0 + beta1*x.temp)
                           ae.temp <- rbinom(mclust[i], 1, ae.eta)
                           
                           dat <- rbind(dat,
                                        data.frame(ae = ae.temp, x = x.temp,
                                                   clust = rep(i, mclust[i])))
                         }
                         
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
                         
                         beta0.post.ae <- unlist(model1.samples.ae[,ncol(model1.samples.ae[[1]]) - 2])
                         beta1.post.ae <- unlist(model1.samples.ae[,ncol(model1.samples.ae[[1]]) - 1])
                         
                         save0.ae <- NULL
                         save1.ae <- NULL
                         mm <- nrow(model1.samples.ae[[1]])
                         for (i in 1:mm){
                           mu1.ae <- expit(beta0.post.ae[i] + beta1.post.ae[i] + model1.samples.ae[[1]][i, 1:nclust])
                           w1 <- rexp(nclust, 1)
                           w1 <- w1/sum(w1)
                           
                           mu0.ae <- expit(beta0.post.ae[i] + model1.samples.ae[[1]][i, 1:nclust])
                           w0 <- rexp(nclust, 1)
                           w0 <- w0/sum(w0)
                           
                           save1.ae[i] <- sum(w1*mu1.ae) 
                           save0.ae[i] <- sum(w0*mu0.ae)
                         }
                         ldiff.ae2 <- save1.ae - save0.ae
                         
                         kd.ae2 <- density(na.omit(ldiff.ae2))
                         np.prob.ae2 <- mean(pnorm(d_ae, ifelse(is.finite(ldiff.ae2), ldiff.ae2, 
                                                                max(na.omit(ldiff.ae2)) + 1), kd.ae2$bw, lower.tail = FALSE))
                         
                         if (np.prob.ae2 < 0.00004){
                           np.prob.ae2 <- pnorm(d_ae, mean(na.omit(ldiff.ae2)), sd(na.omit(ldiff.ae2)), lower.tail = FALSE)
                         } 
                         
                         c(np.prob.ae, np.prob.comp, np.prob.tol, np.prob.ae2)
                       }
    write.csv(sim.res, paste0("scen2",ii,"_log_c_", M, ".csv"), row.names = FALSE)
  }
}

## parameters for high ICC setting
Ms <- cs[3]
Mfs <- ceiling(1.25*Ms)

beta0s <- c(-4.4041, -4.4041, -4.4041, -4.4041)
beta1s <- c(0.0000, 0.4330, 0.9937, 1.2008)
phi0s <- c(1.8075, 1.8075, 1.8075, 1.8075)
phi1s <- c(0.0000, -0.1815, -0.2237, -0.5643)
for (ii in 1:4){
  for (j in 1:length(Ms)){
    M <- Ms[j]
    Mf <- Mfs[j]
    sim.res <- foreach(k=1:m, .packages=c('rjags', 'coda', 'purrr'), .combine=rbind,
                       .options.snow=opts) %dopar% {
                         
                         set.seed(880000 + (ii-1)*j*m + (j-1)*m + k)
                         nclust <- M
                         mclust <- rdunif(nclust, 4, 6)
                         
                         beta0 <- beta0s[ii]
                         beta1 <- beta1s[ii]
                         phi0 <- phi0s[ii]
                         phi1 <- phi1s[ii]
                         
                         sig_u <- sqrt(pi^2/9)
                         sig_v <- sqrt(pi^2/3)
                         rho <- 0.9
                         
                         dat <- data.frame(ae = 0, comp = 0, ntol = 0, 
                                           x = 0, clust = 0)
                         for (i in 1:nclust){
                           ae.b0 <- rnorm(1, 0, sig_u)
                           x.temp <- rep(runif(1)< 1/2, mclust[i])
                           ae.eta <- expit(beta0 + ae.b0 + beta1*x.temp)
                           ae.temp <- rbinom(mclust[i], 1, ae.eta)
                           
                           comp.b0 <- rnorm(1, 0, sig_v)
                           comp.eta <- expit(phi0 + comp.b0 + phi1*x.temp)
                           comp.temp <- rbinom(mclust[i], 1, comp.eta)
                           comp.temp <- ifelse(ae.temp == 1, 0, comp.temp)
                           
                           tol.b0 <- rnorm(1, rho*comp.b0, sqrt(1-rho^2)*sig_v)
                           tol.eta <- expit(phi0 + tol.b0 + phi1*x.temp)
                           tol.temp <- rbinom(mclust[i], 1, tol.eta)
                           ntol.temp <- ifelse(ae.temp == 1, 1, 1 - tol.temp)
                           
                           dat <- rbind(dat,
                                        data.frame(ae = ae.temp, comp = comp.temp,
                                                   ntol = ntol.temp, x = x.temp,
                                                   clust = rep(i, mclust[i])))
                         }
                         dat <- dat[-1,]
                         
                         n.chains = 1
                         n.burnin = 500
                         n.draws = 2000
                         n.thin = 1
                         
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
                         
                         beta0.post.ae <- unlist(model1.samples.ae[,ncol(model1.samples.ae[[1]]) - 2])
                         beta1.post.ae <- unlist(model1.samples.ae[,ncol(model1.samples.ae[[1]]) - 1])
                         
                         ## repeat for completion
                         model1.fit <- jags.model(file=log.int.wd,
                                                  data=list(N=nrow(dat), Y=dat$comp, clust = dat$clust,
                                                            nclust = nclust, X = dat$x,
                                                            p0= 0.0001, p1 = 0.001,
                                                            sig.up = 25),
                                                  n.chains = n.chains)
                         
                         update(model1.fit, n.burnin)
                         
                         model1.samples.comp <- coda.samples(model1.fit, c("beta0", "beta1", "sigma_int", "b0"), 
                                                             n.iter=n.draws, thin=n.thin)
                         
                         beta0.post.comp <- unlist(model1.samples.comp[,ncol(model1.samples.comp[[1]]) - 2])
                         beta1.post.comp <- unlist(model1.samples.comp[,ncol(model1.samples.comp[[1]]) - 1])
                         
                         ## repeat for non-tolerability
                         model1.fit <- jags.model(file=log.int.wd,
                                                  data=list(N=nrow(dat), Y=dat$ntol, clust = dat$clust,
                                                            nclust = nclust, X = dat$x,
                                                            p0= 0.0001, p1 = 0.001,
                                                            sig.up = 25),
                                                  n.chains = n.chains)
                         
                         update(model1.fit, n.burnin)
                         
                         model1.samples.tol <- coda.samples(model1.fit, c("beta0", "beta1", "sigma_int", "b0"), 
                                                            n.iter=n.draws, thin=n.thin)
                         
                         beta0.post.tol <- unlist(model1.samples.tol[,ncol(model1.samples.tol[[1]]) - 2])
                         beta1.post.tol <- unlist(model1.samples.tol[,ncol(model1.samples.tol[[1]]) - 1])
                         
                         save0.ae <- NULL
                         save1.ae <- NULL
                         save0.comp <- NULL
                         save1.comp <- NULL
                         save0.tol <- NULL
                         save1.tol <- NULL
                         mm <- nrow(model1.samples.ae[[1]])
                         for (i in 1:mm){
                           mu1.ae <- expit(beta0.post.ae[i] + beta1.post.ae[i] + model1.samples.ae[[1]][i, 1:nclust])
                           mu1.comp <- expit(beta0.post.comp[i] + beta1.post.comp[i] + model1.samples.comp[[1]][i, 1:nclust])
                           mu1.tol <- expit(beta0.post.tol[i] + beta1.post.tol[i] + model1.samples.tol[[1]][i, 1:nclust])
                           w1 <- rexp(nclust, 1)
                           w1 <- w1/sum(w1)
                           
                           mu0.ae <- expit(beta0.post.ae[i] + model1.samples.ae[[1]][i, 1:nclust])
                           mu0.comp <- expit(beta0.post.comp[i] + model1.samples.comp[[1]][i, 1:nclust])
                           mu0.tol <- expit(beta0.post.tol[i] + model1.samples.tol[[1]][i, 1:nclust])
                           w0 <- rexp(nclust, 1)
                           w0 <- w0/sum(w0)
                           
                           save1.ae[i] <- sum(w1*mu1.ae) 
                           save0.ae[i] <- sum(w0*mu0.ae)
                           save1.comp[i] <- sum(w1*mu1.comp) 
                           save0.comp[i] <- sum(w0*mu0.comp)
                           save1.tol[i] <- sum(w1*mu1.tol) 
                           save0.tol[i] <- sum(w0*mu0.tol)
                         }
                         ldiff.ae <- save1.ae - save0.ae
                         
                         kd.ae <- density(na.omit(ldiff.ae))
                         np.prob.ae <- mean(pnorm(d_ae, ifelse(is.finite(ldiff.ae), ldiff.ae, 
                                                               max(na.omit(ldiff.ae)) + 1), kd.ae$bw, lower.tail = FALSE))
                         
                         if (np.prob.ae < 0.00004){
                           np.prob.ae <- pnorm(d_ae, mean(na.omit(ldiff.ae)), sd(na.omit(ldiff.ae)), lower.tail = FALSE)
                         }
                         
                         ldiff.comp <- save0.comp - save1.comp
                         
                         kd.comp <- density(na.omit(ldiff.comp))
                         np.prob.comp <- mean(pnorm(d_comp, ifelse(is.finite(ldiff.comp), ldiff.comp, 
                                                                   max(na.omit(ldiff.comp)) + 1), kd.comp$bw, lower.tail = FALSE))
                         
                         if (np.prob.comp < 0.00004){
                           np.prob.comp <- pnorm(d_comp, mean(na.omit(ldiff.comp)), sd(na.omit(ldiff.comp)), lower.tail = FALSE)
                         }
                         
                         ldiff.tol <- save1.tol - save0.tol
                         
                         kd.tol <- density(na.omit(ldiff.tol))
                         np.prob.tol <- mean(pnorm(d_ntol, ifelse(is.finite(ldiff.tol), ldiff.tol, 
                                                                  max(na.omit(ldiff.tol)) + 1), kd.tol$bw, lower.tail = FALSE))
                         
                         if (np.prob.tol < 0.00004){
                           np.prob.tol <- pnorm(d_ntol, mean(na.omit(ldiff.tol)), sd(na.omit(ldiff.tol)), lower.tail = FALSE)
                         }
                         
                         
                         dat <- dat[, c(-2, -3)]
                         ## now get the probability for the final analysis
                         nclust <- Mf
                         mclust <- c(mclust, rdunif(Mf - M, 4, 6))
                         
                         for (i in (M+1):Mf){
                           ae.b0 <- rnorm(1, 0, sig_u)
                           x.temp <- rep(runif(1)< 1/2, mclust[i])
                           ae.eta <- expit(beta0 + ae.b0 + beta1*x.temp)
                           ae.temp <- rbinom(mclust[i], 1, ae.eta)
                           
                           dat <- rbind(dat,
                                        data.frame(ae = ae.temp, x = x.temp,
                                                   clust = rep(i, mclust[i])))
                         }
                         
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
                         
                         beta0.post.ae <- unlist(model1.samples.ae[,ncol(model1.samples.ae[[1]]) - 2])
                         beta1.post.ae <- unlist(model1.samples.ae[,ncol(model1.samples.ae[[1]]) - 1])
                         
                         save0.ae <- NULL
                         save1.ae <- NULL
                         mm <- nrow(model1.samples.ae[[1]])
                         for (i in 1:mm){
                           mu1.ae <- expit(beta0.post.ae[i] + beta1.post.ae[i] + model1.samples.ae[[1]][i, 1:nclust])
                           w1 <- rexp(nclust, 1)
                           w1 <- w1/sum(w1)
                           
                           mu0.ae <- expit(beta0.post.ae[i] + model1.samples.ae[[1]][i, 1:nclust])
                           w0 <- rexp(nclust, 1)
                           w0 <- w0/sum(w0)
                           
                           save1.ae[i] <- sum(w1*mu1.ae) 
                           save0.ae[i] <- sum(w0*mu0.ae)
                         }
                         ldiff.ae2 <- save1.ae - save0.ae
                         
                         kd.ae2 <- density(na.omit(ldiff.ae2))
                         np.prob.ae2 <- mean(pnorm(d_ae, ifelse(is.finite(ldiff.ae2), ldiff.ae2, 
                                                                max(na.omit(ldiff.ae2)) + 1), kd.ae2$bw, lower.tail = FALSE))
                         
                         if (np.prob.ae2 < 0.00004){
                           np.prob.ae2 <- pnorm(d_ae, mean(na.omit(ldiff.ae2)), sd(na.omit(ldiff.ae2)), lower.tail = FALSE)
                         } 
                         
                         c(np.prob.ae, np.prob.comp, np.prob.tol, np.prob.ae2)
                       }
    write.csv(sim.res, paste0("scen3",ii,"_log_c_", M, ".csv"), row.names = FALSE)
  }
}

## get the estimated operating characteristics based on the confirmatory simulations

## implement for the interim analysis
res.int.dat <- NULL
for (s in 1:3){
  for (j in 1:4){
    M <- cs[s]
    res_temp <- read.csv(paste0("scen",s,j,"_log_c_", M, ".csv"))
    res.int.dat <- c(res.int.dat, mean(1-res_temp[,1] > gam[1] & 1 - res_temp[,2] > gam[2] & 1 - res_temp[,3] > gam[3]))
  }
}

## format into table
res.int.dat <- matrix(res.int.dat, nrow = 3, byrow = TRUE)

## implement for the final analysis
res.fin.dat <- NULL
for (s in 1:3){
  for (j in 1:4){
    M <- cs[s]
    res_temp <- read.csv(paste0("scen",s,j,"_log_c_", M, ".csv"))
    res.fin.dat <- c(res.fin.dat, mean(1-res_temp[,1] > gam[1] & 1 - res_temp[,2] > gam[2] & 1 - res_temp[,3] > gam[3] & 1- res_temp[,4] > 0.94))
  }
}

## format into table
res.fin.dat <- matrix(res.fin.dat, nrow = 3, byrow = TRUE)

## format into final table
tab.seq <- rbind(res.int.lin, res.int.dat, res.fin.lin, res.fin.dat)
tab.seq <- tab.seq[c(1, 4, 7, 10, 2, 5, 8, 11, 4, 6, 9, 12),]

## save results to .csv file
write.csv(tab.seq, "seq_table.csv", row.names = FALSE)