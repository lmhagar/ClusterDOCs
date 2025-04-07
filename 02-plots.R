## BEGIN SETUP ##

## load necessary packages
## run 01-simulations.R first to get the .csv files used in this code file
require(ggplot2)
require(cowplot)
require(ggpubr)

## load colour palette
cbb <- c("#000000", "#E69F00", "#56B4E9", "#009E73", 
         "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

## functions to convert results from .csv files to logits
## of posterior probabilities
logit <- function(x){log(x) - log(1-x)}
## try to get "power" curves
l1 <- function(x){
  -1*logit(x)
}

## function to obtain vectors to plot
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

## similar function to get bootstrap CIs for design
## operating characteristics (and sample size)
bootstrap_doc <- function(ml, mu, cl, cu, lb, ub, gam, Gam){
  # ml is lower matrix
  # mu is upper matrix
  # cl is lower c
  # cu is upper c
  # lb is lower bound
  # ub is upper bound
  # gam are the thresholds
  
  ll <- apply(ml, 2, l1)
  lu <- apply(mu, 2, l1)
  
  for (j in 1:ncol(ll)){
    ll[,j] <- ifelse(ll[,j] == -Inf, min(subset(ll, is.finite(ll[,j]))[,j]) - 1, ll[,j])
    ll[,j] <- ifelse(ll[,j] == Inf, max(subset(ll, is.finite(ll[,j]))[,j]) + 1, ll[,j])
  }
  
  for (j in 1:ncol(lu)){
    lu[,j] <- ifelse(lu[,j] == -Inf, min(subset(lu, is.finite(lu[,j]))[,j]) - 1, lu[,j])
    lu[,j] <- ifelse(lu[,j] == Inf, max(subset(lu, is.finite(lu[,j]))[,j]) + 1, lu[,j])
  }
  
  ll <- cbind(ll, seq(1, nrow(ll), 1))
  lu <- cbind(lu, seq(1, nrow(lu), 1))
  
  # AE
  ll_s <- ll[order(ll[,1]),1]
  lu_s <- lu[order(lu[,1]),1]
  
  l_slope <- (lu_s - ll_s)/(cu-cl)
  l_int <- ll_s - l_slope*cl
  
  l_slope[ll[order(ll[,1]),4]] <- l_slope 
  l_int[ll[order(ll[,1]),4]] <- l_int 
  
  slopes <- NULL
  slopes <- cbind(slopes, l_slope)
  
  ints <- NULL
  ints <- cbind(ints, l_int)
  
  # comp
  ll_s <- ll[order(ll[,2]),2]
  lu_s <- lu[order(lu[,2]),2]
  
  l_slope <- (lu_s - ll_s)/(cu-cl)
  l_int <- ll_s - l_slope*cl
  
  l_slope[ll[order(ll[,2]),4]] <- l_slope 
  l_int[ll[order(ll[,2]),4]] <- l_int 
  
  slopes <- cbind(slopes, l_slope)
  ints <- cbind(ints, l_int)
  
  # tol
  ll_s <- ll[order(ll[,3]),3]
  lu_s <- lu[order(lu[,3]),3]
  
  l_slope <- (lu_s - ll_s)/(cu-cl)
  l_int <- ll_s - l_slope*cl
  
  l_slope[ll[order(ll[,3]),4]] <- l_slope 
  l_int[ll[order(ll[,3]),4]] <- l_int 
  
  slopes <- cbind(slopes, l_slope)
  ints <- cbind(ints, l_int)
  
  ## compile
  
  cs <- seq(lb, ub,1)
  lae <- matrix(0, nrow = nrow(slopes), ncol=length(cs))
  lcomp <- matrix(0, nrow = nrow(slopes), ncol=length(cs))
  ltol <- matrix(0, nrow = nrow(slopes), ncol=length(cs))
  for (i in 1:length(cs)){
    lae[,i] <- ints[,1] + slopes[,1]*cs[i]
    lcomp[,i] <- ints[,2] + slopes[,2]*cs[i]
    ltol[,i] <- ints[,3] + slopes[,3]*cs[i]
  }
  
  ylae <- lae > logit(gam[1])
  ylcomp <- lcomp > logit(gam[2])
  yltol <- ltol > logit(gam[3])
  
  all3 <- ylae + ylcomp + yltol
  
  out.temp <- cbind(cs, colMeans(all3 >= 3))
  
  return(as.numeric(out.temp[,2]))
  
}

## new simulations
cores=detectCores()
cl <- makeSOCKcluster(cores[1]-1)

m <- 10000
registerDoSNOW(cl)
pb <- txtProgressBar(max = m, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

## get bootstrap CIs and plots for each ICC setting
MM <- 10000
for (s in 1:3){

  Ms <- seq(100, 200, 10)
  gam <- c(0.8, 0.55, 0.55)
  
  # process to implement the bootstrap resampling
  for (j in 1:4){
    
    assign(paste0("res", s, j, "_120"), read.csv(paste0("scen", s, j, "_log_c_120.csv"))[,1:3])
    assign(paste0("res", s, j, "_180"), read.csv(paste0("scen", s, j, "_log_c_180.csv"))[,1:3])
    
    samp1 <- get(paste0("res", s, j, "_120"))
    samp2 <- get(paste0("res", s, j, "_180"))
           
    boot.temp <- foreach(k=1:MM, .packages=c('rjags', 'coda', 'purrr'), .combine=rbind,
                         .options.snow=opts) %dopar% {

                           samp1.temp <- samp1[sample(1:m, m, replace = TRUE),]
                           samp2.temp <- samp2[sample(1:m, m, replace = TRUE),]
                           
                           bootstrap_doc(samp1.temp, samp2.temp, 120, 180, 100, 200,
                                    gam, 0.85)
                         }
    
    write.csv(boot.temp, paste0("CIs_OC", s, j, ".csv"), row.names = FALSE)
    
    # for scenario 2, get bootstrap confidence interval for the sample size recommendation
    if (j == 2){
      boot_c = apply(boot.temp, 1, function(x, Gam, low){which.max(as.numeric(x) >= Gam) + low - 1}, 
                     Gam = 0.85, low = 100)
      write.csv(boot_c, paste0("CI_c", s, j, ".csv"), row.names = FALSE)
      print(quantile(boot_c, c(0.025, 0.975)))
    }
  }
  
  # extract 95% bootstrap confidence intervals for design OCs and estimates from naive simulation
  for (j in 1:4){
    assign(paste0("scen", s, j, ".app"), lin_app(get(paste0("res", s,j,"_120")), get(paste0("res", s,j,"_180")), 
                                               120, 180, 100, 200, gam))
    
    scen.act <- NULL
    for (i in 1:length(Ms)){
      M <- Ms[i]
      res_temp <- read.csv(paste0("scen",s,j,"_log_c_", M, ".csv"))
      scen.act[i] <- mean(1-res_temp[,1] > gam[1] & 1 - res_temp[,2] > gam[2] & 1 - res_temp[,3] > gam[3])
      scen.u95 <- as.numeric(apply(read.csv(paste0("CIs_OC", s, j, ".csv")), 2, quantile, probs = 0.975))
      scen.l95 <- as.numeric(apply(read.csv(paste0("CIs_OC", s, j, ".csv")), 2, quantile, probs = 0.025))
    }
    assign(paste0("scen", s,j, ".act"), scen.act)
    assign(paste0("scen", s,j, ".u95"), scen.u95)
    assign(paste0("scen", s,j, ".l95"), scen.l95)
  }
  
  ## create data frame for each scenario and setting combination (for plotting)
  for (j in 1:4){
    assign(paste0("df", s, j),
           data.frame(c = c(get(paste0("scen", s, j, ".app"))[,1], Ms,
                            get(paste0("scen", s, j, ".app"))[,1], get(paste0("scen", s, j, ".app"))[,1]), 
                      prob = c(get(paste0("scen", s, j, ".app"))[,2], get(paste0("scen", s,j, ".act")),
                               get(paste0("scen", s,j, ".u95")), get(paste0("scen", s,j, ".l95"))),
                      type = c(rep(1, 101), rep(2, 11), rep(3, 101), rep(4, 101))))
  }
  
  # get subplots for each ICC setting (one for each scenario)
  for (j in 1:4){
    assign(paste0("plot", s, j), 
           ggplot(get(paste0("df", s, j)), 
                                         aes(x=c, y=prob, color = as.factor(type), linetype = as.factor(type))) + theme_bw() +
      geom_line() +
      scale_color_manual(name = "", labels = c("Estimated", "Simulated", 
                                               "95% Bootstrap CI", "95% Bootstrap CI"),
                         values = c("black", cbb[6], "black", "black")) +
      scale_linetype_manual(name = "", labels = c("Estimated", "Simulated", 
                                                  "95% Bootstrap CI", "95% Bootstrap CI"),
                            values = c("solid", "solid", "dashed", "dashed")) +
      labs(color  = "", linetype = "") +
      labs(x= bquote(italic(c)), y= bquote('Continuation Probability')) +
      theme(plot.title = element_text(size=20,face="bold",
                                      margin = margin(t = 0, 0, 5, 0))) +
      theme(axis.text=element_text(size=16),
            axis.title=element_text(size=18)) +
      theme(legend.text=element_text(size=18)) +
      theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
      theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
      theme(legend.position="bottom") +
      labs(title=paste0("Scenario ", j)))
  }
  
}

## get simplified legend for final plots
df111 <- subset(df11, df11$type < 4)
s = 1; j = 1
plot111 <- ggplot(df111, 
              aes(x=c, y=prob, color = as.factor(type), linetype = as.factor(type))) + theme_bw() +
         geom_line() +
         scale_color_manual(name = "", labels = c("Estimated   ", "Simulated   ", 
                                                  "95% Bootstrap Confidence Interval"),
                            values = c("black", cbb[6], "black")) +
         scale_linetype_manual(name = "", labels = c("Estimated   ", "Simulated   ", 
                                                     "95% Bootstrap Confidence Interval"),
                               values = c("solid", "solid", "dashed")) +
         labs(color  = "", linetype = "") +
         labs(x= bquote(italic(c)), y= bquote('Continuation Probability')) +
         theme(plot.title = element_text(size=20,face="bold",
                                         margin = margin(t = 0, 0, 5, 0))) +
         theme(axis.text=element_text(size=16),
               axis.title=element_text(size=18)) +
         theme(legend.text=element_text(size=18)) +
         theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
         theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
         theme(legend.position="bottom") +
         labs(title=paste0("Scenario ", j))

## create final plots for each ICC setting
for (s in 1:3){
  figp.row1 <- plot_grid(get(paste0("plot", s, 1)) + theme(plot.margin=unit(c(0.25,0.5,0.25,0.5),"cm")) +
                           theme(legend.position="none"), 
                         get(paste0("plot", s, 2)) + theme(plot.margin=unit(c(0.25,0.5,0.25,0.5),"cm")) +
                           theme(legend.position="none"),
                         rel_widths = c(1, 1))
  
  figp.row2 <- plot_grid(get(paste0("plot", s, 3)) + theme(plot.margin=unit(c(0.25,0.5,0.25,0.5),"cm")) +
                           theme(legend.position="none"), 
                         get(paste0("plot", s, 4)) + theme(plot.margin=unit(c(0.25,0.5,0.25,0.5),"cm")) +
                           theme(legend.position="none") + ylim(0.05, 0.11),
                         rel_widths = c(1, 1))
  
  figp <- plot_grid(figp.row1, figp.row2, nrow = 2)
  
  fig_final <- plot_grid(figp, get_legend(get(paste0("plot111"))), ncol = 1, rel_heights = c(2, .1))
  
  # output as .pdf file for the article
  pdf(file = paste0("Fig", s, "DOC.pdf"),   # The directory you want to save the file in
      width = 10.5, # The width of the plot in inches
      height = 7) # The height of the plot in inches
  
  fig_final
  
  dev.off()
}