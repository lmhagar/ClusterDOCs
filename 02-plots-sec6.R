## BEGIN SETUP ##

## load necessary packages
## run 01-simulations-sec6.R first to get the .csv files used in this code file
require(ggplot2)
require(cowplot)
require(ggpubr)

## load colour palette
cbb <- c("#000000", "#E69F00", "#56B4E9", "#009E73", 
         "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

## functions to convert results from .csv files to logits
## of posterior probabilities
logit <- function(x){log(x) - log(1-x)}
l1 <- function(x){
  -1*logit(x)
}

## this function is used to estimate operating characteristics
## based on our linear approximations 
lin_app <- function(ml, mu, cl, cu, lb, ub, gam){
  # ml is the vector of posterior probabilities for c0
  # mu is the vector of posterior probabilities for c1
  # cl is c0
  # cu is c1
  # lb is lower bound for c (in terms of plotting)
  # ub is upper bound for c (in terms of plotting)
  # gam is the threshold
  
  # get logits for the c0 (ll) and c2 (lu)
  ll <- sapply(ml, l1)
  lu <- sapply(mu, l1)
  
  # adjust any infinite logits
  ll <- ifelse(ll == -Inf, min(subset(ll, is.finite(ll))) - 1, ll)
  ll <- ifelse(ll == Inf, max(subset(ll, is.finite(ll))) + 1, ll)
  
  lu <- ifelse(lu == -Inf, min(subset(lu, is.finite(lu))) - 1, lu)
  lu <- ifelse(lu == Inf, max(subset(lu, is.finite(lu))) + 1, lu)
  
  # sort logits
  ll_s <- sort(ll)
  lu_s <- sort(lu)
  
  # get slopes and intercepts of linear approximations based on order stats
  slopes <- (lu_s - ll_s)/(cu-cl)
  ints <- ll_s - slopes*cl
  
  # get matrices of new posterior probabilities corresponding to
  # linear approximations
  cs <- seq(lb, ub,1)
  lae <- matrix(0, nrow = length(slopes), ncol=length(cs))
  ## each column corresponds to a sample size
  ## the rows represent different "logits" based on linear approximations
  for (i in 1:length(cs)){
    lae[,i] <- ints + slopes*cs[i]
  }
  
  # create indicator to determine if decision criteria satisfied
  ylae <- lae > logit(gam)
  
  ## return matrix with cluster counts and estimated power (or type I error rate)
  return(cbind(cs, colMeans(ylae >= 1)))
  
}

## similar function to get bootstrap CIs for design
## operating characteristics (and sample size)
bootstrap_doc <- function(ml, mu, cl, cu, lb, ub, gam){
  # ml is the vector of posterior probabilities for c0
  # mu is the vector of posterior probabilities for c1
  # cl is c0
  # cu is c1
  # lb is lower bound for c (in terms of plotting)
  # ub is upper bound for c (in terms of plotting)
  # gam is the threshold
  
  # get logits for the c0 (ll) and c2 (lu)
  ll <- sapply(ml, l1)
  lu <- sapply(mu, l1)
  
  # adjust any infinite logits
  ll <- ifelse(ll == -Inf, min(subset(ll, is.finite(ll))) - 1, ll)
  ll <- ifelse(ll == Inf, max(subset(ll, is.finite(ll))) + 1, ll)
  
  lu <- ifelse(lu == -Inf, min(subset(lu, is.finite(lu))) - 1, lu)
  lu <- ifelse(lu == Inf, max(subset(lu, is.finite(lu))) + 1, lu)
  
  # sort logits
  ll_s <- sort(ll)
  lu_s <- sort(lu)
  
  # get slopes and intercepts of linear approximations based on order stats
  slopes <- (lu_s - ll_s)/(cu-cl)
  ints <- ll_s - slopes*cl
  
  # get matrices of new posterior probabilities corresponding to
  # linear approximations
  cs <- seq(lb, ub,1)
  lae <- matrix(0, nrow = length(slopes), ncol=length(cs))
  ## each column corresponds to a sample size
  ## the rows represent different "logits" based on linear approximations
  for (i in 1:length(cs)){
    lae[,i] <- ints + slopes*cs[i]
  }
  
  # create indicator to determine if decision criteria satisfied
  ylae <- lae > logit(gam)
  
  ## output just the operating characteristics this time
  out.temp <- cbind(cs, colMeans(ylae >= 1))
  return(as.numeric(out.temp[,2]))
  
}

## set up parallelization
cores=detectCores()
cl <- makeSOCKcluster(cores[1]-1)

## set up progress bar
m <- 10000
registerDoSNOW(cl)
pb <- txtProgressBar(max = m, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

## get bootstrap CIs and plots for each ICC setting
MM <- 10000

## loop over each ICC setting
for (s in 1:3){

  ## set array of cluster counts and decision thresholds
  Ms <- seq(80, 160, 10)
  gam <- 0.970
  
  # process to implement the bootstrap resampling for each scenario
  for (j in 1:4){
    
    ## load in .csv files of posterior probabilities corresponding to c = 100 and c = 140
    assign(paste0("res", s, j, "_100"), read.csv(paste0("scen", s, j, "_log_c_100.csv"))[,1])
    assign(paste0("res", s, j, "_140"), read.csv(paste0("scen", s, j, "_log_c_140.csv"))[,1])
    
    samp1 <- get(paste0("res", s, j, "_100"))
    samp2 <- get(paste0("res", s, j, "_140"))

    ## to construct the bootstrap CI, resample from the initial samples and use bootstrap_doc() function
    boot.temp <- foreach(k=1:MM, .packages=c('rjags', 'coda', 'purrr'), .combine=rbind,
                         .options.snow=opts) %dopar% {

                           samp1.temp <- samp1[sample(1:m, m, replace = TRUE)]
                           samp2.temp <- samp2[sample(1:m, m, replace = TRUE)]

                           bootstrap_doc(samp1.temp, samp2.temp, 100, 140, 80, 160,
                                    gam)
                         }

    ## save matrices of operating characteristics for later use to construct the CIs
    write.csv(boot.temp, paste0("CIs_OC", s, j, ".csv"), row.names = FALSE)

    # for scenario 1, get bootstrap confidence interval for the cluster recommendation
    if (j == 1){
      ## this function call gets the smallest c for each boostrap replicate where power is large enough
      boot_c = apply(boot.temp, 1, function(x, Gam, low){which.max(as.numeric(x) >= Gam) + low - 1},
                     Gam = 0.8, low = 80)
      ## save and output the results for the boostrap CI
      write.csv(boot_c, paste0("CI_c", s, j, ".csv"), row.names = FALSE)
      print(quantile(boot_c, c(0.025, 0.975)))
    }
  }
  
  ## extract 95% bootstrap confidence intervals for design OCs and estimates from naive simulation
  for (j in 1:4){
    ## get solid black line using our linear approximations
    assign(paste0("scen", s, j, ".app"), lin_app(get(paste0("res", s,j,"_100")), get(paste0("res", s,j,"_140")), 
                                               100, 140, 80, 160, gam))
    
    ## read in the csv files to get operating characteristic estimates with naive simulation
    scen.act <- NULL
    
    ## for each sample size in the array
    for (i in 1:length(Ms)){
      M <- Ms[i]
      res_temp <- read.csv(paste0("scen",s,j,"_log_c_", M, ".csv"))
      ## get estimated power or type I error using empirical average
      scen.act[i] <- mean(1-res_temp[,1] > gam)
      ## get the dotted black lines by taking the percentiles from the bootstrap samples
      scen.u95 <- as.numeric(apply(read.csv(paste0("CIs_OC", s, j, ".csv")), 2, quantile, probs = 0.975))
      scen.l95 <- as.numeric(apply(read.csv(paste0("CIs_OC", s, j, ".csv")), 2, quantile, probs = 0.025))
    }
    ## save results to vectors for use later
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
                      type = c(rep(1, 81), rep(2, 9), rep(3, 81), rep(4, 81))))
  }
  
  ## get titles for the plot
  titles <- c("1: Clearly Acceptable", "2: Acceptable", "3: Barely Acceptable", "4: Unacceptable")
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
      labs(x= bquote(italic(c)), y= bquote('Non-inferiority Probability')) +
      theme(plot.title = element_text(size=20,face="bold",
                                      margin = margin(t = 0, 0, 5, 0))) +
      theme(axis.text=element_text(size=16),
            axis.title=element_text(size=18)) +
      theme(legend.text=element_text(size=18)) +
      theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
      theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
      theme(legend.position="bottom") +
      labs(title=titles[j]))
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
         labs(x= bquote(italic(c)), y= bquote('Non-inferiority Probability')) +
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
                           theme(legend.position="none") + ylim(0, 0.05),
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