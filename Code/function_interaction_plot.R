# the first term of the interaction: cvarname2
# bvarname1: the second term of the interaction
#cval1: the first value of cvarname2
# fd_no: bvarname1 = 0, 

logit_interEffect_binary <- function(ModelResults, cvarname2, bvarname1,
                                   cval1, cval2, n.sim = 1000, data, clusterid, label = c("no bvarname1", "bvarname1")){
    #get a sim objective
    require(arm)
    library(multiwayvcov)
    library(lmtest)
    library(ggplot2)
    library(dplyr)
    cluster <- data[,clusterid]
    vcov_cluster <- cluster.vcov(ModelResults, cluster)
    coef_cluster <- coeftest(ModelResults, vcov = vcov_cluster)
    set.seed(12345)
    draw <- mvrnorm(n= n.sim, coef(ModelResults), vcov_cluster)
  ##set simulation
  df <- array(NA, c(2, 4))
  
  X1 <- model.matrix(ModelResults)
  X2 <- model.matrix(ModelResults)
  X1[, bvarname1] = 0
  X1[, cvarname2] = cval1
  X1[, paste(cvarname2, bvarname1, sep = ":")] = 0*cval1
  
  X2[,bvarname1] = 0
  X2[,cvarname2] = cval2
  X2[, paste(cvarname2, bvarname1, sep = ":")] = 0*cval2
  
  fd_no <- apply(apply(X2, 1, function (x) draw %*% x) - 
                   apply(X1, 1, function (x) draw %*% x), 1, mean) #1 indicates row, 2= columns
  
  X3 <- model.matrix(ModelResults)
  X4 <- model.matrix(ModelResults)
  X3[, bvarname1] = 1
  X3[, cvarname2] = cval1
  X3[, paste(cvarname2, bvarname1, sep = ":")] = 1*cval1
  
  X4[,bvarname1] = 1
  X4[,cvarname2] = cval2
  X4[, paste(cvarname2, bvarname1, sep = ":")] = 1*cval2
  
  fd_yes <- apply(apply(X4, 1, function (x) draw %*% x) - 
                    apply(X3, 1, function (x) draw %*% x), 1, mean) #1 indicates row, 2= columns
  
  df[1, 1] <- 0
  df[1, 2] <- mean(fd_no)
  df[1, 3:4] <- quantile(fd_no, probs = c(.05,.95))
  df[2, 1] <- 1
  df[2, 2] <- mean(fd_yes)
  df[2, 3:4] <- quantile(fd_yes, probs = c(.05,.95))
  
  df_plot <- as.data.frame(df)
  
  
  colnames(df_plot) <- c("type", "mean", "lo", "hi")     
  
  p <- ggplot(df_plot, aes(x=type, y=mean)) +
    geom_hline(aes(yintercept=0), linetype=2, color = "black") + 
    geom_point(size=4) + 
    geom_linerange(aes(ymin=lo, ymax=hi),alpha = 1, size = 1) + 
    theme_bw() + scale_x_discrete(limits = c(0,1), labels = label)+
    xlab("") + ylab("First Differences")+
    theme(legend.position="none",
          legend.title=element_blank(),
          axis.text = element_text(size=14),
          text = element_text(size=14),
          plot.title = element_text(hjust = .5, size = 16, face = "bold"))
  
  return(p)
}


