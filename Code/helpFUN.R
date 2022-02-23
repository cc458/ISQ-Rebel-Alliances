## functions
scaled_coef_cluster <- function(ModelResults, data, clusterid, subvar = FALSE, vars = NULL) {
  library(multiwayvcov)
  library(lmtest)
  library(dplyr)
  require(ggplot2)
  library(gridExtra)
  library(tidyr)
  library(stringr)
  library(purrr)
  
  cluster <- data[, clusterid]
  
  vcov_cluster = lapply(ModelResults, function(x) FUN = cluster.vcov(x, 
                                                                     cluster))
  modcoef = Map(function(x, y) coeftest(x, y), x = ModelResults, y = vcov_cluster)
  modcoef <- lapply(modcoef, function(x) FUN =round(x, digits = 4))  
  modelcoef = lapply(modcoef, function(x) FUN = x[, "Estimate"])
  modelse = lapply(modcoef, function(x) FUN = x[, "Std. Error"])
  varnames = lapply(ModelResults, function(x) FUN = names(x$coefficients))
  dfplot <- list()
  for (i in 1:length(modelcoef)){
    ad = 0.5 / modelse[[i]]
    modelcoef[[i]] = modelcoef[[i]]* ad
    modelse[[i]] = modelse[[i]]*ad 
    
    ylo95 =modelcoef[[i]] - qt(.975, nrow(data))*(modelse[[i]])
    yhi95 =modelcoef[[i]] + qt(.975, nrow(data))*(modelse[[i]])
    ylo90 =modelcoef[[i]] - qt(.95, nrow(data))*(modelse[[i]])
    yhi90 =modelcoef[[i]] + qt(.95, nrow(data))*(modelse[[i]])
    dfplot[[i]] = data.frame(varnames[[i]], modelcoef[[i]], modelse[[i]], ylo95, yhi95,
                             ylo90, yhi90,
                             Model = paste('Model',seq_along(modelcoef)[i]))
  }
  
  coefficientdata = do.call(rbind,dfplot)
  colnames(coefficientdata) <- c("Variables", "Coefficients", "S.E",
                                 "Low95CI", "High95CI", 
                                 "Low90CI", "High90CI",
                                 "Model.Name")
  df <- coefficientdata %>%
    dplyr::mutate(Variables = as.character(Variables)) %>%
    dplyr::arrange(Variables)
  
  if(subvar == TRUE){
    df <- df[grep(paste(vars, collapse = "|"), df$Variables),]
    #re-order the variables: if factor, the order of factor levels is used; if character, an alphabetical order ist used
    df$Variables <- factor(df$Variables, levels = vars)
  } 
  
  p = ggplot(df, aes(colour = Model.Name)) + 
    geom_hline(yintercept = 0, lty = 2) +
    geom_linerange(aes(x = Variables, ymin = Low90CI,
                       ymax = High90CI),
                   lwd = 2, position = position_dodge(width = 1/2)) +
    geom_pointrange(aes(x = Variables, y = Coefficients, ymin = Low95CI,
                        ymax = High95CI, shape = Model.Name),
                    lwd = 1, position = position_dodge(width = 1/2),
                    fill = "WHITE") +
    coord_flip() + 
    theme_bw() + xlab("") + ylab("Rescaled Coefficient") + 
    theme(legend.position="bottom",
          legend.title=element_blank(),
          axis.text = element_text(size=14),
          text = element_text(size=14),
          plot.title = element_text(hjust = .5, size = 16, face = "bold"))
  
  return(p)   
} 


## function for plotting marginal effect for vglm
vglm_sim <- function(mod, var, n_sim){
  library(ggplot2); library(ggthemes); library(viridis)
  library(ggridges); library(purrr); library(dplyr)
  library(tidyr)
  # Simulate from sampling distribution of parameters #
  set.seed(1234)
  sim <- mvrnorm(n_sim, coef(mod), vcov(mod))
  formula <- VGAM::formula(mod)
  
  sim_1 <- sim[, grep(":1", colnames(sim))]
  sim_2 <- sim[, grep(":2", colnames(sim))]
  data0 <-data1 <- as.data.frame(model.matrix(lm(formula, data= data)))
  
  data0[,var] <- 1
  data1[,var] <- 0
  n <- nrow(mod@predictors)
  # Calculate effect of MarxistIdeol on each prob #
  
  #p.1 <- array(NA, c(n_sim, n, length(mod@extra$colnames.y)))
  #p.2 <- array(NA, c(n_sim, length(mod@extra$colnames.y)))
  p.1 <- array(NA, c(n_sim, n, length(mod@extra$colnames.y)-1))
  p.2 <- array(NA, c(n_sim, length(mod@extra$colnames.y)-1))
  
  #empty containers
  X1 <- as.matrix(data0)
  X2 <- as.matrix(data1)
  
  for (i in 1:n_sim){
    
    # Calculate exp(XB) for each of the 2 models
    eXB1 <- array(NA, c(2,n))
    eXB2 <- array(NA, c(2,n))
    
    eXB1[1, ] <- apply(X1, 1, function (x) exp(sim_1[i,] %*% x))
    eXB1[2, ] <- apply(X1, 1, function (x) exp(sim_2[i,] %*% x))
    
    eXB2[1, ] <- apply(X2, 1, function (x) exp(sim_1[i,] %*% x))
    eXB2[2, ] <- apply(X2, 1, function (x) exp(sim_2[i,] %*% x))
    
    # Calculate change in prob of each outcome
    #p.1[i, ,1] <- (1 / (1 + eXB1[1,] + eXB1[2,])) - (1 / (1 + eXB2[1,] + eXB2[2,]))
    p.1[i, ,1] <- (eXB1[1,] / (1 + eXB1[1,] + eXB1[2,])) - (eXB2[1,] / (1 + eXB2[1,] + eXB2[2,]))
    p.1[i, ,2] <- (eXB1[2,] / (1 + eXB1[1,] + eXB1[2,])) - (eXB2[2,] / (1 + eXB2[1,] + eXB2[2,]))
    
    # Average over data
    #p.2[i,1] <- mean(p.1[i, ,1])
    p.2[i,1] <- mean(p.1[i, ,1])
    p.2[i,2] <- mean(p.1[i, ,2])
    
  }
  
  df <- as.data.frame(p.2)
  #colnames(df) <- mod@extra$colnames.y
  colnames(df) <- mod@extra$colnames.y[-1]
  
  df_plot <- tidyr::gather(df, key = "type", value = "Pred")
  
  
  fig <- ggplot(df_plot, aes(x = Pred, y = type, height=..density.., fill = type)) +
    geom_density_ridges(col = "grey70", scale = 1.2, show.legend = F) +
    scale_fill_viridis(discrete = TRUE) +
    geom_vline(xintercept = 0, colour = gray(1/2), lty = 2) +
    theme_ridges(font_size = 13, grid = F, center_axis_labels = T) +
    theme(axis.title.y = element_blank(), 
          text = element_text(size=13),
          plot.title = element_text(hjust = .5, size = 15, face = "bold")) 
  
  return(fig)
  
}


## functions
mlogit_scaled_coef <- function(fit, subvar = FALSE, vars = NULL) {
  library(dplyr)
  require(ggplot2)
  library(gridExtra)
  library(tidyr)
  library(stringr)
  library(purrr)
  
  coef_fit <- lapply(fit, function(x) FUN =  coef(x))
  sd_fit <- lapply(fit, function(x) FUN =  sqrt(diag(vcov(x))))
  
  coef1 <- lapply(coef_fit, function(x) FUN = x[grep(":1", names(x))])
  coef2 <- lapply(coef_fit, function(x) FUN = x[grep(":2", names(x))])
  
  
  sd_1 <- lapply(sd_fit, function(x) FUN = x[grep(":1", names(x))])
  sd_2 <- lapply(sd_fit, function(x) FUN = x[grep(":2", names(x))])
  
  y_names <- lapply(fit, function(x) FUN = x@extra$colnames.y)
  x_names <- lapply(fit, function(X) FUN = colnames(X@x))
  
  df1plot <- list()
  
  for (i in 1:length(fit)){
    ad = 0.5/sd_1[[i]]
    coef1[[i]] = coef1[[i]]*ad
    sd_1[[i]] = sd_1[[i]]*ad
    
    ylo95 =  coef1[[i]] - 1.96*sd_1[[i]]
    yhi95 =  coef1[[i]] + 1.96*sd_1[[i]]
    ylo90 =  coef1[[i]] - 1.64*sd_1[[i]]
    yhi90 =  coef1[[i]] + 1.64*sd_1[[i]]
    
    df1plot[[i]] = data.frame(x_names[[i]], coef1[[i]], sd_1[[i]],
                              ylo95, yhi95, ylo90, yhi90, DV = y_names[[i]][2],
                              Model = paste("Model", seq_along(fit)[i]))
    
  }
  
  df1plot <- do.call(rbind, df1plot)
  colnames(df1plot) <- c("Variables", "Coefficients", "S.E",
                         "Low95CI", "High95CI", 
                         "Low90CI", "High90CI","DV",
                         "Model.Name")
  
  df2plot <- list()
  
  for (i in 1:length(fit)){
    ad = 0.5/sd_2[[i]]
    coef2[[i]] = coef2[[i]]*ad
    sd_2[[i]] = sd_2[[i]]*ad
    
    ylo95 =  coef2[[i]] - 1.96*sd_2[[i]]
    yhi95 =  coef2[[i]] + 1.96*sd_2[[i]]
    ylo90 =  coef2[[i]] - 1.64*sd_2[[i]]
    yhi90 =  coef2[[i]] + 1.64*sd_2[[i]]
    
    df2plot[[i]] = data.frame(x_names[[i]], coef2[[i]], sd_2[[i]],
                              ylo95, yhi95, ylo90, yhi90, DV = y_names[[i]][3],
                              Model = paste("Model", seq_along(fit)[i]))
    
  }
  
  df2plot <- do.call(rbind, df2plot)
  colnames(df2plot) <- c("Variables", "Coefficients", "S.E",
                         "Low95CI", "High95CI", 
                         "Low90CI", "High90CI","DV",
                         "Model.Name")
  
  
  dfplot <- bind_rows(df1plot, df2plot)
  
  
  dfplot <- dfplot %>%
    dplyr::mutate(Variables = as.character(Variables)) %>%
    dplyr::arrange(Variables)
  
  if(subvar == TRUE){
    dfplot <- dfplot[grep(paste(vars, collapse = "|"), dfplot$Variables),]
    #re-order the variables: if factor, the order of factor levels is used; if character, an alphabetical order ist used
    dfplot$Variables <- factor(dfplot$Variables, levels = vars)
  } 
  
  
  dfplot <- dfplot %>% 
    dplyr::mutate(DV = factor(DV, levels = c("formal", "informal")))
  
  levels(dfplot$DV)
  levels(dfplot$DV) <- c("Formal Alliance", "Informal Alliance")
  
  p = ggplot(dfplot, aes(colour = Model.Name)) + 
    geom_hline(yintercept = 0, lty = 2) +
    geom_linerange(aes(x = Variables, ymin = Low90CI,
                       ymax = High90CI),
                   lwd = 2, position = position_dodge(width = 1/2)) +
    geom_pointrange(aes(x = Variables, y = Coefficients, ymin = Low95CI,
                        ymax = High95CI, shape = Model.Name),
                    lwd = 1, position = position_dodge(width = 1/2),
                    fill = "WHITE") +
    coord_flip() + facet_wrap(~DV)+
    theme_bw() + xlab("") + ylab("Rescaled Coefficient (Reference Level: No Alliance)") + 
    theme(legend.position="bottom",
          legend.title=element_blank(),
          axis.text = element_text(size=12),
          text = element_text(size=12),
          plot.title = element_text(hjust = .5, size = 14, face = "bold"),
          strip.text.x = element_text(size = 12, color='white',
                                      angle=0),
          strip.background = element_rect(fill = "#525252", color='#525252'))
  
  return(p)   
} 

# coefficient plot for dynamic logit model
cq_coef_plot <- function(ModelResults, data, subvar = FALSE, vars = NULL) {
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  modelcoef = lapply(ModelResults, function(x) FUN = x$coefficients)
  modelse = lapply(ModelResults, function(x) FUN = x$se)
  varnames = lapply(ModelResults, function(x) FUN = names(x$coefficients))
  
  dfplot <- list()
  for (i in 1:length(ModelResults)){
    ad = 0.5 / modelse[[i]]
    modelcoef[[i]] = modelcoef[[i]]* ad
    modelse[[i]] = modelse[[i]]*ad 
    
    ylo95 =modelcoef[[i]] - qt(.975, nrow(data))*(modelse[[i]])
    yhi95 =modelcoef[[i]] + qt(.975, nrow(data))*(modelse[[i]])
    ylo90 =modelcoef[[i]] - qt(.95, nrow(data))*(modelse[[i]])
    yhi90 =modelcoef[[i]] + qt(.95, nrow(data))*(modelse[[i]])
    dfplot[[i]] = data.frame(varnames[[i]], modelcoef[[i]], modelse[[i]], ylo95, yhi95,
                             ylo90, yhi90,
                             Model = paste('Model',seq_along(modelcoef)[i]))
  }
  coefficientdata = do.call(rbind,dfplot)
  colnames(coefficientdata) <- c("Variables", "Coefficients", "S.E",
                                 "Low95CI", "High95CI", 
                                 "Low90CI", "High90CI",
                                 "Model.Name")
  
  df <- coefficientdata %>%
    dplyr::mutate(Variables = as.character(Variables)) %>%
    dplyr::arrange(Variables)
  if(subvar == TRUE){
    df <- df[grep(paste(vars, collapse = "|"), df$Variables),]
    #re-order the variables: if factor, the order of factor levels is used; if character, an alphabetical order ist used
    df$Variables <- factor(df$Variables, levels = vars)
  } 
  
  p = ggplot(df, aes(colour = Model.Name)) + 
    geom_hline(yintercept = 0, lty = 2) +
    geom_linerange(aes(x = Variables, ymin = Low90CI,
                       ymax = High90CI),
                   lwd = 2, position = position_dodge(width = 1/2)) +
    geom_pointrange(aes(x = Variables, y = Coefficients, ymin = Low95CI,
                        ymax = High95CI, shape = Model.Name),
                    lwd = 1, position = position_dodge(width = 1/2),
                    fill = "WHITE") +
    coord_flip() + 
    theme_bw() + xlab("") + ylab("Rescaled Coefficient") + 
    theme(legend.position="bottom",
          legend.title=element_blank(),
          axis.text = element_text(size=14),
          text = element_text(size=14),
          plot.title = element_text(hjust = .5, size = 16, face = "bold"))
  
  return(p)   
} 


### put them together
##make a plot for beta
amen_coef_plot_df <- function(mod, Vars, lab = "Formal Alliance"){
  vc <-  mod$VC %>% as.data.frame()
  vc$ve <- NULL
  beta_vc <- data.frame(beta = apply(vc, 2, mean),
                        se = apply(vc, 2, sd),
                        var = colnames(vc),
                        Modname = lab,
                        Variables = colnames(vc)) 
  
  beta <- mod$BETA %>%
    as.data.frame()
  df <- data.frame(beta = apply(beta, 2, mean),
                   se = apply(beta, 2, sd),
                   var = colnames(beta),
                   Modname = lab)
  df$Variables <- gsub('.','',gsub('dyad|node','',df$var),fixed=TRUE)
  df <- bind_rows(beta_vc,df)
  # df <- df[grep(paste(Vars, collapse = "|"), df$Variables),]
  df$Variables <- factor(df$Variables, levels = Vars)
 
  
  df <- df %>% 
    dplyr::arrange(Variables)%>%
    dplyr::mutate(ad = .5/se) %>%
    dplyr::mutate(beta = beta*ad,
                  se = se*ad) %>%
    dplyr::mutate(lo95 = beta - 1.96*se,
                  hi95 = beta +1.96*se,
                  lo90 = beta - 1.645*se,
                  hi90 = beta +1.645*se) %>%
    dplyr::mutate(color_sig = ifelse(lo95 < 0 & hi95 > 0, "#ef8a62", "#2c7bb6"))
  
  return(df)
  
}



scaled_coef_polr <- function(ModelResults, data, subvar = FALSE, vars = NULL) {
  library(MASS)
  library(dplyr)
  require(ggplot2)
  library(gridExtra)
  library(tidyr)
  library(stringr)
  library(purrr)
  # paste the coef togeter
  modcoef = Map(function(x, y) c(coef(x), y$zeta), x = ModelResults, y = ModelResults)
  modelcoef <- lapply(modcoef, function(x) FUN =round(x, digits = 4))  
  #get the sd
  modelse = lapply(ModelResults, function(x) FUN = sqrt(diag(vcov(x))))
  #get the names together               
  varnames = Map(function(x, y) c(names(x$coefficients), names(y$zeta)), x = ModelResults, y = ModelResults)
  
  dfplot <- list()
  for (i in 1:length(modelcoef)){
    ad = 0.5 / modelse[[i]]
    modelcoef[[i]] = modelcoef[[i]]* ad
    modelse[[i]] = modelse[[i]]*ad 
    
    ylo95 =modelcoef[[i]] - qt(.975, nrow(data))*(modelse[[i]])
    yhi95 =modelcoef[[i]] + qt(.975, nrow(data))*(modelse[[i]])
    ylo90 =modelcoef[[i]] - qt(.95, nrow(data))*(modelse[[i]])
    yhi90 =modelcoef[[i]] + qt(.95, nrow(data))*(modelse[[i]])
    dfplot[[i]] = data.frame(varnames[[i]], modelcoef[[i]], modelse[[i]], ylo95, yhi95,
                             ylo90, yhi90,
                             Model = paste('Model',seq_along(modelcoef)[i]))
  }
  
  coefficientdata = do.call(rbind,dfplot)
  colnames(coefficientdata) <- c("Variables", "Coefficients", "S.E",
                                 "Low95CI", "High95CI", 
                                 "Low90CI", "High90CI",
                                 "Model.Name")
  df <- coefficientdata %>%
    dplyr::mutate(Variables = as.character(Variables)) %>%
    dplyr::mutate(Variables = ifelse(Variables == "0|1", "cut1", Variables)) %>% 
    dplyr::mutate(Variables = ifelse(Variables == "1|2", "cut2", Variables)) %>% 
    dplyr::mutate(Variables = ifelse(Variables == "2|3", "cut3", Variables)) %>% 
    dplyr::mutate(Variables = ifelse(Variables == "3|4", "cut4", Variables)) %>% 
    dplyr::arrange(Variables) 
  
  
  
  if(subvar == TRUE){
    df <- df[grep(paste(vars, collapse = "|"), df$Variables),]
    #re-order the variables: if factor, the order of factor levels is used; if character, an alphabetical order ist used
    df$Variables <- factor(df$Variables, levels = vars)
  } 
  
  p = ggplot(df, aes(colour = Model.Name)) + 
    geom_hline(yintercept = 0, lty = 2) +
    geom_linerange(aes(x = Variables, ymin = Low90CI,
                       ymax = High90CI),
                   lwd = 2, position = position_dodge(width = 1/2)) +
    geom_pointrange(aes(x = Variables, y = Coefficients, ymin = Low95CI,
                        ymax = High95CI, shape = Model.Name),
                    lwd = 1, position = position_dodge(width = 1/2),
                    fill = "WHITE") +
    coord_flip() + 
    theme_bw() + xlab("") + ylab("Rescaled Coefficient") + 
    theme(legend.position="bottom",
          legend.title=element_blank(),
          axis.text = element_text(size=12),
          text = element_text(size=12),
          plot.title = element_text(hjust = .5, size = 14, face = "bold"))
  
  return(p)   
} 
