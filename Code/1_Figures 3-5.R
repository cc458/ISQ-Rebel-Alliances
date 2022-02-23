rm(list=ls())
# Function to load packages
loadPkg=function(toLoad){
  for(lib in toLoad){
    if(! lib %in% installed.packages()[,1])
    {install.packages(lib, repos='http://cran.rstudio.com/')}
    suppressMessages( library(lib, character.only=TRUE))}}

# Load libraries
packs=c("readr", "readxl", "dplyr", "countrycode", "ggplot2", "amen","statnet","VGAM",
        "amen", "purrr", "texreg", "abind", "ergm","reshape2",
        'parallel','foreach','doParallel',"haven", "readstata13")
loadPkg(packs)

#### set working directory to the folder named "Replication_ISQ"
#setwd("Replication_ISQ")
#####################

source("Code/helpFUN.R")
source("Code/function_interaction_plot.R")
###############################################################################################
### This is the code for the models in the main text: Figure 3-5
###############################################################################################

load("Data/glm_df4.RData")
glm_df <- glm_df4
rm(glm_df4)
## creade dyad anme
glm_df <- glm_df %>% 
  dplyr::mutate(dyad = paste(RBLside_A, RBLside_B, sep = "-"))

# create lag Y
glm_df <- glm_df %>% 
  dplyr::arrange(dyad, year) %>% 
  dplyr::group_by(dyad) %>% 
  dplyr::mutate(lag_y = dplyr::lag(alliance, n=1)) 

library(spduration)

#transform to spdur format
glm_df <- as.data.frame(glm_df)
glm_df <- add_duration(glm_df, "alliance", unitID = "dyad", tID ="year", freq = "year", ongoing = FALSE)

## create no-alliance years polynomial terms
glm_df <- glm_df  %>% 
  dplyr::mutate(t_yr = duration,
                t_yrsq = duration^2,
                t_yrcub = duration^3)


##############################################################################
###################  Figures 3      ###########################
##############################################################################

####  Figure 3 (a)

logit_f2 <- glm(alliance ~ co_ideology + allCOETH +splinter_indirect2+
                  rebels_count + pop_log + gdppc_log + milper + rugg_prop + post_cold + 
                  t_yr + t_yrsq + t_yrcub,
                data=glm_df, family='binomial'(link = "logit"))

logit_f3 <- glm(alliance ~ co_constituent +splinter_indirect2+
                  rebels_count + pop_log + gdppc_log + milper + rugg_prop + post_cold + 
                  t_yr + t_yrsq + t_yrcub,
                data=glm_df, family='binomial'(link = "logit"))

logit_f_coef <- scaled_coef_cluster(ModelResults = list(logit_f2), data = glm_df, 
                                    clusterid = "ccode", subvar = TRUE, vars = c(
                                      "(Intercept)",
                                      #"t_yr", "t_yrsq","t_yrcub",
                                      "post_cold", "pop_log","gdppc_log",
                                      "milper","rebels_count","rugg_prop", "splinter_indirect2", 
                                      "allCOETH","co_ideology"
                                      #,"co_constituent"
                                    ))

logit_f_coef <- logit_f_coef +  
  theme(legend.position = "none") +scale_x_discrete(labels=c("Intercept",
                                                             #  "t",
                                                             # expression(t^2),
                                                             #  expression(t^3),
                                                             'Post-Cold War',
                                                             "Population",
                                                             'GDP per capita',
                                                             'Military personnel',
                                                             "Number of rebel groups",
                                                             "Terrain ruggedness",
                                                             "Splinter",
                                                             "Co-ethnic",
                                                             "Co-ideological"#,
                                                             # "Co-constituent"
  )) + 
  ggtitle("")
ggsave("Figure/figure3_a.eps", plot = logit_f_coef, width = 6.5, height = 3.5, dpi = 400, units = "in")


## Figure 3(b)
###################################################
## To produce the simulation-based marginal effect, you need to use `postregplots`
# package from github
#library(devtools)
#install_github("cc458/postregplots")
#######################################################

library(postregplots)
logit_f2_mar1 <- density_ridgeplot(logit_f2, n.sim = 1000, varname = c('allCOETH'), 
                                   data = glm_df, val1 = 0, val2 = 1, clusterid = "ccode")
logit_f2_mar1 <-logit_f2_mar1 + xlab("Average of First Difference") + 
  scale_y_discrete(expand = c(-0.1, 0),labels=c("")) + 
  ggtitle("Co-ethnic") + 
  theme(axis.text = element_text(size=14),
        text = element_text(size=14),
        plot.title = element_text(size = 18))
mean(logit_f2_mar1$data$value) #0.07076735

## marginal for model2 and model 3
logit_f2_mar <- density_ridgeplot(logit_f2, n.sim = 1000, varname = c('co_ideology'), 
                                  data = glm_df, val1 = 0, val2 = 1, clusterid = "ccode")
mean(logit_f2_mar$data$value) #0.1615785
logit_f2_mar <- logit_f2_mar + xlab("Average of First Difference") + 
  ggtitle("Co-ideological") +
  scale_y_discrete(expand = c(-0.1, 0),labels=c("")) + 
  theme(axis.text = element_text(size=14),text = element_text(size=14),
        plot.title = element_text(size = 18))

library(grid)
library(gridExtra)

f <- grid.arrange(logit_f2_mar,logit_f2_mar1, ncol = 2)
ggsave("Figure/figure3_b.eps", plot = f, width = 6.5, height = 3.5, dpi = 400, units = "in")

## predicted
d = predict_probs(ModelResults = logit_f2, n.sim = 1000, varname = c('co_ideology'), 
                  data = glm_df, val1 = 0, val2 = 1,intervals=1, clusterid = "ccode")


##############################################################################
# Figure 4
##############################################################################

## for multinomial 
data <- glm_df

table(data$dv_nominal)
data$dv_nominal <- factor(data$dv_nominal, labels = c("noalliance", "informal", "formal"))
table(data$dv_nominal)
library(VGAM)

nominal_tfit2 <- vglm(dv_nominal ~ co_ideology + allCOETH + splinter_indirect2+
                        rebels_count + pop_log + gdppc_log + milper + rugg_prop + post_cold + 
                        t_yr + t_yrsq + t_yrcub,
                      data= data, multinomial(refLevel="noalliance"))
nominal_tfit3 <- vglm(dv_nominal ~ co_constituent + splinter_indirect2+
                        rebels_count + pop_log + gdppc_log + milper + rugg_prop + post_cold + 
                        t_yr + t_yrsq + t_yrcub,
                      data= data, multinomial(refLevel="noalliance"))

coef_mlogit_t <- mlogit_scaled_coef(fit = list(nominal_tfit2), subvar = TRUE, vars = c(
  "(Intercept)",
  #"t_yr", "t_yrsq","t_yrcub",
  "post_cold", "pop_log","gdppc_log",
  "milper","rebels_count","rugg_prop","splinter_indirect2", 
  "allCOETH","co_ideology"
  #, "co_constituent"
))

coef_mlogit_t <- coef_mlogit_t + theme(legend.position = "none")+ scale_x_discrete(labels=c("Intercept",
                                                                                            #  "t",
                                                                                            # expression(t^2),
                                                                                            #expression(t^3),
                                                                                            'Post-Cold War',
                                                                                            "Population",
                                                                                            'GDP per capita',
                                                                                            'Military personnel',
                                                                                            "Number of rebel groups",
                                                                                            "Terrain ruggedness",
                                                                                            "Splinter",
                                                                                            "Co-ethnic", 
                                                                                            "Co-ideological"
)) + 
  ggtitle("")

ggsave("Figure/figure4_a.eps", plot = coef_mlogit_t, 
       width = 6.5, height = 3.5, dpi = 400,units = "in")


## FIgure 4-b

p_ideology <- vglm_sim(mod = nominal_tfit2, var = "co_ideology", n_sim = 1000)

p_ideology$data %>% dplyr::group_by(type) %>% dplyr::summarise(mean(Pred)) 
#formal         0.104 
#informal       0.0540

p4 <- p_ideology + scale_y_discrete(expand = c(-0.1, 0),
                                    labels=c("Formal alliance","Informal alliance")) + 
  xlab("Average of First Difference") + 
  ggtitle("Co-ideological") + 
  theme(axis.text = element_text(size=12),
        text = element_text(size=12),
        plot.title = element_text(hjust = .5, size = 14))

p_allCOETH <- vglm_sim(mod = nominal_tfit2, var = "allCOETH", n_sim = 1000)
p_allCOETH$data %>% dplyr::group_by(type) %>% dplyr::summarise(mean(Pred)) 
#1 formal        -0.0230
#2 informal       0.103
p3 <- p_allCOETH + scale_y_discrete(expand = c(-0.1, 0),
                                    labels=c("Formal alliance","Informal alliance")) + 
  xlab("Average of First Difference") + 
  ggtitle("Co-ethnic") + 
  theme(axis.text = element_text(size=12),
        text = element_text(size=12),
        plot.title = element_text(hjust = .5, size = 14)) 
### put them together
f1<- grid.arrange(p4, p3, nrow = 1)
ggsave("Figure/figure4_b.eps", plot = f1, width = 6.5, height = 3.5, dpi = 400, units = "in")


##############################################################################
### Figure 5: AME
##############################################################################

## This will be running over 12 hours.
## Please run the R scrip named as "1_Running_AME.R" in the "Code folder"
##################

### formal alliance
load("Data/ame/Fit_AME47_15.RData")

ame_ally_m2 <- amen_coef_plot_df(mod = Fit_AME47_15$ally_ideo, 
                                 Vars = c("va",  "intercept", "post_cold",
                                          "pop_log","gdppc_log", "milper",
                                          "rebels_count",
                                          "rugg_prop",
                                          "splinter_indirect2",
                                          "allCOETH","co_ideology"),
                                 lab = "Alliance")

### formal alliance
ame_formal_m2 <- amen_coef_plot_df(mod = Fit_AME47_15$formal_ideo, 
                                   Vars = c("va",  "intercept", "post_cold",
                                            "pop_log","gdppc_log", "milper",
                                            "rebels_count",
                                              "rugg_prop",
                                              "splinter_indirect2",
                                              "allCOETH","co_ideology"),
                                   lab = "Formal alliance")

ame_df <- bind_rows(ame_ally_m2, ame_formal_m2)

### change the level 
ame_df$Variables <- factor(ame_df$Variables,
                                levels =  c("va",  "intercept", "post_cold",
                                            "pop_log","gdppc_log", "milper",
                                            "rebels_count",
                                            "rugg_prop",
                                            "splinter_indirect2",
                                            "allCOETH","co_ideology"))

ame_ally_p <-  ggplot(ame_df, aes(colour = Modname)) + 
  geom_hline(yintercept = 0, lty = 2) +
  geom_linerange(aes(x = Variables, ymin = lo90,
                     ymax = hi90),
                 lwd = 2, position = position_dodge(width = 1/2)) +
  geom_pointrange(aes(x = Variables, y = beta, ymin = lo95,
                      ymax = hi95, shape = Modname),
                  lwd = 1, position = position_dodge(width = 1/2),
                  fill = "WHITE") +
  coord_flip() + 
  theme_bw() + xlab("") + ylab("Rescaled Coefficient") + 
  theme(legend.position="bottom",
        legend.title=element_blank(),
        axis.text = element_text(size=10),
        text = element_text(size=10),
        plot.title = element_text(hjust = .5, size = 10, face = "bold"))+ 
  ggtitle("")

ame_ally_p <- ame_ally_p +  scale_x_discrete(labels=c("Sender Receiver Covariance", "Intercept", 
                                                      'Post-Cold War (Country-level)',
                                                      "Population (Country-level)",
                                                      'GDP per capita (Country-level)',
                                                      'Military personnel(Country-level)',
                                                      "Number of rebel groups (Country-level)",
                                                      "Terrain ruggedness (Country-level)",
                                                      "Splinter (Dyad-level)",
                                                      "Co-ethnic (Dyad-level)",
                                                      "Co-ideological (Dyad-level)")) 
ggsave("Figure/figure5.eps", plot = ame_ally_p, width = 6.5, 
       height = 4, dpi = 400, units = "in")

######## end of the code for figures in the main text##
## For robustness checks, go the file "2_Appendix_figures A1-A21.R"
