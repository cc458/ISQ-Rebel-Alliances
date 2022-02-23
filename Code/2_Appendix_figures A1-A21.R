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

##############################################################################
################### Online  Appendix   Figures   ###########################
##############################################################################

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
# A1: Robustness Check: Interaction Effects for Co-constituency Variables and Government Power
##############################################################################
logit_inter <- glm(alliance ~ co_ideology*gdppc_log + allCOETH*gdppc_log +splinter_indirect2+
                     rebels_count + pop_log + gdppc_log + milper + rugg_prop + post_cold + 
                     t_yr + t_yrsq + t_yrcub,
                   data=glm_df, family='binomial'(link = "logit"))

logit_inter2 <- glm(alliance ~ co_ideology*milper + milper*allCOETH +splinter_indirect2+
                      rebels_count + pop_log + gdppc_log + milper + rugg_prop + post_cold + 
                      t_yr + t_yrsq + t_yrcub,
                    data=glm_df, family='binomial'(link = "logit"))


library(postregplots)

logit_f_coef <- scaled_coef_cluster(ModelResults = list(logit_inter,logit_inter2), data = glm_df, 
                                    clusterid = "ccode",  subvar = TRUE, vars = c(
                                      "(Intercept)",
                                      #"t_yr", "t_yrsq","t_yrcub",
                                      "post_cold", "pop_log","gdppc_log",
                                      "milper","rebels_count","rugg_prop", "splinter_indirect2", 
                                      "allCOETH","co_ideology",
                                      "co_ideology:gdppc_log","co_ideology:milper",
                                      "gdppc_log:allCOETH","milper:allCOETH"
                                    ))

logit_f_coef_inter <- logit_f_coef +  scale_x_discrete(labels=c("Intercept",
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
                                                                "Co-ideological",
                                                                "Co-ideological x GDP per capita",
                                                                "Co-ideological x Military personnel",
                                                                "Co-ethnic x GDP per capita",
                                                                "Co-ethnic x Military personnel")) + 
  ggtitle("")
ggsave("Appendix_figures/A1_a.pdf", plot = logit_f_coef_inter, width = 12, height = 7.5)

## for multinomial 
data <- glm_df

table(data$dv_nominal)
data$dv_nominal <- factor(data$dv_nominal, labels = c("noalliance", "informal", "formal"))
table(data$dv_nominal)
library(VGAM)


nominal_inter1 <- vglm(dv_nominal ~ co_ideology*gdppc_log + allCOETH*gdppc_log + splinter_indirect2+
                         rebels_count + pop_log + gdppc_log + milper + rugg_prop + post_cold + 
                         t_yr + t_yrsq + t_yrcub,
                       data= data, multinomial(refLevel="noalliance"))

nominal_inter2 <- vglm(dv_nominal ~ co_ideology*milper + allCOETH*milper + splinter_indirect2+
                         rebels_count + pop_log + gdppc_log + milper + rugg_prop + post_cold + 
                         t_yr + t_yrsq + t_yrcub,
                       data= data, multinomial(refLevel="noalliance"))

coef_mlogit_inter <- mlogit_scaled_coef(fit = list(nominal_inter1, nominal_inter2), subvar = TRUE, vars = c(
  "(Intercept)",
  # "t_yr", "t_yrsq","t_yrcub",
  "post_cold", "pop_log","gdppc_log",
  "milper","rebels_count","rugg_prop","splinter_indirect2", 
  "allCOETH","co_ideology", "co_ideology:gdppc_log","co_ideology:milper",
  "gdppc_log:allCOETH","milper:allCOETH"
))

coef_mlogit_inter  <- coef_mlogit_inter  +  scale_x_discrete(labels=c("Intercept",
                                                                      #   "t",
                                                                      #  expression(t^2),
                                                                      #  expression(t^3),
                                                                      'Post-Cold War',
                                                                      "Population",
                                                                      'GDP per capita',
                                                                      'Military personnel',
                                                                      "Number of rebel groups",
                                                                      "Terrain ruggedness",
                                                                      "Splinter",
                                                                      "Co-ethnic", 
                                                                      "Co-ideological",
                                                                      "Co-ideological x GDP per capita",
                                                                      "Co-ideological x Military personnel",
                                                                      "Co-ethnic x GDP per capita",
                                                                      "Co-ethnic x Military personnel")) + 
  ggtitle("")
ggsave("Appendix_figures/A1_b.pdf", plot = coef_mlogit_inter, width = 12, height = 7.5)


##############################################################################
# A2: Robustness Check: Logistic Regression Results for Rebel Alliance (Controlling for strength of central command)
##############################################################################

logit_f2_1 <- glm(alliance ~ co_ideology + allCOETH +splinter_indirect2+strengthcent_dyad_dum+
                    rebels_count + pop_log + gdppc_log + milper + rugg_prop + post_cold + 
                    t_yr + t_yrsq + t_yrcub,
                  data=glm_df, family='binomial'(link = "logit"))

logit_f2_2 <- glm(alliance ~ co_ideology + allCOETH +splinter_indirect2+strengthcent_dyad+
                     rebels_count + pop_log + gdppc_log + milper + rugg_prop + post_cold + 
                     t_yr + t_yrsq + t_yrcub,
                   data=glm_df, family='binomial'(link = "logit"))


logit_f_coef <- scaled_coef_cluster(ModelResults = list(logit_f2_2 ), data = glm_df, 
                                    clusterid = "ccode", subvar = TRUE, vars = c(
                                      "(Intercept)","post_cold", "pop_log","gdppc_log",#"t_yr", "t_yrsq","t_yrcub",
                                      "milper","rebels_count","rugg_prop", "splinter_indirect2", #"centcontrol_dyad_dum",
                                      "strengthcent_dyad", #"strengthcent_dyad_dum",
                                      "allCOETH","co_ideology"#,"co_constituent"
                                    ))

logit_f_coef <- logit_f_coef + theme(legend.position = "none") +  scale_x_discrete(labels=c("Intercept",
                                                          #"Years since",
                                                          #expression(t^2),
                                                          #expression(t^3),
                                                          'Post-Cold War',
                                                          "Population",
                                                          'GDP per capita',
                                                          'Military personnel',
                                                          "Number of rebel groups",
                                                          "Terrain ruggedness",
                                                          "Splinter",
                                                          "Strength central command",
                                                          #"Strength of central command",
                                                          "Co-ethnic",
                                                          "Co-ideological"#,
                                                         # "Co-constituent"
                                                         )) + 
  ggtitle("")
ggsave("Appendix_figures/A2.pdf", plot = logit_f_coef, width = 12, height = 7.5)

##############################################################################
# A3: Robustness Check: Logistic Regression Results for Rebel Alliance (Controlling for strength of central command)
##############################################################################

nominal_tfit2 <- vglm(dv_nominal ~ co_ideology + allCOETH + splinter_indirect2+strengthcent_dyad_dum+
                        rebels_count + pop_log + gdppc_log + milper + rugg_prop + post_cold + 
                        t_yr + t_yrsq + t_yrcub,
                      data= data, multinomial(refLevel="noalliance"))

nominal_tfit2_2 <- vglm(dv_nominal ~ co_ideology + allCOETH + splinter_indirect2+strengthcent_dyad+
                          rebels_count + pop_log + gdppc_log + milper + rugg_prop + post_cold + 
                          t_yr + t_yrsq + t_yrcub,
                        data= data, multinomial(refLevel="noalliance"))

coef_mlogit_t <- mlogit_scaled_coef(fit = list(nominal_tfit2_2), subvar = TRUE, vars = c(
  "(Intercept)",
  # "t_yr", "t_yrsq","t_yrcub",
  "post_cold", "pop_log","gdppc_log",
  "milper","rebels_count","rugg_prop","splinter_indirect2", 
  "strengthcent_dyad",
  #"strengthcent_dyad_dum",
  "allCOETH","co_ideology"#, "co_constituent"
))

coef_mlogit_t <- coef_mlogit_t +  theme(legend.position = "none") + scale_x_discrete(labels=c("Intercept",
                                                            # "t",
                                                            #  expression(t^2),
                                                            #  expression(t^3),
                                                            'Post-Cold War',
                                                            "Population",
                                                            'GDP per capita',
                                                            'Military personnel',
                                                            "Number of rebel groups",
                                                            "Terrain ruggedness",
                                                            "Splinter",
                                                            "Strength central command",
                                                          #  "Strength of central command",
                                                            "Co-ethnic", 
                                                            "Co-ideological")) + 
  ggtitle("")
ggsave("Appendix_figures/A3.pdf", plot = coef_mlogit_t, width = 12, height = 7.5)


##############################################################################
# A4: Robustness Check: Ordered Logistic Regressions on Rebel Strength of central command (1946-2013)
##############################################################################
df <- glm_df

df <- df %>% 
  dplyr::select(year, RBLside_A, rbl_aname,RBLside_B, rbl_bname, ccode,
                rebels_count, pop_log, gdppc_log, milper, rugg_prop ,post_cold,
                islamista, islamistb, RBLside_A_MarxistIdeol, RBLside_B_MarxistIdeol,
                RBLside_A_ethn, RBLside_B_ethn,RBL_B_strengthcent,RBL_A_strengthcent,
                RBL_B_centcontrol, RBL_A_centcontrol
  )

df <- df %>% 
  dplyr::mutate(RBLside_A_MarxistIdeol = as.integer(RBLside_A_MarxistIdeol),
                RBLside_B_MarxistIdeol = as.integer(RBLside_B_MarxistIdeol),
                RBL_A_centcontrol = ifelse(RBL_A_centcontrol == "no", 0,  RBL_A_centcontrol),
                RBL_B_centcontrol = ifelse(RBL_B_centcontrol == "no", 0,  RBL_B_centcontrol),
                RBL_B_strengthcent = ifelse(RBL_B_strengthcent == "low", 1, RBL_B_strengthcent),
                RBL_A_strengthcent = ifelse(RBL_A_strengthcent == "low", 1, RBL_A_strengthcent)) %>% 
  dplyr::mutate(RBL_A_centcontrol = ifelse(RBL_A_centcontrol == "yes", 1,  RBL_A_centcontrol),
                RBL_B_centcontrol = ifelse(RBL_B_centcontrol == "yes", 1,  RBL_B_centcontrol),
                RBL_B_strengthcent = ifelse(RBL_B_strengthcent == "moderate", 2, RBL_B_strengthcent),
                RBL_A_strengthcent = ifelse(RBL_A_strengthcent == "moderate", 2, RBL_A_strengthcent)) %>% 
  dplyr::mutate(RBL_B_strengthcent = ifelse(RBL_B_strengthcent == "high", 3, RBL_B_strengthcent),
                RBL_A_strengthcent = ifelse(RBL_A_strengthcent == "high", 3, RBL_A_strengthcent))

names(df) 
### collaps to group level
df_a <- df %>% 
  dplyr::select(year, RBLside_A, rbl_aname, ccode,
                rebels_count, pop_log, gdppc_log, milper, rugg_prop ,post_cold,
                islamista, RBLside_A_MarxistIdeol, 
                RBLside_A_ethn, RBL_A_strengthcent,
                RBL_A_centcontrol)
df_b <- df %>% 
  dplyr::select(year, RBLside_B, rbl_bname, ccode,
                rebels_count, pop_log, gdppc_log, milper, rugg_prop ,post_cold,
                islamistb, RBLside_B_MarxistIdeol,
                RBLside_B_ethn,RBL_B_strengthcent,
                RBL_B_centcontrol )

names(df_a) <- c("year", "RBLside", "rbl_name", "ccode",
                 "rebels_count", "pop_log", "gdppc_log", "milper", "rugg_prop" ,"post_cold",
                 "islamist", "RBLside_MarxistIdeol",
                 "RBLside_ethn","RBL_strengthcent",
                 "RBL_centcontrol")
names(df_b) <- c("year", "RBLside", "rbl_name", "ccode",
                 "rebels_count", "pop_log", "gdppc_log", "milper", "rugg_prop" ,"post_cold",
                 "islamist", "RBLside_MarxistIdeol",
                 "RBLside_ethn","RBL_strengthcent",
                 "RBL_centcontrol")
## aggrgate
df_group <- bind_rows(df_a, df_b)
dim(df_group )#[1] 4992   15
##get a uniqe
df_group <- unique(df_group)

## rebel strength only goes to 2011
df_group$RBL_centcontrol <- as.integer(df_group$RBL_centcontrol)
df_group <- df_group %>%  
  dplyr::filter(year < 2012) %>% 
  dplyr::mutate(RBL_strengthcent = ifelse(is.na(RBL_strengthcent), 0, RBL_strengthcent),
                RBL_centcontrol = ifelse(is.na(RBL_centcontrol), 0, RBL_centcontrol)) %>% 
  dplyr::mutate(RBL_strengthcent_dum = ifelse(RBL_strengthcent>=1, 1, 0)) %>% 
  dplyr::arrange(RBLside, year) 

library(MASS)

df_group$RBL_strengthcent = as.factor(df_group$RBL_strengthcent)

m1 <- polr(RBL_strengthcent ~ islamist + RBLside_MarxistIdeol + RBLside_ethn +
             rebels_count + pop_log + gdppc_log + milper + rugg_prop + post_cold,
           data = df_group, method="logistic", Hess=TRUE)

### coef plot
coef <-scaled_coef_polr(ModelResults = list(m1), data = df_group, subvar = TRUE, 
                        vars = c('cut1','cut2','cut3',"post_cold", "pop_log",
                                 "gdppc_log",
                                 "milper","rebels_count","rugg_prop",
                                 "RBLside_ethn",
                                 "RBLside_MarxistIdeol","islamist"))
coef  <- coef  + theme(legend.position = "none") +  scale_x_discrete(labels=c('Cut1','Cut2','Cut3',
                                            'Post-Cold War',
                                            "Population",
                                            'GDP per capita',
                                            'Military personnel',
                                            "Number of rebel groups",
                                            "Terrain ruggedness",
                                            "Ethnic organization",
                                            "Marxist organization",
                                            "Islamist organization")) +  ggtitle("")

ggsave("Appendix_figures/A4.pdf", plot = coef , width = 12, height = 7.5)


##############################################################################
# A5: Robustness Check: Logistic Regressions on Rebel Alliance Onset (1946-2015)
##############################################################################
# create lag Y
glm_df2 <- glm_df %>% 
  dplyr::group_by(dyad) %>% 
  dplyr::arrange(dyad, year) %>% 
  dplyr::select(dyad, year, alliance, failure, ongoing, end.spell, cured, atrisk, censor,duration, everything()) %>% 
  dplyr::group_by(dyad) %>% 
  dplyr::arrange(dyad, year) %>% 
  dplyr::mutate(sumcun = cumsum(alliance)) %>% 
  dplyr::select(dyad, year, alliance, sumcun, failure, ongoing, 
                end.spell, cured, atrisk, censor,duration, everything()) %>% 
  dplyr::filter(alliance==sumcun)

logit_f2 <- glm(alliance ~ co_ideology + allCOETH +splinter_indirect2+
                  rebels_count + pop_log + gdppc_log + milper + rugg_prop + post_cold + 
                  t_yr + t_yrsq + t_yrcub,
                data=glm_df2, family='binomial'(link = "logit"))

logit_f_coef <- scaled_coef_cluster(ModelResults = list(logit_f2), data = glm_df2, 
                                    clusterid = "ccode", subvar = TRUE, vars = c(
                                      "(Intercept)",
                                      #"t_yr", "t_yrsq","t_yrcub",
                                      "post_cold", "pop_log","gdppc_log",
                                      "milper","rebels_count","rugg_prop", "splinter_indirect2", 
                                      "allCOETH","co_ideology"
                                    ))

logit_f_coef <- logit_f_coef + theme(legend.position = "none") +  scale_x_discrete(labels=c("Intercept",
                                                                                            #  "t",
                                                                                            #  expression(t^2),
                                                                                            #  expression(t^3),
                                                                                            'Post-Cold War',
                                                                                            "Population",
                                                                                            'GDP per capita',
                                                                                            'Military personnel',
                                                                                            "Number of rebel groups",
                                                                                            "Terrain ruggedness",
                                                                                            "Splinter",
                                                                                            "Co-ethnic",
                                                                                            "Co-ideological")) + 
  ggtitle("")
ggsave("Appendix_figures/A5.pdf", plot = logit_f_coef, width = 12, height = 7.5)


##############################################################################
# A6: Robustness Check: Multinominal Logit Regressions on Rebel Alliance Onset (1946-2015)
##############################################################################

## for multinomial 
data2 <- glm_df2

table(data2$dv_nominal)
data2$dv_nominal <- factor(data2$dv_nominal, labels = c("noalliance", "informal", "formal"))
table(data2$dv_nominal)
library(VGAM)

nominal_tfit2 <- vglm(dv_nominal ~ co_ideology + allCOETH + splinter_indirect2+
                        rebels_count + pop_log + gdppc_log + milper + rugg_prop + post_cold + 
                        t_yr + t_yrsq + t_yrcub,
                      data= data2, multinomial(refLevel="noalliance"))

coef_mlogit_t <- mlogit_scaled_coef(fit = list(nominal_tfit2), subvar = TRUE, vars = c(
  "(Intercept)",
  #"t_yr", "t_yrsq","t_yrcub",
  "post_cold", "pop_log","gdppc_log",
  "milper","rebels_count","rugg_prop","splinter_indirect2", 
  "allCOETH","co_ideology"
))

coef_mlogit_t <- coef_mlogit_t + theme(legend.position = "none") +
  scale_x_discrete(labels=c("Intercept",
                            # "t",
                            #  expression(t^2),
                            #  expression(t^3),
                            'Post-Cold War',
                            "Population",
                            'GDP per capita',
                            'Military personnel',
                            "Number of rebel groups",
                            "Terrain ruggedness",
                            "Splinter",
                            "Co-ethnic", 
                            "Co-ideological")) + 
  ggtitle("")
ggsave("Appendix_figures/A6.pdf", plot = coef_mlogit_t, width = 12, height = 7.5)



##############################################################################
# A7: Robustness Check: Estimated Coefficients for Alliance from Split-Population Duration Models (1946-2015)
##############################################################################

### SPDM
source("Code/spdur_figure_setup.R")

#transform to spdur format
spdm_df_ally <- as.data.frame(glm_df)


ally_spdur_ideo <- spdur(
  duration ~ co_ideology + allCOETH +splinter_indirect2 + 
    rebels_count + pop_log + gdppc_log + milper + rugg_prop + post_cold,
  atrisk ~ splinter_indirect2 + rebels_count + pop_log + gdppc_log + milper + rugg_prop + post_cold,
  data = spdm_df_ally, distr = "loglog", silent = TRUE)
summary(ally_spdur_ideo)

coef_p2 <- plot_spdm_coef(x = ally_spdur_ideo, subvar = TRUE, variables = c("(Intercept)","post_cold",
                                                                            "pop_log","gdppc_log",
                                                                            "milper", "rebels_count",
                                                                            "rugg_prop","splinter_indirect2",
                                                                            "co_ideology","allCOETH"))
coef_p2_dur <- coef_p2[[1]] + scale_x_discrete(labels=c("Intercept",
                                                        'Post-Cold War',
                                                        "Population",
                                                        'GDP per capita',
                                                        'Military personnel',
                                                        "Number of rebel groups",
                                                        "Terrain ruggedness",
                                                        "Splinter",
                                                        "Co-ideological",
                                                        "Co-ethnic")) + 
  ggtitle("Time to Alliance") #duration 
coef_p2_risk <- coef_p2[[2]] + scale_x_discrete(labels=c("Intercept",
                                                         'Post-Cold War',
                                                         "Population",
                                                         'GDP per capita',
                                                         'Military personnel',
                                                         "Number of rebel groups",
                                                         "Terrain ruggedness",
                                                         "Splinter")) + 
  ggtitle("Alliance Risk")  #risk
coef_p2_shap <- coef_p2[[3]] + scale_x_discrete(labels=expression(log(alpha))) #shape

pdf("Appendix_figures/A7.pdf",width = 12, height = 6.5)
grid.arrange(arrangeGrob(coef_p2_risk,coef_p2_shap, ncol = 1,heights = c(6,1)), coef_p2_dur, ncol = 2)
dev.off() 

##############################################################################
# A8: Robustness Check: Estimated Coefficients for Formal Alliance from Split-Population Duration Models (1946-2015)
##############################################################################
### look at formal alliances

fspdm_df <- glm_df %>% 
  dplyr::select(-failure, -ongoing, -end.spell, -cured,
                -atrisk, -censor, -duration, -t.0)

#transform to spdur format
fspdm_df <- as.data.frame(fspdm_df)
fspdm_df_ally <- add_duration(fspdm_df, "FormalAlliance", unitID = "dyad", tID ="year", freq = "year", ongoing = FALSE)
names(fspdm_df_ally)

f_ally_spdur_ideo <- spdur(
  duration ~ co_ideology + allCOETH + splinter_indirect2 + 
    rebels_count + pop_log + gdppc_log + milper + rugg_prop + post_cold,
  atrisk ~  splinter_indirect2 + rebels_count + pop_log + gdppc_log + milper + rugg_prop + post_cold,
  data = fspdm_df_ally, distr = "loglog", silent = TRUE)
summary(f_ally_spdur_ideo)

f_coef_p2 <- plot_spdm_coef(x = f_ally_spdur_ideo, subvar = TRUE, variables = c("(Intercept)","post_cold",
                                                                                "pop_log","gdppc_log",
                                                                                "milper", "rebels_count",
                                                                                "rugg_prop","splinter_indirect2", 
                                                                                "co_ideology","allCOETH"))
f_coef_p2_dur <- f_coef_p2[[1]] + scale_x_discrete(labels=c("Intercept",
                                                            'Post-Cold War',
                                                            "Population",
                                                            'GDP per capita',
                                                            'Military personnel',
                                                            "Number of rebel groups",
                                                            "Terrain ruggedness",
                                                            "Splinter",
                                                            "Co-ideological",
                                                            "Co-ethnic")) + 
  ggtitle("Time to Formal Alliance") #duration 
f_coef_p2_risk <- f_coef_p2[[2]] + scale_x_discrete(labels=c("Intercept",
                                                             'Post-Cold War',
                                                             "Population",
                                                             'GDP per capita',
                                                             'Military personnel',
                                                             "Number of rebel groups",
                                                             "Terrain ruggedness",
                                                             "Splinter")) + 
  ggtitle("Formal Alliance Risk")  #risk
f_coef_p2_shap <- f_coef_p2[[3]] + scale_x_discrete(labels=expression(log(alpha))) #shape

pdf("Appendix_figures/A8.pdf",width = 12, height = 6.5)
grid.arrange(arrangeGrob(f_coef_p2_risk,f_coef_p2_shap, ncol = 1,heights = c(6,1)), f_coef_p2_dur, ncol = 2)
dev.off() 



##############################################################################
# A9: Robustness Check: Logistic Regression Results for Alliance Controlling for State Co-Sponsor (1975-2009)
##############################################################################
logit_r2 <- glm(alliance ~ co_ideology + allCOETH +splinter_indirect2+state_cosponsor_dummy + 
                  rebels_count + pop_log + gdppc_log + milper + rugg_prop + post_cold + 
                  t_yr + t_yrsq + t_yrcub,
                data=glm_df, family='binomial'(link = "logit"))

logit_sponsor <- scaled_coef_cluster(ModelResults = list(logit_r2), data = glm_df, 
                                     clusterid = "ccode", subvar = TRUE, vars = c(
                                       "(Intercept)",
                                       #"t_yr", "t_yrsq","t_yrcub",
                                       "post_cold", "pop_log","gdppc_log",
                                       "milper","rebels_count","rugg_prop", "state_cosponsor_dummy",
                                       "splinter_indirect2", 
                                       "allCOETH","co_ideology"
                                     ))

logit_sponsor <- logit_sponsor + theme(legend.position = "none")+ scale_x_discrete(labels=c("Intercept",
                                                            # "t",
                                                            #  expression(t^2),
                                                            # expression(t^3),
                                                            'Post-Cold War',
                                                            "Population",
                                                            'GDP per capita',
                                                            'Military personnel',
                                                            "Number of rebel groups",
                                                            "Terrain ruggedness",
                                                            "State co-sponsor",
                                                            "Splinter",
                                                            "Co-ethnic",
                                                            "Co-ideological")) + 
  ggtitle("")
ggsave("Appendix_figures/A9_a.pdf", plot = logit_sponsor, width = 12, height = 7.5)


glm_df3 <- glm_df %>% 
  dplyr::filter(year >= 1975 & year < 2010)


logit_rr2 <- glm(alliance ~ co_ideology + allCOETH +splinter_indirect2+
                   rebels_count + pop_log + gdppc_log + milper + rugg_prop + post_cold + 
                   t_yr + t_yrsq + t_yrcub,
                 data=glm_df3, family='binomial'(link = "logit"))

logit_nosponsor <- scaled_coef_cluster(ModelResults = list(logit_rr2), data = glm_df3, 
                                       clusterid = "ccode", subvar = TRUE, vars = c(
                                         "(Intercept)",
                                         # "t_yr", "t_yrsq","t_yrcub",
                                         "post_cold", "pop_log","gdppc_log",
                                         "milper","rebels_count","rugg_prop", 
                                         "splinter_indirect2", 
                                         "allCOETH","co_ideology"
                                       ))

logit_nosponsor <- logit_nosponsor + theme(legend.position = "none")+  scale_x_discrete(labels=c("Intercept",
                                                                # "t",
                                                                #  expression(t^2),
                                                                # expression(t^3),
                                                                'Post-Cold War',
                                                                "Population",
                                                                'GDP per capita',
                                                                'Military personnel',
                                                                "Number of rebel groups",
                                                                "Terrain ruggedness",
                                                                "Splinter",
                                                                "Co-ethnic",
                                                                "Co-ideological")) + 
  ggtitle("")
ggsave("Appendix_figures/A9_b.pdf", plot = logit_nosponsor, width = 12, height = 7.5)



##############################################################################
# A10: Robustness Check: Multinominal Logit Regressions Results on Rebel Alliance with State Co-sponsor (1975-2009)
##############################################################################

data3 <- data %>% 
  dplyr::filter(year >= 1975 & year < 2010)

nominal_rr1 <- vglm(dv_nominal ~ co_ideology + allCOETH + splinter_indirect2+state_cosponsor_dummy + 
                     rebels_count + pop_log + gdppc_log + milper + rugg_prop + post_cold + 
                     t_yr + t_yrsq + t_yrcub,
                   data= data, multinomial(refLevel="noalliance"))

nominal_rr2 <- vglm(dv_nominal ~ co_ideology + allCOETH + splinter_indirect2+ 
                      rebels_count + pop_log + gdppc_log + milper + rugg_prop + post_cold + 
                      t_yr + t_yrsq + t_yrcub,
                    data= data3, multinomial(refLevel="noalliance"))


mlogit_sponsor <- mlogit_scaled_coef(fit = list(nominal_rr1), subvar = TRUE, vars = c(
  "(Intercept)",
  # "t_yr", "t_yrsq","t_yrcub",
  "post_cold", "pop_log","gdppc_log",
  "milper","rebels_count","rugg_prop","splinter_indirect2", "state_cosponsor_dummy",
  "allCOETH","co_ideology"
))

mlogit_sponsor <- mlogit_sponsor + theme(legend.position = "none")+ scale_x_discrete(labels=c("Intercept",
                                                              # "t",
                                                              #  expression(t^2),
                                                              # expression(t^3),
                                                              'Post-Cold War',
                                                              "Population",
                                                              'GDP per capita',
                                                              'Military personnel',
                                                              "Number of rebel groups",
                                                              "Terrain ruggedness",
                                                              "Splinter",
                                                              "State co-sponsor",
                                                              "Co-ethnic", 
                                                              "Co-ideological")) + 
  ggtitle("")
ggsave("Appendix_figures/A10_a.pdf", plot = mlogit_sponsor, width = 12, height = 7.5)

mlogit_sponsor <- mlogit_scaled_coef(fit = list(nominal_rr2), subvar = TRUE, vars = c(
  "(Intercept)",
  # "t_yr", "t_yrsq","t_yrcub",
  "post_cold", "pop_log","gdppc_log",
  "milper","rebels_count","rugg_prop","splinter_indirect2", 
  "allCOETH","co_ideology"
))

mlogit_sponsor <- mlogit_sponsor + theme(legend.position = "none")+ scale_x_discrete(labels=c("Intercept",
                                                                                              # "t",
                                                                                              #  expression(t^2),
                                                                                              # expression(t^3),
                                                                                              'Post-Cold War',
                                                                                              "Population",
                                                                                              'GDP per capita',
                                                                                              'Military personnel',
                                                                                              "Number of rebel groups",
                                                                                              "Terrain ruggedness",
                                                                                              "Splinter",
                                                                                              "Co-ethnic", 
                                                                                              "Co-ideological")) + 
  ggtitle("")
ggsave("Appendix_figures/A10_b.pdf", plot = mlogit_sponsor, width = 12, height = 7.5)


##############################################################################
# A11: Robustness Check: Estimated Coefficients from the AME Models with State Co-sponsor (1975-2009)
##############################################################################
load("Data/ame/Fit_AME7909.RData")

ame79_alliance_m2 <- amen_coef_plot_df(mod = Fit_AME7909$AllyList, 
                                       Vars = c("va",  "intercept", "post_cold",
                                                "pop_log","gdppc_log", "milper",
                                                "rebels_count",
                                                "rugg_prop",
                                                "splinter_indirect2",
                                                "state_cosponsor_dummy",
                                                  "allCOETH","co_ideology"),
                                       lab = "Alliance")


ame79_formal_m2 <- amen_coef_plot_df(mod = Fit_AME7909$Formal_AllyList, 
                                     Vars = c("va",  "intercept", "post_cold",
                                              "pop_log","gdppc_log", "milper",
                                              "rebels_count",
                                              "rugg_prop",
                                              "splinter_indirect2",
                                              "state_cosponsor_dummy",
                                              "allCOETH","co_ideology"),
                                     lab = "Formal alliance")
ame_df_r <- bind_rows(ame79_alliance_m2, ame79_formal_m2)


### change the level 
ame_df_r$Variables <- factor(ame_df_r$Variables,
                           levels =  c("va",  "intercept", "post_cold",
                                       "pop_log","gdppc_log", "milper",
                                       "rebels_count",
                                       "rugg_prop",
                                       "splinter_indirect2",
                                       "state_cosponsor_dummy",
                                       "allCOETH","co_ideology"))



ame_ally_p_r <-  ggplot(ame_df_r, aes(colour = Modname)) + 
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
        axis.text = element_text(size=14),
        text = element_text(size=14),
        plot.title = element_text(hjust = .5, size = 16, face = "bold"))+ 
  ggtitle("")

ame_ally_p_r <-ame_ally_p_r +  scale_x_discrete(labels=c("Sender Receiver Covariance", "Intercept", 
                                                         'Post-Cold War (Country-level)',
                                                         "Population (Country-level)",
                                                         'GDP per capita (Country-level)',
                                                         'Military personnel(Country-level)',
                                                         "Number of rebel groups (Country-level)",
                                                         "Terrain ruggedness (Country-level)",
                                                         "Splinter (Dyad-level)",
                                                      "State co-sponsor (Dyad-level)",
                                                      "Co-ethnic (Dyad-level)",
                                                      "Co-ideological (Dyad-level)")) 
ggsave("Appendix_figures/A11.pdf", plot = ame_ally_p_r, width = 11, height = 5.5)



##############################################################################
# A12: Robustness Check: Logistic Regressions on Alliance with imbalance of power (1946-2013)
##############################################################################
logit_r2 <- glm(alliance ~ co_ideology + allCOETH +splinter_indirect2+powerimbalance+ 
                  rebels_count + pop_log + gdppc_log + milper + rugg_prop + post_cold + 
                  t_yr + t_yrsq + t_yrcub,
                data=glm_df, family='binomial'(link = "logit"))

logit_r_power <- scaled_coef_cluster(ModelResults = list(logit_r2), data = glm_df, 
                                     clusterid = "ccode", subvar = TRUE, vars = c(
                                       "(Intercept)",
                                       # "t_yr", "t_yrsq","t_yrcub",
                                       "post_cold", "pop_log","gdppc_log",
                                       "milper","rebels_count","rugg_prop", "splinter_indirect2", "powerimbalance",
                                       "allCOETH", "co_ideology"
                                     ))

logit_r_power <- logit_r_power + theme(legend.position = "none")+ scale_x_discrete(labels=c("Intercept",
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
                                                                                            "Rebel imbalance of power",
                                                                                            "Co-ethnic",
                                                                                            "Co-ideological")) + 
  ggtitle("")
ggsave("Appendix_figures/A12.pdf", plot = logit_r_power, width = 12, height = 7.5)


##############################################################################
# A13: Robustness Check: Multinomial Logit Regressions on Alliance with imbalance of power (1946-2013)
##############################################################################
nominal_r2 <- vglm(dv_nominal ~ co_ideology + allCOETH + splinter_indirect2 + powerimbalance+ rebels_count + pop_log + 
                     gdppc_log + milper + rugg_prop + post_cold+ 
                     t_yr + t_yrsq + t_yrcub,
                   data= data, multinomial(refLevel="noalliance"))
coef_mlogit_power <- mlogit_scaled_coef(fit = list(nominal_r2), subvar = TRUE, vars = c(
  "(Intercept)",
  # "t_yr", "t_yrsq","t_yrcub",
  "post_cold", "pop_log","gdppc_log",
  "milper","rebels_count","rugg_prop", "splinter_indirect2","powerimbalance", 
  "allCOETH", "co_ideology"
))

coef_mlogit_power <- coef_mlogit_power + theme(legend.position = "none")+ scale_x_discrete(labels=c("Intercept",
                                                                                                    #   "t",
                                                                                                    #  expression(t^2),
                                                                                                    # expression(t^3),
                                                                                                    'Post-Cold War',
                                                                                                    "Population",
                                                                                                    'GDP per capita',
                                                                                                    'Military personnel',
                                                                                                    "Number of rebel groups",
                                                                                                    "Terrain ruggedness",
                                                                                                    "Splinter",
                                                                                                    "Rebel imbalance of power",
                                                                                                    "Co-ethnic",
                                                                                                    "Co-ideological")) + 
  ggtitle("")
ggsave("Appendix_figures/A13.pdf", plot = coef_mlogit_power, width = 12, height = 7.5)


##############################################################################
# A14: Robustness Check: Logistic Regression on Alliance with Mergers Variable (1946-2013)
##############################################################################

logit_rr2 <- glm(alliance ~ co_ideology + allCOETH +splinter_indirect2 +forge_merger_dummy+
                   rebels_count + pop_log + gdppc_log + milper + rugg_prop + post_cold + 
                   t_yr + t_yrsq + t_yrcub,
                 data=glm_df, family='binomial'(link = "logit"))

logit_rr_merger <- scaled_coef_cluster(ModelResults = list(logit_rr2), data = glm_df, 
                                       clusterid = "ccode", subvar = TRUE, vars = c(
                                         "(Intercept)",
                                         # "t_yr", "t_yrsq","t_yrcub",
                                         "post_cold", "pop_log","gdppc_log",
                                         "milper","rebels_count","rugg_prop", "forge_merger_dummy","splinter_indirect2",
                                         "allCOETH", "co_ideology"
                                       ))


logit_rr_merger <- logit_rr_merger +  theme(legend.position = "none")+ scale_x_discrete(labels=c("Intercept",
                                                                # "t",
                                                                #  expression(t^2),
                                                                # expression(t^3),
                                                                'Post-Cold War',
                                                                "Population",
                                                                'GDP per capita',
                                                                'Military personnel',
                                                                "Number of rebel groups",
                                                                "Terrain ruggedness",
                                                                "Merger organizations",
                                                                "Splinter",
                                                                "Co-ethnic",
                                                                "Co-ideological")) + 
  ggtitle("")
ggsave("Appendix_figures/A14.pdf", plot = logit_rr_merger, width = 12, height = 7.5)


##############################################################################
# A15: Robustness Check: Multinomial  Logistic Regressions on Alliance with Mergers Variable (1946-2013)
##############################################################################
nominal_rr2 <- vglm(dv_nominal ~ co_ideology + allCOETH +splinter_indirect2 +forge_merger_dummy+
                      rebels_count + pop_log + gdppc_log + milper + rugg_prop + post_cold + 
                      t_yr + t_yrsq + t_yrcub,
                    data= data, multinomial(refLevel="noalliance"))

####### coefficient plots for multinomial models

coef_mlogit_merge <- mlogit_scaled_coef(fit = list(nominal_rr2), subvar = TRUE, vars = c(
  "(Intercept)",
  # "t_yr", "t_yrsq","t_yrcub", 
  "post_cold", "pop_log","gdppc_log",
  "milper","rebels_count","rugg_prop", "forge_merger_dummy","splinter_indirect2",  
  "allCOETH","co_ideology"
))

coef_mlogit_merge <- coef_mlogit_merge + theme(legend.position = "none")+ scale_x_discrete(labels=c("Intercept",
                                                                    #   "t",
                                                                    #  expression(t^2),
                                                                    #   expression(t^3),
                                                                    'Post-Cold War',
                                                                    "Population",
                                                                    'GDP per capita',
                                                                    'Military personnel',
                                                                    "Number of rebel groups",
                                                                    "Terrain ruggedness",
                                                                    "Merger organizations",
                                                                    "Splinter",
                                                                    "Co-ethnic",
                                                                    "Co-ideological")) + 
  ggtitle("")
ggsave("Appendix_figures/A15.pdf", plot = coef_mlogit_merge, width = 12, height = 7.5)




##############################################################################
# A16: Robustness Check: Logit and Multinomial  Logistic Regressions on Alliance, Breaking Down Co-ideological Constituencies into Sub-types (1946-2013)
##############################################################################

logit_f1 <- glm(alliance ~ MarxistIdeol + islamist + allCOETH + splinter_indirect2+
                  rebels_count + pop_log + gdppc_log + milper + rugg_prop + post_cold + 
                  t_yr + t_yrsq + t_yrcub,
                data=glm_df, family='binomial'(link = "logit"))
logit_m1 <- scaled_coef_cluster(ModelResults = list(logit_f1), data = glm_df, 
                                clusterid = "ccode", subvar = TRUE, vars = c(
                                  "(Intercept)",
                                  #"t_yr", "t_yrsq","t_yrcub",
                                  "post_cold", "pop_log","gdppc_log",
                                  "milper","rebels_count","rugg_prop", 
                                  "splinter_indirect2",
                                  "allCOETH", "islamist", "MarxistIdeol"
                                ))


logit_m1<- logit_m1 +  scale_x_discrete(labels=c("Intercept",
                                                 #            "t",
                                                 #           expression(t^2),
                                                 #          expression(t^3),
                                                 'Post-Cold War',
                                                 "Population",
                                                 'GDP per capita',
                                                 'Military personnel',
                                                 "Number of rebel groups",
                                                 "Terrain ruggedness",
                                                 "Splinter",
                                                 "Co-ethnic",
                                                 "Co-Islamist",
                                                 "Co-Marxist")) + 
  ggtitle("")
ggsave("Appendix_figures/A16_a.pdf", plot = logit_m1, width = 12, height = 7.5)

nominal_tfit1 <- vglm(dv_nominal ~ MarxistIdeol + islamist + allCOETH + splinter_indirect2+
                        rebels_count + pop_log + gdppc_log + milper + rugg_prop + post_cold + 
                        t_yr + t_yrsq + t_yrcub,
                      data= data, multinomial(refLevel="noalliance"))


coef_mlogit_m1 <- mlogit_scaled_coef(fit = list(nominal_tfit1), subvar = TRUE, vars = c(
  "(Intercept)",
  # "t_yr", "t_yrsq","t_yrcub", 
  "post_cold", "pop_log","gdppc_log",
  "milper","rebels_count","rugg_prop", "splinter_indirect2",  
  "allCOETH", "islamist", "MarxistIdeol"
))

coef_mlogit_m1 <- coef_mlogit_m1 +  scale_x_discrete(labels=c("Intercept",
                                                              #     "t",
                                                              #    expression(t^2),
                                                              #   expression(t^3),
                                                              'Post-Cold War',
                                                              "Population",
                                                              'GDP per capita',
                                                              'Military personnel',
                                                              "Number of rebel groups",
                                                              "Terrain ruggedness",
                                                              "Splinter",
                                                              "Co-ethnic",
                                                              "Co-Islamist",
                                                              "Co-Marxist")) + 
  ggtitle("")
ggsave("Appendix_figures/A16_b.pdf", plot = coef_mlogit_m1, width = 12, height = 7.5)



##############################################################################
# A17: Robustness Check: Logistic Regression Results for Alliance Using an Alternative Shared Ideological Constituency Variable (1946-2015)
##############################################################################
glm_df <- glm_df %>% 
  dplyr::mutate(co_ethnic_nuance1 = ifelse(ideolnationlist ==1 & allCOETH ==1, 1, 0),
                co_ethnic_nuance2 = ifelse(allCOETH == 1 & (RBL_A_forge_ideolnat == 0 & RBL_B_forge_ideolnat == 0), 1, 0)) %>% 
  dplyr::mutate(co_ethnic_alternative = ifelse(co_ethnic_nuance1 == 1 | co_ethnic_nuance2 == 1, 1, 0))
table(glm_df$co_ethnic_nuance1)
table(glm_df$co_ethnic_nuance2)
table(glm_df$co_ethnic_alternative)


logit_f2_1 <- glm(alliance ~ co_ideology_new1 + allCOETH +splinter_indirect2+
                    rebels_count + pop_log + gdppc_log + milper + rugg_prop + post_cold + 
                    t_yr + t_yrsq + t_yrcub,
                  data=glm_df, family='binomial'(link = "logit"))


logit_f_coef <- scaled_coef_cluster(ModelResults = list(logit_f2_1 ), data = glm_df, 
                                    clusterid = "ccode", subvar = TRUE, vars = c(
                                      "(Intercept)",
                                      # "t_yr", "t_yrsq","t_yrcub",
                                      "post_cold", "pop_log","gdppc_log",
                                      "milper","rebels_count","rugg_prop", "splinter_indirect2", 
                                      "allCOETH","co_ideology_new1"
                                    ))

logit_f_coef <- logit_f_coef +  theme(legend.position = "none")+  scale_x_discrete(labels=c("Intercept",
                                                                                            # "t",
                                                                                            #  expression(t^2),
                                                                                            #  expression(t^3),
                                                                                            'Post-Cold War',
                                                                                            "Population",
                                                                                            'GDP per capita',
                                                                                            'Military personnel',
                                                                                            "Number of rebel groups",
                                                                                            "Terrain ruggedness",
                                                                                            "Splinter",
                                                                                            "Co-ethnic",
                                                                                            "Co-ideological")) + 
  ggtitle("")
ggsave("Appendix_figures/A17_a.pdf", plot = logit_f_coef, width = 12, height = 7.5)


library(postregplots)
logit_f2_mar1 <- density_ridgeplot(logit_f2_1, n.sim = 1000, varname = c('allCOETH'), 
                                   data = glm_df, val1 = 0, val2 = 1, clusterid = "ccode")
logit_f2_mar1 <-logit_f2_mar1 + xlab("Average of First Difference") + 
  scale_y_discrete(expand = c(-0.1, 0),labels=c("")) + 
  ggtitle("Co-ethnic") + 
  theme(axis.text = element_text(size=14),
        text = element_text(size=14),
        plot.title = element_text(size = 18))


mean(logit_f2_mar1$data$value) #[1] 0.06407264
## marginal for model2 and model 3
logit_f2_mar <- density_ridgeplot(logit_f2_1, n.sim = 1000, varname = c('co_ideology_new1'), 
                                  data = glm_df, val1 = 0, val2 = 1, clusterid = "ccode")

logit_f2_mar <- logit_f2_mar + xlab("Average of First Difference") + 
  ggtitle("Co-ideological") +
  scale_y_discrete(expand = c(-0.1, 0),labels=c("")) + 
  theme(axis.text = element_text(size=14),text = element_text(size=14),
        plot.title = element_text(size = 18))

mean(logit_f2_mar$data$value) #[1] 0.1517364

library(grid)
library(gridExtra)
pdf(file = "Appendix_figures/A17_b.pdf", width = 8.5, height = 4)
grid.arrange(logit_f2_mar,logit_f2_mar1, ncol = 2)
dev.off()



##############################################################################
# A18: Robustness Check: Multinomial Logistic Regression Results for Formal and Informal Alliance Using an Alternative Shared Ideological Constituency Variable (1946-2015)
##############################################################################
data <- data %>% 
  dplyr::mutate(co_ethnic_nuance1 = ifelse(ideolnationlist ==1 & allCOETH ==1, 1, 0),
                co_ethnic_nuance2 = ifelse(allCOETH == 1 & (RBL_A_forge_ideolnat == 0 & RBL_B_forge_ideolnat == 0), 1, 0)) %>% 
  dplyr::mutate(co_ethnic_alternative = ifelse(co_ethnic_nuance1 == 1 | co_ethnic_nuance2 == 1, 1, 0))
nominal_tfit2_1 <- vglm(dv_nominal ~ co_ideology_new1 + allCOETH + splinter_indirect2+
                          rebels_count + pop_log + gdppc_log + milper + rugg_prop + post_cold + 
                          t_yr + t_yrsq + t_yrcub,
                        data= data, multinomial(refLevel="noalliance"))


coef_mlogit_t <- mlogit_scaled_coef(fit = list(nominal_tfit2_1), subvar = TRUE, vars = c(
  "(Intercept)",
  # "t_yr", "t_yrsq","t_yrcub",
  "post_cold", "pop_log","gdppc_log",
  "milper","rebels_count","rugg_prop","splinter_indirect2", 
  "allCOETH","co_ideology_new1"
))

coef_mlogit_t <- coef_mlogit_t + theme(legend.position = "none")+ scale_x_discrete(labels=c("Intercept",
                                                                                            #   "t",
                                                                                            #  expression(t^2),
                                                                                            #  expression(t^3),
                                                                                            'Post-Cold War',
                                                                                            "Population",
                                                                                            'GDP per capita',
                                                                                            'Military personnel',
                                                                                            "Number of rebel groups",
                                                                                            "Terrain ruggedness",
                                                                                            "Splinter",
                                                                                            "Co-ethnic", 
                                                                                            "Co-ideological")) + 
  ggtitle("")
ggsave("Appendix_figures/A18_a.pdf", plot = coef_mlogit_t, width = 12, height = 7.5)


p_ideology <- vglm_sim(mod = nominal_tfit2_1, var = "co_ideology_new1", n_sim = 1000)

p4 <- p_ideology + scale_y_discrete(expand = c(-0.1, 0),
                                    labels=c("Formal alliance","Informal alliance")) + 
  xlab("Average of First Difference") + 
  ggtitle("Co-ideological") + 
  theme(axis.text = element_text(size=14),
        text = element_text(size=14),
        plot.title = element_text(hjust = .5, size = 16))
p_ideology$data %>% group_by(type) %>% summarise(mean = mean(Pred))

## 
p_allCOETH <- vglm_sim(mod = nominal_tfit2_1, var = "allCOETH", n_sim = 1000)

p3 <- p_allCOETH + scale_y_discrete(expand = c(-0.1, 0),
                                    labels=c("Formal alliance","Informal alliance")) + 
  xlab("Average of First Difference") + 
  ggtitle("Co-ethnic") + 
  theme(axis.text = element_text(size=14),
        text = element_text(size=14),
        plot.title = element_text(hjust = .5, size = 16)) 

p_allCOETH$data %>% group_by(type) %>% summarise(mean = mean(Pred))


pdf(file = "Appendix_figures/A18_b.pdf", width = 12.5, height = 5)
grid.arrange(p4, p3, nrow = 1)
dev.off()


###########################################################
# A19: Robustness Check: Regressions Results (types of shared ideology from FORGE
##############################################################################

glm_df <- glm_df %>% 
  dplyr::mutate(co_ideol_forge = ifelse(leftwing==1 | religion ==1, 1, 0))
table(glm_df$co_ideol_forge)

logit_r0 <- glm(alliance ~ co_ideol_forge + allCOETH +splinter_indirect2+
                  rebels_count + pop_log + gdppc_log + milper + rugg_prop + post_cold + 
                  t_yr + t_yrsq + t_yrcub,
                data=glm_df, family='binomial'(link = "logit"))

logit_r1 <- glm(alliance ~ co_ideol_forge + ideolnationlist + allCOETH +splinter_indirect2+
                  rebels_count + pop_log + gdppc_log + milper + rugg_prop + post_cold + 
                  t_yr + t_yrsq + t_yrcub,
                data=glm_df, family='binomial'(link = "logit"))
logit_r2 <- glm(alliance ~ co_ideol_forge + ideolnationlist +splinter_indirect2+
                  rebels_count + pop_log + gdppc_log + milper + rugg_prop + post_cold + 
                  t_yr + t_yrsq + t_yrcub,
                data=glm_df, family='binomial'(link = "logit"))

logit_r22 <- glm(alliance ~ leftwing + religion+ allCOETH +splinter_indirect2+
                   rebels_count + pop_log + gdppc_log + milper + rugg_prop + post_cold + 
                   t_yr + t_yrsq + t_yrcub,
                 data=glm_df, family='binomial'(link = "logit"))

logit_r3 <- glm(alliance ~  leftwing + religion+ allCOETH + ideolnationlist +splinter_indirect2+
                  rebels_count + pop_log + gdppc_log + milper + rugg_prop + post_cold + 
                  t_yr + t_yrsq + t_yrcub,
                data=glm_df, family='binomial'(link = "logit"))

logit_r4 <- glm(alliance ~  leftwing + religion+  ideolnationlist +splinter_indirect2+
                  rebels_count + pop_log + gdppc_log + milper + rugg_prop + post_cold + 
                  t_yr + t_yrsq + t_yrcub,
                data=glm_df, family='binomial'(link = "logit"))


logit_r_coef <- scaled_coef_cluster(ModelResults = list(logit_r0, logit_r1, 
                                                        logit_r2,logit_r22, logit_r3, logit_r4), data = glm_df, 
                                    clusterid = "ccode", subvar = TRUE, vars = c(
                                      "(Intercept)",
                                      # "t_yr", "t_yrsq","t_yrcub",
                                      "post_cold", "pop_log","gdppc_log",
                                      "milper","rebels_count","rugg_prop", "splinter_indirect2", 
                                      "allCOETH",
                                      "religion",
                                      "ideolnationlist",
                                      "leftwing",
                                      "co_ideol_forge"
                                    ))

logit_r_coef <- logit_r_coef +  scale_x_discrete(labels=c("Intercept",
                                                          #    "t",
                                                          #   expression(t^2),
                                                          #    expression(t^3),
                                                          'Post-Cold War',
                                                          "Population",
                                                          'GDP per capita',
                                                          'Military personnel',
                                                          "Number of rebel groups",
                                                          "Terrain ruggedness",
                                                          "Splinter",
                                                          "Co-ethnic",
                                                          "Co-religious",
                                                          "Co-nationalist",
                                                          "Co-leftwing",
                                                          "Co-ideological"
)) + 
  ggtitle("")
ggsave("Appendix_figures/A19_a.pdf", plot = logit_r_coef, width = 12, height = 7.5)

## for multinomial 


data <- data %>% 
  dplyr::mutate(co_ideol_forge = ifelse(leftwing==1 |  religion ==1, 1, 0))


n_0 <- vglm(dv_nominal ~  co_ideol_forge + allCOETH +splinter_indirect2+
              rebels_count + pop_log + gdppc_log + milper + rugg_prop + post_cold + 
              t_yr + t_yrsq + t_yrcub,
            data= data, multinomial(refLevel="noalliance"))

n_1 <- vglm(dv_nominal ~  co_ideol_forge + ideolnationlist + allCOETH +splinter_indirect2+
              rebels_count + pop_log + gdppc_log + milper + rugg_prop + post_cold + 
              t_yr + t_yrsq + t_yrcub,
            data= data, multinomial(refLevel="noalliance"))

n_2 <- vglm(dv_nominal ~ co_ideol_forge + ideolnationlist +splinter_indirect2+
              rebels_count + pop_log + gdppc_log + milper + rugg_prop + post_cold + 
              t_yr + t_yrsq + t_yrcub,
            data= data, multinomial(refLevel="noalliance"))

n_22 <- vglm(dv_nominal ~ leftwing + religion+ allCOETH +splinter_indirect2+
               rebels_count + pop_log + gdppc_log + milper + rugg_prop + post_cold + 
               t_yr + t_yrsq + t_yrcub,
             data= data, multinomial(refLevel="noalliance"))

n_3 <- vglm(dv_nominal ~ leftwing + religion+ allCOETH + ideolnationlist +splinter_indirect2+
              rebels_count + pop_log + gdppc_log + milper + rugg_prop + post_cold + 
              t_yr + t_yrsq + t_yrcub,
            data= data, multinomial(refLevel="noalliance"))

n_4 <- vglm(dv_nominal ~ leftwing + religion+  ideolnationlist +splinter_indirect2+
              rebels_count + pop_log + gdppc_log + milper + rugg_prop + post_cold + 
              t_yr + t_yrsq + t_yrcub,
            data= data, multinomial(refLevel="noalliance"))


coef_mlogit_t3 <- mlogit_scaled_coef(fit = list(n_0, n_1, n_2,  n_22, n_3, n_4 ), subvar = TRUE, vars = c(
  "(Intercept)",
  # "t_yr", "t_yrsq","t_yrcub",
  "post_cold", "pop_log","gdppc_log",
  "milper","rebels_count","rugg_prop","splinter_indirect2", 
  "allCOETH",
  "religion",
  "ideolnationlist",
  "leftwing",
  "co_ideol_forge"
))

coef_mlogit_t3 <- coef_mlogit_t3 +  scale_x_discrete(labels=c("Intercept",
                                                              #   "t",
                                                              #  expression(t^2),
                                                              #  expression(t^3),
                                                              'Post-Cold War',
                                                              "Population",
                                                              'GDP per capita',
                                                              'Military personnel',
                                                              "Number of rebel groups",
                                                              "Terrain ruggedness",
                                                              "Splinter",
                                                              "Co-ethnic",
                                                              "Co-religious",
                                                              "Co-nationalist",
                                                              "Co-leftwing",
                                                              "Co-ideological")) + 
  ggtitle("")
ggsave("Appendix_figures/A19_b.pdf", plot = coef_mlogit_t3, width = 12, height = 7.5)

##
##############################################################################
# A20: Robustness Check: Logistic and Multinomial Logistic Regression Results Using Proximate Ideology and Proxi-Distant Ideology (1946-2015)
##############################################################################
glm_df <- glm_df %>% 
  dplyr::mutate(proximate_ideo1 = ifelse(co_ideology_new1 == 1, 1, 0),
                proximate_distant_ideo1 = ifelse(co_ideology ==1 & co_ideology_new1 == 0, 1, 0),
                proximate_ideo2 = ifelse(co_ideology_new2 == 1, 1, 0),
                proximate_distant_ideo2 = ifelse(co_ideology ==1 & co_ideology_new2 == 0, 1, 0))



logit_r1c5 <- glm(alliance ~ proximate_distant_ideo1 + proximate_ideo1+
                    allCOETH +splinter_indirect2+
                    rebels_count + pop_log + gdppc_log + milper + rugg_prop + post_cold + 
                    t_yr + t_yrsq + t_yrcub,
                  data=glm_df, family='binomial'(link = "logit"))




logit_f_coef <- scaled_coef_cluster(ModelResults = list(logit_r1c5), data = glm_df, 
                                    clusterid = "ccode", subvar = TRUE, vars = c(
                                      "(Intercept)",
                                      # "t_yr", "t_yrsq","t_yrcub",
                                      "post_cold", "pop_log","gdppc_log",
                                      "milper","rebels_count","rugg_prop", "splinter_indirect2", 
                                      "allCOETH",
                                      "proximate_ideo2",
                                      "proximate_ideo1",
                                      "proximate_distant_ideo1"
                                    ))

logit_f_coef <- logit_f_coef + theme(legend.position = "none")+
  scale_x_discrete(labels=c("Intercept",
                            # "t",
                            #  expression(t^2),
                            #  expression(t^3),
                            'Post-Cold War',
                            "Population",
                            'GDP per capita',
                            'Military personnel',
                            "Number of rebel groups",
                            "Terrain ruggedness",
                            "Splinter",
                            "Co-ethnic",
                            #"Proximate ideology-2",
                            "Proximate ideology",
                            "Proxi-distant ideology")) + 
  ggtitle("")
ggsave("Appendix_figures/A20_a.pdf", plot = logit_f_coef, width = 12, height = 7.5)


data <- data %>% 
  dplyr::mutate(proximate_ideo1 = ifelse(co_ideology_new1 == 1, 1, 0),
                proximate_distant_ideo1 = ifelse(co_ideology ==1 & co_ideology_new1 == 0, 1, 0),
                proximate_ideo2 = ifelse(co_ideology_new2 == 1, 1, 0),
                proximate_distant_ideo2 = ifelse(co_ideology ==1 & co_ideology_new2 == 0, 1, 0))


nominal_r1c5 <- vglm(dv_nominal ~ proximate_ideo1 + proximate_distant_ideo1 + 
                       allCOETH +splinter_indirect2+
                       rebels_count + pop_log + gdppc_log + milper + rugg_prop + post_cold + 
                       t_yr + t_yrsq + t_yrcub,
                     data= data, multinomial(refLevel="noalliance"))


coef_mlogit_t <- mlogit_scaled_coef(fit = list(nominal_r1c5), subvar = TRUE, vars = c(
  "(Intercept)",
  # "t_yr", "t_yrsq","t_yrcub",
  "post_cold", "pop_log","gdppc_log",
  "milper","rebels_count","rugg_prop", "splinter_indirect2", 
  "allCOETH",#"proximate_ideo2",
  "proximate_ideo1","proximate_distant_ideo1"))

coef_mlogit_t <- coef_mlogit_t + theme(legend.position = "none") +
  scale_x_discrete(labels=c("Intercept",
                            # "t",
                            #  expression(t^2),
                            #  expression(t^3),
                            'Post-Cold War',
                            "Population",
                            'GDP per capita',
                            'Military personnel',
                            "Number of rebel groups",
                            "Terrain ruggedness",
                            "Splinter",
                            "Co-ethnic",
                            #"Proximate ideology-2",
                            "Proximate ideology",
                            "Proxi-distant ideology")) + 
  ggtitle("")
ggsave("Appendix_figures/A20_b.pdf", plot = coef_mlogit_t, width = 12, height = 7.5)

###################################################################################### 
#Appendix: Figure A21
###################################################################################### 
 
map_df <- glm_df %>% 
  dplyr::select(year, ccode, rbl_aname, rbl_bname, alliance, FormalAlliance, InformalAlliance, 
                co_constituent, co_ideology, allCOETH) %>% 
  #group_by(ccode, rbl_aname, rbl_bname) %>% 
  dplyr::group_by(ccode) %>% 
  dplyr::summarise(alliance = sum(alliance),
                   FormalAlliance = sum(FormalAlliance), 
                   InformalAlliance = sum(InformalAlliance), 
                   co_constituent = sum(co_constituent),
                   co_ideology = sum(co_ideology),
                   allCOETH = sum(allCOETH),
                   dyads = n())

map_df <- map_df %>% 
  dplyr::mutate(alliance = 100*alliance/dyads,
                FormalAlliance = 100*FormalAlliance/dyads, 
                InformalAlliance = 100*InformalAlliance/dyads, 
                co_constituent = 100*co_constituent/dyads,
                co_ideology = 100*co_ideology/dyads,
                allCOETH = 100*allCOETH/dyads)

map_df <- map_df %>%
  dplyr::mutate_at(c("alliance", "FormalAlliance", "InformalAlliance", 
                     "co_constituent", "co_ideology", "allCOETH"), round)
map_df$country_name <- countrycode::countrycode(map_df$ccode, "cown", "country.name")
library(cshapes)
library(scales)
library(ggmap)
world <- cshp(date = as.Date("2014-12-31"), useGW = FALSE)

map = fortify(world, region="COWCODE")
map$id <- as.numeric(map$id)
#join
map <- left_join(map, map_df, by =c("id" = "ccode"))

map$co_constituent_cut <- cut(map$co_constituent,c(0, 10, 30, 50, 80, 100), 
                              include.lowest = TRUE)
map$co_ideology_cut <- cut(map$co_ideology,c(0, 10, 30, 50, 80, 100), 
                           include.lowest = TRUE)

map$allCOETH_cut <- cut(map$allCOETH,c(0, 10, 30, 50, 80, 100), 
                        include.lowest = TRUE)
levels(map$allCOETH_cut)

map$co_constituent_cut <- factor(map$co_constituent_cut, labels = c("0 ~ 10%",   "10% ~ 30%", "30% ~ 50%",
                                                                    "50% ~ 80%", "80% ~ 100%"))
map$co_ideology_cut <- factor(map$co_ideology_cut, labels = c("0 ~ 10%",   "10% ~ 30%", "30% ~ 50%",
                                                              "50% ~ 80%", "80% ~ 100%"))

map$allCOETH_cut <- factor(map$allCOETH_cut, labels = c("0 ~ 10%",   "10% ~ 30%", "30% ~ 50%",
                                                        "50% ~ 80%", "80% ~ 100%"))



pdf("Appendix_figures/A21_a.pdf", width = 7, height = 4.5)
ggplot() + 
  geom_polygon(data = map, aes(x = long, y = lat, group = group,
                               fill = co_ideology_cut), size = 0.25) +
  #scale_color_viridis(na.value="#FFFFFF00") +
  xlab('') + ylab('') + ggtitle("Co-ideological Rebel Dyads (1946-2015)") +
  theme(
    line = element_blank(),
    rect = element_blank(), #defien the margin line
    legend.position = "right",
    legend.title=element_blank(),
    panel.border=element_blank(),
    panel.grid=element_blank(),
    axis.ticks=element_blank(),
    axis.text=element_blank(),
    plot.title = element_text(hjust = .5, size = 14, face = "bold"))
dev.off()

pdf("Appendix_figures/A21_b.pdf", width = 7, height = 4.5)
ggplot() + 
  geom_polygon(data = map, aes(x = long, y = lat, group = group,
                               fill = allCOETH_cut), size = 0.25) +
  #scale_color_viridis(na.value="#FFFFFF00") +
  xlab('') + ylab('') + ggtitle("Co-ethnic Rebel Dyads (1946-2015)") +
  theme(
    line = element_blank(),
    rect = element_blank(), #defien the margin line
    legend.position = "right",
    legend.title=element_blank(),
    panel.border=element_blank(),
    panel.grid=element_blank(),
    axis.ticks=element_blank(),
    axis.text=element_blank(),
    plot.title = element_text(hjust = .5, size = 14, face = "bold"))
dev.off()

## end of Appendix figures A1-21
