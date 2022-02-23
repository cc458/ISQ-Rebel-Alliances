rm(list=ls())
# Function to load packages
loadPkg=function(toLoad){
  for(lib in toLoad){
    if(! lib %in% installed.packages()[,1])
    {install.packages(lib, repos='http://cran.rstudio.com/')}
    suppressMessages( library(lib, character.only=TRUE))}}

# Load libraries
packs=c("readr", "readxl", "dplyr", "countrycode", "ggplot2", "amen","statnet",
         "amen", "purrr", "texreg", "abind", "ergm","reshape2",
        'parallel','foreach','doParallel',"haven", "readstata13")
loadPkg(packs)

########################################################################
# for the AME package, we use the tailored version of the R amen package
# due to changes in compositions of civil war networks over time
# see Gallop, Max, Shahryar Minhas and Cassy Dorff. 2020.
# “Networks of Violence: Predicting Conflict in Nigeria.” The Journal of Politics 82(2):476–493.
# Parts of the R code are also from Gallop et al 2020
############################################################################


######################################################
### 1947-2015
## Figure 5 is based on the results
######################################################

################# AME Modeling

yrs <- as.character(seq(1947, 2015, by =1))
#################
load("Data/ame/AllyList.RData")
load("Data/ame/Formal_AllyList.RData")
load("Data/ame/xDyadL.RData")
load("Data/ame/xNodeL.RData")
AllyList <- AllyList[yrs]
Formal_AllyList <- Formal_AllyList[yrs]


xNodeL <- xNodeL[yrs]
xDyadL <- xDyadL[yrs]
dimnames(xDyadL[[1]])[3]
dimnames(xNodeL[[1]])[2]
# set up model specs
subListArray = function(lA, vars, dims=2){
  if(dims==2){ return( lapply(lA, function(x){ x[,vars, drop=FALSE] }) ) }
  if(dims==3){ return( lapply(lA, function(x){ x[,,vars,drop=FALSE] }) ) } }

designArrays = list(
  ally_ideo =list(
    dyadCovar=subListArray(xDyadL, c("co_ideology", "allCOETH",
                                     "rugg_prop","rebels_count", "splinter_indirect2",
                                     "milper", "gdppc_log", "pop_log", "post_cold"), 3)
  ),
  ally_coco =list(
    dyadCovar=subListArray(xDyadL, c("co_constituent" ,
                                     "rugg_prop","rebels_count", "splinter_indirect2",
                                     "milper", "gdppc_log", "pop_log", "post_cold"), 3)
  ),
  formal_ideo =list(
    dyadCovar=subListArray(xDyadL, c("co_ideology", "allCOETH",
                                     "rugg_prop","rebels_count", "splinter_indirect2",
                                     "milper", "gdppc_log", "pop_log", "post_cold"), 3)
  ),
  formal_coco =list(
    dyadCovar=subListArray(xDyadL, c("co_constituent" ,
                                     "rugg_prop","rebels_count", "splinter_indirect2",
                                     "milper", "gdppc_log", "pop_log", "post_cold"), 3)
  )
)
rm(xDyadL, xNodeL)

# loop 
# parallelize model run
library(amen)
library(parallel)
library(foreach)
## do both at the same time
Y <- list(AllyList, AllyList, Formal_AllyList,Formal_AllyList)
cores=length(Y) ; cl=makeCluster(cores) ; registerDoParallel(cl)

Fit_AME47_15 = foreach(i = 1:length(Y), .packages=c('amen')) %dopar% {
  fit = ame_repL(
    Y=Y[[i]], Xdyad=designArrays[[i]]$dyadCovar,
    symmetric=TRUE, nvar=TRUE, R=2,  
    model='bin', intercept=TRUE, seed=6886,
    burn= 1000, nscan=10000, odens=25,
    plot=FALSE, gof=TRUE, periodicSave=FALSE )
  return(fit)
} ; stopCluster(cl)

names(Fit_AME47_15) = c("ally_ideo", "ally_coco", "formal_ideo", "formal_coco")
save(Fit_AME47_15, file = "Data/ame/Fit_AME47_15.RData")

######################################################
### AME with state-support only 1979-2009
## Figure A11
######################################################

yrs <- as.character(seq(1979, 2009, by =1))
#################
load("Data/ame/AllyList.RData")
load("Data/ame/Formal_AllyList.RData")
load("Data/ame/xDyadL.RData")
load("Data/ame/xNodeL.RData")
AllyList <- AllyList[yrs]
Formal_AllyList <- Formal_AllyList[yrs]


xNodeL <- xNodeL[yrs]
xDyadL <- xDyadL[yrs]
dimnames(xDyadL[[1]])[3]
dimnames(xNodeL[[1]])[2]
# set up model specs
subListArray = function(lA, vars, dims=2){
  if(dims==2){ return( lapply(lA, function(x){ x[,vars, drop=FALSE] }) ) }
  if(dims==3){ return( lapply(lA, function(x){ x[,,vars,drop=FALSE] }) ) } }

designArrays = list(
  MarxistIdeol=list(
    dyadCovar=subListArray(xDyadL, c("MarxistIdeol", "allCOETH", "islamist" ,
                                     "rugg_prop","rebels_count", "splinter_indirect2", "state_cosponsor_dummy",
                                     "milper", "gdppc_log", "pop_log", "post_cold"), 3)
  ),
  RBL_islamist =list(
    dyadCovar=subListArray(xDyadL, c("co_ideology", "allCOETH",
                                     "rugg_prop","rebels_count", "splinter_indirect2","state_cosponsor_dummy",
                                     "milper", "gdppc_log", "pop_log", "post_cold"), 3)
  ),
  RBL_ethn =list(
    dyadCovar=subListArray(xDyadL, c("co_constituent" ,
                                     "rugg_prop","rebels_count", "splinter_indirect2","state_cosponsor_dummy",
                                     "milper", "gdppc_log", "pop_log", "post_cold"), 3)
  )
)
rm(xDyadL, xNodeL)



loadPkg(c('parallel','foreach'))
## do both at the same time
Y <- list(AllyList, Formal_AllyList)
cores=length(Y) ; cl=makeCluster(cores) ; registerDoParallel(cl)

Fit_AME7909 = foreach(i = 1:length(Y), .packages=c('amen')) %dopar% {
  fit = ame_repL(
    Y=Y[[i]], Xdyad=designArrays[[2]]$dyadCovar,
    symmetric=TRUE, nvar=TRUE, R=2,  
    model='bin', intercept=TRUE, seed=6886,
    burn= 1000, nscan=10000, odens=25,
    plot=FALSE, gof=TRUE, periodicSave=FALSE )
  return(fit)
} ; stopCluster(cl)

names(Fit_AME7909) = c("AllyList", "Formal_AllyList")
save(Fit_AME7909, file = "Data/ame/Fit_AME7909.RData")