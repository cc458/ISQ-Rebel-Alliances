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

#### set working directory to the folder named "Replication_ISQ"
#setwd("Replication_ISQ")
#####################

source("Code/helpFUN.R")


##############################################################################################################################
##############################################################################################################################


##############################################################################
# Figure 1: temporal trends
##############################################################################

## make a figure for temporal pattern
load("Data/glm_df4.RData")
time_df <- glm_df4 %>% 
        dplyr::select(year, alliance, FormalAlliance, InformalAlliance, 
                      co_constituent, co_ideology, allCOETH) %>% 
            dplyr::group_by(year) %>% 
            dplyr::summarise(alliance = sum(alliance),
                      FormalAlliance = sum(FormalAlliance), 
                      InformalAlliance = sum(InformalAlliance), 
                      co_constituent = sum(co_constituent),
                      co_ideology = sum(co_ideology),
                      allCOETH = sum(allCOETH),
                      dyads = n())

time_df <- time_df %>% 
  dplyr::mutate(alliance = 100*alliance/dyads,
                FormalAlliance = 100*FormalAlliance/dyads, 
                InformalAlliance = 100*InformalAlliance/dyads, 
                co_constituent = 100*co_constituent/dyads,
                co_ideology = 100*co_ideology/dyads,
                allCOETH = 100*allCOETH/dyads)

time_df <- time_df %>%
        dplyr::mutate_at(c("alliance", "FormalAlliance", "InformalAlliance", 
                           "co_constituent", "co_ideology", "allCOETH"), round)

library(ggrepel)
library(directlabels)
library(tidyr)
time_ally <- time_df %>% 
  dplyr::select(year, alliance, FormalAlliance, InformalAlliance) %>% 
  ungroup() %>% 
  gather(ally_type, number, alliance: InformalAlliance) %>% 
  dplyr::mutate(ally_type = factor(ally_type))
levels(time_ally$ally_type) <- c("Alliance", "Formal alliance", "Informal alliance")


##############################################################################
###################  Figure 1-a  ###########################
##############################################################################

time_ally <- ggplot(time_ally, aes(x = year, y =  number)) + 
  geom_line() + facet_wrap(~ ally_type) + 
  scale_colour_discrete(guide = 'none') +
  theme_bw() +
  scale_x_continuous(breaks = seq(1945, 2015, 15)) + 
  ylab("% of dyad-years") + xlab("Year") + 
  ggtitle("") + 
  theme(legend.position="bottom",
        legend.title=element_blank(),
        axis.text = element_text(size=12),
        text = element_text(size=12),
        strip.text.x = element_text(size = 12, color='white',
                                    angle=0),
        strip.background = element_rect(fill = "#525252", color='#525252'),
        plot.title = element_text(hjust = .5, size = 13, face = "bold"))
ggsave("Figure/figure1_a.eps", plot =time_ally, width = 6.5, height = 3.5, dpi = 400, units = "in")



##############################################################################
###################  Figure 1-b  ###########################
##############################################################################

time_constituency <- time_df %>% 
  dplyr::select(year,co_constituent, co_ideology, allCOETH) %>% 
  ungroup() %>% 
  gather(constituency_type, number, co_constituent: allCOETH) %>% 
  dplyr::mutate(constituency = factor(constituency_type, levels = c("co_constituent","co_ideology","allCOETH" )))

levels(time_constituency$constituency) <- c("Co-constituent", "Co-ideological", "Co-ethnic")

time_constituency <- ggplot(time_constituency, aes(x = year, y =  number)) + 
  geom_line() + facet_wrap(~ constituency) + 
  scale_colour_discrete(guide = 'none') +
  theme_bw() +
  scale_x_continuous(breaks = seq(1945, 2015, 15)) + 
  ylab("% of dyad-years") + xlab("Year") + 
  ggtitle("") + 
  theme(legend.position="bottom",
        legend.title=element_blank(),
        axis.text = element_text(size=12),
        text = element_text(size=12),
        strip.text.x = element_text(size = 12, color='white',
                                    angle=0),
        strip.background = element_rect(fill = "#525252", color='#525252'),
        plot.title = element_text(hjust = .5, size = 13, face = "bold"))
ggsave("Figure/figure1_b.eps", plot =time_constituency, width = 6.5, height = 3.5, dpi = 400, units = "in")


##############################################################################
###################  Figure 2  ###########################
##############################################################################


## as % of world dyads in each country
# Map
map_df2 <- glm_df4 %>% 
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
                   dyads = n()) %>% 
  ungroup()

map_df2 <- map_df2 %>% 
  dplyr::mutate(alliance_sum = sum(alliance),
                FormalAlliance_sum = sum(FormalAlliance), 
                InformalAlliance_sum = sum(InformalAlliance), 
                co_constituent_sum = sum(co_constituent),
                co_ideology_sum = sum(co_ideology),
                allCOETH_sum = sum(allCOETH),
                dyads_sum = sum(dyads))


map_df2 <- map_df2 %>% 
  dplyr::mutate(alliance = 100*alliance/alliance_sum,
                FormalAlliance = 100*FormalAlliance/FormalAlliance_sum, 
                InformalAlliance = 100*InformalAlliance/InformalAlliance_sum, 
                co_constituent = 100*co_constituent/co_constituent_sum,
                co_ideology = 100*co_ideology/co_ideology_sum,
                allCOETH = 100*allCOETH/allCOETH_sum)

map_df2 <- map_df2 %>%
  dplyr::mutate_at(c("alliance", "FormalAlliance", "InformalAlliance", 
                     "co_constituent", "co_ideology", "allCOETH"), round)
map_df2$country_name <- countrycode::countrycode(map_df2$ccode, "cown", "country.name")
library(cshapes)
library(scales)
library(ggmap)
world <- cshp(date = as.Date("2014-12-31"), useGW = FALSE)

#Taiwan add to China
map = fortify(world, region="GWCODE")
map$id <- as.numeric(map$id)
#join
map <- left_join(map, map_df2, by =c("id" = "ccode"))
#0-23
map$co_constituent_cut <- cut(map$co_constituent,c(0, 0.99, 5, 10, 20, 100), 
                              include.lowest = TRUE)
#0-30
map$co_ideology_cut <- cut(map$co_ideology,c(0, 0.99, 5, 10, 20, 100), 
                           include.lowest = TRUE)
#0-17
map$allCOETH_cut <- cut(map$allCOETH,c(0, 0.99, 5, 10, 15, 100), 
                        include.lowest = TRUE)
levels(map$allCOETH_cut)

map$co_constituent_cut <- factor(map$co_constituent_cut, labels = c("0 ~ 1%",  "1% ~ 5%", "5% ~ 10%",
                                                                    "10% ~ 20%", ">20%"))
map$co_ideology_cut <- factor(map$co_ideology_cut, labels = c("0 ~ 1%",  "1% ~ 5%", "5% ~ 10%",
                                                              "10% ~ 20%", ">20%"))

map$allCOETH_cut <- factor(map$allCOETH_cut, labels = c("0 ~ 1%", "1% ~ 5%", "5% ~ 10%",
                                                        "10% ~ 15%", ">15%"))

##############################################################################
###################  Figure 2-a  ###########################
##############################################################################


f2_a <- ggplot() + 
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
ggsave("Figure/figure2_a.eps", plot = f2_a, width = 6.5, height = 3.5, dpi = 400, units = "in")

##############################################################################
###################  Figure 2-b  ###########################
##############################################################################

fig2_b<- ggplot() + 
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

ggsave("Figure/figure2_b.eps", plot = fig2_b, width = 6.5, height = 3.5, dpi = 400, units = "in")



