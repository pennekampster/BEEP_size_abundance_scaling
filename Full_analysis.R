library(modelsummary)
library(sjPlot)
library(sjmisc)
library(ggplot2)
library(dplyr)
library(ggeffects)
library(interactions)
library(ggthemes)
library(tidyverse)
library(lubridate)
library(here)
library(viridis)
library(lme4)
library(lmerTest)
library(broom.mixed)
library(car)
library(DHARMa)
library(merTools)
library(dotwhisker)
library(dplyr)
library('mgcv')
library(mgcViz)
library(patchwork)
library(viridisLite)

#   dataset
dat <- read_csv("Size_abundance_data.csv")

## Part 1: understanding the effect of temperature and richness on size-abundance scaling
dat <- dat %>% filter(!(is.na(mean_volume) | is.na(N)))

# aggregate on the microcosm level to understand effect of richness and temperature
dat2 <- dat %>% 
  filter(!(richness %in% c(0)) & day %in% seq(1,42, by=2))  %>%
  group_by(temperature, richness, combination, replicate, microcosmID) %>% 
  summarize(mean_major = weighted.mean(mean_major, na.rm=T, w=N), 
            # scale individual volume in microcubicmeter to milicubicmeter dividing by 10^-9
            mean_volume = weighted.mean(mean_volume, na.rm=T, w=N)/1000000000,
            # scale up to total number of individuals per mL
            abundance = sum(N)*13.84/0.6) 

# log-transform size and abundance and encode incubator information
dat3 <- dat2 %>%
  mutate(mean_major = log10(mean_major),
         mean_volume = log10(mean_volume),
         abundance = log10(abundance),
         incubator = case_when(
           as.numeric(microcosmID) > 0 & as.numeric(microcosmID) <= 60 ~ "1",
           as.numeric(microcosmID) > 60 & as.numeric(microcosmID) <= 120 ~ "2",
           as.numeric(microcosmID) > 120 & as.numeric(microcosmID) <= 180 ~ "3",
           as.numeric(microcosmID) > 180 & as.numeric(microcosmID) <= 240 ~ "4",
           as.numeric(microcosmID) > 240 & as.numeric(microcosmID) <= 300 ~ "5",
           as.numeric(microcosmID) > 300 & as.numeric(microcosmID) <= 360 ~ "6",
           as.numeric(microcosmID) > 360 & as.numeric(microcosmID) <= 420 ~ "7",
           as.numeric(microcosmID) > 420 & as.numeric(microcosmID) <= 480 ~ "8",
           as.numeric(microcosmID) > 480 & as.numeric(microcosmID) <= 540 ~ "9",
           as.numeric(microcosmID) > 540 & as.numeric(microcosmID) <= 600 ~ "10",
           as.numeric(microcosmID) > 600 & as.numeric(microcosmID) <= 660 ~ "11",
           as.numeric(microcosmID) > 660 & as.numeric(microcosmID) <= 720 ~ "12"
         ))

# center variables
dat4 <- dat3 %>% ungroup() %>% mutate(abundance_c = abundance - mean(abundance, na.rm=T),
                                      temp_c = temperature-mean(temperature),
                                      rich_c = richness - mean(richness),
                                      major_c = mean_major - mean(mean_major, na.rm=T),
                                      volume_c = mean_volume - mean(mean_volume, na.rm=T))

## Figure 2 and associated analyses

# size-abundance scaling (Figure 2a)
mixed.slope.simple <- lmer(abundance ~ volume_c + (1 |combination) + (1|incubator), dat4)
tidy(mixed.slope.simple, conf.int = T, conf.level = 0.95)
summary(mixed.slope.simple)

# size-abundance scaling with richness and their interaction (Figure 2b)
mixed.slope.rich <- lmer(abundance ~ volume_c * richness  + (1 |combination) + (1|incubator), data = dat4)
summary(mixed.slope.rich)


# figure 2
d.new <- expand.grid(richness= unique(dat4$richness),
                     rich_c = unique(dat4$rich_c),
                     temp_c = unique(dat4$temp_c),
                     major_c = seq(from = min(dat4$major_c), to = max(dat4$major_c), length.out=10),
                     volume_c = seq(from = min(dat4$volume_c), to = max(dat4$volume_c), length.out=10)
)

d.new$abundance_pred_simple <- predict(mixed.slope.simple, newdata = d.new, re.form=NA)
d.new$abundance_pred_rich <- predict(mixed.slope.rich, newdata = d.new, re.form=NA)

cc <- scales::seq_gradient_pal("blue", "red", "Lab")(seq(0,1,length.out=6))


p1 <- ggplot() +
  geom_point(data = dat4, aes(x = volume_c,
                              y = abundance, colour=as.factor(temperature), shape = as.factor(richness)), size = 3) +
  geom_line(data = d.new,aes(x = volume_c,
                             y = abundance_pred_simple), color='black',linetype="dashed") +
  guides(shape=guide_legend(title="richness"),
         colour=guide_legend(title="temperature")) +
  labs(x = expression(paste("Log10 Cell Volume (", mu*l, " ,centered)")), y = "Log10 Abundance (Ind. per mL) ") + theme_few() + scale_colour_manual(values=cc) # + xlim(-0.4, 0.4) + ylim(0.5, 3.5) 

p2 <- ggplot() +
  geom_point(data = dat4, aes(x = volume_c,
                              y = abundance, colour=as.factor(temperature), shape = as.factor(richness)), size = 3) +
  geom_line(data = d.new,aes(x = volume_c,
                             y = abundance_pred_rich, group=richness), linetype="dashed") +
  guides(shape=guide_legend(title="richness"),
         colour=guide_legend(title="temperature")) +
  labs(x = expression(paste("Log10 Cell Volume (", mu*l, " ,centered)")), y = "Log10 Abundance (Ind. per mL) ") + theme_few() + scale_colour_manual(values=cc) + 
  facet_grid(.~richness, labeller = labeller(richness = label_both) ) # + xlim(-0.4, 0.4) + ylim(0.5, 3.5)  

p1 / p2 + patchwork::plot_layout(guides = 'collect', ncol=1) + patchwork::plot_annotation(tag_levels = 'A') + plot_layout(heights = c(4, 2))
ggsave(here("output/figure2.png"), width=8, height=8)


# Table 1 and Figure 3

# size-abundance scaling with temperature and richness effects and all possible interactions
mixed.slope.full <- lmer(abundance ~ volume_c * temp_c * rich_c  + (1 |combination) + (1|incubator), data = dat4)
summary(mixed.slope.full)

modelsummary(mixed.slope.full, 
             estimate = "{estimate} [{conf.low}, {conf.high}]", 
             statistic = NULL,
             coef_map = c('(Intercept)' = 'Intercept', 'volume_c' = 'volume', 'temp_c' = 'temperature', 'rich_c' = 'richness', 
                          'volume_c:temp_c' = 'volume:temperature', 'volume_c:rich_c' = 'volume:richness', 'temp_c:rich_c' = 'richness:temperature',
                          'volume_c:temp_c:rich_c' = 'volume:temperature:richness'),
             gof_map = c("none"),
             output = here("output/Part1_model_table.docx"))


# plot temperature x richness
cc <- scales::seq_gradient_pal("blue", "red", "Lab")(seq(0,1,length.out=6))
p <- plot_model(mixed.slope.full, type = "pred", prefix.labels = c("label"), terms=c("volume_c", "temp_c", "rich_c"), ci.lvl = NA, legend.title = "Temperature", title="") + theme_few() + scale_colour_manual(values=cc, labels=seq(15,25, by=2)) + geom_abline(intercept=4.7, slope=-0.75, linetype="dashed") + labs(x = "Log10 Cell Volume (centered)", y = "Log10 Abundance (Ind. per mL) ")
p

p$data$facet <- factor(case_when(
  p$data$facet == "rich_c = -2.04993" ~ "richness = 1", 
  p$data$facet == "rich_c = -1.04993" ~ "richness = 2", 
  p$data$facet == "rich_c = -0.0499266" ~ "richness = 3", 
  p$data$facet == "rich_c = 0.950073" ~ "richness = 4", 
  p$data$facet == "rich_c = 1.95007" ~ "richness = 5", 
  p$data$facet == "rich_c = 2.95007" ~ "richness = 6", 
), levels=c("richness = 1", "richness = 2", "richness = 3", "richness = 4", "richness = 5", "richness = 6"))
p
ggsave(here("output/figure3.png"), width=7, height=5)


# Figure 3

p <- plot_model(mixed.slope.full, type = "pred", terms = c("volume_c", "rich_c", "temp_c"), ci.lvl = NA, legend.title = "Richness", title="") + theme_few() + 
  scale_colour_manual(values=viridis::viridis(6), labels=seq(1,6, by=1)) + geom_abline(intercept=4.7, slope=-0.75, linetype="dashed") + labs(x = "Log10 Cell Volume (centered)", y = "Log10 Abundance (Ind. per mL) ")

p$data$facet <- factor(case_when(
  p$data$facet == "temp_c = -5.00734" ~ "temperature = 15", 
  p$data$facet == "temp_c = -3.00734" ~ "temperature = 17", 
  p$data$facet == "temp_c = -1.00734" ~ "temperature = 19", 
  p$data$facet == "temp_c = 0.992658" ~ "temperature = 21", 
  p$data$facet == "temp_c = 2.99266" ~ "temperature = 23", 
  p$data$facet == "temp_c = 4.99266" ~ "temperature = 25", 
), levels=c("temperature = 15", "temperature = 17", "temperature = 19", "temperature = 21", "temperature = 23", "temperature = 25"))
p
ggsave(here("output/Figure_S1.png"), width=7, height=5)



## Understanding size-abundance scaling through time

# analysis per day
dat2 <- dat %>% 
  dplyr::filter(., richness > 0)  %>%
  group_by(temperature, combination, richness, day, replicate, microcosmID) %>% 
  summarize(mean_major = weighted.mean(mean_major, na.rm=T, w=N), 
            # scale individual volume in microcubicmeter to milicubicmeter dividing by 10^-9
            mean_volume = weighted.mean(mean_volume, na.rm=T, w=N)/1000000000,
            # scale up to total number of individuals per mL
            abundance = sum(N)*13.84/0.6,
            rel_richness=mean(rel_richness)) 

dat2 <- dat2 %>% filter(mean_volume > 0 & abundance > 0)

dat3 <- dat2 %>%
  filter(day %in% seq(1,42, by=2)) %>%
  mutate(mean_major = log10(mean_major),
         mean_volume = log10(mean_volume),
         abundance = log10(abundance),
         incubator = case_when(
           as.numeric(microcosmID) > 0 & as.numeric(microcosmID) <= 60 ~ "1",
           as.numeric(microcosmID) > 60 & as.numeric(microcosmID) <= 120 ~ "2",
           as.numeric(microcosmID) > 120 & as.numeric(microcosmID) <= 180 ~ "3",
           as.numeric(microcosmID) > 180 & as.numeric(microcosmID) <= 240 ~ "4",
           as.numeric(microcosmID) > 240 & as.numeric(microcosmID) <= 300 ~ "5",
           as.numeric(microcosmID) > 300 & as.numeric(microcosmID) <= 360 ~ "6",
           as.numeric(microcosmID) > 360 & as.numeric(microcosmID) <= 420 ~ "7",
           as.numeric(microcosmID) > 420 & as.numeric(microcosmID) <= 480 ~ "8",
           as.numeric(microcosmID) > 480 & as.numeric(microcosmID) <= 540 ~ "9",
           as.numeric(microcosmID) > 540 & as.numeric(microcosmID) <= 600 ~ "10",
           as.numeric(microcosmID) > 600 & as.numeric(microcosmID) <= 660 ~ "11",
           as.numeric(microcosmID) > 660 & as.numeric(microcosmID) <= 720 ~ "12"
         )) 

dat3 <- dat3 %>% mutate(time = day)

dat4 <- dat3 %>% ungroup() %>% mutate(temp_c = temperature - mean(temperature, na.rm=T),
                                      # uncomment if analysis should be run on realized rather than initial richness
                                      rich_c = richness - mean(richness, na.rm=T),
                                      #rich_c = rel_richness - mean(rel_richness, na.rm=T),
                                      major_c = mean_major - median(mean_major, na.rm=T),
                                      volume_c = mean_volume - mean(mean_volume, na.rm=T),
                                      temp_fac = as.factor(temperature),
                                      rich_fac = as.factor(richness),
                                      day_fac = as.factor(day))




dat5 <- dat4 %>% 
  dplyr::select(temperature, combination, richness, day, temp_c, rich_c, temp_fac, rich_fac, major_c, volume_c, abundance, replicate, microcosmID, incubator) %>%  
  group_by(day) %>% 
  nest()

# error on day 39, need to catch error
posslmer = possibly(.f = lmer, otherwise = NULL)

#lmer(abundance ~ volume_c * temp_c * rich_c  + (1 |combination) + (1|incubator), 
#     data=dat5$data[[2]], 
#     na.action=na.omit)

mods <- dat5 %>% mutate(mod_per_day = map(data, ~ posslmer(abundance ~ volume_c * temp_c * rich_c  + (1 |combination) + (1|incubator), 
                                                           data=., 
                                                           na.action=na.omit)))

# remove failing model on day 39 (row 20 due to only working with every second day)
#mods <- mods[-20,]

fixed_per_day <-  plyr::ldply(map(mods$mod_per_day, ~ fixed.effects(.)))

confid_per_day <-  bind_rows(lapply(1:nrow(mods),
                                    function(x) tibble::rownames_to_column(as.data.frame(confint(mods$mod_per_day[[x]])))))


confid_per_day_t <- confid_per_day %>%
  dplyr::mutate(time = rep(mods$day, each=11)) %>%
  rename(lower=`2.5 %`, upper=`97.5 %`) %>%
  dplyr::select(rowname, lower, upper, time) %>%
  gather(var, value, -rowname, -time) %>%
  mutate(rowname2 = paste0(rowname, "_", var),
         rowname3 = paste0(rowname, "_", var, "_", time)) %>%
  dplyr::select(time, rowname2, value) %>%
  spread(rowname2, value)

mods <- cbind(mods, fixed_per_day)

gg_volume <- ggplot()  +
  geom_ribbon(data= confid_per_day_t, aes(ymin = volume_c_lower, ymax = volume_c_upper, x=time), fill="blue", alpha=.25) +
  xlab("")+ geom_abline(intercept=-0.75, slope=0, linetype="dashed") + geom_line(data=mods, aes(y=volume_c, x=day), colour="blue") + geom_point(data=mods, aes(y=volume_c, x=day), colour="blue") +
  theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 30, unit = "pt")) +
  ylab("V slope") + theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

gg_rich <- ggplot() +
  geom_ribbon(data= confid_per_day_t, aes(ymin = `rich_c_lower`, ymax = `rich_c_upper`, x=time), fill="blue", alpha=.25) +
  geom_line(data=mods, aes(y=rich_c, x=day), colour="blue")  + xlab("") + geom_abline(intercept=0, slope=0, linetype="dashed") + geom_point(data=mods, aes(y=rich_c, x=day), colour="blue")+
  ylab("Richness\nmain effect")+ theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

gg_temp <- ggplot() +
  geom_ribbon(data= confid_per_day_t, aes(ymin = `temp_c_lower`, ymax = `temp_c_upper`, x=time), fill="blue", alpha=.25) +
  geom_line(data=mods, aes(y=temp_c, x=day), colour="blue")  + xlab("") + geom_abline(intercept=0, slope=0, linetype="dashed") +geom_point(data=mods, aes(y=temp_c, x=day), colour="blue") +
  ylab("Temperature\nmain effect")+ theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

gg_tr <- ggplot() +
  geom_ribbon(data= confid_per_day_t, aes(ymin = `temp_c:rich_c_lower`, ymax = `temp_c:rich_c_upper`, x=time), fill="blue", alpha=.25) +
  geom_line(data=mods, aes(y=`temp_c:rich_c`, x=day), colour="blue") + geom_abline(intercept=0, slope=0, linetype="dashed") + geom_point(data=mods, aes(y=`temp_c:rich_c`, x=day), colour="blue") +
  ylab("TxR slope") + xlab("") + theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

gg_mr <- ggplot() +
  geom_ribbon(data= confid_per_day_t, aes(ymin = `volume_c:rich_c_lower`, ymax = `volume_c:rich_c_upper`, x=time), fill="blue", alpha=.25) +
  geom_line(data=mods, aes(y=`volume_c:rich_c`, x=day), colour="blue") + geom_abline(intercept=0, slope=0, linetype="dashed") +geom_point(data=mods, aes(y=`volume_c:rich_c`, x=day), colour="blue") +
  ylab("VxR slope") + xlab("") + theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

gg_mt <- ggplot() +
  geom_ribbon(data= confid_per_day_t, aes(ymin = `volume_c:temp_c_lower`, ymax = `volume_c:temp_c_upper`, x=time), fill="blue", alpha=.25) +
  geom_line(data=mods, aes(y=`volume_c:temp_c`, x=day), colour="blue") + geom_abline(intercept=0, slope=0, linetype="dashed") +geom_point(data=mods, aes(y=`volume_c:temp_c`, x=day), colour="blue") +
  ylab("VxT slope") + xlab("") + theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

gg_mtr <- ggplot() +
  geom_ribbon(data= confid_per_day_t, aes(ymin = `volume_c:temp_c:rich_c_lower`, ymax = `volume_c:temp_c:rich_c_upper`, x=time), fill="blue", alpha=.25) +
  geom_line(data=mods, aes(y=`volume_c:temp_c:rich_c`, x=day), colour="blue") + geom_abline(intercept=0, slope=0, linetype="dashed") +geom_point(data=mods, aes(y=`volume_c:temp_c:rich_c`, x=day), colour="blue") +
  ylab("VxTxR slope") + xlab("Day") + theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

cowplot::plot_grid(gg_volume, gg_mr, gg_mt, gg_mtr, align = "hv",  ncol=1)
#ggsave(here("output/Figure_S2.png"), width=12, height=9)



##model across time selected#

#FIGURE 8
#plot richness x temperature across day 10, 17, 24, 31

##model across time selected#
# model where time (day) can change in a nonlinear fashion
mixed.slope3 <- lmer(abundance ~ volume_c * temp_c * rich_c * day_fac  + (1 |combination) + (1|incubator), data = dat4)
summary(mixed.slope3)

cc <- scales::seq_gradient_pal("blue", "red", "Lab")(seq(0,1,length.out=6))

p1 <- plot_model(mixed.slope3, type = "pred", terms = c("volume_c", "temp_c", "rich_c[1, 6]", "day_fac [11]"), 
      show.values = T, ci.lvl = NA, title = "Day 11", legend.title = "Temperature") + 
  geom_abline(intercept=3.7, slope=-0.75, linetype="dashed")+ 
  scale_colour_manual(values=cc, labels=seq(15,25, by=2)) + 
  theme_few() + labs(x = expression(paste("Log10 Cell Volume (", mu*l, " ,centered)")), y = "Log10 Abundance (Ind. per mL) ") 

p2 <- plot_model(mixed.slope3, type = "pred", terms = c("volume_c", "temp_c", "rich_c[1, 6]", "day_fac [17]"), 
           ci.lvl = NA, title = "Day 17", legend.title = "Temperature") + 
  geom_abline(intercept=3.7, slope=-0.75, linetype="dashed")+ 
  scale_colour_manual(values=cc, labels=seq(15,25, by=2)) + 
  theme_few()  + labs(x = expression(paste("Log10 Cell Volume (", mu*l, " ,centered)")), y = "Log10 Abundance (Ind. per mL) ")

p3 <- plot_model(mixed.slope3, type = "pred", terms = c("volume_c", "temp_c", "rich_c[1, 6]", "day_fac [23]"), 
           ci.lvl = NA, title = "Day 23", legend.title = "Temperature") + 
  geom_abline(intercept=3.7, slope=-0.75, linetype="dashed")+ scale_colour_manual(values=cc, labels=seq(15,25, by=2)) + 
  theme_few() + labs(x = expression(paste("Log10 Cell Volume (", mu*l, " ,centered)")), y = "Log10 Abundance (Ind. per mL) ")

p4 <- plot_model(mixed.slope3, type = "pred", terms = c("volume_c", "temp_c", "rich_c[1, 6]", "day_fac [29]"), 
           ci.lvl = NA, title = "Day 29", legend.title = "Temperature") + 
  geom_abline(intercept=3.7, slope=-0.75, linetype="dashed")+ scale_colour_manual(values=cc, labels=seq(15,25, by=2))  + 
  theme_few() + labs(x = expression(paste("Log10 Cell Volume (", mu*l, " ,centered)")), y = "Log10 Abundance (Ind. per mL) ")


p1$data$facet <- factor(case_when(
  p1$data$facet == "rich_c = 1" ~ "richness = 1", 
  p1$data$facet == "rich_c = 6" ~ "richness = 6", 
), levels=c("richness = 1", "richness = 6"))

p2$data$facet <- factor(case_when(
  p2$data$facet == "rich_c = 1" ~ "richness = 1", 
  p2$data$facet == "rich_c = 6" ~ "richness = 6", 
), levels=c("richness = 1", "richness = 6"))

p3$data$facet <- factor(case_when(
  p3$data$facet == "rich_c = 1" ~ "richness = 1", 
  p3$data$facet == "rich_c = 6" ~ "richness = 6", 
), levels=c("richness = 1", "richness = 6"))

p4$data$facet <- factor(case_when(
  p4$data$facet == "rich_c = 1" ~ "richness = 1", 
  p4$data$facet == "rich_c = 6" ~ "richness = 6", 
), levels=c("richness = 1", "richness = 6"))


p <- (p1 + p3)  / (p2 + p4) + patchwork::plot_layout(guides = 'collect', axis_titles = 'collect') + patchwork::plot_annotation(tag_levels = 'A') & theme(legend.position = "bottom")
p
ggsave(here("output/Figure_S3.png"), width=10, height=8)




#FIGURE 9
##Proceeding data of simulated communities
combinations <- unique(dat$combination)  
combinations2 <- combinations[grepl(",", combinations)]

#Proceeding data#
mono <- dat %>% 
  filter(richness == 1 & predicted_species != "none") %>% group_by(combination, day, temperature) %>% 
  summarize(mean_major = mean(mean_major), mean_volume = mean(mean_volume),  N = mean(N))

simulated_poly <- vector("list", length(combinations2))



for (i in 1:length(combinations2)){
  print(combinations2[i])
  
  species_list <- unlist(str_split(combinations2[i], ","))
  species_list <- gsub(" ", "", species_list)
  
  simulated_poly[[i]] <- mono %>% filter(combination %in% species_list) %>% group_by(day, temperature) %>% 
    summarize(mean_major1 = mean(mean_major),
              mean_major_cwm = weighted.mean(mean_major, w = N), 
              mean_volume_cwm = weighted.mean(mean_volume, na.rm=T, w=N)/1000000000, 
              abundance = sum(N)*13.84/0.6, 
              abundance_corrected = sum(N)/length(species_list)) %>% 
    mutate(log_mean_major1 = log10(mean_major1),
           log_mean_major_cwm = log10(mean_major_cwm),
           log_mean_volume_cwm = log10(mean_volume_cwm),
           log_abundance = log10(abundance),
           log_abundance_corrected = log10(abundance_corrected),
           combination = combinations2[i],
           richness = length(species_list)) 
}

simulated_poly_dd <- bind_rows(simulated_poly)

dat2 <- dat %>% 
  #filter(., day > 9) %>%
  #filter(., day > 21) %>%
  filter(!(richness %in% c(0)))  %>%
  filter(day %in% seq(1,42, by=2)) %>%
  group_by(temperature, combination, richness, day) %>% 
  summarize(mean_major = weighted.mean(mean_major, na.rm=T, w=N), 
            mean_volume = weighted.mean(mean_volume, na.rm=T, w=N)/1000000000,
            abundance = sum(N)*13.84/0.6) 

original_communities <- dat2 %>%
  filter(richness > 1) %>%
  filter(day %in% seq(1,42, by=2)) %>%
  mutate(log_mean_major = log10(mean_major),
         log_mean_volume = log10(mean_volume),
         log_abundance = log10(abundance))

# plot size-abundance scaling across simulated and real communities (not shown in ms)
ggplot() + geom_point(data=original_communities, aes(y=log_abundance, x=log_mean_volume), alpha=.1) +
  geom_point(data=simulated_poly_dd, aes(y=log_abundance, x=log_mean_volume_cwm), alpha=.1, colour="red") + 
  stat_smooth(data=original_communities, aes(y=log_abundance, x=log_mean_volume), method="lm", colour="black")+ 
  stat_smooth(data=simulated_poly_dd, aes(y=log_abundance, x=log_mean_volume_cwm), method="lm", colour="red")


dat4 <- simulated_poly_dd %>% ungroup() %>% 
  filter(day %in% seq(1,42, by=2)) %>% mutate(temp_c = temperature - mean(temperature, na.rm=T),
                                              rich_c = richness - mean(richness, na.rm=T),
                                              major_c = log_mean_major_cwm - mean(log_mean_major_cwm, na.rm=T),
                                              volume_c = log_mean_volume_cwm - mean(log_mean_volume_cwm, na.rm=T),
                                              temp_fac = as.factor(temperature),
                                              rich_fac = as.factor(richness))

dat5 <- dat4 %>% 
  dplyr::select(temperature, combination, richness, day, temp_c, rich_c, temp_fac, rich_fac, major_c, volume_c, log_abundance) %>%  
  group_by(day) %>% 
  nest()

# error on day 39, need to catch error
posslmer = possibly(.f = lmer, otherwise = NULL)

mods_sim <- dat5 %>% mutate(mod_per_day = map(data, ~ posslmer(log_abundance ~ volume_c * temp_c * rich_c  + (1 |combination), 
                                                               data=., 
                                                               na.action=na.omit)))

# remove failing model on day 39
#mods_sim <- mods_sim[-20,]

fixed_per_day_sim <-  plyr::ldply(map(mods_sim$mod_per_day, ~ fixed.effects(.)))

confid_per_day_sim <-  bind_rows(lapply(1:nrow(mods_sim),
                                        function(x) tibble::rownames_to_column(as.data.frame(confint(mods_sim$mod_per_day[[x]])))))


confid_per_day_t_sim <- confid_per_day_sim %>%
  dplyr::mutate(time = rep(mods_sim$day, each=10)) %>%
  rename(lower=`2.5 %`, upper=`97.5 %`) %>%
  dplyr::select(rowname, lower, upper, time) %>%
  gather(var, value, -rowname, -time) %>%
  mutate(rowname2 = paste0(rowname, "_", var),
         rowname3 = paste0(rowname, "_", var, "_", time)) %>%
  dplyr::select(time, rowname2, value) %>%
  spread(rowname2, value)

mods_sim <- cbind(mods_sim, fixed_per_day_sim)


# joint figure comparing observed and simulated community size abundance scaling



cols <- c("simulated" = "red", "real" = "blue")

gg_volume <- ggplot()  +
  geom_ribbon(data= confid_per_day_t, aes(ymin = volume_c_lower, ymax = volume_c_upper, x=time, fill="real"), alpha=.25) +
  xlab("")+ geom_abline(intercept=-0.75, slope=0, linetype="dashed") + geom_line(data=mods, aes(y=volume_c, x=day, colour="real")) + 
  geom_ribbon(data= confid_per_day_t_sim, aes(ymin = volume_c_lower, ymax = volume_c_upper, x=time, fill="simulated"), alpha=.25) +
  geom_line(data=mods_sim, aes(y=volume_c, x=day, colour="simulated")) + 
  theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 30, unit = "pt")) +
  ylab("V slope") + scale_color_manual(values = cols, aesthetics = c("color", "fill")) + 
  guides(fill=guide_legend("Community"), colour=guide_legend("Community")) + theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

gg_rich <- ggplot() +
  geom_ribbon(data= confid_per_day_t, aes(ymin = `rich_c_lower`, ymax = `rich_c_upper`, x=time, fill="real"), alpha=.25) +
  geom_line(data=mods, aes(y=rich_c, x=day, colour="real"))  + 
  geom_ribbon(data= confid_per_day_t_sim, aes(ymin = `rich_c_lower`, ymax = `rich_c_upper`, x=time, fill="simulated"), alpha=.25) +
  geom_line(data=mods_sim, aes(y=rich_c, x=day, colour="simulated"))  + 
  xlab("") + geom_abline(intercept=0, slope=0, linetype="dashed") + 
  ylab("Richness\nmain effect")+ scale_color_manual(values = cols, aesthetics = c("color", "fill")) + guides(fill=guide_legend("Community"), colour=guide_legend("Community"))+ theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

gg_temp <- ggplot() +
  geom_ribbon(data= confid_per_day_t, aes(ymin = `temp_c_lower`, ymax = `temp_c_upper`, x=time, fill="real"), alpha=.25) +
  geom_line(data=mods, aes(y=temp_c, x=day, colour="real"))  + xlab("") +
  geom_ribbon(data= confid_per_day_t_sim, aes(ymin = `temp_c_lower`, ymax = `temp_c_upper`, x=time, fill="simulated"), alpha=.25) +
  geom_line(data=mods_sim, aes(y=temp_c, x=day, colour="simulated"))  + xlab("") +
  geom_abline(intercept=0, slope=0, linetype="dashed") +
  ylab("Temperature\nmain effect") + scale_color_manual(values = cols, aesthetics = c("color", "fill")) + guides(fill=guide_legend("Community"), colour=guide_legend("Community"))+ theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

gg_tr <- ggplot() +
  geom_ribbon(data= confid_per_day_t, aes(ymin = `temp_c:rich_c_lower`, ymax = `temp_c:rich_c_upper`, x=time, fill="real"), alpha=.25) +
  geom_line(data=mods, aes(y=`temp_c:rich_c`, x=day, colour="real")) + 
  geom_ribbon(data= confid_per_day_t_sim, aes(ymin = `temp_c:rich_c_lower`, ymax = `temp_c:rich_c_upper`, x=time, fill="simulated"), alpha=.25) +
  geom_line(data=mods_sim, aes(y=`temp_c:rich_c`, x=day, colour="simulated")) + 
  geom_abline(intercept=0, slope=0, linetype="dashed") +
  ylab("TxR slope") + xlab("Day") + scale_color_manual(values = cols, aesthetics = c("color", "fill")) + guides(fill=guide_legend("Community"), colour=guide_legend("Community"))+ theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

gg_mr <- ggplot() +
  geom_ribbon(data= confid_per_day_t, aes(ymin = `volume_c:rich_c_lower`, ymax = `volume_c:rich_c_upper`, x=time, fill="real"), alpha=.25) +
  geom_line(data=mods, aes(y=`volume_c:rich_c`, x=day, colour="real")) + 
  geom_ribbon(data= confid_per_day_t_sim, aes(ymin = `volume_c:rich_c_lower`, ymax = `volume_c:rich_c_upper`, x=time, fill="simulated"), alpha=.25) +
  geom_line(data=mods_sim, aes(y=`volume_c:rich_c`, x=day, colour="simulated")) + 
  geom_abline(intercept=0, slope=0, linetype="dashed") +
  ylab("VxR slope") + xlab("") + scale_color_manual(values = cols, aesthetics = c("color", "fill")) + guides(fill=guide_legend("Community"), colour=guide_legend("Community"))+ theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

gg_mt <- ggplot() +
  geom_ribbon(data= confid_per_day_t, aes(ymin = `volume_c:temp_c_lower`, ymax = `volume_c:temp_c_upper`, x=time, fill="real"), alpha=.25) +
  geom_line(data=mods, aes(y=`volume_c:temp_c`, x=day, colour="real")) + 
  geom_ribbon(data= confid_per_day_t_sim, aes(ymin = `volume_c:temp_c_lower`, ymax = `volume_c:temp_c_upper`, x=time, fill="simulated"), alpha=.25) +
  geom_line(data=mods_sim, aes(y=`volume_c:temp_c`, x=day, colour="simulated")) + 
  geom_abline(intercept=0, slope=0, linetype="dashed") +
  ylab("VxT slope") + xlab("") + scale_color_manual(values = cols, aesthetics = c("color", "fill")) + guides(fill=guide_legend("Community"), colour=guide_legend("Community"))+ theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

gg_mtr <- ggplot() + 
  geom_ribbon(data= confid_per_day_t, aes(ymin = `volume_c:temp_c:rich_c_lower`, ymax = `volume_c:temp_c:rich_c_upper`, x=time, fill="real"), alpha=.25) +
  geom_line(data=mods, aes(y=`volume_c:temp_c:rich_c`, x=day, colour="real")) + 
  geom_ribbon(data= confid_per_day_t_sim, aes(ymin = `volume_c:temp_c:rich_c_lower`, ymax = `volume_c:temp_c:rich_c_upper`, x=time, fill="simulated"), alpha=.25) +
  geom_line(data=mods_sim, aes(y=`volume_c:temp_c:rich_c`, x=day, colour="simulated")) + 
  geom_abline(intercept=0, slope=0, linetype="dashed") +
  ylab("VxTxR slope") + xlab("Day") + scale_color_manual(values = cols, aesthetics = c("color", "fill")) + guides(fill=guide_legend("Community"), colour=guide_legend("Community"))+ theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

gg_volume + gg_mr + gg_mt + gg_mtr + patchwork::plot_layout(guides = 'collect', ncol=1) + patchwork::plot_annotation(tag_levels = 'A', caption = 'V = volume, R = richness, T = temperature') & theme(legend.position = "bottom", plot.caption = element_text(hjust = 0))
ggsave(here("output/Figure4.png"), width=8, height=10)
