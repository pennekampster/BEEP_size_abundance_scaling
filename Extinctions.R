library(tidyverse)
library(here)
library(viridisLite)
library(readr)

# assess coexistence probability using previously published dataset that contains species-level biomass dynamics:
#https://github.com/pennekampster/Code_and_data_OverallEcosystemStability
dd <- read_csv("https://raw.githubusercontent.com/pennekampster/Code_and_data_OverallEcosystemStability/refs/heads/master/data/species_biomass_BEEP_OES.csv")

max_bm_rich1 <- dd %>% 
  filter(richness == 1) %>% 
  group_by(predicted_species, temperature, day) %>%
  summarize(mean_bm = mean(species_biomass, na.rm=T)) %>%
  group_by(predicted_species, temperature) %>%
  summarize(max_bm = max(mean_bm, na.rm=T))

dd2 <- merge(dd, max_bm_rich1) %>% mutate(bm_scaled = species_biomass / max_bm)

dd2 %>% filter(richness == 2) %>% group_by(richness, combination, temperature, predicted_species, day) %>% summarize(bm_scaled = mean(bm_scaled, na.rm=T)) %>%
ggplot(data=., aes(y=bm_scaled, x=day, group=predicted_species, colour=predicted_species)) + geom_line() + facet_grid(combination~temperature, scales="free")

dd %>% filter(richness == 2 & day > 30) %>% group_by(richness, combination, temperature, predicted_species, day) %>% 
  summarize(species_biomass = (mean(species_biomass, na.rm=T))) %>%
ggplot(data=., aes(y=(species_biomass), x=day, group=predicted_species, colour=predicted_species)) + geom_line() + facet_grid(combination~temperature, scales="free") + geom_hline(yintercept=0.03) + scale_y_log10()


# for species to coexist, they need an average biomass larger than 0.01 over the period > 30 days
dd3 <- dd %>% filter(day > 30) %>% 
  group_by(richness, combination, temperature, predicted_species) %>% 
  summarize(species_biomass = (mean(species_biomass, na.rm=T))) %>% 
  mutate(extinct = case_when(
    species_biomass > 0.01 ~  0,
    species_biomass <= 0.01 ~ 1,
  ))

# extinction per species
cc <- scales::seq_gradient_pal("blue", "red", "Lab")(seq(0,1,length.out=6))
ggplot(data=dd3, aes(x=richness, y=extinct, colour = as.factor(temperature))) + geom_point() + scale_colour_manual(values=cc) + 
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  facet_grid(temperature~predicted_species)  +  
  geom_smooth(method = "glm", method.args = list(family = "binomial"), se=F) + guides(colour="none")
ggsave(here("output/Figure_S4.png"), width=10, height=8)
