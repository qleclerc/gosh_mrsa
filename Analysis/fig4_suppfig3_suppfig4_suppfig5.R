
library(lubridate)
library(dplyr)
library(tidyr)
library(reshape2)
library(ggplot2)
library(scales)
library(cowplot)

mrsa_col = "#EFC000FF"
mssa_col = "#0073C2FF"

staph_isolates = read.csv(here::here("Clean", "staph_isolates.csv")) %>%
  mutate(date = as_date(date))

# #mrsa
mrsa_resistances = staph_isolates %>%
  filter(SpeciesName == "Methicillin-Resistant Staphylococcus aureus") %>%
  select(5:56) %>%
  melt(.,id.vars = "date") %>%
  filter(!is.na(value)) %>%
  filter(value %in% c("S", "R")) %>%
  mutate(date = floor_date(date, "month")) %>%
  arrange(date) %>%
  group_by(date, variable, value) %>%
  summarise(n = n()) %>%
  mutate(prop = n/sum(n)) %>%
  mutate(variable = as.character(variable)) %>%
  mutate(SpeciesName = "Methicillin-Resistant Staphylococcus aureus")
# 
# 
# resistances = unique(mrsa_resistances$variable)
# 
# for(i in c(1,9,17,25,33,41)){
#   
#   pp = mrsa_resistances %>%
#     filter(value == "R") %>%
#     filter(variable %in% resistances[(i):(i+7)]) %>%
#     ggplot() +
#     geom_line(aes(date, prop, colour = variable)) +
#     theme_bw() +
#     scale_x_date(limits = as.Date(c("2000-02-01", "2021-11-01")))
#   
#   plot(pp)
#   
# }
# 
# mrsa_resistances %>%
#   filter(value == "R") %>%
#   filter(variable == "Gent.Cipro") %>%
#   ggplot() +
#   geom_line(aes(date, prop, colour = variable)) +
#   theme_bw() +
#   scale_x_date(limits = as.Date(c("2000-02-01", "2021-11-01")))
# 
# 
# #mssa
mssa_resistances = staph_isolates %>%
  filter(SpeciesName == "Methicillin-Susceptible Staphylococcus aureus") %>%
  select(5:56) %>%
  melt(.,id.vars = "date") %>%
  filter(!is.na(value)) %>%
  filter(value %in% c("S", "R")) %>%
  mutate(date = floor_date(date, "month")) %>%
  arrange(date) %>%
  group_by(date, variable, value) %>%
  summarise(n = n()) %>%
  mutate(prop = n/sum(n)) %>%
  mutate(variable = as.character(variable)) %>%
  mutate(SpeciesName = "Methicillin-Susceptible Staphylococcus aureus")
# 
# 
# resistances = unique(mssa_resistances$variable)
# 
# for(i in c(1,9,17,25,33,41)){
#   
#   pp = mssa_resistances %>%
#     filter(value == "R") %>%
#     filter(variable %in% resistances[(i):(i+7)]) %>%
#     ggplot() +
#     geom_line(aes(date, prop, colour = variable)) +
#     theme_bw() +
#     scale_x_date(limits = as.Date(c("2000-02-01", "2021-11-01")))
#   
#   plot(pp)
#   
# }
# 
# mssa_resistances %>%
#   filter(value == "R") %>%
#   filter(variable == "Gent.Cipro") %>%
#   ggplot() +
#   geom_line(aes(date, prop, colour = variable)) +
#   theme_bw() +
#   scale_x_date(limits = as.Date(c("2000-02-01", "2021-11-01")))

staph_resistances = rbind(mrsa_resistances, mssa_resistances) %>%
  tibble() %>%
  complete(date, variable, value, SpeciesName)

pa = staph_resistances %>%
  filter(value == "R") %>%
  filter(variable %in% c("Amikacin", "Amik.Fluclox", "Gent.Cipro", "Gentamicin", "Rifampicin")) %>%
  ggplot() +
  geom_line(aes(date, prop, colour = SpeciesName)) +
  facet_wrap(~variable, ncol = 1) +
  theme_bw() +
  scale_y_continuous(limits = c(0,1)) +
  labs(x = "Time (months)", y = "Proportion of isolates resistant to antibiotic", colour = "") +
  scale_x_date(limits = as.Date(c("2000-02-01", "2021-11-01"))) +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        strip.text = element_text(size = 12),
        legend.text = element_text(size = 12)) +
  scale_colour_manual(values = c(mrsa_col, mssa_col))

pb = staph_resistances %>%
  filter(value == "R") %>%
  filter(variable %in% c("Ciprofloxacin", "Mupirocin", "Clindamycin")) %>%
  ggplot() +
  geom_line(aes(date, prop, colour = SpeciesName)) +
  facet_wrap(~variable, nrow = 1) +
  theme_bw() +
  scale_y_continuous(limits = c(0,1)) +
  labs(x = "Time (months)", y = "Proportion of isolates resistant to antibiotic") +
  scale_x_date(limits = as.Date(c("2000-02-01", "2021-11-01"))) +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        strip.text = element_text(size = 12),
        legend.text = element_text(size = 12)) +
  scale_colour_manual(values = c(mrsa_col, mssa_col))

pc = staph_resistances %>%
  filter(value == "R") %>%
  filter(variable %in% c("Trimethoprim", "Erythromycin")) %>%
  ggplot() +
  geom_line(aes(date, prop, colour = SpeciesName)) +
  facet_wrap(~variable, nrow = 1) +
  theme_bw() +
  scale_y_continuous(limits = c(0,1)) +
  labs(x = "Time (months)", y = "Proportion of isolates resistant to antibiotic") +
  scale_x_date(limits = as.Date(c("2000-02-01", "2021-11-01"))) +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        strip.text = element_text(size = 12),
        legend.text = element_text(size = 12)) +
  scale_colour_manual(values = c(mrsa_col, mssa_col))

pd = staph_resistances %>%
  filter(value == "R") %>%
  filter(variable %in% c("Fosfomycin", "Cotrimoxazole")) %>%
  ggplot() +
  geom_line(aes(date, prop, colour = SpeciesName)) +
  facet_wrap(~variable, nrow = 1) +
  theme_bw() +
  scale_y_continuous(limits = c(0,1)) +
  labs(x = "Time (months)", y = "Proportion of isolates resistant to antibiotic") +
  scale_x_date(limits = as.Date(c("2000-02-01", "2021-11-01"))) +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        strip.text = element_text(size = 12),
        legend.text = element_text(size = 12)) +
  scale_colour_manual(values = c(mrsa_col, mssa_col))

plot_grid(plot_grid(pa + theme(legend.position = "none"),
                    plot_grid(pb + theme(legend.position = "none"),
                              plot_grid(pc + theme(legend.position = "none"),
                                        pd + theme(legend.position = "none"),
                                        nrow = 1, labels = c("c)", "d)")),
                              nrow = 2, labels = c("b)", "")),
                    ncol = 2, rel_widths = c(0.4, 1.3), labels = c("a)", "")),
          get_legend(pa + theme(legend.position = "bottom")),
          nrow = 2, rel_heights = c(1,0.1))

ggsave(here::here("Figures", "fig4.png"), height = 10, width = 14)

sp1 = staph_resistances %>%
  filter(value == "R") %>%
  ggplot() +
  geom_line(aes(date, prop, colour = SpeciesName)) +
  facet_wrap(~variable) +
  theme_bw() +
  scale_y_continuous(limits = c(0,1)) +
  labs(x = "Time (months)", y = "Proportion of isolates resistant to antibiotic", colour = "") +
  scale_x_date(limits = as.Date(c("2000-02-01", "2021-11-01"))) +
  scale_colour_manual(values = c(mrsa_col, mssa_col)) +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(here::here("Figures", "suppfig5.png"), sp1, height = 12, width = 12)


#supp correlations
#remove resistances with less than 50 resistant isolates over all time period
library(Hmisc)
library(corrplot)

abx_dat = mrsa_resistances %>%
  filter(value == "R") %>%
  group_by(variable) %>%
  summarise(n = sum(n)) %>%
  filter(n > 50) %>%
  select(variable) %>% pull

cor_res = mrsa_resistances %>%
  filter(variable %in% abx_dat) %>%
  filter(value == "R") %>%
  dcast(., date~variable, value.var = "prop", fill = 0) %>%
  select(-"date")

cor_res = rcorr(as.matrix(cor_res))
cor_res$P[is.na(cor_res$P)] = 0.04
cor_res$P[which(cor_res$P < 0.05)] = 1
cor_res$P[which(cor_res$P != 1)] = 0

cor_res$r = cor_res$r * cor_res$P

png(here::here("Figures", "suppfig3.png"), width = 1100, height = 700)
corrplot(cor_res$r,
         method = "number", type = "upper", order = "hclust", tl.col = 'black', cl.cex = 1)
dev.off()


abx_dat = mssa_resistances %>%
  filter(value == "R") %>%
  group_by(variable) %>%
  summarise(n = sum(n)) %>%
  filter(n > 50) %>%
  select(variable) %>% pull

cor_res = mssa_resistances %>%
  filter(variable %in% abx_dat) %>%
  filter(value == "R") %>%
  dcast(., date~variable, value.var = "prop", fill = 0) %>%
  select(-"date")

cor_res = rcorr(as.matrix(cor_res))
cor_res$P[is.na(cor_res$P)] = 0.04
cor_res$P[which(cor_res$P < 0.05)] = 1
cor_res$P[which(cor_res$P != 1)] = 0

cor_res$r = cor_res$r * cor_res$P

png(here::here("Figures", "suppfig4.png"), width = 1100, height = 700)
corrplot(cor_res$r,
         method = "number", type = "upper", order = "hclust", tl.col = 'black', cl.cex = 1)
dev.off()
