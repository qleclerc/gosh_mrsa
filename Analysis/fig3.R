
library(dplyr)
library(lubridate)
library(reshape2)
library(ggplot2)
library(scales)
library(Hmisc)
library(corrplot)
library(cowplot)

staph_isolates = read.csv(here::here("Clean", "staph_isolates.csv"), stringsAsFactors = F) %>%
  mutate(date = as_date(date))

#the line to remove all isolates with NAs only (ie nothing tested at all)
#staph_isolates = staph_isolates[!apply(staph_isolates[,6:59], 1, function(x) all(is.na(x))),]

mrsa_isolates_m = staph_isolates %>%
  filter(SpeciesName == "Methicillin-Resistant Staphylococcus aureus") %>%
  select(-c(1:4)) %>%
  melt(id.vars = "date") %>%
  filter(!is.na(value)) %>%
  mutate(date = floor_date(date, "month")) %>%
  group_by(date, variable) %>%
  summarise(n = n()) %>%
  filter(variable %in% c("Chloramphenicol", "Mupirocin", "Tetracycline",
                         "Teicoplanin", "Vancomycin", "Syncercid", "Septrin",
                         "Linezolid", "Clindamycin", "Fosfomycin", "Cefoxitin",
                         "Co.Trimoxazole..Septrin."))

mssa_isolates_m = staph_isolates %>%
  filter(SpeciesName == "Methicillin-Susceptible Staphylococcus aureus") %>%
  select(-c(1:4)) %>%
  melt(id.vars = "date") %>%
  filter(!is.na(value)) %>%
  mutate(date = floor_date(date, "month")) %>%
  group_by(date, variable) %>%
  summarise(n = n()) %>%
  filter(variable %in% c("Chloramphenicol", "Mupirocin", "Tetracycline",
                         "Teicoplanin", "Vancomycin", "Syncercid", "Septrin",
                         "Linezolid", "Clindamycin", "Fosfomycin", "Cefoxitin",
                         "Co.Trimoxazole..Septrin."))

mrsa_per_day = staph_isolates %>%
  filter(SpeciesName == "Methicillin-Resistant Staphylococcus aureus") %>%
  mutate(date = floor_date(date, "month")) %>%
  group_by(date) %>%
  summarise(total = n())

mrsa_isolates_m = left_join(mrsa_isolates_m, mrsa_per_day, by = "date")
mrsa_isolates_m$species = "Methicillin-Resistant Staphylococcus aureus"

mssa_per_day = staph_isolates %>%
  filter(SpeciesName == "Methicillin-Susceptible Staphylococcus aureus") %>%
  mutate(date = floor_date(date, "month")) %>%
  group_by(date) %>%
  summarise(total = n())

mssa_isolates_m = left_join(mssa_isolates_m, mssa_per_day, by = "date")
mssa_isolates_m$species = "Methicillin-Susceptible Staphylococcus aureus"

all_isolates = rbind(mrsa_isolates_m, mssa_isolates_m)

p1 = ggplot(all_isolates %>% 
              filter(variable %in% c("Teicoplanin", "Vancomycin", "Syncercid", "Linezolid"))) +
  geom_line(aes(date, n/total, colour = species)) +
  scale_y_continuous(limits = c(0,1)) +
  facet_wrap(~variable, ncol = 1) +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        strip.text = element_text(size = 12),
        legend.text = element_text(size = 12)) +
  labs(x = "Time", y = "Proportion of isolates tested for resistance", colour = "") +
  scale_x_date(limits = as.Date(c("2000-02-01", "2021-11-01")))

p2 = ggplot(all_isolates %>%
              filter(variable %in% c("Septrin", "Fosfomycin"))) +
  geom_line(aes(date, n/total, colour = species)) +
  scale_y_continuous(limits = c(0,1)) +
  facet_wrap(~variable, ncol = 1) +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        strip.text = element_text(size = 12),
        legend.text = element_text(size = 12)) +
  labs(x = "Time", y = "Proportion of isolates tested for resistance") +
  scale_x_date(limits = as.Date(c("2000-02-01", "2021-11-01")))

p3 = ggplot(all_isolates %>% 
              filter(variable %in% c("Mupirocin", "Chloramphenicol", "Tetracycline")) %>%
              mutate(variable = factor(variable, levels = c("Mupirocin", "Chloramphenicol", "Tetracycline")))) +
  geom_line(aes(date, n/total, colour = species)) +
  scale_y_continuous(limits = c(0,1)) +
  facet_wrap(~variable) +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        strip.text = element_text(size = 12),
        legend.text = element_text(size = 12)) +
  labs(x = "Time", y = "Proportion of isolates tested for resistance") +
  scale_x_date(limits = as.Date(c("2000-02-01", "2021-11-01")))

p4 = ggplot(all_isolates %>% 
              filter(variable %in% c("Clindamycin", "Cefoxitin", "Co.Trimoxazole..Septrin."))) +
  geom_line(aes(date, n/total, colour = species)) +
  scale_y_continuous(limits = c(0,1)) +
  facet_wrap(~variable) +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        strip.text = element_text(size = 12),
        legend.text = element_text(size = 12)) +
  labs(x = "Time", y = "Proportion of isolates tested for resistance") +
  scale_x_date(limits = as.Date(c("2000-02-01", "2021-11-01")))


plot_grid(plot_grid(plot_grid(p1 + theme(legend.position = "none"),
                              NULL,
                              p2 + theme(legend.position = "none"),
                              ncol = 3, rel_widths = c(1,0.05,1),
                              labels = c("a)", "", "b)"), hjust = 0, label_size = 12),
                    NULL,
                    plot_grid(p3 + theme(legend.position = "none"),
                              NULL,
                              p4 + theme(legend.position = "none"),
                              ncol = 1, rel_heights = c(1, 0.05, 1),
                              labels = c("c)", "", "d)"), hjust = 0, label_size = 12),
                    ncol = 3, rel_widths = c(1,0.05,1.2)),
          get_legend(p1+theme(legend.position = "bottom")),
          ncol = 1, rel_heights = c(1,0.05))

ggsave(here::here("Figures","fig3.png"), height = 10, width = 14)

# #cor
# cor_res = staph_isolates_m %>%
#   dcast(., date~variable, value.var = "n", fill = 0) %>%
#   select(-"date")
# 
# cor_res = rcorr(as.matrix(cor_res))
# cor_res$P[is.na(cor_res$P)] = 0.04
# cor_res$P[which(cor_res$P < 0.05)] = 1
# cor_res$P[which(cor_res$P != 1)] = 0
# 
# cor_res$r = cor_res$r * cor_res$P
# 
# png(here::here("Figures", "suppfig3.png"), width = 1100, height = 701)
# corrplot(cor_res$r,
#          method = "number", type = "upper", order = "hclust", tl.col = 'black', cl.cex = 1)
# dev.off()
# 
# 
# 
# #mrsa
# staph_isolates_m = staph_isolates %>%
#   filter(SpeciesName == "Methicillin-Resistant Staphylococcus aureus") %>%
#   #filter(SpeciesName == "Methicillin-Susceptible Staphylococcus aureus") %>%
#   select(-c(1:4)) %>%
#   melt(id.vars = "date") %>%
#   filter(!is.na(value)) %>%
#   mutate(date = floor_date(date, "month")) %>%
#   group_by(date, variable) %>%
#   summarise(n = n()) %>%
#   filter(variable %in% c("Amikacin", "Amik.Fluclox", "Chloramphenicol", "Erythromycin",
#                          "Ciprofloxacin", "Flucloxacillin", "Fucidin", "Gentamicin",
#                          "Gent.Cipro", "Linezolid", "Mupirocin","Penicillin", "Rifampicin",
#                          "Teicoplanin", "Tetracycline", "Trimethoprim", "Vancomycin"))
# 
# isolates_per_day = staph_isolates %>%
#   filter(SpeciesName == "Methicillin-Resistant Staphylococcus aureus") %>%
#   #filter(SpeciesName == "Methicillin-Susceptible Staphylococcus aureus") %>%
#   mutate(date = floor_date(date, "month")) %>%
#   group_by(date) %>%
#   summarise(total = n())
# 
# staph_isolates_m = left_join(staph_isolates_m, isolates_per_day, by = "date")
# 
# p4 = ggplot(staph_isolates_m %>%
#               filter(variable %in% c("Mupirocin", "Teicoplanin", "Vancomycin", "Linezolid"))) +
#   geom_line(aes(date, n/total, colour = variable)) +
#   scale_y_continuous(limits = c(0,1)) +
#   theme_bw() +
#   labs(x = "Date", y = "Proportion of S aureus isolates tested for resistance", colour = "Antibiotic")
# 
# p2 = ggplot(staph_isolates_m %>% 
#               filter(variable %in% c("Amikacin", "Amik.Fluclox",
#                                      "Trimethoprim","Fucidin", "Gentamicin",
#                                      "Gent.Cipro", "Rifampicin"))) +
#   geom_line(aes(date, n/total, colour = variable)) +
#   scale_y_continuous(limits = c(0,1)) +
#   theme_bw() +
#   labs(x = "Date", y = "Proportion of S aureus isolates tested for resistance", colour = "Antibiotic")
# 
# p3 = ggplot(staph_isolates_m %>% 
#               filter(variable %in% c( "Chloramphenicol","Tetracycline"))) +
#   geom_line(aes(date, n/total, colour = variable)) +
#   scale_y_continuous(limits = c(0,1)) +
#   theme_bw() +
#   labs(x = "Date", y = "Proportion of S aureus isolates tested for resistance", colour = "Antibiotic")
# 
# p1 = ggplot(staph_isolates_m %>% 
#               filter(variable %in% c("Erythromycin",
#                                      "Ciprofloxacin", "Flucloxacillin",
#                                      "Penicillin"))) +
#   geom_line(aes(date, n/total, colour = variable)) +
#   scale_y_continuous(limits = c(0,1)) +
#   theme_bw() +
#   labs(x = "Date", y = "Proportion of S aureus isolates tested for resistance", colour = "Antibiotic")
# 
# plot_grid(p1, p2, p3, p4,
#           ncol = 2)
# 
# ggsave(here::here("Figures","suppfig1.png"), height = 10, width = 12)
# 
# #mssa
# staph_isolates_m = staph_isolates %>%
#   #filter(SpeciesName == "Methicillin-Resistant Staphylococcus aureus") %>%
#   filter(SpeciesName == "Methicillin-Susceptible Staphylococcus aureus") %>%
#   select(-c(1:4)) %>%
#   melt(id.vars = "date") %>%
#   filter(!is.na(value)) %>%
#   mutate(date = floor_date(date, "month")) %>%
#   group_by(date, variable) %>%
#   summarise(n = n()) %>%
#   filter(variable %in% c("Amikacin", "Amik.Fluclox", "Chloramphenicol", "Erythromycin",
#                          "Ciprofloxacin", "Flucloxacillin", "Fucidin", "Gentamicin",
#                          "Gent.Cipro", "Linezolid", "Mupirocin","Penicillin", "Rifampicin",
#                          "Teicoplanin", "Tetracycline", "Trimethoprim", "Vancomycin"))
# 
# isolates_per_day = staph_isolates %>%
#   #filter(SpeciesName == "Methicillin-Resistant Staphylococcus aureus") %>%
#   filter(SpeciesName == "Methicillin-Susceptible Staphylococcus aureus") %>%
#   mutate(date = floor_date(date, "month")) %>%
#   group_by(date) %>%
#   summarise(total = n())
# 
# staph_isolates_m = left_join(staph_isolates_m, isolates_per_day, by = "date")
# 
# p4 = ggplot(staph_isolates_m %>%
#               filter(variable %in% c("Mupirocin", "Teicoplanin", "Vancomycin", "Linezolid"))) +
#   geom_line(aes(date, n/total, colour = variable)) +
#   scale_y_continuous(limits = c(0,1)) +
#   theme_bw() +
#   labs(x = "Date", y = "Proportion of S aureus isolates tested for resistance", colour = "Antibiotic")
# 
# p2 = ggplot(staph_isolates_m %>% 
#               filter(variable %in% c("Amikacin", "Amik.Fluclox",
#                                      "Trimethoprim","Fucidin", "Gentamicin",
#                                      "Gent.Cipro", "Rifampicin"))) +
#   geom_line(aes(date, n/total, colour = variable)) +
#   scale_y_continuous(limits = c(0,1)) +
#   theme_bw() +
#   labs(x = "Date", y = "Proportion of S aureus isolates tested for resistance", colour = "Antibiotic")
# 
# p3 = ggplot(staph_isolates_m %>% 
#               filter(variable %in% c( "Chloramphenicol","Tetracycline"))) +
#   geom_line(aes(date, n/total, colour = variable)) +
#   scale_y_continuous(limits = c(0,1)) +
#   theme_bw() +
#   labs(x = "Date", y = "Proportion of S aureus isolates tested for resistance", colour = "Antibiotic")
# 
# p1 = ggplot(staph_isolates_m %>% 
#               filter(variable %in% c("Erythromycin",
#                                      "Ciprofloxacin", "Flucloxacillin",
#                                      "Penicillin"))) +
#   geom_line(aes(date, n/total, colour = variable)) +
#   scale_y_continuous(limits = c(0,1)) +
#   theme_bw() +
#   labs(x = "Date", y = "Proportion of S aureus isolates tested for resistance", colour = "Antibiotic")
# 
# plot_grid(p1, p2, p3, p4,
#           ncol = 2)
# 
# ggsave(here::here("Figures","suppfig2.png"), height = 10, width = 12)
# 
# 
# #other abx
# staph_isolates_m = staph_isolates %>%
#   #filter(SpeciesName == "Methicillin-Resistant Staphylococcus aureus") %>%
#   #filter(SpeciesName == "Methicillin-Susceptible Staphylococcus aureus") %>%
#   select(-c(1:4)) %>%
#   melt(id.vars = "date") %>%
#   filter(!is.na(value)) %>%
#   mutate(date = floor_date(date, "month")) %>%
#   group_by(date, variable) %>%
#   summarise(n = n()) %>%
#   filter(!(variable %in% c("Amikacin", "Amik.Fluclox", "Chloramphenicol", "Erythromycin",
#                            "Ciprofloxacin", "Flucloxacillin", "Fucidin", "Gentamicin",
#                            "Gent.Cipro", "Linezolid", "Mupirocin","Penicillin", "Rifampicin",
#                            "Teicoplanin", "Tetracycline", "Trimethoprim", "Vancomycin")))
# 
# isolates_per_day = staph_isolates %>%
#   #filter(SpeciesName == "Methicillin-Resistant Staphylococcus aureus") %>%
#   #filter(SpeciesName == "Methicillin-Susceptible Staphylococcus aureus") %>%
#   mutate(date = floor_date(date, "month")) %>%
#   group_by(date) %>%
#   summarise(total = n())
# 
# staph_isolates_m = left_join(staph_isolates_m, isolates_per_day, by = "date")
# 
# ggplot(staph_isolates_m) +
#   geom_line(aes(date, n/total, colour = variable)) +
#   scale_y_continuous(limits = c(0,1)) +
#   theme_bw() +
#   labs(x = "Date", y = "Proportion of S aureus isolates tested for resistance", colour = "Antibiotic")
# 
# ggsave(here::here("Figures", "suppfigx.png"))
# 
# 
# 
