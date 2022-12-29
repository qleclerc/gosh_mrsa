
library(dplyr)
library(lubridate)
library(reshape2)
library(ggplot2)
library(ggtext)
library(Hmisc)
library(corrplot)
library(cowplot)

mrsa_col = "#EFC000FF"
mssa_col = "#0073C2FF"

staph_isolates = read.csv(here::here("Clean", "staph_isolates.csv")) %>%
  mutate(date = as_date(date))

#number of resistances per isolate
staph_isolates$n_res = apply(staph_isolates[,-c(1:5)], 1, function(x) sum(x == "R", na.rm = T))
#number of resistance tested for per isolate
staph_isolates$n_test = apply(staph_isolates[,-c(1:5, 57)], 1, function(x) sum(!is.na(x)))

staph_isolates = staph_isolates[,c(1:5, 57, 58)]

#link between number of tests and resistances
# p1 = staph_isolates %>%
#   ggplot() +
#   geom_jitter(aes(n_test, n_res, colour = SpeciesName), alpha = 0.3) +
#   geom_smooth(data = staph_isolates %>% 
#                 filter(SpeciesName == "Methicillin-Resistant Staphylococcus aureus"),
#               aes(n_test, n_res), method = "lm", colour = "yellow4") +
#   geom_smooth(data = staph_isolates %>%
#                 filter(SpeciesName == "Methicillin-Susceptible Staphylococcus aureus"),
#               aes(n_test, n_res), method = "lm", colour = "blue4") +
#   geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
#   theme_bw() +
#   labs(x = "Number of susceptibility tests conducted", y = "Number of antibiotic resistances\ndetected",
#        colour = "") +
#   theme(axis.text = element_text(size = 12),
#         axis.title = element_text(size = 12),
#         legend.text = element_text(size = 12),
#         legend.position = "bottom") +
#   scale_x_continuous(breaks = seq(0,30,5)) +
#   coord_cartesian(ylim = c(0,20), xlim = c(0,30)) +
#   scale_colour_manual(values = c(mrsa_col, mssa_col))

p2 = staph_isolates %>%
  mutate(date = floor_date(date, "year")) %>%
  group_by(date, SpeciesName) %>%
  mutate(q25_test = quantile(n_test, 0.25)) %>%
  mutate(q75_test = quantile(n_test, 0.75)) %>%
  mutate(median_test = median(n_test)) %>%
  group_by(date, SpeciesName) %>%
  summarise(q25_test = mean(q25_test),
            q75_test = mean(q75_test),
            median_test = mean(median_test)) %>%
  ggplot() +
  geom_point(aes(date, median_test, colour = SpeciesName)) +
  geom_errorbar(aes(x = date, ymin = q25_test, ymax = q75_test, colour = SpeciesName),
                alpha = 0.5, size = 1) +
  theme_bw() +
  labs(x = "Time (years)", y = "Number of susceptibility tests\nconducted per isolate", colour = "") +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  scale_colour_manual(values = c(mrsa_col, mssa_col),
                      breaks = c("Methicillin-Resistant Staphylococcus aureus",
                                 "Methicillin-Susceptible Staphylococcus aureus"),
                      labels = c("Methicillin-Resistant *Staphylococcus aureus*",
                                 "Methicillin-Susceptible *Staphylococcus aureus*")) +
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=12),
        legend.text = element_markdown(size=12),
        legend.position = "bottom")

p3 = staph_isolates %>%
  mutate(date = floor_date(date, "year")) %>%
  group_by(date, SpeciesName) %>%
  mutate(q25_res = quantile(n_res, 0.25)) %>%
  mutate(q75_res = quantile(n_res, 0.75)) %>%
  mutate(median_res = median(n_res)) %>%
  group_by(date, SpeciesName) %>%
  summarise(q25_res = mean(q25_res),
            q75_res = mean(q75_res),
            median_res = mean(median_res)) %>%
  ggplot() +
  geom_point(aes(date, median_res, colour = SpeciesName)) +
  geom_errorbar(aes(x = date, ymin = q25_res, ymax = q75_res, colour = SpeciesName),
                alpha = 0.5, size = 1) +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.position = "bottom") +
  labs(x = "Time (years)", y = "Number of antibiotic resistances\ndetected per isolate") +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  scale_colour_manual(values = c(mrsa_col, mssa_col))


# #suppfig
# sp1 = staph_isolates %>%
#   filter(date < as.Date("2010-12-31")) %>%
#   ggplot() +
#   geom_jitter(aes(n_test, n_res, colour = SpeciesName), alpha = 0.3) +
#   geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
#   theme_bw() +
#   labs(x = "Number of susceptibility tests conducted", y = "Number of antibiotic resistances detected",
#        colour = "") +
#   scale_x_continuous(breaks = seq(0,30,5)) +
#   coord_cartesian(ylim = c(0,20), xlim = c(0,30)) +
#   scale_colour_manual(values = c(mrsa_col, mssa_col),
#                       breaks = c("Methicillin-Resistant Staphylococcus aureus",
#                                  "Methicillin-Susceptible Staphylococcus aureus"),
#                       labels = c("Methicillin-Resistant *Staphylococcus aureus*",
#                                  "Methicillin-Susceptible *Staphylococcus aureus*")) +
#   theme(axis.text = element_text(size=12),
#         axis.title = element_text(size=12),
#         legend.text = element_markdown(size=12),
#         legend.position = "bottom")
# 
# sp2 = staph_isolates %>%
#   filter(date > as.Date("2010-12-31")) %>%
#   ggplot() +
#   geom_jitter(aes(n_test, n_res, colour = SpeciesName), alpha = 0.3) +
#   geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
#   theme_bw() +
#   labs(x = "Number of susceptibility tests conducted", y = "Number of antibiotic resistances detected",
#        colour = "") +
#   theme(axis.text = element_text(size = 12),
#         axis.title = element_text(size = 12),
#         legend.text = element_text(size = 12),
#         legend.position = "bottom") +
#   scale_x_continuous(breaks = seq(0,30,5)) +
#   coord_cartesian(ylim = c(0,20), xlim = c(0,30)) +
#   scale_colour_manual(values = c(mrsa_col, mssa_col))
# 
# plot_grid(plot_grid(sp1+theme(legend.position = "none"), sp2+theme(legend.position = "none"),
#                     nrow = 1, labels = c("a)", "b)"), hjust = 0),
#           get_legend(sp1), nrow = 2, rel_heights = c(1,0.1))
# 
# ggsave(here::here("Figures", "suppfig3.png"), height = 5, width = 10, dpi = 300)




#descriptive stats

#how many isolates with no reportedresistance?
sum(staph_isolates$n_test == 0)
#how many patients?
(staph_isolates %>%
    select(project_id) %>%
    pull %>%
    unique %>%
    length) - 
  (staph_isolates %>%
     filter(n_test != 0) %>%
     select(project_id) %>%
     pull %>%
     unique %>%
     length)
sum(staph_isolates$n_test[staph_isolates$SpeciesName == "Methicillin-Resistant Staphylococcus aureus"] == 0)
sum(staph_isolates$n_test[staph_isolates$SpeciesName == "Methicillin-Susceptible Staphylococcus aureus"] == 0)

nosusnose = staph_isolates %>% filter(n_test == 0) %>% select(SpecimenType) %>% pull
nosusnose = length((grep("ose|hroat|asal", nosusnose, value = T)))
nosusnonose = 10029-nosusnose
nosusnose/(nosusnonose+nosusnose)

susnose = staph_isolates %>% filter(n_test > 0) %>% select(SpecimenType) %>% pull
susnose = length((grep("ose|hroat|asal", susnose, value = T)))
susnonose = 62178 - susnose
susnose/(susnonose+susnose)

#isolates with/without sus test data x isolates from nose,throat,skin/not
chisq.test(matrix(c(susnose,nosusnose,
                    susnonose,nosusnonose),
                  byrow = T, nrow = 2))


#how many isolates with n reported resistance == n tested resistance
sum(staph_isolates$n_test == staph_isolates$n_res)
sum(staph_isolates$n_test[staph_isolates$SpeciesName == "Methicillin-Resistant Staphylococcus aureus"] == staph_isolates$n_res[staph_isolates$SpeciesName == "Methicillin-Resistant Staphylococcus aureus"])
sum(staph_isolates$n_test[staph_isolates$SpeciesName == "Methicillin-Susceptible Staphylococcus aureus"] == staph_isolates$n_res[staph_isolates$SpeciesName == "Methicillin-Susceptible Staphylococcus aureus"])


#linear regression
summary(lm(data = staph_isolates,
           n_res ~ n_test + SpeciesName))

summary(lm(data = staph_isolates %>%
             mutate(date = floor_date(date, "year")) %>%
             filter(date < as.Date("2011-01-01")),
           n_test ~ date + SpeciesName))
summary(lm(data = staph_isolates %>%
             mutate(date = floor_date(date, "year")) %>%
             filter(date >= as.Date("2011-01-01")),
           n_test ~ date + SpeciesName))

summary(lm(data = staph_isolates %>%
             mutate(date = floor_date(date, "year")) %>%
             filter(date >= as.Date("2011-01-01")),
           n_res ~ date + SpeciesName))




#top resistances for MSSA and MRSA
#second part of fig2, compile and save
staph_isolates = read.csv(here::here("Clean", "staph_isolates.csv")) %>%
  mutate(date = as_date(date))

all_data = data.frame()

d1 = staph_isolates %>%
  filter(SpeciesName == "Methicillin-Susceptible Staphylococcus aureus") %>%
  select(-c(1:5)) %>%
  apply(., 2, function(x) round(sum(!is.na(x))/51020, 2)) %>%
  sort(decreasing = T) %>%
  .[1:15]

all_data = rbind(all_data,
                 data.frame(abx = names(d1),
                            val = d1,
                            spec = "MSSA",
                            var = "test"))

d1 = staph_isolates %>%
  filter(SpeciesName == "Methicillin-Susceptible Staphylococcus aureus") %>%
  select(-c(1:5)) %>%
  apply(., 2, function(x) round(sum(x == "R", na.rm = T)/51020, 2)) %>%
  .[all_data %>% filter(spec == "MSSA" & var == "test") %>% select(abx) %>% pull]

all_data = rbind(all_data,
                 data.frame(abx = names(d1),
                            val = d1,
                            spec = "MSSA",
                            var = "res"))

d1 = staph_isolates %>%
  filter(SpeciesName == "Methicillin-Resistant Staphylococcus aureus") %>%
  select(-c(1:5)) %>%
  apply(., 2, function(x) round(sum(!is.na(x))/21187, 2)) %>%
  sort(decreasing = T) %>%
  .[1:15]

all_data = rbind(all_data,
                 data.frame(abx = names(d1),
                            val = d1,
                            spec = "MRSA",
                            var = "test"))

d1 = staph_isolates %>%
  filter(SpeciesName == "Methicillin-Resistant Staphylococcus aureus") %>%
  select(-c(1:5)) %>%
  apply(., 2, function(x) round(sum(x == "R", na.rm = T)/21187, 2)) %>%
  .[all_data %>% filter(spec == "MRSA" & var == "test") %>% select(abx) %>% pull]

all_data = rbind(all_data,
                 data.frame(abx = names(d1),
                            val = d1,
                            spec = "MRSA",
                            var = "res"))


p4 = ggplot() +
  geom_col(data = all_data %>%
             filter(spec == "MRSA" & var == "test") %>%
             mutate(abx = factor(abx, levels = abx[rev(order(val))])),
           aes(abx, val), fill = mrsa_col, alpha = 0.5) +
  geom_col(data = all_data %>%
             filter(spec == "MRSA" & var == "res"),
           aes(abx, val), fill = mrsa_col) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        axis.title = element_text(size = 12)) +
  labs(x = "", y = "Proportion of MRSA\nisolates tested and resistant") +
  scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.2))

p5 = ggplot() +
  geom_col(data = all_data %>%
             filter(spec == "MSSA" & var == "test") %>%
             mutate(abx = factor(abx, levels = abx[rev(order(val))])),
           aes(abx, val), fill = mssa_col, alpha = 0.5) +
  geom_col(data = all_data %>%
             filter(spec == "MSSA" & var == "res"),
           aes(abx, val), fill = mssa_col) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        axis.title = element_text(size = 12)) +
  labs(x = "", y = "Proportion of MSSA\nisolates tested and resistant") +
  scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.2))


pall = plot_grid(plot_grid(p2 + theme(legend.position = "none"),
                           NULL,
                           p3 + theme(legend.position = "none"),
                           ncol = 3,
                           labels = c("a)", "", "b)"),
                           rel_widths = c(1,0.05,1),
                           label_size = 12,
                           hjust = 0),
                 get_legend(p2),
                 p4,
                 p5,
                 rel_heights = c(1.1,0.1,1,1),
                 labels = c("", "", "c)", "d)"),
                 nrow = 4,
                 label_size = 12,
                 hjust = 0,
                 vjust = 0)

ggsave(here::here("Figures", "fig2.png"), pall, height = 10, width = 12, dpi = 300)

