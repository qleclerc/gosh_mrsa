
library(dplyr)
library(lubridate)
library(reshape2)
library(ggplot2)
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

#significant correlation, but r^2 = 0.572 suggesting not all var in n res explained by n tests 
summary(lm(data = staph_isolates,
           n_res ~ n_test + SpeciesName))

#link between number of tests and resistances
p1 = staph_isolates %>%
  ggplot() +
  geom_jitter(aes(n_test, n_res, colour = SpeciesName), alpha = 0.3) +
  geom_smooth(data = staph_isolates %>% 
                filter(SpeciesName == "Methicillin-Resistant Staphylococcus aureus"),
              aes(n_test, n_res), method = "lm", colour = "yellow4") +
  geom_smooth(data = staph_isolates %>%
                filter(SpeciesName == "Methicillin-Susceptible Staphylococcus aureus"),
              aes(n_test, n_res), method = "lm", colour = "blue4") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  theme_bw() +
  labs(x = "Number of susceptibility tests conducted", y = "Number of antibiotic resistances detected",
       colour = "") +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.position = "bottom") +
  scale_x_continuous(breaks = seq(0,30,5)) +
  coord_cartesian(ylim = c(0,20), xlim = c(0,30)) +
  scale_colour_manual(values = c(mrsa_col, mssa_col))


#how many isolates with no reported resistance?
sum(staph_isolates$n_test == 0)
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
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.position = "bottom") +
  labs(x = "Time (years)", y = "Number of susceptibility tests\nconducted per isolate", colour = "") +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  scale_colour_manual(values = c(mrsa_col, mssa_col))

#significant correlation, but r^2 = 0.572 suggesting not all var in n res explained by n tests 
summary(lm(data = staph_isolates %>%
             mutate(date = floor_date(date, "year")) %>%
             filter(date < as.Date("2011-01-01")),
           n_test ~ date + SpeciesName))

summary(lm(data = staph_isolates %>%
             mutate(date = floor_date(date, "year")) %>%
             filter(date >= as.Date("2011-01-01")),
           n_test ~ date + SpeciesName))


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

summary(lm(data = staph_isolates %>%
             mutate(date = floor_date(date, "year")) %>%
             filter(date >= as.Date("2011-01-01")),
           n_res ~ date + SpeciesName))

pall = plot_grid(p1 + theme(legend.position = "none"),
                 NULL,
                 plot_grid(p2 + theme(legend.position = "none"),
                           NULL,
                           p3 + theme(legend.position = "none"),
                           ncol = 3,
                           labels = c("b)", "", "c)"),
                           rel_widths = c(1,0.05,1),
                           label_size = 12,
                           hjust = 0),
                 get_legend(p2),
                 rel_heights = c(1.2,0.05,1.1,0.1),
                 labels = c("a)"),
                 nrow = 4,
                 label_size = 12,
                 hjust = 0)

ggsave(here::here("Figures", "fig2.png"), pall, height = 9, width = 12, dpi = 300)



sp1 = staph_isolates %>%
  filter(date < as.Date("2010-12-31")) %>%
  ggplot() +
  geom_jitter(aes(n_test, n_res, colour = SpeciesName), alpha = 0.3) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  theme_bw() +
  labs(x = "Number of susceptibility tests conducted", y = "Number of antibiotic resistances detected",
       colour = "") +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.position = "bottom") +
  scale_x_continuous(breaks = seq(0,30,5)) +
  coord_cartesian(ylim = c(0,20), xlim = c(0,30)) +
  scale_colour_manual(values = c(mrsa_col, mssa_col))

sp2 = staph_isolates %>%
  filter(date > as.Date("2010-12-31")) %>%
  ggplot() +
  geom_jitter(aes(n_test, n_res, colour = SpeciesName), alpha = 0.3) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  theme_bw() +
  labs(x = "Number of susceptibility tests conducted", y = "Number of antibiotic resistances detected",
       colour = "") +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.position = "bottom") +
  scale_x_continuous(breaks = seq(0,30,5)) +
  coord_cartesian(ylim = c(0,20), xlim = c(0,30)) +
  scale_colour_manual(values = c(mrsa_col, mssa_col))

plot_grid(plot_grid(sp1+theme(legend.position = "none"), sp2+theme(legend.position = "none"),
                    nrow = 1, labels = c("a)", "b)"), hjust = 0),
          get_legend(sp1), nrow = 2, rel_heights = c(1,0.1))

ggsave(here::here("Figures", "suppfig2.png"), height = 5, width = 10, dpi = 300)
