
library(dplyr)
library(lubridate)
library(ggplot2)
library(scales)
library(cowplot)

mrsa_col = "#EFC000FF"
mssa_col = "#0073C2FF"

staph_isolates = read.csv(here::here("Clean", "staph_isolates.csv")) %>%
  mutate(date = as_date(date))

all_data = data.frame()

#top resistances for MSSA and MRSA
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


p1 = ggplot() +
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
  labs(x = "", y = "Proportion of MRSA isolates\ntested and resistant") +
  scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.2))

p2 = ggplot() +
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
  labs(x = "", y = "Proportion of MSSA isolates\ntested and resistant") +
  scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.2))

plot_grid(p1, p2,
          rel_heights = c(1,1),
          ncol = 1,
          labels = c("a)","b)"))

ggsave(here::here("Figures", "fig3.png"), height = 8, width = 10)
