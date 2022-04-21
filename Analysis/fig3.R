
library(dplyr)
library(lubridate)
library(ggplot2)
library(scales)
library(cowplot)

staph_isolates = read.csv(here::here("Clean", "staph_isolates.csv")) %>%
  mutate(date = as_date(date))


#top resistances for MSSA and MRSA
d1 = staph_isolates %>%
  filter(SpeciesName == "Methicillin-Susceptible Staphylococcus aureus") %>%
  select(-c(1:5)) %>%
  apply(., 2, function(x) round(sum(!is.na(x))/51020*100, 2)) %>%
  sort(decreasing = T) %>%
  .[1:10]

p1 = data.frame(abx = factor(names(d1), levels = names(d1)),
                val = d1) %>%
  ggplot() +
  geom_col(aes(abx, val), fill = "#00BFC4") +
  theme_bw() +
  theme(axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1,
                                   face = c("bold", "plain", "plain", "plain",
                                            "bold", "bold", "plain", "plain",
                                            "plain", "plain")),
        axis.title = element_text(size = 12)) +
  labs(x = "", y = "Percentage of MSSA\nisolates tested") +
  scale_y_continuous(limits = c(0,100), breaks = seq(0,100,20))

d1 = staph_isolates %>%
  filter(SpeciesName == "Methicillin-Susceptible Staphylococcus aureus") %>%
  select(-c(1:5)) %>%
  apply(., 2, function(x) round(sum(x == "R", na.rm = T)/51020*100, 2)) %>%
  sort(decreasing = T) %>%
  .[1:10]

p2 = data.frame(abx = factor(names(d1), levels = names(d1)),
                val = d1) %>%
  ggplot() +
  geom_col(aes(abx, val), fill = "#00BFC4") +
  theme_bw() +
  theme(axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1,
                                   face = c("plain", "plain", "plain", "plain",
                                            "plain", "bold", "plain", "bold",
                                            "plain", "bold")),
        axis.title = element_text(size = 12)) +
  labs(x = "", y = "Percentage of MSSA\nisolates resistant") +
  scale_y_continuous(limits = c(0,100), breaks = seq(0,100,20))

d1 = staph_isolates %>%
  filter(SpeciesName == "Methicillin-Resistant Staphylococcus aureus") %>%
  select(-c(1:5)) %>%
  apply(., 2, function(x) round(sum(!is.na(x))/21187*100, 2)) %>%
  sort(decreasing = T) %>%
  .[1:10]

p3 = data.frame(abx = factor(names(d1), levels = names(d1)),
                val = d1) %>%
  ggplot() +
  geom_col(aes(abx, val), fill = "#F8766D") +
  theme_bw() +
  theme(axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1,
                                   face = c("plain", "plain", "plain", "plain",
                                            "bold", "plain", "plain", "plain",
                                            "plain", "bold")),
        axis.title = element_text(size = 12)) +
  labs(x = "", y = "Percentage of MRSA\nisolates tested") +
  scale_y_continuous(limits = c(0,100), breaks = seq(0,100,20))

d1 = staph_isolates %>%
  filter(SpeciesName == "Methicillin-Resistant Staphylococcus aureus") %>%
  select(-c(1:5)) %>%
  apply(., 2, function(x) round(sum(x == "R", na.rm = T)/21187*100, 2)) %>%
  sort(decreasing = T) %>%
  .[1:10]

p4 = data.frame(abx = factor(names(d1), levels = names(d1)),
                val = d1) %>%
  ggplot() +
  geom_col(aes(abx, val), fill = "#F8766D") +
  theme_bw() +
  theme(axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1,
                                   face = c("plain", "plain", "plain", "plain",
                                            "bold", "plain", "plain", "plain",
                                            "bold", "plain")),
        axis.title = element_text(size = 12)) +
  labs(x = "", y = "Percentage of MRSA\nisolates resistant") +
  scale_y_continuous(limits = c(0,100), breaks = seq(0,100,20))

plot_grid(p1, NULL, p2,
          NULL, NULL, NULL, 
          p3, NULL, p4,
          rel_heights = c(1,0.05,1),
          rel_widths = c(1,0.05,1),
          ncol = 3,
          labels = c("a)", "", "b)",
                     "", "", "",
                     "c)","", "d)"))

ggsave(here::here("Figures", "fig3.png"), height = 8, width = 10)
