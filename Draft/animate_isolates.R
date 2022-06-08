
library(dplyr)
library(lubridate)
library(reshape2)
library(ggplot2)
library(gganimate)

mrsa_col = "#EFC000FF"
mssa_col = "#0073C2FF"

staph_isolates = read.csv(here::here("Clean", "staph_isolates.csv")) %>%
  mutate(date = as_date(date))

#number of resistances per isolate
staph_isolates$n_res = apply(staph_isolates[,-c(1:5)], 1, function(x) sum(x == "R", na.rm = T))
#number of resistance tested for per isolate
staph_isolates$n_test = apply(staph_isolates[,-c(1:5, 57)], 1, function(x) sum(!is.na(x)))

staph_isolates = staph_isolates[,c(1:5, 57, 58)]

staph_isolates$year = as.numeric(year(staph_isolates$date))

p = ggplot(staph_isolates) +
  geom_jitter(aes(n_test, n_res, colour = SpeciesName), alpha = 0.3) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  theme_bw() +
  labs(x = "Number of susceptibility tests conducted", y = "Number of antibiotic resistances detected",
       colour = "", title = "Year: {frame_time}") +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.position = "bottom") +
  scale_x_continuous(breaks = seq(0,30,5)) +
  coord_cartesian(ylim = c(0,20), xlim = c(0,30)) +
  scale_colour_manual(values = c(mrsa_col, mssa_col)) +
  transition_time(year)

animate(p, duration = 30, fps = 10)

