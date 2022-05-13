#time-series analysis of s aureus cases
#tldr: no seasonality

library(tidyr)
library(dplyr)
library(lubridate)
library(ggplot2)
library(scales)
library(cowplot)
library(forecast)
library(zoo)

mrsa_col = "#EFC000FF"
mssa_col = "#0073C2FF"

staph_isolates = read.csv(here::here("Clean", "staph_isolates.csv")) %>% 
  mutate(date = as_date(date))

#some general plots
#baseline incidence
p1 = staph_isolates %>%
  mutate(date = floor_date(date, "month")) %>%
  # group_by(project_id, date, SpeciesName) %>%
  # summarise(n = 1) %>%
  # ungroup() %>%
  group_by(date) %>%
  summarise(n = n()) %>%
  ggplot() +
  geom_line(aes(date, n), size = 0.8) +
  theme_bw() +
  labs(y = "Number of S. aureus isolates", x = "Time (months)") +
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=12),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12)) +
  scale_x_date(limits = as.Date(c("2000-02-01", "2021-11-01")),
               date_breaks = "2 years", date_labels = "%Y") +
  geom_vline(xintercept = as.Date("2020-03-26"), linetype = 2, colour = "green4", size = 1)



#prop MRSA/MSSA
p2 = staph_isolates %>%
  mutate(date = floor_date(date, "month")) %>%
  # group_by(project_id, date, SpeciesName) %>%
  # summarise(n = 1) %>%
  # ungroup() %>%
  group_by(date, SpeciesName) %>%
  summarise(n = n()) %>%
  mutate(prop = n/sum(n)) %>%
  ggplot() +
  geom_line(aes(date, prop, colour = SpeciesName), size = 0.8) +
  theme_bw() +
  labs(y = "Proportion of MSSA and MRSA isolates", x = "Time (months)", colour = "") +
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=12),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12)) +
  scale_x_date(limits = as.Date(c("2000-02-01", "2021-11-01")),
               date_breaks = "2 years", date_labels = "%Y") +
  geom_vline(xintercept = as.Date("2020-03-26"), linetype = 2, colour = "green4", size = 1) +
  scale_colour_manual(values = c(mrsa_col, mssa_col))

p3 = staph_isolates %>%
  mutate(date = floor_date(date, "year")) %>%
  group_by(date, SpeciesName, project_id) %>%
  summarise(n = n()) %>%
  ggplot() +
  geom_boxplot(aes(date, n, group = interaction(date, SpeciesName), colour = SpeciesName),
               outlier.shape = NA, size = 0.8) +
  theme_bw() +
  scale_y_continuous(limits = c(0,10), breaks = seq(0,10,2)) +
  labs(x = "Time (years)", y = "Number of isolates per patient", colour = "") +
  theme(legend.position = "bottom",
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 12)) +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  scale_colour_manual(values = c(mrsa_col, mssa_col))

plot_grid(plot_grid(p1 + theme(legend.position = "none"),
                    NULL,
                    p2 + theme(legend.position = "none"),
                    NULL,
                    p3 + theme(legend.position = "none"),
                    nrow = 5, labels = c("a)", "", "b)", "", "c)"), label_size = 12,
                    rel_heights = c(0.7, 0.05, 1, 0.05, 0.5), vjust = c(1,0,0)),
          get_legend(p2 + theme(legend.position = "bottom")),
          nrow = 2, rel_heights = c(1,0.05))

ggsave(here::here("Figures", "fig1.png"), heigh = 14, width = 10)


# #descriptive stats
cat(length(unique(staph_isolates$project_id)), "patients with a positive S aureus isolate between",
    as.character(min(staph_isolates$date)), "and", as.character(max(staph_isolates$date)))
cat(staph_isolates %>%
      filter(SpeciesName == "Methicillin-Resistant Staphylococcus aureus") %>%
      select(project_id) %>% pull %>% unique %>% length,
    "patients with a positive MRSA isolate")
cat(staph_isolates %>%
      filter(SpeciesName == "Methicillin-Susceptible Staphylococcus aureus") %>%
      select(project_id) %>% pull %>% unique %>% length,
    "patients with a positive MSSA isolate")
(19777+3506-22206)
19777 - 1077
3506 - 1077

(19780+3519-22206)/22206
(19777 - 1077)/22206
(3506 - 1077)/22206

table(staph_isolates$SpeciesName)

#lm
lm_dat = staph_isolates %>%
  mutate(date = floor_date(date, "year")) %>%
  group_by(date, SpeciesName, project_id) %>%
  summarise(n = n()) %>%
  select(date, SpeciesName, n)

lm_res = lm(n ~ date + SpeciesName, data = lm_dat)
summary(lm_res)


# 
# #periods in msts: the period defines how often any trend should repeat.
# #in our case, because our data is weekly, and we think the seasonality is annual, then
# # the period should be 52 (because there are 52 weeks in 1 year)
# #if instead we had monthly data, then the seasonality would be 12 (12 months per year)
# mrsa_ts = staph_isolates %>%
#   filter(SpeciesName == "Methicillin-Resistant Staphylococcus aureus") %>%
#   mutate(date = floor_date(date, "month")) %>%
#   group_by(date, SpeciesName) %>%
#   summarise(n = n()) 
# 
# msts(mrsa_ts$n, c(12), start = c(2000,2,1)) %>%
#   mstl() %>%
#   autoplot() + theme_bw()
# 
# mrsa_ts = msts(mrsa_ts$n, c(12), start = c(2000,2,1))
# mrsa_ts_decom = mstl(stats::window(mrsa_ts, end = c(2007,1)))
# mrsa_ts = data.frame(data = as.matrix(mrsa_ts), date = as_date(as.yearmon(time(mrsa_ts))))
# 
# arima_global = auto.arima(mrsa_ts_decom[,2]+mrsa_ts_decom[,3])
# arima_trend = auto.arima(mrsa_ts_decom[,2])
# arima_season = auto.arima(mrsa_ts_decom[,3])
# 
# plot(arima_global$fitted)
# lines(arima_global$x, col = "red")
# plot(arima_trend$fitted)
# lines(arima_trend$x, col = "red")
# plot(arima_season$fitted)
# lines(arima_season$x, col = "red")
# 
# forecast_trend = forecast(arima_trend, h = 142)
# forecast_season = forecast(arima_season, h = 142)
# forecast_combine = forecast_trend[["mean"]] + forecast_season[["mean"]]
# forecast_upper = forecast_trend[["upper"]] + forecast_season[["upper"]]
# forecast_lower = forecast_trend[["lower"]] + forecast_season[["lower"]]
# forecast_interval = cbind(forecast_upper, forecast_lower)
# colnames(forecast_interval) = c("upper80", "upper95", "lower80", "lower95")
# 
# forecast_interval = cbind(as.data.frame(forecast_interval),
#                           date = as_date(as.yearmon(time(forecast_interval))))
# forecast_interval$mean = as.matrix(forecast_combine)
# 
# forecast_global = forecast(arima_global, h = 142)
# forecast_interval_global = cbind(as.data.frame(forecast_global),
#                                  date = forecast_interval$date)
# colnames(forecast_interval_global) = c("mean", "lower80","upper80", "lower95", "upper95", "date")
# 
# ggplot() +
#   geom_ribbon(data = forecast_interval_global, aes(ymin = lower95, ymax = upper95,
#                                                    x = date), alpha = 0.3, fill = "blue") +
#   geom_line(data = forecast_interval_global, aes(x = date, y = mean), colour = "blue") +
#   geom_ribbon(data = forecast_interval, aes(ymin = lower95, ymax = upper95,
#                                             x = date), alpha = 0.3, fill = "red") +
#   geom_line(data = forecast_interval, aes(x = date, y = mean), colour = "red") +
#   geom_line(data = mrsa_ts, aes(x = date, y = data)) +
#   theme_bw()
