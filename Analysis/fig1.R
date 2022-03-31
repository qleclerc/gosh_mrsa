#time-series analysis of s aureus cases
#tldr: no seasonality

library(tidyr)
library(dplyr)
library(lubridate)
library(ggplot2)
library(cowplot)
library(forecast)
library(zoo)

staph_isolates = read.csv(here::here("Clean", "staph_isolates.csv")) %>% 
  mutate(date = as_date(date))

#some general plots
#baseline incidence
p1 = staph_isolates %>%
  mutate(date = floor_date(date, "month")) %>%
  # group_by(project_id, date, SpeciesName) %>%
  # summarise(n = 1) %>%
  # ungroup() %>%
  group_by(date, SpeciesName) %>%
  summarise(n = n()) %>%
  ggplot() +
  geom_line(aes(date, n, colour = SpeciesName)) +
  theme_bw() +
  labs(y = "Number of isolates", x = "Time (years)", colour = "Group:") +
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=12),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12))


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
  geom_line(aes(date, prop, colour = SpeciesName)) +
  theme_bw() +
  labs(y = "Proportion of isolates", x = "Time (years)", colour = "Group:") +
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=12),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12))

plot_grid(plot_grid(p1 + theme(legend.position = "none"),
                    p2 + theme(legend.position = "none"),
                    ncol = 2, labels = c("a)", "b)")),
          get_legend(p1 + theme(legend.position = "bottom")),
          nrow = 2,
          rel_heights = c(1,0.1))

ggsave(here::here("Figures", "fig1.png"))

#descriptive stats
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
(19780+3519-22206)/22206
(19780-1093)/22206
(3514-1093)/22206

#periods in msts: the period defines how often any trend should repeat.
#in our case, because our data is weekly, and we think the seasonality is annual, then
# the period should be 52 (because there are 52 weeks in 1 year)
#if instead we had monthly data, then the seasonality would be 12 (12 months per year)
mrsa_ts = staph_isolates %>%
  filter(SpeciesName == "Methicillin-Resistant Staphylococcus aureus") %>%
  mutate(date = floor_date(date, "month")) %>%
  group_by(project_id, date, SpeciesName) %>%
  summarise(n = 1) %>%
  ungroup() %>%
  group_by(date, SpeciesName) %>%
  summarise(n = sum(n)) 

msts(mrsa_ts$n, c(12), start = c(2000,1,1)) %>%
  mstl() %>%
  autoplot() + theme_bw()

mrsa_ts = msts(mrsa_ts$n, c(12), start = c(2000,1,1))
mrsa_ts_decom = mstl(stats::window(mrsa_ts, end = c(2009,1)))
mrsa_ts = data.frame(data = as.matrix(mrsa_ts), date = as_date(as.yearmon(time(mrsa_ts))))

arima_global = auto.arima(mrsa_ts_decom[,2]+mrsa_ts_decom[,3])
arima_trend = auto.arima(mrsa_ts_decom[,2])
arima_season = auto.arima(mrsa_ts_decom[,3])

plot(arima_global$fitted)
lines(arima_global$x, col = "red")
plot(arima_trend$fitted)
lines(arima_trend$x, col = "red")
plot(arima_season$fitted)
lines(arima_season$x, col = "red")

forecast_trend = forecast(arima_trend, h = 48)
forecast_season = forecast(arima_season, h = 48)
forecast_combine = forecast_trend[["mean"]] + forecast_season[["mean"]]
forecast_upper = forecast_trend[["upper"]] + forecast_season[["upper"]]
forecast_lower = forecast_trend[["lower"]] + forecast_season[["lower"]]
forecast_interval = cbind(forecast_upper, forecast_lower)
colnames(forecast_interval) = c("upper80", "upper95", "lower80", "lower95")

forecast_interval = cbind(as.data.frame(forecast_interval),
                          date = as_date(as.yearmon(time(forecast_interval))))
forecast_interval$mean = as.matrix(forecast_combine)

forecast_global = forecast(arima_global, h = 48)
forecast_interval_global = cbind(as.data.frame(forecast_global),
                                 date = forecast_interval$date)
colnames(forecast_interval_global) = c("mean", "lower80","upper80", "lower95", "upper95", "date")

ggplot() +
  geom_ribbon(data = forecast_interval_global, aes(ymin = lower95, ymax = upper95,
                                                   x = date), alpha = 0.3, fill = "blue") +
  geom_line(data = forecast_interval_global, aes(x = date, y = mean), colour = "blue") +
  geom_ribbon(data = forecast_interval, aes(ymin = lower95, ymax = upper95,
                                            x = date), alpha = 0.3, fill = "red") +
  geom_line(data = forecast_interval, aes(x = date, y = mean), colour = "red") +
  geom_line(data = mrsa_ts, aes(x = date, y = data)) +
  theme_bw()
