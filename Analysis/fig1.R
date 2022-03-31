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

#is the peak because of under 5 patients?
admissions = read.csv(here::here("Data", "combined_patient_hospital_admissions.csv")) %>%
  mutate(start_datetime = as_date(start_datetime)) %>%
  select(project_id, start_datetime)
demo = read.csv(here::here("Data", "caboodle_patient_demographics.csv")) %>%
  mutate(birth_date = as_date(birth_date)) %>%
  select(project_id, birth_date)

admissions = left_join(admissions, demo, "project_id") %>%
  mutate(admin_age = floor(interval(birth_date, start_datetime)/years(1))) %>%
  filter(start_datetime > "2000-02-01") %>%
  filter(!is.na(admin_age)) %>%
  filter(admin_age >= 0) #some patients were admitted before born, not sure if that's preterm?

admissions %>%
  mutate(young = (admin_age < 5)) %>%
  mutate(date = floor_date(start_datetime, "month")) %>%
  group_by(date) %>%
  #summarise(admin_age = median(admin_age)) %>%
  summarise(young = sum(young)) %>%
  ggplot() +
  geom_line(aes(date, young)) +
  theme_bw()


#is the peak because of bone marrow transplant surgery?
theatre1 = read.csv(here::here("Data", "caboodle_patient_theatre_list.csv")) %>%
  mutate(date = as_date(start_datetime)) %>%
  filter(grepl("MARROW", procedure_name)) %>%
  select(date)
theatre2 = read.csv(here::here("Data", "datalake_patient_theatre_list.csv")) %>%
  mutate(date = as_date(start_datetime)) %>%
  filter(grepl("arrow", procedure_name)) %>%
  select(date)

rbind(theatre1, theatre2) %>%
  mutate(date = floor_date(date, "year")) %>%
  group_by(date) %>%
  summarise(n = n()) %>%
  ggplot() +
  geom_line(aes(date, n)) +
  theme_bw()


#is the peak because of any surgery?

theatre1 = read.csv(here::here("Data", "caboodle_patient_theatre_list.csv")) %>%
  mutate(date = as_date(start_datetime)) %>%
  select(date)
theatre2 = read.csv(here::here("Data", "datalake_patient_theatre_list.csv")) %>%
  mutate(date = as_date(start_datetime)) %>%
  select(date)

rbind(theatre1, theatre2) %>%
  mutate(date = floor_date(date, "month")) %>%
  group_by(date) %>%
  summarise(n = n()) %>%
  ggplot() +
  geom_line(aes(date, n)) +
  theme_bw()


#is the peak because of ICU?
icu_patients = read.csv(here::here("Data", "combined_patient_ward_stays.csv")) %>%
  filter(icu_ward_stay == 1) %>%
  mutate(date = as_date(start_datetime))

icu_patients %>%
  mutate(date = floor_date(date, "month")) %>%
  group_by(date) %>%
  summarise(n = n()) %>%
  ggplot()+
  geom_line(aes(date, n)) +
  theme_bw()


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
  group_by(date, SpeciesName) %>%
  summarise(n = n()) 

msts(mrsa_ts$n, c(12), start = c(2000,2,1)) %>%
  mstl() %>%
  autoplot() + theme_bw()

mrsa_ts = msts(mrsa_ts$n, c(12), start = c(2000,2,1))
mrsa_ts_decom = mstl(stats::window(mrsa_ts, end = c(2007,1)))
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

forecast_trend = forecast(arima_trend, h = 142)
forecast_season = forecast(arima_season, h = 142)
forecast_combine = forecast_trend[["mean"]] + forecast_season[["mean"]]
forecast_upper = forecast_trend[["upper"]] + forecast_season[["upper"]]
forecast_lower = forecast_trend[["lower"]] + forecast_season[["lower"]]
forecast_interval = cbind(forecast_upper, forecast_lower)
colnames(forecast_interval) = c("upper80", "upper95", "lower80", "lower95")

forecast_interval = cbind(as.data.frame(forecast_interval),
                          date = as_date(as.yearmon(time(forecast_interval))))
forecast_interval$mean = as.matrix(forecast_combine)

forecast_global = forecast(arima_global, h = 142)
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
