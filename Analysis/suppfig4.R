#time-series analysis of s aureus cases
#tldr: no seasonality

library(tidyr)
library(dplyr)
library(lubridate)
library(ggplot2)
library(RColorBrewer)
library(ggtext)
library(scales)
library(cowplot)
library(forecast)
library(zoo)

mrsa_col = "#EFC000FF"
mssa_col = "#0073C2FF"

cpalette = c(RColorBrewer::brewer.pal(12, "Paired"), "black")

staph_isolates = read.csv(here::here("Clean", "staph_isolates.csv")) %>% 
  mutate(date = as_date(date))

typing_data = read.csv(here::here("Clean", "typing_data.csv")) %>%
  mutate(across(matches("datetime"), as_date))

typing_data %>%
  mutate(CC = replace(CC, grepl("ST", CC), "Other")) %>%
  filter(CC != "None") %>%
  select(project_id, start_datetime, CC) %>%
  distinct() %>%
  mutate(start_datetime = floor_date(start_datetime, "year")) %>%
  filter(start_datetime > as.Date("2010-01-01")) %>%
  group_by(start_datetime, CC) %>%
  summarise(n = n()) %>%
  ggplot() +
  geom_col(aes(start_datetime, n, fill = CC), position = "fill") +
  theme_bw() +
  labs(y = "Proportion of typed isolates", x = "Year") +
  scale_fill_manual(values = cpalette) +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12)) +
  scale_x_date(breaks = as.Date(c("2011-01-01", "2013-01-01", "2015-01-01",
                                  "2017-01-01", "2019-01-01", "2021-01-01")), date_labels = "%Y",
               limits = as.Date(c("2010-07-01", "2021-07-01")))

ggsave(here::here("figures", "new_suppfig.png"))
