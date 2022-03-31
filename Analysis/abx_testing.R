
library(dplyr)
library(lubridate)
library(reshape2)
library(ggplot2)
library(Hmisc)
library(corrplot)

staph_isolates = read.csv(here::here("Clean", "staph_isolates.csv")) %>%
  mutate(date = as_date(date))

staph_isolates_m = staph_isolates %>%
  #filter(SpeciesName == "Methicillin-Resistant Staphylococcus aureus") %>%
  #filter(SpeciesName == "Methicillin-Susceptible Staphylococcus aureus") %>%
  select(-c(1:4)) %>%
  melt(id.vars = "date") %>%
  filter(!is.na(value)) %>%
  mutate(date = floor_date(date, "month")) %>%
  group_by(date, variable) %>%
  summarise(n = n()) %>%
  filter(n > 10) %>%
  filter(variable %in%c("Amikacin", "Amik.Fluclox", "Chloramphenicol", "Erythromycin",
                        "Ciprofloxacin", "Flucloxacillin", "Fucidin", "Gentamicin",
                        "Gent.Cipro", "Linezolid", "Mupirocin","Penicillin", "Rifampicin",
                        "Teicoplanin", "Tetracycline", "Trimethoprim", "Vancomycin"))

isolates_per_day = staph_isolates %>%
  mutate(date = floor_date(date, "month")) %>%
  group_by(date) %>%
  summarise(total = n())

staph_isolates_m = left_join(staph_isolates_m, isolates_per_day, by = "date")

ggplot(staph_isolates_m ) +
  geom_line(aes(date, n/total, colour = variable)) +
  theme_bw()


cor_res = staph_isolates_m %>%
  dcast(., date~variable, value.var = "n", fill = 0) %>%
  select(-"date")

cor_res = rcorr(as.matrix(cor_res))
cor_res$P[is.na(cor_res$P)] = 0.04
cor_res$P[which(cor_res$P < 0.05)] = 1
cor_res$P[which(cor_res$P != 1)] = 0

cor_res$r = cor_res$r * cor_res$P

png(here::here("Figures", "suppfig2.png"), width = 1100, height = 701)
corrplot(cor_res$r,
         method = "number", type = "upper", order = "hclust", tl.col = 'black', cl.cex = 1)
dev.off()

