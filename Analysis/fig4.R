
library(dplyr)
library(lubridate)
library(reshape2)
library(ggplot2)
library(cowplot)
library(Hmisc)
library(corrplot)

staph_isolates = read.csv(here::here("Clean", "staph_isolates.csv")) %>%
  mutate(date = as_date(date))

staph_resistances = staph_isolates %>%
  select(5:59) %>%
  melt(.,id.vars = "date") %>%
  filter(!is.na(value)) %>%
  filter(value %in% c("S", "R")) %>%
  mutate(date = floor_date(date, "month")) %>%
  arrange(date) %>%
  group_by(date, variable, value) %>%
  summarise(n = n()) %>%
  mutate(prop = n/sum(n)) %>%
  mutate(variable = as.character(variable))


resistances = unique(staph_resistances$variable)

for(i in c(1,10,19,28,37,44)){
  
  pp = staph_resistances %>%
    filter(value == "R") %>%
    filter(variable %in% resistances[(i):(i+8)]) %>%
    ggplot() +
    geom_line(aes(date, n, colour = variable)) +
    theme_bw()
  
  plot(pp)
  
}

good_resistances = c("Amikacin", "Amik.Fluclox", "Chloramphenicol", "Erythromycin",
                     "Ciprofloxacin", "Flucloxacillin", "Fucidin", "Gentamicin",
                     "Gent.Cipro", "Linezolid", "Mupirocin",
                     "Neomycin", "Nitrofurantoin", "Penicillin", "Rifampicin",
                     "Teicoplanin", "Tetracycline", "Trimethoprim", "Vancomycin")

staph_resistances = staph_resistances %>%
  filter(variable %in% good_resistances)


pa = staph_resistances %>%
  filter(value == "R") %>%
  filter(variable %in% c("Amikacin", "Amik.Fluclox", "Flucloxacillin", "Ciprofloxacin")) %>%
  ggplot() +
  geom_line(aes(date, prop, colour = variable)) +
  theme_bw() +
  scale_y_continuous(limits = c(0,1)) +
  labs(x = "Time (years)", y = "Proportion of resistant isolates", colour = "Antibiotic:")

pb = staph_resistances %>%
  filter(value == "R") %>%
  filter(variable %in% c("Gentamicin", "Gent.Cipro", "Rifampicin", "Chloramphenicol")) %>%
  ggplot() +
  geom_line(aes(date, prop, colour = variable)) +
  theme_bw() +
  scale_y_continuous(limits = c(0,1)) +
  labs(x = "Time (years)", y = "Proportion of resistant isolates", colour = "Antibiotic:")

pc = staph_resistances %>%
  filter(value == "R") %>%
  filter(variable %in% c("Fucidin", "Erythromycin", "Trimethoprim")) %>%
  ggplot() +
  geom_line(aes(date, prop, colour = variable)) +
  theme_bw() +
  scale_y_continuous(limits = c(0,1)) +
  labs(x = "Time (years)", y = "Proportion of resistant isolates", colour = "Antibiotic:")

pd = staph_resistances %>%
  filter(value == "R") %>%
  filter(variable %in% c("Mupirocin", "Penicillin", "Tetracycline")) %>%
  ggplot() +
  geom_line(aes(date, prop, colour = variable)) +
  theme_bw() +
  scale_y_continuous(limits = c(0,1)) +
  labs(x = "Time (years)", y = "Proportion of resistant isolates", colour = "Antibiotic:")

plot_grid(pa, pb, pc, pd,
          ncol = 2, labels = c("a)", "b)", "c)", "d)"))

ggsave(here::here("Figures", "fig4.png"))

#supp
staph_resistances %>%
  filter(value == "R") %>%
  filter(variable %in% c("Teicoplanin", "Vancomycin","Linezolid")) %>%
  ggplot() +
  geom_line(aes(date, prop, colour = variable)) +
  theme_bw() +
  scale_y_continuous(limits = c(0,0.1)) +
  labs(x = "Time (years)", y = "Proportion of resistant isolates", colour = "Antibiotic:")

ggsave(here::here("Figures", "suppfig3.png"))


#correlations
cor_res = staph_resistances %>%
  filter(value == "R") %>%
  dcast(., date~variable, value.var = "prop", fill = 0) %>%
  select(-"date")

cor_res = rcorr(as.matrix(cor_res))
cor_res$P[is.na(cor_res$P)] = 0.04
cor_res$P[which(cor_res$P < 0.05)] = 1
cor_res$P[which(cor_res$P != 1)] = 0

cor_res$r = cor_res$r * cor_res$P

png(here::here("Figures", "suppfig4.png"), width = 1100, height = 701)
corrplot(cor_res$r,
         method = "number", type = "upper", order = "hclust", tl.col = 'black', cl.cex = 1)
dev.off()

#investigate high risk of ARG transduction in phage therapy 
#in of risk phage therapy
#investigate high ARG transduction

## antibiotic classes ####

# abx_classes = data.frame(variable = c("Amik.Fluclox", "Amikacin", "Ampicillin", 
#                                       "Augmentin", "Cefotaxime", "Ceftazidime","Cephradine",
#                                       "Chloramphenicol", "Ciprofloxacin", "Colistin",
#                                       "Erythromycin", "Flucloxacillin", 
#                                       "Fucidin", "Gent.Ceftaz", "Gent.Cipro", "Gent.Pip.Taz",
#                                       "Gentamicin", "Mupirocin", "Naladixic.Acid",
#                                       "Neomycin", "Nitrofurantoin", "Penicillin",
#                                       "Piperacillin...Tazobactam", "Pip.Taz.Amik",
#                                       "Pip.Taz.Cipro", "Rifampicin", "Teicoplanin",
#                                       "Tetracycline", "Trimethoprim", "Vancomycin",
#                                       "Sulphonamide", "Oxacillin", "Pristinomycin",
#                                       "Tobramycin", "Mecillinam", "Syncercid",
#                                       "Ceftriaxone", "Septrin", "Cefuroxime",
#                                       "Linezolid", "Clindamycin", "Timentin",
#                                       "Imipenem", "Meropenem", "Tigecycline",
#                                       "Daptomycin", "Doxycycline", "Minocycline",
#                                       "Fosfomycin", "Metronidazole", "Cefoxitin",
#                                       "Co.Trimoxazole..Septrin."),
#                          class = c("Mixed", 
#                                    "Aminoglycosides", "Penicillins", "Penicillins",
#                                    "Cephalosporins", "Cephalosporins", "Cephalosporins",
#                                    "Chloramphenicol", "Fluoroquinolones", "Polypetides",
#                                    "Macrolides", "Penicillins", "Fucidin", "Mixed", 
#                                    "Mixed", "Mixed", "Aminoglycosides", "Mupirocin",
#                                    "Fluoroquinolones", "Aminoglycosides",
#                                    "Nitrofurans", "Penicillins", "Mixed", "Mixed", 
#                                    "Mixed", "Antimycobacterials", "Glycopeptides",
#                                    "Tetracyclines", "Trimethoprim", "Glycopeptides",
#                                    "Sulphonamides", "Penicillins", "Pristinomycin", 
#                                    "Aminoglycosides", "Penicillins", "Syncercid",
#                                    "Cephalosporins", "Sulfonamides", "Cephalosporins",
#                                    "Oxazolidinones", "Lincosamides", "Mixed",
#                                    "Carbapenems", "Carbapenems", "Tigecycline",
#                                    "Lipopetides", "Tetracyclines", "Tetracyclines",
#                                    "Fosfomycin", "Metronidazole", "Cephalosporins", 
#                                    "Sulfonamides"))
# 
# staph_resistances = left_join(staph_resistances, abx_classes, by = "variable")
# 
# classes = unique(staph_resistances$class)
# 
# for(i in c(1,4,7,10,12)){
#   
#   pp = staph_resistances %>%
#     group_by(date, class, value) %>%
#     summarise(n = sum(n)) %>%
#     mutate(prop = n/sum(n)) %>%
#     filter(value == "R") %>%
#     filter(class %in% classes[i:(i+2)]) %>%
#     ggplot() +
#     geom_line(aes(date, prop, colour = class)) +
#     theme_bw() +
#     scale_y_continuous(limits = c(0,1))
#   
#   plot(pp)
#   
# }
# 
