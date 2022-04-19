library(ppcor)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(lubridate)
library(reshape2)
library(FSA)


staph_isolates = read.csv(here::here("Clean", "staph_isolates.csv")) %>%
  mutate(date = as_date(date))

#replace S by 1, R by 2. This way, a mean value of 2 between two samples for a given abx implies a change in resistance
#that approach allows me to ignore NA values in the comparison
staph_isolates_sr = staph_isolates[,-c(1:5)]
staph_isolates_sr[staph_isolates_sr == "S"] = 1
staph_isolates_sr[staph_isolates_sr == "R"] = 3

staph_isolates_sr = apply(staph_isolates_sr, c(1,2), as.numeric)

staph_isolates[,-c(1:5)] = staph_isolates_sr

#assumption: only keep antibiotics with more than 50 results (S or I)
good_cols = which(colSums(staph_isolates_sr, na.rm = T) > 50)
good_rows = which(rowSums(staph_isolates_sr, na.rm = T) > 0)

staph_isolates_profiles = staph_isolates %>%
  select(c(1:5,(good_cols+5)))

staph_isolates_profiles = staph_isolates_profiles %>%
  filter(row_number() %in% good_rows)

#only keep patients with more than one sample
good_ids = staph_isolates_profiles %>%
  group_by(project_id) %>%
  summarise(n = n()) %>%
  filter(n > 1) %>%
  select(project_id) %>%
  pull

staph_isolates_profiles = staph_isolates_profiles %>%
  filter(project_id %in% good_ids)

diffs_data = data.frame(project_id = staph_isolates_profiles$project_id[1],
                        date = staph_isolates_profiles$date[1],
                        species = "staph",
                        samples = 0,
                        #diversity_score = 0,
                        #diffs = 0,
                        n_profiles = 0)
abx_diffs = c()
#heavy lifting here
#for each row, compare the resistance profiles with all other profiles from same patient, day and species
for(i in unique(staph_isolates_profiles$project_id)){
  
  compare_i = staph_isolates_profiles %>%
    filter(project_id == i)
  
  for(j in unique(compare_i$date)){
    
    compare_j = compare_i %>%
      filter(date == j)
    
    for(k in unique(compare_j$SpeciesName)){
      
      compare_k = compare_j %>%
        filter(SpeciesName == k) %>%
        select(-c(1:5)) 
      
      #skip if just one sample
      if(nrow(compare_k) == 1) next
      
      #remove duplicates
      compare_k_c = distinct(compare_k)
      
      #initially assume all isolates are unique
      profiles_k = nrow(compare_k_c)
      
      if(profiles_k > 1){
        
        #go through each isolate
        for(l in 1:(nrow(compare_k_c)-1)){
          
          compare_l = compare_k_c[l,]
          
          #compare with all subsequent isolates
          for(m in (l+1):nrow(compare_k_c)){
            
            #if there's a difference, flag isolate l as non unique (reducing number of resistance
            #profiles by 1), and move on to next isolate
            if(!(any(apply(rbind(compare_l, compare_k_c[m,]),
                           2, mean, na.rm = T) == 2, na.rm = T))){
              
              profiles_k = profiles_k-1
              break
              
            }
          }
        }
      }
      
      # #compare all profiles
      # #a mean between 1 and 3 will indicate a difference
      # #na.rm means we're ignoring NAs instead of assuming they're resistant or susceptible
      # tt = apply(compare_k, 2, mean, na.rm = T)
      # #standardise
      # tt[is.nan(tt)] = 1
      # tt[tt <= 2] = tt[tt <= 2] - 1
      # tt[tt > 2] = 3 - tt[tt > 2]
      # #this way, diffs values are between 0 and 1
      # #a diff of 1 means that there is maximum diversity (half of samples are S/R, other half are R/S)
      # #a diff of 0 means no difference across samples
      # 
      # abx_diffs = c(abx_diffs, names(tt[tt!=0]))
      # 
      # # #diversity assuming NA = R
      # # div_r_assum = compare_k %>%
      # #   mutate(across(everything(), ~replace(., is.na(.), 3))) %>%
      # #   distinct %>%
      # #   nrow
      # # 
      # # #diversity assuming NA = S
      # # div_s_assum = compare_k %>%
      # #   mutate(across(everything(), ~replace(., is.na(.), 1))) %>%
      # #   distinct %>%
      # #   nrow
      # 
      
      diffs_data_k = data.frame(project_id = i,
                                date = as_date(j),
                                species = k,
                                samples = nrow(compare_k),
                                #diversity_score = mean(tt),
                                #diffs = sum(tt!=0),
                                n_profiles = profiles_k)
      # upper_div = max(div_r_assum, div_s_assum),
      # lower_div = min(div_r_assum, div_s_assum))
      
      diffs_data = rbind(diffs_data, diffs_data_k)
      
    }
  }
}

diffs_data = diffs_data[-1,]

write.csv(diffs_data, here::here("Clean", "isolates_diversity.csv"), row.names = F)

table(diffs_data$n_profiles)/nrow(diffs_data)

ranked_diffs = diffs_data %>%
  arrange(-n_profiles) %>%
  mutate(cumn = cumsum(n_profiles)) %>%
  mutate(cumsamples = cumsum(samples)) %>%
  mutate(cumn = cumn/max(cumn)) %>%
  mutate(cumsamples = cumsamples/max(cumsamples))

min(which(ranked_diffs$cumn > 0.8))/nrow(ranked_diffs)
#80% of diversity in 10% of patients with more than 1 sample taken
ranked_diffs$cumsamples[min(which(ranked_diffs$cumn > 0.8))]
#These 10% of patients represent 30% of all samples

ggplot(diffs_data) +
  geom_jitter(aes(x = samples, y = n_profiles, colour = species), alpha = 0.3) +
  geom_abline(intercept = 0, slope = 1) +
  facet_wrap(~species)+
  scale_y_continuous(limits = c(0,6), breaks  = seq(1,6,2)) +
  scale_x_continuous(limits = c(0,14), breaks  = seq(1,14,2)) +
  theme_bw() +
  theme(legend.position = "bottom")

diffs_data %>%
  mutate(date = floor_date(date, "year")) %>%
  ggplot() +
  geom_boxplot(aes(y = n_profiles, group = species, colour = species)) +
  theme_bw()


boxplot(diffs_data$n_profiles~diffs_data$samples)

kruskal.test(n_profiles~samples, data = diffs_data)
dunnTest(n_profiles ~ as.factor(samples), data = diffs_data, method = "bonferroni")


ggplot(ranked_diffs) +
  geom_line(aes(x = 1:nrow(ranked_diffs), cumsamples, colour = "Cumulative samples taken")) +
  geom_line(aes(x = 1:nrow(ranked_diffs), cumn, colour = "Cumulative isolates")) +
  theme_bw() +
  labs(x = "Patients", y = "Proportion", colour = "Variable:")

ggsave(here::here("Figures", "diversity_fig.png"))
