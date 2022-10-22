#Reading file
#Load file

#Phenotype folder
#LASKY PHENOTYPE_ filtered for africa
lasky_africa <- as_tibble(read.csv("lasky_africa.csv")) %>% relocate(hapmap_id)
head(lasky_africa)

#HARVEST DATA
planting_data <- as_tibble(read.csv("africa_harvest.csv", header = FALSE)) %>% janitor::row_to_names(1,remove_rows_above = FALSE) 
names(planting_data)[1] <- "country"
head(planting_data)


#Replacing planting data with temperature columns in 
plant <- planting_data %>% 
  mutate(Plant = str_split(c(Plant), ",")) %>%
  unnest(Plant) %>%
  mutate(Group = 
           case_when(Plant == 1 ~ "tmax_1",
                     Plant == 2 ~ "tmax_2",
                     Plant == 3 ~ "tmax_3",
                     Plant == 4 ~ "tmax_4",
                     Plant == 5 ~ "tmax_5",
                     Plant == 6 ~ "tmax_6",
                     Plant == 7 ~ "tmax_7",
                     Plant == 8 ~ "tmax_8",
                     Plant == 9 ~ "tmax_9",
                     Plant == 10 ~ "tmax_10",
                     Plant == 11 ~ "tmax_11",
                     Plant == 12 ~ "tmax_12"
           ))

mid_season <- planting_data %>% 
  mutate(`Mid,Season` = str_split(c(`Mid,Season`), ",")) %>%
  unnest(`Mid,Season`) %>%
  mutate(Group = 
           case_when(`Mid,Season` == 1 ~ "tmax_1",
                     `Mid,Season` == 2 ~ "tmax_2",
                     `Mid,Season` == 3 ~ "tmax_3",
                     `Mid,Season` == 4 ~ "tmax_4",
                     `Mid,Season` == 5 ~ "tmax_5",
                     `Mid,Season` == 6 ~ "tmax_6",
                     `Mid,Season` == 7 ~ "tmax_7",
                     `Mid,Season` == 8 ~ "tmax_8",
                     `Mid,Season` == 9 ~ "tmax_9",
                     `Mid,Season` == 10 ~ "tmax_10",
                     `Mid,Season` == 11 ~ "tmax_11",
                     `Mid,Season` == 12 ~ "tmax_12"
           ))

Harvest <- planting_data %>% 
  mutate(Harvest = str_split(c(Harvest), ",")) %>%
  unnest(Harvest) %>%
  mutate(Group = 
           case_when(Harvest == 1 ~ "tmax_1",
                     Harvest == 2 ~ "tmax_2",
                     Harvest == 3 ~ "tmax_3",
                     Harvest == 4 ~ "tmax_4",
                     Harvest == 5 ~ "tmax_5",
                     Harvest == 6 ~ "tmax_6",
                     Harvest == 7 ~ "tmax_7",
                     Harvest == 8 ~ "tmax_8",
                     Harvest == 9 ~ "tmax_9",
                     Harvest == 10 ~ "tmax_10",
                     Harvest == 11 ~ "tmax_11",
                     Harvest == 12 ~ "tmax_12"
           ))



# Combining both data set by country
lasky_africa <- left_join(lasky_africa,planting_data)
head(lasky_africa)







#Separating the planting seasons
lasky_africa[c("planting_start","planting_end")] <- str_split_fixed(lasky_africa$Plant,',',2)
lasky_africa[c("mid_season_start","mid_season_end")] <- str_split_fixed(lasky_africa$`Mid,Season`,',',2)
lasky_africa[c("harvest_start","harvest_end")] <- str_split_fixed(lasky_africa$Harvest,',',2)



