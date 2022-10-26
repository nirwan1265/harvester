# Setting up table
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Research/Data/Phenotype")


#Reading file
#Load file

#Phenotype folder
#LASKY PHENOTYPE_ filtered for africa
lasky_africa <- as_tibble(read.csv("lasky_africa.csv")) %>% relocate(Taxa) %>% rename(Taxa = "hapmap_id")
head(lasky_africa)

#HARVEST DATA
planting_data <- as_tibble(read.csv("africa_harvest.csv", header = FALSE)) %>% janitor::row_to_names(1,remove_rows_above = FALSE) 
names(planting_data)[1] <- "country"
head(planting_data)


# Combining both data set by country
lasky_africa <- left_join(lasky_africa,planting_data)
head(lasky_africa)



#Planting data

#Replacing planting data with temperature columns in 
plant <- lasky_africa %>% 
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

# Subsetting max temp
plant <- as.data.frame(plant[,c(1,21:32,52)])
plant[is.na(plant)] <- 0
plant$maxtemp <- 0

#Replacing values
for (i in 2:13){
  for (j in 1:nrow(plant)){
    if (plant[j,14] == names(plant)[i]) {
      plant[j,15] <- plant[j,i]
    } else {
      next
    }
  }
}

#plant <- plant[,c(1,15)]
plant <- plant %>%
  dplyr::select("hapmap_id", "maxtemp") %>%
  dplyr::group_by(hapmap_id) %>%
  dplyr::mutate(rid = row_number()) %>%
  pivot_wider(
    id_cols = hapmap_id,
    names_from = rid,
    values_from = maxtemp
  ) %>% 
  `row.names<-`(., NULL) %>%
  tibble::column_to_rownames(var = "hapmap_id") %>% 
  dplyr::mutate(mean_maxtemp = rowMeans(., na.rm = TRUE)) %>%
  dplyr::select(mean_maxtemp)


plant_missing <-  plant %>% 
  filter(mean_maxtemp == 0)
  
 
  



#Mid Season data

mid_season <- lasky_africa %>% 
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



# Subsetting max temp
mid_season <- as.data.frame(mid_season[,c(1,21:32,52)])
mid_season[is.na(mid_season)] <- 0
mid_season$maxtemp <- 0

#Replacing values
for (i in 2:13){
  for (j in 1:nrow(mid_season)){
    if (mid_season[j,14] == names(mid_season)[i]) {
      mid_season[j,15] <- mid_season[j,i]
    } else {
      next
    }
  }
}

#mid_season <- mid_season[,c(1,15)]
mid_season <- mid_season %>%
  select(hapmap_id, maxtemp) %>%
  group_by(hapmap_id) %>%
  mutate(rid = row_number()) %>%
  pivot_wider(
    id_cols = hapmap_id,
    names_from = rid,
    values_from = maxtemp
  ) %>% `row.names<-`(., NULL) %>%
  column_to_rownames(var = "hapmap_id") %>% 
  mutate(mean_maxtemp = rowMeans(., na.rm = TRUE)) %>%
  select(mean_maxtemp)


mid_season_missing <-  mid_season %>% 
  filter(mean_maxtemp == 0)




harvest <- lasky_africa %>% 
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

#harvesting data

#Replacing harvesting data with temperature columns in 
harvest <- lasky_africa %>% 
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

# Subsetting max temp
harvest <- as.data.frame(harvest[,c(1,21:32,52)])
harvest[is.na(harvest)] <- 0
harvest$maxtemp <- 0

#Replacing values
for (i in 2:13){
  for (j in 1:nrow(harvest)){
    if (harvest[j,14] == names(harvest)[i]) {
      harvest[j,15] <- harvest[j,i]
    } else {
      next
    }
  }
}

#harvest <- harvest[,c(1,15)]
harvest <- harvest %>%
  select(hapmap_id, maxtemp) %>%
  group_by(hapmap_id) %>%
  mutate(rid = row_number()) %>%
  pivot_wider(
    id_cols = hapmap_id,
    names_from = rid,
    values_from = maxtemp
  ) %>% `row.names<-`(., NULL) %>%
  column_to_rownames(var = "hapmap_id") %>% 
  mutate(mean_maxtemp = rowMeans(., na.rm = TRUE)) %>%
  select(mean_maxtemp)



harvest_missing <-  harvest %>% 
  filter(mean_maxtemp == 0)



####################################################################################################
# Finding empty rows 
####################################################################################################

#LASKY PHENOTYPE filtered for africa
lasky_africa <- as_tibble(read.csv("lasky_africa.csv")) %>% relocate(Taxa) %>% rename(Taxa = "hapmap_id")
head(lasky_africa)

plant_missing <- plant_missing %>%
  rownames_to_column() %>%
  rename("rowname" = "hapmap_id")
head(plant_missing)  


missing_data <- left_join(lasky_africa, plant_missing)

missing_data <- left_join(missing_data, planting_data, by = "country")

unique(missing_data$hapmap_id)



plant <- lasky_africa %>% 
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

# Subsetting max temp
plant <- as.data.frame(plant[,c(1,21:32,52)])
plant[is.na(plant)] <- 0
plant$maxtemp <- 0
