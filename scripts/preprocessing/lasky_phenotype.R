#Reading file
#Load file

#Phenotype folder
#LASKY PHENOTYPE_ filtered for africa
lasky_africa <- as_tibble(read.csv("lasky_africa.csv")) %>% relocate(hapmap_id)
head(lasky_africa)



#HARVEST DATA
planting_data <- as_tibble(read.csv("africa_harvest.csv", header = TRUE)) %>% janitor::row_to_names(1,remove_rows_above = FALSE) 
names(planting_data)[1] <- "country"
head(planting_data)

# Combining both dataset by country
lasky_africa <- left_join(lasky_africa,planting_data)
head(lasky_africa)





