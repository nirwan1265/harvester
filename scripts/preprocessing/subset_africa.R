# Setting up table
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Research/Data/Phenotype")
long_lat <- read.table("long_lat.txt", sep = "")

x <- long_lat %>%
  column_to_rownames("hapmap_id") %>%
  dplyr::select(Longitude, Latitude)
  #mutate(country = y)
  #coords2country()


y <- coords2country(x)

x <- x %>%
  mutate(country = y)

x <- x %>%
  !filter(x, country == c("Afganistan","Bangladesh","China","India","Indonesia","Iraq",
                          "Japan","Myanmar","Pakistan","Philippines","Sri Lanka","Syria","Taiwan","Turkey",
                          "United States of America","Yemen"))
  
  



