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

others <-c("Afganistan","Bangladesh","China","India","Indonesia","Iraq",
             "Japan","Myanmar","Pakistan","Philippines","Sri Lanka","Syria","Taiwan","Turkey",
             "United States of America","Yemen")

z <- filter(x, country != "Afghanistan" & country != "Bangladesh" & country != "China" & country != "India" & country != "Indonesia"
            & country != "Iraq" & country != "Japan"& country != "Myanmar"& country != "Pakistan" & country != "Philippines"& country != "Sri Lanka"& country != "Syria"
            & country != "Taiwan"& country != "Turkey"& country != "United States of America" & country != "Yemen" & country != "El Salvador" & country != "Saudi Arabia" & country != "Nicaragua" & country != "Barbados" & country != "Dominican Republic")

# Combining countries latitude and temp
head(z)
z$hapmap_id <- rownames(z)
head(z)
head(geo_hap)

c <- inner_join(z, geo_hap, by = "hapmap_id")

geo_temp <- c %>%
  dplyr::select(hapmap_id,country,avgt_min,avgt_max)

write.csv(geo_temp,"geo_temp.csv", row.names = F)

unique(z$country)




# Crop Map
# https://ipad.fas.usda.gov/ogamaps/cropcalendar.aspx


# Combining countries temp and alt
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

others <-c("Afganistan","Bangladesh","China","India","Indonesia","Iraq",
           "Japan","Myanmar","Pakistan","Philippines","Sri Lanka","Syria","Taiwan","Turkey",
           "United States of America","Yemen")

z <- filter(x, country != "Afghanistan" & country != "Bangladesh" & country != "China" & country != "India" & country != "Indonesia"
            & country != "Iraq" & country != "Japan"& country != "Myanmar"& country != "Pakistan" & country != "Philippines"& country != "Sri Lanka"& country != "Syria"
            & country != "Taiwan"& country != "Turkey"& country != "United States of America" & country != "Yemen" & country != "El Salvador" & country != "Saudi Arabia" & country != "Nicaragua" & country != "Barbados" & country != "Dominican Republic")

# Combining countries latitude and temp
head(z)
z$hapmap_id <- rownames(z)
head(z)
head(geo_hap)

c <- inner_join(z, alt.ph, by = "hapmap_id")

geo_alt.ph <- c %>%
  dplyr::select(hapmap_id,country,Alt, topsoil_pH)


temp.alt.ph <- inner_join(geo_temp,geo_alt.ph)

write.csv(temp.alt.ph,"temp.alt.ph.csv", row.names = F)

