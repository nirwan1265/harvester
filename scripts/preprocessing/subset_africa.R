# Setting up table
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Research/Data/Phenotype")
long_lat <- read.table("long_lat.txt", sep = "")
#Change first column to row names or read the file from map_id

x <- long_lat %>%
  column_to_rownames("hapmap_id") %>%
  dplyr::select(Longitude, Latitude)
  


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

#Filtering only Africa

accession_africa_filtered <- geo_temp[,1]
write.table(accession_africa_filtered, "accession_africa_filtered.txt", quote = F, row.names = F, col.names = F, eol = "\t")





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
           "Japan","Myanmar"x,"Pakistan","Philippines","Sri Lanka","Syria","Taiwan","Turkey",
           "United States of America","Yemen")

z <- filter(x, country != "Afghanistan" & country != "Bangladesh" & country != "China" & country != "India" & country != "Indonesia"
            & country != "Iraq" & country != "Japan"& country != "Myanmar"& country != "Pakistan" & country != "Philippines"& country != "Sri Lanka"& country != "Syria"
            & country != "Taiwan"& country != "Turkey"& country != "United States of America" & country != "Yemen" & country != "El Salvador" & country != "Saudi Arabia" & country != "Nicaragua" & country != "Barbados" & country != "Dominican Republic")

head(z)
z$hapmap_id <- rownames(z)
head(z)
head(geo_hap)

c <- inner_join(z, alt.ph, by = "hapmap_id")

geo_alt.ph <- c %>%
  dplyr::select(hapmap_id,country,Alt, topsoil_pH)


temp.alt.ph <- inner_join(geo_temp,geo_alt.ph)

write.csv(temp.alt.ph,"temp.alt.ph.csv", row.names = F)



#########################################################################################################
## LASKY phenotypes
#########################################################################################################


# Combining countries latitude and temp
head(z)
z$hapmap_id <- rownames(z)
head(z)
head(lasky)

c <- inner_join(z, geo_hap, by = "hapmap_id")



write.csv(c,"lasky_africa.csv", row.names = F)




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
           "Japan","Myanmar"x,"Pakistan","Philippines","Sri Lanka","Syria","Taiwan","Turkey",
           "United States of America","Yemen")

z <- filter(x, country != "Afghanistan" & country != "Bangladesh" & country != "China" & country != "India" & country != "Indonesia"
            & country != "Iraq" & country != "Japan"& country != "Myanmar"& country != "Pakistan" & country != "Philippines"& country != "Sri Lanka"& country != "Syria"
            & country != "Taiwan"& country != "Turkey"& country != "United States of America" & country != "Yemen" & country != "El Salvador" & country != "Saudi Arabia" & country != "Nicaragua" & country != "Barbados" & country != "Dominican Republic")

head(z)
z$hapmap_id <- rownames(z)
head(z)
head(geo_hap)

c <- inner_join(z, alt.ph, by = "hapmap_id")

geo_alt.ph <- c %>%
  dplyr::select(hapmap_id,country,Alt, topsoil_pH)


temp.alt.ph <- inner_join(geo_temp,geo_alt.ph)

write.csv(temp.alt.ph,"temp.alt.ph.csv", row.names = F)
(hapmap_id, Latitude, Longitude, Alt, AridityIndex, prec_1,prec_2,prec_3,prec_4,
  prec_5,prec_6,prec_7,prec_8,prec_9,prec_10,prec_11,prec_12, tmax_1,tmax_2,tmax_3,
  tmax_4,tmax_5,tmax_6,tmax_7,tmax_8,tmax_9,tmax_10,tmax_11,tmax_12
  