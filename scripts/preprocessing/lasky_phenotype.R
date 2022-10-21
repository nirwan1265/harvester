#Reading file
#Load file

#Phenotype folder
#LASKY PHENOTYPE
"lasky_africa.csv"
lasky_africa <- c
head(lasky_africa)



#HARVEST DATA
planting_data <- read.csv("africa_harvest.csv", header = TRUE)
names(planting_data) <- planting_data[1,]
planting_data <- planting_data[-1,]

#planting_data <- data.frame(lapply(planting_data, function(x){
  gsub("1","")
}))


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





# Getting country and region information
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Research/Data/Phenotype")
long_lat <- read.table("long_lat.txt", sep = "")
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

head(z)


# Joining Country/region with Planting Data

africa_planting_pheno <- inner_join(z, alt.ph, by = "hapmap_id")

geo_alt.ph <- c %>%
  dplyr::select(hapmap_id,country,Alt, topsoil_pH)


temp.alt.ph <- inner_join(geo_temp,geo_alt.ph)

write.csv(temp.alt.ph,"temp.alt.ph.csv", row.names = F)



