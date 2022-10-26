# Working directory
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Research/Data/Phenotype")

# country_code.csv is produced from the python script

# Country code:
coco <- read.csv("country_code.csv")
coco <- toupper(coco[,1])
coco

# Change to country name
countryname <- as.data.frame(countrycode(coco, 'iso2c', "country.name"))


# co-ordinates:
long_lat <- read.table("long_lat.txt", sep = "") %>% `row.names<-`(., NULL) %>%
  column_to_rownames(var = "hapmap_id")
long_lat


# Combining longitude and latitude data with country code
long_lat_country <- cbind(long_lat, countryname)
names(long_lat_country)[3] <- "country"

# African countries

africa <- c("Algeria","Angola","Benin","Botswana","Burkina Faso","Burundi",
            "Cabo Verde","Cameroon","Central African Republic","Chad",
            "Comoros","Democratic Republic of the Congo","Congo - Kinshasa",
            "Republic of the Congo","Cote d'Ivoire","Djibouti","Egypt",
            "Equatorial Guinea","Eritrea","Eswatini","Ethiopia","Gabon","Gambia"
            ,"Ghana","Guinea","Guinea-Bissau","Kenya","Lesotho","Liberia",
            "Libya","Madagascar","Malawi","Mali","Mauritania","Mauritius",
            "Morocco","Mozambique","Namibia","Niger","Nigeria","Rwanda",
            "Sao Tome and Principe","Senegal","Seychelles","Sierra Leone",
            "Somalia","South Africa","South Sudan","Sudan","Tanzania","Togo",
            "Tunisia","Uganda","Zambia","Zimbabwe")


acc_africa <- long_lat_country[long_lat_country$country %in% africa, ] %>% tibble::rownames_to_column("Taxa") %>% relocate(Taxa)


#Phosphorus data
phospho <- as.data.frame(read.csv("taxa_geoloc_pheno.csv"))

#Subsetting only african accesions
acc_pheno_africa <- left_join(acc_africa,phospho) 
acc_pheno_africa <- acc_pheno_africa[,(-2:-10)]

#Saving
write.csv(acc_pheno_africa, "Sorghum_allphospho_africa.csv", row.names = FALSE)



# Subsetting maize
head(phospho)
mays_phospho <- phospho[which(phospho$sp == "Zea mays"), ]
mays_phospho <- mays_phospho[,c(-1,-3:-7)]


#Saving
write.csv(mays_phospho, "Maize_allphospho.csv", row.names = FALSE)

