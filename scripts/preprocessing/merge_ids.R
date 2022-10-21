setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Research/Data/Phenotype")

georef <- read.csv(file = "temp_minmax.csv", na.strings = "NA")
georef

hapmap <- read.table("genotype_ids.txt", header = F)
colnames(hapmap)[1] = "hapmap_id"
hapmap

colnames(georef)
by_IS <- georef %>%
  inner_join(hapmap, by = c(is_no = "V2")) %>% print()

by_PI <- georef %>%
  inner_join(hapmap, by = c(pi = "V2")) %>% print()

nrow(by_IS) + nrow(by_PI)
nrow(georef)

geo_hap <-  rbind(
  georef %>%
    inner_join(hapmap, by = c(is_no = "V2")),
  georef %>%
    inner_join(hapmap, by = c(pi = "V2")) 
) %>%
  group_by(hapmap_id, avgt_max, avgt_min) %>%
  summarise(count = length(hapmap_id)) %>%
  arrange(-count)  %>%
  dplyr::select(hapmap_id, avgt_min, avgt_max)

long_lat <-  rbind(
  georef %>%
    inner_join(hapmap, by = c(is_no = "V2")),
  georef %>%
    inner_join(hapmap, by = c(pi = "V2")) 
) %>%
  group_by(hapmap_id, Latitude, Longitude) %>%
  summarise(count = length(hapmap_id)) %>%
  arrange(-count)  %>%
  dplyr::select(hapmap_id, Latitude, Longitude)


write.table(long_lat, "long_lat.txt", quote = FALSE)

quartz()




#################################################################################
# Altitude and pH
#################################################################################



setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Research/Data/Phenotype")

georef <- read.csv(file = "alt.ph.csv", na.strings = "NA")
georef

hapmap <- read.table("genotype_ids.txt", header = F)
colnames(hapmap)[1] = "hapmap_id"
hapmap

colnames(georef)
by_IS <- georef %>%
  inner_join(hapmap, by = c(is_no = "V2")) %>% print()

by_PI <- georef %>%
  inner_join(hapmap, by = c(pi = "V2")) %>% print()

nrow(by_IS) + nrow(by_PI)
nrow(georef)

geo_hap <-  rbind(
  georef %>%
    inner_join(hapmap, by = c(is_no = "V2")),
  georef %>%
    inner_join(hapmap, by = c(pi = "V2")) 
) %>%
  group_by(hapmap_id, Alt, topsoil_pH) %>%
  summarise(count = length(hapmap_id)) %>%
  arrange(-count)  %>%
  dplyr::select(hapmap_id, Alt, topsoil_pH)

alt.ph <-  rbind(
  georef %>%
    inner_join(hapmap, by = c(is_no = "V2")),
  georef %>%
    inner_join(hapmap, by = c(pi = "V2")) 
) %>%
  group_by(hapmap_id,Alt, topsoil_pH) %>%
  summarise(count = length(hapmap_id)) %>%
  arrange(-count)  %>%
  dplyr::select(hapmap_id, Alt, topsoil_pH)


write.table(long_lat, "long_lat.txt", quote = FALSE)

quartz()



#############################################################################################################################
#    Laskey Phenotype
#############################################################################################################################



setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Research/Data/Phenotype")

lasky <- read.csv(file = "lasky_phenotype.csv", na.strings = "NA")
head(lasky)

hapmap <- read.table("genotype_ids.txt", header = F)
colnames(hapmap)[1] = "hapmap_id"
hapmap

colnames(lasky)
by_IS <- lasky %>%
  inner_join(hapmap, by = c(is_no = "V2")) %>% print()

by_PI <- georef %>%
  inner_join(hapmap, by = c(pi = "V2")) %>% print()

nrow(by_IS) + nrow(by_PI)
nrow(lasky)



lasky <-  rbind(
  lasky %>%
    inner_join(hapmap, by = c(is_no = "V2")),
  lasky %>%
    inner_join(hapmap, by = c(pi = "V2")) 
) %>%
  group_by(hapmap_id, Latitude, Longitude, Alt, AridityIndex, prec_1,prec_2,prec_3,prec_4,
           prec_5,prec_6,prec_7,prec_8,prec_9,prec_10,prec_11,prec_12, tmax_1,tmax_2,tmax_3,
           tmax_4,tmax_5,tmax_6,tmax_7,tmax_8,tmax_9,tmax_10,tmax_11,tmax_12,tmin_1,tmin_2,
           tmin_3,tmin_4,tmin_5,tmin_6,tmin_7,tmin_8,tmin_9,tmin_10,tmin_11,tmin_12,soil_m,
           topsoil_pH
           ) %>%
  summarise(count = length(hapmap_id)) %>%
  arrange(-count)  %>%
  dplyr::select(hapmap_id, Latitude, Longitude, Alt, AridityIndex, prec_1,prec_2,prec_3,prec_4,
                prec_5,prec_6,prec_7,prec_8,prec_9,prec_10,prec_11,prec_12, tmax_1,tmax_2,tmax_3,
                tmax_4,tmax_5,tmax_6,tmax_7,tmax_8,tmax_9,tmax_10,tmax_11,tmax_12,tmin_1,tmin_2,
                tmin_3,tmin_4,tmin_5,tmin_6,tmin_7,tmin_8,tmin_9,tmin_10,tmin_11,tmin_12,soil_m,
                topsoil_pH)



write.table(lasky, "lasky_phenotype.txt", quote = FALSE)

