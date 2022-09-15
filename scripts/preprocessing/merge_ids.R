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
