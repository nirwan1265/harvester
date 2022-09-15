setwd("~/Desktop")
library(tidyverse)
gantt <- read.csv("africa_harvest.time.csv", header = T)
acts <- c("Botswana",
          "Democratic Republic of Congo",
          "Eswatini (Swaziland)",
          "Lesotho",
          "Malawi",
          "Mozambique",
          "Namibia",
          "South Africa",
          "United Republic of Tanzania(Bimodal, Masika)",
          "United Republic of Tanzania(Unimodal, Msimu)",
          "Zambia",
          "Zimbabwe",
          "Burundi(Sorghum A)",
          "Burundi(Sorghum B)",
          "Ethiopia",
          "Eritrea",
          "Kenya",
          "Somalia(Deyr)",
          "Somalia(Gu)",
          "South Sudan(South, Main)",
          "South Sudan(South, Second)",
          "South Sudan(Unimodal)",
          "Sudan",
          "Uganda",
          "Benin",
          "Burkina Faso",
          "Chad",
          "Central African Republic(North)",
          "Gambia",
          "Ghana",
          "Guinea",
          "Mali",
          "Nigeria",
          "Niger",
          "Sierra Leone",
          "Senegal",
          "Togo")
els <- c("Plant","Mid-Season","Harvest")
reg <- c("Southern Africa", "East Africa","West Africa")
head(gantt)

g.gantt <- gather(gantt, "state", "date", 4:5) %>% mutate(date = as.Date(date, "%Y.%m.%d"), Country=factor(Country, acts[length(acts):1]), Timeline=factor(Timeline, els), Region=factor(Region,reg))
head(g.gantt)
dev.off()
tiff("test.tiff", units="in", width=20, height=10, res=300)
ggplot(g.gantt, aes(date, Country, color = Timeline, group=Item)) +
  geom_line(size = 5) +
  facet_grid(Region~., scales = "free")+
  scale_color_manual(values=c("#1e847f","grey","#ecc19c"), name="Timeline") +
  labs(x="Planting year", y=NULL, title="Planting timeline")
