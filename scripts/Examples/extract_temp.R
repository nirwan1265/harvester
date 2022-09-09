#website:https://www.benjaminbell.co.uk/2018/01/extracting-data-and-making-climate-maps.html
setwd("~/Downloads/wc2.1_2.5m_tmin_2000-2009/")
system("ls")
temp1 <- raster("wc2.1_2.5m_tmin_2009-11.tif")
temp2 <- raster("wc2.1_2.5m_tmin_2009-12.tif")
temp1

# Create a data.frame with sample site coordinates
site <- c("Manchester", "Liverpool", "Oxford", "London")
lon <- c(-2.24, -2.98, -1.25, -0.11)
lat <- c(53.47, 53.4, 51.74, 51.49)
samples <- data.frame(site, lon, lat, row.names="site")
samples



# Extract data from WorldClim for your sites
temp.data <- samples 
temp.data$Nov18 <- extract(temp1, samples)
temp.data$Dec18 <- extract(temp2, samples)

tempcol <- colorRampPalette(c("purple", "blue", "skyblue", "green", "lightgreen", "yellow", "orange", "red", "darkred"))
plot(temp1, col=tempcol(100))

plot(temp1, xlim=c(-12, 4), ylim=c(48, 64), zlim=c(-10,30), col=tempcol(100))
