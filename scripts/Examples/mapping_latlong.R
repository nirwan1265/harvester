
# The single argument to this function, points, is a data.frame in which:
#   - column 1 contains the longitude in degrees
#   - column 2 contains the latitude in degrees

# set up some points to test
points = data.frame(lon=c(0, 5, 10, 15, 20), lat=c(51.5, 50, 48.5, 47, 44.5))
# plot them on a map
points(points$lon, points$lat, col="red")
points
# get a list of country names
coords2country(points)
# returns [1] "UK"         "Belgium"    "Germany"    "Austria"    "Yugoslavia"
# number 5 should probably be in Serbia...