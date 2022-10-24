
# Packages
from geopy.geocoders import Nominatim
import os
import pandas as pd
import numpy as np
import country_converter as coco

# Set working directory
os.chdir("/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Research/Data/Phenotype")


# Open latitude longitude file
long_lat = pd.read_table("long_lat.txt", delim_whitespace=True)

# Separating lat and long and changing to np arrays
lat = long_lat[['Latitude']]
long = long_lat[['Longitude']]
lat = lat.to_numpy()
lat = np.concatenate(lat)
long = long.to_numpy()
long = np.concatenate(long)

# initialize Nominatim API
geolocator = Nominatim(user_agent="geoapiExercises")


#Long and lat inputs 
c = []
cd = []
code = []
for i,j in zip(long,lat):
    i = str(i)
    j = str(j)
    c = geolocator.reverse(j+","+i)
    address = c.raw['address']
    code = address.get('country_code')
    cd.append(code)
    
    
#Saving
np.savetxt('/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Research/Data/Phenotype/country_code.csv', cd, delimiter = ',')

    

