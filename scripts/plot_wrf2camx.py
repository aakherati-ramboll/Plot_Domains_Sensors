import xarray as xr
import folium
from folium import plugins
import numpy as np
import geopandas as gpd

# Open the netCDF file using xarray
nc_file = '../inputs/franklin/franklin_wrf2camx.2d.1000m.2024101100.nc'
ds = xr.open_dataset(nc_file)

# Extract latitude and longitude
lats = ds.latitude.values
lons = ds.longitude.values

# Calculate the center of the map
center_lat = np.mean(lats)
center_lon = np.mean(lons)

# Create a map centered on the domain
m = folium.Map(location=[center_lat, center_lon], zoom_start=8)

# Add the grid to the map
for i in range(lats.shape[0] - 1):
    for j in range(lats.shape[1] - 1):
        coordinates = [
            [lats[i, j], lons[i, j]],
            [lats[i+1, j], lons[i+1, j]],
            [lats[i+1, j+1], lons[i+1, j+1]],
            [lats[i, j+1], lons[i, j+1]]
        ]
        folium.Polygon(
            locations=coordinates,
            color='red',
            weight=1,
            fill=False,
        ).add_to(m)

# Add Ohio county lines
ohio_counties = gpd.read_file('https://raw.githubusercontent.com/plotly/datasets/master/geojson-counties-fips.json')
ohio_counties = ohio_counties[ohio_counties['STATE'] == '39']  # FIPS code for Ohio is 39

folium.GeoJson(
    ohio_counties,
    style_function=lambda feature: {
        'fillColor': 'transparent',
        'color': 'blue',
        'weight': 2,
    }
).add_to(m)

# Add a minimap
minimap = plugins.MiniMap()
m.add_child(minimap)

# Add layer control
folium.LayerControl().add_to(m)

# Add a title (as text on the map, since Folium doesn't have a built-in title function)
title_html = '''
             <h3 align="center" style="font-size:16px"><b>Domain Grid and Ohio Counties - {}</b></h3>
             '''.format(ds.SDATEC.item())

m.get_root().html.add_child(folium.Element(title_html))

# Save the map
m.save("../figures/domain_grid_map_with_counties.html")

print("Map saved as domain_grid_map_with_counties.html")

# Close the dataset
ds.close()
