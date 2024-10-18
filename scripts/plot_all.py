import xarray as xr
import folium
from folium import plugins
import numpy as np
import geopandas as gpd
import pandas as pd

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

# Add Ohio county lines (now in gray)
ohio_counties = gpd.read_file('https://raw.githubusercontent.com/plotly/datasets/master/geojson-counties-fips.json')
ohio_counties = ohio_counties[ohio_counties['STATE'] == '39']  # FIPS code for Ohio is 39

folium.GeoJson(
    ohio_counties,
    style_function=lambda feature: {
        'fillColor': 'transparent',
        'color': 'gray',
        'weight': 2,
    }
).add_to(m)

# Function to add WPS grid to the map
def add_wps_grid(file_path, color):
    ds_wps = xr.open_dataset(file_path)
    lats_wps = ds_wps.XLAT_M.values[0]
    lons_wps = ds_wps.XLONG_M.values[0]
    
    # Create a list to hold all coordinates
    all_coordinates = []
    
    for i in range(lats_wps.shape[0]):
        row_coords = [[lats_wps[i, j], lons_wps[i, j]] for j in range(lats_wps.shape[1])]
        all_coordinates.append(row_coords)
    
    for j in range(lats_wps.shape[1]):
        col_coords = [[lats_wps[i, j], lons_wps[i, j]] for i in range(lats_wps.shape[0])]
        all_coordinates.append(col_coords)
    
    # Add all lines at once
    folium.PolyLine(all_coordinates, color=color, weight=2, opacity=0.8).add_to(m)
    
    ds_wps.close()

# Add WPS grids
add_wps_grid('../inputs/franklin/geo_em.d01.nc', 'lightcoral')  # Light red
add_wps_grid('../inputs/franklin/geo_em.d02.nc', 'lightblue')   # Light blue

# Add the CAMx2WRF grid to the map (now in green)
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
            color='green',
            weight=1,
            fill=False,
        ).add_to(m)

# Add AirNow stations
airnow_stations = pd.read_csv('../inputs/franklin/franklin_airnow_list.csv')
for _, station in airnow_stations.iterrows():
    folium.Marker(
        location=[station['lat'], station['lon']],
        popup=station['name'],
        icon=folium.Icon(color='orange', icon='info-sign')
    ).add_to(m)

# Add a minimap
minimap = plugins.MiniMap()
m.add_child(minimap)

# Add layer control
folium.LayerControl().add_to(m)

# Add a title (as text on the map, since Folium doesn't have a built-in title function)
title_html = '''
             <h3 align="center" style="font-size:16px"><b>Ohio Counties, WPS Grids, CAMx2WRF Grid, and AirNow Stations - {}</b></h3>
             '''.format(ds.SDATEC.item())

m.get_root().html.add_child(folium.Element(title_html))

# Save the map
m.save("../figures/domain_grid_map_with_counties_wps_and_airnow.html")

print("Map saved as domain_grid_map_with_counties_wps_and_airnow.html")

# Close the datasets
ds.close()
