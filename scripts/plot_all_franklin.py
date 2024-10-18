import xarray as xr
import folium
from folium import plugins
import numpy as np
import geopandas as gpd
import pandas as pd
from shapely.geometry import Polygon, Point, LineString
import pyproj

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

# Function to add WPS grid to the map and return GeoDataFrame
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
    
    # Create GeoDataFrame for WPS grid
    lines = [LineString([(lon, lat) for lat, lon in coords]) for coords in all_coordinates]
    wps_gdf = gpd.GeoDataFrame(geometry=lines, crs="EPSG:4326")
    wps_gdf['grid_type'] = 'WPS'
    
    ds_wps.close()
    return wps_gdf

# Add WPS grids
wps_d01_gdf = add_wps_grid('../inputs/franklin/geo_em.d01.nc', 'lightcoral')  # Light red
wps_d02_gdf = add_wps_grid('../inputs/franklin/geo_em.d02.nc', 'lightblue')   # Light blue

# Add the CAMx2WRF grid to the map (now in green) and create GeoDataFrame
camx2wrf_polygons = []
for i in range(lats.shape[0] - 1):
    for j in range(lats.shape[1] - 1):
        coordinates = [
            (lons[i, j], lats[i, j]),
            (lons[i+1, j], lats[i+1, j]),
            (lons[i+1, j+1], lats[i+1, j+1]),
            (lons[i, j+1], lats[i, j+1]),
            (lons[i, j], lats[i, j])
        ]
        polygon = Polygon(coordinates)
        camx2wrf_polygons.append(polygon)
        folium.Polygon(
            locations=[(lat, lon) for lon, lat in coordinates],
            color='green',
            weight=1,
            fill=False,
        ).add_to(m)

camx2wrf_gdf = gpd.GeoDataFrame(geometry=camx2wrf_polygons, crs="EPSG:4326")
camx2wrf_gdf['grid_type'] = 'CAMx2WRF'

# Add AirNow stations
airnow_stations = pd.read_csv('../inputs/franklin/franklin_airnow_list.csv')
for _, station in airnow_stations.iterrows():
    folium.Marker(
        location=[station['lat'], station['lon']],
        popup=station['name'],
        icon=folium.Icon(color='orange', icon='info-sign')
    ).add_to(m)

# Create GeoDataFrame for AirNow stations
airnow_gdf = gpd.GeoDataFrame(
    airnow_stations,
    geometry=gpd.points_from_xy(airnow_stations.lon, airnow_stations.lat),
    crs="EPSG:4326"
)
airnow_gdf['type'] = 'AirNow Station'

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

# Export data as GeoPackage
output_gpkg = "../outputs/franklin_data.gpkg"
ohio_counties.to_file(output_gpkg, layer='ohio_counties', driver="GPKG")
wps_d01_gdf.to_file(output_gpkg, layer='wps_grid_d01', driver="GPKG")
wps_d02_gdf.to_file(output_gpkg, layer='wps_grid_d02', driver="GPKG")
camx2wrf_gdf.to_file(output_gpkg, layer='camx2wrf_grid', driver="GPKG")
airnow_gdf.to_file(output_gpkg, layer='airnow_stations', driver="GPKG")

print(f"Data exported as GeoPackage: {output_gpkg}")

# Close the datasets
ds.close()
