import xarray as xr
import folium
from folium import plugins
import numpy as np
import geopandas as gpd
import pandas as pd
from shapely.geometry import Polygon, Point, LineString
import pyproj
from typing import Tuple, List
from dotenv import load_dotenv
import os

load_dotenv()

# Input files and parameters
NC_FILE: str = '../inputs/franklin/franklin_wrf2camx.2d.1000m.2024101100.nc'
WPS_D01_FILE: str = '../inputs/franklin/geo_em.d01.nc'
WPS_D02_FILE: str = '../inputs/franklin/geo_em.d02.nc'
AIRNOW_STATIONS_FILE: str = '../inputs/franklin/franklin_airnow_list.csv'
OHIO_COUNTIES_URL: str = 'https://raw.githubusercontent.com/plotly/datasets/master/geojson-counties-fips.json'
OUTPUT_HTML: str = "../figures/domain_grid_map_with_counties_wps_and_airnow.html"
OUTPUT_GPKG: str = "../outputs/franklin_data.gpkg"
MAPBOX_TOKEN: str = os.getenv("MAPBOX_TOKEN")

def load_data() -> Tuple[xr.Dataset, gpd.GeoDataFrame, pd.DataFrame]:
    """
    Load necessary data from files.

    Returns:
        Tuple[xr.Dataset, gpd.GeoDataFrame, pd.DataFrame]: A tuple containing:
            - xarray Dataset with WRF2CAMx data
            - GeoDataFrame with Ohio counties
            - DataFrame with AirNow stations
    """
    ds: xr.Dataset = xr.open_dataset(NC_FILE)
    ohio_counties: gpd.GeoDataFrame = gpd.read_file(OHIO_COUNTIES_URL)
    ohio_counties = ohio_counties[ohio_counties['STATE'] == '39']  # FIPS code for Ohio is 39
    airnow_stations: pd.DataFrame = pd.read_csv(AIRNOW_STATIONS_FILE)
    return ds, ohio_counties, airnow_stations

def create_map(ds: xr.Dataset) -> folium.Map:
    """
    Create a Folium map centered on the dataset's domain.

    Args:
        ds (xr.Dataset): The dataset containing latitude and longitude information.

    Returns:
        folium.Map: A Folium map object with various tile layers added.
    """
    lats, lons = ds.latitude.values, ds.longitude.values
    center_lat, center_lon = np.mean(lats), np.mean(lons)
    m = folium.Map(location=[center_lat, center_lon], zoom_start=8)

    # Add various tile layers
    folium.TileLayer('openstreetmap').add_to(m)
    folium.TileLayer('cartodbdark_matter', name='Dark Map').add_to(m)
    folium.TileLayer('cartodbpositron', name='Light Map').add_to(m)
    folium.TileLayer(
        tiles=f'https://api.mapbox.com/styles/v1/mapbox/streets-v11/tiles/{{z}}/{{x}}/{{y}}?access_token={MAPBOX_TOKEN}',
        attr='Mapbox',
        name='Mapbox Streets',
    ).add_to(m)
    folium.TileLayer(
        tiles=f'https://api.mapbox.com/styles/v1/mapbox/outdoors-v11/tiles/{{z}}/{{x}}/{{y}}?access_token={MAPBOX_TOKEN}',
        attr='Mapbox',
        name='Mapbox Outdoors',
    ).add_to(m)
    folium.TileLayer(
        tiles=f'https://api.mapbox.com/styles/v1/mapbox/light-v10/tiles/{{z}}/{{x}}/{{y}}?access_token={MAPBOX_TOKEN}',
        attr='Mapbox',
        name='Mapbox Light',
    ).add_to(m)
    folium.TileLayer(
        tiles=f'https://api.mapbox.com/styles/v1/mapbox/dark-v10/tiles/{{z}}/{{x}}/{{y}}?access_token={MAPBOX_TOKEN}',
        attr='Mapbox',
        name='Mapbox Dark',
    ).add_to(m)
    folium.TileLayer(
        tiles=f'https://api.mapbox.com/styles/v1/mapbox/satellite-v9/tiles/{{z}}/{{x}}/{{y}}?access_token={MAPBOX_TOKEN}',
        attr='Mapbox',
        name='Satellite',
    ).add_to(m)

    # Add measure control
    plugins.MeasureControl(position='topright', primary_length_unit='kilometers', secondary_length_unit='miles').add_to(m)

    return m

def add_ohio_counties(m: folium.Map, ohio_counties: gpd.GeoDataFrame) -> None:
    """
    Add Ohio county boundaries to the map.

    Args:
        m (folium.Map): The Folium map object to add the counties to.
        ohio_counties (gpd.GeoDataFrame): GeoDataFrame containing Ohio county geometries.
    """
    folium.GeoJson(
        ohio_counties,
        style_function=lambda feature: {
            'fillColor': 'transparent',
            'color': 'gray',
            'weight': 2,
        }
    ).add_to(m)

def add_wps_grid(m: folium.Map, file_path: str, color: str) -> gpd.GeoDataFrame:
    """
    Add WPS grid to the map and return it as a GeoDataFrame.

    Args:
        m (folium.Map): The Folium map object to add the grid to.
        file_path (str): Path to the WPS file.
        color (str): Color of the grid lines.

    Returns:
        gpd.GeoDataFrame: A GeoDataFrame containing the WPS grid lines.
    """
    ds_wps: xr.Dataset = xr.open_dataset(file_path)
    lats_wps, lons_wps = ds_wps.XLAT_M.values[0], ds_wps.XLONG_M.values[0]
    
    all_coordinates: List[List[List[float]]] = []
    for i in range(lats_wps.shape[0]):
        all_coordinates.append([[lats_wps[i, j], lons_wps[i, j]] for j in range(lats_wps.shape[1])])
    for j in range(lats_wps.shape[1]):
        all_coordinates.append([[lats_wps[i, j], lons_wps[i, j]] for i in range(lats_wps.shape[0])])
    
    folium.PolyLine(all_coordinates, color=color, weight=2, opacity=0.8).add_to(m)
    
    lines: List[LineString] = [LineString([(lon, lat) for lat, lon in coords]) for coords in all_coordinates]
    wps_gdf: gpd.GeoDataFrame = gpd.GeoDataFrame(geometry=lines, crs="EPSG:4326")
    wps_gdf['grid_type'] = 'WPS'
    
    ds_wps.close()
    return wps_gdf

def add_camx2wrf_grid(m: folium.Map, lats: np.ndarray, lons: np.ndarray) -> gpd.GeoDataFrame:
    """
    Add CAMx2WRF grid to the map and return it as a GeoDataFrame.

    Args:
        m (folium.Map): The Folium map object to add the grid to.
        lats (np.ndarray): Array of latitude values.
        lons (np.ndarray): Array of longitude values.

    Returns:
        gpd.GeoDataFrame: A GeoDataFrame containing the CAMx2WRF grid polygons.
    """
    camx2wrf_polygons: List[Polygon] = []
    for i in range(lats.shape[0] - 1):
        for j in range(lats.shape[1] - 1):
            coordinates: List[Tuple[float, float]] = [
                (lons[i, j], lats[i, j]),
                (lons[i+1, j], lats[i+1, j]),
                (lons[i+1, j+1], lats[i+1, j+1]),
                (lons[i, j+1], lats[i, j+1]),
                (lons[i, j], lats[i, j])
            ]
            polygon: Polygon = Polygon(coordinates)
            camx2wrf_polygons.append(polygon)
            folium.Polygon(
                locations=[(lat, lon) for lon, lat in coordinates],
                color='green',
                weight=1,
                fill=False,
            ).add_to(m)
    
    camx2wrf_gdf: gpd.GeoDataFrame = gpd.GeoDataFrame(geometry=camx2wrf_polygons, crs="EPSG:4326")
    camx2wrf_gdf['grid_type'] = 'CAMx2WRF'
    return camx2wrf_gdf

def add_airnow_stations(m: folium.Map, airnow_stations: pd.DataFrame) -> gpd.GeoDataFrame:
    """
    Add AirNow stations to the map and return them as a GeoDataFrame.

    Args:
        m (folium.Map): The Folium map object to add the stations to.
        airnow_stations (pd.DataFrame): DataFrame containing AirNow station information.

    Returns:
        gpd.GeoDataFrame: A GeoDataFrame containing the AirNow stations.
    """
    for _, station in airnow_stations.iterrows():
        folium.Marker(
            location=[station['lat'], station['lon']],
            popup=station['name'],
            icon=folium.Icon(color='orange', icon='info-sign')
        ).add_to(m)
    
    airnow_gdf: gpd.GeoDataFrame = gpd.GeoDataFrame(
        airnow_stations,
        geometry=gpd.points_from_xy(airnow_stations.lon, airnow_stations.lat),
        crs="EPSG:4326"
    )
    airnow_gdf['type'] = 'AirNow Station'
    return airnow_gdf

def add_map_features(m: folium.Map, ds: xr.Dataset) -> None:
    """
    Add additional features to the map such as minimap, layer control, and title.

    Args:
        m (folium.Map): The Folium map object to add features to.
        ds (xr.Dataset): The dataset containing date information for the title.
    """
    minimap: plugins.MiniMap = plugins.MiniMap()
    m.add_child(minimap)
    folium.LayerControl().add_to(m)
    title_html: str = f'''
                 <h3 align="center" style="font-size:16px"><b>Ohio Counties, WPS Grids, CAMx2WRF Grid, and AirNow Stations - {ds.SDATEC.item()}</b></h3>
                 '''
    m.get_root().html.add_child(folium.Element(title_html))

def export_data(ohio_counties: gpd.GeoDataFrame, wps_d01_gdf: gpd.GeoDataFrame, wps_d02_gdf: gpd.GeoDataFrame, camx2wrf_gdf: gpd.GeoDataFrame, airnow_gdf: gpd.GeoDataFrame) -> None:
    """
    Export all data to a GeoPackage file.

    Args:
        ohio_counties (gpd.GeoDataFrame): GeoDataFrame containing Ohio county geometries.
        wps_d01_gdf (gpd.GeoDataFrame): GeoDataFrame containing WPS domain 1 grid.
        wps_d02_gdf (gpd.GeoDataFrame): GeoDataFrame containing WPS domain 2 grid.
        camx2wrf_gdf (gpd.GeoDataFrame): GeoDataFrame containing CAMx2WRF grid.
        airnow_gdf (gpd.GeoDataFrame): GeoDataFrame containing AirNow stations.
    """
    ohio_counties.to_file(OUTPUT_GPKG, layer='ohio_counties', driver="GPKG")
    wps_d01_gdf.to_file(OUTPUT_GPKG, layer='wps_grid_d01', driver="GPKG")
    wps_d02_gdf.to_file(OUTPUT_GPKG, layer='wps_grid_d02', driver="GPKG")
    camx2wrf_gdf.to_file(OUTPUT_GPKG, layer='camx2wrf_grid', driver="GPKG")
    airnow_gdf.to_file(OUTPUT_GPKG, layer='airnow_stations', driver="GPKG")
    print(f"Data exported as GeoPackage: {OUTPUT_GPKG}")

def main() -> None:
    """
    Main function to orchestrate the creation of the map and data export.
    """
    ds, ohio_counties, airnow_stations = load_data()
    m: folium.Map = create_map(ds)
    add_ohio_counties(m, ohio_counties)
    wps_d01_gdf: gpd.GeoDataFrame = add_wps_grid(m, WPS_D01_FILE, 'lightcoral')
    wps_d02_gdf: gpd.GeoDataFrame = add_wps_grid(m, WPS_D02_FILE, 'lightblue')
    camx2wrf_gdf: gpd.GeoDataFrame = add_camx2wrf_grid(m, ds.latitude.values, ds.longitude.values)
    airnow_gdf: gpd.GeoDataFrame = add_airnow_stations(m, airnow_stations)
    add_map_features(m, ds)
    m.save(OUTPUT_HTML)
    print(f"Map saved as {OUTPUT_HTML}")
    export_data(ohio_counties, wps_d01_gdf, wps_d02_gdf, camx2wrf_gdf, airnow_gdf)
    ds.close()

if __name__ == "__main__":
    main()
