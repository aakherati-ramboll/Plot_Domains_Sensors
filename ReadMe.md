
# Geospatial Data Visualization for Atmospheric Modeling

This repository contains Python scripts for visualizing geospatial data related to atmospheric modeling and air quality monitoring. The primary focus is on plotting data from various sources, including:

- Weather Research and Forecasting (WRF) model
- Comprehensive Air Quality Model with Extensions (CAMx)
- Other Chemical Transport Models (CTMs)
- Air quality sensor data (e.g., AirNow, PurpleAir, and other desired sensors)

## Features

- Visualization of model domains and grids (WRF, CAMx, etc.)
- Plotting of meteorological and air quality data
- Integration of various data sources (NetCDF, CSV, etc.)
- Customizable map layouts using libraries like Folium
- Support for different map projections and coordinate systems
- Export capabilities for further analysis and presentation

## Scripts

The repository includes several Python scripts for different visualization tasks:

1. `plot_all_franklin.py`: Comprehensive script for plotting multiple data layers, including:
   - Ohio county boundaries
   - WPS (WRF Preprocessing System) grids
   - CAMx to WRF grid mapping
   - AirNow station locations

2. `plot_all.py`: A simpler version for basic plotting of model domains and data

## Data Handling

The scripts are designed to work with various data formats:

- NetCDF files for model outputs (WRF, CAMx)
- CSV files for sensor data (AirNow, PurpleAir)
- GeoJSON for geographical boundaries

## Dependencies

Main libraries used in this project include:

- xarray
- folium
- geopandas
- pandas
- numpy
- shapely

## Usage

To use these scripts, ensure you have the required dependencies installed and the necessary input data. Modify the input and output paths in the scripts as needed for your specific use case.

## Contributing

Contributions to enhance the functionality, add new visualization capabilities, or improve existing scripts are welcome. Please submit pull requests or open issues for any bugs or feature requests.


