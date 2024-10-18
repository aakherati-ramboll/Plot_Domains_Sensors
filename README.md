
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

## Setup

To set up the project environment and install the required dependencies, follow these steps:

1. Ensure you have Python 3.7+ installed on your system.

2. Clone this repository to your local machine:
   ```
   git clone <repohttps://github.com/your-username/geospatial-data-visualization.git
   cd geospatial-data-visualization
   ```

3. Create a virtual environment:
   ```
   python -m venv .venv
   ```

4. Activate the virtual environment:
   - On Windows:
     ```
     .venv\Scripts\activate
     ```
   - On macOS and Linux:
     ```
     source .venv/bin/activate
     ```

5. Install the required dependencies:
   ```
   pip install -r requirements.txt
   ```

6. Set up the Mapbox token:
   - Create a `.env` file in the root directory of the project if it doesn't exist already.
   - Add your Mapbox token to the `.env` file in the following format:
     ```
     MAPBOX_TOKEN=your_mapbox_token_here
     ```
   - Replace `your_mapbox_token_here` with your actual Mapbox token.

   Note: The `.env` file is included in the `.gitignore` to prevent accidentally sharing your private token. Never commit this file to version control.


This will set up a isolated Python environment with all the necessary libraries installed.

## Contributing

Contributions to enhance the functionality, add new visualization capabilities, or improve existing scripts are welcome. Please submit pull requests or open issues for any bugs or feature requests.


