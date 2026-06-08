"""
==============================================================================
SCRIPT: Solar Radiation Modeling using GRASS GIS (r.sun)
PROJECT: Gargano Forest Dieback Monitoring
DESCRIPTION: This script automates the GRASS GIS 'r.sun' module to calculate 
             clear-sky global solar radiation for specific days of the year (DOY).
             It imports a high-resolution DEM, computes slope and aspect, and 
             batch-processes the irradiance models.
             
USAGE INSTRUCTIONS & REQUIREMENTS:
    1. GRASS GIS must be installed on the local machine (OSGeo4W recommended for Windows).
    2. The 'grass-session' Python package must be installed (`pip install grass-session`).
    3. Update the 'GRASSBIN' environment variable to point to your local GRASS executable.
    4. Update the 'DEM_PATH' and 'OUTPUT_FOLDER' paths according to your local directories.
==============================================================================
"""

import os
import sys

# 1. CONFIGURATION 
# Set the path to the GRASS GIS executable (Modify if necessary based on your OS/Installation)
os.environ['GRASSBIN'] = r'C:\OSGeo4W\bin\grass84.bat'  

# Input/Output paths, modify with your directories
DEM_PATH = r'C:\DroughtWoods\Gargano\DEM\DEM_100M.tif'
OUTPUT_FOLDER = r'C:\DroughtWoods\Gargano\DEM\solar_radiation_days_pre'

# GRASS Database structure variables
GRASSDATA_DIR = os.path.join(os.path.dirname(OUTPUT_FOLDER), 'grassdata_pre')
LOCATION_NAME = 'radiation_location'
MAPSET_NAME = 'PERMANENT'

from tqdm import tqdm
from grass_session import Session
import grass.script as gs


def check_file_exists(path):
    """Utility function to check if the input DEM exists."""
    if not os.path.isfile(path):
        print(f"File not found: {path}")
        sys.exit(1)

def setup_grass():
    """Initializes the GRASS GIS Location based on the input DEM."""
    location_path = os.path.join(GRASSDATA_DIR, LOCATION_NAME)
    if not os.path.exists(location_path):
        print(f"Creating new GRASS location in {location_path}")
        try:
            gs.create_location(
                dbase=GRASSDATA_DIR,
                location=LOCATION_NAME,
                filename=DEM_PATH,
                epsg='32633'  # UTM Zone 33N projection (matches the study area), replace with the CRS that best suits your study area
            )
        except Exception as e:
            print(f"Error creating GRASS location: {e}")
            sys.exit(1)
    else:
        print("GRASS location already exists.")

def calculate_slope_aspect():
    """Imports the DEM into GRASS and computes topographical parameters (slope and aspect)."""
    with Session(gisdb=GRASSDATA_DIR, location=LOCATION_NAME, mapset=MAPSET_NAME):
        print("Importing DEM and calculating slope/aspect...")
        try:
            gs.run_command('r.in.gdal', input=DEM_PATH, output='dem', overwrite=True)
            info = gs.raster_info('dem')
            print("DEM info:", info)
            
            if info["min"] is None or info["max"] is None:
                print("The DEM is empty or contains only null values.")
                sys.exit(1)

            gs.run_command('r.slope.aspect', elevation='dem', slope='slope', aspect='aspect', overwrite=True)
        except Exception as e:
            print(f"Error calculating slope/aspect: {e}")
            sys.exit(1)

def calculate_radiation_for_day(day_of_year):
    """Executes the r.sun algorithm to compute solar radiation for a specific Day of Year (DOY)."""
    output_raster = os.path.join(OUTPUT_FOLDER, f'solar_rad_day_pre{day_of_year:03d}.tif')
    
    if os.path.exists(output_raster):
        print(f"Day {day_of_year}: file already exists. Skipping.")
        return True

    try:
        with Session(gisdb=GRASSDATA_DIR, location=LOCATION_NAME, mapset=MAPSET_NAME):
            gs.run_command('g.remove', type='raster', name='MASK', flags='f')  # Remove any existing masks
            gs.run_command('g.region', raster='dem')

            region = gs.parse_command('g.region', flags='p')
            print(f"Computational region for day {day_of_year}:", region)

            rad_raster = f'global_radiation_{day_of_year}'

            # Run the r.sun solar irradiance model
            gs.run_command('r.sun',
                elevation='dem',
                slope='slope',
                aspect='aspect',
                day=day_of_year,
                step=1.0,
                dist=100,
                glob_rad=rad_raster,
                overwrite=True)

            info = gs.raster_info(rad_raster)
            print(f" Solar radiation day {day_of_year}: min={info.get('min')}, max={info.get('max')}")

            if info["min"] is None or info["max"] is None:
                print(f"Warning: raster {rad_raster} is empty.")
                return False

            # Export the result to an external GeoTIFF
            gs.run_command('r.out.gdal',
                input=rad_raster,
                output=output_raster,
                format='GTiff',
                createopt='COMPRESS=DEFLATE',
                overwrite=True)

            gs.run_command('g.remove', type='raster', name=rad_raster, flags='f')
            print(f"Completed day {day_of_year}")
            return True
            
    except Exception as e:
        print(f"Error during calculation for day {day_of_year}: {e}")
        return False


def main():
    check_file_exists(DEM_PATH)
    os.makedirs(OUTPUT_FOLDER, exist_ok=True)
    os.makedirs(GRASSDATA_DIR, exist_ok=True)

    setup_grass()
    calculate_slope_aspect()

    print("Enter the days of the year (1-365), separated by commas:")
    try:
        user_input = input("Example: 80, 120, 180\n> ")
        days_to_process = [int(day.strip()) for day in user_input.split(',') if day.strip().isdigit() and 1 <= int(day.strip()) <= 366]
        
        if not days_to_process:
            raise ValueError("No valid day entered.")

        for day in tqdm(days_to_process, desc="Calculating solar radiation"):
            calculate_radiation_for_day(day)

        print("All calculations completed!")

    except ValueError as e:
        print(f"Input error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
