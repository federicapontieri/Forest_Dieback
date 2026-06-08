"""
==============================================================================
SCRIPT: MODIS MOD16A2 Evapotranspiration (ET) Data Extraction (Python/Colab)
PROJECT: Gargano Forest Dieback Monitoring
DESCRIPTION: This script authenticates with Google Earth Engine (through Google Colab), loads the 
             MODIS ET 8-day composite collection, applies the official 
             scaling factor, clips the data to the Gargano Region of Interest 
             (roi), and exports the time-series mosaics to Google Drive 
             projected in UTM 33N.
==============================================================================
"""

# 0. ENVIRONMENT SETUP (For Google Colab)
# Installation of the Earth Engine API
!pip install earthengine-api

# 1. IMPORT REQUIRED LIBRARIES
import ee
import folium
from folium import plugins

# 2. AUTHENTICATION AND INITIALIZATION 
# Trigger the authentication flow and initialize the GEE project
ee.Authenticate()
ee.Initialize(project='ee-xxxxxxx') # This must be replaced with your GEE project

# 3. REGION OF INTEREST (ROI) DEFINITION
# Load the pre-defined bounding box asset for the Gargano Promontory
aoi = ee.FeatureCollection('projects/ee-xxxxxxx/assets/roi') # This asset must be uploaded in your GEE projects; replace it with your directory

# --- 4. DATASET LOADING AND FILTERING ---
# Load MOD16A2 Version 6.1 (8-Day Evapotranspiration composite)
# Filter by the period and select the core ET band
collection = ee.ImageCollection('MODIS/061/MOD16A2') \
    .filterDate('2024-03-01', '2025-05-31') \
    .select('ET')  

# 5. BATCH EXPORT ITERATION
# Convert the Earth Engine ImageCollection into a List to iterate over it client-side
image_list = collection.toList(collection.size())

# Iterate through each image in the collection
for i in range(image_list.size().getInfo()):
    img = ee.Image(image_list.get(i))
    
    # Extract the acquisition date to use in the output filename
    date = img.date().format('YYYY_MM_dd').getInfo()

    # --- 6. PRE-PROCESSING ---
    # Apply the official MODIS scaling factor for ET (0.1) and clip to the AOI
    processed_img = img.multiply(0.1).clip(aoi)  

    # --- 7. EXPORT TO GOOGLE DRIVE ---
    # Configure the export task ensuring proper CRS (EPSG:32633 - UTM Zone 33N) 
    # and spatial resolution (500m) to match the project's analytical framework
    task = ee.batch.Export.image.toDrive(
        image=processed_img,
        description=f'ET_{date}_Gargano',
        folder='GEE_ET_post',
        fileNamePrefix=f'ET_{date}_UTM33N',
        region=aoi.geometry(),
        scale=500,
        crs='EPSG:32633',
        maxPixels=1e13,
        fileFormat='GeoTIFF'
    )
    
    # Start the export task
    task.start()
    
    # Log the progress in the console
    print(f'Export started: ET_{date} | Task ID: {task.id}')
