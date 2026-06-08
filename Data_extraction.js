/**
 * ==============================================================================
 * SCRIPT: Landsat 8 Surface Temperature & Surface Reflectance Data Extraction
 * PROJECT: Gargano Forest Dieback Monitoring
 * DESCRIPTION: Filters, cloud-masks, calibrates, scales, and exports 100m daily 
 * Landsat 8 mosaics. Generates a metadata CSV for sun elevation.
 * ==============================================================================
 * NOTE FOR USERS: Before running this script, please import your Region of Interest 
 * (ROI) vector/geometry and name the variable 'table'.
 */

// 1. SETUP PARAMETERS
var targetCRS = 'EPSG:32633'; // Use the CRS that best suits your region of interest (roi)
var exportFolder = 'Landsat8_img_100m'; // This folder will appear in your Google Drive as soon as you download the images

// 2. LOAD AND FILTER LANDSAT 8 COLLECTION (TIER 1, LEVEL 2)
var dataset = ee.ImageCollection("LANDSAT/LC08/C02/T1_L2")
  .filterBounds(table) //your roi
  .filterDate("2024-09-01", "2025-04-30")
  .filterMetadata("CLOUD_COVER", "LESS_THAN", 80)
  .map(function(img) {
    // Select core operational bands and preserve native metadata
    var selected = img.select([
      "SR_B1", "SR_B2", "SR_B3", 
      "SR_B4", "SR_B5", "SR_B6", 
      "SR_B7", "ST_B10", "QA_PIXEL"
    ]);
    // Retain solar geometry property for downstream analysis
    return selected.set('SUN_ELEVATION', img.get('SUN_ELEVATION'));
  });
  
// 3. CLOUD AND SHADOW MASKING FUNCTION ---
var cloudMaskUsingQAPixel = function(image) {
  var qaPixel = image.select('QA_PIXEL');
  
  // Isolate cloud and cloud shadow bits using bitwise operations
  var cloud = qaPixel.bitwiseAnd(1 << 3).neq(0);        // Bit 3: Cloud
  var cloudShadow = qaPixel.bitwiseAnd(1 << 4).neq(0);  // Bit 4: Cloud Shadow
  
  // Combine masks to preserve clear-sky pixels only
  var mask = cloud.or(cloudShadow).not(); 
  return image.updateMask(mask);
};
  
// 4. RADIOMETRIC CALIBRATION, SCALING, AND RESAMPLING
var preprocessImage = function(image) {
  // Apply official USGS scaling factors for optical bands (SR)
  var opticalBands = image.select(['SR_B1', 'SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B6', 'SR_B7'])
    .multiply(0.0000275)
    .add(-0.2);
  
  // Convert Thermal Infrared Band (ST_B10) from DN to Celsius
  var thermalBandCelsius = image.select('ST_B10')
    .multiply(0.00341802)
    .add(149.0)
    .subtract(273.15);
  
  // Reassemble image components, clip to roi, and resample to target resolution
  return image
    .addBands(opticalBands, null, true)
    .addBands(thermalBandCelsius, null, true)
    .clip(table)
    .resample('bilinear');
};
  
// Apply preprocessing workflow across the collection
var processedDataset = dataset
  .map(cloudMaskUsingQAPixel) 
  .map(preprocessImage);
  
// 5. DAILY MOSAIC GENERATION (SWATH INTERSECTION RESOLUTION)
var dailyMosaic = ee.ImageCollection(
  processedDataset
    .map(function(image) {
      // Inject standard date string attribute
      return image.set('date', ee.Date(image.get('system:time_start')).format('YYYY-MM-dd'));
    })
    .distinct('date')
    .map(function(dateImage) {
      var date = dateImage.get('date');
      var dailyImages = processedDataset.filter(ee.Filter.date(ee.Date(date), ee.Date(date).advance(1, 'day')));
      // Compute spatial median for overlapping paths/swaths on the same day
      return dailyImages.median().set('system:time_start', ee.Date(date).millis());
    })
);
  
// 6. MAP VISUALIZATION (optional)
var visParam = {bands: ["SR_B4", "SR_B3", "SR_B2"], min: 0, max: 0.4};
Map.addLayer(dailyMosaic, visParam, "Merged Daily Mosaic");
  
// 7. BATCH IMAGE EXPORT (TO GOOGLE DRIVE)
var imageList = dailyMosaic.toList(dailyMosaic.size());
var numImages = imageList.size().getInfo(); 

for (var i = 0; i < numImages; i++) {
  var image = ee.Image(imageList.get(i));
  var date = ee.Date(image.get('system:time_start')).format('YYYY-MM-dd').getInfo();
  
  print('Exporting raster task generated for date:', date);
  Export.image.toDrive({
    image: image,
    description: 'Landsat8_100_' + date,
    scale: 100,
    crs: targetCRS,
    region: table,
    maxPixels: 1e13,
    folder: exportFolder,
    fileFormat: 'GeoTIFF' 
  });
}

// 8. METADATA EXPORT: SOLAR GEOMETRY METRICS (TO CSV)
var imageList_1 = dataset.toList(dataset.size());
var features = [];

for (var i = 0; i < imageList_1.size().getInfo(); i++) {
  var image = ee.Image(imageList_1.get(i));
  var date = ee.Date(image.get('system:time_start')).format('YYYY-MM-dd').getInfo();
  var sunElev = image.get('SUN_ELEVATION').getInfo();
  features.push(ee.Feature(null, {'date': date, 'SUN_ELEVATION': sunElev}));
}

var fc = ee.FeatureCollection(features);

Export.table.toDrive({
  collection: fc,
  description: 'SunElevation_per_image',
  fileFormat: 'CSV',
  folder: exportFolder
});
