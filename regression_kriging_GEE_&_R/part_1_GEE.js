/**
 * ==============================================================================
 * SCRIPT: Random Forest Regression for Evapotranspiration (ET) Modeling
 * PROJECT: Gargano Forest Dieback Monitoring
 * DESCRIPTION:
 *   This script implements a Random Forest regression pipeline in Google Earth
 *   Engine (GEE) to predict Evapotranspiration (ET) across the Gargano forest
 *   area. The model uses three predictor variables:
 *     - LST : Land Surface Temperature
 *     - CLC : Corine Land Cover classification (static, 2018)
 *     - SR  : Solar Radiation
 *
 *   The workflow covers:
 *     1. Temporal joining of multi-source image collections
 *     2. Training/validation dataset preparation (70/30 split)
 *     3. Random Forest model training (500 trees)
 *     4. Accuracy assessment: RMSE, R², MAPE on the validation set
 *     5. Spatial prediction and residual computation per image
 *     6. Batch export of predicted ET images and residual point tables to Drive
 *        (residuals are intended for subsequent Kriging interpolation in GIS)
 * WORKFLOW POSITION:
 *   This is PART 1 of a two-step pipeline. The residuals exported here feed
 *   into a separate geostatistical analysis (Kriging) to correct spatial
 *   autocorrelation in the RF predictions. 
 *   The second part will be carried out in R environment.
 *
 * ASSETS REQUIRED (GEE), the assets must be uploaded in GEE, replace the directories with the ones from your GEE projects:
 *   - ET ImageCollection: projects/ee-xxxxxxxx/assets/ET
 *   - LST ImageCollection: projects/ee-xxxxxxxx/assets/LST
 *   - SR ImageCollection: projects/ee-xxxxxxxx/assets/SR
 *   - CLC Image: projects/ee-xxxxxxxx/assets/CLC_2018_Gargano
 *   - Area of Interest: imported as 'table' (FeatureCollection)
 *
 * OUTPUTS (exported to Google Drive, folder: GEE_Output):
 *   - accuracy_rf_model.csv: RMSE, R², MAPE metrics
 *   - pred_rf_<YYYY_MM_DD>.tif: Predicted ET raster per date (EPSG:32633, 100m)
 *   - residuals_<YYYY_MM_DD>.csv: Sampled residual points per date (for Kriging)
 * ==============================================================================
 */

// 1. AREA OF INTEREST & ASSET CONFIGURATION
// The Area of Interest (AOI) is imported as a GEE asset named 'table'.
// It should be a FeatureCollection defining the Gargano study region.
var areaOfInterest = table;

// Asset IDs for the three temporal ImageCollections and the static CLC layer.
// Update these paths if assets are moved or renamed.
var ET_COLLECTION_ID  = 'projects/ee-xxxxxxxx/assets/ET';
var LST_COLLECTION_ID = 'projects/ee-xxxxxxxx/assets/LST';
var SR_COLLECTION_ID  = 'projects/ee-xxxxxxxx/assets/SR';
var CLC_ASSET_ID      = 'projects/ee-xxxxxxxx/assets/CLC_2018_Gargano';


// 2. LOADING COLLECTIONS & TEMPORAL JOINING

// Load the three temporal ImageCollections and the static CLC image.
var etCollection  = ee.ImageCollection(ET_COLLECTION_ID);
var lstCollection = ee.ImageCollection(LST_COLLECTION_ID);
var srCollection  = ee.ImageCollection(SR_COLLECTION_ID);

// CLC is a single static image; rename the band for clarity in downstream steps.
var clc = ee.Image(CLC_ASSET_ID).rename('CLC');

// Define a maximum temporal difference filter: images from different collections
// are considered a valid match only if acquired within 5 days of each other.
// The threshold is expressed in milliseconds.
var maxTimeDiff = ee.Filter.maxDifference({
  difference: 5 * 24 * 60 * 60 * 1000, // 5 days in milliseconds
  leftField: 'system:time_start',
  rightField: 'system:time_start'
});

// Use saveBest joins so that each ET image is paired with the temporally
// closest LST and SR images (within the 5-day window).
// The matched images are stored as properties ('lstMatch', 'srMatch') on the
// resulting ET features, allowing them to be retrieved later.
var saveBestLstJoin = ee.Join.saveBest({ matchKey: 'lstMatch', measureKey: 'timeDiff' });
var saveBestSrJoin  = ee.Join.saveBest({ matchKey: 'srMatch',  measureKey: 'timeDiff' });

// Join LST to ET.
var etWithLst = saveBestLstJoin.apply(etCollection, lstCollection, maxTimeDiff);

// Join SR to the already-joined ET+LST collection.
var etWithLstSr = saveBestSrJoin.apply(etWithLst, srCollection, maxTimeDiff);


// 3. TRAINING & VALIDATION DATASET PREPARATION

// Keep only ET images that successfully matched both an LST and an SR image.
// This removes any edge-case entries where the temporal filter found no pair.
var successfulMatches = etWithLstSr.filter(
  ee.Filter.notNull(['lstMatch', 'srMatch'])
);

// Debug prints — visible in the GEE Console.
print('Matched triplets (ET + LST + SR):', successfulMatches);
print('Number of temporal matches found:', successfulMatches.size());

/**
 * createTrainingImage
 * -------------------
 * Converts a single matched ET feature (with embedded lstMatch and srMatch
 * properties) into a multi-band image suitable for sampling training points.
 *
 * Output bands:
 *   LST  — Land Surface Temperature (band b8 from the LST asset)
 *   CLC  — Corine Land Cover class (static, reprojected to LST grid)
 *   SR   — Solar Radiation
 *   ET   — Evapotranspiration (resampled to 100 m, bilinear)
 *
 * All bands are co-registered to the LST native projection at 100 m scale.
 *
 * @param {ee.Feature} feature - A feature from successfulMatches.
 * @returns {ee.Image} - A 4-band image (LST, CLC, SR, ET).
 */
var createTrainingImage = function(feature) {
  var et  = ee.Image(feature);
  var lst = ee.Image(feature.get('lstMatch')).select('b8').rename('LST');
  var sr  = ee.Image(feature.get('srMatch')).rename('SR');

  // Use LST's projection as the common spatial reference for all layers.
  var targetProjection = lst.projection();
  var targetScale = 100; // metres

  // Resample ET to 100 m using bilinear interpolation, then reproject to
  // match LST's coordinate system exactly.
  var et_resampled = et
    .resample('bilinear')
    .reproject({ crs: targetProjection, scale: targetScale })
    .rename('ET');

  // Stack all predictor bands and the target variable into one image.
  return ee.Image.cat([lst, clc, sr, et_resampled]);
};

/**
 * sampleRegions
 * -------------
 * Randomly samples up to 2000 pixels from a given stacked image within the
 * AOI. Each pixel becomes one training/validation sample (a GEE Feature with
 * band values as properties).
 *
 * @param {ee.Image} image - A 4-band image from createTrainingImage.
 * @returns {ee.FeatureCollection} - Sampled point features.
 */
var sampleRegions = function(image) {
  return image.sample({
    region: areaOfInterest.geometry(),
    scale: 100,
    numPixels: 2000,
    geometries: true
  });
};

// Build the full sample dataset by applying both functions across the collection
// and flattening the result into a single FeatureCollection.
var fullData = ee.ImageCollection(successfulMatches.map(createTrainingImage))
                 .map(sampleRegions)
                 .flatten();

// Append a uniformly distributed random number [0, 1) to each sample.
// This is used to perform a reproducible random train/validation split.
fullData = fullData.randomColumn();

// 70 % training / 30 % validation split based on the random column.
var trainingData   = fullData.filter(ee.Filter.lt('random',  0.7));
var validationData = fullData.filter(ee.Filter.gte('random', 0.7));

print('Training samples:',   trainingData.size());
print('Validation samples:',  validationData.size());


// 4. RANDOM FOREST MODEL TRAINING & VALIDATION

//Model Training

print('Starting Random Forest model training...');

// Instantiate a Random Forest regressor with 500 trees.
// REGRESSION output mode returns continuous numeric predictions (not class labels).
var rfRegressor = ee.Classifier.smileRandomForest({ numberOfTrees: 500 })
  .setOutputMode('REGRESSION');

// Train the model on the 70 % training partition.
// Predictor variables: LST, CLC (land cover class code), SR.
// Target variable:     ET.
rfRegressor = rfRegressor.train({
  features: trainingData,
  classProperty: 'ET',
  inputProperties: ['LST', 'CLC', 'SR']
});

print('Model training complete.');


// Accuracy Assessment on the Validation Set

print('Computing accuracy metrics on the validation set...');

// Apply the trained model to the 30 % validation partition.
// Predictions are stored in a new property called 'classification'.
var validationPredictions = validationData.classify(rfRegressor);


// Root Mean Squared Error (RMSE)
// Compute the squared difference between observed and predicted ET for each
// validation sample, then take the mean and its square root.
var calculateSquaredError = function(feature) {
  var et_real = ee.Number(feature.get('ET'));
  var et_pred = ee.Number(feature.get('classification'));
  var error   = et_real.subtract(et_pred);
  return feature.set('squared_error', error.pow(2));
};

var validationWithErrors = validationPredictions.map(calculateSquaredError);

var meanSquaredError = ee.Number(
  validationWithErrors.reduceColumns({
    reducer: ee.Reducer.mean(),
    selectors: ['squared_error']
  }).get('mean')
);

var rmse = meanSquaredError.sqrt();
print('RF Model RMSE (on validation set):', rmse);


// Coefficient of Determination (R²)
// R² = 1 - (SS_res / SS_tot)
// SS_res: sum of squared residuals (= MSE × n)
// SS_tot: total sum of squares around the observed ET mean
var ss_res = meanSquaredError.multiply(validationWithErrors.size());

var et_mean = ee.Number(
  validationData.reduceColumns(ee.Reducer.mean(), ['ET']).get('mean')
);

var sstot = ee.Number(
  validationData.map(function(feature) {
    var et_real = ee.Number(feature.get('ET'));
    var diff    = et_real.subtract(et_mean);
    return feature.set('total_error_sq', diff.pow(2));
  }).reduceColumns(ee.Reducer.sum(), ['total_error_sq']).get('sum')
);

var r2 = ee.Number(1).subtract(ss_res.divide(sstot));
print('RF Model R² (on validation set):', r2);


// Mean Absolute Percentage Error (MAPE)
// MAPE = mean(|ET_obs - ET_pred| / ET_obs) × 100
// When ET_obs = 0 the percentage error is set to 0 to avoid division by zero.
var calculateMAPE = function(feature) {
  var et_real   = ee.Number(feature.get('ET'));
  var et_pred   = ee.Number(feature.get('classification'));
  var abs_error = et_real.subtract(et_pred).abs();

  var percentage_error = ee.Algorithms.If(
    et_real.eq(0),
    0,                           // Guard against division by zero
    abs_error.divide(et_real)
  );
  return feature.set('abs_perc_error', percentage_error);
};

var validationWithMAPE = validationPredictions.map(calculateMAPE);

var mape = ee.Number(
  validationWithMAPE.reduceColumns(ee.Reducer.mean(), ['abs_perc_error']).get('mean')
).multiply(100); // Convert to percentage

print('RF Model MAPE % (on validation set):', mape);


// Export Accuracy Metrics to Google Drive
// Saves RMSE, R², and MAPE as a single-row CSV for record-keeping.
var accuracyFeature = ee.Feature(null, {
  'RMSE': rmse,
  'R_squared': r2,
  'MAPE': mape
});

Export.table.toDrive({
  collection: ee.FeatureCollection([accuracyFeature]),
  description: 'export_accuracy_rf_model',
  folder: 'GEE_Output',
  fileNamePrefix: 'accuracy_rf_model',
  fileFormat: 'CSV'
});



// 5. SPATIAL PREDICTION & RESIDUAL COMPUTATION

/**
 * createPrediction
 * ----------------
 * Applies the trained RF model to a single matched ET triplet to produce:
 *   - ET_pred: RF-predicted ET image (band name: 'ET_pred')
 *   - residual: Observed ET − Predicted ET (stored as an image property
 *                 'residuals_image' and later sampled for Kriging input)
 *
 * The original acquisition date is preserved as a metadata property so it can
 * be used to construct meaningful export file names in Section 6.
 *
 * @param {ee.Feature} feature - A feature from successfulMatches.
 * @returns {ee.Image} - Predicted ET image with metadata and residuals attached.
 */
var createPrediction = function(feature) {
  var et  = ee.Image(feature);
  var lst = ee.Image(feature.get('lstMatch')).select('b8').rename('LST');
  var sr  = ee.Image(feature.get('srMatch')).rename('SR');

  var targetProjection = lst.projection();
  var targetScale = 100; // metres

  // Stack the predictor bands (no ET needed here, this is the prediction step).
  var predictors = ee.Image.cat([lst, clc, sr]);

  // Resample observed ET to 100 m for residual computation.
  var et_resampled = et
    .resample('bilinear')
    .reproject({ crs: targetProjection, scale: targetScale });

  // Apply the RF model to obtain the predicted ET surface.
  var et_predicted_rf = predictors.classify(rfRegressor).rename('ET_pred');

  // Compute residuals: positive values = model under-predicts ET,
  // negative values = model over-predicts ET.
  var residuals = et_resampled.subtract(et_predicted_rf).rename('residuals');

  // Attach the date string and residual image as properties for use in export.
  return et_predicted_rf
    .set('system:time_start', et.get('system:time_start'))
    .set('date_string',       et.date().format('YYYY_MM_dd'))
    .set('residuals_image',   residuals);
};

// Map the prediction function over all successfully matched triplets.
var predictions    = successfulMatches.map(createPrediction);
var predictionsList = predictions.toList(predictions.size());
var numImages      = predictionsList.size();

print('Number of prediction images to export:', numImages);


// 6. BATCH EXPORT OF PREDICTIONS AND RESIDUALS

/**
 * exportSinglePrediction
 * ----------------------
 * For a given list index, exports:
 *   1. The predicted ET raster to Google Drive as a GeoTIFF (EPSG:32633, 100 m).
 *   2. A CSV of sampled residual points (up to 10 000 per image) including
 *      pixel coordinates and the residual value, tagged with the image date.
 *      These point tables are the direct input for external Kriging workflows.
 *
 * Note: The .evaluate() call is required because the date string is a
 * server-side object; it must be pulled to the client before being used
 * in export descriptions (which are plain JS strings).
 *
 * @param {number} index - Zero-based index into predictionsList.
 */
var exportSinglePrediction = function(index) {
  var prediction   = ee.Image(predictionsList.get(index));
  var residualsImg = ee.Image(prediction.get('residuals_image'));
  var dateString   = prediction.get('date_string');

  // Evaluate brings the server-side date string to the client so it can be
  // used as a plain JS string in export task descriptions and file names.
  dateString.evaluate(function(dateStr) {

    // Export 1: Predicted ET raster 
    Export.image.toDrive({
      image: prediction.toFloat().clip(areaOfInterest.geometry()),
      description: 'pred_rf_' + dateStr,
      folder: 'GEE_Output',
      fileNamePrefix: 'pred_rf_' + dateStr,
      region: areaOfInterest.geometry(),
      crs: 'EPSG:32633',  // UTM Zone 33N — appropriate for the Gargano area
      scale: 100,
      maxPixels: 1e9
    });

    // Export 2: Residual point samples (for Kriging)
    // Randomly sample up to 10 000 pixels from the residual image within the AOI.
    // Geometries are included so that spatial coordinates are preserved in the CSV.
    var residualPoints = residualsImg.sample({
      region: areaOfInterest.geometry(),
      scale: 100,
      numPixels: 10000,
      projection: 'EPSG:32633',
      geometries: true
    }).map(function(f) {
      // Tag each point with its acquisition date for traceability.
      return f.set('date', dateStr);
    });

    Export.table.toDrive({
      collection: residualPoints,
      description: 'residuals_' + dateStr,
      folder: 'GEE_Output',
      fileNamePrefix: 'residuals_' + dateStr,
      fileFormat: 'CSV'
    });
  });
};

// Evaluate the total image count on the client side, then loop over all images
// to schedule individual export tasks. Tasks appear in the GEE Tasks panel.
numImages.evaluate(function(count) {
  print('Scheduling exports for', count, 'images...');
  for (var i = 0; i < count; i++) {
    exportSinglePrediction(i);
  }
  print('All export tasks submitted. Check the Tasks tab to monitor progress.');
});


// Optionally display the AOI boundary (hidden by default — toggle in the Layers panel).
Map.addLayer(areaOfInterest, { color: 'red' }, 'Area of Interest', false);
