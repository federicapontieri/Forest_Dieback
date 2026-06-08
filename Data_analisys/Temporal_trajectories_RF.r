  # ==============================================================================
  # SCRIPT: Random Forest Modelling of NDVI Trend Drivers
  # PROJECT: Gargano Forest Dieback Monitoring
  #
  # DESCRIPTION:
  #   This script performs a complete spatial analysis pipeline to identify and
  #   model the environmental drivers of NDVI decline in the Gargano forest.
  #   It focuses on two forest types: Holm oak (Quercus ilex) and Aleppo pine
  #   (Pinus halepensis) mapped at high Tree Cover Density (TCD ≥ 80 %).
  #
  #   The pipeline covers:
  #     1. Output folder setup
  #     2. Reference raster definition (100 m, EPSG:32633)
  #     3. Landsat 8 processing: NDVI computation and land cover masking
  #     4. Pixel-wise NDVI temporal trend analysis (slope + significance)
  #     5. LST summary metrics: skewness, heat stress (≥ 35 °C), optimum (20–28 °C)
  #     6. ET summary metrics: mean and skewness
  #     7. Multicollinearity check (VIF)
  #     8. Random Forest model training with 10-fold cross-validation (caret + ranger)
  #     9. Model validation plot (observed vs predicted)
  #    10. Variable importance chart
  #    11. Partial Dependence Plots (PDP) per land-use class via iml
  #
  # INPUTS (remembre to replace it with your paths):
  #   - Landsat 8 multi-band rasters: Images/Landsat8_img_100m_pre/*.tif
  #   - Hybrid ET rasters (100 m): Images/Hybrid_Analysis/Hybrid_Output_Final/ET_100m_FINAL/*.tif
  #   - Land cover mask (TCD ≥ 80 %): DATA/land_cover_final_TCD80.tiff
  #   - Digital Elevation Model: DEM/DEM_Gargano.tif
  #
  # OUTPUTS (all written under High_TCD/Outputs_hightcd_NDVI_t1/):
  #   - NDVI rasters (full + masked): NDVI/ and NDVI/CUT/
  #   - NDVI trend maps (slope + p<0.05): NDVI/TREND/
  #   - LST summary rasters: LST/ and LST/CUT/
  #   - ET summary rasters: ET/  and ET/CUT/
  #   - Variogram PNGs (from Part 2): (not produced here)
  #   - Model validation scatter plot: Plots/
  #   - Variable importance chart: Plots/
  #   - PDP plots per variable (.png/.svg): Plots/
  #   - Trained RF model (.rds): Models/
  #   - Training data snapshot (.rds): Models/
  #
  # NOTES:
  #   - All rasters are co-registered to the 100 m LST reference grid (EPSG:32633), you can use the
  #     the CRS that best suits your study area;
  #   - Trend analysis uses linear regression (lm);
  # ==============================================================================
  
  
  #### 1. PACKAGE LOADING ####
  
  # terra      — raster read/write, spatial operations, terrain analysis
  # pbapply    — progress-bar lapply for batch loops
  # parallel   — multi-core processing (makeCluster, clusterExport, etc.)
  # ggplot2    — all plots (trends, validation, variable importance, PDP)
  # ggtext     — element_markdown() for italic species labels in PDP legends
  # ranger     — fast Random Forest engine (called via caret method = "ranger")
  # caret      — unified model training, cross-validation, and variable importance
  # doParallel — parallel backend for caret's cross-validation
  # usdm       — Variance Inflation Factor (VIF) collinearity screening
  # iml        — model-agnostic interpretability: Predictor + FeatureEffect (PDP)
  # e1071      — skewness() for LST and ET summary statistics
  # svglite    — SVG graphics device for ggsave() vector output
  
  library(terra)
  library(pbapply)
  library(parallel)
  library(ggplot2)
  library(ggtext)
  library(ranger)
  library(caret)
  library(doParallel)
  library(usdm)
  library(iml)
  library(e1071)
  library(svglite)
  
  

  #### 2. WORKING DIRECTORY & OUTPUT FOLDER SETUP ####

  # Set the working directory, replace it with your own
  
  setwd("~/Scrivania/DroughtWoods/Gargano")
  
  # Output folder for this analysis run (High TCD, NDVI, during-disturbance period T1).
  output_base <- "High_TCD/Outputs_hightcd_NDVI_t1"
  
  # Sub-folders for each index and its masked version.
  subfolders <- c(
    "MASK",
    "NDVI", "NDVI/CUT", "NDVI/TREND",
    "LST",  "LST/CUT",
    "ET",   "ET/CUT"
  )
  
  # Create all output directories in one loop (no error if already present).
  for (folder in file.path(output_base, subfolders)) {
    if (!dir.exists(folder)) dir.create(folder, recursive = TRUE)
  }
  
  # Additional top-level folders for diagnostic plots and saved models.
  plots_folder  <- file.path(output_base, "Plots")
  models_folder <- file.path(output_base, "Models")
  if (!dir.exists(plots_folder))  dir.create(plots_folder,  recursive = TRUE)
  if (!dir.exists(models_folder)) dir.create(models_folder, recursive = TRUE)
  
  

  #### 3. REFERENCE RASTER DEFINITION ####

  # All input layers are resampled to a common spatial reference defined here:
  #   - CRS: EPSG:32633 (UTM Zone 33N)
  #   - Resolution: ~100 m
  #   - Extent: matched to the 16th Landsat 8 scene (band 8 = LST band)
  #
  # Using a single reference raster ensures pixel-perfect alignment across all layers.
  
  lst_folder <- "./Images/Landsat8_img_100m_pre/"
  lst_files  <- list.files(lst_folder, pattern = "\\.tif$", full.names = TRUE)
  
  if (length(lst_files) == 0) {
    stop("No LST rasters found in Images/Landsat8_img_100m_pre/")
  }
  
  # Load band 8 (LST) from the 16th file as the spatial template.
  reference_raster <- rast(lst_files[15])
  reference_raster <- reference_raster[[8]]
  
  # Reproject to EPSG:32633 if not already in that CRS.
  if (crs(reference_raster) != "EPSG:32633") {
    reference_raster <- project(reference_raster, "EPSG:32633")
  }
  
  

  #### 4. LANDSAT 8 PROCESSING, NDVI COMPUTATION ####

  process_sentinel_file <- function(file_path, ref_raster, vector_mask_path, output_base) {
    
    # Load and resample to 100 m reference grid (bilinear for continuous data).
    r <- rast(file_path)
    r <- resample(r, ref_raster, method = "bilinear")
    
    # Extract spectral bands.
    # Expected band order (check the raster before to assess the order of the bands):
    #                      B1(Aero) B2(Blue) B3(Green) B4(Red) B5(NIR)
    #                      B6(SWIR1) B7(SWIR2) ST_B10(LST) QA_PIXEL
    B4 <- r[[4]]   # Red
    B5 <- r[[5]]   # Near-Infrared (NIR)
    
    # Compute NDVI, mask physically impossible values.
    NDVI <- (B5 - B4) / (B5 + B4)
    NDVI <- ifel(NDVI >= -1 & NDVI <= 1, NDVI, NA)
    
    # Save full-area NDVI.
    ndvi_name <- file.path(output_base, "NDVI", paste0("NDVI_", basename(file_path)))
    writeRaster(NDVI, ndvi_name, overwrite = TRUE)
    
    # Apply the raster land cover mask (TCD ≥ 80 %): crop to extent, then mask.
    land_use_eroded <- rast("DATA/land_cover_final_TCD80.tiff")
    NDVI_cut <- mask(crop(NDVI, land_use_eroded), land_use_eroded)
    ndvi_cut_name <- file.path(output_base, "NDVI", "CUT", paste0("NDVI_CUT_", basename(file_path)))
    writeRaster(NDVI_cut, ndvi_cut_name, overwrite = TRUE)
    
    return(basename(file_path))
  }
  
  # Parallel execution of Landsat 8 processing
  sentinel_folder <- "Images/Landsat8_img_100m_pre/"
  sentinel_files  <- list.files(sentinel_folder, pattern = "\\.tif$", full.names = TRUE)
  
  # Reserve 2 cores for the OS; each worker handles one image independently.
  num_cores <- detectCores() - 2
  cl <- makeCluster(num_cores)
  clusterEvalQ(cl, { library(terra) })
  clusterExport(cl,
                varlist = c("process_sentinel_file", "reference_raster", "output_base"),
                envir   = environment())
  
  results <- pblapply(
    sentinel_files,
    function(f) process_sentinel_file(f, reference_raster, NULL, output_base)
  )
  stopCluster(cl)
  
  

  #### 5. PIXEL-WISE NDVI TEMPORAL TREND ANALYSIS ####

  
  # calculate_significant_slope
  #
  # Fits a linear model to the time series of a single pixel and returns the
  # standardised slope and a significance flag.
  #
  # A linear regression is used (lm). 
  # The independent variable is scaled (standardised) so that slopes 
  # are directly comparable across pixels
  # regardless of the time span.
  
  calculate_significant_slope <- function(pv) {
    time_indices <- as.Date(colnames(pv))
    valid_values <- pv[!is.na(pv)]
    valid_indices <- time_indices[!is.na(pv)]
    
    # Express time as days elapsed since the first valid observation.
    sel          <- valid_indices[1]
    valid_indices <- sapply(valid_indices, function(x) x - sel)
    
    # Return NA/0 if fewer than 2 valid observations or if the series is constant
    # (a constant series has zero variance and cannot be modelled).
    if (length(valid_values) < 2 || length(unique(valid_values)) == 1) {
      return(list(slope = NA, significant = 0))
    }
    
    DATA <- data.frame(valid_indices, valid_values)
    
    # Linear model 
    # scale() standardises valid_indices so the slope is in units of NDVI per SD.
    trend     <- lm(valid_values ~ as.numeric(scale(valid_indices)), data = DATA)
    slope     <- coef(trend)[2]
    p_value   <- summary(trend)$coefficients[2, 4]
    significant <- ifelse(p_value < 0.05, 1, 0)
    
    return(list(slope = slope, significant = significant))
  }
  
  # extract_date
  #
  # Parses the acquisition date from a filename following the convention:
  # <prefix>_YYYY-MM-DD.tif
  # @param filename  Character. File basename or full path.
  # @return          Date object.
  
  extract_date <- function(filename) {
    date_str <- sub(".*_(\\d{4}-\\d{2}-\\d{2})\\.tif$", "\\1", filename)
    as.Date(date_str, format = "%Y-%m-%d")
  }
  
  # plot_raster
  #
  # Creates a ggplot2 map of a single-band SpatRaster using a diverging colour
  # scale centred on zero (purple = negative trend, orange = positive trend).
  # The plot is saved as a 300 dpi PNG to the Plots folder.
  #
  # @param raster      SpatRaster. Single-band raster to map.
  # @param plot_title  Character. Title string for the plot.
  # @param out_name    Character. Output filename stem (no extension).
  
  plot_raster <- function(raster, plot_title, out_name) {
    df <- as.data.frame(raster, xy = TRUE)
    colnames(df)[3] <- "trend"
    
    p <- ggplot(df, aes(x = x, y = y, fill = trend)) +
      geom_raster() +
      scale_fill_gradient2(
        low      = "#8f00ff",
        mid      = "#AEAEAE",
        high     = "#DB7E17",
        midpoint = 0,
        na.value = "white",
        name     = "Trend",
        limits   = c(min(df$trend, na.rm = TRUE), max(df$trend, na.rm = TRUE))
      ) +
      coord_fixed() +
      labs(title = plot_title, x = "Longitude", y = "Latitude") +
      theme_bw() +
      theme(
        plot.background  = element_rect(fill = "white", color = NA),
        plot.title       = element_text(size = 25, hjust = 0.5),
        axis.title.x     = element_text(size = 20),
        axis.title.y     = element_text(size = 20),
        axis.text.x      = element_text(size = 17, color = "black"),
        axis.text.y      = element_text(size = 17, color = "black"),
        panel.grid.minor = element_blank()
      )
    
    ggsave(file.path(plots_folder, paste0(out_name, ".png")), p, width = 8, height = 6)
  }
  
  # process_raster_trend_parallel
  #
  # Loads all NDVI rasters in a folder, fits calculate_significant_slope() to
  # each pixel's time series in parallel, and writes two output rasters:
  #   - slope_values.tif: standardised linear slope per pixel
  #   - significance_map.tif: 1 where p < 0.05, 0 elsewhere
  #
  # @param raster_folder  Character. Folder containing time-stamped NDVI .tif files.
  # @param output_folder  Character. Folder where output rasters are saved.
  # @param plot_title     Character. Title prefix for the diagnostic plots.
  # @param num_cores      Integer. Number of parallel workers (default: nCores - 1).
  # @return               List of two output file paths (slope, significance).
  
  process_raster_trend_parallel <- function(
      raster_folder,
      output_folder,
      plot_title,
      num_cores = detectCores() - 1
  ) {
    # Load and sort files chronologically.
    r_files <- list.files(raster_folder, pattern = "\\.tif$", full.names = TRUE)
    r_files <- r_files[order(sapply(r_files, extract_date))]
    
    raster_list     <- lapply(r_files, rast)
    raster_data_list <- lapply(raster_list, function(r) as.data.frame(r, xy = TRUE, na.rm = FALSE))
    
    # Merge all time steps into one wide data.frame (rows = pixels, columns = dates).
    merged_df    <- Reduce(function(x, y) merge(x, y, by = c("x", "y"), all = TRUE), raster_data_list)
    pixel_values <- merged_df[, -(1:2)]
    
    # Name columns with the acquisition date string (used inside calculate_significant_slope).
    col_names <- gsub(".tif", "", basename(r_files))
    col_names <- gsub("NDVI_CUT_Landsat8_100_", "", col_names)
    colnames(pixel_values) <- col_names
    
    # Run pixel-wise trend computation in parallel.
    cl <- makeCluster(num_cores)
    clusterExport(cl, c("pixel_values", "calculate_significant_slope"), envir = environment())
    # mgcv is loaded on workers in case the GAM branch is reactivated.
    clusterEvalQ(cl, { library(mgcv) })
    
    slope_results <- pblapply(seq_len(nrow(pixel_values)), function(i) {
      calculate_significant_slope(pv = pixel_values[i, ])
    }, cl = cl)
    
    stopCluster(cl)
    
    # Extract slopes and significance flags into separate vectors.
    slopes       <- sapply(slope_results, function(res) res$slope)
    significance <- sapply(slope_results, function(res) res$significant)
    
    # Reconstruct SpatRasters from the result vectors.
    slope_df        <- cbind(merged_df[, 1:2], slope       = slopes)
    significance_df <- cbind(merged_df[, 1:2], significant = significance)
    
    slope_raster        <- rast(slope_df,        type = "xyz", crs = crs(raster_list[[1]]))
    significance_raster <- rast(significance_df, type = "xyz", crs = crs(raster_list[[1]]))
    
    if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)
    
    slope_file        <- file.path(output_folder, "slope_values.tif")
    significance_file <- file.path(output_folder, "significance_map.tif")
    writeRaster(slope_raster,        slope_file,        overwrite = TRUE)
    writeRaster(significance_raster, significance_file, overwrite = TRUE)
    
    plot_raster(slope_raster,        paste0(plot_title, " - Slopes"),       paste0(plot_title, "_Slopes"))
    plot_raster(significance_raster, paste0(plot_title, " - Significance"), paste0(plot_title, "_Significance"))
    
    return(list(slope_file, significance_file))
  }
  
  # Run the NDVI trend analysis.
  ndvi_result <- process_raster_trend_parallel(
    file.path(output_base, "NDVI", "CUT"),
    file.path(output_base, "NDVI", "TREND"),
    "NDVI Trend"
  )

  #### 6. LST SUMMARY METRICS ####

  # Load all LST images, reproject to EPSG:32633, and resample to the 100 m grid.
  lst_list <- lapply(lst_files, function(file) {
    r <- rast(file)
    if (crs(r) != "EPSG:32633") {
      r <- project(r, "EPSG:32633")
      r <- r[[8]]  # Band 8 = ST_B10 (surface temperature)
    }
    resample(r, reference_raster, method = "bilinear")
  })
  
  lst_stack <- rast(lst_list)
  
  # Note: digital-number to Celsius conversion is commented out below.
  # Uncomment and adjust if raw DN values have not yet been converted upstream.
  
  # lst_stack <- (lst_stack * 0.00341802) + 149.0 - 273.15
  
  # Skewness of the LST time series per pixel 
  # Skewness captures asymmetry in temperature distribution
  
  lst_sk <- app(lst_stack, fun = e1071::skewness, na.rm = TRUE)
  writeRaster(lst_sk, file.path(output_base, "LST", "LST_SK.tif"), overwrite = TRUE)
  
  cnat_mask  <- rast("DATA/land_cover_final_TCD80.tiff")
  lst_sk_cut <- mask(crop(lst_sk, cnat_mask), cnat_mask)
  writeRaster(lst_sk_cut, file.path(output_base, "LST", "CUT", "LST_SK_CUT.tif"), overwrite = TRUE)
  
  # Percentage of days exceeding 35 °C (heat stress threshold, see manuscript for bibliographic references) 

  exceed_binary_stack <- ifel(lst_stack >= 35, 1, 0)
  exceed_count        <- sum(exceed_binary_stack, na.rm = TRUE)
  
  valid_binary_stack <- ifel(is.na(lst_stack), 0, 1)
  valid_count        <- sum(valid_binary_stack, na.rm = TRUE)
  
  exceed_percentage <- ifel(valid_count == 0, NA, (exceed_count / valid_count) * 100)
  
  writeRaster(exceed_percentage,
              file.path(output_base, "LST", "LST_exceed_35C_percent.tif"),
              overwrite = TRUE)
  
  exceed_percentage_cut <- mask(crop(exceed_percentage, cnat_mask), cnat_mask)
  writeRaster(exceed_percentage_cut,
              file.path(output_base, "LST", "CUT", "LST_exceed_35C_percent_cut.tif"),
              overwrite = TRUE)
  
  # Percentage of days within the photosynthetic optimum (20–28 °C, check manuscript for bibliographic references) 
  
  optimum_binary_stack <- ifel(lst_stack >= 20 & lst_stack <= 28, 1, 0)
  optimum_count        <- sum(optimum_binary_stack, na.rm = TRUE)
  
  optimum_percentage <- ifel(valid_count == 0, NA, (optimum_count / valid_count) * 100)
  
  writeRaster(optimum_percentage,
              file.path(output_base, "LST", "LST_optimum_percent.tif"),
              overwrite = TRUE)
  
  optimum_percentage_cut <- mask(crop(optimum_percentage, cnat_mask), cnat_mask)
  writeRaster(optimum_percentage_cut,
              file.path(output_base, "LST", "CUT", "LST_optimum_percent_cut.tif"),
              overwrite = TRUE)
  
  #### 7. ET SUMMARY METRICS ####
  
  # Load the final hybrid ET rasters (100 m) produced by the Kriging script (Part 2).
  et_folder <- "./Images/Hybrid_Analysis/Hybrid_Output_Final/ET_100m_FINAL"
  et_files  <- list.files(et_folder, pattern = "\\.tif$", full.names = TRUE)
  
  if (length(et_files) == 0) {
    stop("No ET rasters found in Images/Hybrid_Analysis/Hybrid_Output_Final/ET_100m_FINAL/")
  }
  
  et_list <- lapply(et_files, function(file) {
    r <- rast(file)
    if (crs(r) != "EPSG:32633") r <- project(r, "EPSG:32633")
    return(r)
  })
  
  et_stack <- rast(et_list)
  
  # Mean ET across the time series 
  et_mean <- app(et_stack, fun = mean, na.rm = TRUE)
  writeRaster(et_mean, file.path(output_base, "ET", "ET_MEAN.tif"), overwrite = TRUE)
  
  # Skewness of the ET time series
  et_sk <- app(et_stack, fun = e1071::skewness, na.rm = TRUE)
  writeRaster(et_sk, file.path(output_base, "ET", "ET_SK.tif"), overwrite = TRUE)
  
  # Apply land cover mask and save clipped versions.
  cnat_mask   <- rast("DATA/land_cover_final_TCD80.tiff")
  et_mean_cut <- mask(crop(et_mean, cnat_mask), cnat_mask)
  et_sk_cut   <- mask(crop(et_sk,   cnat_mask), cnat_mask)
  
  writeRaster(et_mean_cut, file.path(output_base, "ET", "CUT", "ET_MEAN_CUT.tif"), overwrite = TRUE)
  writeRaster(et_sk_cut,   file.path(output_base, "ET", "CUT", "ET_SK_CUT.tif"),   overwrite = TRUE)
  
  #### 8. PREDICTOR STACK ASSEMBLY ####

  # Define paths for all pre-computed predictor layers, replace them with your paths.
  ndvi_trend_path  <- file.path(output_base, "NDVI", "TREND", "slope_values.tif")
  significance_path <- file.path(output_base, "NDVI", "TREND", "significance_map.tif")
  lst_sk_path      <- file.path(output_base, "LST",  "CUT",   "LST_SK_CUT.tif")
  lst_exceed_path  <- file.path(output_base, "LST",  "CUT",   "LST_exceed_35C_percent_cut.tif")
  lst_optimum_path <- file.path(output_base, "LST",  "CUT",   "LST_optimum_percent_cut.tif")
  et_sk_path       <- file.path(output_base, "ET",   "CUT",   "ET_SK_CUT.tif")
  et_mean_path     <- file.path(output_base, "ET",   "CUT",   "ET_MEAN_CUT.tif")
  dem_path         <- "DEM/DEM_Gargano.tif"
  
  # Load each layer and resample to the 100 m reference grid.
  # Categorical layers use nearest-neighbour, while continuous layers use bilinear.
  
  ndvi_trend <- rast(ndvi_trend_path)
  if (is.na(crs(ndvi_trend))) crs(ndvi_trend) <- "EPSG:32633"
  ndvi_trend <- resample(ndvi_trend, reference_raster, method = "bilinear")
  
  dem <- rast(dem_path)
  if (crs(dem) != "EPSG:32633") dem <- project(dem, "EPSG:32633")
  dem <- resample(dem, reference_raster, method = "bilinear")
  
  # Land cover (categorical: 10 = Holm oak, 20 = Aleppo pine) -> nearest-neighbour.
  land_use <- rast("DATA/land_cover_final_TCD80.tiff")
  land_use <- resample(land_use, reference_raster, method = "near")
  
  significance <- rast(significance_path)
  significance <- resample(significance, reference_raster, method = "near")
  
  lst_sk_raster  <- resample(rast(lst_sk_path),      reference_raster, method = "bilinear")
  lst_exceed     <- resample(rast(lst_exceed_path),   reference_raster, method = "bilinear")
  lst_optimum    <- resample(rast(lst_optimum_path),  reference_raster, method = "bilinear")
  et_sk_raster   <- resample(rast(et_sk_path),        reference_raster, method = "bilinear")
  et_mean_raster <- resample(rast(et_mean_path),      reference_raster, method = "bilinear")
  
  # Derive morphometric variables from the DEM.
  slope     <- terrain(dem, v = "slope",     unit = "degrees")
  aspect    <- terrain(dem, v = "aspect",    unit = "degrees")
  tpi       <- terrain(dem, v = "TPI")
  tri       <- terrain(dem, v = "TRI")
  roughness <- terrain(dem, v = "roughness")
  
  # Stack all layers into a single multi-band SpatRaster.
  raster_stack <- c(
    ndvi_trend,
    dem, slope, aspect, tpi, tri, roughness,
    land_use, significance,
    lst_sk_raster, lst_exceed, lst_optimum,
    et_sk_raster, et_mean_raster
  )
  names(raster_stack) <- c(
    "NDVI_Trend",
    "DEM", "Slope", "Aspect", "TPI", "TRI", "Roughness",
    "Land_Use", "Significance",
    "LST_SK", "LST_Exceed", "LST_Optimum",
    "ET_SK", "ET_MEAN"
  )
  
  
  # Remove pixels that have NA in ANY layer (complete-case filtering).
 
  na_count           <- app(raster_stack, fun = function(x) sum(is.na(x)))
  complete_case_mask <- ifel(na_count == 0, 1, NA)
  raster_stack_clean <- mask(raster_stack, complete_case_mask)
  
  # Convert the clean stack to a data.frame for modelling (one row per pixel).
  data_matrix <- as.data.frame(raster_stack_clean, na.rm = TRUE)
  
  #### 9. VARIANCE INFLATION FACTOR (VIF) SCREENING ####

  # VIF > 5 indicates problematic multicollinearity. 

  vif_matrix <- data_matrix
  vif_matrix$NDVI_Trend   <- NULL
  vif_matrix$Land_Use     <- NULL
  vif_matrix$Significance <- NULL
  
  vif_results <- vifstep(vif_matrix, th = 5)
  print(vif_results)
  
  #### 10. RANDOM FOREST MODEL TRAINING ####
  
  # Prepare modelling data 
  # Keep only pixels with a statistically significant NDVI trend (p < 0.05).
  # These are the pixels where a real signal exists to model.
  
  model_data <- data_matrix[data_matrix$Significance == 1, ]
  model_data$Significance <- NULL  # No longer needed as a column
  
  # Remove variables identified as collinear by VIF analysis (TRI, Roughness in our specific case, that depends on the VIF results).
  model_data$TRI       <- NULL
  model_data$Roughness <- NULL
  
  # Land_Use must be a factor so caret/ranger treats it as categorical.
  model_data$Land_Use <- as.factor(model_data$Land_Use)
  
  # Custom summary function for caret (adds SMAPE to default RMSE + R²) 
  #
  # SMAPE (Symmetric Mean Absolute Percentage Error) is included because it is
  # scale-independent and handles both over- and under-predictions symmetrically,
  # complementing RMSE for comparing performance across folds.
  custom_summary <- function(data, lev = NULL, model = NULL) {
    rmse      <- sqrt(mean((data$obs - data$pred)^2))
    smape     <- mean(2 * abs(data$obs - data$pred) /
                        (abs(data$obs) + abs(data$pred))) * 100
    r_squared <- cor(data$obs, data$pred)^2
    return(c(RMSE = rmse, SMAPE = smape, Rsquared = r_squared))
  }
  
  # Set up parallel backend for caret cross-validation
  num_cores <- detectCores() - 1
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  
  # 10-fold cross-validation configuration
  train_control <- trainControl(
    method          = "cv",
    number          = 10,
    verboseIter     = TRUE,
    allowParallel   = TRUE,
    summaryFunction = custom_summary,
    savePredictions = "final"   # Save fold-level predictions for the scatter plot
  )
  
  # Train the Random Forest
  
  set.seed(666)
  rf_model_cv <- train(
    NDVI_Trend ~ .,
    data       = model_data,
    method     = "ranger",
    trControl  = train_control,
    tuneLength = 5,
    importance = "permutation"  # Permutation importance: unbiased for mixed variable types
  )
  
  print(rf_model_cv)
  
  stopCluster(cl)
  registerDoSEQ()
  
  # Save the trained model and training data snapshot for PDP analysis.
  saveRDS(rf_model_cv, file = file.path(models_folder, "rf_model_CDNedit_LST_100_hightcd_t1.rds"))
  saveRDS(model_data,  file = file.path(models_folder, "training_data_for_pdp.rds"))
  
  # Report best performance metrics.
  results <- rf_model_cv$results
  cat("Best SMAPE : ", min(results$SMAPE),     "\n")
  cat("Best R²    : ", max(results$Rsquared),  "\n")
  
  #### 11. MODEL VALIDATION PLOT (OBSERVED vs PREDICTED) ####

  cv_predictions <- rf_model_cv$pred
  
  mae_val  <- round(mean(abs(cv_predictions$obs - cv_predictions$pred)), 4)
  bias_val <- round(mean(cv_predictions$pred - cv_predictions$obs),      5)
  r2_val   <- round(cor(cv_predictions$obs, cv_predictions$pred)^2,      3)
  
  stats_text <- paste0("R² = ", r2_val, "\nMAE = ", mae_val, "\nBias = ", bias_val)
  
  p_val <- ggplot(cv_predictions, aes(x = obs, y = pred)) +
    geom_bin2d(bins = 100) +
    scale_fill_gradient(low = "#e0f2f1", high = "#004d40", name = "Pixel Count") +
    geom_abline(slope = 1, intercept = 0, color = "#8c301b",
                linetype = "dashed", linewidth = 0.9) +
    geom_smooth(method = "lm", color = "gray30", linewidth = 1) +
    annotate("label",
             x = min(cv_predictions$obs), y = max(cv_predictions$pred),
             label = stats_text, hjust = 0, vjust = 1,
             size = 5, fontface = "bold", alpha = 0.8) +
    labs(
      title    = "Model Validation: Stress Phase (T1)",
      subtitle = "Cross-Validation Predictions (10-fold)",
      x        = "Observed NDVI Trend",
      y        = "Predicted NDVI Trend"
    ) +
    theme_bw() +
    theme(
      plot.title    = element_text(hjust = 0.5, size = 20, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 12),
      axis.title    = element_text(size = 15),
      axis.text     = element_text(size = 13, color = "black"),
      legend.text   = element_text(size = 15),
      legend.box.background = element_rect(color = "black", linewidth = 0.3),
      legend.box.margin     = margin(2, 2, 2, 2)
    )
  
  ggsave(file.path(plots_folder, "Model_Validation_Obs_vs_Pred_T1.png"),
         p_val, width = 8, height = 7, dpi = 300)
  
  #### 12. VARIABLE IMPORTANCE PLOT ####

  vi_raw <- varImp(rf_model_cv)
  vi_df  <- data.frame(
    Variable   = rownames(vi_raw$importance),
    Importance = vi_raw$importance$Overall
  )
  vi_df <- vi_df[order(vi_df$Importance, decreasing = TRUE), ]
  
  # Identify the top 7 variables for use in the PDP section.
  top_vars <- vi_df$Variable[1:7]
  
  # Aggregate all Land_Use dummy variables into a single "Land_Use" entry.
  # caret expands factors into dummies internally; summing their importances
  # gives the aggregate contribution of land cover class.
  vi_df_plot      <- vi_df
  land_use_rows   <- grepl("^Land_Use", vi_df_plot$Variable)
  if (any(land_use_rows)) {
    land_use_imp  <- sum(vi_df_plot$Importance[land_use_rows])
    vi_df_plot    <- vi_df_plot[!land_use_rows, ]
    vi_df_plot    <- rbind(vi_df_plot, data.frame(Variable = "Land_Use", Importance = land_use_imp))
  }
  
  # Convert to percentage of total importance for interpretability.
  vi_df_plot$ImportancePercent <- (vi_df_plot$Importance / sum(vi_df$Importance)) * 100
  
  max_limit <- 60
  vi_plot <- ggplot(vi_df_plot,
                    aes(x = reorder(gsub("_", " ", Variable), -Importance),
                        y = ImportancePercent)) +
    geom_bar(stat = "identity", fill = "#DB7E17", color = "black") +
    coord_flip() +
    scale_y_continuous(
      limits = c(0, max_limit),
      breaks = seq(0, max_limit, by = 20),
      expand = expansion(mult = c(0, 0.05))
    ) +
    labs(
      title = "Variable Importance (%)",
      x     = "Variable",
      y     = "Importance (%)"
    ) +
    theme_bw(base_size = 14) +
    theme(
      plot.background  = element_rect(fill = "white", color = NA),
      plot.title       = element_text(size = 25, hjust = 0.5),
      axis.title.x     = element_text(size = 20),
      axis.title.y     = element_text(size = 20),
      axis.text.x      = element_text(size = 17, color = "black"),
      axis.text.y      = element_text(size = 17, color = "black"),
      panel.grid.minor = element_blank()
    )
  
  ggsave(file.path(plots_folder, "Variable_Importance.tiff"),
         vi_plot, width = 8, height = 6)
  
  #### 13. PARTIAL DEPENDENCE PLOTS (PDP) VIA iml ####

  # Reload the saved model and training data if running this section independently.
  rf_model  <- readRDS(file.path(models_folder, "rf_model_CDNedit_LST_100_hightcd_t1.rds"))
  # model_data <- readRDS(file.path(models_folder, "training_data_for_pdp.rds"))
  
  # Sub-sample 60 % of the training data to reduce PDP computation time.
  # Increasing this fraction improves precision at the cost of runtime.
  set.seed(666)
  sampled_data <- model_data[sample(nrow(model_data), size = 0.60 * nrow(model_data)), ]
  
  # Ensure factor levels in the sample match those seen during training.
  training_data <- rf_model$trainingData
  for (var in colnames(training_data)) {
    if (var %in% colnames(sampled_data)) {
      train_class <- class(training_data[[var]])
      if (train_class == "factor") {
        sampled_data[[var]] <- as.factor(sampled_data[[var]])
        levels(sampled_data[[var]]) <- levels(training_data[[var]])
      } else if (train_class %in% c("numeric", "integer")) {
        sampled_data[[var]] <- as.numeric(sampled_data[[var]])
      }
    }
  }
  
  # Custom prediction function required by iml's Predictor object.
  predict_function <- function(model, newdata) predict(model, newdata = newdata, type = "raw")
  
  # Build the iml Predictor wrapper around the caret model.
  predictor <- Predictor$new(
    model       = rf_model,
    data        = sampled_data,
    y           = sampled_data$NDVI_Trend,
    predict.fun = predict_function
  )
  
  # Generate PDP for each of the top 7 variables, split by land-use class
  # Each plot shows how predicted NDVI Trend responds to one predictor,
  # separately for Holm oak (Qi = class 10) and Aleppo pine (Ph = class 20).
  # The raw per-sample lines (dashed) are overlaid with a LOESS smoother.
  for (var in top_vars) {
    if (!var %in% colnames(sampled_data)) {
      warning(paste("Variable", var, "not found in sampled data — skipping PDP."))
      next
    }
    
    feature_effect <- FeatureEffect$new(
      predictor,
      feature = c(var, "Land_Use"),
      method  = "pdp"
    )
    
    pdp_data <- feature_effect$results
    
    # Map numeric land-use codes to italic species abbreviations for the legend.
    # Qi = Quercus ilex (Holm oak, code 10)
    # Ph = Pinus halepensis (Aleppo pine, code 20)
    pdp_data$Land_Use <- factor(
      pdp_data$Land_Use,
      levels = c(10, 20),
      labels = c("*Qi*", "*Ph*")
    )
    
    pdp_plot <- ggplot(pdp_data, aes(x = .data[[var]], y = .value,
                                     group = Land_Use, color = Land_Use)) +
      geom_line(linetype = "dashed", linewidth = 1, color = "gray30") +
      geom_smooth(method = "loess", se = TRUE) +
      scale_color_manual(values = c(
        "*Qi*" = "#df7f25",   
        "*Ph*" = "#8f09e9"  
      )) +
      labs(
        x     = gsub("_", " ", var),
        y     = "NDVI Trend",
        title = paste("Partial Dependence for", gsub("_", " ", var)),
        color = "Land Use"
      ) +
      theme_bw(base_size = 14) +
      theme(
        plot.title    = element_markdown(size = 25, hjust = 0.5),
        axis.title.x  = element_markdown(size = 20),
        axis.title.y  = element_markdown(size = 20, margin = margin(r = 10)),
        axis.text.x   = element_markdown(size = 20, color = "black"),
        axis.text.y   = element_markdown(size = 20, color = "black"),
        legend.position  = "bottom",
        legend.text      = element_markdown(size = 25),
        legend.title     = element_markdown(size = 25)
      )
    
    # Save as high-resolution TIFF and SVG (vector, for editing).
    tiff(
      filename = file.path(plots_folder, paste0("PDP_", var, ".tiff")),
      width = 30, height = 20, units = "cm", res = 300
    )
    print(pdp_plot)
    dev.off()
    
    ggsave(
      filename = file.path(plots_folder, paste0("PDP_", var, ".svg")),
      plot     = pdp_plot,
      width    = 30, height = 20, units = "cm"
    )
  }
