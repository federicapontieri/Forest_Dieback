# ==============================================================================
# SCRIPT: Kriging Correction of RF-Predicted Evapotranspiration
# PROJECT: Gargano Forest Dieback Monitoring
# DESCRIPTION:
#   This script is PART 2 of a two-step hybrid modelling pipeline for
#   Evapotranspiration (ET) estimation at 100 m resolution over the Gargano
#   forest area.
#
#   It reads the outputs produced by the GEE Random Forest script (Part 1):
#     - Predicted ET rasters  (pred_rf_<YYYY_MM_DD>.tif)
#     - Residual point tables (residuals_<YYYY_MM_DD>.csv)
#   and applies Ordinary Kriging to interpolate the RF residuals spatially.
#   The interpolated residual surface is then added back to the RF prediction
#   to produce a bias-corrected final ET map.
#
#   The hybrid approach (RF + Kriging) is designed to:
#     1. Capture the non-linear relationships between ET and its drivers (RF step)
#     2. Correct the spatially structured error that RF leaves behind (Kriging step)
#
# WORKFLOW POSITION:
#   Requires outputs from Part 1 (GEE Random Forest script) before running.
#
# INPUTS (per date):
#   - pred_rf_<YYYY_MM_DD>.tif: RF-predicted ET raster (GEE output)
#   - residuals_<YYYY_MM_DD>.csv: Sampled residual points with .geo column
#                                      containing GeoJSON coordinates (GEE output)
#   - Original ET raster: Reference ET at 500 m for accuracy validation
#
# OUTPUTS (per date, written to Hybrid_Output_Final/):
#   - ET_final_<date>_RF_only.tif: RF prediction raster (copy)
#   - ET_final_<date>_resid_kriged.tif: Kriged residual surface
#   - ET_100m_FINAL/ET_final_<date>_ET_final.tif: Final bias-corrected ET raster
#   - Variograms/variogram_<date>.png: Fitted variogram plot (300 dpi)
#   - accuracy_results.csv / .xlsx: Per-date accuracy metrics
#                                          (RMSE, MAE, R² for RF-only and RF+Kriging)
#
# NOTES:
#   - Dates with missing files or empty residual CSVs (e.g. cloudy scenes) are
#     skipped automatically and reported in the console.
#   - Kriging is skipped for dates with fewer than 50 residual points; in that
#     case the RF prediction is used as-is for the final ET map.
#   - All rasters are projected to EPSG:32633 (UTM Zone 33N).
# ==============================================================================


#### 1_PACKAGE INSTALLATION & LOADING #####

# Define all required packages. Each is installed automatically if not present.
required_packages <- c(
  "terra",      # Raster read/write and spatial operations
  "gstat",      # Variogram modelling and Kriging interpolation
  "sp",         # Legacy spatial classes required by gstat
  "jsonlite",   # Parsing GeoJSON coordinate strings from the residuals CSV
  "openxlsx",   # Exporting the accuracy summary table as .xlsx
  "pbapply",    # lapply with a progress bar for batch processing
  "ggplot2"     # Variogram plots
)

for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}


#### 2. PATH CONFIGURATION ####

# Set the working directory to the project root.
# All relative paths below are resolved from here.
setwd("C:/DroughtWoods/Gargano/Images/Hybrid_Analysis_Post/") # Remember to replace it with your specific path

# Input paths
# Folder containing GEE outputs: RF prediction rasters and residual CSVs.
input_folder_gee <- "./GEE_Output/"

# Folder containing the original 500 m ET rasters used for accuracy validation.
input_folder_et_original <- "C:/DroughtWoods/Gargano/Images/et_500m_post/" # Remember to replace it with your specific path

# Output paths 
# Root output directory; all results are written under subdirectories here.
output_dir     <- file.path(getwd(), "Hybrid_Output_Final")

# Final bias-corrected ET rasters (100 m, one per date).
et_final_dir   <- file.path(output_dir, "ET_100m_FINAL")

# Variogram diagnostic plots (one PNG per date).
variogram_dir  <- file.path(output_dir, "Variograms")

# Create output directories (silently, no error if they already exist).
dir.create(output_dir,    showWarnings = FALSE, recursive = TRUE)
dir.create(variogram_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(et_final_dir,  showWarnings = FALSE, recursive = TRUE)


#### 3. CORE PROCESSING FUNCTION (PER-DATE)####

# process_single_date
#
# Processes one date of the RF+Kriging hybrid pipeline:
#   1. Loads the RF prediction raster and residual point CSV;
#   2. Parses GeoJSON coordinates from the '.geo' column (GEE export format);
#   3. Fits an empirical variogram to the residuals and selects the best model
#      (Spherical, Exponential, or Gaussian);
#   4. Performs Ordinary Kriging to interpolate residuals onto the RF grid.
#   5. Adds the kriged residuals to the RF prediction → final ET surface;
#   6. Saves intermediate and final rasters plus a variogram plot;
#   7. Computes and returns per-date accuracy metrics against the original ET.
#
# @param pred_rf_file  Character. Full path to a pred_rf_<YYYY_MM_DD>.tif file.
# @return A one-row data.frame with accuracy statistics, or NULL if the date
#         is skipped due to missing/empty input files.
process_single_date <- function(pred_rf_file) {

  #
  # File Existence & Validity Checks

  # Extract the date string (YYYY_MM_DD) from the RF prediction filename.
  date_str_underscore <- sub(
    ".*pred_rf_(\\d{4}_\\d{2}_\\d{2})\\.tif", "\\1",
    basename(pred_rf_file)
  )

  # Build the expected residuals CSV filename from the same date string.
  residuals_csv_file <- file.path(
    input_folder_gee,
    paste0("residuals_", date_str_underscore, ".csv")
  )

  # Search for the matching original ET raster (any file containing the date).
  et_original_file <- list.files(
    input_folder_et_original,
    pattern   = date_str_underscore,
    full.names = TRUE
  )

  # Skip this date if the residuals CSV or original ET raster is missing.
  if (!file.exists(residuals_csv_file) || length(et_original_file) == 0) {
    cat("\nWARNING: Missing files for date", date_str_underscore,
        "— skipping.\n")
    return(NULL)
  }

  # Skip if the residuals CSV is empty (can happen for fully cloud-covered scenes).
  if (file.info(residuals_csv_file)$size == 0) {
    cat("\nWARNING: Residuals file for date", date_str_underscore,
        "is empty (likely due to cloud cover) — skipping.\n")
    return(NULL)
  }

  cat("\nProcessing date:", date_str_underscore, "\n")

  # Load Input Data

  # RF prediction raster (100 m, EPSG:32633, single band).
  pred_rf        <- rast(pred_rf_file)
  names(pred_rf) <- "ET_RF"

  # Original ET raster resampled to match the RF prediction grid (for validation).
  # Using the first match in case multiple files share the date string.
  et_original    <- rast(et_original_file[1])
  et_resampled   <- resample(et_original, pred_rf)
  names(et_resampled) <- "ET_original"

  # Load residual points exported from GEE.
  # The '.geo' column contains GeoJSON strings with WGS84 coordinates.
  residual_points_raw <- read.csv(residuals_csv_file)

  # Parse GeoJSON Coordinates

  # GEE exports geometry as a JSON string in the '.geo' column, e.g.:
  # '{"type":"Point","coordinates":[-15.5, 41.2]}'
  # We extract the [longitude, latitude] pair from each row.
  coords_list <- lapply(residual_points_raw$.geo, function(json_string) {
    return(jsonlite::fromJSON(json_string)$coordinates)
  })
  coords <- do.call(rbind, coords_list)
  colnames(coords) <- c("x", "y")

  # Build a plain data.frame of coordinates + residual values.
  df <- data.frame(
    x         = coords[, 1],
    y         = coords[, 2],
    residuals = residual_points_raw$residuals
  )

  # Kriging Interpolation of Residuals

  # Initialise the kriged residuals raster to zero.
  # This default is used if Kriging is skipped (too few points or flat variogram).
  target_crs        <- "EPSG:32633"
  residuals_kriged  <- rast(pred_rf)
  values(residuals_kriged) <- 0
  names(residuals_kriged) <- "residuals_kriged"

  # Require at least 50 residual points to fit a meaningful variogram.
  if (nrow(df) >= 50) {

    # Convert the data.frame to an sp SpatialPointsDataFrame in WGS84,
    # then reproject to UTM Zone 33N (matching the RF raster CRS).
    sp::coordinates(df) <- ~x + y
    sp::proj4string(df) <- sp::CRS("+proj=longlat +datum=WGS84")
    df_krig <- sp::spTransform(
      df,
      sp::CRS("+proj=utm +zone=33 +datum=WGS84 +units=m +no_defs")
    )

    # Fit an empirical variogram to the residuals (assumes stationarity).
    # gstat tries Spherical, Exponential, and Gaussian models; the best fit
    # is selected automatically based on weighted least squares.
    vario <- variogram(residuals ~ 1, data = df_krig)
    fit   <- fit.variogram(vario, vgm(c("Sph", "Exp", "Gau")))

    # Only proceed with Kriging if the fitted model has a meaningful sill
    # (psill > 1e-9). A near-zero sill indicates no spatial structure.
    if (fit$psill[nrow(fit)] > 1e-9) {

      # Variogram Diagnostic Plot
      
      fitted_line <- variogramLine(fit, maxdist = max(vario$dist))

      variogram_plot <- ggplot(data = vario, aes(x = dist, y = gamma)) +
        geom_point(color = '#0072B2', alpha = 0.8, size = 2.5) +
        geom_line(
          data     = fitted_line,
          aes(x = dist, y = gamma),
          color    = '#D55E00',
          linewidth = 1
        ) +
        labs(
          title    = paste("Variogram of RF Residuals:", date_str_underscore),
          subtitle = paste(
            "Model:", fit$model[2],
            "(psill =", round(fit$psill[2], 3),
            "| range =", round(fit$range[2], 1),
            ") | Nugget =", round(fit$psill[1], 3)
          ),
          x = "Lag Distance (m)",
          y = expression("Semivariance (" * gamma * ")")
        ) +
        theme_bw()

      # Save the variogram plot as a 300 dpi PNG.
      ggsave(
        filename = file.path(variogram_dir, paste0("variogram_", date_str_underscore, ".png")),
        plot     = variogram_plot,
        width    = 8, height = 6, dpi = 300
      )

      # Build Prediction Grid & Run Kriging
      
      # Convert the RF raster to an sp SpatialPixelsDataFrame to serve as the
      # Kriging prediction grid (gstat requires the legacy sp format here).
      grd_df <- as.data.frame(pred_rf, xy = TRUE, na.rm = FALSE)
      coordinates(grd_df) <- ~x + y
      gridded(grd_df)     <- TRUE
      proj4string(grd_df) <- CRS(target_crs)

      # Perform Ordinary Kriging: interpolate residuals onto the prediction grid.
      krige_res <- krige(residuals ~ 1, df_krig, grd_df, fit)

      # Convert the Kriging output back to a terra SpatRaster.
      kriged_rast_attempt <- rast(krige_res["var1.pred"])
      crs(kriged_rast_attempt) <- target_crs

      # Replace the zero-initialised residuals raster only if Kriging succeeded.
      if (terra::hasValues(kriged_rast_attempt)) {
        residuals_kriged <- kriged_rast_attempt
      }
    }
  } else {
    cat("  INFO: Fewer than 50 residual points for", date_str_underscore,
        "— Kriging skipped. Using RF prediction only.\n")
  }

  # Combine RF Prediction + Kriged Residuals → Final ET

  # Final ET = RF prediction + spatially interpolated bias correction.
  # Where Kriging was skipped, the residuals layer is zero so ET_final == ET_RF.
  et_final        <- pred_rf + residuals_kriged
  names(et_final) <- "ET_final"

  
  # Write Output Rasters

  out_name_prefix <- paste0("ET_final_", date_str_underscore)

  # RF-only prediction (intermediate; useful for direct comparison).
  writeRaster(
    pred_rf,
    file.path(output_dir, paste0(out_name_prefix, "_RF_only.tif")),
    overwrite = TRUE
  )

  # Kriged residual surface (intermediate; useful for diagnosing the correction).
  writeRaster(
    residuals_kriged,
    file.path(output_dir, paste0(out_name_prefix, "_resid_kriged.tif")),
    overwrite = TRUE
  )

  # Final bias-corrected ET raster (primary deliverable).
  writeRaster(
    et_final,
    file.path(et_final_dir, paste0(out_name_prefix, "_ET_final.tif")),
    overwrite = TRUE
  )

  # Accuracy Assessment Against Original ET

  # Stack all three layers and extract pixel values for comparison.
  # na.rm = TRUE discards pixels where any layer has NoData.
  validation_df <- as.data.frame(
    c(pred_rf, et_final, et_resampled),
    xy    = TRUE,
    na.rm = TRUE
  )

  # Compute residuals for both the RF-only and the final (RF + Kriging) maps.
  res_rf    <- validation_df$ET_original - validation_df$ET_RF
  res_final <- validation_df$ET_original - validation_df$ET_final

  # Helper: total sum of squares for R² denominator.
  ss_tot <- sum((validation_df$ET_original - mean(validation_df$ET_original))^2)

  # Assemble a one-row accuracy summary for this date.
  stats <- data.frame(
    Date                    = date_str_underscore,
    RMSE_RF_only            = sqrt(mean(res_rf^2)),
    MAE_RF_only             = mean(abs(res_rf)),
    R2_RF_only              = 1 - (sum(res_rf^2)    / ss_tot),
    RMSE_RF_plus_Kriging    = sqrt(mean(res_final^2)),
    MAE_RF_plus_Kriging     = mean(abs(res_final)),
    R2_RF_plus_Kriging      = 1 - (sum(res_final^2) / ss_tot),
    Num_Validation_Points   = nrow(validation_df)
  )

  return(stats)
}


#### 4. BATCH EXECUTION & RESULTS EXPORT ####

# Discover all RF prediction rasters in the GEE output folder.
# Only files matching the naming convention pred_rf_*.tif are selected.
pred_files <- list.files(
  input_folder_gee,
  pattern    = "^pred_rf_.*\\.tif$",
  full.names = TRUE
)

cat("Starting automated processing for", length(pred_files), "dates...\n")

# Apply process_single_date to every prediction file.
# pblapply is a drop-in replacement for lapply that shows a progress bar —
# useful when processing dozens of dates.
all_stats_list <- pblapply(pred_files, process_single_date)

# Combine per-date results into a single data.frame, discarding NULL entries
# (dates that were skipped due to missing or invalid inputs).
final_results_table <- do.call(
  rbind,
  all_stats_list[!sapply(all_stats_list, is.null)]
)

# Save the accuracy summary table if at least one date was processed successfully.
if (!is.null(final_results_table) && nrow(final_results_table) > 0) {

  write.csv(
    final_results_table,
    file.path(output_dir, "accuracy_results.csv"),
    row.names = FALSE
  )

  write.xlsx(
    final_results_table,
    file.path(output_dir, "accuracy_results.xlsx"),
    rowNames = FALSE
  )

  cat("\n\nAccuracy metrics table saved successfully.\n")
  print(final_results_table)

} else {
  cat("\n\nNo dates were processed successfully. Check file names and paths.\n")
}

cat(
  "\nProcessing complete for all dates!\n",
  "Final rasters, variogram plots, and accuracy table saved to:", output_dir, "\n"
)
