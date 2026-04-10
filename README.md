# bathymetric-filtering-and-flowline-analysis
MATLAB workflow for bathymetric analysis combining spectral methods, multi-scale filtering, and flow routing (TopoToolbox). Identifies dominant wavelengths and extracts scale-dependent flow networks from DEMs for geomorphic analysis.


# Bathymetric Spectral and Flow Analysis (MATLAB)

## Overview

This repository provides a MATLAB-based workflow for analyzing bathymetric digital elevation models (DEMs) to identify dominant spatial wavelengths and evaluate geomorphic structure across scales.

The workflow integrates spectral analysis, spatial filtering, and flow routing to identify and isolate wavelengths that dominate seafloor morphology. It is designed for reproducibility and easy adaptation to other bathymetric datasets.

## Methods

### Spectral Analysis

A 2D Fast Fourier Transform (FFT) is applied to the DEM to compute the power spectrum. Radial frequency is used to convert results into wavelength space, allowing identification of dominant spatial scales.

A power-law relationship is fit in log-log space and removed to produce a normalized (detrended) spectrum, highlighting deviations from background scaling behavior. Essentially removing the long wavelength signal. 

### Radial Averaging

The detrended power spectrum is radially averaged using logarithmically spaced bins. Mean power and standard deviation are calculated to identify statistically significant wavelength peaks. The data is binned and averaged in order to improve computational efficiency. 

### Spatial Filtering

Gaussian high-pass filters are applied to isolate bathymetric features at different spatial scales (e.g., 500 m, 1 km, 5 km). Each filtered DEM is preserved and exported independently.

### Flow Analysis

Filtered DEMs are analyzed using TopoToolbox. Bathymetry is inverted to enable flow routing, and flow accumulation is computed for each spatial scale.

### Flow Network Extraction

Flow networks are extracted using a drainage area threshold and exported as shapefiles for analysis in GIS platforms such as QGIS.

## Resolution Handling

The DEM is downsampled **only for spectral (FFT) analysis** to improve computational efficiency. All filtering and flow routing are performed on the original full-resolution DEM to preserve geomorphic detail.

## Inputs

* Bathymetric DEM (GeoTIFF format)
* Grid resolution (meters)

## Outputs

For each filter scale:

* Filtered DEM (GeoTIFF)
* Flow accumulation grid (GeoTIFF)
* Flow line network (shapefile)

Additional outputs:

* Raw and normalized power spectra
* Radially averaged spectra

## Requirements

* MATLAB
* Image Processing Toolbox
* TopoToolbox

## Usage

1. Update the input DEM file path.
2. Set parameters:

   * Downsample factor (FFT only)
   * Filter cutoff wavelengths
   * Flow accumulation threshold
    
3. Run the script or live script.

## Notes

* Input DEMs should be projected in a consistent coordinate system (e.g., UTM).
* Outputs are written using the specified EPSG code.
* Flow network density depends on the chosen drainage area threshold.

## Author

Sarah Rysanek
PhD Dissertation
