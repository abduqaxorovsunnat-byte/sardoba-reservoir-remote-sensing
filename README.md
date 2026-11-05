# sardoba-reservoir-remote-sensing

Google Earth Engine scripts for analysing surface-water loss at the Sardoba Reservoir (Uzbekistan) between April and May 2020, using Sentinel-1 radar and Sentinel-2 optical imagery, with 2019 as a control year.

Author: Sunnat Abdukaxorov

## What is here

`earth_engine_scripts.js` contains six separate Earth Engine scripts, each in its own numbered section:

1. Radar, simple threshold (20 km box)
2. Radar, change detection with patch filter (40 km box, orbit 151)
3. Sensitivity to the water thresholds (radar and optical)
4. Optical (Sentinel-2) lost water at NDWI > 0
5. Figure export to Google Drive
6. Figure 5 map and Sentinel-1 scene dates

## How to run

1. Open https://code.earthengine.google.com (a free Earth Engine account is needed).
2. Clear the editor, paste ONE numbered section, and press Run.
3. Read the results in the Console tab.

## Data

- Sentinel-1 GRD (COPERNICUS/S1_GRD), interferometric wide mode, VV, descending, relative orbit 151
- Sentinel-2 surface reflectance (COPERNICUS/S2_SR_HARMONIZED), cloud cover below 30%

Both are open data from the Copernicus programme, available through Google Earth Engine.

## Notes

The thresholds (-16 dB for radar, NDWI above 0 for optical) are assumed and tested for sensitivity in section 3. Results are for research use and have no independent field validation.
