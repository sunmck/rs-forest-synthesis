# Forest Structure Dynamics in Germany - A Synthesis of Remote Sensing Products

## Overview
This repository contains the code and data preprocessing pipelines used in my Master’s thesis, which investigates forest structure changes across Germany from 2017 to 2023 by synthesising three selected remote sensing (RS) data products. 

The analysis begins with an introductory examination of the individual, univariate RS products to identify hotspots with significant changes of forests between 2017 and 2023. The RS products
are then combined into multivariate data structures. Initially, exploratory bivariate correlations are analysed, followed by a national investigation of the effects of different disturbance patterns on forest structure dynamics, within areas previously identified as having undergone significant changes. For this, forest disturbances are classified into seven distinct disturbance patterns, based on the RS products, and unique characteristics are explored throughout the disturbance process, from pre-disturbance to post-disturbance phases.

Subsequent regional analyses allow for investigations into forest structure dynamics after specific disturbance events. For this purpose, two regions of interest (ROI) with available groundtruth reference data on forest disturbances and their management are selected, allowing the disturbance patterns to be defined using these datasets. These ROIs will serve the purpose to assess whether the effects of different disturbance patterns on forest structure dynamics, as observed at the national level, can be confirmed and to test the applicability of the RS products for smaller-scale analyses.

## Data
The RS products used in this analysis comprise three key products, provided by the German Aerospace Centre (DLR):
* Forest Structure (Kacic et al. 2023), 
* Forest Canopy Cover Loss (FCCL) (Thonfeld et al. 2022) and 
* Fractional Cover of Standing Deadwood (FCSD) (Schiefer et al. 2024). 

These products cover the entire study area of Germany, offering high temporal (monthly or yearly) and spatial (10 m) resolution for the period between 2017 and 2023. The forest structure attributes analysed are: Total Canopy Cover (TCC), Canopy Height (CH), Above-ground Biomass Density (AGBD), Foliage-height Diversity Index (FHDI) and FCSD.

To distinguish different change agents, data from the European Foreest Disturbance Atlas (Viana-Soto and Senf 2025) is used. For the regional analyses, data from the Copernicus Emergency Management Service (EMS) and the Bavarian Forest National Park are used.

## Installation / Requirements

To run the code in this repository, the following Python (and R) libraries are required:

### Python
- rioxarray
- xarray
- rasterio
- geopandas
- pandas
- matplotlib
- numpy
- dask
- scipy
- shapely
- seaborn
- lexcube
- ipywidgets
- pyvista
- jupyterlab

To install the required Python packages, run:

Run `pip install -r requirements.txt` to install the dependencies.

### R
- gdalUtilities
- terra
- sf
- dplyr
- tidyverse

## Repository Structure

```
├── data/ # Placeholder for raw/processed data (not uploaded)
├── results/ # Output tables, plots, etc.
├── scripts/
│ ├── polygon_based/
│ │ ├── reproject.R
│ │ └── poly_calculatestats_hexagons.R
│ │ └── poly_calculatestats_otherpolygons.R
│ └── pixel_based/
│ └── pixel_datacube_mainanalysis.py
└── README.md
└── requirements.txt # Python requirements
```

> **Note:** The `data/` folder is included for structure, but raw and processed datasets are not uploaded due to size and licensing restrictions. Please contact the data providers or me for access.
> 

## Key Findings
* Stand-replacing disturbances (e.g., windthrow, severe bark beetle outbreaks) cause the most significant structural losses
* Structural declines often precede visible disturbances—especially in aboveground biomass
* Post-disturbance recovery varies; canopy cover recovers fastest
* Remote sensing shows strong potential for early detection, forest inventory updates, and post-disturbance monitoring

## Polygon-based Analysis
* Aggregation of pixel-wise information to polygons (hexagons or other polygons such as forest growth regions) by means of calculating descriptive statistics of the single RS products (`poly_calculatestats_hexagons.R` and `poly_calculatestats_otherpolygons.R`)
* Computation of time series and bivariate correlations between the RS product attributes (`TODO`)

## Pixel-based Analysis (`pixel_datacube_mainanalysis.py`)
* Analysis of a timeseries of descriptive statistics of the forest structure attributes and FCSD for distinct disturbance patterns at 10 m resolution
* Definition and mapping of seven disturbance patterns based on RS-derived classification of forest condition
* Two regions with available ground-truth disturbance data and known management histories used to define and validate disturbance patterns and test robustness of RS-derived classifications at finer spatial scales

### Disturbance Pattern Classification
Forest disturbances are classified into seven distinct patterns based on the intensity, scale and primary change agent observed in the RS products. 

![distpatterns](https://github.com/user-attachments/assets/f3a18c51-ecc7-4398-b28d-1b184526c210)

### Disturbance Phase Classification
For all pixels that fall into one of the disturbance patterns in a specific year, a timeseries of descriptive statistics of the forest structure attributes and FCSD is calculated. This approach enables the analysis of forest structure dynamics across the entire available time series of disturbed pixels. Specific phases of a disturbance are analysed, including:

<img src="https://github.com/user-attachments/assets/3fa12329-4297-4e38-94df-ef25ad7cb223" alt="methods_prepostdist" width="400"/>

## License

This repository is licensed under the MIT License.

*Note:* Remote sensing products used in this analysis may have their own license restrictions—please refer to the data providers for usage guidelines.

## References

- Kacic, P.; Thonfeld, F.;  Gessner, U.; Kuenzer, C. Forest  Structure Characterization in  Germany: Novel Products and  Analysis Based on GEDI, Sentinel-1  and Sentinel-2 Data. Remote Sens.  2023, 15, 1969. https://doi.org/  10.3390/rs15081969S. 
- Thonfeld, F.; Zhu, Z.; Gessner, U.; Dech, S.; et al. A First Assessment of Canopy Cover Loss in Germany’s Forests after the 2018–2020 Drought Years. Remote Sens. 2022, 14(3), 562. https://doi.org/10.3390/rs14030562.
- Schiefer, F.; Schmidtlein, S.; Frick, A.; Frey, J.; Klinke, R.; Zielewska-Büttner, K.; Junttila, S.; Uhl, A.; Kattenborn, T. UAV-Based Reference Data for the Prediction of Fractional Cover of Standing Deadwood from Sentinel Time Series. Karlsruhe Institute of Technology (KIT), University of Freiburg, FVA Baden-Württemberg, University of Eastern Finland, Leipzig University, iDiv.
- Viana-Soto, A.; Senf, C. European Forest Disturbance Atlas. Version 2.1.1. Aug. 20, 2024. DOI: 10.5281/zenodo.13333034. Available: https://zenodo.org/records/13333034.

## Acknowledgements

Thank you for the supporting and insightful supervision and all resources provided by the Chair of Remote Sensing of the University of Würzburg and the DLR to make this thesis possible!

