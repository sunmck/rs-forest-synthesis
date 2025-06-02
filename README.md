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

To distinguish different change agents, data from the European Foreest Disturbance Atlas (Viana-Soto and Senf 2025) is used.

## Key Findings
* Stand-replacing disturbances (e.g., windthrow, severe bark beetle outbreaks) cause the most significant structural losses
* Structural declines often precede visible disturbances—especially in aboveground biomass
* Post-disturbance recovery varies; canopy cover recovers fastest
* Remote sensing shows strong potential for early detection, forest inventory updates, and post-disturbance monitoring

## Polygon-based Analysis
* Computation of hexagonal maps of descriptive statistics of the single RS products
* Bivariate correlations between the RS product attributes

## Pixel-based Analysis
* Analysis of a timeseries of descriptive statistics of the forest structure attributes and FCSD for distinct disturbance patterns at 10 m resolution
* Definition and mapping of seven disturbance patterns based on RS-derived classification of forest condition
* Two regions with available ground-truth disturbance data and known management histories used to define and validate disturbance patterns and test robustness of RS-derived classifications at finer spatial scales

### Disturbance Pattern Classification
Forest disturbances are classified into seven distinct patterns based on the intensity, scale and primary change agent observed in the RS products. 

![distpatterns](https://github.com/user-attachments/assets/f3a18c51-ecc7-4398-b28d-1b184526c210)




