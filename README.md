# Built-up Surface Ensemble Model for Romania Based on OSM, MSBF, and GHS Data Sources Using Triple Collocation Analysis
*Authored by Zsolt Magyari-Sáska and Ionel Haidu*

The repository contains two R scripts:

- **build_ensemble.R** – creates the ensemble model at the settlement level  
- **morris.R** – computes the Morris sensitivity-analysis index for the parameters used in the creation script  

---

### Parameters to set or modify in the scripts

- **shp_path** – path to the polygon layer of settlements  
- **osmWeight** – confidence level for the OSM data source when Triple Collocation Analysis cannot be performed  
- **GHSdiv** – weight-reducing factor for the GHS data source when it is the only source  
- **MSBFdiv** – weight-reducing factor for the MSBF data source when it is the only source  
- **cell_area** – area of a raster cell  

---

By default, the scripts produce the ensemble model for settlements in **Covasna County**.  
Before running them, unzip **GHS_CV.zip**, **MSBF_CV.zip**, and **OSM_CV.zip** from the **DATA** folder. The results will be generated in the **ENS** folder. 

---

The ensemble model covering the entire territory of Romania is available on Zenodo: [https://doi.org/10.5281/zenodo.16742094](https://doi.org/10.5281/zenodo.16742094)


