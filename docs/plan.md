# Karakoram Highway corridor: multi-hazard risk, and safer alignments

Aim: map flood, landslide and debris flow, glacial lake outburst flood (GLOF), erosion and earthquake hazard, exposure
and vulnerability along the Karakoram Highway (N-35) and its feeder roads, then test realignments, bridges and
tunnels that lower the overall risk to people, the economy and ecosystems. Every number is traceable to open data and
code; every step is documented and repeatable.

## Corridor (config/corridor.toml)
- N-35 Karakoram Highway: Hasan Abdal / Havelian - Abbottabad - Mansehra - Thakot - Besham - Dasu - Chilas - Raikot -
  Gilgit - Hunza (Attabad tunnels) - Sost - Khunjerab Pass (and the China G314 side to Tashkurgan as an extension)
- N-15 Mansehra - Naran - Babusar Pass - Chilas
- Astore road (Jaglot / Raikot bridge - Astore - Rama, on to Deosai and Chillum)
- Gilgit - Ghizer - Gahkuch - Phander - Shandur (to Chitral)
- S-1 Jaglot - Skardu (Strategic Road), Skardu - Khaplu, Skardu - Deosai
- Shahrah-e-Kaghan, Nagar and Hispar roads, Chapursan and Misgar roads (second pass)
- Existing and planned works: Attabad and Lowari-type tunnels, Thakot - Raikot CPEC realignment (planned), Diamer-Bhasha
  dam (reservoir will drown ~100 km of the KKH: the planned realignment is part of the comparison), Dasu, Tarbela,
  Satpara, and run-of-river schemes in OSM / WAPDA lists

## Data, best open first
| Need | Source | Res. | Notes |
|---|---|---|---|
| Surface (DSM) | Copernicus DEM GLO-30 (TanDEM-X derived, 2011-15) | 30 m | best open global DSM; GLO-90 for checks |
| Surface (DSM) | ALOS AW3D30 v4.1, NASADEM, SRTM | 30 m | inter-comparison, voids |
| Terrain (DTM) | FABDEM (forest/building removed GLO-30) | 30 m | for hydrology and slope where trees bias the DSM |
| Higher resolution | High Mountain Asia 8 m DEM (NSIDC HMA_DEM8m, from WorldView/GeoEye stereo) | 8 m | patchy; mosaic where present |
| SAR DEM (built here) | Sentinel-1 InSAR pairs (ASF HyP3, 6-12 day, short baseline) | 20-40 m | tutorial: make a DEM, and show why S-1 is weak for DEMs (short baselines, layover) |
| SAR change | Sentinel-1 InSAR / SBAS velocity (HyP3, LiCSAR, EGMS-like) | 40-80 m | active landslides, creeping slopes |
| Optical | Sentinel-2, Landsat, PlanetScope NICFI (no: not tropics) | 10-30 m | debris flows, scars, lakes |
| Land cover | kashmir-landcover v2 (this workspace) | 30 m | classes from config/classes.toml |
| Rain | CHIRPS, IMERG, ERA5-Land | 5-10 km | triggers, return levels |
| Rivers, floods | MERIT Hydro, HydroSHEDS, GloFAS reanalysis, JRC surface water, 2010 and 2022 flood extents | 90 m-0.1 deg | |
| Glacial lakes | Wang et al. HMA glacial lakes, ICIMOD GLOF inventory, RGI 7.0 glaciers | | GLOF runout with a simple flow model |
| Landslide inventory | NASA COOLR, ICIMOD, literature inventories for KKH (e.g. Hunza, Attabad 2010) | points | training and validation |
| Earthquake | GEM global hazard (PGA 475 yr), USGS catalogue, GEM active faults (MMT, MKT, Raikot fault) | | |
| Geology | USGS/ICIMOD geological map of Pakistan, lithology | 1:1M-1:250k | |
| Exposure | OSM roads, bridges, tunnels, buildings; Open Buildings; WorldPop / GHS-POP | | |
| Ecology | WDPA protected areas (Khunjerab NP, Deosai NP, Central Karakoram NP), KBAs, habitat | | |

## Method
1. Hazard susceptibility per hazard: frequency ratio and weights of evidence (baseline), logistic regression and random
   forest (data driven, spatially blocked cross-validation), and AHP / weighted overlay where inventories are thin
   (expert weights, consistency ratio < 0.1). Report AUC, spatially blocked; compare methods (as in tree comparison).
2. GLOF: lake inventory, dam type, volume (Cook and Quincey 2015 area-volume), peak discharge (Huggel), runout path on
   the DTM (flow routing, modified single flow), road crossings hit.
3. Flood: HAND (height above nearest drainage), return-period discharge from GloFAS, bridge and causeway exposure.
4. Earthquake: PGA and distance to active faults; co-seismic landslide susceptibility (Newmark-style).
5. Risk = hazard x exposure x vulnerability, per 100 m road segment, per hazard and combined (ISO 31000 / UNDRR terms).
6. Realignment: least-cost path over a cost surface (risk, slope, construction cost per km for road on cut, bridge,
   tunnel; ecology as cost and as no-go), with tunnel candidates where the surface route crosses high-risk cliffs;
   output new segments as polylines with length, cost, risk avoided. Multi-objective: Pareto set of alignments, not
   one answer.
7. Costs: unit costs from published Pakistani and regional contracts (NHA, CPEC KKH Phase II, ADB/World Bank road
   appraisals), tunnel cost per km, maintenance and closure-days cost; risk in expected annual loss and lives.
8. EIA screening of each option (IFC Performance Standards, World Bank ESF, ADB SPS, Pakistan EPA 1997 / GB EPA,
   IEE/EIA Regulations 2000); engineering checks (AASHTO Green Book geometry, mountain grade limits, PIARC tunnel and
   road safety, ISO 31000 risk, Sendai Framework).

## Outputs
Per step: script in src/, data in data/, figures in figures/, a record in docs/record.tex (LaTeX, one source), the
public tutorial in tutorial/ (SAR DEM tutorial first), and map layers for the Third Pole atlas.
