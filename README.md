# KKH hazard corridor

Multi-hazard susceptibility, risk and safer alignment options along the Karakoram Highway (N-35) and its feeder roads:
Babusar (N-15), Astore, Gilgit - Ghizer - Phander - Shandur, and Jaglot - Skardu (S-1). Floods, landslides and debris
flows, glacial lake outburst floods, erosion and earthquakes; exposure of roads, bridges, tunnels and dams; candidate
realignments and tunnels with cost, risk and environmental screening.

Open data only (OpenStreetMap, Copernicus DEM, Sentinel-1, Sentinel-2, Landsat, GEM, USGS). See `docs/plan.md`.

| Step | Script | State |
|---|---|---|
| Roads, bridges, tunnels, dams, 100 m segments | `src/01_roads.py` | done |
| Copernicus GLO-30 elevation (30 m surface) | `src/02_dem.py` | done |
| Google GEM-15 elevation (15 m surface and bare-earth terrain) | `src/02b_gem15.py` | downloading |
| Terrain predictors (slope, curvature, TWI, SPI, HAND, relief) | `src/03_terrain.py` | running |
| Faults, earthquakes, landslide inventories | `src/05_inventories.py` | faults, quakes done |
| SAR elevation tutorial (Sentinel-1 InSAR) | `tutorial/sar_dem/` | pairs found |
