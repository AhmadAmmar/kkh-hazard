# Ground data: soil, rock, faults, minerals (status 2026-11-16)

| layer | source | resolution | state |
|---|---|---|---|
| soil properties (clay, sand, silt, coarse fragments, bulk density, organic carbon, pH; 0-5, 15-30, 60-100 cm) | ISRIC SoilGrids 2.0 via Earth Engine | 250 m | downloading (src/06, --only soil) |
| soil texture class (USDA) | OpenLandMap v02 | 250 m | downloading |
| mineral indices: AlOH, MgOH/carbonate, ferric, ferrous, quartz | ASTER L1T 2000-2008, snow-free summer medians, computed in Earth Engine | 15-30 m (90 m TIR) | downloading (--only aster) |
| surface mineralogy (spectral identification, band depth) | NASA EMIT L2B MIN (ISS imaging spectrometer) | 60 m | 1,532 granules over the box 2022-08 to 2026-10; newest 60 downloaded (data/raw/geology/emit) |
| active faults, sutures | HimaTibetMap (Styron, Taylor, Yin et al.) + GEM Global Active Faults | 1:1M-ish lines | done (inventories.gpkg: htm_HimaTibetMap 11, htm_sutures 6; faults 28) |
| mineral occurrences | USGS MRDS | points | done, 25 in box (thin for the north) |
| geology of South Asia (geo8ag, Wandrey 1998) | USGS, ScienceBase | 1:5M-ish polygons | MANUAL: ScienceBase blocks scripts; download in a browser from https://www.sciencebase.gov/catalog/item/60c3b89fd34e86b93897ef19 into data/raw/geology/usgs_geo8ag/ |
| Geological map of Pakistan 1964 (Abu Bakr and Jackson) | National Library of Australia scan, https://nla.gov.au/nla.obj-233602228/view | 1:2,000,000 | to fetch (TIF download from the viewer), georeference |
| Geology of trans-Himalayan Gilgit and Baltistan (Abu Bakr 1965, GSP Records XI/3) | GSP scan | report | downloaded; the Plate 1 fold-out map is NOT in the scan (blank page): need another copy |
| Geology and geochemistry of part of Hunza (GSP) | GSP scan | 1:250,000-ish | downloaded; Fig. 2 geological map of upper Hunza (B/W hatching, graticule ticks) and Fig. 3 section Gulmit - Khunjerab; Table 4 rock chemistry (SiO2, Fe2O3, Al2O3, CaO, MgO) to transcribe |
| Geology and slope stability along the Gilgit - Skardu road (GSP) | GSP report page | report | link page found; PDF to locate |
| Geological map of North Pakistan (Searle and Khan 1996) | Oxford / GSP | 1:650,000 | not open online: needs a scan (library) |
| GSP KKH road log (Highways geological map series no. 1, 2000) | GSP | strip map | not found online |
| Kohistan preliminary geological map (Peshawar Univ. Geol. Bull. 11, 1979) | University of Peshawar | 1:250,000 | not found online |
| global lithology GLiM (Hartmann and Moosdorf 2012) | PANGAEA (0.5 deg grid only) / CCGM gdb | 1:3.75M | coarse: only as a fallback |

## Georeferencing and digitising (src/07_georeference.py)
- Georeference from ground control points (graticule crossings, confluences, bridges, peaks) with leave-one-out RMS
  reported; first or second order for printed maps, thin-plate spline only for distorted paper.
- Digitise colour maps by nearest legend colour in CIELAB with a majority filter. Black-and-white hatched maps (the
  GSP Hunza map) need pattern classification or hand tracing in QGIS on the georeferenced scan; done by hand where
  the hatching is ambiguous, recorded as such.
- Expected accuracy is set by the source scale: 1:250,000 drawn to 0.5 mm is ~125 m at best; 1:650,000 ~325 m;
  1:2M ~1 km. These layers are used as categorical covariates at that scale, not at 15 m.

## Why the finest layers are spectral, not maps
At 15-60 m the best rock information comes from ASTER and EMIT spectra (clays and micas, carbonates, iron oxides,
quartz) and from texture in GEM-15; published maps add the names and ages of units at 1:250,000 and coarser.
