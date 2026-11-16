# The mountain network: roads, valleys, passes, settlements, and a railway

Extends docs/plan.md (2026-11-16) from the KKH corridor to the network that serves the people of Gilgit-Baltistan,
Chitral, Hazara, Kaghan, Neelum and the Murree hills. Three questions drive it:

1. **Access as a public duty.** For every village: how far, over which road, in which season, to the nearest school,
   clinic, hospital, market, power and network; how often that road is cut, and how long for. A people's state plans
   from the most cut-off household outward, not from the trunk road in.
2. **Hazard on every road, and what to do about it.** Landslide, rockfall, debris flow, flood, GLOF, avalanche,
   snow closure and earthquake susceptibility per 100 m segment; where a segment fails repeatedly (the Astore road
   below Dashkin and Gorikot, Tatta Pani and Shatial on the KKH, the Kaghan and Neelum slides), compare options:
   stabilise (drainage, rock bolts, nets, retaining walls, shelters), realign on the opposite bank, bridge, or tunnel,
   each with cost, risk removed, construction risk and ecological cost.
3. **A railway for the mountains.** A planning-level alignment for a Havelian - Khunjerab railway (the CPEC
   proposal, ~680 km) and for links to Skardu, Gilgit - Chitral and Muzaffarabad, under railway geometry limits, with
   tunnels and bridges where the terrain demands, energy from hydro along the route, and the same hazard and ecology
   costs as the roads.

## Network (config/corridor.toml), routed through OpenStreetMap town to town

| group | roads |
|---|---|
| trunk | KKH N-35 Hasan Abdal - Khunjerab; M-15 (E-35) Burhan - Thakot; GT Road N-5 Islamabad - Hasan Abdal; M-1; Islamabad - Swat Expressway - Dir - Lowari Tunnel - Chitral (N-45); old Murree Road (N-75) and Murree Expressway to Kohala; N-15 Mansehra - Kaghan - Naran - Babusar - Chilas; S-1 Jaglot - Skardu |
| access | Gilgit - Ghizer - Phander - Shandur - Mastuj - Chitral; Astore - Chillum; Neelum: Muzaffarabad - Athmuqam - Keran - Sharda - Kel - Taobat; Shimshal; Chapursan; Nagar - Hoper; Naltar; Yasin; Ishkoman; Shigar - Askole; Khaplu - Hushe; Kalash valleys: Bumburet, Rumbur, Birir |
| passes | Burzil (Chillum - Minimarg); Deosai (Skardu - Sadpara - Deosai - Chillum); Noori Top (Sharda - Jalkhad) |

Still to place: "Baboon" (name to confirm: Babusar is in; if it means another pass or valley, add it to the config);
Hunza side valleys beyond Shimshal and Chapursan (Misgar, Hispar), Bagrot, Haramosh, Rondu villages, Kharmang,
Basho, Thalle, the Kalash summer pastures; the China side (G314 to Tashkurgan) as an extension.

## Settlements and services
- All OSM places in the box (6,401 points: city, town, village, hamlet), with district (OSM admin level 6) and
  population from WorldPop 100 m / GHS-POP; buildings from Google Open Buildings v3 and OSM.
- Services: OSM amenities (school, clinic, hospital, pharmacy, police, bank), health facility lists (Pakistan
  health facility registry where open), power lines and plants (OSM, GRW solar), mobile coverage (OpenCelliD).
- Accessibility: travel time over the road network (speeds by class and gradient, seasonal closures from snow and
  slide records) and off-road walking (Tobler's hiking function on GEM-15 slope), to the nearest facility of each
  type; Malaria Atlas Project friction surface as a check.
- Kalash villages: all OSM places in Bumburet, Rumbur and Birir (Krakal, Brun, Anish, Batrik, Balanguru, Grom,
  Guru, ...), named as their residents name them, with the threats they face (floods of 2015 and 2022, road cuts).

## Data at the highest resolution available
| layer | best open source | resolution |
|---|---|---|
| terrain | Google GEM-15 dtm (2022-24 stereo), gaps from Copernicus GLO-30 | 15 m (30 m) |
| local high resolution | High Mountain Asia 8 m DEM (NSIDC, WorldView stereo) where it exists | 8 m |
| imagery | Sentinel-2 L2A (10 m), Landsat 8/9 (30 m, 15 m pan), PlanetScope where licensed | 3-10 m |
| radar | NISAR L-band (GUNW, GOFF, GCOV), Sentinel-1 C-band (HyP3) | 20-80 m |
| land cover | kashmir-landcover v2 scheme and maps, ESA WorldCover 10 m, Dynamic World | 10-30 m |
| buildings | Google Open Buildings v3, OSM | footprints |
| glaciers, lakes | RGI 7.0, ICIMOD glacial lake inventory, Sentinel-2 lake mapping | 10 m |
| snow | MODIS snow cover (uib-snow-cover), Sentinel-2 snow | 10-500 m |
| rain | CHIRPS, IMERG, ERA5-Land | 5-10 km |
| geology, faults | GEM active faults, geological map of Pakistan (GSP), MMT/MKT/Raikot fault traces | 1:250k-1:1M |

## Astore road landslide study (first detailed site)
1. Inventory: map every scar and fan along the Astore road from Sentinel-2 (2016-now) and Landsat (1990-now)
   change, NASA/ICIMOD catalogues, news reports of closures (dated, cited), field reports where available.
2. Motion: NISAR GUNW stacks (12-day, 2025-now) and GOFF offsets for fast movers; Sentinel-1 SBAS (HyP3) for
   2017-now; velocity per slope unit.
3. Factors per slope unit: slope, relief, curvature, aspect, lithology, distance to faults and river, road cut
   height (DSM minus DTM and slope break), land cover, snowmelt timing, rainfall extremes, PGA.
4. Models: frequency ratio, logistic regression, random forest, with spatially blocked validation; AHP as the
   expert check. Output: susceptibility and expected closures per 100 m.
5. Options per hotspot: drainage and slope protection; rockfall barriers and galleries; realignment to the other
   bank (cost surface on GEM-15); tunnels where the valley wall is unstable for > 500 m. Cost per option from
   regional unit costs, risk removed in closure-days and lives, ecological cost.

## Railway planning (eco-socialist Third Pole)
- Geometry: ruling gradient 1.5-2.5 % (mountain line; 4 % only with rack or for short ramps), minimum curve radius
  300-500 m (mountain, 80-120 km/h design), vertical curves, tunnels where needed; spiral and loop tunnels to gain
  height (as on the Albula and Lötschberg lines and the Qinghai-Tibet and Sichuan-Tibet railways).
- Method: least-cost corridor on GEM-15 with a gradient-constrained path search (graph on a 3D lattice with
  gradient and curvature limits), costs = construction (open line, bridge, tunnel per km) + hazard + ecology, and
  settlements served as attractors; output polylines with chainage, gradient profile, tunnel and bridge lists.
- Comparisons: the 2016-2018 CPEC Havelian - Khunjerab feasibility notes (published figures only), the existing ML-1
  terminus at Havelian, the Diamer-Bhasha reservoir (the line must climb around it).
- Ecology and people: protected areas (Khunjerab NP, Deosai NP, Central Karakoram NP) as costs or no-go; villages
  within reach of a station; freight for valleys (food, fuel, materials) against tourism; electricity demand met
  by run-of-river hydro along the line.

## Standards
EIA: IFC Performance Standards, World Bank ESF, ADB SPS 2009, Pakistan EPA 1997 and GB EPA, IEE/EIA Regulations 2000.
Roads: AASHTO Green Book geometry, PIARC road tunnel and road safety manuals, Pakistan NHA specifications.
Railway: UIC codes and leaflets for geometry and tunnels, EN 1991-2 loads, Chinese TB 10001 (mountain railways)
for comparison. Risk: ISO 31000, Sendai Framework.
