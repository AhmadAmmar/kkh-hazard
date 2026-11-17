# Mountain railways for Gilgit-Baltistan and Hazara: first planning sketch (2026-11-17)

**What this is.** A reproducible, planning-level first answer: where a railway could run, how steep it must be, how
much of it would be in tunnel or on bridge, and roughly what it would cost. Made with open terrain data (Copernicus
GLO-30 at 90 m), the OSM road network, and stated assumptions (config/railway.toml). It is **not** a feasibility
study: geology, hazard, hydrology, construction method, operations and costs all need survey-level work.
Code: src/40_railway.py. Outputs: outputs/railway_<line>_{profile,structures,summary}.csv, data/processed/railway_<line>.gpkg,
figures/railway_<line>_{map,profile}.png.

## Method in brief
1. **Plan line.** Least-cost path on a 90 m grid within a band around the existing road (valleys are the only
   routes). The cost grows with slope and with height above the road's own level at that point, so the line keeps
   to the valley on benches along its sides instead of crossing ridges. Ground below the full pool of the future
   **Diamer-Bhasha (1,160 m)** and **Dasu (950 m)** reservoirs, plus 15 m freeboard, is costed as lake, so the line
   is pushed onto the valley sides above them. The line is then smoothed to a minimum curve radius (500 m main line,
   400 m branches, measured over 250 m chords).
2. **Profile.** Dynamic programming along the line in 100 m steps:
   - **Grade limit:** the rail may rise or fall by at most the ruling grade (2.5 %) per step.
   - **Open line:** rail within 12 m of the ground is open line on cut or fill, costing more the deeper the cut or
     higher the fill.
   - **Tunnel:** rail more than 12 m below the ground is tunnel, with a mild charge for depth.
   - **Bridge:** rail more than 12 m above the ground is bridge, costing more as piers get taller, and never over
     150 m.
   - **Summit tunnel:** the Khunjerab line may end in a summit tunnel portal up to 1,000 m below the pass.
   The cheapest profile is kept.
3. **Cost.** Lengths x unit cost ranges (US$ million per km, 2024): open line 8-20, tunnel 30-70, bridge 25-50.
   These are indicative ranges for mountain lines, not quotes, and should be replaced with regional contract data.

## Results (first run)
| line | length | rail heights | tunnel | bridge | longest tunnel | indicative cost |
|---|---|---|---|---|---|---|
| Havelian - Khunjerab (CPEC proposal) | 634 km | 860 m -> summit tunnel portal ~4,110 m under the 4,700 m pass | 231 km (37 %), 199 sections | 63 km | 28.6 km | US$ 11-26 bn |
| Jaglot - Skardu | 159 km | 1,420 -> 2,250 m | 7.6 km | 7.2 km | 2.5 km | US$ 1.6-3.8 bn |
| Gilgit - Shandur - Chitral | 325 km | 1,470 -> 3,290 m (Shandur) -> 1,510 m | 86 km | 10 km | 42.2 km | US$ 4.7-11.1 bn |

What the profile of the main line shows (figures/railway_havelian_khunjerab_profile.png):
- **Havelian to Thakot (km 0-95):** the line crosses the Hazara uplands near 1,200 m and drops to the Indus at
  about 520 m.
- **Reservoir reach (km 135-335):** along the future Dasu and Diamer-Bhasha lakes it has to hold 965-1,175 m on the
  valley sides. That is a side-hill line of many short tunnels and bridges, and it is the most expensive and most
  hazardous part of the route (rockfall, debris fans, the Raikot and MMT fault zones).
- **Gilgit and Hunza (km 335-560):** a steady 0.5-1 % climb to 2,500-2,800 m.
- **Khunjerab (last ~60 km):** the ground rises at about 3.4 %, steeper than a 2.5 % railway can follow. The line
  needs either extra length (loops or spirals, as on the Albula or Bernina lines) or a long summit (base) tunnel.
  The model chose a climbing approach in tunnel to a portal about 600 m below the pass.

## What is missing before this can be argued in public
- **Hazard and ecology:** hazard costs per cell (the susceptibility models of this project) and protected areas
  (Khunjerab National Park) as costs or no-go zones.
- **Geology:** faults, rock quality for tunnels, and squeezing ground (the MMT and MKT zones). Earthquake design.
- **Plan geometry:** a true 3D alignment search with curvature, so loops and spirals can be found automatically.
  Grade compensation in curves and in long tunnels.
- **Operations:** single track with passing loops, electrification from run-of-river hydro along the line, and
  freight against passenger demand (valley supply, minerals, the China trade). Who the line serves: the stations
  within reach of villages, against a corridor built only for transit.
- **Costs:** unit costs from published contracts in the region, and an uncertainty analysis.
- **Comparison:** the published CPEC Havelian - Khunjerab studies (figures only where they are publicly documented).
