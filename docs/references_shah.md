# Muddasir Shah's work used as a reference

Credit: Muddasir Shah (ex-colleague). Used for ideas and methods only; nothing copied without credit.

1. "SAR interferometry, phase unwrapping with SNAPHU and DEM generation" (LinkedIn article):
   https://www.linkedin.com/pulse/sar-interferometry-phase-unwrapping-snaphu-dem-generation-shah-fgxsf/
   Read (PDF supplied by the user, 2026-11-16). Taken into tutorial/sar_dem:
   - the phase budget dphi = flat + elevation + displacement + atmosphere + noise + error; only the first three carry
     signal, and which one is wanted decides a DEM job versus a deformation job;
   - DEM job: short time gap so displacement ~ 0; z = -lambda R0 sin(theta) / (4 pi B_n) dphi_topo;
     height of ambiguity h_a = lambda R0 sin(theta) / (2 B_n): larger baseline, finer height, denser fringes;
   - deformation job: remove topography with an external DEM; d_LOS = lambda / (4 pi) dphi, v = d_LOS / dT; no
     baseline term; DEM error and atmosphere leak straight into apparent motion;
   - SNAPHU: fewest tiles memory allows (2x2 on 8 GB), overlap 100-200 px (0 gives seams or empty output), tile cost
     threshold 500, processors = physical cores, initial method MCF not MST, cost mode TOPO for DEM, DEFO for motion.
2. "From VV/VH dB to covariance matrices: first lessons from NISAR L2" (LinkedIn article):
   https://www.linkedin.com/pulse/from-vvvh-db-covariance-matrices-first-lessons-nisar-l2-muddasir-shah-vohxf/
   Read (PDF supplied by the user). Taken in:
   - GCOV is geocoded, radiometrically terrain-corrected (gamma0) polarimetric covariance; diagonal terms are power
     (HHHH ~ HH, HVHV ~ HV, VVVV, VHVH), dB = 10 log10; sigma0 = gamma0 x rtcGammaToSigmaFactor when needed;
   - one file can hold frequencyA and frequencyB with different bandwidths and polarisations (e.g. 20+5 MHz;
     QP on A, DH on B): full-pol decompositions (Cloude-Pottier H/A/alpha, Freeman-Durden) only where off-diagonals
     exist; otherwise dual-pol ratios and texture;
   - mask 1-5 valid (sub-swath), 0 invalid, 255 fill: respect it, do not fill, keep it the same across dates;
   - very bright patches: corner-reflector geometry, double bounce, rough water; here: bridges, gabions, cut faces;
   - reference: Shiroma (2024) NISAR L2 GCOV product specification, JPL D-102274.
   Over the KKH our tracks are DHDH (dual-pol HH/HV) mostly, with QPDH on descending track 063 some dates: full-pol
   decomposition possible on those, for debris, scree, vegetation and wet ground classes.
3. Colab notebook (shared link; local copy in data/raw/reference_shah/, not redistributed):
   NISAR L2 GCOV quicklooks with h5py. Points taken:
   - read HDF5 with striding (downsample on read) to keep memory low;
   - path science/LSAR/GCOV/grids/frequencyA|B/<term>; plot the real diagonal terms (HHHH, HVHV, ...) only;
   - apply the GCOV mask: keep values 1-5 (valid, by number of subswaths), drop 0 (invalid) and 255 (fill);
   - diagonal terms are power: dB = 10 log10(power); display stretch by 2-98 percentiles;
   - parse polarisation mode from the file name (e.g. QPDH: frequency A quad-pol, B dual-pol).

## What this changes here
NISAR (L band, 24 cm) covers the corridor every 12 days from Nov 2025 to now, ascending (tracks 055, 156) and
descending (063), as GCOV, GSLC, GUNW (unwrapped interferograms) and GOFF (pixel offsets). L band keeps coherence on
steep, vegetated and debris slopes much better than Sentinel-1 C band, so:
- tutorial/sar_dem part B: NISAR GUNW, residual phase to height and to motion, compared with Sentinel-1 (part A);
- landslides: GUNW stacks for slow slopes, GOFF for fast movers (cm-m per 12 days);
- floods, debris flows, new scars: GCOV backscatter change (Shah's quicklook method, then change detection).
