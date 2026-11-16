# Muddasir Shah's work used as a reference

Credit: Muddasir Shah (ex-colleague). Used for ideas and methods only; nothing copied without credit.

1. "SAR interferometry, phase unwrapping with SNAPHU and DEM generation" (LinkedIn article):
   https://www.linkedin.com/pulse/sar-interferometry-phase-unwrapping-snaphu-dem-generation-shah-fgxsf/
   LinkedIn blocks automated reading; the text is still to be read (pasted in by hand) and folded into tutorial/sar_dem.
2. "From VV/VH dB to covariance matrices: first lessons from NISAR L2" (LinkedIn article):
   https://www.linkedin.com/pulse/from-vvvh-db-covariance-matrices-first-lessons-nisar-l2-muddasir-shah-vohxf/
   Same: text to be read.
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
