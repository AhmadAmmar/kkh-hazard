#!/bin/bash
# Rebuild everything that depends on the corridor after the roads changed (2026-11-16). Old outputs archived, not deleted.
cd ~/dev/kkh-hazard
export GDAL_DATA=C:/Users/Owner/anaconda3/envs/eo/Library/share/gdal PYTHONIOENCODING=utf-8
PY=/c/Users/Owner/anaconda3/envs/eo/python.exe
A=data/processed/_archive_2026-11-16; mkdir -p $A
step(){ echo "== $1 $(date +%T)"; shift; "$@" 2>&1 | grep -v -i warn | tail -4; echo "exit ${PIPESTATUS[0]}"; }
# step roads (done 19:48) $PY src/01_roads.py
step dem $PY src/02_dem.py
for f in dem_gem15_dtm_utm.vrt dem_gem15_dsm_utm.vrt dem_gem15_dsm_std_utm.vrt; do [ -f data/processed/$f ] && mv data/processed/$f $A/ ; done
step gem15 $PY src/02b_gem15.py
ls data/raw/dem/gem15/tiles/tile_*.tif > data/processed/gem15_list.txt
i=1; for b in dtm dsm dsm_std; do /c/Users/Owner/anaconda3/envs/eo/Library/bin/gdalbuildvrt.exe -q -b $i -vrtnodata -9999 -input_file_list data/processed/gem15_list.txt data/processed/dem_gem15_${b}_utm.vrt; i=$((i+1)); done
[ -d data/processed/terrain ] && mv data/processed/terrain $A/terrain_v1
step terrain $PY src/03_terrain.py data/processed/dem_glo30_utm_full.tif
# inventories.gpkg kept: it now also holds ground layers (src/06); src/05 overwrites only its own layers
step inventories $PY src/05_inventories.py
step maps_all $PY src/90_maps.py --area every
echo "== rebuild done $(date +%T)"
