cd ~/dev/kkh-hazard
export GDAL_DATA=C:/Users/Owner/anaconda3/envs/eo/Library/share/gdal PYTHONIOENCODING=utf-8
/c/Users/Owner/anaconda3/envs/eo/python.exe -u src/02b_gem15.py 2>&1 | grep -v -i warn | tail -3
ls data/raw/dem/gem15/tiles/tile_E*.tif > data/processed/gem15_list.txt
i=1; for b in dtm dsm dsm_std; do /c/Users/Owner/anaconda3/envs/eo/Library/bin/gdalbuildvrt.exe -q -b $i -vrtnodata -9999 -input_file_list data/processed/gem15_list.txt data/processed/dem_gem15_${b}_utm.vrt; i=$((i+1)); done
/c/Users/Owner/anaconda3/envs/eo/python.exe -u src/90_maps.py --area every 2>&1 | grep -v -i warn | tail -2
