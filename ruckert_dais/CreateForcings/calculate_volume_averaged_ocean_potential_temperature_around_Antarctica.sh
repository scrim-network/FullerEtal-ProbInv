#!/bin/bash

# calculate volume averaged ocean potential temperature anomaly (calculated with respect
# to 1961-1990) around Antarctica - 
# 52 to 70 degree South, for depth between 200 and 800 m.
# Last modified by Varada Vaidya on May 7, 2015

#vol_calculation="/woju/s0/vvv3/Kelsey_project/vol_calculation"
vol_calculation="/home/scrim/klr324/Ruckertetal_DAIS_codes/CreateForcings/vol_calculation"
mkdir -p ${vol_calculation} |exit

# Input file names 
MODEL="CNRM-CM5"
cellarea="/woju/s0/common/ethz_rsync/cmip5/historical/fx/areacello/${MODEL}/r0i0p0/areacello_fx_${MODEL}_historical_r0i0p0.nc"

thetao_hist="/woju/s0/data_archive/CMIP5/devel_CMIP5/3D_variables/3D_ocean_variables_giant_merged_files/thetao.mon.${MODEL}.historical.r1i1p1.185001-200512.nc"
thetao_rcp85="/woju/s0/data_archive/CMIP5/devel_CMIP5/3D_variables/3D_ocean_variables_giant_merged_files/thetao.mon.${MODEL}.rcp85.r1i1p1.200601-230012.nc"

LLBOX="0,360,-70,-52"

# Intermediate products file names
thetao="${vol_calculation}/thetao.mon.${MODEL}.historical+rcp85.around_Antarctica.185001-230012.nc"
antarctic_cellarea="${vol_calculation}/cellarea.${MODEL}_around_Antarctica.nc"
int_thetao="${vol_calculation}/interpolated_thetao.mon.${MODEL}.historical+rcp85.r1i1p1.around_Antarctica.185001-230012.nc"
thetao_anom="${vol_calculation}/anom.interpolated_thetao.mon.${MODEL}.historical+rcp85.r1i1p1.around_Antarctica.185001-230012.nc"
grid_depth="${vol_calculation}/depth_at_each_grid_box_${MODEL}.nc"
total_depth="${vol_calculation}/total_depth_at_each_grid_column_${MODEL}.nc"
grid_volume="${vol_calculation}/volume_at_each_grid_box_${MODEL}.nc"
total_thetao="${vol_calculation}/total_thetao_${MODEL}.nc"
total_volume="${vol_calculation}/total_volume_${MODEL}.nc"
lat_lon_dim_removed="${vol_calculation}/thetao.mon.anom.${MODEL}.no_lat-lon-dim.nc"
volume_avg_thetao="${vol_calculation}/thetao.mon.anom.${MODEL}.with-lat-lon-dim.nc"
volume_avg_thetao_anom="thetao.year.anom.${MODEL}.historical_rcp85.r1i1p1.volume-averaged.around_Antarctica.185001-230012.nc"

# select thetao and cellarea in the LLBOX
cdo -a sellonlatbox,$LLBOX ${thetao_hist} ${vol_calculation}/thetao.mon.${MODEL}.historical.r1i1p1.around_Antarctica.nc
cdo -a sellonlatbox,$LLBOX ${thetao_rcp85} ${vol_calculation}/thetao.mon.${MODEL}.rcp85.r1i1p1.around_Antarctica.nc

cdo -a mergetime ${vol_calculation}/thetao.mon.${MODEL}.historical.r1i1p1.around_Antarctica.nc ${vol_calculation}/thetao.mon.${MODEL}.rcp85.r1i1p1.around_Antarctica.nc ${thetao}

cdo -a sellonlatbox,$LLBOX ${cellarea} ${antarctic_cellarea}

# interpolate thetao for the desired depths
cdo -a intlevel,250,350,450,550,650,750 ${thetao} ${int_thetao}

# Calculate anomalies based on 1961 to 1990 climatology
cdo -a ymonsub ${int_thetao} -ymonmean -seldate,1961-01-01,1990-12-31 ${int_thetao} ${thetao_anom}

# Calculate bathimetry based volume for each water column
# Each layer after interpolation is 100 m apart.
# Subtracting the interpolated thetao at 100 m apart layers with itself, leaves N/A 
# values as it is and turns all other values to 0. Add 100 to it, and verically add 
# to calculate total depth at each grid.
# Calculate volume at each grid box to weight the temperatures
cdo -a -addc,100 -sub ${int_thetao} ${int_thetao} ${grid_depth}
cdo -a mul ${grid_depth} ${antarctic_cellarea} ${grid_volume}

# calculate volume weighted thetao at each grid box by multipying interpolated thetao
# at each grid box and volume of each grid box and then sum it over the region around
# Antarctica between depths of 200m and 800m.
cdo -a -fldsum -vertsum -mul ${thetao_anom} ${grid_volume} ${total_thetao}

# Calculate total depth of each grid column
#cdo -a -vertsum -addc,100 -sub ${int_thetao} ${int_thetao} ${total_depth}
cdo -a -vertsum ${grid_depth} ${total_depth} 

# Multiply total depth by cell area to calculate volume at each grid
# Then sum volume of each grid column to find the total ocean volume
cdo -a -fldsum -mul ${total_depth} ${antarctic_cellarea} ${total_volume}

# vertically add interpolated thetao values at all layers between 200m and 800m 
# and then sum all the thetao values to get total thetao between 200 and 800m for each year 
#cdo -a -fldsum -vertsum ${thetao_anom} ${total_thetao}

# divide total thetao by total volume and then calculate yearly mean
cdo -a -yearmean -div ${total_thetao} ${total_volume} ${volume_avg_thetao}

# Remove lat, lon and bnds dimentions from the file
ncwa -a lat,lon,bnds ${volume_avg_thetao} ${lat_lon_dim_removed}

# Remove lat,lon and time_bnds variables from the file
ncks -O -x -v lat,lon,time_bnds ${lat_lon_dim_removed} ${volume_avg_thetao_anom}

# change attributes for thetao
ncatted -O -a standard_name,thetao,o,c,"volume_averaged_sea_water_potential_temperature" ${volume_avg_thetao_anom}
ncatted -O -a long_name,thetao,o,c,"Volume Averaged Sea Water Potential Temperature Anomaly" ${volume_avg_thetao_anom}
