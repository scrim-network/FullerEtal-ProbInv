#!/bin/bash

# Create a 1 degree x 1 degree land mask for Antarctica with missing values
# over the ocean and an add altitude adjustment for the reduction of surface air
# temperature to sea level based on a dry adiabatic lapse rate of 9.8 K/km.
# 
# Produces surface air temperature on and around Antarctica that is adjusted for 
# altitude and dry adiabatic lapse rate.
#
# For land masks and elevation fields at lots of different resolutions, see
# http://jisao.washington.edu/datasets/elevation/.
#
# Modified Fri 14 Nov 2014 by Robert Nicholas.
# Modified Tue 21 Apr 2015 by Varada Vaidya
# Last modified Wed 22 Apr 2015 by Robert Nicholas.
# -----------------------------
#
# In CDO The weights are exclusive computed from the grid cell area. The netCDF 
# variable for the cell area is defined by the attribute cell_measures:
# var:cell_measures = "area: cell_area" ;
# If the cell area is not available then it will be computed from the geographical 
# coordinates via spherical triangles. This is only possible if the geographical 
# coordinates of the grid cell corners are available or derivable. Otherwise 
# fldmean gives a warning message and uses constant area weights. You can 
# overwrite/set the grid cell area with the undocumented CDO function setgridarea. 
# The parameter for this function is a data file with one field. This field must 
# have the same grid size as the input file:
# cdo -fldmean -setgridarea,gridareafile  ifile ofile
# -----------------------------

MODEL="CNRM-CM5"

# Read files for elevation, temperature, cell area and land area fraction. 
ELEVATION="/woju/s0/common/ethz_rsync/cmip5/historical/fx/orog/${MODEL}/r0i0p0/orog_fx_${MODEL}_historical_r0i0p0.nc"
TEMPERATURE="/woju/s0/data_archive/CMIP5/devel_CMIP5/historical+rcp85_merged/tas.mon.${MODEL}.historical+rcp85.r1i1p1.185001-230012.nc"
CELL_AREA="/woju/s0/common/ethz_rsync/cmip5/historical/fx/areacella/${MODEL}/r0i0p0/areacella_fx_${MODEL}_historical_r0i0p0.nc"
LAND_AREA_FRACTION="/woju/s0/common/ethz_rsync/cmip5/historical/fx/sftlf/${MODEL}/r0i0p0/sftlf_fx_${MODEL}_historical_r0i0p0.nc"

# Specify output file name.
OUTFILE_NAME="tas.mon_anom.${MODEL}.historical+rcp85.r1i1p1.Antarctica_mean_with_adjustment_to_sea_level.185001-230012.nc"


# function Generate_Antarctica_Adjusted_Mean()
#
# This function calculates the sea-level adjusted monthly mean temperature
# anomaly over Antarctica relative to a 1961-1990 climatology.
#
# Arguments passed to the function are as follows:
# 
#   1. name to use for output file $OUTFILE_NAME
#   2. input elevation file $ELEVATION
#   3. input land area fraction file $LAND_AREA_FRACTION
#   4. input surface air temperature file $TEMPERATURE
#   5. input cell area file $CELL_AREA
#
function Generate_Antarctica_Adjusted_Mean {

   TMP="`mktemp -d tmp_XXXXXXXXXX`"

   # Longitude and latitude bounds for mask.
   LLBOX="0,360,-90,-60"
   OUTFILE=$1

   # Chop tepmperature, land area fraction and cell area files for LLBOX
   # coordinates.  Use elevation file for lapse rate adjustment.
   cdo -a -mulc,0.0098 -sellonlatbox,$LLBOX $2 ${TMP}/lra_orog.nc
   cdo -a sellonlatbox,$LLBOX $3 ${TMP}/sftlf.nc
   cdo -a sellonlatbox,$LLBOX $4 ${TMP}/tas.nc
   cdo -a sellonlatbox,$LLBOX $5 ${TMP}/areacella.nc

   # Add lapse rate ajustment to temperature.
   cdo -a  add ${TMP}/tas.nc ${TMP}/lra_orog.nc ${TMP}/tas+lra_orog.nc

   # Calculate anomalies relative to a 1961-1990 climatology
   cdo -a ymonsub ${TMP}/tas+lra_orog.nc -ymonmean -seldate,1961-01-01,1990-12-31 ${TMP}/tas+lra_orog.nc ${TMP}/tas+lra_orog_anom.nc

   # Calculate Ti*Ai*Li and Ai*Li, which can be used as a weight to
   # calculate field mean, where Ti is lapse rate adjusted temperature, Li
   # is land area fraction, and Ai is cell area.
   cdo -a  mul ${TMP}/tas+lra_orog_anom.nc ${TMP}/sftlf.nc ${TMP}/TiLi.nc
   cdo -a  mul ${TMP}/TiLi.nc ${TMP}/areacella.nc ${TMP}/TiAiLi.nc
   cdo -a  mul ${TMP}/sftlf.nc ${TMP}/areacella.nc ${TMP}/AiLi.nc
   cdo -a fldsum ${TMP}/TiAiLi.nc ${TMP}/sum_TiAiLi.nc
   cdo -a fldsum ${TMP}/AiLi.nc ${TMP}/sum_AiLi.nc
   cdo -a div ${TMP}/sum_TiAiLi.nc ${TMP}/sum_AiLi.nc ${OUTFILE}

   rm -rf $TMP
   echo $OUTFILE;
   return
}   
   

output_file=$(Generate_Antarctica_Adjusted_Mean ${OUTFILE_NAME} ${ELEVATION} ${LAND_AREA_FRACTION} ${TEMPERATURE} ${CELL_AREA})
echo $output_file
 
