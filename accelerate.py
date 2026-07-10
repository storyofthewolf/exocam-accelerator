#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# accelerate.py
# 
# AUTHOR: Wolf, E.T
#
# PURPOSE Asynchronously  advance the restart files trends to speed up convergence.
#         Use with caution with attention to time-evolving trends (see trend.py).
#
# situation 1: (option --aqua_ice)
#     Snowball simulations and cold tidally locked simulations can take hundreds
#     of years to reach true equilibrium due to slow cooling and continued ice sheet
#     growth.  Here we asynchronously increase the ice sheet volume and energy.
#
# Use with care and caution!  Model behavior is not always straight forward.
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

import netCDF4 as nc
import numpy   as np
import os
import sys
import exocampy_tools  as exo
import trend_utils as trend
import argparse
import sys
import subprocess

parser = argparse.ArgumentParser()
parser.add_argument('case_id'       , type=str,     nargs=1, default=' ',  help='Set simulation time series case name')
parser.add_argument('date'          , type=str,     nargs=1, default='0001-01-01', help='Date of restart files')

# modes of operation
parser.add_argument('--aqua_ice',      action='store_true', help='option to accelerate ice sheet energy')
parser.add_argument('--land_temp',     action='store_true', help='option to accelerate land surface temperatures')
parser.add_argument('--land_hydro',    action='store_true', help='option to accelerate land hydrology')
parser.add_argument('--atm_temp',      action='store_true', help='option to accelerate atmosphere temperatures')

# scaling factors
# if these start to get more complicated, move these to an input file
parser.add_argument('-fice'          , type=float,   default='1.0', help='scaling factor')
parser.add_argument('-flnd'          , type=float,   default='1.0', help='scaling factor')
parser.add_argument('-fatm'          , type=float,   default='1.0', help='scaling factor')

args = parser.parse_args()

case_id     = str(args.case_id[0])
date        = str(args.date[0])
path        = "/discover/nobackup/etwolf/cesm_scratch/rundir/"
root        = path + case_id + "/run/"
tsec        = "-00000"

# constrcut restart file names
file_cicer  = root + case_id  + '.cice.r.'  + date + tsec + '.nc' 
file_camr   = root + case_id  + '.cam.r.'   + date + tsec + '.nc' 
file_camrs  = root + case_id  + '.cam.rs.'  + date + tsec + '.nc' 
file_docnr  = root + case_id  + '.docn.r.'  + date + tsec + '.nc'
file_clmr   = root + case_id  + '.clm2.r.'  + date + tsec + '.nc' 

file_cicer_backup  = root + case_id  + '.cice.r.'  + date + tsec + '_backup.nc' 
file_camr_backup   = root + case_id  + '.cam.r.'   + date + tsec + '_backup.nc' 
file_camrs_backup  = root + case_id  + '.cam.rs.'  + date + tsec + '_backup.nc' 
file_docnr_backup  = root + case_id  + '.docn.r.'  + date + tsec + '_backup.nc' 
file_clmr_backup   = root + case_id  + '.clm2.r.'  + date + tsec + '_backup.nc' 

# rpointer files
rp_atm = "rpointer.atm"
rp_drv = "rpointer.drv"
rp_ice = "rpointer.ice"
rp_ocn = "rpointer.ocn"
rp_lnd = "rpointer.lnd"


if (args.aqua_ice == True):
    print("------- aqua sea ice energy acceleration -------")
    print("   scaling factor: ", args.fice)
    file_present = os.path.isfile(file_cicer_backup)
    if file_present == False:
        print("   no backup file exists for ",file_cicer)
        print("   creating backup file      ",file_cicer_backup)
        command = ['cp', file_cicer, file_cicer_backup]
        subprocess.run(command)
        print("------------------------------------------------")
    if file_present == True:
        print("   back up file exists ",file_cicer_backup)
        print("   copying to          ",file_cicer)
        command = ['cp', file_cicer_backup, file_cicer]      
        subprocess.run(command)
        print("------------------------------------------------")

    # open file for read in
    # read from backup file
    ncid = nc.Dataset(file_cicer_backup, 'r')
    Tsfcn_in = ncid['Tsfcn'][:]
    eicen_in = ncid['eicen'][:]
    esnon_in = ncid['esnon'][:]
    vicen_in = ncid['vicen'][:]
    vsnon_in = ncid['vsnon'][:]
    ncid.close()

    # scale sea ice and snow
    eicen_out = np.zeros((eicen_in.shape[:]), dtype=float)
    esnon_out = np.zeros((esnon_in.shape[:]), dtype=float)
    vicen_out = np.zeros((vicen_in.shape[:]), dtype=float)
    vsnon_out = np.zeros((vsnon_in.shape[:]), dtype=float)

    eicen_out[:,:,:] = eicen_in[:,:,:] * args.fice
    esnon_out[:,:,:] = esnon_in[:,:,:] * args.fice
    vicen_out[:,:,:] = vicen_in[:,:,:] * args.fice
    vsnon_out[:,:,:] = vsnon_in[:,:,:] * args.fice
   
    # open file for writing
    ncid = nc.Dataset(file_cicer, 'r+')
    eicen = ncid.variables['eicen']
    eicen[:,:,:] = eicen_out[:,:,:] 
    vicen = ncid.variables['vicen']
    vicen[:,:,:] = vicen_out[:,:,:] 
    esnon = ncid.variables['esnon']
    esnon[:,:,:] = esnon_out[:,:,:] 
    vsnon = ncid.variables['vsnon']
    vsnon[:,:,:] = vsnon_out[:,:,:] 
    ncid.close()

    print("   cice.r. sea ice and snow scaling complete")


# Future work
# There is no "T" or "TS" in the cam.r file, although it seems to have everything else.  Hmmm....
if (args.atm_temp == True):
    print("------- atmosphere temperature acceleration -------")
    print("   scaling factor: ", args.fatm)
    file_present = os.path.isfile(file_camr_backup)
    if file_present == False:
        print("   no backup file exists for ",file_camr)
        print("   creating backup file      ",file_camr_backup)
        command = ['cp', file_camr, file_camr_backup]
        subprocess.run(command)
        print("------------------------------------------------")
    if file_present == True:
        print("   back up file exists ",file_camr_backup)
        print("   copying to          ",file_camr)
        command = ['cp', file_camr_backup, file_camr]
        subprocess.run(command)
        print("------------------------------------------------")

    # open file for read in
    # read from backup file
    ncid = nc.Dataset(file_camr_backup, 'r')
    T_in  = ncid['T'][:]
    TS_in = ncid['TS'][:]
    ncid.close()

    # scale surface and atmosphere temperatures
    T_out  = np.zeros((T_in.shape[:]), dtype=float)
    TS_out = np.zeros((TS_in.shape[:]), dtype=float)

    T_out[:,:,:]  = T_in[:,:,:] * args.fatm
    TS_out[:,:,:] = TS_in[:,:,:] * args.fatm

    # open file for writing
    ncid = nc.Dataset(file_camr, 'r+')
    T = ncid.variables['T']
    T[:,:,:] = T_out[:,:,:]
    TS = ncid.variables['TS']
    TS[:,:,:] = TS_out[:,:,:]
    ncid.close()

    print("   cam.r. atmosphere and surfaxe temperature scaling complete")




if (args.land_temp == True):
    print("------- land temperature acceleration -------")
    print("   scaling factor: ", args.flnd)
    file_present = os.path.isfile(file_clmr_backup)
    if file_present == False:
        print("   no backup file exists for ",file_clmr)
        print("   creating backup file      ",file_clmr_backup)
        command = ['cp', file_clmr, file_clmr_backup]
        subprocess.run(command)
        print("------------------------------------------------")
    if file_present == True:
        print("   back up file exists ",file_clmr_backup)
        print("   copying to          ",file_clmr)
        command = ['cp', file_clmr_backup, file_clmr]
        subprocess.run(command)
        print("------------------------------------------------")

    # open file for read in                                                                                                                
    # read from backup file                                                                                                                
    ncid         = nc.Dataset(file_clmr_backup, 'r')
    T_GRND_in    = ncid['T_GRND'][:]
    T_GRND_U_in  = ncid['T_GRND_U'][:]
    T_SOISNO_in  = ncid['T_SOISNO'][:] # soil-snow temperature
    T_LAKE_in    = ncid['T_LAKE'][:]   # soil-snow temperature
    T_VEG_in     = ncid['T_VEG'][:]    # vegetation temperature
    T_REF2M_in   = ncid['T_REF2M'][:]  # 2m height surface air temperature
    # there are a number of other 2-m reference air temperature variables
    ncid.close()

    # scale land surface temperature fields
    T_GRND_out   = np.zeros((T_GRND_in.shape[:]), dtype=float)
    T_GRND_U_out = np.zeros((T_GRND_U_in.shape[:]), dtype=float)
    T_SOISNO_out = np.zeros((T_SOISNO_in.shape[:]), dtype=float)
    T_LAKE_out   = np.zeros((T_LAKE_in.shape[:]), dtype=float)
    T_VEG_out    = np.zeros((T_VEG_in.shape[:]), dtype=float)
    T_REF2M_out  = np.zeros((T_REF2M_in.shape[:]), dtype=float)

    T_GRND_out   =  T_GRND_in[:,:,:] * args.flnd
    T_GRND_U_out =  T_GRND_U_in[:,:,:] * args.flnd
    T_SOISNO_out =  T_SOISNO_in[:,:,:] * args.flnd
    T_LAKE_out   =  T_LAKE_in[:,:,:] * args.flnd
    T_VEG_out    =  T_VEG_in[:,:,:] * args.flnd
    T_REF2M_out  =  T_REF2M_in[:,:,:] * args.flnd


    # open file for writing                                                                                                                
    ncid = nc.Dataset(file_clmr, 'r+')
    T_GRND = ncid.variables['T_GRND']
    T_GRND[:,:,:] = T_GRND_out[:,:,:]
    T_GRND_U = ncid.variables['T_GRND_U']
    T_GRND_U[:,:,:] = T_GRND_U_out[:,:,:]
    T_SOISNO = ncid.variables['T_SOISNO']
    T_SOISNO[:,:,:] = T_SOISNO_out[:,:,:]
    T_LAKE = ncid.variables['T_LAKE']
    T_LAKE[:,:,:] = T_LAKE_out[:,:,:]
    T_VEG = ncid.variables['T_VEG']
    T_VEG[:,:,:] = T_VEG_out[:,:,:]
    T_REF2M = ncid.variables['T_REF2M']
    T_REF2M[:,:,:] = T_REF2M_out[:,:,:]
    ncid.close()

    print("   clm.r. surface temperature scaling complete")



if (args.land_hydro == True):
    print('land hydrology acceleration not currently set up')
    pass



