#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
-----------------------------------------------------------------
@purpose : Run rttov 13.4 with realistic profiles 
@author  : Vito Galligani
@email   : vito.galligani@gmail.com
@condaenv: conda activate RTTOV14
-----------------------------------------------------------------
@main    : Runs clear-sky RTTOV v13.4 with a few IFS profiles
           using run_example_fwd_VITO.sh
          
           remember use for reference: https://github.com/lkugler/RTTOV-WRF/blob/master/rttov_wrf.py
           and also to open remote session:
               ssh -L 6000:yakaira:22 vito.galligani@portal.cima.fcen.uba.ar
               nohup python -m spyder_kernels.console -matplotlib='inline' -f=./remotemachine.json &
        
-----------------------------------------------------------------
"""
#------------------------------------------------------------------------------
# to ensure that pyrttov is importable
import sys
#import pyrttov
#import glob
import numpy as np
from netCDF4 import Dataset
import os 
import typhon  as ty
import matplotlib.pyplot as plt

#------------------------------------------------------------------------------
# This script reads nc file of IFS profiles
def ReadIFS(ncfile): 
    

    # Open ncfile
    data = Dataset(ncfile,'r')   
    
    # Hydrometeor content
    A = {}

    A['swc']  = data.variables['snow'][:]      # kg/m3
    A['iwc']  = data.variables['ice'][:]       # kg/m3
    A['gwc']  = data.variables['graupel'][:]   # kg/m3
    A['rwc']  = data.variables['rain'][:]      # kg/m3
    A['cc']   = data.variables['cc'][:]        # cloud cover
    
    A['h20']  = data.variables['q'][:]         # kg/kg 
    A['o3']   = data.variables['o3'][:]        # kg/kg
    
    A['pf']   = data.variables['pressure_full'][:]   # Pa, n = 137
    A['ph']   = data.variables['pressure_half'][:]   # Pa, n = 138
    A['t']    = data.variables['t'][:]               # temperature
    A['tsfc'] = data.variables['tsfc'][:]            # surface temperature
    A['z0']   = data.variables['orography'][:]       # orography
    
    print( 'IFS variables in data: ' )
    print( data.variables.keys() )
    print( '-------------------------------------------------------' )
    return A  

#------------------------------------------------------------------------------
# This script outputs a specific profile of the IFS global profiles
def SelectProf(A,iprof): 

    Ai = {}
    F1 = list(A.keys())
    for k in range(len(F1)):
        fname = F1[k]
        if len(A[fname].shape) == 1:
            Ai[fname] = A[fname][iprof]
        else:
            Ai[fname] = A[fname][iprof,:]

    return Ai  

#------------------------------------------------------------------------------
# This script crate the input clear-sky and all-sky profiles for rttov v13.4
def rttov13_prof_realistic(Ai, fid, fid_clear): 
        
    # vmr q units
    mh2o       = 18.01528
    mair       = 28.9644
    mr2vmr_h2o = mair/mh2o

    # Pressure levels (hPa): space ---> surface in FULL LEVELS n=137
    full_p = Ai['pf']                   # 'Full levels' n = 137     
    full_t = Ai['t']                    # Temperature 
    full_q = mr2vmr_h2o * Ai['h20']     # Water vapor # kg/kg  == D.vmr_field.data(3,:)  
    half_p = Ai['ph']                   # 'half levels' nlevels = n+1     
        
    # We need to convert to kg/kg
    rho_air = ty.physics.density(full_p, full_t)

    # Surface 
    rte_pos = np.abs( Ai['z0']/1e3 ) 
    t_sfc = Ai['tsfc']
    
    #--------------------------------------------------------------------------
    # Hydrometeor info: (qx_mass_grid[:,0] == CLW) in full levels
    qx_mass_full      = np.zeros(( len(full_p), 5))
    qx_mass_full[:,1] = Ai['iwc'] /rho_air  # iwc kg/kg
    qx_mass_full[:,2] = Ai['rwc'] /rho_air  # rwc kg/kg
    qx_mass_full[:,3] = Ai['swc'] /rho_air  # swc kg/kg
    qx_mass_full[:,4] = Ai['gwc'] /rho_air  # gwc kg/kg

    #--------------------------------------------------------------------------
    # For cloud fraction    
    #id1 = qx_mass_full[:,0]>0   #clw is 0
    id2 = qx_mass_full[:,1]>0
    id3 = qx_mass_full[:,2]>0
    id4 = qx_mass_full[:,3]>0
    id5 = qx_mass_full[:,4]>0

    indices_id2 = np.where(id2 > 0)[0]
    indices_id3 = np.where(id3 > 0)[0]
    indices_id4 = np.where(id4 > 0)[0]
    indices_id5 = np.where(id5 > 0)[0]

    # Concatenate the indices (lwc is zero)
    cc = np.concatenate((indices_id2, indices_id3, indices_id4, indices_id5))
    cc = np.unique(cc)

    cf = np.zeros(( len(full_p) )) 
    cf[cc] = 1
    
    #--------------------------------------------------------------------------    
    # Separate profiles
    fid.write('! --- Profile ---\n')
    
    # Write size of profile    
    fid.write('! Size of profile\n')
    fid.write(f'{len(full_p)}\n')
    
    fid.write('! Gas units (must be same for all profiles)\n')
    fid.write('2\n')                                # Writing integer value

    # Write pressure header
    fid.write('! Full pressure (hPa)    Half-pressure (hPa)    Temp.      full_q      cloud fraction      clw / iwc / rwc / swc / gwc \n')
    
    # Loop over each element in 'full_p' (assuming it's a list or NumPy array)
    for nn in range(len(full_p)):
        # Extract values for the current iteration
        full_p_nn       = full_p[nn]
        half_p_nn       = half_p[nn]
        full_t_nn       = full_t[nn]
        full_q_nn       = full_q[nn]
        cf_nn           = cf[nn]
        qx_mass_full_nn = qx_mass_full[nn]
    
        # Write formatted data to filefull_qn
        fid.write(f'{full_p_nn/100:.6f} {half_p_nn/100:.6f} {full_t_nn:.10f} {full_q_nn*1e6:.10f} '
                  f'{cf_nn} {qx_mass_full_nn[0]:.10f} {qx_mass_full_nn[1]:.10f} {qx_mass_full_nn[2]:.10f} '
                  f'{qx_mass_full_nn[3]:.10f} {qx_mass_full_nn[4]:.10f}\n')
    
    # Write near-surface variables header
    fid.write('! Near-surface variables:\n')
    fid.write('! 2m T (K)    2m q (ppmv) 2m p (hPa) 10m wind u (m/s)  10m wind v (m/s)\n')
    
    # Write near-surface variables data
    fid.write(f'{full_t[-1]:.2f}   {full_q[-1]*1e6:.2f}   {half_p[-1]/100:.2f}   0.0   0.0\n')
    
    # Write skin variables header
    fid.write('! Skin variables:\n')
    fid.write('! Skin T (K)  FASTEM parameters for land surfaces\n')
    
    # Write skin variables data
    fid.write(f'{t_sfc:.2f}  3.0  5.0  15.0  0.1  0.3\n')
    
    # Write elevation, latitude, and longitude header
    fid.write('! Elevation (km), latitude and longitude (degrees)\n')
    fid.write(f'{rte_pos:.3f}\n')
    
    # Write zenith and azimuth angles header
    #fid.write('! Sat. zenith and azimuth angles, solar zenith and azimuth angles (degrees)\n')
    #fid.write(f'{zenith:.2f}   0.0\n')
    
    # Write end of profile marker
    fid.write('!\n')
    fid.write('! --- End of profile ---\n')
    
    #--------------------------------------------------------------------------    
    # Save profile info to file
    # Separate profiles
    fid_clear.write('! --- Profile ---\n')
    
    # Write size of profile    
    fid_clear.write('! Size of profile\n')
    fid_clear.write(f'{len(full_p)}\n')

    # Write gas units
    fid_clear.write('! Gas units\n') 
    fid_clear.write(f'{2}\n')  # Writing integer value (assuming gas units is 2)

    # Write pressure
    fid_clear.write('! Pressure levels (hPa)\n')     
    # Loop over each element in 'full_p' (assuming it's a list or NumPy array) 
    for nn in range(len(full_p)): 
        # Extract pressure level value for the current iteration
        pressure_hPa = full_p[nn] / 100.0   # Convert pressure to hPa
        fid_clear.write(f'{pressure_hPa:.6f}\n')  # Writing pressure level in hPa format with 6 decimal places

    # Write temp
    fid_clear.write('! Temperature levels (K)\n')
    # Loop over each element in 'full_t' (assuming it's a list or NumPy array)
    for temperature in full_t:
        fid_clear.write(f'{temperature:.6f}\n')  # Writing temperature level in Kelvin with 6 decimal places
   
    # Write water vapor profiles header
    fid.write('! Water vapour profiles (ppmv)\n')
    # Loop over each element in 'full_q' (assuming it's a list or NumPy array)
    for water_vapor in full_q:
        # Convert water vapor level to ppmv and write to file
        water_vapor_ppmv = water_vapor * 1e6  # Convert to ppmv
        fid_clear.write(f'{water_vapor_ppmv:.6f}\n')  # Writing water vapor level in ppmv with 6 decimal places

    
    # Write near-surface variables header
    fid_clear.write('! Near-surface variables:\n')
    
    # Write near-surface variable descriptions
    fid_clear.write('! 2m T (K)    2m q (ppmv) 2m p (hPa) 10m wind u (m/s)  10m wind v (m/s)\n')
    
    # Write near-surface variables data
    fid_clear.write(f'{full_t[-1]:.2f}   {full_q[-1]*1e6:.2f}   {half_p[-1]/100:.2f}   0.00   0.00\n')
    
    # Write skin variables header 
    fid_clear.write('! Skin variables:\n')
    fid_clear.write('! Skin T (K)  FASTEM parameters for land surfaces\n')

    # Write skin variables data 
    fid_clear.write(f'{t_sfc:.2f}  3.00  5.00  15.00  0.10  0.30\n')
    
    # Write elevation, latitude, and longitude header
    fid_clear.write('! Elevation (km), latitude and longitude (degrees)\n')
    fid_clear.write(f'{rte_pos:.3f}\n')
    
    # Write zenith and azimuth angles header
    #fid.write('! Sat. zenith and azimuth angles, solar zenith and azimuth angles (degrees)\n')
    #fid.write(f'{zenith:.2f}   0.0\n')
    
    # Write end of profile marker
    fid_clear.write('!\n')
    fid_clear.write('! --- End of profile ---\n')   
    
    return 

#------------------------------------------------------------------------------


# This script crate the input clear-sky and all-sky profiles for rttov v4
def rttov14_prof_realistic(Ai, fid, fid_hydro): 
        
    # vmr q units
    mh2o       = 18.01528
    mair       = 28.9644
    mr2vmr_h2o = mair/mh2o

    # Pressure levels (hPa): space ---> surface in FULL LEVELS n=137
    full_p = Ai['pf']                   # 'Full levels' n = 137     
    full_t = Ai['t']                    # Temperature 
    full_q = mr2vmr_h2o * Ai['h20']     # Water vapor # kg/kg  == D.vmr_field.data(3,:)  
    half_p = Ai['ph']                   # 'half levels' nlevels = n+1     
    
    # We need to convert to kg/kg
    rho_air = ty.physics.density(full_p, full_t)

    # Surface 
    rte_pos = np.abs( Ai['z0']/1e3 ) 
    t_sfc = Ai['tsfc']
    
    #--------------------------------------------------------------------------
    # Hydrometeor info: (qx_mass_grid[:,0] == CLW) in full levels
    qx_mass_full      = np.zeros(( len(full_p), 5))
    qx_mass_full[:,1] = Ai['iwc'] /rho_air  # iwc kg/kg
    qx_mass_full[:,2] = Ai['rwc'] /rho_air  # rwc kg/kg
    qx_mass_full[:,3] = Ai['swc'] /rho_air  # swc kg/kg
    qx_mass_full[:,4] = Ai['gwc'] /rho_air  # gwc kg/kg

    #--------------------------------------------------------------------------
    # For cloud fraction    
    #id1 = qx_mass_full[:,0]>0   #clw is 0
    id2 = qx_mass_full[:,1]>0
    id3 = qx_mass_full[:,2]>0
    id4 = qx_mass_full[:,3]>0
    id5 = qx_mass_full[:,4]>0

    indices_id2 = np.where(id2 > 0)[0]
    indices_id3 = np.where(id3 > 0)[0]
    indices_id4 = np.where(id4 > 0)[0]
    indices_id5 = np.where(id5 > 0)[0]

    # Concatenate the indices (lwc is zero)
    cc = np.concatenate((indices_id2, indices_id3, indices_id4, indices_id5))
    cc = np.unique(cc)

    cf = np.zeros(( len(full_p) )) 
    cf[cc] = 1
    
    #--------------------------------------------------------------------------    
    # Construct hydro file
        
    # Save profile info to file
    # Separate profiles
    fid_hydro.write('! --- Profile ---\n')
    
    # Write cloud fraction and hydromteor concentration profiles
    fid_hydro.write('! cloud fraction and hydromteor concentration profiles \n')     
    # Loop over each element in 'full_p' (assuming it's a list or NumPy array) 
    for nn in range(len(full_p)): 
        # Extract eachlevel value 
        fid_hydro.write(f'{cf[nn]:.6f}\n')  # 6 decimal places
        qx_mass_full_nn = qx_mass_full[nn]        
        fid_hydro.write(f'{qx_mass_full_nn[0]:.10f}\n')  # 6 decimal places
        fid_hydro.write(f'{qx_mass_full_nn[1]:.10f}\n')  # 6 decimal places
        fid_hydro.write(f'{qx_mass_full_nn[2]:.10f}\n')  # 6 decimal places
        fid_hydro.write(f'{qx_mass_full_nn[3]:.10f}\n')  # 6 decimal places
        fid_hydro.write(f'{qx_mass_full_nn[4]:.10f}\n')  # 6 decimal places
       
    # Write end of profile marker
    fid_hydro.write('!\n')
    fid_hydro.write('! --- End of profile ---\n')   
           
    
    # Clear-Sky profiles
    fid.write('! --- Profile ---\n')
    
    # Write gas units
    fid.write('! Gas units\n') 
    fid.write(f'{2}\n')  # Writing integer value (assuming gas units is 2)

    # Write pressure
    fid.write('! half pressure levels (hPa)\n')     
    # Loop over each element in 'full_p' (assuming it's a list or NumPy array) 
    for nn in range(len(half_p)): 
        # Extract pressure level value for the current iteration
        pressure_hPa = half_p[nn] / 100.0   # Convert pressure to hPa
        fid.write(f'{pressure_hPa:.6f}\n')  # Writing pressure level in hPa format with 6 decimal places
    
    # Write temp
    fid.write('! Temperature levels (K)\n')
    # Loop over each element in 'full_t' (assuming it's a list or NumPy array)
    for temperature in full_t:
        fid.write(f'{temperature:.6f}\n')  # Writing temperature level in Kelvin with 6 decimal places
   
    # Write water vapor profiles header
    fid.write('! Water vapour profiles (ppmv)\n')
    # Loop over each element in 'full_q' (assuming it's a list or NumPy array)
    for water_vapor in full_q:
        # Convert water vapor level to ppmv and write to file
        water_vapor_ppmv = water_vapor * 1e6  # Convert to ppmv
        fid.write(f'{water_vapor_ppmv:.6f}\n')  # Writing water vapor level in ppmv with 6 decimal places

    # Write near-surface variables header
    fid.write('! Near-surface variables:\n')
    
    # Write near-surface variable descriptions (no longer readl 2mp as == lowest half pressure level)
    fid.write('! 2m T (K)    2m q (ppmv)   10m wind u (m/s)      10m wind v (m/s)\n')
    
    # Write near-surface variables data
    fid.write(f'{full_t[-1]:.2f}   {full_q[-1]*1e6:.2f}   0.00   0.00\n')
    
    # Write skin variables header 
    fid.write('! Skin variables:\n')
    fid.write('! Skin T (K)  FASTEM parameters for land surfaces\n')

    # Write skin variables data 
    fid.write(f'{t_sfc:.2f}  3.00  5.00  15.00  0.10  0.30\n')
    
    # Write elevation, latitude, and longitude header
    fid.write('! Elevation (km), latitude and longitude (degrees)\n')
    fid.write(f'{rte_pos:.3f}\n')
    
    # Write zenith and azimuth angles header
    #fid.write('! Sat. zenith and azimuth angles, solar zenith and azimuth angles (degrees)\n')
    #fid.write(f'{zenith:.2f}   0.0\n')
    
    # Write end of profile marker
    fid.write('!\n')
    fid.write('! --- End of profile ---\n')   
    

    return 

#------------------------------------------------------------------------------
# This is the main rttov 13.4 version set-up 
def run_IFS(ipnd):

    paths_rttov = '/home/vito.galligani/Work/RTTOV/rttov13.2/'
    sys.path.append(paths_rttov+'/wrapper')     
    
    # flag-name (rttov13 and IFS_W1)
    ncfolder  = '/home/vito.galligani/datosmunin3/RTTOV_betaTest/data/'
    ncfile    = ncfolder+'AIRS_processed_W1.nc'
    outfilename = 'rttov132_1File_' +str( ncfile[-5:-3] )
    
    # rttov main output path 
    general_out = '/home/vito.galligani/datosmunin3/RTTOV_betaTest/rttov13/'
         
    # Load all profiles
    A = ReadIFS(ncfile)
    
    print('-------------------------------------------------------------------')    

    # Open outputfilename All-sky
    surface_emissivity = 1
    filename           = f"atm_em{surface_emissivity:.1f}{outfilename}.dat"    
    print('filename: '+filename)
    file_path = os.path.join(general_out, f"InputAll/1FileProfs/TestAllSky/{filename}")   
    
    # Open outputfilename Clear-sky
    # %- Correction for clear sky, the last level was the same as at the surface
    filename   = f"atm_em{surface_emissivity:.1f}{outfilename}.dat"    
    print('filename: '+filename)
    file_path_clear = os.path.join(general_out, f"InputAll/1FileProfs/TestClearSky/{filename}")    
    
    
    # Open file for writing
    with open(file_path, 'w') as fid:
        fid.write('! Number of profiles\n')
        fid.write(f'{22000}\n')
        with open(file_path_clear, 'w') as fid_cl:
            fid_cl.write('! Number of profiles\n')
            fid_cl.write(f'{22000}\n')
            # Write metadata to file
            for iprof in range(22000):
                print('iprof: '+ str(iprof))
                Ai    = SelectProf(A,iprof)
                
                # In case there is ice hydrometeor above 290 K, set the content to zero from 285 to 280
                Ai['swc'][Ai['t']>280]=0 
                Ai['iwc'][Ai['t']>280]=0 
                Ai['gwc'][Ai['t']>280]=0
          
                rttov13_prof_realistic(Ai, fid, fid_cl)
        
                print('----------------------------------------------------------- iprof: '+ str(iprof))
            fid_cl.close()
        fid.close()
    
    return A

#------------------------------------------------------------------------------
# This is the main rttov 13.4 version set-up 
def run_IFS_rttov14(ipnd):
    
    paths_rttov = '/home/vito.galligani/Work/RTTOV/rttov14.0_beta'
    sys.path.append(paths_rttov+'/wrapper')     

    # flag-name (rttov13 and IFS_W1)
    ncfolder  = '/home/vito.galligani/datosmunin3/RTTOV_betaTest/data/'
    ncfile    = ncfolder+'AIRS_processed_W1.nc'
    outfilename = 'rttov14_1File_' +str( ncfile[-5:-3] )
    
    # rttov main output path 
    general_path = '/home/vito.galligani/datosmunin3/RTTOV_betaTest/rttov14'
    A = ReadIFS(ncfile)
    
    print('-------------------------------------------------------------------')    

    # Open outputfilename All-sky
    surface_emissivity = 1
    filename           = f"atm_em{surface_emissivity:.1f}{outfilename}.dat"    
    print('filename: '+filename)
    file_path = os.path.join(general_path, f"InputAll/1FileProfs/{filename}")   
    
    # Open outputfilename Clear-sky
    # %- Correction for clear sky, the last level was the same as at the surface
    filename   = f"atm_em{surface_emissivity:.1f}{outfilename}_HYDRO.dat"    
    print('filename: '+filename)
    file_path_hydro = os.path.join(general_path, f"InputAll/1FileProfs/{filename}")    
    
    # Open file for writing
    with open(file_path, 'w') as fid:
        fid.write('! Number of profiles\n')
        fid.write(f'{22000}\n')
        with open(file_path_hydro, 'w') as fid_hydro:
            for iprof in range(22000):
                
                Ai    = SelectProf(A,iprof)
                # In case there is ice hydrometeor above 290 K, set the content to zero from 285 to 280
                Ai['swc'][Ai['t']>280]=0 
                Ai['iwc'][Ai['t']>280]=0 
                Ai['gwc'][Ai['t']>280]=0
          
                rttov14_prof_realistic(Ai, fid, fid_hydro)        

                print('----------------------------------------------------------- iprof: '+ str(iprof))
            fid_hydro.close()
        fid.close()

    return A

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
# RTTOV 13:
#A = run_IFS(ipnd=1)

# RTTOV 14:
run_IFS_rttov14(ipnd=1)

