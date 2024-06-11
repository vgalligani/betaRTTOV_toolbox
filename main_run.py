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
def rttov13_prof_realistic(Ai, A2filename, paths_rttov): 
        
    # vmr q units
    mh2o       = 18.01528
    mair       = 28.9644
    mr2vmr_h2o = mair/mh2o

    # Pressure levels (hPa): space ---> surface in FULL LEVELS n=137
    full_p = Ai['pf']                   # 'Full levels' n = 137     
    full_t = Ai['t']                    # Temperature 
    full_q = mr2vmr_h2o * Ai['h20']     # Water vapor # kg/kg  == D.vmr_field.data(3,:)  
    half_p = Ai['ph']                   # 'half levels' nlevels = n+1     
        
    # Plot the levels 
    # plt.axhline(np.nan, color = 'r', linestyle = '-', label='full levels') 
    # plt.axhline(np.nan, color = 'k', linestyle = '-', label='half levels') 
    # plt.legend(ncol=2)
    # for i in range(len(full_p)):
    #     plt.axhline(y = full_p[i]/100, color = 'r', linestyle = '-') 
    # for i in range(len(half_p)):
    #     plt.axhline(y = half_p[i]/100, color = 'k', linestyle = '-')         
    # plt.ylim([1000, 1035])
    # plt.ylabel('Pressure [hPa]')    
    # plt.title('Pressure levels')    
    
    # We need to convert to kg/kg
    rho_air = ty.physics.density(full_p, full_t)

    # Surface 
    rte_pos = np.abs( Ai['z0']/1e3 ) 
    surface_emissivity = 1
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
    # Construct filename string All-sky
    #ifname          = '{:02d}'.format(zenith)  
    #filename = f"atm_zenith{ifname}em{surface_emissivity:.1f}{A2filename}.dat"
    filename = f"atm_em{surface_emissivity:.1f}{A2filename}.dat"    
    print('filename: '+filename)
    
    # Construct full file path
    file_path = os.path.join(paths_rttov, f"InputAll/TestAllSky/{filename}")   
    
    # Open file for writing
    with open(file_path, 'w') as fid:
        # Write metadata to file
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
        
        # Close the file
        fid.close()
    
    # Clear-Sky profiles
    # %- Correction for clear sky, the last level was the same as at the surface
    print('filename: '+filename)
    file_path = os.path.join(paths_rttov, f"InputAll/TestClearSky/{filename}")
    
    # Open file for writing
    with open(file_path, 'w') as fid:
        
        # Write metadata to file    
        fid.write('! Number of profiles\n')
        fid.write(f'{1}\n')

        # Write size of profile
        fid.write('! Size of profile\n')
        fid.write(f'{len(full_p)}\n')

        # Save profile info to file
        # Separate profiles
        fid.write('! --- Profile ---\n')
        
        # Write gas units
        fid.write('! Gas units\n') 
        fid.write(f'{2}\n')  # Writing integer value (assuming gas units is 2)

        # Write pressure
        fid.write('! Pressure levels (hPa)\n')     
        # Loop over each element in 'full_p' (assuming it's a list or NumPy array) 
        for nn in range(len(full_p)): 
            # Extract pressure level value for the current iteration
            pressure_hPa = full_p[nn] / 100.0   # Convert pressure to hPa
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
        
        # Write near-surface variable descriptions
        fid.write('! 2m T (K)    2m q (ppmv) 2m p (hPa) 10m wind u (m/s)  10m wind v (m/s)\n')
        
        # Write near-surface variables data
        fid.write(f'{full_t[-1]:.2f}   {full_q[-1]*1e6:.2f}   {half_p[-1]/100:.2f}   0.00   0.00\n')
        
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
        
        # Close the file
        fid.close()


    return filename

#------------------------------------------------------------------------------


# This script crate the input clear-sky and all-sky profiles for rttov v4
#def rttov14_prof_realistic(Ai, A2filename, zenith, paths_rttov): 
def rttov14_prof_realistic(Ai, A2filename, paths_rttov): 
        
    # vmr q units
    mh2o       = 18.01528
    mair       = 28.9644
    mr2vmr_h2o = mair/mh2o

    # Pressure levels (hPa): space ---> surface in FULL LEVELS n=137
    full_p = Ai['pf']                   # 'Full levels' n = 137     
    full_t = Ai['t']                    # Temperature 
    full_q = mr2vmr_h2o * Ai['h20']     # Water vapor # kg/kg  == D.vmr_field.data(3,:)  
    half_p = Ai['ph']                   # 'half levels' nlevels = n+1     
        
    # Plot the levels 
    # plt.axhline(np.nan, color = 'r', linestyle = '-', label='full levels') 
    # plt.axhline(np.nan, color = 'k', linestyle = '-', label='half levels') 
    # plt.legend(ncol=2)
    # for i in range(len(full_p)):
    #     plt.axhline(y = full_p[i]/100, color = 'r', linestyle = '-') 
    # for i in range(len(half_p)):
    #     plt.axhline(y = half_p[i]/100, color = 'k', linestyle = '-')         
    # plt.ylim([1000, 1035])
    # plt.ylabel('Pressure [hPa]')    
    # plt.title('Pressure levels')    
    
    # We need to convert to kg/kg
    rho_air = ty.physics.density(full_p, full_t)

    # Surface 
    rte_pos = np.abs( Ai['z0']/1e3 ) 
    surface_emissivity = 1
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
    # Construct filename string All-sky
    #ifname          = '{:02d}'.format(zenith)  
    #filename = f"atm_zenith{ifname}em{surface_emissivity:.1f}{A2filename}.dat"
    #filename_hydro = f"atm_zenith{ifname}em{surface_emissivity:.1f}{A2filename}_hydro.dat"
    filename = f"atm_em{surface_emissivity:.1f}{A2filename}.dat"
    filename_hydro = f"atm_em{surface_emissivity:.1f}{A2filename}_hydro_REAL.dat"
    print('filename: '+filename)
    print('filename: '+filename_hydro)
    
    # # Construct full file path
    file_path_hydro = os.path.join(paths_rttov, f"InputAll/{filename_hydro}")
    with open(file_path_hydro, 'w') as fid:
        
        # Save profile info to file
        # Separate profiles
        fid.write('! --- Profile ---\n')
        
        # Write cloud fraction and hydromteor concentration profiles
        fid.write('! cloud fraction and hydromteor concentration profiles \n')     
        # Loop over each element in 'full_p' (assuming it's a list or NumPy array) 
        for nn in range(len(full_p)): 
            # Extract eachlevel value 
            fid.write(f'{cf[nn]:.6f}\n')  # 6 decimal places
            qx_mass_full_nn = qx_mass_full[nn]        
            fid.write(f'{qx_mass_full_nn[0]:.10f}\n')  # 6 decimal places
            fid.write(f'{qx_mass_full_nn[1]:.10f}\n')  # 6 decimal places
            fid.write(f'{qx_mass_full_nn[2]:.10f}\n')  # 6 decimal places
            fid.write(f'{qx_mass_full_nn[3]:.10f}\n')  # 6 decimal places
            fid.write(f'{qx_mass_full_nn[4]:.10f}\n')  # 6 decimal places
           
        # Write end of profile marker
        fid.write('!\n')
        fid.write('! --- End of profile ---\n')   
        
        # Close the file
        fid.close()            
    
    # Clear-Sky profiles
    # %- Correction for clear sky, the last level was the same as at the surface
    file_path = os.path.join(paths_rttov, f"InputAll/{filename}")
    print('filename: '+file_path)
    
    # Open file for writing
    with open(file_path, 'w') as fid:
        
        # Write metadata to file    
        #fid.write('! Number of profiles\n')
        #fid.write(f'{1}\n')

        # Write size of profile
        #fid.write('! Size of profile\n')
        #fid.write(f'{len(half_p)}\n')

        # Save profile info to file
        # Separate profiles
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
        #print('CHECK READ ALL HALF PRESSURE LEVELS (n): ', str(len(half_p)))
        
        # Write temp
        fid.write('! Temperature levels (K)\n')
        # Loop over each element in 'full_t' (assuming it's a list or NumPy array)
        for temperature in full_t:
            fid.write(f'{temperature:.6f}\n')  # Writing temperature level in Kelvin with 6 decimal places
        #print('CHECK t,q in FULL LEVELS nn (n-1): ', str(len(full_t)))
   
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
        
        # Close the file
        fid.close()


    return filename

#------------------------------------------------------------------------------




# This is the main rttov 13.4 version set-up 
def run_IFS(ipnd):

    paths_rttov = '/home/vito.galligani/Work/RTTOV/rttov13.2/'
    sys.path.append(paths_rttov+'/wrapper')     
    
    # flag-name (rttov13 and IFS_W1)
    ncfolder  = '/home/vito.galligani/datosmunin3/RTTOV_betaTest/data/'
    ncfile    = ncfolder+'AIRS_processed_W1.nc'
    flag_name = 'rttov132_Profile_' +str( ncfile[-5:-3] )
    
    # rttov main output path 
    general_out = '/home/vito.galligani/datosmunin3/RTTOV_betaTest/rttov13/'
    #paths_outfolder  = general_out+flag_name

    # Creates main output folder     
    #if not os.path.exists(paths_outfolder):
    #    os.makedirs(paths_outfolder)
         
    # Load all profiles
    A = ReadIFS(ncfile)
    
    print('-------------------------------------------------------------------')    

    
    for iprof in range(22000):
        print('iprof: '+ str(iprof))
        Ai    = SelectProf(A,iprof)
        
        # In case there is ice hydrometeor above 290 K, set the content to zero from 285 to 280
        Ai['swc'][Ai['t']>280]=0 
        Ai['iwc'][Ai['t']>280]=0 
        Ai['gwc'][Ai['t']>280]=0
  
        # Filenames following flag
        flag2name = flag_name+'_'+str(iprof)  
    
        #- Files containing the input files for clear- (cs) and all-sky (as) simulations
        #fas_path = f"{paths_rttov}/rttov_test/TestAllSky/{flag2name}.txt"
        #fcs_path = f"{paths_rttov}/rttov_test/TestClearSky/{flag2name}.txt"

        #fcs = open(fcs_path, 'w') 
        #fas = open(fas_path, 'w') 
    
        prof_filename = rttov13_prof_realistic(Ai, flag2name, general_out)
        #fas.write(f'{prof_filename} {ipnd:.2f}\n')
        #fcs.write(f'{prof_filename}\n')

        # Close the files
        #fas.close()
        #fcs.close()


        print('-----------------------------------------------------------')

    return A

#------------------------------------------------------------------------------
# This is the main rttov 13.4 version set-up 
def run_IFS_rttov14(ipnd):
    
    paths_rttov = '/home/vito.galligani/Work/RTTOV/rttov14.0_beta'
    sys.path.append(paths_rttov+'/wrapper')     

    # flag-name (rttov13 and IFS_W1)
    ncfolder  = '/home/vito.galligani/datosmunin3/RTTOV_betaTest/data/'
    ncfile    = ncfolder+'AIRS_processed_W1.nc'
    flag_name = 'rttov14_Profile_' +str( ncfile[-5:-3] )
    
    # rttov main output path 
    general_path = '/home/vito.galligani/datosmunin3/RTTOV_betaTest/rttov14'
    #paths_outfolder  = general_out+flag_name

    # Creates main output folder     
    #if not os.path.exists(paths_outfolder):
    #    os.makedirs(paths_outfolder)
         
    # Load all profiles
    A = ReadIFS(ncfile)
    
    print('-------------------------------------------------------------------')    

    zenith = np.arange(0,60,5)
    
    for iprof in range(22000): 
        print('iprof: '+ str(iprof))
        Ai    = SelectProf(A,iprof)
        
        # In case there is ice hydrometeor above 290 K, set the content to zero from 285 to 280
        Ai['swc'][Ai['t']>280]=0 
        Ai['iwc'][Ai['t']>280]=0 
        Ai['gwc'][Ai['t']>280]=0
  
        # Filenames following flag
        flag2name = flag_name+'_'+str(iprof)  
    
        #- Files containing the input files for clear- (cs) and all-sky (as) simulations
        #fas_path = f"{paths_rttov}/rttov_test/TestAllSky/{flag2name}.txt"
        #fcs_path = f"{paths_rttov}/rttov_test/TestClearSky/{flag2name}.txt"
        #fas_path = f"{general_path}/Input/TestAllSky/{flag2name}.txt"
        #fcs_path = f"{general_path}/Input/TestClearSky/{flag2name}.txt"
        
        #fcs = open(fcs_path, 'w') 
        #fas = open(fas_path, 'w') 
    
        #for iz in range(len(zenith)): 
        #prof_filename = rttov14_prof_realistic(Ai, flag2name, zenith[iz], general_path)
        prof_filename = rttov14_prof_realistic(Ai, flag2name, general_path)        
        print(prof_filename)
        #fas.write(f'{prof_filename} {ipnd:.2f}\n')
        #fcs.write(f'{prof_filename}\n')

        # Close the files
        #fas.close()
        #fcs.close()


        print('-----------------------------------------------------------')


    return A

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
# RTTOV 13:
#A = run_IFS(ipnd=1)
# RTTOV 14:
run_IFS_rttov14(ipnd=1)


#import main_run
#A = main_run.run_IFS(sat="ssmis", ipnd=1)