#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
-----------------------------------------------------------------
@purpose : Analyse rttov 13.4 outputs  
@author  : Vito Galligani
@email   : vito.galligani@gmail.com
@condaenv: conda activate RTTOV14
-----------------------------------------------------------------
@main    : First simple rttov 13.2 and 14.0 comparison plots for atms
            - BTs, tramisssions, tau ... 
          
           Remember to open remote session:
               ssh -L 6000:yakaira:22 vito.galligani@portal.cima.fcen.uba.ar
               nohup python -m spyder_kernels.console -matplotlib='inline' -f=./remotemachine.json &
        
-----------------------------------------------------------------
"""
# Load libraries
import glob
import numpy as np
import os
from datetime import datetime
import pandas as pd
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import typhon  as ty
import Plots4Analysis as P4A
import Tools2Read as T2R

#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
def main(sat, Aiprofile):
    
    plt.matplotlib.rc('font', family='serif', size = 12)

    # half pressure levels in IFS == nlevels 
    nlev = 138  
    # T, q, etc in rttov is now nlevels-1  == nlayers 

    # Satellite zenith angle    
    zenith =  np.arange(0,60,5)
    
    # Sensor-related 
    if sat == "ici":
        sat   = 'ici'
        nchan = 13
        print("sat == ici")
    elif sat == "ssmis":
        sat   = 'ssmis'
        nchan = 3
        print("sat == ssmis")
        
    elif sat == "atms":
        sat   = 'atms'
        nchan = 22
        opt = {'nrows':2,'ncols':3,'width':11,'height':8}
        # [CH1: 23.800, CH2:31.400, CH9: 55.500, CH16: 88.2, CH17: 165.5, CH18:183.31±7.0, 
        # CH20: 183.31±3.0, CH22: 183.31±1.0]
        #nchan_plots = [1, 2, 9, 16, 17, 18, 20, 22]
        nchan_plots = [2, 16, 17, 18, 20, 22]
        chantitles = ['23.8', '31.4', '50.3', '51.76', '52.8', 
                      '53.596', '54.4', '54.94','55.5', '57.29', 
                       '57.29$\pm$0.5329', '57.29$\pm$0.3702', '57.29$\pm$0.3442', 
                       '57.29$\pm$0.3322', '57.29$\pm$0.3267', '88.2', '165.5',  
                       '183.31$\pm$7', '183.31$\pm$4.5', '183.31$\pm$3', 
                       '183.31$\pm$1.8', '183.31$\pm$1']
        print("sat == atms")
        
    elif sat == "mhs":
        sat   = 'mhs'
        nchan = 5
        opt = {'nrows':2,'ncols':3,'width':8,'height':8}
        # [CH1: 23.800, CH2:31.400, CH9: 55.500, CH16: 88.2, CH17: 165.5, CH18:183.31±7.0, 
        # CH20: 183.31±3.0, CH22: 183.31±1.0]
        nchan_plots = [1, 2, 3, 4, 5]
        chantitles = ['89', '157', '183.31$\pm$3', '183.31$\pm$1', '190']
        print("sat == mhs")  
        
    # Profiles related profile file names    
    ncfolder  = '/home/vito.galligani/datosmunin3/RTTOV_betaTest/data/'
    ncfile    = 'AIRS_processed_W1.nc'
    data      = Dataset(ncfolder+ncfile,'r')   

    Ap_full   = data.variables['pressure_full'][:]   # Pa, n = 137
    p_full    = Ap_full[Aiprofile]
    t_full    = data.variables['t'][:][Aiprofile]

    plevels = {}
    plevels['rttov14_halflevels'] = data.variables['pressure_half'][:][Aiprofile].copy()
    plevels['rttov13_fulllevels'] = p_full.copy()
        
    # Re-allocate pressure levels between half pressure levels == plevels['rttov13_fulllevels']  
    pressure = np.empty(len(plevels['rttov14_halflevels'])-1)
    for ilayer in range(len(plevels['rttov14_halflevels'])-1):
        pressure[ilayer] = (plevels['rttov14_halflevels'][ilayer] + plevels['rttov14_halflevels'][ilayer+1])/2
    plevels['betwee_halflevels'] = pressure 

    # Re-allocate pressure levels between full pressure levels == plevels['between_fulllevels_136']  
    pressure = np.empty(len(plevels['rttov13_fulllevels'])-1)
    for ilayer in range(len(plevels['rttov13_fulllevels'])-1):
        pressure[ilayer] = (plevels['rttov13_fulllevels'][ilayer] + plevels['rttov13_fulllevels'][ilayer+1])/2
    plevels['between_fulllevels_136'] = pressure 
    
    # Output Plot dir (Within this folder, define the name of a sub-folder according to date)
    plot_dir   = '/home/vito.galligani/Work/RTTOV/Plots/'
    path_plots = os.path.join(plot_dir,datetime.now().strftime('%d%m%Y'))

    # overview plot config
    opt_lev = { 'nrows':opt['nrows'],'ncols':opt['ncols'],'width':opt['width'],
                'height':opt['height'],'path':path_plots+'/'+sat+'/'+str(Aiprofile), 'nlev':nlev, 'nchan':nchan, 
                'nchan_plots':nchan_plots, 'chan_plots_titles':chantitles}
    
    # Plot pressure levels 
    P4A.plot_pressuregrid(plevels, opt_lev)   
    
    
    # READ FILES   
    Tb, Trans = T2R.read_outouts(ncfile, Aiprofile, nlev, nchan, sat)

    #-------------------------------------------------------------------------
    # 1/ tb (brightness temperature)
    P4A.ClearSkySensor(Tb, 'BTs ATMS Clearsky', opt_lev)
    P4A.ClearSkySensor_coeffsComp(Tb,'BTs: Clear-sky Coefficients Test', opt_lev)

    # 2. Transmission [transmission % tau_total: Transmittance from surface to TOA]
    # and tau_level transmission
    # OD tau == - ln(Transmittance)  
    #-------------------------------------------------------------------------
    # print totals
    #for iz in range(len(zenith)): 
    #    print('------- zenith: '+str(zenith[iz]))
    #    print('rttov total rttov 13. transmittance () ', np.round(Trans['rttov13_cs_Surface2spaceTr'][iz,:],3))
    #    print('rttov total rttov 14. transmittance () ', np.round(Trans['rttov14_cs_Surface2spaceTr'][iz,:],3))
    P4A.TransSensor(Trans, 'Surf2TOA Total Transmittance ATMS', opt_lev, plevels['rttov14_halflevels'], plevels['rttov13_fulllevels'])        
    P4A.ODSensor(Trans, plevels, opt_lev)        
    P4A.TransSensor_coeffs(Trans, 'Surf2TOA total Transmittance (Clearsky) Coefficients Test', 
                             opt_lev, plevels['rttov14_halflevels'], plevels['rttov13_fulllevels'])  
    P4A.ODSensor_coeffs(Trans, plevels, opt_lev) 
        
    # # NOW ADD TRANSMSSION IN RTTOV-SCTT RTTOV 13 Y OPTICAL DEPTHS DIFFERENTES
    P4A.NadirOD_angle_RTTOVSCATT(Trans['rttov13_as_extPro_angles'],'rttov13', opt_lev)
    P4A.NadirOD_angle_RTTOVSCATT(Trans['rttov14_aux_extClear'].swapaxes(1, 2), 'rttov14', opt_lev)
    
        
    # # Try to understund how is slant path zenith angle handled in rttov 14 and if numerial issues
    # # were solved when tau is turned back to tau nadir in rttov-scatt 
    P4A.scatt_aux(Trans, plevels, opt_lev)
    P4A.scatt_aux_ext_compare(Trans['rttov13_as_extPro_angles'], Trans['rttov14_aux_extClear'].swapaxes(1, 2), plevels, opt_lev)

    



    return Tb, Trans


#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
def plots_4SingleProfiles(): 
    
    # Get the two profiles of interest for the single profile analysis: 
    ncfolder  = '/home/vito.galligani/datosmunin3/RTTOV_betaTest/data/'
    ncfile    = ncfolder+'AIRS_processed_W1.nc'
    plot_dir  = '/home/vito.galligani/Work/RTTOV/Plots/'
    savepath  = os.path.join(plot_dir,datetime.now().strftime('%d%m%Y'))
    
    
    # Open ncfile
    data = Dataset(ncfile,'r')   
        
    # Hydrometeor content
    A = {}
    A['swc']  = data.variables['snow'][:]      # kg/m3
    A['iwc']  = data.variables['ice'][:]       # kg/m3
    A['gwc']  = data.variables['graupel'][:]   # kg/m3
    A['rwc']  = data.variables['rain'][:]      # kg/m3
    A['h20']  = data.variables['q'][:]         # kg/kg  
    A['pf']   = data.variables['pressure_full'][:]   # Pa, n = 137
    A['ph']   = data.variables['pressure_half'][:]   # Pa, n = 138
    A['t']    = data.variables['t'][:]               # temperature
    A['gwp'] = data.variables['gwp'][:]       # orography
    A['swp'] = data.variables['swp'][:]       # orography
    A['rwp'] = data.variables['rwp'][:]       # orography
    A['iwp'] = data.variables['iwp'][:]       # orography


    A['longitude'] = data.variables['longitude'][:]       # orography
    A['lat'] = data.variables['lats'][:]       # orography
    
    # Plot IFS cloudy:
    index = 19304

    print(A['longitude'][index])
    print(A['lat'][index])
    
    plt.matplotlib.rc('font', family='serif', size = 8)
    fig = plt.figure(figsize=(6,4))
    ax1 = fig.add_subplot(1,2,1)        
    ax1.plot( A['rwc'][index,:], A['pf'][index,:]/100, linewidth=1.2, linestyle='-', color='darkblue',label='Rain')
    ax1.plot( A['swc'][index,:], A['pf'][index,:]/100, linewidth=1.2, linestyle='-', color='darkred',label='Snow')
    ax1.plot( A['gwc'][index,:], A['pf'][index,:]/100, linewidth=1.2, linestyle='-', color='darkgreen',label='Graupel')
    ax1.plot( A['iwc'][index,:], A['pf'][index,:]/100, linewidth=1.2, linestyle='-', color='darkviolet',label='Ice')
    ax1.set_xlabel(r'Hydrometeor Contents (kg/m${^3}$)') #, color='darkblue')
    ax1.set_ylabel('Pressure (hPa)')
    ax1.set_ylim([np.max(A['pf'][index,:])/100, 100])
    plt.legend()
    plt.grid(True)


    ax2 = fig.add_subplot(1,2,2)        
    line1, = ax2.plot( A['t'][index,:], A['pf'][index,:]/100, linewidth=1.2, linestyle='-', color='darkblue')
    ax2.set_xlabel('Temperature (K)', color='darkblue')
    ax2.set_ylim([np.max(A['pf'][index,:])/100, 600])
    ax2.set_xlim([270, 302])
    # Set color of spines and tick labels for the first axis
    ax2.spines['bottom'].set_color(line1.get_color())
    ax2.tick_params(axis='x', colors=line1.get_color())
    ax2.xaxis.label.set_color(line1.get_color())
    plt.grid(True)
    
    # Creating a second y-axis sharing the same y-axis
    ax3 = ax2.twiny()
    line2, = ax3.plot( A['h20'][index,:], A['pf'][index,:]/100, linewidth=1.2, linestyle='-', color='darkred')
    ax3.set_xlabel('q (kg/kg)', color='darkred')
    ax3.set_xscale('log')
    ax3.set_ylim([np.max(A['pf'][index,:])/100, 600])    
    # Set color of spines and tick labels for the first axis
    ax3.spines['top'].set_color(line2.get_color())    
    ax3.tick_params(axis='x', colors=line2.get_color())
    ax3.xaxis.label.set_color(line2.get_color())
    ax3.set_xlim([0.0009, 0.02])
    plt.grid(True)

    plt.suptitle('Test SCATT Profile', fontsize=10) #(-58.9,-30.3)
    plt.tight_layout()
    plt.savefig(savepath+'/'+'ExampleProf_scatteringandATM.png', bbox_inches='tight', dpi=300 )

    
    # I will be working with iprof1 and iprof2  
    if 1:
        iprof1, iprof2 = P4A.run_IFS_plots(savepath, ncfile)
    
    
    if 0:
        Tb1, Trans1 = main('atms', iprof1)    # MidLats
        Tb2, Trans2 = main('atms', iprof2)   # Tropics

        # iprof1:: MidLat and iprof2:: Tropics    
        Tb1, Trans1 = main('atms', iprof1)    # MidLats
        Tb2, Trans2 = main('atms', iprof2)   # Tropics

    return

#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
def main_scatt(sat, Aiprofile):
    
    plt.matplotlib.rc('font', family='serif', size = 12)

    # half pressure levels in IFS == nlevels 
    nlev = 138  
    # T, q, etc in rttov is now nlevels-1  == nlayers 

    # Satellite zenith angle    
    zenith =  np.arange(0,55,5)
    
    # Sensor-related 
    if sat == "ici":
        sat   = 'ici'
        nchan = 13
        print("sat == ici")
    elif sat == "ssmis":
        sat   = 'ssmis'
        nchan = 3
        print("sat == ssmis")
        
    elif sat == "atms":
        sat   = 'atms'
        nchan = 22
        opt = {'nrows':2,'ncols':3,'width':11,'height':8}
        # [CH1: 23.800, CH2:31.400, CH9: 55.500, CH16: 88.2, CH17: 165.5, CH18:183.31±7.0, 
        # CH20: 183.31±3.0, CH22: 183.31±1.0]
        #nchan_plots = [1, 2, 9, 16, 17, 18, 20, 22]
        nchan_plots = [2, 16, 17, 18, 20, 22]
        chantitles = ['23.8', '31.4', '50.3', '51.76', '52.8', 
                      '53.596', '54.4', '54.94','55.5', '57.29', 
                       '57.29$\pm$0.5329', '57.29$\pm$0.3702', '57.29$\pm$0.3442', 
                       '57.29$\pm$0.3322', '57.29$\pm$0.3267', '88.2', '165.5',  
                       '183.31$\pm$7', '183.31$\pm$4.5', '183.31$\pm$3', 
                       '183.31$\pm$1.8', '183.31$\pm$1']
        print("sat == atms")
        
    elif sat == "mhs":
        sat   = 'mhs'
        nchan = 5
        opt = {'nrows':2,'ncols':3,'width':8,'height':8}
        # [CH1: 23.800, CH2:31.400, CH9: 55.500, CH16: 88.2, CH17: 165.5, CH18:183.31±7.0, 
        # CH20: 183.31±3.0, CH22: 183.31±1.0]
        nchan_plots = [1, 2, 3, 4, 5]
        chantitles = ['89', '157', '183.31$\pm$3', '183.31$\pm$1', '190']
        print("sat == mhs")  
        
    # Profiles related profile file names    
    ncfolder  = '/home/vito.galligani/datosmunin3/RTTOV_betaTest/data/'
    ncfile    = 'AIRS_processed_W1.nc'
    data      = Dataset(ncfolder+ncfile,'r')   

    Ap_full   = data.variables['pressure_full'][:]   # Pa, n = 137
    p_full    = Ap_full[Aiprofile]
    t_full    = data.variables['t'][:][Aiprofile]

    plevels = {}
    plevels['rttov14_halflevels'] = data.variables['pressure_half'][:][Aiprofile].copy()
    plevels['rttov13_fulllevels'] = p_full.copy()
        
    # Re-allocate pressure levels between half pressure levels == plevels['rttov13_fulllevels']  
    pressure = np.empty(len(plevels['rttov14_halflevels'])-1)
    for ilayer in range(len(plevels['rttov14_halflevels'])-1):
        pressure[ilayer] = (plevels['rttov14_halflevels'][ilayer] + plevels['rttov14_halflevels'][ilayer+1])/2
    plevels['betwee_halflevels'] = pressure 

    # Re-allocate pressure levels between full pressure levels == plevels['between_fulllevels_136']  
    pressure = np.empty(len(plevels['rttov13_fulllevels'])-1)
    for ilayer in range(len(plevels['rttov13_fulllevels'])-1):
        pressure[ilayer] = (plevels['rttov13_fulllevels'][ilayer] + plevels['rttov13_fulllevels'][ilayer+1])/2
    plevels['between_fulllevels_136'] = pressure 
    
    # Output Plot dir (Within this folder, define the name of a sub-folder according to date)
    plot_dir   = '/home/vito.galligani/Work/RTTOV/Plots/'
    path_plots = os.path.join(plot_dir,datetime.now().strftime('%d%m%Y'))

    # overview plot config
    opt_lev = { 'nrows':opt['nrows'],'ncols':opt['ncols'],'width':opt['width'],
                'height':opt['height'],'path':path_plots+'/'+sat+'/'+str(Aiprofile), 'nlev':nlev, 'nchan':nchan, 
                'nchan_plots':nchan_plots, 'chan_plots_titles':chantitles}
    

    Tb, dTb, SCATT, HYDROSCATT = T2R.read_outouts_SCATT(ncfile, Aiprofile, nlev, nchan, sat)
    if 0:
        P4A.AllSkySensor(dTb, Tb, opt_lev)



    opt_lev = { 'nrows':opt['nrows'],'ncols':opt['ncols'],'width':opt['width'],
                'height':opt['height'],'path':path_plots+'/'+sat+'/'+str(Aiprofile), 'nlev':nlev, 'nchan':nchan, 
                'nchan_plots':nchan_plots, 'chan_plots_titles':chantitles, 'xlabel':r'SSA ($\omega$)', 'title':'Single scattering albedo (gas+hydro)', 
                'name':'compSSA'}
    P4A.BulkCompSensor(SCATT['rttov13_ssa'], SCATT['rttov14_ssa'],  plevels['betwee_halflevels']/100, 
                       1000, 800, opt_lev, ['rttov13','rttov14'])

    opt_lev = { 'nrows':opt['nrows'],'ncols':opt['ncols'],'width':opt['width'],
                'height':opt['height'],'path':path_plots+'/'+sat+'/'+str(Aiprofile), 'nlev':nlev, 'nchan':nchan, 
                'nchan_plots':nchan_plots, 'chan_plots_titles':chantitles, 'xlabel':'Asym. parameter (g)', 'title':'asymmetry parameter (gas+hydro)', 
                'name':'compASM'}
    P4A.BulkCompSensor(SCATT['rttov13_asm'], SCATT['rttov14_asm'],  plevels['betwee_halflevels']/100, 
                       1000, 600, opt_lev, ['rttov13','rttov14'])

    opt_lev = { 'nrows':opt['nrows'],'ncols':opt['ncols'],'width':opt['width'],
                'height':opt['height'],'path':path_plots+'/'+sat+'/'+str(Aiprofile), 'nlev':nlev, 'nchan':nchan, 
                'nchan_plots':nchan_plots, 'chan_plots_titles':chantitles, 'xlabel':r'Ext. coeff. (k$_{ext}$)', 'title':'Extinction coefficient (gas+hydro)', 
                'name':'compEXT'}
    P4A.BulkCompSensor(SCATT['rttov13_ext'], SCATT['rttov14_ext'],  plevels['betwee_halflevels']/100, 
                       1000, 900, opt_lev, ['rttov13','rttov14'])


    
    if 0: 
        opt_lev = { 'nrows':opt['nrows'],'ncols':opt['ncols'],'width':opt['width'],
                    'height':opt['height'],'path':path_plots+'/'+sat+'/'+str(Aiprofile), 'nlev':nlev, 'nchan':nchan, 
                    'nchan_plots':nchan_plots, 'chan_plots_titles':chantitles, 'xlabel':r'Ext. coeff. (k$_{ext}$)', 'title':'Hydro-only Extinction coefficient', 
                    'name':'kext'}    
        P4A.HYDRO_CompSensor(HYDROSCATT['rttov14_ext'], HYDROSCATT['rttov13_ext'], plevels['betwee_halflevels']/100,  
                           1000, 10, opt_lev)
                             
        opt_lev = { 'nrows':opt['nrows'],'ncols':opt['ncols'],'width':opt['width'],
                    'height':opt['height'],'path':path_plots+'/'+sat+'/'+str(Aiprofile), 'nlev':nlev, 'nchan':nchan, 
                    'nchan_plots':nchan_plots, 'chan_plots_titles':chantitles, 'xlabel':'Asym. parameter (g)', 'title':'Hydro-only asymmetry parameter', 
                    'name':'asym'}    
        P4A.HYDRO_CompSensor(HYDROSCATT['rttov14_asm'], HYDROSCATT['rttov13_asm'], plevels['betwee_halflevels']/100,  
                           1000, 10, opt_lev)        
    
        opt_lev = { 'nrows':opt['nrows'],'ncols':opt['ncols'],'width':opt['width'],
                    'height':opt['height'],'path':path_plots+'/'+sat+'/'+str(Aiprofile), 'nlev':nlev, 'nchan':nchan, 
                    'nchan_plots':nchan_plots, 'chan_plots_titles':chantitles, 'xlabel':r'SSA ($\omega$)', 'title':'Hydro-only Single scattering albedo ', 
                    'name':'ssa'}    
        P4A.HYDRO_CompSensor(HYDROSCATT['rttov14_ssa'], HYDROSCATT['rttov13_ssa'], plevels['betwee_halflevels']/100,  
                           1000, 10, opt_lev)        




    return Tb, dTb, SCATT, HYDROSCATT



#------------------------------------------------------------------------------
plot_dir  = '/home/vito.galligani/Work/RTTOV/Plots/'
savepath  = os.path.join(plot_dir,datetime.now().strftime('%d%m%Y'))
ncfolder  = '/home/vito.galligani/datosmunin3/RTTOV_betaTest/data/'
ncfile    = ncfolder+'AIRS_processed_W1.nc'
        


iprof1, iprof2 = P4A.run_IFS_plots(savepath, ncfile)
Tb2, Trans2 = main('atms', iprof2)   # Tropics    
#Tb1, Trans1 = main('atms', iprof1)   # midlat    


# data.variables['pressure_half'][iprof1][136]/100
# data.variables['pressure_full'][iprof1][136]/100
# data.variables['orography'][iprof1]/100

file = 'Level2SpaceTransmitt_CS_minus1_0.0.txt'
path1 = '/home/vito.galligani/Work/RTTOV/rttov14.0_beta/rttov_test/beta_test/atms/TestClearSky_Nprof6577/Prop/'
path2 = '/home/vito.galligani/Work/RTTOV/rttov14.0_beta/rttov_test/beta_test/atms/TestClearSky_Nprof6577_option1/Prop/' 
path3 = '/home/vito.galligani/Work/RTTOV/rttov14.0_beta/rttov_test/beta_test/atms/TestClearSky_Nprof6577_option2b/Prop/' 
test1 = T2R.SingleCol(path1+file)
test2 = T2R.SingleCol(path2+file)
test3 = T2R.SingleCol(path3+file)


# plots_4SingleProfiles()

# Aiprofile = 19304
# Tb, dTb, SCATT, HYDROSCATT = main_scatt('atms', Aiprofile)    # MidLats

      
# y tambien quiero very hydro fractions









