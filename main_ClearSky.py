#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
-----------------------------------------------------------------
@purpose : Analyse rttov 13.4 outputs  
@author  : Vito Galligani
@email   : vito.galligani@gmail.com
@condaenv: conda activate RTTOV14
-----------------------------------------------------------------
@main    : Analysis of clear-sky RTTOV v13.4 with a few IFS profiles
           using run_example_fwd_VITO.sh and looks at tbs, tau_levels, tau_tot
           for runs using both ClearSky and rttov-scatt with no clouds 
          
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

#-----------------------------------------------------------------------------
def MultiRow(fdat):
    """
    -------------------------------------------------------------
    Loads multiple ascii files stored in a single row. Example:
    brightness temperature from the modified fortran scripts. Each
    file contains one row, one value per channel.
    -------------------------------------------------------------
    OUT   dset    Dataset
    IN    fdat    Filename
    -------------------------------------------------------------
    """
    dset = []
    for fi in fdat:
        dsi = np.genfromtxt(fi)
        dset.append(dsi)
    dset = np.array(dset)
    return dset

#-----------------------------------------------------------------------------
def SingleCol(fdat):
    """
    -------------------------------------------------------------
    Loads a single column ascii file. Note here that these files
    are stored as follows:

    Do ichan = 1, nchannels
      Do ilayer = 1, nlevels
      Enddo
    Enddo
    -------------------------------------------------------------
    OUT   dset    Dataset
    IN    fdat    Filename
    -------------------------------------------------------------
    """
    with open(fdat,"r") as f :
        lines   = f.readlines()
        dset    = pd.DataFrame(l.strip().split() for l in lines)
    dset        = np.array(dset)
    dset        = dset[:,0].astype(float)
    f.close()
    return dset

#-----------------------------------------------------------------------------
def Level2SpeceTau(nchan, data_path, options):
    """
    -------------------------------------------------------------
    Loads total tau per surface zenith angle and returns the full
    dataset
    -------------------------------------------------------------
    OUT   orig_lev   Level to space tau at RTTOV's original 53
                     levels
          user_lev   Corresponding tau at user levels
    IN    nchan      Number of channels
          data_path  Path of the related data
          path       Path for storing the plots (date-wise)
    HIN   [-]        Corresponding ascii files
    -------------------------------------------------------------
    """
    # Calculated by rttov_opdep.F90
    orig_lev = SingleCol(str(data_path)+'Prop/Level2SpaceOpticalDepth54_CS.txt')
    user_lev = SingleCol(str(data_path)+'Prop/Level2SpaceOpticalDepthUserLevels_CS.txt')


    orig_lev = np.reshape(orig_lev, (nchan,53))
    user_lev = np.reshape(user_lev, (nchan,int(options['nlev'])))
    
    P4A.Plot4Level2SpaceTauSensor(-orig_lev, r'RTTOV levels [-]','Level2SpaceOpticalDepth54', options)
    P4A.Plot4Level2SpaceTauSensor(-user_lev, r'User levels [-]','Level2SpaceOpticalDepthUserLevels', options)
    
    return -orig_lev,-user_lev

#-----------------------------------------------------------------------------
def MultiCol(nlev,nchan,fdat):
    """
    -------------------------------------------------------------
    For all earth incident angles, loads and appends multi-column
    output from RTTOV. Example: transmission. It has nchan + 1
    columns. First column contains the number of model layers.
    Final row contains the # of channels.
    -------------------------------------------------------------
    OUT   out    Full dataset
    IN    fdat   List of filenames per sza
          nlev   The number of levels to read in
          nchan  Number of channels
    -------------------------------------------------------------
    """
    out = np.empty((0,nlev,nchan))
    for fi in fdat:
        with open(fi,"r") as f :
            lines   = f.readlines()
            data    = pd.DataFrame(l.strip().split() for l in lines)
            data    = np.array(data)
            data    = data.astype(float)
            data    = data[0:nlev,1:]    # From 2nd column
            data    = np.expand_dims(data,0)
            out     = np.append(out,data,axis=0)
            f.close()
    out = np.array(out)
    return out

#-----------------------------------------------------------------------------
def CumulTau(nlev,nchan,ds1,ds2):
    """
    -------------------------------------------------------------
    Loads tau data from RTTOV
    -------------------------------------------------------------
    OUT   out    Full dataset of cumulative tau
    IN    ds1    Column tau cumul data (nchan,nlev-1)
          ds2    Surface cumulative tau (nchan)
          nlev   The number of levels to read in. 
          nchan  Number of channels
    -------------------------------------------------------------
    """
    out = np.empty((0,nchan,nlev))
    
    for i,j in zip(ds1,ds2):
        dt  = np.empty((nchan,nlev))
        d1 = SingleCol(i)
        d1 = np.reshape(d1,(nchan,nlev-1))
        d1 = - d1

        d2 = SingleCol(j)
        d2 = np.reshape(d2,(nchan,1))
        d2 = - d2
        
        dt[:,:-1] = d1
        dt[:,-1]  = d2[:,0] + d1[:,-1]

        dt = np.expand_dims(dt,0)
        out= np.append(out,dt,axis=0)

    return out

#-----------------------------------------------------------------------------
def ToTau(nchan,fdat):
    """
    -------------------------------------------------------------
    Loads total tau per surface zenith angle and returns the full
    dataset
    -------------------------------------------------------------
    OUT   out    Full dataset of total tau
    IN    fdat   Dataset containing total taul for all zeniths
          nchan  Number of channels
    -------------------------------------------------------------
    """
    out = np.empty((0,nchan))
    for i in fdat:
        dsi = SingleCol(i)
        dsi = np.expand_dims(dsi,0)
        out = np.append(out,dsi,axis=0)
    return -out

#-----------------------------------------------------------------------------
def SingleLayerTau_ltr(nlev,nchan,ds2):                         
    """
    -------------------------------------------------------------
    Loads total tau per surface zenith angle and returns the full
    dataset
    -------------------------------------------------------------
    OUT   ds_u   Layer tau at RTTOV levels
    IN    ds2    Datset for layer tau at RTTOV layers (files)
          nlev   The number of levels to read in.
          nchan  Number of channels
    -------------------------------------------------------------
    """
    out = np.empty((0,nchan,nlev))
    for i in ds2:
        d2 = SingleCol(i)
        d2 = np.reshape(d2,(nchan,nlev))
        d2 = np.expand_dims(d2,0)
        out = np.append(out,d2,axis=0)
    return out

#-----------------------------------------------------------------------------
def LoadColAndResh(fname,nlev,nchan):
    """
    -------------------------------------------------------------
    Loads data stored column-wise by RTTOV-SCATT and reshapes it
    to (nlev,nchan)
    -------------------------------------------------------------
    OUT   var    Variable in nlev,nchan
    IN    fname  Variable in nlev x nchan
          nlev   Number of levels to read in. Based on ARTS
          nchan  Number of channels
    -------------------------------------------------------------
    """
    fname = sorted(glob.glob(str(fname)))
    print('paok',fname)
    var   = SingleCol(fname[0])
    var   = np.reshape(var,(nlev,nchan))
    return var

#-----------------------------------------------------------------------------
def RetrieveTau(exti,dlev,nchan):
    """
    -------------------------------------------------------------
    Compute layer optical depth from extinction and layer depth
    -------------------------------------------------------------
    OUT   taui   Layer optical depth
          taum   Total optical depth
    IN    exti   Layer extinction
          dlev   Layer depth
          nchan  Number of channels
    -------------------------------------------------------------
    """
    taui = np.zeros_like(exti)
    for i in np.arange(nchan):
        taui[:,i] = exti[:,i] * dlev[:]
    taum = np.sum(taui,axis = 0)
    return taui,taum

#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
def main(sat, Aiprofile):
    
    plt.matplotlib.rc('font', family='serif', size = 12)

    
    # Full pressure levels in IFS
    nlev = 137  

    # Satellite zenith angle    
    zenith =  np.arange(0,60,5)

    
    # Sensor-related 
    if sat == "ici":
        sat   = 'ici'
        nchan = 13
        nfreq = 2
        #tdim = (nchan-9)*nfreq + (nchan-4)*nfreq*2        
        print("sat == ici")
    elif sat == "ssmis":
        sat   = 'ssmis'
        nchan = 3
        nfreq = 2
        #tdim = (nchan-2)*nfreq + (nchan-1)*nfreq*2
        print("sat == ssmis")
    elif sat == "atms":
        sat   = 'atms'
        nchan = 22
        nfreq = 2
        opt = {'nrows':3,'ncols':3,'width':8,'height':8}
        # [CH1: 23.800, CH2:31.400, CH9: 55.500, CH16: 88.2, CH17: 165.5, CH18:183.31±7.0, 
        # CH20: 183.31±3.0, CH22: 183.31±1.0]
        nchan_plots = [1, 2, 9, 16, 17, 18, 20, 22]
        chantitles = ['23.8', '31.4', '50.3', '51.76', '52.8', 
                      '53.596', '54.4', '54.94','55.5', '57.29', 
                      '57.29$\pm$0.5329', '57.29$\pm$0.3702', '57.29$\pm$0.3442', 
                      '57.29$\pm$0.3322', '57.29$\pm$0.3267', '88.2', '165.5',  
                      '183.31$\pm$7', '183.31$\pm$4.5', '183.31$\pm$3', 
                      '183.31$\pm$1.8', '183.31$\pm$1']
        #tdim = (nchan-2)*nfreq + (nchan-1)*nfreq*2;
        print("sat == atms")
        
    # Path to main dir 
    parent_dir = '/home/vito.galligani/Work/RTTOV/rttov13.2/'
    # Path of rttov clear sky (cs) output
    prttov_cs  = parent_dir+'rttov_test/beta_test/atms/TestClearSky/'
    # Path of rttov all sky (as) output
    prttov_as  = parent_dir+'rttov_test/beta_test/atms/TestAllSky/'
    
    # Profiles related profile file names    
    ncfolder  = '/home/vito.galligani/datosmunin3/RTTOV_betaTest/data/'
    ncfile    = 'AIRS_processed_W1.nc'
    flag_name = 'rttov132_Profile_' +str( ncfile[-5:-3] )+'_'+str(Aiprofile)
    
    data = Dataset(ncfolder+ncfile,'r')   
    Ap_full   = data.variables['pressure_full'][:]   # Pa, n = 137
    p_full    =  Ap_full[Aiprofile]
    t_full    =  data.variables['t'][:][Aiprofile]
    
    # Output Plos
    # Define the folder name to store the plots
    plot_dir   = 'Plots/'

    # Within this folder, define the name of a sub-folder according to date
    path_plots = os.path.join(parent_dir,plot_dir,datetime.now().strftime('%d%m%Y'))
    
    # Load tau at original levels (_t53) and at given user levels (_tu) and give an overview plot
    opt_lev = { 'nrows':opt['nrows'],'ncols':opt['ncols'],'width':opt['width'],
                'height':opt['height'],'path':path_plots, 'nlev':nlev, 'nchan':nchan, 
                'nchan_plots':nchan_plots, 'chan_plots_titles':chantitles}
    
    
    # RTTOV 13.4 Clear Sky data
    #-------------------------------------------------------------------------
    # 1/ tb (brightness temperature)
    frttov_tb  = sorted(glob.glob(str(prttov_cs)+'output_tb_atm_zenith*em1.0'+str(flag_name)+'.dat'))
    rttov_tb   = MultiRow(frttov_tb)

    frttov_tb      = sorted(glob.glob(str(prttov_as)+'output_atm_zenith*em1.0'+str(flag_name)+'.dat'))
    rttov_tb_scatt = MultiRow(frttov_tb)

    #-------------------------------------------------------------------------
    # 2/ tr (transmittance) (In old-code == output_trans )
    frttov_tr  = sorted(glob.glob(str(prttov_cs)+'output_taulevel_atm_zenith*em1.0'+str(flag_name)+'.dat'))
    rttov_tr   = MultiCol(nlev,nchan,frttov_tr)
    
    #-------------------------------------------------------------------------
    # 3/ Cumulative tau == optical depth
    # ct (cumulative tau, [nchan,136] ) 
    frttov_ct  = sorted(glob.glob(str(prttov_cs)+'Prop/Level2SpaceTauUserLevels_CS_minus1_*.txt'))[::-1]
    # cts (Last cumulative tau level, [nchan,1]
    frttov_cts = sorted(glob.glob(str(prttov_cs)+'Prop/Level2SpaceTauUserLevels_CS_R_*.txt'))[::-1]
    # Cumulative tau
    rttov_ct = CumulTau(nlev,nchan,frttov_ct,frttov_cts)  

    #-------------------------------------------------------------------------
    # tt (Total tau == total optical depth)
    frttov_tt  = sorted(glob.glob(str(prttov_cs)+'Prop/Surface2SpaceOpticalDepth_CS_*.txt'))[::-1]
    rttov_tt = ToTau(nchan,frttov_tt)         
    for iz in range(len(zenith)): 
        print('------- zenith: '+str(zenith[iz]))
        print('rttov total od: ', np.round(rttov_tt[iz,:],3))

    #-------------------------------------------------------------------------
    # 4/ Cumulative tau too == Level to space optical depth (sin tomar en cuenta surface altitude)
    # Load tau at original levels (_t53) and at given user levels (_tu) and give an overview plot
    # HERE I LOOK AT Level2SpaceOpticalDepth54. La diferencia con Level2SpaceTauUserLevels_CS_minus1_*.txt 
    # es que en los de rttov_trans se toma en cuenta la superficie. sino son iguales ... 
    rttov_t53, rttov_tu = Level2SpeceTau(nchan, prttov_cs, opt_lev)

    #-------------------------------------------------------------------------
    # 5/ Plot transmittance per earth incident angle
    # Overview transmission plot. One plot per surface zenithal angle
    # TowardsTrasnmOver == Overview plot
    alt = ty.physics.pressure2height(np.flipud(p_full.data), np.flipud(t_full.data))
    alt = np.flipud(alt)
    for iz in range(len(zenith)): 
        zenith_in = str(zenith[iz])
        suptitle = r'Transmission - Sat. zenith angle ('+zenith_in+r'$^{\circ}$)'
        filename = r'Transmission_SatZenithAngle_'+zenith_in
        P4A.OverviewSensor(alt,rttov_tr[iz,:,:],-0.05,1.05,0,14,suptitle,filename,opt_lev,r'Transmission [-]')    

    #-------------------------------------------------------------------------
    # 6/ Plot TBs
    P4A.ClearSkySensorOverMulti(rttov_tb,rttov_tb_scatt,'BT (Clearsky)', opt_lev)
    
    #------------------------ RTTOV-SCATT TEST CLEAR -------------------------
    #-------------------------------------------------------------------------

    #-------------------------------------------------------------------------
    # 7/ OD analysis (RTTOV 13.2 and RTTOV-SCATT 13.2 for clear-sky)
    #-------------------------------------------------------------------------
    # 7.1 Read RTTOV-SCATT OD (w/ RTTOV calculated od) after rttov_mieproc.F90. 
    # note that RTTOV_"ltr" is the RTTOV levels optical depth at full pressure levels
    # BUT with nlevels+1 levels, where the lowest level (nlevels+1) is the difference
    # between total tau and tau_level(nlevels). 
    # (rttov_iniscatt)
    rttov_ltr_beEdd = sorted(glob.glob(prttov_as+'atm_em1.0'+flag_name+'_Prop/'+'OpticalDepthLayerFromTOA_BeEdd_*.txt'))[::-1]                      
    rttov_beEdd = SingleLayerTau_ltr(nlev,nchan,rttov_ltr_beEdd) 
    
    rttov_ltr_AfterEdd = sorted(glob.glob(prttov_as+'atm_em1.0'+flag_name+'_Prop/'+'OpticalDepthLayerFromTOA_AfterEdd_*.txt'))[::-1]                      
    rttov_AfterEdd = SingleLayerTau_ltr(nlev,nchan,rttov_ltr_AfterEdd) 

    # - Read RTTOV-SCATT OD (below we read rttov"LTU" == single layer user levels)
    # note that rttov_ltu is the od Between half levels,i.e., also full levels (see rttov_iniscatt.F90). 
    # rttov_ltu (od in rttov_iniscatt.F90) is the od used in rttov-scatt!! 
    rttov_ltu  = SingleCol(prttov_as+'atm_em1.0'+flag_name+'_Prop/'+'SingleLayerOpticalDepthsOurLevels_cs_180.0.txt')
    rttov_ltu = np.reshape(rttov_ltu, (nchan, nlev))   # This is 're-allocation in 1/2 levels'! nlevels    
    
    # This are ltr == single layer od using transmission % tau_levels on full levels 
    frttov_ltr = SingleCol(prttov_as+'atm_em1.0'+flag_name+'_Prop/'+'SingleLayerOpticalDepthsRTTOVLevels_cs_180.0.txt')
    rttov_ltr  = np.reshape(frttov_ltr, (nchan, 138))  
       
    #--------------
    # - Read RTTOV-SCATT SSA (SCATT-AUX%SSA)
    # (mw_scatt/rttov_mieproc.F90)
    rttov_ssa  = np.reshape(SingleCol(prttov_as+'atm_em1.0'+flag_name+'_Prop/'+'ssa_full_180.0.txt'),(nchan,nlev))
    
    #- Read RTTOV-SCATT EXT (SCATT-AUX%EXT)
    # (mw_scatt/rttov_mieproc.F90)
    rttov_ext  = np.reshape(SingleCol(prttov_as+'atm_em1.0'+flag_name+'_Prop/'+'ext_full_180.0.txt'),(nchan,nlev))
          
    # Ext from mw/rttov_iniscatt before hydro-mieprox-iniedd.
    extProNadir = LoadColAndResh(prttov_as+'atm_em1.0'+flag_name+'_Prop/'+'ExtProHydro_180.0.txt',nlev,nchan)
    
    # Ext from mw/rttov_iniscatt before hydro-mieprox-iniedd.
    fextPro_angles = sorted(glob.glob(prttov_as+'atm_em1.0'+flag_name+'_Prop/'+'ExtProHydro_*'))[::-1]  
    extPro_angles = np.empty((0,nlev,nchan))
    for i in fextPro_angles:
        d = LoadColAndResh(i,nlev, nchan)
        d = np.expand_dims(d,0)
        extPro_angles = np.append(extPro_angles,d,axis=0)
    extPro_angles = np.array(extPro_angles)   
        
    # Ext from  mw_scatt/rttov_mieproc.F90 (KPP HYDRO ONLY)
    rttov_files = sorted(glob.glob(prttov_as+'atm_em1.0'+flag_name+'_Prop/'+'ext4angle_*.txt'))[::-1]                      
    scatt_ext = np.empty((0,nlev,nchan))
    for i in rttov_files:
        d_ext = LoadColAndResh(i, nlev, nchan)
        d_ext = np.expand_dims(d_ext,0)
        scatt_ext = np.append(scatt_ext,d_ext,axis=0)
    scatt_ext = np.array(scatt_ext)


    # Level depth (mw_scatt/rttov_iniscatt.F90) 
    rttov_files_dz  = sorted(glob.glob(prttov_as+'atm_em1.0'+flag_name+'_Prop/'+'HalfLevelDepth_*.txt'))
    scatt_dz = np.empty((0,nlev))
    for i in rttov_files_dz:
        d = SingleCol(i)
        d = np.expand_dims(d,0)
        scatt_dz = np.append(scatt_dz,d,axis=0)
    scatt_dz = np.array(scatt_dz)

    # Derive the layer and the total optical depth for each zenith angle:
    scatt_tau = np.empty((0,nlev,nchan)) 
    scatt_taum = np.empty((0,nchan)) 
    for i in range(scatt_ext.shape[0]):
        dscatt_tau, dscatt_taum = RetrieveTau(scatt_ext[i,:,:],scatt_dz[i,:],nchan)
        dscatt_tau  = np.expand_dims(dscatt_tau,0)
        dscatt_taum = np.expand_dims(dscatt_taum,0)
        scatt_tau   = np.append(scatt_tau,dscatt_tau,axis=0)
        scatt_taum  = np.append(scatt_taum,dscatt_taum,axis=0)
    scatt_tau  = np.array(scatt_tau)
    scatt_taum = np.array(scatt_taum)
    
    # Single scattering albedo. The ssa4angle_180.0.txt is multiplied by ext
    # (mw_scatt/rttov_mieproc.F90)
    rttov_files = sorted(glob.glob(prttov_as+'atm_em1.0'+flag_name+'_Prop/'+'ssa4angleT_*.txt'))[::-1]                      
    scatt_ssa = np.empty((0,nlev,nchan))
    for i in rttov_files:
        d_ssa = LoadColAndResh(i, nlev, nchan)
        d_ssa = np.expand_dims(d_ssa,0)
        scatt_ssa = np.append(scatt_ssa,d_ssa,axis=0)
    scatt_ssa = np.array(scatt_ssa)

    # Needed for the plot
    scatt_ssa[np.where(scatt_ssa == 0)] = np.nan

    # Asymmetry parameter
    rttov_files = sorted(glob.glob(prttov_as+'atm_em1.0'+flag_name+'_Prop/'+'ass4angleT_*.txt'))[::-1]                      
    scatt_ass = np.empty((0,nlev,nchan))
    for i in rttov_files:
        d_ass = LoadColAndResh(i, nlev, nchan)
        d_ass = np.expand_dims(d_ass,0)
        scatt_ass = np.append(scatt_ass,d_ass,axis=0)
    scatt_ass = np.array(scatt_ass)

    # In rttov_transmit check this is od clear sky rttov ? 
    # frttov  = sorted(glob.glob(prttov_cs+'Prop/'+'Level2SpaceTransmitt_CS_minus1_*.0.txt'))[::-1]
    # clear_od = np.empty((0,nchan,nlev-1))
    # for i in frttov:
    #     d = LoadColAndResh(i, nchan, nlev-1)
    #     d = np.expand_dims(d,0)
    #     clear_od = np.append(clear_od,d,axis=0)
    # clear_od = np.array(clear_od)

    # # Add the surface level:
    # rttov_clearODs = np.empty((len(zenith),nchan,nlev))
    # rttov_clearODs[:,:,:-1] = clear_od
    # for j in range(len(frttov_cts)):     
    #     d2 = SingleCol(frttov_cts[j])
    #     d2 = np.reshape(d2,(nchan,1))
    #     d2 = - d2
    #     rttov_clearODs[j,:,-1] = d2[:,0]

    # COMPARE SINGLE LAYER OD in rttov and rttov_scatt
    # rttov_ltr is the one used in the clear sky rttov BUT nlev+1
    # rttov_ltr[0]      == Top RTTOV level to space
    # rttov_ltr[:]      == log( transmission%tau_levels(ilayer-1,ichan) ) - log( transmission%tau_levels (ilayer,ichan) )
    # rttov_ltr[nlev+1] == Surface to bottom RTTOV (full pressure) level
    # rttov_clearODs from rttov_transmit ? 
    #-----------------------------------------------------------------------------
    # rttov_ltu is the od between half levels == rttov_scatt y donde entra HE dry-air g(z)==constant
    #-----------------------------------------------------------------------------
    # Esto es lo que  me intersa mostrar que se arregla en v 14
    fname   = 'atm_em1.0'+flag_name
    opt_tau = {'nchan':nchan,'nchan_plots':nchan_plots,'title':str(sat)+' level optical depth: '+str(flag_name[9:])+'(nadir)',
               'nrows':opt_lev['nrows'],'ncols':opt_lev['ncols'],'width':2*opt_lev['width'],'height':2*opt_lev['height'],
               'name':'BulkLevelTau'+str(fname)+'Nadir','path':path_plots,'chan_plots_titles':chantitles,
               'xlabel':r'$\tau$ [-]','x_rd':r'$\Delta \tau$ [%]'}   
    # rttov single layer od: ltr vs. ltu () at 55degrees and nadir
    print(rttov_ltr.shape)
    P4A.Plot4SingleLayerTauSensor(rttov_ltr, rttov_ltu, rttov_beEdd[0,:,:], rttov_AfterEdd[0,:,:], 
                                  0, [0.005, 0.003, 0.4, 0.01, 0.05, 0.15, 0.4, 0.6], 
                                  0, 50, r'$\tau$(0$^o$) clear in',
                                 'rttovcs_LayerTau_zenith0', opt_tau, ['rttov_cs', 'rttovscatt_in', 'rttovscatt_beEdd', 
                                                                        'rttovscatt_AfterEdd'])    

    P4A.Plot4SingleLayerTauSensor(rttov_ltr, rttov_ltu, rttov_beEdd[-1,:,:], rttov_AfterEdd[-1,:,:], 
                                  0, [0.005, 0.003, 0.4, 0.01, 0.05, 0.15, 0.4, 0.6], 
                                  0, 50, r'$\tau$(55$^o$) clear in',
                                 'rttovcs_LayerTau_zenith55', opt_tau, ['rttov_cs', 'rttovscatt_in', 'rttovscatt_beEdd', 
                                                                        'rttovscatt_AfterEdd'])    

    #-------------------------------------------------------------------------------------------------------
    # total od compare
    #plot total od per zenith angle for clear and cloudy:
    frttov_tt  = sorted(glob.glob(str(prttov_cs)+'Prop/Surface2SpaceOpticalDepth_CS_*.txt'))[::-1]
    rttov_tt = ToTau(nchan,frttov_tt)         
    frttov_tt  = sorted(glob.glob(str(prttov_as)+'Prop/TotalOpticalDepthsRTTOVLevels_cs_*.txt'))[::-1]
    rttovscatt_tt = ToTau(nchan,frttov_tt)            

    
    fig = plt.figure(figsize=(opt_lev['width'],opt_lev['height']))
    for iz in range(len(zenith)): 
        ax1 = fig.add_subplot(4,3,iz+1)
        
        total_      = np.round(rttov_tt[iz,:],3)
        totalscatt_ = -1*np.round(rttovscatt_tt[iz,:],3)
       
        if iz==0:
            ax1.plot(np.arange(total_.shape[0]), total_, linestyle='None', color='darkred', marker='o', markersize=2, label='rttov cs')
            ax1.plot(np.arange(total_.shape[0]), totalscatt_, linestyle='None', color='darkblue', marker='o', markersize=2, label='rttov-scatt (clear)')
        else:             
            ax1.plot(np.arange(total_.shape[0]), total_, linestyle='None', color='darkred', marker='o', markersize=2)
            ax1.plot(np.arange(total_.shape[0]), totalscatt_, linestyle='None', color='darkblue', marker='o', markersize=2)
            
        ax1.set_xlabel(r'CH', color='k')
        ax1.set_ylabel(r'total $\tau$', color='k')
        ax1.grid('true')
        plt.title( r'zenith='+str(zenith[iz])+'$^o$', fontsize='12', fontweight='bold')
        
    fig.legend(   loc     = "lower right")

    fig.suptitle(r'Total $\tau$ (rttovscatt vs rttov)', fontweight='bold' )
    plt.tight_layout()
    plt.subplots_adjust(top=0.899)
    plt.savefig(opt_lev['path']+'/'+str('od_nadir_interpolation_error_rttovscatt_check')+'.png')
    plt.close()    
    
    fig = plt.figure(figsize=(opt_lev['width'],opt_lev['height']))
    for iz in range(len(zenith)): 
        ax1 = fig.add_subplot(4,3,iz+1)
        
        total_      = np.round(rttov_tt[iz,:],3)
        totalscatt_ = -1*np.round(rttovscatt_tt[iz,:],3)
       
        ax1.plot(np.arange(total_.shape[0]), 100*(total_-totalscatt_)/total_, linestyle='None', color='k', marker='o', markersize=2)

        ax1.set_xlabel(r'CH', color='k')
        ax1.set_ylabel('rttov-rttovscatt \n /rttov [%]', color='k')
        ax1.grid('true')
        plt.title( r'zenith='+str(zenith[iz])+'$^o$', fontsize='12', fontweight='bold')
        
    fig.suptitle(r'Total $\tau$ (rttov-rttovscatt/rttov) [%]', fontweight='bold' )
    plt.tight_layout()
    plt.subplots_adjust(top=0.899)
    plt.savefig(opt_lev['path']+'/'+str('DIFF_REL_od_nadir_interpolation_error_rttovscatt_check')+'.png')
    plt.close()    
    


    


    ############################################################# 
    # -----   HARD CODED FOR NADIR ONLY FOR NOW         --------#   
    #      OJO QUE TODO ESTO ES PARA HYDRO SCAT PROP            #
    #############################################################     
    #fname   = 'atm_em1.0'+flag_name
    #opt_ext = {'nchan':nchan,'nchan_plots':nchan_plots,'title':str(sat)+' bulk extinction: '+str(flag_name[9:])+'(nadir)',
    #           'nrows':opt_lev['nrows'],'ncols':opt_lev['ncols'],'width':opt_lev['width'],'height':opt_lev['height'],
    #           'name':'BulkExt'+fname+'Nadir','path':path_plots,'chan_plots_titles':chantitles,
    #           'xlabel':r'$\beta_\mathrm{ext}$ [km$^{-1}$]','x_rd':r'$\Delta \beta_\mathrm{ext}$ [%]'}
    #P4A.BulkCompSensor(np.squeeze(scatt_ext[0,:,:]), alt, -0.005, 0.05, 6, 12, opt_ext)     # key = 'mext'  (hydro-only ok)
    # 
    # Bulk level tau
    #opt_tau = {'nchan':nchan,'nchan_plots':nchan_plots,'title':str(sat)+' bulk level optical depth: '+str(flag_name[9:])+'(nadir)',
    #           'nrows':opt_lev['nrows'],'ncols':opt_lev['ncols'],'width':opt_lev['width'],'height':opt_lev['height'],
    #           'name':'BulkLevelTau'+str(fname)+'Nadir','path':path_plots,'chan_plots_titles':chantitles,
    #           'xlabel':r'$\tau$ [-]','x_rd':r'$\Delta \tau$ [%]'}   
    #P4A.BulkCompSensor(np.squeeze(scatt_tau[0,:,:]), alt, 0, 0.005, -0.2, 12, opt_tau)   # key = 'mdtau'
    #
    # Bulk single scattering albedo
    #opt_w = {'nchan':nchan,'nchan_plots':nchan_plots,'title':str(sat)+' bulk SSA: '+str(flag_name[9:])+'(nadir)',
    #           'nrows':opt_lev['nrows'],'ncols':opt_lev['ncols'],'width':opt_lev['width'],'height':opt_lev['height'],
    #           'name':'BulkSSA'+fname+'Nadir','path':path_plots,'chan_plots_titles':chantitles,
    #           'xlabel':r'$\omega_\mathrm{0}$ [-]','x_rd':r'$\Delta \omega_\mathrm{0}$ [%]'}   
    #P4A.BulkCompSensor(scatt_ssa[0,:,:], alt, 0.6, 1, 6, 12, opt_w)      # key = 'mw0'
    #
    # Bulk assymetry parameter
    #opt_g = {'nchan':nchan,'nchan_plots':nchan_plots,'title':str(sat)+' bulk asym.: '+str(flag_name[9:])+'(nadir)',
    #           'nrows':opt_lev['nrows'],'ncols':opt_lev['ncols'],'width':opt_lev['width'],'height':opt_lev['height'],
    #           'name':'BulkSSA'+fname+'Nadir','path':path_plots,'chan_plots_titles':chantitles,
    #           'xlabel':r'$g$ [-]','x_rd':r'$\Delta g$ [%]'}   
    #scatt_asy[np.where(scatt_asy == 0.0)] = np.nan
    #P4A.BulkCompSensor(abulk,scatt_asy,'asy',-0.005,0.5,6,12,opt_g)


    return
       
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
main('atms', 0)
    
    
    
