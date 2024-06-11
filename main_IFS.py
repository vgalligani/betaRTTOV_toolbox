#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
-----------------------------------------------------------------
@purpose : Analyse rttov 13.4 outputs of REALISITC IFS 
@author  : Vito Galligani
@email   : vito.galligani@gmail.com
@condaenv: conda activate RTTOV14
-----------------------------------------------------------------
@main    :  Plots for IFS simulations 
          
           Remember to open remote session:
               ssh -L 6000:yakaira:22 vito.galligani@portal.cima.fcen.uba.ar
               nohup python -m spyder_kernels.console -matplotlib='inline' -f=./remotemachine.json &
"""
#################################################
import numpy          as np
import Tools2Read     as T2R
import xarray         as xr
import Plots4IFS      as P4I
import glob
import os
from pathlib import Path
from datetime import datetime

#################################################


#-----------------------------------------------------------------------------
def Metrics(ds): #,dt_t):
    """
    ---------------------------------------------------------------------------
    Derives the metrics describing the agreement of RTMs per channel and per
    frequency.
    ---------------------------------------------------------------------------
    ---------------------------------------------------------------------------
    """
    dset = {}
    dset['dt_clearsky']  = ds.rttov13_cs - ds.rttov14_cs
    dset['dt_allsky']    = ds.rttov13_as - ds.rttov14_as
    dset['dt_rttov13']   = ds.rttov13_as - ds.rttov13_cs
    dset['dt_rttov14']   = ds.rttov14_as - ds.rttov14_cs
    dset['dt_CLEAR0']    = ds.rttov13_cs - ds.rttov14_cs
    dset['dt_CLEAR1']    = ds.rttov13_cs - ds.rttov14_cs_opt1
    dset['dt_CLEAR2']    = ds.rttov13_cs -  ds.rttov14_cs_opt2b
    dset['dt_d']         = dset['dt_rttov14'] - dset['dt_rttov13']
    dset['dt_rd']  = 100*(dset['dt_rttov13']-dset['dt_rttov14'])/dset['dt_rttov14']  # 100*(rttov13-rttov14)/rttov14)
    #dset['rmse']   = InCloudRMSE(dset['dt_a'].data,dset['dt_r'].data,dt_t)
    dset['MBE'], dset['std'] = InCloudMetrics(dset['dt_d'].data, dset['dt_rttov14'])
    dset['cs_MBE0'], dset['cs_std0'] = ClearMetrics(dset['dt_CLEAR0'].data)
    dset['cs_MBE1'], dset['cs_std1'] = ClearMetrics(dset['dt_CLEAR1'].data)
    dset['cs_MBE2'], dset['cs_std2'] = ClearMetrics(dset['dt_CLEAR2'].data)

    return dset

#-----------------------------------------------------------------------------
def InCloudMetrics(v,t):
    """
    ---------------------------------------------------------------------------
    dtb is to keep only cloud sceces 
    ---------------------------------------------------------------------------
    ---------------------------------------------------------------------------
    """
    mbe = np.zeros((v.shape[0],v.shape[2]))
    std = np.zeros((v.shape[0],v.shape[2]))
    for j in np.arange(v.shape[0]):           # nzen
        diff = v[j,:,:]
        dtb  = t[j,:,:]
        for i in np.arange(v.shape[2]):         # nchan
            mbe[j,i] = np.mean(diff[dtb[:,i]<-1,i])
            std[j,i] = np.std( diff[dtb[:,i]<-1,i])

    return mbe,std

#-----------------------------------------------------------------------------
def ClearMetrics(v):
    """
    ---------------------------------------------------------------------------
    dtb is to keep only cloud sceces 
    ---------------------------------------------------------------------------
    ---------------------------------------------------------------------------
    """
    mbe = np.zeros((v.shape[0],v.shape[2]))
    std = np.zeros((v.shape[0],v.shape[2]))
    for j in np.arange(v.shape[0]):           # nzen
        diff = v[j,:,:]
        for i in np.arange(v.shape[2]):         # nchan
            mbe[j,i] = np.mean(diff[:,i])
            std[j,i] = np.std( diff[:,i])

    return mbe,std


# #################################################
# def InCloudRMSE(dta,dtr,dt_t):
#     rmse  = np.zeros((dta.shape[1],dta.shape[2]))
#     for j in np.arange(dta.shape[1]):           # nzen
#       a = dta[:,j,:]
#       r = dtr[:,j,:]
#       for i in np.arange(dta.shape[2]):         # nchan
#           rmse[j,i] = sqrt(mean_squared_error( a[a[:,i]<-1,i],r[a[:,i]<-1,i] ))
#     return rmse



#-----------------------------------------------------------------------------

if __name__ == '__main__':

    # Output Plot dir (Within this folder, define the name of a sub-folder according to date)
    plot_dir   = '/home/vito.galligani/Work/RTTOV/Plots/'
    path_plots = os.path.join(plot_dir,datetime.now().strftime('%d%m%Y'))
    # IFS rttov simulations file 
    full_IFS   = '/home/vito.galligani/Work/RTTOV/EXPncfiles/full_experiments.nc'
    ds         = xr.open_dataset(full_IFS)

    # And create a netcdf file
    # "rttov13_cs"     
    # "rttov13_as"     
    # "rttov14_as"      
    # "rttov14_cs"
    # "rttov14_cs_opt1"
    # "rttov14_cs_opt2b"

    #-----     
    sat = 'ATMS'
    nchan = 22
    opt = {'nrows':2,'ncols':3,'width':11,'height':8}
    nchan_plots = [2, 16, 17, 18, 20, 22]
    chantitles = ['23.8', '31.4', '50.3', '51.76', '52.8', 
                          '53.596', '54.4', '54.94','55.5', '57.29', 
                           '57.29$\pm$0.5329', '57.29$\pm$0.3702', '57.29$\pm$0.3442', 
                           '57.29$\pm$0.3322', '57.29$\pm$0.3267', '88.2', '165.5',  
                           '183.31$\pm$7', '183.31$\pm$4.5', '183.31$\pm$3', 
                           '183.31$\pm$1.8', '183.31$\pm$1']
    
    # overview plot config
    optplot = { 'nrows':opt['nrows'],'ncols':opt['ncols'],'width':opt['width'],
                    'height':opt['height'],'path':path_plots+'/IFS', 'nchan':nchan, 
                    'nchan_plots':nchan_plots, 'chan_plots_titles':chantitles}


    #-----    
    # CLOUD METRICS and other calcs
    stats = Metrics(ds)
    rttov14_dTB      = ds.rttov14_as.data[:]-ds.rttov14_cs.data[:]
    rttov13_dTB      = ds.rttov13_as.data[:]-ds.rttov13_cs.data[:]
        
    #-----        
    if 0:
        
        #-----    
        # scatter plot clear sky 
        P4I.TBTB_csplots(ds.rttov13_cs.data[:], ds.rttov14_cs.data[:], ds.rttov14_cs_opt1.data[:], ds.rttov14_cs_opt2b.data[:], optplot)
        
        #-----         
        # scatter plot all sky BTs 
        P4I.TBTB_asplots(ds.rttov13_as.data[:], ds.rttov14_as.data[:], optplot, 'AllSkyTB', '', 0)    
        
    #-----  
    # scatter plot dBTs        
    P4I.TBTB_asplots(rttov13_dTB, rttov14_dTB, optplot, 'dTB', '$\Delta$', 1)
                
    #-----     
    # Mean dTb with cloud impacts (no clear sky! )
    P4I.stats_plots1(stats, optplot, '', '')    # stats['MBE']
    P4I.stats_plots2(stats, optplot, '', '')
    
    #-----        
    if 0:    
        #-----     
        # Mean cleat Tbs for different coeefs
        P4I.stats_plots3(stats, optplot, '', '')         # Clear sky as a function of ATMs CHANNEL for nadir  (and subplot with 50degrees)
        P4I.stats_plots4_clear(stats, optplot, '', '')   # Clear sky as a function of zenith angle for selected channels





        
        


    





