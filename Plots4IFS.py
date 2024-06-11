"""
-----------------------------------------------------------------
@purpose : Analyse rttov 13.4 outputs  
@author  : Vito Galligani
@email   : vito.galligani@gmail.com
@condaenv: conda activate RTTOV14
-----------------------------------------------------------------
@main    : Plotting functions for IFS simulations
          
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

import sys
import numpy as np
from netCDF4 import Dataset
import os 
import typhon  as ty
import matplotlib.pyplot as plt
import cartopy
import matplotlib.pyplot as plt
import matplotlib
import cartopy.feature as cfeature
import cartopy.crs as ccrs
import matplotlib as mpl
import matplotlib.ticker as mticker
from cartopy.mpl.ticker import (LongitudeFormatter, LatitudeFormatter,
                            LatitudeLocator)
    




#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
def TBTB_csplots(rttov13, rttov14, rttov14_1, rttov14_2b, options):
    """
    -------------------------------------------------------------
    scatter plots of TBrttov13 vs TBrtto14 (for both clear and all sky)
    For rttov14 add other coeff runs. 
    -------------------------------------------------------------
    OUT    
    IN     
    -------------------------------------------------------------
    """
    plt.matplotlib.rc('font', family='serif', size = 12)

    unity = np.arange(230,300,1)    

    freq_title = []
    for i in options['nchan_plots']:
        freq_title.append(r'CH'+str(i))

    fig = plt.figure(figsize=(options['width'],options['height']))
    for i, item in enumerate(options['nchan_plots']):
        ax1 = fig.add_subplot(options['nrows'],options['ncols'],i+1)
        ax1.plot(rttov14[0,:,item-1], rttov13[0,:,item-1], linestyle='none',color='k', marker='o',markersize = 2, mfc='none')
        ax1.plot(unity, unity, linestyle='-',color='red', linewidth=1.2)        
        plt.title( freq_title[i]+'('+str(options['chan_plots_titles'][item-1])+' GHz)', fontsize='12', fontweight='bold')
        ax1.set_xlabel(r'RTTOV14 $T_\mathrm{B}$ [K]', color='k')
        ax1.set_ylabel(r'RTTOV13 $T_\mathrm{B}$ [K]', color='k')
        ax1.grid('true')
        # agrego mean y std 
        diff = rttov14[0,:,item-1] - rttov13[0,:,item-1] 
        mbe = 'mean: '+str( round(np.mean(diff),2)) + '\n'
        std = 'std: ' + str( round(np.std(diff), 2))
        plt.text(240, 280, mbe+std)
        
        
    fig.suptitle('Nadir ClearSky ATMS (rtcoeff13)' ,fontweight='bold' )
    plt.tight_layout()
    plt.subplots_adjust(bottom=0.2, wspace=0.33) # top=0.899)
    plt.savefig(options['path']+'/'+'IFS_cs_TBscatternadir_rtcoef13.png', bbox_inches='tight')

    fig = plt.figure(figsize=(options['width'],options['height']))
    for i, item in enumerate(options['nchan_plots']):
        ax1 = fig.add_subplot(options['nrows'],options['ncols'],i+1)
        ax1.plot(rttov14_1[0,:,item-1], rttov13[0,:,item-1], linestyle='none',color='k', marker='o',markersize = 2, mfc='none')
        ax1.plot(unity, unity, linestyle='-',color='red', linewidth=1.2)        
        plt.title( freq_title[i]+'('+str(options['chan_plots_titles'][item-1])+' GHz)', fontsize='12', fontweight='bold')
        ax1.set_xlabel(r'RTTOV14 $T_\mathrm{B}$ [K]', color='k')
        ax1.set_ylabel(r'RTTOV13 $T_\mathrm{B}$ [K]', color='k')
        ax1.grid('true')
        # agrego mean y std 
        diff = rttov14_1[0,:,item-1] - rttov13[0,:,item-1] 
        mbe = 'mean: '+str( round(np.mean(diff),2)) + '\n'
        std = 'std: ' + str( round(np.std(diff), 2))
        plt.text(240, 280, mbe+std)
        
    fig.suptitle('Nadir ClearSky ATMS (rttov14 option1 coeffs)' ,fontweight='bold' )
    plt.tight_layout()
    plt.subplots_adjust(bottom=0.2, wspace=0.33) # top=0.899)
    plt.savefig(options['path']+'/'+'IFS_cs_TBscatternadir_rtcoef14option1.png', bbox_inches='tight')

    fig = plt.figure(figsize=(options['width'],options['height']))
    for i, item in enumerate(options['nchan_plots']):
        ax1 = fig.add_subplot(options['nrows'],options['ncols'],i+1)
        ax1.plot(rttov14_2b[0,:,item-1], rttov13[0,:,item-1], linestyle='none',color='k', marker='o',markersize = 2, mfc='none')
        ax1.plot(unity, unity, linestyle='-',color='red', linewidth=1.2)        
        plt.title( freq_title[i]+'('+str(options['chan_plots_titles'][item-1])+' GHz)', fontsize='12', fontweight='bold')
        ax1.set_xlabel(r'RTTOV14 $T_\mathrm{B}$ [K]', color='k')
        ax1.set_ylabel(r'RTTOV13 $T_\mathrm{B}$ [K]', color='k')
        ax1.grid('true')
        # agrego mean y std 
        diff = rttov14_2b[0,:,item-1] - rttov13[0,:,item-1] 
        mbe = 'mean: '+str( round(np.mean(diff),2)) + '\n'
        std = 'std: ' + str( round(np.std(diff), 2))
        plt.text(240, 280, mbe+std)
        
        fig.suptitle('Nadir ClearSky ATMS (rttov14 option2b coeffs)' ,fontweight='bold' )
    
    plt.tight_layout()
    plt.subplots_adjust(bottom=0.2, wspace=0.33) # top=0.899)
    plt.savefig(options['path']+'/'+'IFS_cs_TBscatternadir_rtcoef14option2b.png', bbox_inches='tight')
    
    
    return

#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
def TBTB_asplots(rttov13, rttov14, options, title, axistitle, depression):
    """
    -------------------------------------------------------------
    scatter plots of TBrttov13 vs TBrtto14 (for both clear and all sky)
    For rttov14 add other coeff runs. 
    -------------------------------------------------------------
    OUT    
    IN     
    -------------------------------------------------------------
    """
    plt.matplotlib.rc('font', family='serif', size = 12)

    unity = np.arange(50,300,1)    

    freq_title = []
    for i in options['nchan_plots']:
        freq_title.append(r'CH'+str(i))

    fig = plt.figure(figsize=(options['width'],options['height']))
    for i, item in enumerate(options['nchan_plots']):
        ax1 = fig.add_subplot(options['nrows'],options['ncols'],i+1)
        pl1 = ax1.plot(rttov14[0,:,item-1], rttov13[0,:,item-1], linestyle='none',color='k', marker='o',markersize = 2, mfc='none')
        
        ax1.plot(unity, unity, linestyle='-',color='red', linewidth=1.2)

        plt.title( freq_title[i]+'('+str(options['chan_plots_titles'][item-1])+' GHz)', fontsize='12', fontweight='bold')

        
        ax1.set_xlabel(r'RTTOV14 '+ axistitle +' $T_\mathrm{B}$ [K]', color='k')
        ax1.set_ylabel(r'RTTOV13 '+ axistitle +' $T_\mathrm{B}$ [K]', color='k')
        ax1.grid('true')
        
        if depression == 1:
            ax1.set_xlim([-100, 0])
            ax1.set_ylim([-100, 0])
            # agrego mean y std 
            diff = rttov14[0,:,item-1] - rttov13[0,:,item-1] 
            mbe = 'mean: '+str( round(np.mean(diff[ rttov14[0,:,item-1]<-1]) ,2)) + '\n'
            std = 'std: ' + str( round(np.std(diff[ rttov14[0,:,item-1]<-1]), 2))
            plt.text(-90, -40, mbe+std)
            
        else:
            ax1.set_xlim([100, 300])
            ax1.set_ylim([100, 300])            
            # agrego mean y std pero para cloud impacts:
            diff = rttov14[0,:,item-1] - rttov13[0,:,item-1] 
            mbe = 'mean: '+str( round(np.mean(diff),2)) + '\n'
            std = 'std: ' + str( round(np.std(diff), 2))
            plt.text(120, 250, mbe+std)
            

    # External legend; line labels

    fig.suptitle('Nadir AllSky '+axistitle + '$T_\mathrm{B}$ ATMS' ,fontweight='bold' )
    plt.tight_layout()
    plt.savefig(options['path']+'/'+'IFS_as_'+title+'TBscatternadir.png', bbox_inches='tight',)

    
    return

#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
def stats_plots1(stats, options, title, axistitle):
    """
    -------------------------------------------------------------
    scatter plots of TBrttov13 vs TBrtto14 (for both clear and all sky)
    For rttov14 add other coeff runs. 
    -------------------------------------------------------------
    OUT    
    IN     
    -------------------------------------------------------------
    """
    plt.matplotlib.rc('font', family='serif', size = 12)


    freq_title = []
    for i in options['nchan_plots']:
        freq_title.append(r'CH'+str(i))

    fig = plt.figure(figsize=(options['width'],options['height']))
    for i, item in enumerate(options['nchan_plots']):
        ax1 = fig.add_subplot(options['nrows'],options['ncols'],i+1)
        pl1 = ax1.plot(np.arange(0,55,5),stats['MBE'][:,item-1], linestyle='none',color='darkred', marker='o',markersize = 2, mfc='none')
        ax1.fill_between(np.arange(0,55,5), stats['MBE'][:,item-1] - stats['std'][:,item-1], stats['MBE'][:,item-1] + stats['std'][:,item-1], 
                         color='darkred', alpha=0.2, label='1 Std Dev')
        
        plt.title( freq_title[i]+'('+str(options['chan_plots_titles'][item-1])+' GHz)', fontsize='12', fontweight='bold')

        ax1.set_xlabel(r'Satellite zenith angle', color='k')
        ax1.set_ylabel(r'Average deviation in $\Delta$ $T_\mathrm{B}$ [K]', color='k')
        ax1.grid('true')


    fig.suptitle('ATMS RTTOV14-RTTOV13' ,fontweight='bold' )
    plt.tight_layout()
    plt.savefig(options['path']+'/'+'IFS_ATMS_MEANTB_zenithangle', bbox_inches='tight',)
    
    return


#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
def stats_plots2(stats, options, title, axistitle):
    """
    -------------------------------------------------------------
    -------------------------------------------------------------
    OUT    
    IN     
    -------------------------------------------------------------
    """
    plt.matplotlib.rc('font', family='serif', size = 12)


    freq_title = []
    for i in options['nchan_plots']:
        freq_title.append(r'CH'+str(i))

    fig = plt.figure(figsize=(8.3,5))
    ax1 = fig.add_subplot(1,1,1)

    ax1.plot(np.arange(1,23,1), stats['MBE'][0,:], linestyle='none',color='darkred', marker='o',markersize = 2, mfc='none')
    ax1.fill_between(np.arange(1,23,1), stats['MBE'][0,:]-stats['std'][0,:], stats['MBE'][0,:] + stats['std'][0,:], color='darkred', alpha=0.2, label='1 Std Dev')

    ax1.plot(np.arange(1,23,1), stats['MBE'][-1,:], linestyle='none',color='darkblue', marker='o',markersize = 2, mfc='none')
    ax1.fill_between(np.arange(1,23,1), stats['MBE'][-1,:]-stats['std'][-1,:], stats['MBE'][-1,:] + stats['std'][-1,:], color='darkblue', alpha=0.2, label='1 Std Dev')

    ticks = np.arange(1, 23, 1)  # x-ticks from 0 to 10 with step of 1
    tick_labels = [f'CH{i}' for i in ticks]  # Custom tick labels
    plt.xticks(ticks, tick_labels, rotation=45)

    plt.title(r'ATMS $\Delta$ TB$_{RTTOV14}$ − $\Delta$ TB$_{RTTOV13}$', fontsize='12', fontweight='bold')

    ax1.set_xlabel(r'ATMS channel no.', color='k')
    ax1.set_ylabel(r'Average deviation in $\Delta$ $T_\mathrm{B}$ [K]', color='k')
    ax1.grid('true')
    #plt.tight_layout()
    plt.savefig(options['path']+'/'+'IFS_ATMS_MEANdTB_zenithangle', bbox_inches='tight',)
    
    return

#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
def stats_plots3(stats, options, title, axistitle):
    """
    -------------------------------------------------------------
    -------------------------------------------------------------
    OUT    
    IN      FOR CLEAR SKY! ONLY NO CLOUD THRESHOLD
    -------------------------------------------------------------
    """
    plt.matplotlib.rc('font', family='serif', size = 12)


    freq_title = []
    for i in options['nchan_plots']:
        freq_title.append(r'CH'+str(i))

    fig = plt.figure(figsize=(16,5))
    ax1 = fig.add_subplot(1,2,1)
    pl1=ax1.plot(np.arange(1,23,1), stats['cs_MBE0'][0,:], linestyle='none',color='darkred', marker='o',markersize = 2, mfc='none', label= 'rtcoef v13')
    pl2=ax1.plot(np.arange(1,23,1), stats['cs_MBE1'][0,:], linestyle='none',color='darkblue', marker='o',markersize = 2, mfc='none', label= 'rtcoef v14 option1')
    pl3=ax1.plot(np.arange(1,23,1), stats['cs_MBE2'][0,:], linestyle='none',color='darkgreen', marker='o',markersize = 2, mfc='none', label= 'rtcoef v14 option2b')

    # Collect handles and labels for the first legend
    handles1, labels1 = ax1.get_legend_handles_labels()
    legend1 = ax1.legend(handles1, labels1, loc="upper left")

    ax1.fill_between(np.arange(1,23,1), stats['cs_MBE0'][0,:]-stats['cs_std0'][0,:], stats['cs_MBE0'][0,:] + stats['cs_std0'][0,:], color='darkred', alpha=0.2, label='1 Std Dev')
    ax1.fill_between(np.arange(1,23,1), stats['cs_MBE1'][0,:]-stats['cs_std1'][0,:], stats['cs_MBE1'][0,:] + stats['cs_std1'][0,:], color='darkblue', alpha=0.2, label='1 Std Dev')
    ax1.fill_between(np.arange(1,23,1), stats['cs_MBE2'][0,:]-stats['cs_std2'][0,:], stats['cs_MBE2'][0,:] + stats['cs_std2'][0,:], color='darkgreen', alpha=0.2, label='1 Std Dev')

    ticks = np.arange(1, 23, 1)  # x-ticks from 0 to 10 with step of 1
    tick_labels = [f'CH{i}' for i in ticks]  # Custom tick labels
    plt.xticks(ticks, tick_labels, rotation=45)
    plt.title(r'nadir', fontsize='12', fontweight='bold')
    ax1.set_xlabel(r'ATMS channel no.', color='k')
    ax1.set_ylabel(r'Average deviation in $\Delta$ $T_\mathrm{B}$ [K]', color='k')
    ax1.grid('true')
    ax1.set_ylim([-0.4, 0.4])

    ax1 = fig.add_subplot(1,2,2)
    ax1.plot(np.arange(1,23,1), stats['cs_MBE0'][-1,:], linestyle='none',color='darkred', marker='o',markersize = 2, mfc='none', label= 'rtcoef v13')
    ax1.plot(np.arange(1,23,1), stats['cs_MBE1'][-1,:], linestyle='none',color='darkblue', marker='o',markersize = 2, mfc='none', label= 'rtcoef v14 option1')
    ax1.plot(np.arange(1,23,1), stats['cs_MBE2'][-1,:], linestyle='none',color='darkgreen', marker='o',markersize = 2, mfc='none', label= 'rtcoef v14 option2b')
    
    ax1.fill_between(np.arange(1,23,1), stats['cs_MBE0'][-1,:]-stats['cs_std0'][-1,:], stats['cs_MBE0'][-1,:] + stats['cs_std0'][-1,:], color='darkred', alpha=0.2, label='1 Std Dev')
    ax1.fill_between(np.arange(1,23,1), stats['cs_MBE1'][-1,:]-stats['cs_std1'][-1,:], stats['cs_MBE1'][-1,:] + stats['cs_std1'][-1,:], color='darkblue', alpha=0.2, label='1 Std Dev')
    ax1.fill_between(np.arange(1,23,1), stats['cs_MBE2'][-1,:]-stats['cs_std2'][-1,:], stats['cs_MBE2'][-1,:] + stats['cs_std2'][-1,:], color='darkgreen', alpha=0.2, label='1 Std Dev')

    ticks = np.arange(1, 23, 1)  # x-ticks from 0 to 10 with step of 1
    tick_labels = [f'CH{i}' for i in ticks]  # Custom tick labels
    plt.xticks(ticks, tick_labels, rotation=45)
    plt.title(r'50$^{o}$', fontsize='12', fontweight='bold')
    ax1.set_ylim([-0.4, 0.4])

    ax1.set_xlabel(r'ATMS channel no.', color='k')
    ax1.set_ylabel(r'Average deviation in $\Delta$ $T_\mathrm{B}$ [K]', color='k')
    ax1.grid('true')
    #plt.tight_layout()
    plt.suptitle(r'ATMS ClearSky TB (RTTOV13 -RTTOV14)', fontsize='12', fontweight='bold')
    plt.savefig(options['path']+'/'+'IFS_ATMS_MEAN_CLEARTBs_zenithangle_NADIRAND50', bbox_inches='tight',)

    
    return

#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
def stats_plots4_clear(stats, options, title, axistitle):
    """
    -------------------------------------------------------------
    scatter plots of TBrttov13 vs TBrtto14 (for both clear and all sky)
    For rttov14 add other coeff runs. 
    -------------------------------------------------------------
    OUT    
    IN     
    -------------------------------------------------------------
    """
    plt.matplotlib.rc('font', family='serif', size = 12)


    freq_title = []
    for i in options['nchan_plots']:
        freq_title.append(r'CH'+str(i))

    fig = plt.figure(figsize=(options['width'],options['height']))
    for i, item in enumerate(options['nchan_plots']):
        ax1 = fig.add_subplot(options['nrows'],options['ncols'],i+1)
        pl1 = ax1.plot(np.arange(0,55,5),stats['cs_MBE0'][:,item-1], linestyle='none',color='darkred', marker='o',markersize = 2, mfc='none', label= 'rtcoef v13')
        pl3 = ax1.plot(np.arange(0,55,5),stats['cs_MBE1'][:,item-1], linestyle='none',color='darkblue', marker='o',markersize = 2, mfc='none', label= 'rtcoef v14 option1')
        pl2 = ax1.plot(np.arange(0,55,5),stats['cs_MBE2'][:,item-1], linestyle='none',color='darkgreen', marker='o',markersize = 2, mfc='none', label= 'rtcoef v14 option2b')

        # Collect handles and labels for the first legend
        handles1, labels1 = ax1.get_legend_handles_labels()
        legend1 = ax1.legend(handles1, labels1, loc="upper left")
    
        ax1.fill_between(np.arange(0,55,5), stats['cs_MBE0'][:,item-1] - stats['cs_std0'][:,item-1], stats['cs_MBE0'][:,item-1] + stats['cs_std0'][:,item-1], 
                         color='darkred', alpha=0.2)
        ax1.fill_between(np.arange(0,55,5), stats['cs_MBE1'][:,item-1] - stats['cs_std1'][:,item-1], stats['cs_MBE1'][:,item-1] + stats['cs_std1'][:,item-1], 
                         color='darkblue', alpha=0.2)        
        ax1.fill_between(np.arange(0,55,5), stats['cs_MBE2'][:,item-1] - stats['cs_std2'][:,item-1], stats['cs_MBE2'][:,item-1] + stats['cs_std2'][:,item-1], 
                         color='darkblue', alpha=0.2) 
        
        plt.title( freq_title[i]+'('+str(options['chan_plots_titles'][item-1])+' GHz)', fontsize='12', fontweight='bold')

        ax1.set_xlabel(r'Satellite zenith angle', color='k')
        ax1.set_ylabel(r'Average deviation in ClearSky $T_\mathrm{B}$ [K]', color='k')
        ax1.grid('true')


    fig.suptitle('ATMS RTTOV14-RTTOV13' ,fontweight='bold' )
    plt.tight_layout()
    plt.savefig(options['path']+'/'+'IFS_ATMS_MEANclearTB_zenithangle', bbox_inches='tight',)
    
    return

