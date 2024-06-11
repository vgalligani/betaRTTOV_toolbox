#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
-----------------------------------------------------------------
@purpose : Analyse rttov 13.4 outputs  
@author  : Vito Galligani
@email   : vito.galligani@gmail.com
@condaenv: conda activate RTTOV14
-----------------------------------------------------------------
@main    : Plotting functions for single profile simulation analysis
          
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
def Plot4Level2SpaceTauSensor(ltau,ylab,name,options):
    """
    -------------------------------------------------------------
    Level to space optical depth at original RTTOV levels and
    interpolated to the user levels
    -------------------------------------------------------------
    OUT    name.png  Plot stored at the given path
    IN     ltau      Level-to-space tau
           ylab      ylabel
           name      Name to store the figure
           path      Path to store the figure
    -------------------------------------------------------------
    """
    plt.matplotlib.rc('font', family='serif', size = 12)
    
    freq_title = []
    for i in options['nchan_plots']:
        freq_title.append('CH'+str(i))

    fig = plt.figure(figsize=(options['width'],options['height']))

    for i,item in enumerate(options['nchan_plots']):
        ax1 = fig.add_subplot(options['nrows'],options['ncols'],i+1)
        ax1.plot(np.flipud(ltau[item-1,:]),np.arange(ltau.shape[1]),linewidth=0.05,color='darkblue' , marker='o',markersize = 2, mfc='none')

        ax1.plot(np.flipud(ltau[item-1,:]),np.arange(ltau.shape[1]),linewidth=0.05,color='darkred' , marker='o',markersize = 2, mfc='none')


        #plt.title( freq_title[i], fontsize='12', fontweight='bold')
        plt.title( freq_title[i]+'('+str(options['chan_plots_titles'][item-1])+' GHz)', fontsize='12', fontweight='bold')

        ax1.set_ylabel(str(ylab), color='k')
        ax1.set_xlabel('Optical depth [-]', color='k')
        ax1.grid('true')
    fig.suptitle('Level-to-space od ('+ylab[:-5]+')' ,fontweight='bold' )
    plt.tight_layout()
    plt.subplots_adjust(bottom=0.2, wspace=0.33)
    plt.savefig(options['path']+'/'+str(name)+'.png')
    plt.close()

    

    return

#-----------------------------------------------------------------------------
#def OverviewSensor(alt,d1,d2,x1,x2,y1,y2,ititle,name,options):
def OverviewSensor(alt,d1,x1,x2,y1,y2,ititle,name,options,xlabel):
    """
    -------------------------------------------------------------
    Comparison between RTTOV and ARTS for clear sky
    -------------------------------------------------------------
    OUT    name.png   Plot stored at the given path
    IN     alt        Altitude
           d1         RTTOV clear sky data (v13.2)
           d2         RTTOV-SCATT clear sky data (v13.2)
           x1         Lower limit at x axis
           x2         Upper limit at x axis
           y1         Lower limit at y axis
           y2         Upper limit at y axis
           xlabel     xlabel
           ititle     Title of figure
           name       Name to store the figure
           path       Path to store the figure
    -------------------------------------------------------------
    """
    plt.matplotlib.rc('font', family='serif', size = 12)
    
    freq_title = []
    for i in options['nchan_plots']:
        freq_title.append('CH'+str(i))

    fig = plt.figure(figsize=(options['width'],options['height']))
    for i, item in enumerate(options['nchan_plots']):
        ax1 = fig.add_subplot(options['nrows'],options['ncols'],i+1)
        pl1 = ax1.plot(d1[:,item-1],alt/1e3,linewidth=0.05,color='darkblue', marker='o',markersize = 2, mfc='none')
        # pl2 = ax1.plot(np.flipud(d2[:,i]),alt,linewidth=0.05,color='darkred'  , marker='x',markersize = 2)
        # plt.xscale('log')
        
        #plt.title( freq_title[i], fontsize='12', fontweight='bold')
        plt.title( freq_title[i]+'('+str(options['chan_plots_titles'][item-1])+' GHz)', fontsize='12', fontweight='bold')
        
        
        ax1.axhline(y=5   ,ls='--',color='k')
        ax1.set_ylabel(r'Altitude [km]', color='k')
        ax1.set_xlabel(xlabel, color='k')
        ax1.grid('true')

        if type(x2) == float:
            plt.xlim(x1,x2)
        else:
            plt.xlim(x1,x2[i])
        plt.ylim(y1,y2)
    
    # External legend; line labels
    llabels = ['RTTOV']#,'RTTOV']
    fig.legend([pl1], #, pl2],                   # The line objects
               labels        = llabels,      # The labels for each line
               loc           = "lower right",# Position of legend
               borderaxespad = 0.5,          # Small spacing around legend box
               title         = "RTTOV version:")     # Title for the legend
    
    fig.suptitle(str(ititle) ,fontweight='bold' )
    plt.tight_layout()
    plt.subplots_adjust(bottom=0.2, wspace=0.33)
    plt.savefig(options['path']+'/'+str(name)+'.png')
    plt.close()

    return

#-----------------------------------------------------------------------------
def ClearSkySensorOverMulti(d1,d2,title,options):
    """
    -------------------------------------------------------------
    As in ClearSkySensorOver but for 5 differemt types and ARTS
    or RTTOV
    -------------------------------------------------------------
    OUT    name.png   Plot stored at the given path
    IN     d1i        RTTOV  i = 1,...,5
           y1         Lower limit at y axis
           y2         Upper limit at y axis
           ybins      yticks, e.g., np.arange(y_lo,y_up,increment)
           options    Figure dependent options
    -------------------------------------------------------------
    """
    plt.matplotlib.rc('font', family='serif', size = 12)

    #yl1 = [290,280,240,255,265]
    #yl2 = [300,290,255,265,280]
    
    freq_title = []
    for i in options['nchan_plots']:
        freq_title.append(r'CH'+str(i))

    fig = plt.figure(figsize=(options['width'],options['height']))
    for i, item in enumerate(options['nchan_plots']):
        ax1 = fig.add_subplot(options['nrows'],options['ncols'],i+1)
        pl1 = ax1.plot(np.arange(0,55,5),d1[:,item-1],linewidth=1.2,color='black' , marker='o',markersize = 0, mfc='none')
        pl1 = ax1.plot(np.arange(0,55,5),d2[:,item-1],linewidth=1.2,color='black' , marker='o',markersize = 0, mfc='none')
        #pl2 = ax1.plot(np.arange(0,50,5),d1[:,i],linewidth=1.2,color='darkblue' , marker='o',markersize = 0, mfc='none')
        #pl3 = ax1.plot(np.arange(0,50,5),d1[:,i],linewidth=1.2,color='darkgreen' , marker='o',markersize = 0, mfc='none')
        #pl4 = ax1.plot(np.arange(0,50,5),d1[:,i],linewidth=1.2,color='darkred' , marker='o',markersize = 0, mfc='none')
        #pl5 = ax1.plot(np.arange(0,50,5),d1[:,i],linewidth=1.2,color='grey' , marker='o',markersize = 0, mfc='none')
        
        #plt.title( freq_title[i]+'('+str(options['chan_plots_titles'][item-1])+' GHz)', fontsize='12', fontweight='bold')
        plt.title( freq_title[i]+'('+str(options['chan_plots_titles'][item-1])+' GHz)', fontsize='12', fontweight='bold')

        ax1.set_xticks(np.arange(0,60,10))
        ax1.set_ylabel(r'$T_\mathrm{B}$ [K]', color='k')
        ax1.set_xlabel(r'Sat. zenith angle [$^\circ$]', color='k')
        ax1.grid('true')

        #ax1.axhline(y=0.1 ,ls='-',color='k')
        #ax1.axhline(y=0   ,ls='--',color='k')
        #ax1.axhline(y=-0.1,ls='-',color='k')
        
        #print(yl1[i],yl2[i])
        #plt.ylim(yl1[i],yl2[i])
        #ax1.set_yticks(ybins)

    # External legend; line labels
    llabels = ['RTTOV 13.2', 'RTTOV-SCATT 13.2']

    fig.legend([pl1], #, pl2, pl3, pl4, pl5],                   # The line objects
                labels        = llabels,      # The labels for each line
                loc           = "lower right",# Position of legend
                borderaxespad = 0.5,          # Small spacing around legend box
                title         = "RTTOV version:")     # Title for the legend

    fig.suptitle(title ,fontweight='bold' )
    plt.tight_layout()
    plt.subplots_adjust(bottom=0.2, wspace=0.33)

    plt.savefig(options['path']+'/'+'BT_overview.png')
    plt.close()
    
    return

#-----------------------------------------------------------------------------
def ClearSkySensor(d,title,options):
    """
    -------------------------------------------------------------
    As in ClearSkySensorOver but for 5 differemt types and ARTS
    or RTTOV
    -------------------------------------------------------------
    OUT    name.png   Plot stored at the given path
    IN     d          RTTOV  i = 1,...,5
           y1         Lower limit at y axis
           y2         Upper limit at y axis
           ybins      yticks, e.g., np.arange(y_lo,y_up,increment)
           options    Figure dependent options
    -------------------------------------------------------------
    """
    plt.matplotlib.rc('font', family='serif', size = 12)

    #yl1 = [290,280,240,255,265]
    #yl2 = [300,290,255,265,280]
    
    freq_title = []
    for i in options['nchan_plots']:
        freq_title.append(r'CH'+str(i))

    fig = plt.figure(figsize=(options['width'],options['height']))
    for i, item in enumerate(options['nchan_plots']):
        ax1 = fig.add_subplot(options['nrows'],options['ncols'],i+1)
        pl1 = ax1.plot(np.arange(0,55,5),d['rttov14'][:-1,item-1],linewidth=1.2,color='darkred', marker='o',markersize = 2, mfc='none')
        pl2 = ax1.plot(np.arange(0,55,5),d['rttov14_dummyscattTRUE'][:-1,item-1],linewidth=1.2,color='red', marker='o',markersize = 2, mfc='none')
        pl3 = ax1.plot(np.arange(0,55,5),d['rttov14_dummyscattFALSE'][:-1,item-1],linewidth=1.2,color='magenta', marker='o',markersize = 2, mfc='none')
        pl4 = ax1.plot(np.arange(0,55,5),d['rttov13'][:-1,item-1],linewidth=1.2,color='darkblue', marker='o',markersize = 2, mfc='none')
        pl5 = ax1.plot(np.arange(0,55,5),d['rttov13_dummyscatt'][:-1,item-1],linewidth=1.2,color='blue', marker='o',markersize = 2, mfc='none')

        plt.title( freq_title[i]+'('+str(options['chan_plots_titles'][item-1])+' GHz)', fontsize='12', fontweight='bold')

        ax1.set_xticks(np.arange(0,60,10))
        ax1.set_ylabel(r'$T_\mathrm{B}$ [K]', color='k')
        ax1.set_xlabel(r'Sat. zenith angle [$^\circ$]', color='k')
        ax1.grid('true')

    # External legend; line labels
    llabels = ['RTTOV 14', 'RTTOV 14 (scatt.T.)',  'RTTOV 14 (scatt.F.)', 'RTTOV 13.2', 'RTTOV 13.2 (scatt dummy)']

    fig.suptitle(title ,fontweight='bold' )
    plt.tight_layout()
    plt.subplots_adjust(bottom=0.2, wspace=0.33) # top=0.899)

    fig.legend([pl1, pl2, pl3, pl4, pl5], #, pl2, pl3, pl4, pl5],                   # The line objects
                labels        = llabels,      # The labels for each line
                loc           = "lower center",# Position of legend
                ncol          = 3, 
                borderaxespad = 0.5,          # Small spacing around legend box
                title         = "RTTOV version:", 
                fancybox=False, shadow=False)     # Title for the legend

    plt.savefig(options['path']+'/'+'BT_overview.png', bbox_inches='tight',)
    #plt.close()
    
    
    fig = plt.figure(figsize=(options['width'],options['height']))
    for i, item in enumerate(options['nchan_plots']):
        ax1 = fig.add_subplot(options['nrows'],options['ncols'],i+1)
        pl1 = ax1.plot(np.arange(0,55,5),d['rttov13'][:-1,item-1]-d['rttov13_dummyscatt'][:-1,item-1],linewidth=1.2,color='blue', marker='o',markersize = 2, mfc='none')
        pl2 = ax1.plot(np.arange(0,55,5),d['rttov14'][:-1,item-1]-d['rttov14_dummyscattTRUE'][:-1,item-1],linewidth=1.2,color='red', marker='o',markersize = 2, mfc='none')
        pl3 = ax1.plot(np.arange(0,55,5),d['rttov14'][:-1,item-1]-d['rttov14_dummyscattFALSE'][:-1,item-1],linewidth=1.2,color='magenta', marker='o',markersize = 2, mfc='none')

        plt.title( freq_title[i]+'('+str(options['chan_plots_titles'][item-1])+' GHz)', fontsize='12', fontweight='bold')

        ax1.set_xticks(np.arange(0,60,10))
        ax1.set_ylabel(r'$\Delta$ $T_\mathrm{B}$ [K]', color='k')
        ax1.set_xlabel(r'Sat. zenith angle [$^\circ$]', color='k')
        ax1.grid('true')

        #ax1.axhline(y=0.1 ,ls='-', color='gray')
        #ax1.axhline(y=0   ,ls='--',color='gray')
        #ax1.axhline(y=-0.1,ls='-', color='gray')
        
    # External legend; line labels
    llabels = ['RTTOV13.2(clear-scatt)', 'RTTOV14(clear-scatt.T.)', 'RTTOV14(clear-scatt.F.)']

    fig.suptitle(r'$\Delta$ '+title ,fontweight='bold' )
    plt.tight_layout()
    plt.subplots_adjust(bottom=0.2, wspace=0.33) # top=0.899)

    fig.legend([pl1, pl2, pl3], #, pl2, pl3, pl4, pl5],                   # The line objects
                labels        = llabels,      # The labels for each line
                loc           = "lower center",# Position of legend
                borderaxespad = 0.5,          # Small spacing around legend box
                title         = "RTTOV version:", 
                ncol          = 3, 
                fancybox=False, shadow=False)     # Title for the legend
    plt.savefig(options['path']+'/'+'BT_DIFFoverview.png')
    #plt.close()
            
    fig = plt.figure(figsize=(options['width'],options['height']))
    for i, item in enumerate(options['nchan_plots']):
        ax1 = fig.add_subplot(options['nrows'],options['ncols'],i+1)
        pl1 = ax1.plot(np.arange(0,55,5), 100*(d['rttov13'][:-1,item-1]-d['rttov13_dummyscatt'][:-1,item-1])/d['rttov13'][:-1,item-1] ,linewidth=1.2,color='blue', marker='o',markersize = 2, mfc='none')
        pl2 = ax1.plot(np.arange(0,55,5), 100*(d['rttov14'][:-1,item-1]-d['rttov14_dummyscattTRUE'][:-1,item-1])/d['rttov14'][:-1,item-1] ,linewidth=1.2,color='red', marker='o',markersize = 2, mfc='none')
        pl2 = ax1.plot(np.arange(0,55,5), 100*(d['rttov14'][:-1,item-1]-d['rttov14_dummyscattFALSE'][:-1,item-1])/d['rttov14'][:-1,item-1] ,linewidth=1.2,color='magenta', marker='o',markersize = 2, mfc='none')
        
        #pl3 = ax1.plot(np.arange(0,55,5), 100*(d['rttov14'][:,item-1]-d['rttov13'][:,item-1])/d['rttov14'][:,item-1] , linewidth=1.2, color='k', marker='o',markersize = 2, mfc='none')

        plt.title( freq_title[i]+'('+str(options['chan_plots_titles'][item-1])+' GHz)', fontsize='12', fontweight='bold')

        ax1.set_xticks(np.arange(0,60,10))
        ax1.set_ylabel(r'$\Delta$ $T_\mathrm{B}$ [%]', color='k')
        ax1.set_xlabel(r'Sat. zenith angle [$^\circ$]', color='k')
        ax1.grid('true')

        #ax1.axhline(y=0.1 ,ls='-', color='gray')
        #ax1.axhline(y=0   ,ls='--',color='gray')
        #ax1.axhline(y=-0.1,ls='-', color='gray')
        
    # External legend; line labels
    llabels = ['RTTOV13.2(clear-scatt)', 'RTTOV14(clear-scatt.T.)', 'RTTOV14(clear-scatt.F.)'] #'RTTOV14-RTTOV13.2(clear)']

    fig.suptitle(r'$\Delta$ Rel. '+title ,fontweight='bold' )
    plt.tight_layout()
    plt.subplots_adjust(bottom=0.2, wspace=0.33) # top=0.899)
    fig.legend([pl1, pl2, pl3], #, pl2, pl3, pl4, pl5],                   # The line objects
                labels        = llabels,      # The labels for each line
                loc           = "lower center",# Position of legend
                borderaxespad = 0.5,          # Small spacing around legend box
                title         = "RTTOV version:",
                fancybox=False, shadow=False, ncol=3)     # Title for the legend

    plt.savefig(options['path']+'/'+'BT_RELDIFFoverview.png')
    #plt.close()
    
    return


#-----------------------------------------------------------------------------
def AllSkySensor(diff,d,options):
    """
    -------------------------------------------------------------
    As in ClearSkySensorOver but for 5 differemt types and ARTS
    or RTTOV
    -------------------------------------------------------------
    OUT    name.png   Plot stored at the given path
    IN     d          RTTOV  i = 1,...,5
           y1         Lower limit at y axis
           y2         Upper limit at y axis
           ybins      yticks, e.g., np.arange(y_lo,y_up,increment)
           options    Figure dependent options
    -------------------------------------------------------------
    """
    plt.matplotlib.rc('font', family='serif', size = 12)

    #yl1 = [290,280,240,255,265]
    #yl2 = [300,290,255,265,280]
    
    freq_title = []
    for i in options['nchan_plots']:
        freq_title.append(r'CH'+str(i))


    fig = plt.figure(figsize=(options['width'],options['height']))
    for i, item in enumerate(options['nchan_plots']):
        ax1 = fig.add_subplot(options['nrows'],options['ncols'],i+1)
        pl1 = ax1.plot(np.arange(0,55,5),d['rttov13_scatt'][:-1,item-1],linewidth=1.2,color='darkred', marker='o',markersize = 2, mfc='none')
        pl2 = ax1.plot(np.arange(0,55,5),d['rttov14_scatt'][:-1,item-1],linewidth=1.2,color='darkblue', marker='o',markersize = 2, mfc='none')

        #pl1 = ax1.plot(np.arange(0,55,5),d['rttov13'][:-1,item-1],linewidth=1.2,color='darkred', linestyle='--', marker='o',markersize = 2, mfc='none')
        #pl2 = ax1.plot(np.arange(0,55,5),d['rttov14'][:-1,item-1],linewidth=1.2,color='darkblue', linestyle='--', marker='o',markersize = 2, mfc='none')
        plt.title( freq_title[i]+'('+str(options['chan_plots_titles'][item-1])+' GHz)', fontsize='12', fontweight='bold')

        ax1.set_xticks(np.arange(0,60,10))
        ax1.set_ylabel(r'$T_\mathrm{B}$ [K]', color='k')
        ax1.set_xlabel(r'Sat. zenith angle [$^\circ$]', color='k')
        ax1.grid('true')

        
    # External legend; line labels
    llabels = ['RTTOV13', 'RTTOV14']

    fig.suptitle(r'AllSky $T_\mathrm{B}$ ATMS' ,fontweight='bold' )
    plt.tight_layout()
    plt.subplots_adjust(bottom=0.2, wspace=0.33) # top=0.899)

    fig.legend([pl1, pl2],              # The line objects
                labels        = llabels,      # The labels for each line
                loc           = "lower center",# Position of legend
                borderaxespad = 0.5,          # Small spacing around legend box
                title         = "RTTOV version:", 
                ncol          = 3, 
                fancybox=False, shadow=False)     # Title for the legend
    plt.savefig(options['path']+'/'+'AllSky_BT_overview.png')
    



    fig = plt.figure(figsize=(options['width'],options['height']))
    for i, item in enumerate(options['nchan_plots']):
        ax1 = fig.add_subplot(options['nrows'],options['ncols'],i+1)
        pl1 = ax1.plot(np.arange(0,55,5),d['rttov14_scatt'][:-1,item-1]-d['rttov13_scatt'][:-1,item-1],linewidth=1.2,color='k', marker='o',markersize = 2, mfc='none')
        pl2 = ax1.plot(np.arange(0,55,5),d['rttov14'][:-1,item-1]-d['rttov13'][:-1,item-1],linewidth=1.2,linestyle='--', color='k', marker='o',markersize = 2, mfc='none')

        plt.title( freq_title[i]+'('+str(options['chan_plots_titles'][item-1])+' GHz)', fontsize='12', fontweight='bold')

        ax1.set_xticks(np.arange(0,60,10))
        ax1.set_ylabel(r'$T_\mathrm{B}$ [K]', color='k')
        ax1.set_xlabel(r'Sat. zenith angle [$^\circ$]', color='k')
        ax1.grid('true')

    fig.suptitle(r'AllSky $T_\mathrm{B}$(RTTOV14-RTTOV13) ATMS' ,fontweight='bold' )
    plt.tight_layout()
    plt.subplots_adjust(bottom=0.2, wspace=0.33) # top=0.899)

    fig.legend([pl1, pl2], #, pl2, pl3, pl4, pl5],                   # The line objects
                labels        = ['AllSky', 'Clear'],      # The labels for each line
                loc           = "lower center",# Position of legend
                ncol          = 3, 
                borderaxespad = 0.5,          # Small spacing around legend box
                title         = "RTTOV version:", 
                fancybox=False, shadow=False)     # Title for the legend
    
    plt.savefig(options['path']+'/'+'AllSky_diffBT_overview.png')



    fig = plt.figure(figsize=(options['width'],options['height']))
    for i, item in enumerate(options['nchan_plots']):
        ax1 = fig.add_subplot(options['nrows'],options['ncols'],i+1)
        pl1 = ax1.plot(np.arange(0,55,5),diff['rttov13'][:-1,item-1],linewidth=1.2,color='darkred', marker='o',markersize = 2, mfc='none')
        pl2 = ax1.plot(np.arange(0,55,5),diff['rttov14'][:-1,item-1],linewidth=1.2,color='darkblue', marker='o',markersize = 2, mfc='none')

        plt.title( freq_title[i]+'('+str(options['chan_plots_titles'][item-1])+' GHz)', fontsize='12', fontweight='bold')

        ax1.set_xticks(np.arange(0,60,10))
        ax1.set_ylabel(r'$\Delta$ T$_\mathrm{B}$ [K]', color='k')
        ax1.set_xlabel(r'Sat. zenith angle [$^\circ$]', color='k')
        ax1.grid('true')

    # External legend; line labels
    llabels = ['RTTOV 13', 'RTTOV 14']

    fig.suptitle(r'T$_\mathrm{B}$(AllSky) - T$_\mathrm{B}$(clear) ATMS',fontweight='bold' )
    plt.tight_layout()
    plt.subplots_adjust(bottom=0.2, wspace=0.33) # top=0.899)

    fig.legend([pl1, pl2], #, pl2, pl3, pl4, pl5],                   # The line objects
                labels        = llabels,      # The labels for each line
                loc           = "lower center",# Position of legend
                ncol          = 3, 
                borderaxespad = 0.5,          # Small spacing around legend box
                title         = "RTTOV version:", 
                fancybox=False, shadow=False)     # Title for the legend

    plt.savefig(options['path']+'/'+'AllSky_dTB_overview.png', bbox_inches='tight',)
    #plt.close()
    
    fig = plt.figure(figsize=(options['width'],options['height']))
    for i, item in enumerate(options['nchan_plots']):
        ax1 = fig.add_subplot(options['nrows'],options['ncols'],i+1)
        pl1 = ax1.plot(np.arange(0,55,5),100*(diff['rttov14'][:-1,item-1]-diff['rttov13'][:-1,item-1])/diff['rttov14'][:-1,item-1],
                       linewidth=1.2,color='k', marker='o',markersize = 2, mfc='none')

        plt.title( freq_title[i]+'('+str(options['chan_plots_titles'][item-1])+' GHz)', fontsize='12', fontweight='bold')

        ax1.set_xticks(np.arange(0,60,10))
        ax1.set_ylabel(r'$\Delta$ T$_\mathrm{B}$(RTTOV14)-$\Delta$ T$_\mathrm{B}$(RTTOV13) [%]', color='k')
        ax1.set_xlabel(r'Sat. zenith angle [$^\circ$]', color='k')
        ax1.grid('true')

    fig.suptitle(r'$\Delta$ T$_\mathrm{B}$(RTTOV14-RTTOV13) ATMS [%]',fontweight='bold' )
    plt.tight_layout()
    plt.savefig(options['path']+'/'+'AllSky_dTB_relDIFF.png', bbox_inches='tight',)




    fig = plt.figure(figsize=(options['width'],options['height']))
    for i, item in enumerate(options['nchan_plots']):
        ax1 = fig.add_subplot(options['nrows'],options['ncols'],i+1)
        pl1 = ax1.plot(np.arange(0,55,5),(diff['rttov14'][:-1,item-1]-diff['rttov13'][:-1,item-1]),
                       linewidth=1.2,color='k', marker='o',markersize = 2, mfc='none')

        plt.title( freq_title[i]+'('+str(options['chan_plots_titles'][item-1])+' GHz)', fontsize='12', fontweight='bold')

        ax1.set_xticks(np.arange(0,60,10))
        ax1.set_ylabel(r'$\Delta$ T$_\mathrm{B}$(RTTOV14)-$\Delta$ T$_\mathrm{B}$(RTTOV13) [K]', color='k')
        ax1.set_xlabel(r'Sat. zenith angle [$^\circ$]', color='k')
        ax1.grid('true')

    fig.suptitle(r'$\Delta$ T$_\mathrm{B}$(RTTOV14-RTTOV13) ATMS',fontweight='bold' )
    plt.tight_layout()
    plt.savefig(options['path']+'/'+'AllSky_dTB_DIFF.png', bbox_inches='tight',)
    
    
    return



    
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
def ClearSkySensor_coeffsComp(d,title,options):
    """
    -------------------------------------------------------------
    As in ClearSkySensorOver but for 5 differemt types and ARTS
    or RTTOV
    -------------------------------------------------------------
    OUT    name.png   Plot stored at the given path
    IN     d          RTTOV  i = 1,...,5
           y1         Lower limit at y axis
           y2         Upper limit at y axis
           ybins      yticks, e.g., np.arange(y_lo,y_up,increment)
           options    Figure dependent options
    -------------------------------------------------------------
    """
    plt.matplotlib.rc('font', family='serif', size = 12)

    #yl1 = [290,280,240,255,265]
    #yl2 = [300,290,255,265,280]
    
    freq_title = []
    for i in options['nchan_plots']:
        freq_title.append(r'CH'+str(i))

    fig = plt.figure(figsize=(options['width'],options['height']))
    for i, item in enumerate(options['nchan_plots']):
        ax1 = fig.add_subplot(options['nrows'],options['ncols'],i+1)
        pl1 = ax1.plot(np.arange(0,55,5),d['rttov14'][:-1,item-1],linewidth=1.2,color='k', marker='o',markersize = 2, mfc='none')
        pl2 = ax1.plot(np.arange(0,55,5),d['rttov13'][:-1,item-1],linewidth=1.2,color='darkblue', marker='o',markersize = 2, mfc='none')
        pl3 = ax1.plot(np.arange(0,55,5),d['rttov14_option1'][:-1,item-1],linewidth=1.2,color='red', marker='o',markersize = 2, mfc='none')
        #pl4 = ax1.plot(np.arange(0,55,5),d['rttov14_option2b'][:-1,item-1],linewidth=1.2,color='darkred', marker='o',markersize = 2, mfc='none')

        plt.title( freq_title[i]+'('+str(options['chan_plots_titles'][item-1])+' GHz)', fontsize='12', fontweight='bold')

        ax1.set_xticks(np.arange(0,60,10))
        ax1.set_ylabel(r'$T_\mathrm{B}$ [K]', color='k')
        ax1.set_xlabel(r'Sat. zenith angle [$^\circ$]', color='k')
        ax1.grid('true')

    # External legend; line labels
    llabels = ['RTTOV 14', 'RTTOV 13.2', 'RTTOV 14 (option1)', 'RTTOV 14 (option2b)']

    fig.legend([pl1, pl2, pl3], #, pl2, pl3, pl4, pl5],                   # The line objects
                labels        = llabels,      # The labels for each line
                loc           = "lower right",# Position of legend
                borderaxespad = 0.5,          # Small spacing around legend box
                title         = "RTTOV version:")     # Title for the legend

    fig.suptitle(title ,fontweight='bold' )
    plt.tight_layout()
    plt.subplots_adjust(bottom=0.2, wspace=0.33)

    plt.savefig(options['path']+'/'+'BT_overview_coeffsComparison.png')
    #plt.close()
    
    
    fig = plt.figure(figsize=(options['width'],options['height']))
    for i, item in enumerate(options['nchan_plots']):
        ax1 = fig.add_subplot(options['nrows'],options['ncols'],i+1)
        pl1 = ax1.plot(np.arange(0,55,5),d['rttov14'][:-1,item-1]-d['rttov13'][:-1,item-1],linewidth=1.2,color='k', marker='o',markersize = 2, mfc='none')
        pl2 = ax1.plot(np.arange(0,55,5),d['rttov14'][:-1,item-1]-d['rttov14_option1'][:-1,item-1],linewidth=1.2,color='red', marker='o',markersize = 2, mfc='none')
        pl3 = ax1.plot(np.arange(0,55,5),d['rttov14'][:-1,item-1]-d['rttov14_option2b'][:-1,item-1],linewidth=1.2,color='darkred', marker='o',markersize = 2, mfc='none')
        #pl4 = ax1.plot(np.arange(0,55,5),d['rttov13'][:,item-1]-d['rttov14_option2b'][:,item-1],linewidth=1.2,linestyle='--',color='k', marker='o',markersize = 2, mfc='none')
        #pl4 = ax1.plot(np.arange(0,55,5),d['rttov14_option2b'][:-1,item-1]-d['rttov13'][:-1,item-1],linewidth=1.2,linestyle='--',color='k', marker='o',markersize = 2, mfc='none')


        plt.title( freq_title[i]+'('+str(options['chan_plots_titles'][item-1])+' GHz)', fontsize='12', fontweight='bold')

        ax1.set_xticks(np.arange(0,60,10))
        ax1.set_ylabel(r'$\Delta$ $T_\mathrm{B}$ [K]', color='k')
        ax1.set_xlabel(r'Sat. zenith angle [$^\circ$]', color='k')
        ax1.grid('true')

        
    # External legend; line labels
    llabels = ['RTTOV14-RTTOV13.2: coefsv13', 'RTTOV14(coefs13-option1)', 'RTTOV14(coefs13-option2b)']

    fig.suptitle(r'$\Delta$ '+title ,fontweight='bold' )
    plt.tight_layout()
    plt.subplots_adjust(bottom=0.2, wspace=0.33)
    
    
    fig.legend([pl1, pl2, pl3], #, pl2, pl3, pl4, pl5],                   # The line objects
                labels        = llabels,      # The labels for each line
                loc           = "lower right",# Position of legend
                borderaxespad = 0.5,          # Small spacing around legend box
                title         = "RTTOV version:", ncol=3)     # Title for the legend

    plt.savefig(options['path']+'/'+'BT_DIFFoverview_coefs.png')
    #plt.close()
            
    fig = plt.figure(figsize=(options['width'],options['height']))
    for i, item in enumerate(options['nchan_plots']):
        ax1 = fig.add_subplot(options['nrows'],options['ncols'],i+1)
        pl1 = ax1.plot(np.arange(0,55,5),100*(d['rttov14'][:-1,item-1]-d['rttov13'][:-1,item-1])/d['rttov14'][:-1,item-1],linewidth=1.2,color='k', marker='o',markersize = 2, mfc='none')
        pl2 = ax1.plot(np.arange(0,55,5),100*(d['rttov14'][:-1,item-1]-d['rttov14_option1'][:-1,item-1])/d['rttov14'][:-1,item-1],linewidth=1.2,color='red', marker='o',markersize = 2, mfc='none')
        pl3 = ax1.plot(np.arange(0,55,5),100*(d['rttov14'][:-1,item-1]-d['rttov14_option2b'][:-1,item-1])/d['rttov14'][:-1,item-1],linewidth=1.2,color='darkred', marker='o',markersize = 2, mfc='none')

        plt.title( freq_title[i]+'('+str(options['chan_plots_titles'][item-1])+' GHz)', fontsize='12', fontweight='bold')

        ax1.set_xticks(np.arange(0,60,10))
        ax1.set_ylabel(r'$\Delta$ $T_\mathrm{B}$ [%]', color='k')
        ax1.set_xlabel(r'Sat. zenith angle [$^\circ$]', color='k')
        ax1.grid('true')

        
    # External legend; line labels
    llabels = ['RTTOV14-RTTOV13.2: coefsv13', 'RTTOV14(coefs13-option1)', 'RTTOV14(coefs13-option2b)']

    fig.legend([pl1, pl2, pl3],           
                labels        = llabels,      # The labels for each line
                loc           = "lower right",# Position of legend
                borderaxespad = 0.5,          # Small spacing around legend box
                title         = "RTTOV version:")     # Title for the legend

    fig.suptitle(r'$\Delta$ Rel. '+title ,fontweight='bold' )
    plt.tight_layout()
    plt.subplots_adjust(bottom=0.2, wspace=0.33)

    plt.savefig(options['path']+'/'+'BT_RELDIFFoverview_coeffsRELdiff.png')
    #plt.close()
    
    return


    
#-----------------------------------------------------------------------------




def NadirOD_angle_RTTOVSCATT(extPro_angles, title, options):
    """
    -------------------------------------------------------------
    As in ClearSkySensorOver but for 5 differemt types and ARTS
    or RTTOV
    -------------------------------------------------------------
    OUT    name.png   Plot stored at the given path
    IN     d          RTTOV  i = 1,...,5
           y1         Lower limit at y axis
           y2         Upper limit at y axis
           ybins      yticks, e.g., np.arange(y_lo,y_up,increment)
           options    Figure dependent options
    -------------------------------------------------------------
    """
    plt.matplotlib.rc('font', family='serif', size = 12)

    zenith =  np.arange(0,60,5)

    # check if extPro_angles == nadir od for all angle   (interpolation errors!)
    freq_title = []
    for i in options['nchan_plots']:
        freq_title.append(r'CH'+str(i))

    fig = plt.figure(figsize=(options['width'],options['height']))
    for i, item in enumerate(options['nchan_plots']):
        ax1 = fig.add_subplot(options['nrows'],options['ncols'],i+1)
        for iz in range(extPro_angles.shape[0]-1):
            diffrel =  100*(np.flipud(extPro_angles[0,:,item-1]) -  np.flipud(extPro_angles[iz,:,item-1]))/np.flipud(extPro_angles[0,:,item-1])
            if iz == 0:
                if i==0:
                    ax1.plot(diffrel, np.arange(extPro_angles.shape[1]), '-k', label=r'$\theta$='+str(zenith[iz])+'$^o$')
                else:
                    ax1.plot(diffrel, np.arange(extPro_angles.shape[1]), '-k', markersize=2)
            else:
                if i==0:
                    ax1.plot(diffrel, np.arange(extPro_angles.shape[1]), linestyle='None', marker='o', markersize=2, label=r'$\theta$='+str(zenith[iz])+'$^o$')
                else:
                    ax1.plot(diffrel, np.arange(extPro_angles.shape[1]), linestyle='None', marker='o', markersize=2)
                        
        plt.title( freq_title[i]+'('+str(options['chan_plots_titles'][item-1])+' GHz)', fontsize='12', fontweight='bold')

        ax1.set_ylabel(r'model level', color='k')
        ax1.set_xlabel(r'$\tau_{0^o}$-$\tau_{\theta,0^o}$ [%]', color='k')
        ax1.grid('true')


    fig.suptitle('nadir ext in mw scatt ('+title+')', fontweight='bold' )
    plt.tight_layout()
    plt.subplots_adjust(bottom=0.2, wspace=0.33)
    
    fig.legend(loc  = "lower center", ncol=6, borderaxespad = 0.5)
    print('ok')
    
    plt.savefig(options['path']+'/'+str('od_nadir_interpolation_error_rttovscatt_checkreldiff')+title+'.png')
    
    
    fig = plt.figure(figsize=(options['width'], options['height']))
    for i, item in enumerate(options['nchan_plots']):
        ax1 = fig.add_subplot(options['nrows'],options['ncols'],i+1)
        for iz in range(extPro_angles.shape[0]):
            diffrel =  np.flipud(extPro_angles[0,:,item-1]) -  np.flipud(extPro_angles[iz,:,item-1])
            if iz == 0:
                if i==0:
                    ax1.plot(diffrel, np.arange(extPro_angles.shape[1]), '-k', label=r'$\theta$='+str(zenith[iz])+'$^o$')
                else:
                    ax1.plot(diffrel, np.arange(extPro_angles.shape[1]), '-k', markersize=2)
            else:
                if i==0:
                    ax1.plot(diffrel, np.arange(extPro_angles.shape[1]), linestyle='None', marker='o', markersize=2, label=r'$\theta$='+str(zenith[iz])+'$^o$')
                else:
                    ax1.plot(diffrel, np.arange(extPro_angles.shape[1]), linestyle='None', marker='o', markersize=2)
                              
        plt.title( freq_title[i]+'('+str(options['chan_plots_titles'][item-1])+' GHz)', fontsize='12', fontweight='bold')
        ax1.set_ylim([0, 80])
        ax1.set_ylabel(r'model level', color='k')
        ax1.set_xlabel(r'$\tau_{0^o}$-$\tau_{\theta,0^o}$ ', color='k')
        ax1.grid('true')

    # External legend; line labels
    fig.suptitle(r'$\tau$(nadir) for scatttering ('+title+')', fontweight='bold' )
    plt.tight_layout()
    plt.subplots_adjust(bottom=0.2, wspace=0.33)

    fig.legend(loc  = "lower center", ncol=6, borderaxespad = 0.5)
    print('ok')

    plt.savefig(options['path']+'/'+str('od_nadir_interpolation_error_rttovscatt_checkdiff')+title+'.png')
    
    return
    
    
#-----------------------------------------------------------------------------
def scatt_aux_ext_compare(rttov13, rttov14, plevels, options):
    """
    -------------------------------------------------------------
    As in ClearSkySensorOver but for 5 differemt types and ARTS
    or RTTOV
    -------------------------------------------------------------
    OUT    name.png   Plot stored at the given path
    IN     d          RTTOV  i = 1,...,5
           y1         Lower limit at y axis
           y2         Upper limit at y axis
           ybins      yticks, e.g., np.arange(y_lo,y_up,increment)
           options    Figure dependent options
    -------------------------------------------------------------
    """
    plt.matplotlib.rc('font', family='serif', size = 12)

    zenith =  np.arange(0,55,5)

    # check if extPro_angles == nadir od for all angle   (interpolation errors!)
    freq_title = []
    for i in options['nchan_plots']:
        freq_title.append(r'CH'+str(i))

    fig = plt.figure(figsize=(2*options['width'],2*options['height']))
    for i, item in enumerate(options['nchan_plots']):
        ax1 = fig.add_subplot(options['nrows'],options['ncols'],i+1)
        pl1=ax1.plot(rttov13[0,:,item-1], plevels['betwee_halflevels']/100, linestyle='-', color='darkblue', marker='o', markersize=2)
        pl2=ax1.plot(rttov14[0,:,item-1], plevels['betwee_halflevels']/100, linestyle='-', color='darkred', marker='o', markersize=2)
                
        plt.title( freq_title[i]+'('+str(options['chan_plots_titles'][item-1])+' GHz)', fontsize='12', fontweight='bold')

        ax1.set_ylabel(r'Pressure [hPa]', color='k')
        ax1.set_xlabel(r'scatt_aux%ext $\tau$ [km-1]', color='k')
        ax1.grid('true')

    # External legend; line labels
    llabels = ['RTTOV13.2', 'RTTOV14'] #'RTTOV14-RTTOV13.2(clear)']

    fig.legend([pl1, pl2], #, pl2, pl3, pl4, pl5],                   # The line objects
                labels        = llabels,      # The labels for each line
                loc           = "lower right",# Position of legend
                borderaxespad = 0.5,          # Small spacing around legend box
                title         = "RTTOV version:")     # Title for the legend

    fig.suptitle(r'nadir scatt_aux%ext', fontweight='bold' )
    plt.tight_layout()
    plt.subplots_adjust(bottom=0.2, wspace=0.33)
    plt.savefig(options['path']+'/'+str('nadir_rttov_scattaux_compare')+'.png')    
    
    fig = plt.figure(figsize=(2*options['width'],2*options['height']))
    for i, item in enumerate(options['nchan_plots']):
        ax1 = fig.add_subplot(options['nrows'],options['ncols'],i+1)
        pl1=ax1.plot(rttov13[-1,:,item-1], plevels['betwee_halflevels']/100, linestyle='-', color='darkblue', marker='o', markersize=2)
        pl2=ax1.plot(rttov14[-1,:,item-1], plevels['betwee_halflevels']/100, linestyle='-', color='darkred', marker='o', markersize=2)
                
        plt.title( freq_title[i]+'('+str(options['chan_plots_titles'][item-1])+' GHz)', fontsize='12', fontweight='bold')

        ax1.set_ylabel(r'Pressure [hPa]', color='k')
        ax1.set_xlabel(r'scatt_aux%ext $\tau$ [km-1]', color='k')
        ax1.grid('true')

    # External legend; line labels
    llabels = ['RTTOV13.2', 'RTTOV14'] #'RTTOV14-RTTOV13.2(clear)']

    fig.legend([pl1, pl2], #, pl2, pl3, pl4, pl5],                   # The line objects
                labels        = llabels,      # The labels for each line
                loc           = "lower right",# Position of legend
                borderaxespad = 0.5,          # Small spacing around legend box
                title         = "RTTOV version:")     # Title for the legend

    fig.suptitle(r'55$^o$ scatt_aux%ext', fontweight='bold' )
    plt.tight_layout()
    plt.subplots_adjust(bottom=0.2, wspace=0.33)
    plt.savefig(options['path']+'/'+str('55degree_rttov_scattaux_compare')+'.png')    
    
    return


#-----------------------------------------------------------------------------
def ODSensor(d, plevels, options):
    """
    -------------------------------------------------------------
    As in ClearSkySensorOver but for 5 differemt types and ARTS
    or RTTOV
    -------------------------------------------------------------
    OUT    name.png   Plot stored at the given path
    IN     d          RTTOV  i = 1,...,5
           y1         Lower limit at y axis
           y2         Upper limit at y axis
           ybins      yticks, e.g., np.arange(y_lo,y_up,increment)
           options    Figure dependent options
    -------------------------------------------------------------
    """
    plt.matplotlib.rc('font', family='serif', size = 12)

    
    freq_title = []
    for i in options['nchan_plots']:
        freq_title.append(r'CH'+str(i))

    fig = plt.figure(figsize=(options['width'],options['height']))
    for i, item in enumerate(options['nchan_plots']):
        ax1 = fig.add_subplot(options['nrows'],options['ncols'],i+1)
        pl1 = ax1.plot( d['rttov14_cs_tauUser'][item-1,:], plevels['rttov13_fulllevels']/100, linewidth=1.2,color='darkred', marker='o',markersize = 2, mfc='none')
        pl2 = ax1.plot( d['rttov13_cs_tauUser'][item-1,:], plevels['between_fulllevels_136']/100, linewidth=1.2,color='cyan', marker='x',markersize = 2, mfc='none')
        pl3 = ax1.plot( d['rttov13_cs_layerOD_nlev+1_fullp'][item-1,:], plevels['rttov14_halflevels']/100, linewidth=1.2,color='darkblue', marker='o',markersize = 2, mfc='none')
        pl4 = ax1.plot( d['rttov13_as_layerOD_nlev_halfp'][item-1,:], plevels['betwee_halflevels']/100, linewidth=1.2,color='blue', marker='o',markersize = 2, mfc='none')

        plt.title( freq_title[i]+'('+str(options['chan_plots_titles'][item-1])+' GHz)', fontsize='12', fontweight='bold')

        ax1.set_xlabel(r'optical depth $\tau$', color='k')
        ax1.set_ylabel(r'Pressure [hPa]', color='k')
        ax1.grid('true')

    fig.suptitle(r'Layer $\tau$ (nadir)', fontweight='bold' ) 
    plt.tight_layout()  
    plt.subplots_adjust(bottom=0.2, wspace=0.33)
    
    # External legend; line labels
    llabels = ['RTTOV 14 (transmit)', 'RTTOV 13.2(transmit)', 'RTTOV 13.2 (cs-mieproc)', 'RTTOV 13.2 (as-mieproc)'] #, 'RTTOV 13.2 (as after Edd)']
    #llabels = ['RTTOV 13.2 (cs)', 'RTTOV 13.2 (as)'] #, 'RTTOV 13.2 (as after Edd)']
    fig.legend([pl1, pl2, pl3, pl4], #, pl3], #, pl2, pl3, pl4, pl5],                   # The line objects
                labels        = llabels,      # The labels for each line
                loc           = "lower center",# Position of legend
                borderaxespad = 0.5,          # Small spacing around legend box
                title         = "RTTOV version:", 
                fancybox=False, shadow=False, ncol=3)     # Title for the legend

    plt.savefig(options['path']+'/'+'LayerTau_overview.png')
    #plt.close()
    
    print('rttov 14: '+ str(d['rttov14_cs_tauUser'][:,-1]))
    print('rttov 13: '+ str(d['rttov13_cs_tauUser'][:,-1]))    
    print('rttov-scatt 13: '+ str(d['rttov13_as_layerOD_nlev_halfp'][:,-1]))    
    

    fig = plt.figure(figsize=(options['width'],options['height']))
    for i, item in enumerate(options['nchan_plots']):
        ax1 = fig.add_subplot(options['nrows'],options['ncols'],i+1)
        pl1 = ax1.plot( d['rttov14_cs_tauUser'][item-1,:], plevels['rttov13_fulllevels']/100, linewidth=1.2, linestyle='none', color='darkred', marker='o',markersize = 3, mfc='none')
        pl2 = ax1.plot( d['rttov13_cs_tauUser'][item-1,:], plevels['between_fulllevels_136']/100, linewidth=1.2,linestyle='none',color='red', marker='x',markersize = 3, mfc='none')
        #pl3 = ax1.plot( d['rttov13_cs_layerOD_nlev+1_fullp'][item-1,:], plevels['rttov14_halflevels']/100, linewidth=1.2,linestyle='none',color='blue', marker='x',markersize = 3, mfc='none')
        pl4 = ax1.plot( d['rttov13_as_layerOD_nlev_halfp'][item-1,:], plevels['betwee_halflevels']/100, linewidth=1.2,linestyle='none',color='blue', marker='x',markersize = 3, mfc='none')

        plt.title( freq_title[i]+'('+str(options['chan_plots_titles'][item-1])+' GHz)', fontsize='12', fontweight='bold')

        ax1.set_xlabel(r'optical depth $\tau$', color='k')
        ax1.set_ylabel(r'Pressure [hPa]', color='k')
        ax1.grid('true')

        plt.ylim(1000,980)

    fig.suptitle(r'Layer $\tau$ (nadir)', fontweight='bold' ) 
    plt.tight_layout()  
    plt.subplots_adjust(bottom=0.2, wspace=0.33)


    # External legend; line labels (si no hay parentesis biene de rttov_transmit !)
    #llabels = ['RTTOV 14', 'RTTOV 13.2', 'RTTOV 13.2 (od_rttov: rttov_iniscatt)', 'RTTOV 13.2 (od: rttov_iniscatt)'] #, 'RTTOV 13.2 (as after Edd)']
    llabels = ['RTTOV 14', 'RTTOV 13.2',  'RTTOV 13.2 (RTTOV-SCATT)'] #, 'RTTOV 13.2 (as after Edd)']   
    fig.legend([pl1, pl2, pl3, pl4], #, pl3], #, pl2, pl3, pl4, pl5],                   # The line objects
                labels        = llabels,      # The labels for each line
                loc           = "lower center",# Position of legend
                borderaxespad = 0.5,          # Small spacing around legend box
                title         = "RTTOV version:",      # Title for the legend
                fancybox=False, shadow=False, ncol=3)     # Title for the legend

    plt.savefig(options['path']+'/'+'LayerTau_transmittOD.png')
    #plt.close()
    
    
    return
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
def ODSensor_coeffs(d, plevels, options):
    """
    -------------------------------------------------------------
    As in ClearSkySensorOver but for 5 differemt types and ARTS
    or RTTOV
    -------------------------------------------------------------
    OUT    name.png   Plot stored at the given path
    IN     d          RTTOV  i = 1,...,5
           y1         Lower limit at y axis
           y2         Upper limit at y axis
           ybins      yticks, e.g., np.arange(y_lo,y_up,increment)
           options    Figure dependent options
    -------------------------------------------------------------
    """
    plt.matplotlib.rc('font', family='serif', size = 12)

    #yl1 = [290,280,240,255,265]
    #yl2 = [300,290,255,265,280]
    
    #d['rttov14_cs_tauUser_option2b'][:,136] = np.nan
    #d['rttov14_cs_tauUser_option1'][:,136] = np.nan
    
    print('opt2b: '+ str(d['rttov14_cs_tauUser_option2b'][:,136]))
    print('opt1: '+ str(d['rttov14_cs_tauUser_option1'][:,136]))
    
    print('opt2b: '+ str(d['rttov14_cs_tauUser_option2b'][:,135]))
    print('opt1: '+ str(d['rttov14_cs_tauUser_option1'][:,135]))    
    
    freq_title = []
    for i in options['nchan_plots']:
        freq_title.append(r'CH'+str(i))

    fig = plt.figure(figsize=(options['width'],options['height']))
    for i, item in enumerate(options['nchan_plots']):
        ax1 = fig.add_subplot(options['nrows'],options['ncols'],i+1)
        pl1 = ax1.plot( d['rttov14_cs_tauUser'][item-1,:], plevels['rttov13_fulllevels']/100, linewidth=1.2, linestyle='none', color='k', marker='o',markersize = 5, mfc='none')
        pl2 = ax1.plot( d['rttov13_cs_tauUser'][item-1,:], plevels['between_fulllevels_136']/100, linewidth=1.2,linestyle='none',color='blue', marker='x',markersize = 3, mfc='none')
        pl3 = ax1.plot( d['rttov14_cs_tauUser_option1'][item-1,:], plevels['rttov13_fulllevels']/100, linewidth=1.2, linestyle='none', color='magenta', marker='x',markersize = 2, mfc='none')
        pl4 = ax1.plot( d['rttov14_cs_tauUser_option2b'][item-1,:], plevels['rttov13_fulllevels']/100, linewidth=1.2, linestyle='none', color='darkred', marker='o',markersize = 3, mfc='none')

        plt.title( freq_title[i]+'('+str(options['chan_plots_titles'][item-1])+' GHz)', fontsize='12', fontweight='bold')

        #ax1.set_xticks(np.arange(0,60,10))
        ax1.set_xlabel(r'optical depth $\tau$', color='k')
        ax1.set_ylabel(r'Pressure [hPa]', color='k')
        ax1.grid('true')


        plt.ylim(1050,800)

    fig.suptitle(r'Layer $\tau$ (nadir) Coefficients Test', fontweight='bold' ) 
    plt.tight_layout()  
    plt.subplots_adjust(bottom=0.2, wspace=0.33)
    
    # External legend; line labels
    llabels = ['RTTOV14', 'RTTOV13.2', 'RTTOV14 (option1)', 'RTTOV14 (option2b)'] #, 'RTTOV 13.2 (as after Edd)']
    #llabels = ['RTTOV 13.2 (cs)', 'RTTOV 13.2 (as)'] #, 'RTTOV 13.2 (as after Edd)']
    fig.legend([pl1, pl2, pl3], #, pl3], #, pl2, pl3, pl4, pl5],                   # The line objects
                labels        = llabels,      # The labels for each line
                loc           = "lower center",# Position of legend
                borderaxespad = 0.5,          # Small spacing around legend box
                title         = "RTTOV version:", ncol=4)     # Title for the legend


    #plt.savefig(options['path']+'/'+'LayerTau_transmittOD_coeffs_WITHoutLOWESTOPTION2LEVEL.png')
    plt.savefig(options['path']+'/'+'LayerTau_transmittOD_coeffs_WITHLOWESTOPTION2LEVEL.png')
    
    #plt.close()
    
    
    
    #plt.close()

    
    
    return
#-----------------------------------------------------------------------------
def scatt_aux(d, plevels, options):
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

    fig = plt.figure(figsize=(options['width'],options['height']))
    for i, item in enumerate(options['nchan_plots']):
        ax1 = fig.add_subplot(options['nrows'],options['ncols'],i+1)
        pl1 = ax1.plot(d['rttov14_aux_extAll'][0,item-1,:], plevels['betwee_halflevels']/100, linewidth=1.2,linestyle='-', color='darkred', marker='o',markersize = 0, mfc='none')
        pl2 = ax1.plot(d['rttov14_aux_extClear'][0,item-1,:], plevels['betwee_halflevels']/100, linewidth=1.2,linestyle='none',color='darkblue', marker='x',markersize = 4, mfc='none')

        plt.title( freq_title[i]+'('+str(options['chan_plots_titles'][item-1])+' GHz)', fontsize='12', fontweight='bold')

        ax1.set_xlabel(r'$\tau$ [km-1 and at nadir]', color='k')
        ax1.set_ylabel(r'Pressure [hPa]', color='k')
        ax1.grid('true')

    # External legend; line labels
    llabels = ['extClear', 'extAll'] #, 'RTTOV 13.2 (as)']

    fig.legend([pl1, pl2],                    # The line objects
                labels        = llabels,      # The labels for each line
                loc           = "lower right",# Position of legend
                borderaxespad = 0.5,          # Small spacing around legend box
                title         = "RTTOV version:")     # Title for the legend

    fig.suptitle('RTTOV-14 scatt_aux consistency (clear/scatt dummy)', fontweight='bold' )
    plt.tight_layout()
    plt.subplots_adjust(bottom=0.2, wspace=0.33)

    plt.savefig(options['path']+'/'+'rttov14_scatt_aux_check.png')
    #plt.close()
    

    return

#-----------------------------------------------------------------------------
def TransSensor_coeffs(d, title, options, rttov14_halflevels, rttov13_fulllevels):
    """
    -------------------------------------------------------------
    As in ClearSkySensorOver but for 5 differemt types and ARTS
    or RTTOV
    -------------------------------------------------------------
    OUT    name.png   Plot stored at the given path
    IN     d          RTTOV  i = 1,...,5
           y1         Lower limit at y axis
           y2         Upper limit at y axis
           ybins      yticks, e.g., np.arange(y_lo,y_up,increment)
           options    Figure dependent options
    -------------------------------------------------------------
    """
    plt.matplotlib.rc('font', family='serif', size = 12)
    
    freq_title = []
    for i in options['nchan_plots']:
        freq_title.append(r'CH'+str(i))

    fig = plt.figure(figsize=(options['width'],options['height']))
    for i, item in enumerate(options['nchan_plots']):
        ax1 = fig.add_subplot(options['nrows'],options['ncols'],i+1)
        pl1 = ax1.plot(np.arange(0,55,5),d['rttov14_cs_Surface2spaceTr'][:-1,item-1],linewidth=1.2,color='k', marker='o',markersize = 2, mfc='none')
        pl2 = ax1.plot(np.arange(0,55,5),d['rttov13_cs_Surface2spaceTr'][:-1,item-1],linewidth=1.2,color='darkblue', marker='o',markersize = 2, mfc='none')
        pl3 = ax1.plot(np.arange(0,55,5),d['rttov14_cs_option1_Surface2spaceTr'][:-1,item-1],linewidth=1.2,color='red', marker='o',markersize = 2, mfc='none')
        pl4 = ax1.plot(np.arange(0,55,5),d['rttov14_cs_option2b_Surface2spaceTr'][:-1,item-1],linewidth=1.2,color='darkred', marker='o',markersize = 2, mfc='none')
        plt.title( freq_title[i]+'('+str(options['chan_plots_titles'][item-1])+' GHz)', fontsize='12', fontweight='bold')

        ax1.set_xticks(np.arange(0,60,10))
        ax1.set_ylabel(r'Transmittance', color='k')
        ax1.set_xlabel(r'Sat. zenith angle [$^\circ$]', color='k')
        ax1.grid('true')

    # External legend; line labels
    llabels = ['RTTOV14', 'RTTOV13.2', 'RTTOV14 (option1)', 'RTTOV14 (option2b)']

    # fig.legend([pl1,pl2, pl3, pl4],             # The line objects
    #             labels        = llabels,      # The labels for each line
    #             loc           = "upper center",# Position of legend 'lower right'
    #             borderaxespad = 0.5,          # Small spacing around legend box
    #             title         = "RTTOV version:")     # Title for the legend


    fig.suptitle(title, fontweight='bold' )
    plt.tight_layout()
    plt.subplots_adjust(bottom=0.3, wspace=0.33) # top=0.899)

    fig.legend([pl1,pl2, pl3, pl4],             # The line objects
                labels        = llabels,      # The labels for each line
                loc           = "upper center",# Position of legend 'lower right'
                borderaxespad = 0.5,          # Small spacing around legend box
                title         = "RTTOV version:")     # Title for the legend

    plt.savefig(options['path']+'/'+'Surf2TOA_TotalTransmittance_overview_coeffs.png')
    #plt.close()
    
    
    # Check if  rttov14_as_Surface2spaceTr_FALSE == rttov14_as_Surface2spaceTr_TRUE
    fig = plt.figure(figsize=(options['width'],options['height']))
    for i, item in enumerate(options['nchan_plots']):
        ax1 = fig.add_subplot(options['nrows'],options['ncols'],i+1)
        pl1 = ax1.plot(np.arange(0,55,5),d['rttov14_cs_Surface2spaceTr'][:-1,item-1]-d['rttov14_cs_option1_Surface2spaceTr'][:-1,item-1],
                       linewidth=1.2,color='red', marker='x',markersize = 4, mfc='none')
        pl2 = ax1.plot(np.arange(0,55,5),d['rttov14_cs_Surface2spaceTr'][:-1,item-1]-d['rttov14_cs_option2b_Surface2spaceTr'][:-1,item-1],
                       linewidth=1.2,color='darkred', marker='o',markersize = 2, mfc='none')
        pl3 = ax1.plot(np.arange(0,55,5),d['rttov14_cs_option1_Surface2spaceTr'][:-1,item-1]-d['rttov14_cs_option2b_Surface2spaceTr'][:-1,item-1],
                       linewidth=1.2,color='magenta', marker='o',markersize = 2, mfc='none')
        
        plt.title( freq_title[i]+'('+str(options['chan_plots_titles'][item-1])+' GHz)', fontsize='12', fontweight='bold')

        ax1.set_xticks(np.arange(0,60,10))
        ax1.set_ylabel(r'Transmittance', color='k')
        ax1.set_xlabel(r'Sat. zenith angle [$^\circ$]', color='k')
        ax1.grid('true')

    # External legend; line labels
    llabels = ['RTTOV14 (coef13-option1)','RTTOV14 (coef13-option2b)', 'RTTOV14 (option1-option2b)']

    fig.legend([pl1, pl2, pl3], #, pl3], #, pl2, pl3, pl4, pl5],                   # The line objects
                labels        = llabels,      # The labels for each line
                loc           = "lower right",# Position of legend
                borderaxespad = 0.5,          # Small spacing around legend box
                title         = "RTTOV version:")     # Title for the legend
    
    fig.suptitle(title+' difference',fontweight='bold' )
    plt.tight_layout()
    plt.subplots_adjust(bottom=0.2, wspace=0.33)

    plt.savefig(options['path']+'/'+'Surf2TOA_TotalTransmittance_coeffsDifference.png')
    
    return
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
def TransSensor(d, title, options, rttov14_halflevels, rttov13_fulllevels):
    """
    -------------------------------------------------------------
    As in ClearSkySensorOver but for 5 differemt types and ARTS
    or RTTOV
    -------------------------------------------------------------
    OUT    name.png   Plot stored at the given path
    IN     d          RTTOV  i = 1,...,5
           y1         Lower limit at y axis
           y2         Upper limit at y axis
           ybins      yticks, e.g., np.arange(y_lo,y_up,increment)
           options    Figure dependent options
    -------------------------------------------------------------
    """
    plt.matplotlib.rc('font', family='serif', size = 12)

    #yl1 = [290,280,240,255,265]
    #yl2 = [300,290,255,265,280]
    
    freq_title = []
    for i in options['nchan_plots']:
        freq_title.append(r'CH'+str(i))

    fig = plt.figure(figsize=(options['width'],options['height']))
    for i, item in enumerate(options['nchan_plots']):
        ax1 = fig.add_subplot(options['nrows'],options['ncols'],i+1)
        pl1 = ax1.plot(np.arange(0,55,5),d['rttov14_cs_Surface2spaceTr'][:-1,item-1],linewidth=1.2,color='darkred', marker='o',markersize = 2, mfc='none')
        #pl2 = ax1.plot(np.arange(0,55,5),d['rttov14_as_Surface2spaceTr_FALSE'][:,item-1],linewidth=1.2,color='red', marker='o',markersize = 2, mfc='none')
        #pl3 = ax1.plot(np.arange(0,55,5),d['rttov14_as_Surface2spaceTr_TRUE'][:,item-1],linewidth=1.2,color='magenta', marker='o',markersize = 2, mfc='none')
        pl4 = ax1.plot(np.arange(0,55,5),d['rttov13_cs_Surface2spaceTr'][:-1,item-1],linewidth=1.2,color='darkblue', marker='o',markersize = 2, mfc='none')
        pl5 = ax1.plot(np.arange(0,55,5),(d['rttov13_as_Surface2spaceTr'][:-1,item-1]),linewidth=1.2,color='blue', marker='o',markersize = 2, mfc='none')

        plt.title( freq_title[i]+'('+str(options['chan_plots_titles'][item-1])+' GHz)', fontsize='12', fontweight='bold')

        ax1.set_xticks(np.arange(0,60,10))
        ax1.set_ylabel(r'Transmittance', color='k')
        ax1.set_xlabel(r'Sat. zenith angle [$^\circ$]', color='k')
        ax1.grid('true')

    # External legend; line labels
    llabels = ['RTTOV14', 'RTTOV13.2', 'RTTOV13.2 (RTTOV-SCATT)']#, 'RTTOV 13.2 (as)']

    fig.suptitle(title, fontweight='bold' )
    plt.tight_layout()
    plt.subplots_adjust(bottom=0.2, wspace=0.33) # top=0.899)

    fig.legend([pl1, pl4, pl5], #, pl3], #, pl2, pl3, pl4, pl5],                   # The line objects
                labels        = llabels,      # The labels for each line
                loc           = "lower center",# Position of legend
                ncol          = 3, 
                borderaxespad = 0.5,          # Small spacing around legend box
                title         = "RTTOV version:", 
                fancybox=False, shadow=False)     # Title for the legend

    plt.savefig(options['path']+'/'+'Surf2TOA_TotalTransmittance_overview.png')
    #plt.close()
    
    
    # Check if  rttov14_as_Surface2spaceTr_FALSE == rttov14_as_Surface2spaceTr_TRUE
    fig = plt.figure(figsize=(options['width'],options['height']))
    for i, item in enumerate(options['nchan_plots']):
        ax1 = fig.add_subplot(options['nrows'],options['ncols'],i+1)
        pl1 = ax1.plot(np.arange(0,55,5),d['rttov14_as_Surface2spaceTr_FALSE'][:-1,item-1]-d['rttov14_as_Surface2spaceTr_TRUE'][:-1,item-1],
                       linewidth=1.2,color='darkred', marker='x',markersize = 4, mfc='none')
        pl2 = ax1.plot(np.arange(0,55,5),d['rttov14_cs_Surface2spaceTr'][:-1,item-1]-d['rttov14_as_Surface2spaceTr_TRUE'][:-1,item-1],
                       linewidth=1.2,color='red', marker='o',markersize = 2, mfc='none')
        plt.title( freq_title[i]+'('+str(options['chan_plots_titles'][item-1])+' GHz)', fontsize='12', fontweight='bold')

        ax1.set_xticks(np.arange(0,60,10))
        ax1.set_ylabel(r'Transmittance', color='k')
        ax1.set_xlabel(r'Sat. zenith angle [$^\circ$]', color='k')
        ax1.grid('true')

    # External legend; line labels
    llabels = ['RTTOV14(scatt.F) - RTTOV14(scatt.True)', 'RTTOV14(cs) - RTTOV14(scatt.True)']#, 'RTTOV 13.2 (as)']

    fig.suptitle(title+' RTTOV14 internal consistency' ,fontweight='bold' )
    plt.tight_layout()
    plt.subplots_adjust(bottom=0.2, wspace=0.33)

    fig.legend([pl1, pl2], #, pl3], #, pl2, pl3, pl4, pl5],                   # The line objects
                labels        = llabels,      # The labels for each line
                loc           = "lower center",# Position of legend
                borderaxespad = 0.5,          # Small spacing around legend box
                title         = "RTTOV version:",
                fancybox=False, shadow=False, ncol=3)     # Title for the legend
    

    plt.savefig(options['path']+'/'+'Surf2TOA_TotalTransmittance_checkRTTOV14consistency.png')
    
    return
#-----------------------------------------------------------------------------
def Plot4SingleLayerTauSensor(dcs,d1,d2,d3,x1,x2,y1,y2,ititle,name,options,llabels):
    """
    -------------------------------------------------------------
    Comparison between RTTOV and ARTS for clear sky
    -------------------------------------------------------------
    OUT    name.png   Plot stored at the given path
    IN     
           d2         RTTOV clear sky data
           d2         RTTOV clear sky data
           x1         Lower limit at x axis
           x2         Upper limit at x axis
           y1         Lower limit at y axis
           y2         Upper limit at y axis
           xlabel     xlabel
           ititle     Title of figure
           name       Name to store the figure
           path       Path to store the figure
    -------------------------------------------------------------
    """
    plt.matplotlib.rc('font', family='serif', size = 12)

    freq_title = []
    for i in options['nchan_plots']:
        freq_title.append(r'CH'+str(i))

    fig = plt.figure(figsize=(options['width'],options['height']))
    for i, item in enumerate(options['nchan_plots']):
        ax1 = fig.add_subplot(options['nrows'],options['ncols'],i+1)
               
        pl1cs = ax1.plot( np.flipud(dcs[item-1,:]), np.arange(dcs.shape[1]), linewidth=1.0, color='magenta', marker='o', markersize = 2 )
        pl1 = ax1.plot( np.flipud(d1[item-1,:]), np.arange(d1.shape[1]), linewidth=1.0, color='black', marker='o', markersize = 2 )
        pl2 = ax1.plot( np.flipud(d2[item-1,:]), np.arange(d2.shape[1]), linewidth=0.5, color='darkred', marker='x',markersize = 2)
        pl3 = ax1.plot( np.flipud(d3[item-1,:]), np.arange(d3.shape[1]), linewidth=0.5, color='darkblue', marker='x',markersize = 2)

        #plt.title( freq_title[i], fontsize='12', fontweight='bold')
        plt.title( freq_title[i]+'('+str(options['chan_plots_titles'][item-1])+' GHz)', fontsize='12', fontweight='bold')

        ax1.set_ylabel(r'model level', color='k')
        ax1.set_xlabel(options['xlabel'], color='k')
        ax1.grid('true')
        ax1.set_ylim(y1,y2)

        if len(x2) > 1:
            ax1.set_xlim(x1,x2[i])
        else:
            ax1.set_xlim(x1,x2)
        
    # External legend; line labels
    fig.legend([pl1cs,pl1, pl2, pl3],                   # The line objects
               labels        = llabels,      # The labels for each line
               loc           = "lower right",# Position of legend
               borderaxespad = 0.5,          # Small spacing around legend box
               title         = "RTTOV version:")     # Title for the legend

    fig.suptitle(str(ititle), fontweight='bold' )
    plt.tight_layout()
    plt.subplots_adjust(bottom=0.2, wspace=0.33)
    plt.savefig(options['path']+'/'+str(name)+'.png')
    #plt.close()
    return

#-----------------------------------------------------------------------------
def HYDRO_CompSensor(d1,d2,alt,y1,y2,options):
    """
    -------------------------------------------------------------
    Comparison between RTTOV and ARTS in terms of their bulk
    properties
    -------------------------------------------------------------
    OUT    name.png   Plot stored at the given path
    IN     d1         ARTS clear sky data
           d2         RTTOV clear sky data
           ikey       Corresponding bulk property
           x1         Lower limit at x axis
           x2         Upper limit at x axis
           y1         Lower limit at y axis
           y2         Upper limit at y axis
           options    Figure dependent options
    -------------------------------------------------------------
    """
    plt.matplotlib.rc('font', family='serif', size = 18)
    
    freq_title = []
    for i in options['nchan_plots']:
        freq_title.append(r'CH'+str(i))

    hydro = ['Rain', 'Snow', 'Graupel', 'Cloud', 'Ice']
    for ihydro in [0,1,2,4]:
        
        fig = plt.figure(figsize=(options['width'],options['height']))
        for i, item in enumerate(options['nchan_plots']):
            ax1 = fig.add_subplot(options['nrows'],options['ncols'],i+1)        
    
            pl1 = ax1.plot(d1[item-1,:,ihydro], alt,  color='darkblue', linestyle='-', linewidth = 1.2)
            pl2 = ax1.plot(d2[:, item-1,ihydro], alt,  color='darkred', marker='o', linestyle='none', markersize = 5)

            #pl2 = ax1.plot(d1[item-1,:,1], alt,  color='darkred', linestyle='-', linewidth = 1.2)
            #pl3 = ax1.plot(d1[item-1,:,2], alt,  color='darkgreen', linestyle='-', linewidth = 1.2)
            ##pl4 = ax1.plot(d1[item-1,:,3], alt,  color='darkviolet', linestyle='-', linewidth = 1.2) IGNORE CLOUD
            #pl5 = ax1.plot(d1[item-1,:,4], alt, color='darkviolet', linestyle='-', linewidth = 1.2)
    
            #pl2 = ax1.plot(d2[item-1,:], alt, linewidth=0.75, color='darkred', marker='x', markersize = 5)
    
            ax1.grid('true')
            plt.title( freq_title[i]+'('+str(options['chan_plots_titles'][item-1])+' GHz)', fontsize='12', fontweight='bold')
            plt.xlabel(options['xlabel'])
            plt.ylabel(r'Pressure [hPa]')
            plt.ylim(y1,y2)
    
        fig.suptitle(options['title']+' ('+hydro[ihydro]+')' ,fontweight='bold' )
        plt.tight_layout()
        plt.subplots_adjust(bottom=0.25, wspace=0.5)
        
        # External legend; line labels
        llabels = ['RTTOV14','RTTOV13']
        fig.legend([pl1, pl2],            # The line objects
                    labels        = llabels,      # The labels for each line
                   loc           = "lower center",# Position of legend
                   borderaxespad = 0.5,          # Small spacing around legend box
                    title         = "RTTOV version:",     # Title for the legend
                    fancybox=False, shadow=False, ncol=2)     # Title for the legend
    
        plt.savefig(options['path']+'/'+options['name']+hydro[ihydro]+'.png')
        
    


    return


#-----------------------------------------------------------------------------
def BulkCompSensor(d1,d2,alt,y1,y2,options,llabels):
    """
    -------------------------------------------------------------
    Comparison between RTTOV and ARTS in terms of their bulk
    properties
    -------------------------------------------------------------
    OUT    name.png   Plot stored at the given path
    IN     d1         ARTS clear sky data
           d2         RTTOV clear sky data
           ikey       Corresponding bulk property
           x1         Lower limit at x axis
           x2         Upper limit at x axis
           y1         Lower limit at y axis
           y2         Upper limit at y axis
           options    Figure dependent options
    -------------------------------------------------------------
    """
    plt.matplotlib.rc('font', family='serif', size = 18)
    
    freq_title = []
    for i in options['nchan_plots']:
        freq_title.append(r'CH'+str(i))

    fig = plt.figure(figsize=(options['width'],options['height']))
    for i, item in enumerate(options['nchan_plots']):
        ax1 = fig.add_subplot(options['nrows'],options['ncols'],i+1)        

        pl1 = ax1.plot(d1[item-1,:], alt, linewidth=0.75, color='darkblue', marker='o', markersize = 5)
        pl2 = ax1.plot(d2[item-1,:], alt, linewidth=0.75, color='darkred', marker='x', markersize = 5)

        ax1.grid('true')
        plt.title( freq_title[i]+'('+str(options['chan_plots_titles'][item-1])+' GHz)', fontsize='12', fontweight='bold')
        plt.xlabel(options['xlabel'])
        plt.ylabel(r'Pressure [hPa]')
        plt.ylim(y1,y2)

    fig.suptitle(options['title'] ,fontweight='bold' )
    plt.tight_layout()
    plt.subplots_adjust(bottom=0.25, wspace=0.5)
    
    # External legend; line labels
    #llabels = ['RTTOV'] #,'RTTOV 3.2']
    fig.legend([pl1, pl2],            # The line objects
                labels        = llabels,      # The labels for each line
                loc           = "lower center",# Position of legend
                borderaxespad = 0.5,          # Small spacing around legend box
                title         = "RTTOV version:", ncol=2,     # Title for the legend
                fancybox=False, shadow=False)     # Title for the legend

    plt.savefig(options['path']+'/'+options['name']+'.png')
    
    


    return



#-----------------------------------------------------------------------------
def plot_pressuregrid(plevels, option):
    
    fig, axs = plt.subplots(1, 1, figsize=(8, 10))  # 2 rows, 1 column of subplots
    # nan legends
    plt.axhline(np.nan, color='darkred', linestyle='-', label='IFS half level')
    plt.axhline(np.nan, color='darkblue', linestyle='-', label='IFS half levels average')
    plt.axhline(np.nan, color='k',  marker='o', linestyle='None', label='IFS full levels')
    plt.legend()
    for level in plevels['rttov14_halflevels']:
        plt.axhline(level/100, color='darkred', linestyle='-', label='IFS half level')
    for level in plevels['betwee_halflevels']:
        plt.axhline(level/100, color='darkblue', linestyle='-', label='IFS half levels average')
    for level in plevels['rttov13_fulllevels']:
        plt.axhline(level/100, color='k', marker='o', linestyle='None', label='IFS full levels')
    plt.ylabel('Pressure levels (hPa)')
    plt.ylim(900,1030)
    plt.savefig(option['path']+'/'+'PressureLevels.png')
    
    return
    
#-----------------------------------------------------------------------------
#-------------------------------------------------------------------------
def run_IFS_plots(options_path, ncfile):
    
    
    #plot_dir   = '/home/vito.galligani/mountpoint/Work/RTTOV/Plots/'
    #options_path = os.path.join(plot_dir,datetime.now().strftime('%d%m%Y'))
    #
    #ncfile='/home/vito.galligani/Work/RTTOV/main/AIRS_processed_W1.nc'
 
    #------------------------------------------------------------------------------
    #------------------------------------------------------------------------------
    
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

    A['longitude'] = data.variables['longitude'][:]       # orography
    A['lats']      = data.variables['lats'][:]       # orography

    A['gwp'] = data.variables['gwp'][:]       # orography
    A['swp'] = data.variables['swp'][:]       # orography
    A['rwp'] = data.variables['rwp'][:]       # orography
    A['iwp'] = data.variables['iwp'][:]       # orography
    
    
    crs_latlon = ccrs.PlateCarree()

    # FIND CLEAR SKY PROFILE! == A['xwp'][36] in (350.53,17.05524 ) AND (38.61, -5.9 ) 
    lonlon = data.variables['longitude'][:]  
    latlat = data.variables['lats'][:]  
    LON  = lonlon[np.where( (A['rwp'][:] == 0) & (lonlon > 300) & (latlat > -40) & (latlat < -20))]     
    LAT  = latlat[np.where( (A['rwp'][:] == 0) & (lonlon > 300) & (latlat > -40) & (latlat < -20))]     
    
    fig, ax1 = plt.subplots(figsize=(5, 5), subplot_kw={'projection': ccrs.PlateCarree()})
    plt.plot( LON, LAT, 'ok', transform=crs_latlon) # vmax=1E1)
    ax1.coastlines(resolution='10m', edgecolor='g',linewidth=0.5, zorder=2)
    fig = plt.figure(figsize=(5,5))
    clr_ = ['k','r','b','darkred','darkblue','darkgreen','magenta', 'cyan']
    for i in range(len(LON)):
        plt.plot( LON[i], LAT[i],'o', color=clr_[i], label='index: '+str(i))
    plt.legend()
    # USE EITHER PROFILES 3
    ProfN = 5
    index = [np.where( (lonlon == LON[ProfN]) & (latlat == LAT[ProfN])) ][0][0][0]    


    # FIND CLEAR SKY PROFILE! in the tropics too TRY INDEX 3
    lonlon = data.variables['longitude'][:]  
    latlat = data.variables['lats'][:]  
    LON2  = lonlon[np.where( (A['rwp'][:] == 0) & (lonlon > 300) & (latlat > -15) & (latlat < 15))]     
    LAT2  = latlat[np.where( (A['rwp'][:] == 0) & (lonlon > 300) & (latlat > -15) & (latlat < 15))]     
    fig, ax1 = plt.subplots(figsize=(5, 5), subplot_kw={'projection': ccrs.PlateCarree()})
    plt.plot( LON2, LAT2, 'ok', transform=crs_latlon) # vmax=1E1)
    ax1.coastlines(resolution='10m', edgecolor='g',linewidth=0.5, zorder=2)
    ax1.set_extent([-60,0,-20,20], crs=crs_latlon) 
 
    fig = plt.figure(figsize=(5,5))
    clr_ = ['k','r','b','darkred','darkblue','darkgreen','magenta', 'cyan']
    for i in range(len(LON)):
        plt.plot( LON2[i], LAT2[i],'o', color=clr_[i], label='index: '+str(i))    
    plt.legend()
    # USE EITHER PROFILES 3
    ProfN = 3
    index2 = [np.where( (lonlon == LON2[ProfN]) & (latlat == LAT2[ProfN])) ][0][0][0]       

    fig = plt.figure(figsize = (8,4))  #[width, height]    
    ax1 = plt.subplot(111, projection=ccrs.PlateCarree())
    ax1.plot(lonlon[np.where( A['rwp'][:] == 0)], latlat[np.where( A['rwp'][:] == 0)], 'ok', markersize=2, transform=crs_latlon)
    ax1.set_extent([-180,0,-70,70], crs=crs_latlon) 
    ax1.coastlines(resolution='10m', edgecolor='g',linewidth=0.5, zorder=2)
    ax1.plot(lonlon[index], latlat[index], 'xr',  transform=crs_latlon)
    ax1.plot(lonlon[index2], latlat[index2], 'xb',  transform=crs_latlon)
    
     
    plt.matplotlib.rc('font', family='serif', size = 8)
    fig, ax1 = plt.subplots(figsize=(2,3)) 	#8.3 x 11.7

    #fig, ax1 = plt.subplots(figsize=(4,8)) 	#8.3 x 11.7
    line1, = ax1.plot( A['t'][index,:], A['pf'][index,:]/100, linewidth=1.2, linestyle='-', color='darkblue')
    ax1.plot( A['t'][index2,:], A['pf'][index2,:]/100, linewidth=1.2, linestyle='--', color='darkblue')
    ax1.set_xlabel('Temperature (K)', color='darkblue')
    ax1.set_ylabel('Pressure (hPa)')
    ax1.set_ylim([np.max(A['pf'][index,:])/100, 600])
    ax1.set_xlim([270, 302])
    # Set color of spines and tick labels for the first axis
    ax1.spines['bottom'].set_color(line1.get_color())
    ax1.tick_params(axis='x', colors=line1.get_color())
    ax1.xaxis.label.set_color(line1.get_color())
    plt.grid(True)
    
    # Creating a second y-axis sharing the same y-axis
    ax2 = ax1.twiny()
    line2, = ax2.plot( A['h20'][index,:], A['pf'][index,:]/100, linewidth=1.2, linestyle='-', color='darkred')
    ax2.plot( A['h20'][index2,:], A['pf'][index2,:]/100, linewidth=1.2, linestyle='--', color='darkred')
    ax2.set_xlabel('q (kg/kg)', color='darkred')
    ax2.set_xscale('log')
    ax2.set_ylim([np.max(A['pf'][index,:])/100, 600])    
    # Set color of spines and tick labels for the first axis
    ax2.spines['top'].set_color(line2.get_color())    
    ax2.tick_params(axis='x', colors=line2.get_color())
    ax2.xaxis.label.set_color(line2.get_color())
    ax2.set_xlim([0.0009, 0.02])
    plt.grid(True)

    plt.title('Test Profiles', fontsize=10) #(-58.9,-30.3)
    plt.tight_layout()
    plt.savefig(options_path+'/'+'ExampleProf_2profiles.png', bbox_inches='tight', dpi=300 )

    plt.matplotlib.rc('font', family='serif', size = 8)
    fig, ax1 = plt.subplots(figsize=(3,4)) 	#8.3 x 11.7
    ax1.plot( A['t'][index2,:], A['pf'][index2,:]/100, linewidth=1.2, linestyle='-', color='darkblue')
    ax1.set_xlabel('Temperature (K)', color='darkblue')
    ax1.set_ylabel('Pressure (hPa)')
    ax1.set_ylim([np.max(A['pf'][index,:])/100, 600])
    ax1.set_xlim([270, 302])
    # Set color of spines and tick labels for the first axis
    ax1.spines['bottom'].set_color(line1.get_color())
    ax1.tick_params(axis='x', colors=line1.get_color())
    ax1.xaxis.label.set_color(line1.get_color())
    plt.grid(True)
    
    # Creating a  y-axis sharing the same y-axis
    ax2 = ax1.twiny()
    ax2.plot( A['h20'][index2,:], A['pf'][index2,:]/100, linewidth=1.2, linestyle='-', color='darkred')
    ax2.set_xlabel('q (kg/kg)', color='darkred')
    ax2.set_xscale('log')
    ax2.set_ylim([np.max(A['pf'][index,:])/100, 600])    
    # Set color of spines and tick labels for the first axis
    ax2.spines['top'].set_color(line2.get_color())    
    ax2.tick_params(axis='x', colors=line2.get_color())
    ax2.xaxis.label.set_color(line2.get_color())
    ax2.set_xlim([0.0009, 0.02])
    plt.grid(True)

    #plt.title('Test Profiles', fontsize=10) #(-58.9,-30.3)
    plt.tight_layout()
    plt.savefig(options_path+'/'+'ExampleProf_1profiles.png', bbox_inches='tight', dpi=300 )
    

    #-----------------------------------------------------------------------
    # Graficamos
    threshold = 1e1

    crs_latlon = ccrs.PlateCarree()
    countries = cfeature.NaturalEarthFeature(category='cultural', name='admin_0_countries',scale='10m',facecolor='none')
    xticks = [-180, -150, -120,  -90,  -60,  -30,    0]
    yticks = [-70, -60, -40, -20,  -0,  20,  40,  60, 70]

    plt.matplotlib.rc('font', family='serif', size = 10)

    #-------
    #fig, ax1 = plt.subplots(figsize=(10, 5), subplot_kw={'projection': ccrs.PlateCarree()})
    fig = plt.figure(figsize = (8,4))  #[width, height]
    
    ax1 = plt.subplot(221, projection=ccrs.PlateCarree())
    ax1.set_facecolor('gray')
    ax1.scatter(A['longitude'], A['lats'], color='white', edgecolor='white', marker='o', s=0.1,
            alpha=1, transform=ccrs.PlateCarree(), zorder=0)
    pc = ax1.scatter(A['longitude'], A['lats'], c=A['gwp'], cmap='rainbow', marker='o', s=0.1,alpha=1,
                 norm=matplotlib.colors.LogNorm(vmin=1e-6, vmax=1e1),  transform=crs_latlon, zorder=1) # vmax=1E1)
    cbar = plt.colorbar(pc, ax=ax1, orientation='vertical', shrink=0.9,extend='max', spacing='proportional')#, pad=0.05, aspect=50)
    cbar.set_label('Graupel (kgm$^{-2}$)')
    cbar.set_ticks([1e-6,1e-5,1e-4,1e-3,1e-2,1e-1,1e0,1e1])
    cbar.ax.tick_params(labelsize=8) 
    ax1.coastlines(resolution='10m', edgecolor='g',linewidth=0.5, zorder=2)
    ax1.set_extent([-180,0,-70,70], crs=crs_latlon) 
    gl = ax1.gridlines(draw_labels=True,linewidth=0.5, color='gray', alpha=1, linestyle='--', xlocs=xticks, ylocs=yticks, zorder=0)
    gl.xformatter = LongitudeFormatter()
    gl.yformatter = LatitudeFormatter()
    gl.top_labels   = False
    gl.right_labels = False
    gl.xlines       = True
    gl.xlocator = mticker.FixedLocator([-180, -150, -120,  -90,  -60,  -30,    0])
    gl.xlabel_style =  {'family': 'monospace', 'size': 8}
    gl.ylabel_style = {'family': 'monospace', 'size': 8}

    ax1 = plt.subplot(222, projection=ccrs.PlateCarree())
    ax1.set_facecolor('gray')
    ax1.scatter(A['longitude'], A['lats'], color='white', edgecolor='white', marker='o', s=0.1,
            alpha=1, transform=ccrs.PlateCarree(), zorder=0)
    pc = ax1.scatter(A['longitude'], A['lats'], c=A['swp'], cmap='rainbow', marker='o', s=0.1,
                 norm=matplotlib.colors.LogNorm(vmin=1e-6, vmax=1e1),  transform=crs_latlon, zorder=1) # vmin=1E-6, vmax=1E1)
    cbar = plt.colorbar(pc, ax=ax1, orientation='vertical', shrink=0.9,extend='max', spacing='proportional')#, pad=0.05, aspect=50)
    cbar.set_label('Snow (kgm$^{-2}$)')
    cbar.set_ticks([1e-6,1e-5,1e-4,1e-3,1e-2,1e-1,1e0,1e1])
    cbar.ax.tick_params(labelsize=8) 
    ax1.coastlines(resolution='10m', edgecolor='g',linewidth=0.5, zorder=2)
    ax1.set_extent([-180,0,-70,70], crs=crs_latlon) 
    gl = ax1.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                       linewidth=1.2, color='gray', alpha=0.5, linestyle='--', xlocs=xticks, ylocs=yticks)
    gl.top_labels   = False
    gl.right_labels = False
    gl.xlines       = True
    gl.xlocator = mticker.FixedLocator([-180, -150, -120,  -90,  -60,  -30,    0])
    gl.xformatter = LongitudeFormatter()
    gl.yformatter = LatitudeFormatter()
    gl.xlabel_style =  {'family': 'monospace', 'size': 8}
    gl.ylabel_style = {'family': 'monospace', 'size': 8}
    
    
    ax1 = plt.subplot(223, projection=ccrs.PlateCarree())
    ax1.set_facecolor('gray')
    ax1.scatter(A['longitude'], A['lats'], color='white', edgecolor='white', marker='o', s=0.1,
            alpha=1, transform=ccrs.PlateCarree(), zorder=0)
    pc = ax1.scatter(A['longitude'], A['lats'], c=A['rwp'], cmap='rainbow', marker='o', s=0.1,
                 norm=matplotlib.colors.LogNorm(vmin=1e-6, vmax=1e1),  transform=crs_latlon, zorder=1) # vmin=1E-6, vmax=1E1)
    cbar = plt.colorbar(pc, ax=ax1, orientation='vertical', shrink=0.9,extend='max', spacing='proportional')#, pad=0.05, aspect=50)
    cbar.set_ticks([1e-6,1e-5,1e-4,1e-3,1e-2,1e-1,1e0,1e1])
    cbar.set_label('Rain (kgm$^{-2}$)')
    cbar.ax.tick_params(labelsize=8) 
    ax1.coastlines(resolution='10m', edgecolor='g',linewidth=0.5, zorder=2)
    ax1.set_extent([-180,0,-70,70], crs=crs_latlon) 
    gl = ax1.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                       linewidth=1.2, color='gray', alpha=0.5, linestyle='--', xlocs=xticks, ylocs=yticks)
    gl.top_labels   = False
    gl.right_labels = False
    gl.xlines       = True
    gl.xlocator = mticker.FixedLocator([-180, -150, -120,  -90,  -60,  -30,    0])
    gl.xformatter = LongitudeFormatter()
    gl.yformatter = LatitudeFormatter()
    gl.xlabel_style =  {'family': 'monospace', 'size': 8}
    gl.ylabel_style = {'family': 'monospace', 'size': 8}
    
    
    ax1 = plt.subplot(224, projection=ccrs.PlateCarree())
    ax1.set_facecolor('gray')
    ax1.scatter(A['longitude'], A['lats'], color='white', edgecolor='white', marker='o', s=0.1,
            alpha=1, transform=ccrs.PlateCarree(), zorder=0)
    pc = ax1.scatter(A['longitude'], A['lats'], c=A['iwp'], cmap='rainbow', marker='o', s=0.1,
                 norm=matplotlib.colors.LogNorm(vmin=1e-6, vmax=1e1),  transform=crs_latlon, zorder=1) # vmin=1E-6, vmax=1E1)
    cbar = plt.colorbar(pc, ax=ax1, orientation='vertical', shrink=0.9,extend='max', spacing='proportional')#, pad=0.05, aspect=50)
    cbar.set_ticks([1e-6,1e-5,1e-4,1e-3,1e-2,1e-1,1e0,1e1])
    cbar.set_label('Ice (kgm$^{-2}$)')
    cbar.ax.tick_params(labelsize=8) 
    ax1.coastlines(resolution='10m', edgecolor='g',linewidth=0.5, zorder=2)
    ax1.set_extent([-180,0,-70,70], crs=crs_latlon) 
    gl = ax1.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                       linewidth=1.2, color='gray', alpha=0.5, linestyle='--', xlocs=xticks, ylocs=yticks)
    gl.top_labels   = False
    gl.right_labels = False
    gl.xlines       = True
    gl.xlocator = mticker.FixedLocator([-180, -150, -120,  -90,  -60,  -30,    0])
    gl.xformatter = LongitudeFormatter()
    gl.yformatter = LatitudeFormatter()
    gl.xlabel_style =  {'family': 'monospace', 'size': 8}
    gl.ylabel_style = {'family': 'monospace', 'size': 8}

       
    fig.suptitle('IFS Hydrometeor water paths' ,fontweight='bold' )
    plt.savefig(options_path+'/'+'IFSwaterpaths.png', bbox_inches='tight', dpi=300 )


    return index, index2

#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
