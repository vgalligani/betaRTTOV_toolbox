#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
-----------------------------------------------------------------
@purpose : Analyse rttov 13.4 outputs of REALISITC IFS 
@author  : Vito Galligani
@email   : vito.galligani@gmail.com
@condaenv: conda activate RTTOV14
-----------------------------------------------------------------
@main    : Combine into netcdf all simulations 
          
           Remember to open remote session:
               ssh -L 6000:yakaira:22 vito.galligani@portal.cima.fcen.uba.ar
               nohup python -m spyder_kernels.console -matplotlib='inline' -f=./remotemachine.json &
"""
#################################################
import numpy          as np
import Tools2Read     as T2R
import xarray         as xr
import time
#################################################
def CombineRTTOV(sList):
    
    #------------ RTTOV14
    out_path = '/home/vito.galligani/datosmunin3/RTTOV_betaTest/rttov14/Output/'
    
    prttov_as = out_path+'TestAllSky/atms'
    prttov_cs   = out_path+'TestClearSky/atms/TestClearSky'
    prttov_cs1  = out_path+'TestClearSky/atms/TestClearSky_option1'
    prttov_cs2b = out_path+'TestClearSky/atms/TestClearSky_option2b'
    
    r14_as = []
    r14_cs = []
    r14_cs_option1 = []
    r14_cs_option2b = []
    
    start_time = time.time()
    print('start: as rttov14')
    for i in sList:  
        #print('RTTOV14 as: '+str(i))
        rttov_asfiles   = [f"{prttov_as}/output_tb_{x}_atm_em1.0rttov14_Profile_W1_{i}.dat" for x in range(0, 55, 5) ]
        f_scatt_01 = rttov_asfiles
        scatt_01   = T2R.MultiRow(f_scatt_01)
        r14_as.append(scatt_01)
    #-----------------------------------------------------------------------------
    #exp = 'TestAllSky'
    ds  = xr.Dataset( {  "rttov14_tb_as": ( ("observation", "theta", "chan"), r14_as, {'units':'K'})  }   )
    ds.to_netcdf('/home/vito.galligani/Work/RTTOV/EXPncfiles/RTTOV14_ATMS_asky.nc', 'w')
    del r14_as
    
    end_time = time.time()
    # Calculate elapsed time
    elapsed_time = end_time - start_time
    print("Elapsed time:", elapsed_time, "seconds")    
    
    start_time = time.time()
    print('start: cs rttov14')
    for i in sList:
        #print('RTTOV14 cs: '+str(i))
        rttov_csfiles   = [f"{prttov_cs}/output_tb_{x}_atm_em1.0rttov14_Profile_W1_{i}.dat" for x in range(0, 55, 5) ]
        frttov_tb  = rttov_csfiles
        rttov_cs   = T2R.MultiRow(frttov_tb)
        r14_cs.append(rttov_cs)
    #-----------------------------------------------------------------------------
    #exp = 'TestClearSky'
    ds  = xr.Dataset( {  "rttov14_tb_cs": ( ("observation", "theta", "chan"), r14_cs, {'units':'K'})  }   )
    ds.to_netcdf('/home/vito.galligani/Work/RTTOV/EXPncfiles/RTTOV14_ATMS_csky.nc', 'w')
    del r14_cs
    end_time = time.time()
    # Calculate elapsed time
    elapsed_time = end_time - start_time
    print("Elapsed time:", elapsed_time, "seconds")
    

    
    start_time = time.time()
    print('start: cs opt1 rttov14')
    for i in sList:
        #print('RTTOV14 cs 1: '+str(i))
        rttov_cs1files  = [f"{prttov_cs1}/output_tb_{x}_atm_em1.0rttov14_Profile_W1_{i}.dat" for x in range(0, 55, 5) ]
        frttov_tb  = rttov_cs1files
        rttov_cs1   = T2R.MultiRow(frttov_tb)
        r14_cs_option1.append(rttov_cs1)    
    #-----------------------------------------------------------------------------
    #exp = 'TestClearSky_option1'
    ds  = xr.Dataset( {  "rttov14_tb_cs_option1": ( ("observation", "theta", "chan"), r14_cs_option1, {'units':'K'})  }   )
    ds.to_netcdf('/home/vito.galligani/Work/RTTOV/EXPncfiles/RTTOV14_ATMS_csky_option1.nc', 'w')
    del r14_cs_option1
    end_time = time.time()
    # Calculate elapsed time
    elapsed_time = end_time - start_time
    print("Elapsed time:", elapsed_time, "seconds")

    start_time = time.time()
    print('start: cs opt2 rttov14')
    for i in sList:
        #print('RTTOV14 cs 2b: '+str(i))
        rttov_cs2bfiles = [f"{prttov_cs2b}/output_tb_{x}_atm_em1.0rttov14_Profile_W1_{i}.dat" for x in range(0, 55, 5) ]
        frttov_tb  = rttov_cs2bfiles
        rttov_cs2   = T2R.MultiRow(frttov_tb)
        r14_cs_option2b.append(rttov_cs2)  
    #-----------------------------------------------------------------------------
    #exp = 'TestClearSky_option2n'
    ds  = xr.Dataset( {  "rttov14_tb_cs_option2b": ( ("observation", "theta", "chan"), r14_cs_option2b, {'units':'K'})  }   )
    ds.to_netcdf('/home/vito.galligani/Work/RTTOV/EXPncfiles/RTTOV14_ATMS_csky_option2b.nc', 'w')
    del r14_cs_option2b
    end_time = time.time()
    # Calculate elapsed time
    elapsed_time = end_time - start_time
    print("Elapsed time:", elapsed_time, "seconds")

    #------------ RTTOV13
    out_path = '/home/vito.galligani/datosmunin3/RTTOV_betaTest/rttov13/Output/'
    prttov_as = out_path+'TestAllSky/atms'
    prttov_cs   = out_path+'TestClearSky/atms'        
        
    r13_as = []
    r13_cs = []

    start_time = time.time()
    print('start: as rttov13')
    for i in sList:
        #print('RTTOV13 as: '+str(i))
        rttov_asfiles   = [f"{prttov_as}/output_tb_{x}atm_em1.0rttov132_Profile_W1_{i}.dat" for x in range(0, 55, 5) ]
        f_scatt_01 = rttov_asfiles
        scatt_01   = T2R.MultiRow(f_scatt_01)
        r13_as.append(scatt_01)
    #-----------------------------------------------------------------------------
    #exp = RTTOV13 'TestAllSky'
    ds  = xr.Dataset( {  "rttov13_tb_as": ( ("observation", "theta", "chan"), r13_as, {'units':'K'})  }   )
    ds.to_netcdf('/home/vito.galligani/Work/RTTOV/EXPncfiles/RTTOV13_ATMS_asky.nc', 'w')
    del r13_as
    end_time = time.time()
    # Calculate elapsed time
    elapsed_time = end_time - start_time
    print("Elapsed time:", elapsed_time, "seconds")

    start_time = time.time()
    print('start: cs rttov13')
    for i in sList:
        #print('RTTOV13 cs: '+str(i))
        rttov_csfiles   = [f"{prttov_cs}/output_tb_{x}_atm_em1.0rttov132_Profile_W1_{i}.dat" for x in range(0, 55, 5) ]
        frttov_tb  = rttov_csfiles
        rttov_cs   = T2R.MultiRow(frttov_tb)
        r13_cs.append(rttov_cs)
    #-----------------------------------------------------------------------------
    #exp = RTTOV13 'TestClearSky'
    ds  = xr.Dataset( {  "rttov13_tb_cs": ( ("observation", "theta", "chan"), r13_cs, {'units':'K'})  }   )
    ds.to_netcdf('/home/vito.galligani/Work/RTTOV/EXPncfiles/RTTOV13_ATMS_csky.nc', 'w')
    del r13_cs 
    end_time = time.time()
    # Calculate elapsed time
    elapsed_time = end_time - start_time
    print("Elapsed time:", elapsed_time, "seconds")
        
                
    return 

def CombineRTTOV14(sList, prttov, dataoutputfile):
        
    r14 = []
    
    start_time = time.time()
    print('start: as rttov14')
    for i in sList:  
        rttov_files   = [f"{prttov_as}/output_tb_{x}_atm_em1.0rttov14_Profile_W1_{i}.dat" for x in range(0, 55, 5) ]
        data   = T2R.MultiRow(rttov_files)
        r14.append(data)
        
    ds  = xr.Dataset( {  "r14": ( ("observation", "theta", "chan"), r14, {'units':'K'})  }   )
    ds.to_netcdf('/home/vito.galligani/Work/RTTOV/EXPncfiles/'+dataoutputfile, 'w')
    del r14
    
    elapsed_time = time.time() - start_time
    print("Elapsed time:", elapsed_time, "seconds")    
                
    return 

def CombineRTTOV13a(sList, prttov, dataoutputfile):
        
    r14 = []
    
    start_time = time.time()
    print('start: as rttov14')
    for i in sList:  
        rttov_files   = [f"{prttov_as}/output_tb_{x}_atm_em1.0rttov132_Profile_W1_{i}.dat" for x in range(0, 55, 5) ]
        data   = T2R.MultiRow(rttov_files)
        r14.append(data)
        
    ds  = xr.Dataset( {  "r14": ( ("observation", "theta", "chan"), r14, {'units':'K'})  }   )
    ds.to_netcdf('/home/vito.galligani/Work/RTTOV/EXPncfiles/'+dataoutputfile, 'w')
    del r14
    
    elapsed_time = time.time() - start_time
    print("Elapsed time:", elapsed_time, "seconds")    
                
    return 

def CombineRTTOV13b(sList, prttov, dataoutputfile):
        
    r14 = []
    
    start_time = time.time()
    print('start: as rttov14')
    for i in sList:  
        rttov_files   = [f"{prttov_as}/output_tb_{x}atm_em1.0rttov132_Profile_W1_{i}.dat" for x in range(0, 55, 5) ]
        data   = T2R.MultiRow(rttov_files)
        r14.append(data)
        
    ds  = xr.Dataset( {  "r14": ( ("observation", "theta", "chan"), r14, {'units':'K'})  }   )
    ds.to_netcdf('/home/vito.galligani/Work/RTTOV/EXPncfiles/'+dataoutputfile, 'w')
    del r14
    
    elapsed_time = time.time() - start_time
    print("Elapsed time:", elapsed_time, "seconds")    
                
    return 



#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
nchan      = 22
nzenith    = 12
dummy      = np.zeros((nzenith,nchan))          
specList1  = [str(i) for i in np.arange(1,21998,1)]   # here the last profile is 22000

#- RTTOV14
out_path = '/home/vito.galligani/datosmunin3/RTTOV_betaTest/rttov14/Output/'

prttov_as = out_path+'TestAllSky/atms'
CombineRTTOV14(specList1, prttov_as, 'RTTOV14_ATMS_asky.nc') &

prttov_cs   = out_path+'TestClearSky/atms/TestClearSky'
CombineRTTOV14(specList1, prttov_cs, 'RTTOV14_ATMS_csky.nc') &

prttov_cs1  = out_path+'TestClearSky/atms/TestClearSky_option1'
CombineRTTOV14(specList1, prttov_cs1, 'RTTOV14_ATMS_csky_option1.nc') &

prttov_cs2b = out_path+'TestClearSky/atms/TestClearSky_option2b'
CombineRTTOV14(specList1, prttov_cs2b, 'RTTOV14_ATMS_csky_option2b.nc') &

#- RTTOV13
out_path = '/home/vito.galligani/datosmunin3/RTTOV_betaTest/rttov13/Output/'

prttov_as = out_path+'TestAllSky/atms'
CombineRTTOV13b(specList1, prttov_as, 'RTTOV13_ATMS_asky.nc') &

prttov_cs   = out_path+'TestClearSky/atms'        
CombineRTTOV13a(specList1, prttov_cs, 'RTTOV13_ATMS_csky.nc') & ^
  








