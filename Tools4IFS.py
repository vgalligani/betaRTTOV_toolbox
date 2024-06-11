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
import multiprocessing
import re

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

#-----------------------------------------------------------------------------
def CombineRTTOV_fast_clear2b(i):
    
    
    """Function to run in parallel"""
    print(f'Worker {i} /11 started')


    out_path = '/home/vito.galligani/datosmunin3/RTTOV_betaTest/rttov14/Output/'
    prttov = out_path+'TestClearSky/atms/TestClearSky_option2b'
    dataoutputfile = 'rttov14_cs_option2b'
    angles    = np.arange(0,55,5)
    x = angles[i]

    print('start: as rttov14')
    start_time = time.time()
    rttov_asfiles   = T2R.MultiRow(sorted([f"{prttov}/output_tb_{x}_atm_em1.0rttov14_Profile_W1_{i}.dat" for i in range(1, 21998, 1) ]))  #21998
    elapsed_time = time.time() - start_time
    print("Elapsed time:", elapsed_time, "seconds")        
    ds  = xr.Dataset( {  "rttov14": ( ("observation", "chan"), rttov_asfiles, {'units':'K'})  }   )
    ds.to_netcdf('/home/vito.galligani/Work/RTTOV/EXPncfiles/'+dataoutputfile+str(x)+'.nc', 'w')
    
    print(f'Worker {i} /11 finished')
    
    return

#-----------------------------------------------------------------------------
def CombineRTTOV_fast_clear1(i):
    
    
    """Function to run in parallel"""
    print(f'Worker {i} /11 started')


    out_path = '/home/vito.galligani/datosmunin3/RTTOV_betaTest/rttov14/Output/'
    prttov = out_path+'TestClearSky/atms/TestClearSky_option1'
    dataoutputfile = 'rttov14_cs_option1'
    angles    = np.arange(0,55,5)
    x = angles[i]

    print('start: as rttov14')
    start_time = time.time()
    rttov_asfiles   = T2R.MultiRow(sorted([f"{prttov}/output_tb_{x}_atm_em1.0rttov14_Profile_W1_{i}.dat" for i in range(1, 21998, 1) ]))  #21998
    elapsed_time = time.time() - start_time
    print("Elapsed time:", elapsed_time, "seconds")        
    ds  = xr.Dataset( {  "rttov14": ( ("observation", "chan"), rttov_asfiles, {'units':'K'})  }   )
    ds.to_netcdf('/home/vito.galligani/Work/RTTOV/EXPncfiles/'+dataoutputfile+str(x)+'.nc', 'w')
    
    print(f'Worker {i} /11 finished')
    
    return

#-----------------------------------------------------------------------------
def CombineRTTOV_fast_clear(i):
    
    
    """Function to run in parallel"""
    print(f'Worker {i} /11 started')


    out_path = '/home/vito.galligani/datosmunin3/RTTOV_betaTest/rttov14/Output/'
    prttov = out_path+'TestClearSky/atms/TestClearSky'
    dataoutputfile = 'rttov14_cs_'
    angles    = np.arange(0,55,5)
    x = angles[i]

    print('start: as rttov14')
    start_time = time.time()
    rttov_asfiles   = T2R.MultiRow(sorted([f"{prttov}/output_tb_{x}_atm_em1.0rttov14_Profile_W1_{i}.dat" for i in range(1, 21998, 1) ]))  #21998
    elapsed_time = time.time() - start_time
    print("Elapsed time:", elapsed_time, "seconds")        
    ds  = xr.Dataset( {  "rttov14": ( ("observation", "chan"), rttov_asfiles, {'units':'K'})  }   )
    ds.to_netcdf('/home/vito.galligani/Work/RTTOV/EXPncfiles/'+dataoutputfile+str(x)+'.nc', 'w')
    
    print(f'Worker {i} /11 finished')
    
    return
#-----------------------------------------------------------------------------
def CombineRTTOV_fast_clear13(i):
    
    
    """Function to run in parallel"""
    print(f'Worker {i} /11 started')

    #------------ RTTOV13
    out_path = '/home/vito.galligani/datosmunin3/RTTOV_betaTest/rttov13/Output/'
    prttov   = out_path+'TestClearSky/atms'        
    dataoutputfile = 'rttov13_cs_'
    angles    = np.arange(0,55,5)
    x = angles[i]

    print('start: cs rttov13')
    start_time = time.time()
    rttov_asfiles   = T2R.MultiRow(sorted([f"{prttov}/output_tb_{x}_atm_em1.0rttov132_Profile_W1_{i}.dat" for i in range(1, 21998, 1) ]))  #21998
    elapsed_time = time.time() - start_time
    print("Elapsed time:", elapsed_time, "seconds")        
    ds  = xr.Dataset( {  "rttov14": ( ("observation", "chan"), rttov_asfiles, {'units':'K'})  }   )
    ds.to_netcdf('/home/vito.galligani/Work/RTTOV/EXPncfiles/'+dataoutputfile+str(x)+'.nc', 'w')
    
    print(f'Worker {i} /11 finished')
    

    return

#-----------------------------------------------------------------------------
def CombineRTTOV_fast_scatt13(i):
    
    
    """Function to run in parallel"""
    print(f'Worker {i} /11 started')

    #------------ RTTOV13
    out_path = '/home/vito.galligani/datosmunin3/RTTOV_betaTest/rttov13/Output/'
    prttov   = out_path+'TestAllSky/atms'        
    dataoutputfile = 'rttov13_as_'
    angles    = np.arange(0,55,5)
    x = angles[i]

    print('start: as rttov13')
    start_time = time.time()
    rttov_asfiles   = T2R.MultiRow(sorted([f"{prttov}/output_tb_{x}atm_em1.0rttov132_Profile_W1_{i}.dat" for i in range(1, 21998, 1) ]))  #21998
    elapsed_time = time.time() - start_time
    print("Elapsed time:", elapsed_time, "seconds")        
    ds  = xr.Dataset( {  "rttov14": ( ("observation", "chan"), rttov_asfiles, {'units':'K'})  }   )
    ds.to_netcdf('/home/vito.galligani/Work/RTTOV/EXPncfiles/'+dataoutputfile+str(x)+'.nc', 'w')
    
    print(f'Worker {i} /11 finished')
    

    return


#-----------------------------------------------------------------------------
def leer_parallel(): 
        
    nchan      = 22
    nzenith    = 12
    dummy      = np.zeros((nzenith,nchan))          
    specList1  = [str(i) for i in np.arange(1,21998,1)]   # here the last profile is 22000
    angles     = np.arange(0,55,5)


    # Create process objects
    processes = []
    
    for i in range(len(angles)):
        # p2 = multiprocessing.Process(target=CombineRTTOV_fast_clear2b, args=(i, ))
        # processes.append(p2)
        # p2.start()

        # p1 = multiprocessing.Process(target=CombineRTTOV_fast_clear1, args=(i, ))
        # processes.append(p1)
        # p1.start()

        # p = multiprocessing.Process(target=CombineRTTOV_fast_clear, args=(i, ))
        # processes.append(p)
        # p.start()        

        p1 = multiprocessing.Process(target=CombineRTTOV_fast_scatt13, args=(i, ))
        processes.append(p1)
        p1.start()

        p2 = multiprocessing.Process(target=CombineRTTOV_fast_clear13, args=(i, ))
        processes.append(p2)
        p2.start()    
        
    # Wait for all processes to finish
    for p in processes:
        p.join()
        p2.join()

    print('All workers finished')
    
    return 
    
# Function to extract the numeric part of the file name
def extract_number(filename):
    match = re.search(r'_(\d+)\.nc', filename)
    return int(match.group(1)) if match else 0
    
    
if __name__ == '__main__':

    path = '/home/vito.galligani/Work/RTTOV/EXPncfiles'
    
    print('1/6')
    # rttov-13 clear sky simulations: 
    files_rttov13_cs = sorted([f"{path}/rttov13_cs_{x}.nc" for x in range(0, 55, 5) ], key=extract_number)[::-1]
    rttov13_cs = np.empty( (11, 21997, 22))     
    for i in range(len(files_rttov13_cs)):  
        ds = xr.open_dataset(files_rttov13_cs[i])
        rttov13_cs[i,:,:] = ds.rttov14.data[:]

    print('2/6')        
    # rttov-13 all sky simulations:         
    files_rttov13_as = sorted([f"{path}/rttov13_as_{x}.nc" for x in range(0, 55, 5) ], key=extract_number)[::-1]
    rttov13_as = np.empty( (11, 21997, 22)) 
    for i in range(len(files_rttov13_as)):  
        ds = xr.open_dataset(files_rttov13_as[i])
        rttov13_as[i,:,:] = ds.rttov14.data[:]    
    
    print('3/6')
    # rttov-14 all sky simulations:         
    files_rttov14_as = sorted([f"{path}/rttov14_as_{x}.nc" for x in range(0, 55, 5) ], key=extract_number)[::-1]
    rttov14_as = np.empty( (11, 21997, 22)) 
    for i in range(len(files_rttov14_as)):  
        ds = xr.open_dataset(files_rttov14_as[i])
        rttov14_as[i,:,:] = ds.rttov14.data[:]  

    print('4/6')
    # rttov-14 clear sky simulations:         
    files_rttov14_cs = sorted([f"{path}/rttov14_cs_{x}.nc" for x in range(0, 55, 5) ], key=extract_number)[::-1]
    rttov14_cs = np.empty( (11, 21997, 22)) 
    for i in range(len(files_rttov14_cs)):  
        ds = xr.open_dataset(files_rttov14_cs[i])
        rttov14_cs[i,:,:] = ds.rttov14.data[:]       
    
    print('5/6')
    # rttov-14 clear sky simulations (option1 coeffs):         
    files_rttov14_cs = sorted([f"{path}/rttov14_cs_option1{x}.nc" for x in range(0, 55, 5) ], key=extract_number)[::-1]
    rttov14_cs_opt1 = np.empty( (11, 21997, 22)) 
    for i in range(len(files_rttov14_cs)):  
        ds = xr.open_dataset(files_rttov14_cs[i])
        rttov14_cs_opt1[i,:,:] = ds.rttov14.data[:]       
            
    print('6/6')
    # rttov-14 clear sky simulations (option2b coeffs):         
    files_rttov14_cs = sorted([f"{path}/rttov14_cs_option2b{x}.nc" for x in range(0, 55, 5) ], key=extract_number)[::-1]
    rttov14_cs_opt2b = np.empty( (11, 21997, 22)) 
    for i in range(len(files_rttov14_cs)):  
        ds = xr.open_dataset(files_rttov14_cs[i])
        rttov14_cs_opt2b[i,:,:] = ds.rttov14.data[:]        
    
    print('saving')
    # And create a netcdf file
    ds = xr.Dataset( {
                    "rttov13_cs":       ( ( "theta", "observation", "chan"), rttov13_cs, {'units':'K'}),
                    "rttov13_as":       ( ( "theta", "observation", "chan"), rttov13_as, {'units':'K'}),
                    "rttov14_as":       ( ( "theta", "observation", "chan"), rttov14_as, {'units':'K'}),
                    "rttov14_cs":       ( ( "theta", "observation", "chan"), rttov14_cs, {'units':'K'}),
                    "rttov14_cs_opt1":  ( ( "theta", "observation", "chan"), rttov14_cs_opt1, {'units':'K'}),
                    "rttov14_cs_opt2b": ( ( "theta", "observation", "chan"), rttov14_cs_opt2b, {'units':'K'})
                    }   )
    ds.to_netcdf('/home/vito.galligani/Work/RTTOV/EXPncfiles/full_experiments.nc', 'w')

















