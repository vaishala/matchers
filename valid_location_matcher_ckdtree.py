#! /home/mas001/msh/python_environments/base_env/latest/bin/python

'''

This script take the input of the pairs of MODIS and GEM-MACH files that match as far as date and time
This script finds all the indices that have lat/lon coincidences under 0.15 degrees 
This script then indexes into the lat,lon,aod,mass of these geo coincidences
Then the script indexes into a subset of above where the AOD from MODIS is valid and not a fill value
The script finds the indices of the AOD where the values are valid
The script then saves the AOD+mass from gemmach and AOD from modis of those locations and their respective lats and lons and the PM2.5 estimate
The script then exports the list of values in a long pandas dataframe to an nc file for NAPS comparison
'''
import numpy as np
import netCDF4 as nc
import pyhdf
from pyhdf.SD import SD, SDC
from scipy.spatial import cKDTree
import pandas as pd
import matplotlib.pyplot as plt
import xarray as xr
from pathlib import Path  

#Get lat lon boundaries from gem-mach data
f1 = "/home/vat000/msh/model/gem-mach/earthcare/2014-15/June_July_2015/nc/000-000nc/2015060100-000-000_000.nc"
ds1 = nc.Dataset(f1) #gem-mach
#Lat Lon 
lat_gem=ds1['lat_1'][:] #gemmach lat in Deg North
lon_gem=ds1['lon_1'][:]-360 #gemmach lon originally in Degrees East, but changes to negative since we are interested in western hemisphere

pairs_csv_filepath='/home/vat000/vaishala_test/data_and_plots/MODIS_3km_GM_Stephen_15km_File_Pairs_2015/file_matched_datetime.csv'
df_pairs=pd.read_csv(pairs_csv_filepath)
print(np.shape(df_pairs))
df_pairs=df_pairs.iloc[:,1:] #Drops the fist column of indices that is unecessary
df_gemmach_modis=pd.DataFrame({})
print(np.shape(lat_gem))

for i in range(np.shape(df_pairs)[0]):
#Index into first pair of gem-mach and modis matched by date and time in the list of 9000-something
    pair=df_pairs.iloc[i,:]

    gem_filepath=pair['Gemmach_FilePath']
    modis_filepath=pair['MODIS_FilePath']

#Index into lat lon of modis to check location coincidences
    hdf=SD(modis_filepath,SDC.READ)
    lat_modis=hdf.select('Latitude')[:,:]
    lon_modis=hdf.select('Longitude')[:,:]
    print(np.shape(lon_modis))

#Index into other relevant values of modis and gemmach by using an array (to prep toward indexing and PM2.5 product)
    aod1=np.multiply(hdf.select('Corrected_Optical_Depth_Land')[1,:,:],0.001) #inedxing into all lat and lon, applying scaling factor 
    aod1=aod1.flatten()
    lat_modis=lat_modis.flatten()
    lon_modis=lon_modis.flatten()

    aod2=(nc.Dataset(gem_filepath))['AODC'][0,:,:] #time 0
    aod2=aod2.flatten()
    mass2=(nc.Dataset(gem_filepath))['AF'][0,-1,:,:] #time 0 and lowest level (altitude is negative z direction) in ug/m3
    mass2=mass2.flatten()

    mass_secondlast=(nc.Dataset(gem_filepath))['AF'][0,-2,:,:] #AF Concentration of fine particulate (Diam. < 2.5 Microns) µg/m³
                                                               #time 0 and second lowest level (altitude is negative z direction)
    mass_secondlast=mass_secondlast.flatten()

    rh2=(nc.Dataset(gem_filepath))['HR'][0,-1,11:-11,11:-11] #surface relative humidity as a fraction
    #last one is lowest because level3 is going down
    rh2=rh2.flatten()

    bl2=(nc.Dataset(gem_filepath))['H'][0,9:-9,9:-9] #boundary layer height in meters (rlat3 rlon3)

    bl2=bl2.flatten()
    # time = UNLIMITED ; // (1 currently)
    # level1 = 23 ;
    # rlat1 = 500 ;
    # rlon1 = 524 ;
    # level2 = 46 ;
    # rlat2 = 522 ;
    # rlon2 = 546 ;
    # rlat3 = 518 ;
    # rlon3 = 542 ;
    # level3 = 23 ;

    lat_gem=lat_gem.flatten()
    lon_gem=lon_gem.flatten()

    lat_lon_data1 = np.column_stack((lat_modis,lon_modis)) #modis (more points)
    lat_lon_data2 = np.column_stack((lat_gem,lon_gem )) #gemmach (fewer points)
    print(np.shape(lat_lon_data1))
    print(np.shape(lat_lon_data2))
    
    tree_2 = cKDTree(lat_lon_data2)
    closest_distances, closest_indices_2= tree_2.query(lat_lon_data1,k=1,distance_upper_bound=0.15) #gives us indices of items in tree that is closest in lat/lon to data1

    valid_indices_1= np.where(closest_indices_2 < lat_gem.shape[0])[0]
    valid_indices_2 = closest_indices_2[closest_indices_2 < lat_gem.shape[0]]  #Index into valid indices only
    
    array_gemmach_variables=np.stack((aod2,mass2,lat_gem,lon_gem,rh2,bl2,mass_secondlast))
    array_modis_variables=np.stack((aod1,lat_modis,lon_modis))

    matched_array_modis_variables=array_modis_variables[:,valid_indices_1] #index into the gm indices

    if np.shape(matched_array_modis_variables[1])!=(0,): #and np.shape(matched_array_modis_variables[1])[1]>=20:
        matched_array_gemmach_variables=array_gemmach_variables[:,valid_indices_2]
        valid_indices_level2=np.where(matched_array_modis_variables[0] > -9) #Where MODIS AOD is valid
        
        if np.shape(valid_indices_level2)!=(1, 0):
            valid_matched_array_modis_variables=(matched_array_modis_variables[:,valid_indices_level2])[:,0]
            valid_matched_array_gemmach_variables=(matched_array_gemmach_variables[:,valid_indices_level2])[:,0]
            gemmach_ratio=np.divide(valid_matched_array_gemmach_variables[1],valid_matched_array_gemmach_variables[0])
            pm_25_estimate=np.multiply(gemmach_ratio,valid_matched_array_modis_variables[0])

            df_pair=pd.DataFrame({"aod_modis":valid_matched_array_modis_variables[0],"aod_gemmach":valid_matched_array_gemmach_variables[0],\
                "mass_gemmach":valid_matched_array_gemmach_variables[1],"relative_humidity":valid_matched_array_gemmach_variables[4],\
                "boundary_layer_height":valid_matched_array_gemmach_variables[5],\
                "lat_gemmach":valid_matched_array_gemmach_variables[2],\
                "lon_gemmach":valid_matched_array_gemmach_variables[3],"lat_modis":valid_matched_array_modis_variables[1],\
                "lon_modis":valid_matched_array_modis_variables[2],\
                "gemmach_ratio":gemmach_ratio,"pm25estimate":pm_25_estimate,"mass_gemmach_secondlowest":valid_matched_array_gemmach_variables[6]})
            df_pair["gemmach_filepath"]=gem_filepath
            df_pair["modis_filepath"]=modis_filepath    
            df_gemmach_modis=pd.concat([df_gemmach_modis,df_pair],ignore_index=True)
            print("Index of pair",i,"/",np.shape(df_pairs)[0])
            print(valid_matched_array_gemmach_variables[1],valid_matched_array_gemmach_variables[-1])
           
ds_gemmach_modis = xr.Dataset.from_dataframe(df_gemmach_modis)
filepath = Path('/home/vat000/vaishala_test/data_and_plots/MODIS_3km_GM_Stephen_15km_File_Pairs_2015/MODIS_GEMMACH_PM25estimate_pairs.nc')  
ds_gemmach_modis.to_netcdf(filepath)
