#! /home/mas001/msh/python_environments/base_env/latest/bin/python
'''
This script take an input of PM2.5 estimates in an nc file calculated from MODIS and GEMMACH data that is matched to the same 
    hour, day, time, and location to 0.15 deg lat and 0.3 deg lon from /home/vat000/projects/vaishala_test/PM25/nc_reader/MODIS_vs_GEMMACH_2km/2valid_geo_coincidence_file_matcher_exporter.py
This script takes an input of the NAPS 2018 PM2.5 observations in a csv file that has UTC offsets available in the file from /home/vat000/projects/vaishala_test/data_and_plots/naps_data/pm25_pm25_cont_naps_2018_ts_with_OFFSETS.csv

This script determines, for each NAPS station location, the MODIS/GEMMACH pairs that match within 0.3 deg lat lon
The script then checks the datetime of the MODIS/GEMMACH pair
The script then finds the PM2.5 observation from NAPS hourly data at the same datetime of the MODIS/GEMMACH pair (specific to a datetime)
The script then exports the pairs of PM2.5 estimate and PM2.5 surface observation to an nc file

'''
import numpy as np
import netCDF4 as nc
import pandas as pd
from pathlib import Path  
import datetime
import xarray as xr

#GEMMACH MODIS PRODUCT
f1 = '/home/vat000/vaishala_test/data_and_plots/MODIS_3km_GM_Stephen_15km_File_Pairs_2015/MODIS_GEMMACH_PM25estimate_pairs.nc'

ds1 = nc.Dataset(f1)

lat_modis=ds1['lat_modis'][:]
print(np.shape(lat_modis))

lon_modis=ds1['lon_modis'][:]
lat_gemmach=ds1['lat_gemmach'][:]
lon_gemmach=ds1['lon_gemmach'][:]
gemmach_mass_lowest=ds1['mass_gemmach'][:] #only used for calculating vertical gradient
gemmach_mass_2lowest=ds1['mass_gemmach_secondlowest'][:] #only used for calculating vertical gradient
AF_diff=gemmach_mass_lowest-gemmach_mass_2lowest
HR=ds1["relative_humidity"][:]
gemmach_ratio=ds1['gemmach_ratio'][:]
mass_gemmach=gemmach_mass_lowest #The mass value used for the numerator of the gemmach ratio is this lowest mass value
gem_filepath=ds1['gemmach_filepath'][:]
modis_filepaths=ds1['modis_filepath'][:]
aod_modis=ds1['aod_modis'][:] #with scaling already removed
aod_gemmach=ds1['aod_gemmach'][:]
H=ds1["boundary_layer_height"][:]

pm_est=ds1['pm25estimate'][:] #V1; original pm2.5 estimate = GEM-MACH_ratio * MODIS_AOD
# V2_HR_corrected_pm_est=pm_est[:]+0.0067-3.75587*HR[:] #HR corrected PM2.5 estimate
# V3_HR_AF_Corrected_pm_est=V2_HR_corrected_pm_est- 15.4406* AF_diff + 3.2698 #Robust vert gradient corrected PM2.5 estimate
# V4_HR_Ratio_Corected_pm_est=V2_HR_corrected_pm_est- 0.0335* gemmach_ratio + 3.8851  #V4 Robust ratio corrected PM2.5 estimate

print("max gm aod",np.amax(aod_gemmach))
print("max modis aod",np.amax(aod_modis))
print("min gm aod",np.amin(aod_gemmach))
print("min modis aod",np.amin(aod_modis))
print("MODIS")
print("Lat",np.amax(lat_modis))
print(np.amin(lat_modis),"\n")

print("Lon",np.amax(lon_modis))
print(np.amin(lon_modis),"\n")

print("GM")
print("Lat",np.amax(lat_gemmach))
print(np.amin(lat_gemmach))

print("Lon",np.amax(lon_gemmach))
print(np.amin(lon_gemmach))
pm25_estimate=pm_est #Which version of the PM2.5 estimate we are using moving forward
print(np.shape(pm25_estimate))

#NAPS
f = '/home/vat000/vaishala_test/data_and_plots/naps_data/PM25_NAPS_2015_with_OFFSETS.csv'
df=pd.read_csv(f)
df = df.reset_index()
df[pd.to_numeric(df['observation'],errors='coerce').notnull()]
i=0
m=0
n=0
indices_list=[]
indices_list2=[]
offset_list=[]
df_final=pd.DataFrame({})
station_list=[]
oil_sands_naps_stations=['NAPS00070301','NAPS00080111','NAPS00080211','NAPS00080402'\
'NAPS00080502','NAPS00081001','NAPS00082002','NAPS00090120'\
'NAPS00090121','NAPS00090130','NAPS00090132','NAPS00090133'\
'NAPS00090134','NAPS00090135','NAPS00090136','NAPS00090229'\
'NAPS00090230','NAPS00090235','NAPS00090250','NAPS00090302'\
'NAPS00090304','NAPS00090402','NAPS00090502','NAPS00090601'\
'NAPS00090607','NAPS00090608','NAPS00090609','NAPS00090701'\
'NAPS00090702','NAPS00090801','NAPS00090805','NAPS00090806'\
'NAPS00090807','NAPS00090808','NAPS00091101','NAPS00091301'\
'NAPS00091501','NAPS00091801','NAPS00091901','NAPS00092001'\
'NAPS00092201','NAPS00092801','NAPS00092901','NAPS00093001'\
'NAPS00093101','NAPS00093202','NAPS00093901','NAPS00094001'\
'NAPS00094301','NAPS00094401','NAPS00094601','NAPS00100701'\
'NAPS00103202','NAPS00103502','NAPS00104003','NAPS00104101'\
'NAPS00105504','NAPS00106101','NAPS00129601']


for index, row in df.iterrows():
    
    matched=False
    if (index-1)%8762==0:
        print("Station index:",i)
        i+=1
        naps_station_id=row['level_0'] 
        naps_station_id=naps_station_id[-13:-1] #a string to use later finding the pm observation if there is a match
        #print(naps_station_id)
        # if naps_station_id in oil_sands_naps_stations: #To see what we can compare with Mahtab's
            #print("It is in oil sands")

        offset=row['start date'] #always negative in Canada and can go up to UTC -7
        lonn=float(row['level_3'])
        latt=float(row['level_2'])
        abs_lon_diff_array=abs(lon_modis-lonn) 
        if np.amin(abs_lon_diff_array)<0.15: #If the minimum of the diff array is within 0.15 degrees lat lon, then there is a match between the longitudes
            ind=np.where(abs_lon_diff_array<0.15)[0] 

            lat_matched=lat_modis[ind]
            lon_matched=lon_modis[ind]
            pm25_estimate_matched=pm25_estimate[ind]
            mass_gemmach_matched=mass_gemmach[ind]
            gemmach_mass_2lowest_matched=gemmach_mass_2lowest[ind]
            modis_filepaths_matched=modis_filepaths[ind]
            aod_modis_matched=aod_modis[ind]
            aod_gemmach_matched=aod_gemmach[ind]
            lat_gemmach_matched=lat_gemmach[ind]
            lon_gemmach_matched=lon_gemmach[ind]
            gemmach_filepath_matched=gem_filepath[ind]
            HR_matched=HR[ind]
            H_matched=H[ind]
            abs_lat_diff_array=abs(lat_matched-latt)

            if np.amin(abs_lat_diff_array)<0.15: 
                ind=np.where(abs_lat_diff_array<0.15)[0]
                m+=1
                indices_list.append(index)
                #print(ind)

                #if len(ind)>=20:

                for loc_match in ind:
                
                    lat_matched_matched=lat_matched[loc_match]

                    lon_matched_matched=lon_matched[loc_match]
                    lat_gemmach_matched_matched=lat_gemmach_matched[loc_match]
                    lon_gemmach_matched_matched=lon_gemmach_matched[loc_match]
                    pm25_estimate_matched_matched=pm25_estimate_matched[loc_match]
                    mass_gemmach_matched_matched=mass_gemmach_matched[loc_match]
                    aod_modis_matched_matched=aod_modis_matched[loc_match]
                    aod_gemmach_matched_matched=aod_gemmach_matched[loc_match]
                    filepaths_matched_matched=modis_filepaths_matched[loc_match]  #Determine the datetime corresponding to where the geolocation of MODIS occurs from the filename
                    gemmach_filepath_matched_matched=gemmach_filepath_matched[loc_match]
                    HR_matched_matched=HR_matched[loc_match]
                    H_matched_matched=H_matched[loc_match]
                    gemmach_mass_2lowest_matched_matched=gemmach_mass_2lowest_matched[loc_match]
                    indices_list2.append(loc_match)

                    day=(filepaths_matched_matched)[-30:-27]
                    hour=(filepaths_matched_matched)[-26:-24] #2 digits
                    minute=(filepaths_matched_matched)[-24:-22]

                    out=datetime.datetime(2015, 1, 1) + datetime.timedelta(days=int(day) - 1,minutes=int(minute),hours=float(hour))

                    out_updated=out + datetime.timedelta(hours=float(offset))
                    def hour_rounder(t):
                        return (t.replace(minute=0, hour=t.hour)+datetime.timedelta(hours=t.minute//30))
                
                    out_updated_rounded=hour_rounder(out_updated)

                    NAPS_datetime=out_updated_rounded.strftime("%Y-%m-%dT%H:%M:00")

                    if offset[-2:]=='50': 
                        NAPS_offset="-0"+offset[2:3]+":"+"30" 
                    elif offset[-2:]=='00':
                        NAPS_offset="-0"+offset[2:3]+":"+"00"
                    offset_list.append(offset) #The unique offsets in the 2.5km run is -6, -7, -8
                    NAPS_level_1=out_updated_rounded.strftime("%Y-%m-%dT%H:%M:00")+NAPS_offset
                    # print(NAPS_level_1)
                    row_of_matched=df.loc[(df['level_1']==" "+NAPS_level_1) & (df['level_0']==naps_station_id)]

                    pm_observed=list(row_of_matched['level_3'])[0]
                    if pm_observed != " nan":
                        n+=1
                        df_pair=pd.DataFrame({"new_pm25_estimate":pm25_estimate_matched_matched,"pm25_naps":float(pm_observed),\
                        "naps_lat":latt,"naps_lon":lonn,"modis_lat":lat_matched_matched,"modis_lon":lon_matched_matched,\
                        "modis_filename":filepaths_matched_matched,"gemmach_filepath":gemmach_filepath_matched_matched,"station_id":naps_station_id,\
                        'gemmach_mass':mass_gemmach_matched_matched,"gemmach_2nd_lowest_mass":gemmach_mass_2lowest_matched_matched,'gemmach_aod':aod_gemmach_matched_matched,\
                        'modis_aod':aod_modis_matched_matched,"lat_gemmach":lat_gemmach_matched_matched,\
                        "lon_gemmach":lon_gemmach_matched_matched,"relative_humidity":HR_matched_matched,"boundary_layer":H_matched_matched,"naps_datetime":NAPS_level_1},index=[0])
                        df_final=pd.concat([df_final,df_pair],ignore_index=True)
print("Number of matches",n)
print("Number of stations that have a match",m)
print("Unique offsets are:",np.unique(offset_list))

ds = xr.Dataset.from_dataframe(df_final)
filepath = Path('/home/vat000/vaishala_test/data_and_plots/MODIS_3km_GM_Stephen_15km_File_Pairs_2015/matched_pairs_arctic_domain_naps_coincident.nc')  
ds.to_netcdf(filepath)
