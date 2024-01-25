# -*- coding: utf-8 -*-
"""
Created on Wed Nov 22 16:25:22 2023

@author: christophe
"""
import numpy as np
import os
import json
import pathlib
import requests
import pandas as pd
from datetime import datetime
from logging import getLogger
from datetime import datetime
from matplotlib.dates import date2num, num2date
from matplotlib.ticker import Formatter,FuncFormatter, NullLocator, FixedLocator
from matplotlib import pyplot as plt
from matplotlib.pyplot import figure, show,savefig,clabel,contour,setp,gcf,getp,gca,close
from osgeo import gdal
from osgeo.osr import SpatialReference, CoordinateTransformation




logger = getLogger(__name__)

      
def select_stations_from_circle(piezometer_coordinates=[140802.697,456881.771],rmax = 1e3, uurgegevens= False):
 

    if uurgegevens == False:
        weatherstations_data_file = 'weerstations_daggegevens_input.json'
    else:
        weatherstations_data_file = 'weerstations_uurgegevens_input.json'
    
    with open(weatherstations_data_file) as json_file:
        weatherstations_dict = json.load(json_file)
       
    
    keys = list(weatherstations_dict.keys())
    
    weatherstations_structured_array= np.empty(len(keys), dtype=[
        ('name', 'S20'),
        ('lon', 'f4'),
        ('lat', 'f4'),
        ('alt', 'f4'),
        ('STN', 'S10'),
        ('X', 'f4'),
        ('Y', 'f4'),
        ('Z', 'f4'),
        ('dist', 'f4')])
    
    Xc = piezometer_coordinates[0]
    Yc = piezometer_coordinates[1]
    
    for i in range(0,len(keys)):
        key = keys[i]
        STN = weatherstations_dict[key]['STN']
        X = weatherstations_dict[key]['X']
        Y = weatherstations_dict[key]['Y']
        Z = weatherstations_dict[key]['Z']
        lat = weatherstations_dict[key]['lat']
        lon = weatherstations_dict[key]['lon']        
        alt = weatherstations_dict[key]['alt']
    
        dist = np.sqrt((X-Xc)**2+(Y-Yc)**2)    
       
        weatherstations_structured_array[i][0] = key
        weatherstations_structured_array[i][1] = lon
        weatherstations_structured_array[i][2] = lat
        weatherstations_structured_array[i][3] = alt
        weatherstations_structured_array[i][4] = STN
        weatherstations_structured_array[i][5] = X
        weatherstations_structured_array[i][6] = Y
        weatherstations_structured_array[i][7] = Z
        weatherstations_structured_array[i][8] = dist
       
    weatherstations_structured_array = np.sort(weatherstations_structured_array, kind='mergesort', order=['dist'])
    # weatherstations_structured_array = np.flip(weatherstations_structured_array) # to reverse the order
    
    
    list_of_dict = []
    
    encoding = 'utf-8'
    for i in range(0,len(weatherstations_structured_array)):
        key = weatherstations_structured_array[i][0].decode(encoding)
        
        list_of_dict.append(weatherstations_dict[key])
        
    
    
    return list_of_dict







  
def make_url(weatherstation_dict = None, tminstr = '01-01-1980', tmaxstr = '31-10-2023', uurgegevens = False, knmiparams = 'PRCP'):


    """
       zie  https://www.knmi.nl/kennis-en-datacentrum/achtergrond/data-ophalen-vanuit-een-script
    
    The text below comes from the KNMI website, the Dutch Royal Meteorological Institute (in Dutch)
    "vars is een lijst van gewenste variabelen in willekeurige volgorde, aangeduid met hun acroniemen 
    (zoals op de selectiepagina) gescheiden door ':', bijvoorbeeld 'TG:TN:EV24'. 
    Hierin zijn de volgende acroniemen gedefiniëerd om groepen van variabelen aan te duiden:
    
    WIND = DDVEC:FG:FHX:FHX:FX wind
    TEMP = TG:TN:TX:T10N temperatuur
    SUNR = SQ:SP:Q Zonneschijnduur en globale straling
    PRCP = DR:RH:EV24 neerslag en potentiële verdamping
    PRES = PG:PGX:PGN druk op zeeniveau
    VICL = VVN:VVX:NG zicht en bewolking
    MSTR = UG:UX:UN luchtvochtigheid
    ALL alle variabelen
    Default is ALL."
    
    
    
    
    KNMI parameters
    DDVEC: Vectorgemiddelde windrichting in graden (360=noord, 90=oost, 180=zuid, 270=west, 0=windstil/variabel). Meer info
    FHVEC: Vectorgemiddelde windsnelheid (in 0.1 m/s). Meer info
    FG: Etmaalgemiddelde windsnelheid (in 0.1 m/s)
    FHX: Hoogste uurgemiddelde windsnelheid (in 0.1 m/s)
    FHXH: Uurvak waarin FHX is gemeten
    FHN: Laagste uurgemiddelde windsnelheid (in 0.1 m/s)
    FHNH: Uurvak waarin FHN is gemeten
    FXX: Hoogste windstoot (in 0.1 m/s)
    FXXH: Uurvak waarin FXX is gemeten
    TG: Etmaalgemiddelde temperatuur (in 0.1 graden Celsius)
    TN: Minimum temperatuur (in 0.1 graden Celsius)
    TNH: Uurvak waarin TN is gemeten
    TX: Maximum temperatuur (in 0.1 graden Celsius)
    TXH: Uurvak waarin TX is gemeten
    T10N: Minimum temperatuur op 10 cm hoogte (in 0.1 graden Celsius)
    T10NH: 6-uurs tijdvak waarin T10N is gemeten 6=0-6 UT, 12=6-12 UT, 18=12-18 UT, 24=18-24 UT
    SQ: Zonneschijnduur (in 0.1 uur) berekend uit de globale straling (-1 voor <0.05 uur)
    SP: Percentage van de langst mogelijke zonneschijnduur
    Q: Globale straling (in J/cm2)
    DR: Duur van de neerslag (in 0.1 uur)
    RH: Etmaalsom van de neerslag (in 0.1 mm) (-1 voor <0.05 mm)
    RHX: Hoogste uursom van de neerslag (in 0.1 mm) (-1 voor <0.05 mm)
    RHXH: Uurvak waarin RHX is gemeten
    PG: Etmaalgemiddelde luchtdruk herleid tot zeeniveau (in 0.1 hPa) berekend uit 24 uurwaarden
    PX: Hoogste uurwaarde van de luchtdruk herleid tot zeeniveau (in 0.1 hPa)
    PXH: Uurvak waarin PX is gemeten / Hourly division in which PX was measured
    PN: Laagste uurwaarde van de luchtdruk herleid tot zeeniveau (in 0.1 hPa)
    PNH: Uurvak waarin PN is gemeten / Hourly division in which PN was measured
    VVN: Minimum opgetreden zicht; 0: <100 m, 1:100-200 m, 2:200-300 m,..., 49:4900-5000 m, 50:5-6 km, 56:6-7 km, 57:7-8 km,..., 79:29-30 km, 80:30-35 km, 81:35-40 km,..., 89: >70 km
    VVNH: Uurvak waarin VVN is gemeten / Hourly division in which VVN was measured
    VVX: Maximum opgetreden zicht; 0: <100 m, 1:100-200 m, 2:200-300 m,..., 49:4900-5000 m, 50:5-6 km, 56:6-7 km, 57:7-8 km,..., 79:29-30 km, 80:30-35 km, 81:35-40 km,..., 89: >70 km)
    VVXH: Uurvak waarin VVX is gemeten
    NG: Etmaalgemiddelde bewolking (bedekkingsgraad van de bovenlucht in achtsten, 9=bovenlucht onzichtbaar)
    UG: Etmaalgemiddelde relatieve vochtigheid (in procenten)
    UX: Maximale relatieve vochtigheid (in procenten)
    UXH: Uurvak waarin UX is gemeten
    UN: Minimale relatieve vochtigheid (in procenten)
    UNH: Uurvak waarin UN is gemeten
    EV24: Referentiegewasverdamping (Makkink) (in 0.1 mm)
    
        
    Example:
    
    url = " https://www.daggegevens.knmi.nl/klimatologie/daggegevens?stns=235:280:260&vars=VICL:PRCP&start=19700101&end=20090818"
    
    
    """
    
       
    
    if weatherstation_dict is None:
        weatherstation_dict = {'name': 'DeBilt',
          'lon': 5.18,
          'lat': 52.1,
          'alt': 1.9,
          'STN': '260',
          'X': 140802.69756360428,
          'Y': 456881.7714955277,
          'Z': 1.9}
    

    
    url = " https://www.daggegevens.knmi.nl/klimatologie/"
       
    if uurgegevens == False:
        url += 'daggegevens?'
    else:
        url += 'uurgegevens?'   
        
    
    stn = weatherstation_dict['STN']
    url += 'stns='+stn+'&'
    
    knmiparams = knmiparams
    
    url += 'vars='+knmiparams+'&'
        
    dd,mm,yyyy = tminstr.split('-')
    start = yyyy+mm+dd
    url += 'start='+start+'&'
        
        
    dd,mm,yyyy = tmaxstr.split('-')
    end = yyyy+mm+dd    
    url += 'end='+end
    
    
    return url            



def fetchdata(weatherstationname='DeBilt', tminstr = '01-01-1980', tmaxstr = '31-10-2023', uurgegevens = False, 
              knmiparams = 'PRCP',url = None):


    try:
        if uurgegevens == False:
            weatherstations_data_file = 'weerstations_daggegevens_input.json'
        else:
            weatherstations_data_file = 'weerstations_uurgegevens_input.json'
        
        with open(weatherstations_data_file) as json_file:
            weatherstations_dict = json.load(json_file)
        
        weatherstation_dict = weatherstations_dict[weatherstationname]
        
    except:
        logger.warning("Weather station not found, standard station in De Bilt is used")
        weatherstation_dict = {'name': 'DeBilt',
          'lon': 5.18,
          'lat': 52.1,
          'alt': 1.9,
          'STN': '260',
          'X': 140802.70,
          'Y': 456881.77,
          'Z': 1.9}
        
        
    
    path_to_meteo_download_data_folder = get_path_to_meteo_download_data_folder()


    station_code = weatherstation_dict['STN']
    lon = weatherstation_dict['lon']
    lat = weatherstation_dict['lat']
    alt = weatherstation_dict['alt']
    X = weatherstation_dict['X']
    Y = weatherstation_dict['Y']
    

    dd,mm,yyyy = tminstr.split('-')
    start = yyyy+mm+dd
    
    dd,mm,yyyy = tmaxstr.split('-')
    end = yyyy+mm+dd            
        
        

    all_lines = []
    

    if url is None: 
        url = make_url(weatherstation_dict = weatherstation_dict, tminstr = tminstr, tmaxstr = tmaxstr, uurgegevens = uurgegevens, knmiparams = knmiparams)
    else:
        url = url
    
    myobj = {} #keep this for possible later development
    x = requests.post(url, json = myobj)
    #print(x.text[:10000])
    all_lines = x.text.split('\n')
    
    fname = station_code+'_'+weatherstationname+'_'+start+'_'+end+'.txt'
    
    
    if os.path.exists(path_to_meteo_download_data_folder):
        path_to_export_file = os.path.join(path_to_meteo_download_data_folder,fname)
    else:
        curdir = os.getcwd()
        path_to_export_file = curdir+'//'+fname


    F=open(path_to_export_file,'w')
    F.write(x.text)
    F.close()

 

    
    return all_lines

def get_path_to_resources_folder():
    # first check if the download directory for meteorology data exists
    
    #curdir=os.getcwd()
    abs_path = os.path.dirname(os.path.abspath(__file__))
    abs_path_splitted = abs_path.split('\\')
    path_to_parent_folder_elements = abs_path_splitted[:-1]
    path_to_parent_folder = os.path.join(*path_to_parent_folder_elements)
    splitted = path_to_parent_folder.split(':')#trick  to  repair  path
    path_to_parent_folder = splitted[0]+':\\'+splitted[1]#trick  to  repair  path
    path_to_resources_folder = os.path.join(path_to_parent_folder,'resources') 

     
    if not os.path.exists(path_to_resources_folder):
        os.mkdir(path_to_resources_folder)

    return path_to_resources_folder



def get_path_to_meteo_download_data_folder():

    path_to_resources_folder = get_path_to_resources_folder()
    path_to_meteo_download_data_folder = os.path.join(path_to_resources_folder,'downloaded_meteorology') 
    
       
    if not os.path.exists(path_to_meteo_download_data_folder):
        os.mkdir(path_to_meteo_download_data_folder)
        # message = (f'Download directory for meteorology data created at {path_to_meteo_download_data_folder}')
        # print(message) 
        
    return path_to_meteo_download_data_folder
    
    
    

    

def parse_knmi_datafile(datafilename=None, tminstr = '01-01-1980', tmaxstr = '31-10-2023', knmiparams = 'PRCP'):

    t = lambda t_str : date2num(datetime.strptime(t_str,'%d-%m-%Y')) # lamda function to transform time string into time num
    
    if tminstr is not None:
        tmin = t(tminstr)       
    
    if tmaxstr is not None:
        tmax = t(tmaxstr)                

    path_to_meteo_download_data_folder = get_path_to_meteo_download_data_folder()
    

              
    
    datafilepath = os.path.join(path_to_meteo_download_data_folder,datafilename)
    F=open(datafilepath)
    all_lines = F.readlines()
    F.close()  

    
     # Define the Rijksdriehoek projection system (EPSG 28992)
    epsg28992 = SpatialReference()
    epsg28992.ImportFromEPSG(28992)
     
    # correct the towgs84 # chob asks: ok, but where do these 7 parameters come from?
    epsg28992.SetTOWGS84(565.237,50.0087,465.658,-0.406857,0.350733,-1.87035,4.0812)
     
    # Define the wgs84 system (EPSG 4326)
    epsg4326 = SpatialReference()
    epsg4326.ImportFromEPSG(4326)
     
    rd2latlon = CoordinateTransformation(epsg28992, epsg4326)
    latlon2rd = CoordinateTransformation(epsg4326, epsg28992)
    
    
    weatherstation_dict = {}
    
    metadataline = all_lines[6] 
    #metadataline_splitted = metadataline.rstrip().split('\t')
    metadataline_splitted_prior = metadataline.rstrip().split('  ')
    metadataline_splitted = []
    for item in metadataline_splitted_prior:
        if item != '':
            metadataline_splitted.append(item)

    name = metadataline_splitted[4]
    namesplitted = name.split(' ')
    if len(namesplitted) > 1:
        name = namesplitted[0]
        for i in range(1,len(namesplitted)):
            name += namesplitted[i]
    
    lon = float(metadataline_splitted[1])
    lat = float(metadataline_splitted[2])
    alt = float(metadataline_splitted[3])
    STN = metadataline_splitted[0].split(' ')[1]
   
    XYZ = latlon2rd.TransformPoint(lat,lon,alt)
    X = XYZ[0]
    Y = XYZ[1]
    Z = XYZ[2]    

    
    weatherstation_dict['name'] = name
    weatherstation_dict['lon'] = lon
    weatherstation_dict['lat'] = lat
    weatherstation_dict['alt'] = alt   
    weatherstation_dict['STN'] = STN
    weatherstation_dict['X'] = X
    weatherstation_dict['Y'] = Y  
    weatherstation_dict['Z'] = Z 
    
    
    dd,mm,yyyy = tminstr.split('-')
    start = yyyy+mm+dd
    
    dd,mm,yyyy = tmaxstr.split('-')
    end = yyyy+mm+dd      
    
        

    

                
    if len(all_lines) > 0:    
        variables = {}    
        
        variables['DDVEC'] = 'Vectorgemiddelde windrichting in graden (360=noord, 90=oost, 180=zuid, 270=west, 0=windstil/variabel). Meer info'
        variables['FHVEC'] = 'Vectorgemiddelde windsnelheid (in 0.1 m/s). Meer info'
        variables['FG'] = 'Etmaalgemiddelde windsnelheid (in 0.1 m/s)'
        variables['FHX'] = 'Hoogste uurgemiddelde windsnelheid (in 0.1 m/s)'
        variables['FHXH'] = 'Uurvak waarin FHX is gemeten'
        variables['FHN'] = 'Laagste uurgemiddelde windsnelheid (in 0.1 m/s)'
        variables['FHNH'] = 'Uurvak waarin FHN is gemeten'
        variables['FXX'] = 'Hoogste windstoot (in 0.1 m/s)'
        variables['FXXH'] = 'Uurvak waarin FXX is gemeten'
        variables['TG'] = 'Etmaalgemiddelde temperatuur (in 0.1 graden Celsius)'
        variables['TN'] = 'Minimum temperatuur (in 0.1 graden Celsius)'
        variables['TNH'] = 'Uurvak waarin TN is gemeten'
        variables['TX'] = 'Maximum temperatuur (in 0.1 graden Celsius)'
        variables['TXH'] = 'Uurvak waarin TX is gemeten'
        variables['T10N'] = 'Minimum temperatuur op 10 cm hoogte (in 0.1 graden Celsius)'
        variables['T10NH'] = '6-uurs tijdvak waarin T10N is gemeten 6=0-6 UT, 12=6-12 UT, 18=12-18 UT, 24=18-24 UT'
        variables['SQ'] = 'Zonneschijnduur (in 0.1 uur) berekend uit de globale straling (-1 voor <0.05 uur)'
        variables['SP'] = 'Percentage van de langst mogelijke zonneschijnduur'
        variables['Q'] = 'Globale straling (in J/cm2)'
        variables['DR'] = 'Duur van de neerslag (in 0.1 uur)'
        variables['RH'] = 'Etmaalsom van de neerslag (in 0.1 mm) (-1 voor <0.05 mm)'
        variables['RHX'] = 'Hoogste uursom van de neerslag (in 0.1 mm) (-1 voor <0.05 mm)'
        variables['RHXH'] = 'Uurvak waarin RHX is gemeten'
        variables['PG'] = 'Etmaalgemiddelde luchtdruk herleid tot zeeniveau (in 0.1 hPa) berekend uit 24 uurwaarden'
        variables['PX'] = 'Hoogste uurwaarde van de luchtdruk herleid tot zeeniveau (in 0.1 hPa)'
        variables['PXH'] = 'Uurvak waarin PX is gemeten / Hourly division in which PX was measured'
        variables['PN'] = 'Laagste uurwaarde van de luchtdruk herleid tot zeeniveau (in 0.1 hPa)'
        variables['PNH'] = 'Uurvak waarin PN is gemeten / Hourly division in which PN was measured'
        variables['VVN'] = 'Minimum opgetreden zicht; 0 :<100 m, 1:100-200 m, 2:200-300 m,..., 49:4900-5000 m, 50:5-6 km, 56:6-7 km, 57:7-8 km,..., 79:29-30 km, 80:30-35 km, 81:35-40 km,..., 89 :>70 km'
        variables['VVNH'] = 'Uurvak waarin VVN is gemeten / Hourly division in which VVN was measured'
        variables['VVX'] = 'Maximum opgetreden zicht; 0: <100 m, 1:100-200 m, 2:200-300 m,..., 49:4900-5000 m, 50:5-6 km, 56:6-7 km, 57:7-8 km,..., 79:29-30 km, 80:30-35 km, 81:35-40 km,..., 89: >70 km)'
        variables['VVXH'] = 'Uurvak waarin VVX is gemeten'
        variables['NG'] = 'Etmaalgemiddelde bewolking (bedekkingsgraad van de bovenlucht in achtsten, 9=bovenlucht onzichtbaar)'
        variables['UG'] = 'Etmaalgemiddelde relatieve vochtigheid (in procenten)'
        variables['UX'] = 'Maximale relatieve vochtigheid (in procenten)'
        variables['UXH'] = 'Uurvak waarin UX is gemeten'
        variables['UN'] = 'Minimale relatieve vochtigheid (in procenten)'
        variables['UNH'] = 'Uurvak waarin UN is gemeten'
        variables['EV24'] = 'Referentiegewasverdamping (Makkink) (in 0.1 mm)'    
        
        
        headerline = all_lines[10]
        headerline_items = headerline.rstrip().split(',')
        
        datalines = all_lines[11:]
        
        columns = {}
        for key in variables:
            for j in range(0,len(headerline_items)):
                if key == headerline_items[j].strip():
                    columns[key]  = j
            
           
        weatherstation_dict['tseries'] = {}
        for key in columns:
            j = columns[key]
            
        
            #First, read the data as given in knmi file
            records_list = []
            for i in range(0,len(datalines)):
                dataline = datalines[i]
                datalinelist = dataline.rstrip().split(',')
                datestring = datalinelist[1]
                yyyy = datestring[:4]
                mm = datestring[4:6]
                dd = datestring[6:8]
                datestringnew = dd+'-'+mm+'-'+yyyy
                dt = datetime.strptime(datestringnew,'%d-%m-%Y') #date numbers begin at mednight
                timenum = date2num(dt)+0.99 # add 0.99 to point to the end of the day (so we have the cumulated precipitation at the end of the day)
                value = float(datalinelist[j])
                    
              
                if key in ['RH','EV24']:
                    if value == -1:
                        value = 0 # because per convention in KNMO data, if precipitation <0.05mm per day is attributed a value of -1
                    value /= 10000. # convert from 0.1mm per day to meter per day
        
                records_list.append([timenum,value])
                
            tseries = np.empty((len(records_list),2))
            for i in range(0,len(records_list)):
                tseries[i,0] = records_list[i][0]
                tseries[i,1] = records_list[i][1]
                
            weatherstation_dict['tseries'][key] = tseries
        
            if key == 'RH':
                parametername = 'PREC'
            elif key == 'EV24':
                parametername = 'EVAP'
            else:
                parametername = key
                
            fname = STN
            fname +='_'+name
            fname +='_'+parametername
            fname +='_'+start
            fname +='_'+end
            fname +='.csv'
            
            
            header1 = 'Meteorological_data_type,time_serie_name,lon,lat,altitude,station_code,X,Y,comment'
            
            
            header2 = parametername
            header2 +=','+fname.split('.')[0]
            header2 +=','+str('%4.2f' %lon)
            header2 +=','+str('%4.2f' %lat)
            header2 +=','+str('%4.2f' %alt)
            header2 +=','+STN
            header2 +=','+str('%10.2f' %X)
            header2 +=','+str('%10.2f' %Y)
            if key in ['RH','EV24']:
                comment = 'Units: meters - Note that the original KNMI parameter is '+variables[key]
            else:
                comment = variables[key]
                
            header2 +=','+comment  
            
            path_to_resources_folder = get_path_to_resources_folder()
            exportfilepath = os.path.join(path_to_resources_folder,fname)
            
            F=open(exportfilepath,'w')
            F.write(header1+'\n')
            F.write(header2+'\n')
            F.write('\n')
            F.write('Date(dd-mm-yyyy HH:MM),value (unit: see comment)\n')
            
            for i in range(0,len(tseries)):
                t = num2date(tseries[i,0])
                value = tseries[i,1]
                record=str(t.day)+'-'+str(t.month)+'-'+str(t.year)+'  '+str(t.hour).zfill(2)+':'+str(t.minute).zfill(2)+','+str('%10.6f' %value)
                F.write(record+'\n')
            F.close()    

            weatherstation_dict['tseries'][key] = tseries
    
    
    
    return all_lines








weatherstations_dict_dag = {
    "IJmond": {
        "name": "IJmond",
        "lon": 4.518,
        "lat": 52.465,
        "alt": 0.0,
        "STN": "209",
        "X": 95931.093433345,
        "Y": 497826.96719725226,
        "Z": 0.0
    },
    "ValkenburgZh": {
        "name": "ValkenburgZh",
        "lon": 4.43,
        "lat": 52.171,
        "alt": -0.2,
        "STN": "210",
        "X": 89519.11346407089,
        "Y": 465192.84233140404,
        "Z": -0.2
    },
    "Voorschoten": {
        "name": "Voorschoten",
        "lon": 4.437,
        "lat": 52.141,
        "alt": -1.1,
        "STN": "215",
        "X": 89954.24818105088,
        "Y": 461848.9366045671,
        "Z": -1.1
    },
    "IJmuiden": {
        "name": "IJmuiden",
        "lon": 4.555,
        "lat": 52.463,
        "alt": 4.4,
        "STN": "225",
        "X": 98442.88140817755,
        "Y": 497574.90517506644,
        "Z": 4.4
    },
    "DeKooy": {
        "name": "DeKooy",
        "lon": 4.781,
        "lat": 52.928,
        "alt": 1.2,
        "STN": "235",
        "X": 114235.45147115004,
        "Y": 549163.1251914556,
        "Z": 1.2
    },
    "Schiphol": {
        "name": "Schiphol",
        "lon": 4.79,
        "lat": 52.318,
        "alt": -3.3,
        "STN": "240",
        "X": 114280.3459602095,
        "Y": 481284.2071902523,
        "Z": -3.3
    },
    "Vlieland": {
        "name": "Vlieland",
        "lon": 4.921,
        "lat": 53.241,
        "alt": 10.8,
        "STN": "242",
        "X": 123874.83750250352,
        "Y": 583924.9767709821,
        "Z": 10.8
    },
    "Wijdenes": {
        "name": "Wijdenes",
        "lon": 5.174,
        "lat": 52.634,
        "alt": 0.8,
        "STN": "248",
        "X": 140566.29430495645,
        "Y": 516298.26080261456,
        "Z": 0.8
    },
    "Berkhout": {
        "name": "Berkhout",
        "lon": 4.979,
        "lat": 52.644,
        "alt": -2.4,
        "STN": "249",
        "X": 127371.41512666266,
        "Y": 517467.6485436297,
        "Z": -2.4
    },
    "HoornTerschelling": {
        "name": "HoornTerschelling",
        "lon": 5.346,
        "lat": 53.392,
        "alt": 0.7,
        "STN": "251",
        "X": 152258.4461636239,
        "Y": 600630.4349474516,
        "Z": 0.7
    },
    "WijkaanZee": {
        "name": "WijkaanZee",
        "lon": 4.603,
        "lat": 52.506,
        "alt": 8.5,
        "STN": "257",
        "X": 101756.71505533597,
        "Y": 502322.80483447673,
        "Z": 8.5
    },
    "Houtribdijk": {
        "name": "Houtribdijk",
        "lon": 5.401,
        "lat": 52.649,
        "alt": 7.3,
        "STN": "258",
        "X": 155933.56143313224,
        "Y": 517946.1824305308,
        "Z": 7.3
    },
    "DeBilt": {
        "name": "DeBilt",
        "lon": 5.18,
        "lat": 52.1,
        "alt": 1.9,
        "STN": "260",
        "X": 140802.69756360428,
        "Y": 456881.7714955277,
        "Z": 1.9
    },
    "Soesterberg": {
        "name": "Soesterberg",
        "lon": 5.274,
        "lat": 52.13,
        "alt": 13.9,
        "STN": "265",
        "X": 147248.64963213712,
        "Y": 460205.3553571091,
        "Z": 13.9
    },
    "Stavoren": {
        "name": "Stavoren",
        "lon": 5.384,
        "lat": 52.898,
        "alt": -1.3,
        "STN": "267",
        "X": 154784.24164376254,
        "Y": 545653.648154086,
        "Z": -1.3
    },
    "Lelystad": {
        "name": "Lelystad",
        "lon": 5.52,
        "lat": 52.458,
        "alt": -3.7,
        "STN": "269",
        "X": 164026.00734818866,
        "Y": 496701.96147179976,
        "Z": -3.7
    },
    "Leeuwarden": {
        "name": "Leeuwarden",
        "lon": 5.752,
        "lat": 53.224,
        "alt": 1.2,
        "STN": "270",
        "X": 179364.06632110506,
        "Y": 581994.1253899828,
        "Z": 1.2
    },
    "Marknesse": {
        "name": "Marknesse",
        "lon": 5.888,
        "lat": 52.703,
        "alt": -3.3,
        "STN": "273",
        "X": 188849.76174235396,
        "Y": 524072.1381718632,
        "Z": -3.3
    },
    "Deelen": {
        "name": "Deelen",
        "lon": 5.873,
        "lat": 52.056,
        "alt": 48.2,
        "STN": "275",
        "X": 188318.79075868436,
        "Y": 452077.68963571027,
        "Z": 48.2
    },
    "Lauwersoog": {
        "name": "Lauwersoog",
        "lon": 6.2,
        "lat": 53.413,
        "alt": 2.9,
        "STN": "277",
        "X": 209047.27058330286,
        "Y": 603272.335235052,
        "Z": 2.9
    },
    "Heino": {
        "name": "Heino",
        "lon": 6.259,
        "lat": 52.435,
        "alt": 3.6,
        "STN": "278",
        "X": 214285.15091157693,
        "Y": 494491.5238856772,
        "Z": 3.6
    },
    "Hoogeveen": {
        "name": "Hoogeveen",
        "lon": 6.574,
        "lat": 52.75,
        "alt": 15.8,
        "STN": "279",
        "X": 235129.9864781872,
        "Y": 529842.8776996406,
        "Z": 15.8
    },
    "Eelde": {
        "name": "Eelde",
        "lon": 6.585,
        "lat": 53.125,
        "alt": 5.2,
        "STN": "280",
        "X": 235179.97007018456,
        "Y": 571581.2231981044,
        "Z": 5.2
    },
    "Hupsel": {
        "name": "Hupsel",
        "lon": 6.657,
        "lat": 52.069,
        "alt": 29.1,
        "STN": "283",
        "X": 242062.1564910458,
        "Y": 454174.13306507206,
        "Z": 29.1
    },
    "Huibertgat": {
        "name": "Huibertgat",
        "lon": 6.399,
        "lat": 53.575,
        "alt": 0.0,
        "STN": "285",
        "X": 222025.617372105,
        "Y": 621469.6201671187,
        "Z": 0.0
    },
    "NieuwBeerta": {
        "name": "NieuwBeerta",
        "lon": 7.15,
        "lat": 53.196,
        "alt": -0.2,
        "STN": "286",
        "X": 272803.15050124796,
        "Y": 580257.8167607762,
        "Z": -0.2
    },
    "Twenthe": {
        "name": "Twenthe",
        "lon": 6.891,
        "lat": 52.274,
        "alt": 34.8,
        "STN": "290",
        "X": 257631.77844409097,
        "Y": 477285.68626925506,
        "Z": 34.8
    },
    "Cadzand": {
        "name": "Cadzand",
        "lon": 3.379,
        "lat": 51.381,
        "alt": 0.0,
        "STN": "308",
        "X": 15204.913169225387,
        "Y": 378794.06751899014,
        "Z": 0.0
    },
    "Vlissingen": {
        "name": "Vlissingen",
        "lon": 3.596,
        "lat": 51.442,
        "alt": 8.0,
        "STN": "310",
        "X": 30475.202405982665,
        "Y": 385185.512857629,
        "Z": 8.0
    },
    "Hoofdplaat": {
        "name": "Hoofdplaat",
        "lon": 3.672,
        "lat": 51.379,
        "alt": 0.0,
        "STN": "311",
        "X": 35593.21776889522,
        "Y": 378050.99033211346,
        "Z": 0.0
    },
    "Oosterschelde": {
        "name": "Oosterschelde",
        "lon": 3.622,
        "lat": 51.768,
        "alt": 0.0,
        "STN": "312",
        "X": 33161.506387629226,
        "Y": 421402.3911734065,
        "Z": 0.0
    },
    "VlaktevanDeRaan": {
        "name": "VlaktevanDeRaan",
        "lon": 3.242,
        "lat": 51.505,
        "alt": 0.0,
        "STN": "313",
        "X": 6075.775561416813,
        "Y": 392856.63421542535,
        "Z": 0.0
    },
    "Hansweert": {
        "name": "Hansweert",
        "lon": 3.998,
        "lat": 51.447,
        "alt": 0.0,
        "STN": "315",
        "X": 58430.28524182568,
        "Y": 385132.1483284156,
        "Z": 0.0
    },
    "Schaar": {
        "name": "Schaar",
        "lon": 3.694,
        "lat": 51.657,
        "alt": 0.0,
        "STN": "316",
        "X": 37843.182266523465,
        "Y": 408937.42539045145,
        "Z": 0.0
    },
    "Westdorpe": {
        "name": "Westdorpe",
        "lon": 3.861,
        "lat": 51.226,
        "alt": 1.7,
        "STN": "319",
        "X": 48393.44633539117,
        "Y": 360739.8987551535,
        "Z": 1.7
    },
    "Wilhelminadorp": {
        "name": "Wilhelminadorp",
        "lon": 3.884,
        "lat": 51.527,
        "alt": 1.4,
        "STN": "323",
        "X": 50689.90913567814,
        "Y": 394188.1480315051,
        "Z": 1.4
    },
    "Stavenisse": {
        "name": "Stavenisse",
        "lon": 4.006,
        "lat": 51.596,
        "alt": 0.0,
        "STN": "324",
        "X": 59300.52685939122,
        "Y": 401696.0347009504,
        "Z": 0.0
    },
    "HoekvanHolland": {
        "name": "HoekvanHolland",
        "lon": 4.122,
        "lat": 51.992,
        "alt": 11.9,
        "STN": "330",
        "X": 68103.53034220003,
        "Y": 445602.26883290964,
        "Z": 11.9
    },
    "Tholen": {
        "name": "Tholen",
        "lon": 4.193,
        "lat": 51.48,
        "alt": 0.0,
        "STN": "331",
        "X": 72044.85907365358,
        "Y": 388562.9503640573,
        "Z": 0.0
    },
    "Woensdrecht": {
        "name": "Woensdrecht",
        "lon": 4.342,
        "lat": 51.449,
        "alt": 19.2,
        "STN": "340",
        "X": 82345.17359185278,
        "Y": 384955.4236763549,
        "Z": 19.2
    },
    "RotterdamGeulhaven": {
        "name": "RotterdamGeulhaven",
        "lon": 4.313,
        "lat": 51.893,
        "alt": 3.5,
        "STN": "343",
        "X": 81058.22539368055,
        "Y": 434377.45888780034,
        "Z": 3.5
    },
    "Rotterdam": {
        "name": "Rotterdam",
        "lon": 4.447,
        "lat": 51.962,
        "alt": -4.3,
        "STN": "344",
        "X": 90380.92096120957,
        "Y": 441925.9109169654,
        "Z": -4.3
    },
    "CabauwMast": {
        "name": "CabauwMast",
        "lon": 4.926,
        "lat": 51.97,
        "alt": -0.7,
        "STN": "348",
        "X": 123307.23288106002,
        "Y": 442498.5025450533,
        "Z": -0.7
    },
    "Gilze-Rijen": {
        "name": "Gilze-Rijen",
        "lon": 4.936,
        "lat": 51.566,
        "alt": 14.9,
        "STN": "350",
        "X": 123715.72551105871,
        "Y": 397548.013806552,
        "Z": 14.9
    },
    "Herwijnen": {
        "name": "Herwijnen",
        "lon": 5.146,
        "lat": 51.859,
        "alt": 0.7,
        "STN": "356",
        "X": 138384.09667411912,
        "Y": 430076.04933944705,
        "Z": 0.7
    },
    "Eindhoven": {
        "name": "Eindhoven",
        "lon": 5.377,
        "lat": 51.451,
        "alt": 22.6,
        "STN": "370",
        "X": 154290.9184795222,
        "Y": 384657.31016766746,
        "Z": 22.6
    },
    "Volkel": {
        "name": "Volkel",
        "lon": 5.707,
        "lat": 51.659,
        "alt": 22.0,
        "STN": "375",
        "X": 177127.99005000386,
        "Y": 407846.336823779,
        "Z": 22.0
    },
    "Ell": {
        "name": "Ell",
        "lon": 5.763,
        "lat": 51.198,
        "alt": 30.0,
        "STN": "377",
        "X": 181267.25444988074,
        "Y": 356578.54363042983,
        "Z": 30.0
    },
    "Maastricht": {
        "name": "Maastricht",
        "lon": 5.762,
        "lat": 50.906,
        "alt": 114.3,
        "STN": "380",
        "X": 181363.98355052367,
        "Y": 324093.5619647858,
        "Z": 114.3
    },
    "Arcen": {
        "name": "Arcen",
        "lon": 6.197,
        "lat": 51.498,
        "alt": 19.5,
        "STN": "391",
        "X": 211231.34266846423,
        "Y": 390198.61134588905,
        "Z": 19.5
    }
}

weatherstations_dict_uur = {
    "IJmond": {
        "name": "IJmond",
        "lon": 4.518,
        "lat": 52.465,
        "alt": 0.0,
        "STN": "209",
        "X": 95931.093433345,
        "Y": 497826.96719725226,
        "Z": 0.0
    },
    "ValkenburgZh": {
        "name": "ValkenburgZh",
        "lon": 4.43,
        "lat": 52.171,
        "alt": -0.2,
        "STN": "210",
        "X": 89519.11346407089,
        "Y": 465192.84233140404,
        "Z": -0.2
    },
    "Voorschoten": {
        "name": "Voorschoten",
        "lon": 4.437,
        "lat": 52.141,
        "alt": -1.1,
        "STN": "215",
        "X": 89954.24818105088,
        "Y": 461848.9366045671,
        "Z": -1.1
    },
    "IJmuiden": {
        "name": "IJmuiden",
        "lon": 4.555,
        "lat": 52.463,
        "alt": 4.4,
        "STN": "225",
        "X": 98442.88140817755,
        "Y": 497574.90517506644,
        "Z": 4.4
    },
    "DeKooy": {
        "name": "DeKooy",
        "lon": 4.781,
        "lat": 52.928,
        "alt": 1.2,
        "STN": "235",
        "X": 114235.45147115004,
        "Y": 549163.1251914556,
        "Z": 1.2
    },
    "Schiphol": {
        "name": "Schiphol",
        "lon": 4.79,
        "lat": 52.318,
        "alt": -3.3,
        "STN": "240",
        "X": 114280.3459602095,
        "Y": 481284.2071902523,
        "Z": -3.3
    },
    "Vlieland": {
        "name": "Vlieland",
        "lon": 4.921,
        "lat": 53.241,
        "alt": 10.8,
        "STN": "242",
        "X": 123874.83750250352,
        "Y": 583924.9767709821,
        "Z": 10.8
    },
    "Wijdenes": {
        "name": "Wijdenes",
        "lon": 5.174,
        "lat": 52.634,
        "alt": 0.8,
        "STN": "248",
        "X": 140566.29430495645,
        "Y": 516298.26080261456,
        "Z": 0.8
    },
    "Berkhout": {
        "name": "Berkhout",
        "lon": 4.979,
        "lat": 52.644,
        "alt": -2.4,
        "STN": "249",
        "X": 127371.41512666266,
        "Y": 517467.6485436297,
        "Z": -2.4
    },
    "HoornTerschelling": {
        "name": "HoornTerschelling",
        "lon": 5.346,
        "lat": 53.392,
        "alt": 0.7,
        "STN": "251",
        "X": 152258.4461636239,
        "Y": 600630.4349474516,
        "Z": 0.7
    },
    "WijkaanZee": {
        "name": "WijkaanZee",
        "lon": 4.603,
        "lat": 52.506,
        "alt": 8.5,
        "STN": "257",
        "X": 101756.71505533597,
        "Y": 502322.80483447673,
        "Z": 8.5
    },
    "Houtribdijk": {
        "name": "Houtribdijk",
        "lon": 5.401,
        "lat": 52.649,
        "alt": 7.3,
        "STN": "258",
        "X": 155933.56143313224,
        "Y": 517946.1824305308,
        "Z": 7.3
    },
    "DeBilt": {
        "name": "DeBilt",
        "lon": 5.18,
        "lat": 52.1,
        "alt": 1.9,
        "STN": "260",
        "X": 140802.69756360428,
        "Y": 456881.7714955277,
        "Z": 1.9
    },
    "Soesterberg": {
        "name": "Soesterberg",
        "lon": 5.274,
        "lat": 52.13,
        "alt": 13.9,
        "STN": "265",
        "X": 147248.64963213712,
        "Y": 460205.3553571091,
        "Z": 13.9
    },
    "Stavoren": {
        "name": "Stavoren",
        "lon": 5.384,
        "lat": 52.898,
        "alt": -1.3,
        "STN": "267",
        "X": 154784.24164376254,
        "Y": 545653.648154086,
        "Z": -1.3
    },
    "Lelystad": {
        "name": "Lelystad",
        "lon": 5.52,
        "lat": 52.458,
        "alt": -3.7,
        "STN": "269",
        "X": 164026.00734818866,
        "Y": 496701.96147179976,
        "Z": -3.7
    },
    "Leeuwarden": {
        "name": "Leeuwarden",
        "lon": 5.752,
        "lat": 53.224,
        "alt": 1.2,
        "STN": "270",
        "X": 179364.06632110506,
        "Y": 581994.1253899828,
        "Z": 1.2
    },
    "Marknesse": {
        "name": "Marknesse",
        "lon": 5.888,
        "lat": 52.703,
        "alt": -3.3,
        "STN": "273",
        "X": 188849.76174235396,
        "Y": 524072.1381718632,
        "Z": -3.3
    },
    "Deelen": {
        "name": "Deelen",
        "lon": 5.873,
        "lat": 52.056,
        "alt": 48.2,
        "STN": "275",
        "X": 188318.79075868436,
        "Y": 452077.68963571027,
        "Z": 48.2
    },
    "Lauwersoog": {
        "name": "Lauwersoog",
        "lon": 6.2,
        "lat": 53.413,
        "alt": 2.9,
        "STN": "277",
        "X": 209047.27058330286,
        "Y": 603272.335235052,
        "Z": 2.9
    },
    "Heino": {
        "name": "Heino",
        "lon": 6.259,
        "lat": 52.435,
        "alt": 3.6,
        "STN": "278",
        "X": 214285.15091157693,
        "Y": 494491.5238856772,
        "Z": 3.6
    },
    "Hoogeveen": {
        "name": "Hoogeveen",
        "lon": 6.574,
        "lat": 52.75,
        "alt": 15.8,
        "STN": "279",
        "X": 235129.9864781872,
        "Y": 529842.8776996406,
        "Z": 15.8
    },
    "Eelde": {
        "name": "Eelde",
        "lon": 6.585,
        "lat": 53.125,
        "alt": 5.2,
        "STN": "280",
        "X": 235179.97007018456,
        "Y": 571581.2231981044,
        "Z": 5.2
    },
    "Hupsel": {
        "name": "Hupsel",
        "lon": 6.657,
        "lat": 52.069,
        "alt": 29.1,
        "STN": "283",
        "X": 242062.1564910458,
        "Y": 454174.13306507206,
        "Z": 29.1
    },
    "Huibertgat": {
        "name": "Huibertgat",
        "lon": 6.399,
        "lat": 53.575,
        "alt": 0.0,
        "STN": "285",
        "X": 222025.617372105,
        "Y": 621469.6201671187,
        "Z": 0.0
    },
    "NieuwBeerta": {
        "name": "NieuwBeerta",
        "lon": 7.15,
        "lat": 53.196,
        "alt": -0.2,
        "STN": "286",
        "X": 272803.15050124796,
        "Y": 580257.8167607762,
        "Z": -0.2
    },
    "Twenthe": {
        "name": "Twenthe",
        "lon": 6.891,
        "lat": 52.274,
        "alt": 34.8,
        "STN": "290",
        "X": 257631.77844409097,
        "Y": 477285.68626925506,
        "Z": 34.8
    },
    "Cadzand": {
        "name": "Cadzand",
        "lon": 3.379,
        "lat": 51.381,
        "alt": 0.0,
        "STN": "308",
        "X": 15204.913169225387,
        "Y": 378794.06751899014,
        "Z": 0.0
    },
    "Vlissingen": {
        "name": "Vlissingen",
        "lon": 3.596,
        "lat": 51.442,
        "alt": 8.0,
        "STN": "310",
        "X": 30475.202405982665,
        "Y": 385185.512857629,
        "Z": 8.0
    },
    "Hoofdplaat": {
        "name": "Hoofdplaat",
        "lon": 3.672,
        "lat": 51.379,
        "alt": 0.0,
        "STN": "311",
        "X": 35593.21776889522,
        "Y": 378050.99033211346,
        "Z": 0.0
    },
    "Oosterschelde": {
        "name": "Oosterschelde",
        "lon": 3.622,
        "lat": 51.768,
        "alt": 0.0,
        "STN": "312",
        "X": 33161.506387629226,
        "Y": 421402.3911734065,
        "Z": 0.0
    },
    "VlaktevanDeRaan": {
        "name": "VlaktevanDeRaan",
        "lon": 3.242,
        "lat": 51.505,
        "alt": 0.0,
        "STN": "313",
        "X": 6075.775561416813,
        "Y": 392856.63421542535,
        "Z": 0.0
    },
    "Hansweert": {
        "name": "Hansweert",
        "lon": 3.998,
        "lat": 51.447,
        "alt": 0.0,
        "STN": "315",
        "X": 58430.28524182568,
        "Y": 385132.1483284156,
        "Z": 0.0
    },
    "Schaar": {
        "name": "Schaar",
        "lon": 3.694,
        "lat": 51.657,
        "alt": 0.0,
        "STN": "316",
        "X": 37843.182266523465,
        "Y": 408937.42539045145,
        "Z": 0.0
    },
    "Westdorpe": {
        "name": "Westdorpe",
        "lon": 3.861,
        "lat": 51.226,
        "alt": 1.7,
        "STN": "319",
        "X": 48393.44633539117,
        "Y": 360739.8987551535,
        "Z": 1.7
    },
    "Wilhelminadorp": {
        "name": "Wilhelminadorp",
        "lon": 3.884,
        "lat": 51.527,
        "alt": 1.4,
        "STN": "323",
        "X": 50689.90913567814,
        "Y": 394188.1480315051,
        "Z": 1.4
    },
    "Stavenisse": {
        "name": "Stavenisse",
        "lon": 4.006,
        "lat": 51.596,
        "alt": 0.0,
        "STN": "324",
        "X": 59300.52685939122,
        "Y": 401696.0347009504,
        "Z": 0.0
    },
    "HoekvanHolland": {
        "name": "HoekvanHolland",
        "lon": 4.122,
        "lat": 51.992,
        "alt": 11.9,
        "STN": "330",
        "X": 68103.53034220003,
        "Y": 445602.26883290964,
        "Z": 11.9
    },
    "Tholen": {
        "name": "Tholen",
        "lon": 4.193,
        "lat": 51.48,
        "alt": 0.0,
        "STN": "331",
        "X": 72044.85907365358,
        "Y": 388562.9503640573,
        "Z": 0.0
    },
    "Woensdrecht": {
        "name": "Woensdrecht",
        "lon": 4.342,
        "lat": 51.449,
        "alt": 19.2,
        "STN": "340",
        "X": 82345.17359185278,
        "Y": 384955.4236763549,
        "Z": 19.2
    },
    "RotterdamGeulhaven": {
        "name": "RotterdamGeulhaven",
        "lon": 4.313,
        "lat": 51.893,
        "alt": 3.5,
        "STN": "343",
        "X": 81058.22539368055,
        "Y": 434377.45888780034,
        "Z": 3.5
    },
    "Rotterdam": {
        "name": "Rotterdam",
        "lon": 4.447,
        "lat": 51.962,
        "alt": -4.3,
        "STN": "344",
        "X": 90380.92096120957,
        "Y": 441925.9109169654,
        "Z": -4.3
    },
    "CabauwMast": {
        "name": "CabauwMast",
        "lon": 4.926,
        "lat": 51.97,
        "alt": -0.7,
        "STN": "348",
        "X": 123307.23288106002,
        "Y": 442498.5025450533,
        "Z": -0.7
    },
    "Gilze-Rijen": {
        "name": "Gilze-Rijen",
        "lon": 4.936,
        "lat": 51.566,
        "alt": 14.9,
        "STN": "350",
        "X": 123715.72551105871,
        "Y": 397548.013806552,
        "Z": 14.9
    },
    "Herwijnen": {
        "name": "Herwijnen",
        "lon": 5.146,
        "lat": 51.859,
        "alt": 0.7,
        "STN": "356",
        "X": 138384.09667411912,
        "Y": 430076.04933944705,
        "Z": 0.7
    },
    "Eindhoven": {
        "name": "Eindhoven",
        "lon": 5.377,
        "lat": 51.451,
        "alt": 22.6,
        "STN": "370",
        "X": 154290.9184795222,
        "Y": 384657.31016766746,
        "Z": 22.6
    },
    "Volkel": {
        "name": "Volkel",
        "lon": 5.707,
        "lat": 51.659,
        "alt": 22.0,
        "STN": "375",
        "X": 177127.99005000386,
        "Y": 407846.336823779,
        "Z": 22.0
    },
    "Ell": {
        "name": "Ell",
        "lon": 5.763,
        "lat": 51.198,
        "alt": 30.0,
        "STN": "377",
        "X": 181267.25444988074,
        "Y": 356578.54363042983,
        "Z": 30.0
    },
    "Maastricht": {
        "name": "Maastricht",
        "lon": 5.762,
        "lat": 50.906,
        "alt": 114.3,
        "STN": "380",
        "X": 181363.98355052367,
        "Y": 324093.5619647858,
        "Z": 114.3
    },
    "Arcen": {
        "name": "Arcen",
        "lon": 6.197,
        "lat": 51.498,
        "alt": 19.5,
        "STN": "391",
        "X": 211231.34266846423,
        "Y": 390198.61134588905,
        "Z": 19.5
    }
}

if __name__ == "__main__":
    weatherstations_daggegevens_data_file = 'weerstations_daggegevens_input.json'
    weatherstations_uurgegevens_data_file = 'weerstations_uurgegevens_input.json'
            
    with open(weatherstations_daggegevens_data_file) as json_file:
        weatherstations_daggegevens_dict = json.load(json_file)
        
    with open(weatherstations_uurgegevens_data_file) as json_file:
        weatherstations_uurgegevens_dict = json.load(json_file)        

    #fetchdata(weatherstationname='DeBilt', tminstr = '01-01-1980', tmaxstr = '31-10-2023', uurgegevens = False, knmiparams = 'PRCP')
    
     
    parse_knmi_datafile(datafilename='260_DeBilt_19800101_20231031.txt', tminstr = '01-01-1980', tmaxstr = '31-10-2023', knmiparams = 'PRCP')
        
    