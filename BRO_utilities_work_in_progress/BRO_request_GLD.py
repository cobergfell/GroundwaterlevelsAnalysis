# -*- coding: utf-8 -*-

"""Main module."""
import json
import pathlib
import logging

import requests
import pandas as pd
import pytz
import tqdm
import dateutil
from datetime import datetime

from logging import getLogger
import numpy as np
import os
from datetime import datetime
from matplotlib.dates import date2num, num2date
from matplotlib.ticker import Formatter,FuncFormatter, NullLocator, FixedLocator
from matplotlib import pyplot as plt
from matplotlib.pyplot import figure, show,savefig,clabel,contour,setp,gcf,getp,gca,close
# from utilities import orderOfMagnitude,generate_ticks,datestring2num,stdev_back_from_log
# from plotfunctions import simpleplot
from wellmetadata import WellMetaData
from bs4 import BeautifulSoup
logger = getLogger(__name__)



broId = "GLD000000002303"
# broId = "GLD000000014660"

url = f'https://publiek.broservices.nl//gm/gld//v1//objects//{broId}'
x = requests.get(url)
gldsoup = BeautifulSoup(x.text, "xml")

try:
    monitoringPoint = gldsoup.find_all("ns4:monitoringPoint", limit=1)[0].text  
except:
    monitoringPoint = gldsoup.find_all("monitoringPoint", limit=1)[0].text
    
brocom:broId    
tubeNumber = gldsoup.find_all("tubeNumber", limit=1)[0].text 


try:
    gldbroId = gldsoup.find_all("broId", limit=1)[0].text
except:
    gldbroId = gldsoup.find_all("brocom:broId", limit=1)[0].text 


gmwbroId = gldsoup.find_all("gldcommon:broId", limit=1)[0].text 


url = "https://publiek.broservices.nl//gm/gmw//v1//objects//"+gmwbroId
	

gmw = requests.get(url)
gmwsoup = BeautifulSoup(gmw.text)
print(gmwsoup.prettify())


nitgcode = gmwsoup.find_all("nitgcode", limit=1)[0].text
gmwbroid = gmwsoup.find_all("brocom:broid", limit=1)[0].text
xy = gmwsoup.find_all("gml:pos", limit=1)[0].text
Xstr,Ystr = xy.split(' ')


crs = gmwsoup.find_all("gmwcommon:location", limit=100)
srsname = crs[0].get("srsname")
epsg = srsname.split('::')[-1]


X = float(Xstr)
Y =  float(Ystr)


ind = int(tubeNumber)-1

groundlevelposition = gmwsoup.find_all("gmwcommon:groundlevelposition", limit=1)[ind].text

screenlength = gmwsoup.find_all("screenlength", limit=1)[ind].text
screentopposition = gmwsoup.find_all("screentopposition", limit=1)[ind].text
screenbottomposition = gmwsoup.find_all("screenbottomposition", limit=1)[ind].text
plaintubepartlength = gmwsoup.find_all("gmwcommon:plaintubepartlength", limit=1)[ind].text
tubetopposition = gmwsoup.find_all("tubetopposition", limit=1)[ind].text
verticaldatum = gmwsoup.find_all("gmwcommon:verticaldatum", limit=1)[ind].text
localverticalreferencepoint = gmwsoup.find_all("gmwcommon:localverticalreferencepoint", limit=1)[ind].text


screentopposition = float(screentopposition)
screenbottomposition = float(screenbottomposition)
groundlevelposition = float(groundlevelposition)

screen_top = screentopposition
screen_bottom = screenbottomposition
Z = 0.5 * (screenbottomposition + screenbottomposition)        
from wellmetadata import WellMetaData              
metadata = WellMetaData(name = nitgcode, X = X, Y = Y, Z = Z, surface_level = groundlevelposition, 
             screen_top = screentopposition, screen_bottom = screenbottomposition , screen_number= tubeNumber)                  
           

points = gldsoup.find_all("point")
N_data = len(points)  



observations_sequences_list = gldsoup.find_all("observation")
observations_sequences_indexes = np.empty((len(observations_sequences_list),2))
observations_sequences_indexes[:,0]=np.arange(0,len(observations_sequences_list))
j = -1
for observations_sequence in observations_sequences_list:
    j = j+1
    firstpoint = observations_sequence.find("point")
    firstpointtimestr = observations_sequence.find("time").text
    datestr, timestr = firstpointtimestr.split('T')
    timestr = timestr.split('+')[0]
    reshaped_datetimestr = datestr+' '+timestr
    dt = datetime.strptime(reshaped_datetimestr,'%Y-%m-%d %H:%M:%S')
    timenum = date2num(dt)  
    observations_sequences_indexes[j,1] = timenum
    
b = np.argsort(observations_sequences_indexes[:,1])#sort from low to high
#b = b[::-1]#sort from high to low
sorted_observations_sequences_indexes = observations_sequences_indexes[b,:]

    



records_list = []


c=-1
for i in range(0,len(observations_sequences_list)):#read the records time series records from ascii file
    index = int(sorted_observations_sequences_indexes[i,0])
    observations_sequence = observations_sequences_list[index]
    
    points = observations_sequence.find_all("point")
    for j in range(0,len(points)):
        point = points[j]
        c += 1
        if c%5000 == 0:
            print('158 c',c)
    
        datatimestr = point.find("time").text
        #2021-01-10T00:00:00+01:00
        value_str = point.find("value").text     
            
        datestr, timestr = datatimestr.split('T')
        timestr = timestr.split('+')[0]
        reshaped_datetimestr = datestr+' '+timestr
        dt = datetime.strptime(reshaped_datetimestr,'%Y-%m-%d %H:%M:%S')
        time_num = date2num(dt)  
        
        value = float(value_str)
        
        records_list.append([time_num,value])    
    
tseries = np.array(records_list)

# mask = tseries[:,1] < 1e9
# tseries = tseries[mask]

# mask = (tseries[:,0] > tmin) & (tseries[:,0] < tmax)
# tseries = tseries[mask]





    
#export to my own format for time series analysis

Xstr = str('%8.2f' %X)
Ystr = str('%8.2f' %Y)


begin_date = num2date(tseries[0,0])
begin_date_string= str(begin_date.day)+'-'+str(begin_date.month)+'-'+str(begin_date.year)
end_date=num2date(tseries[-1,0])
end_date_string=str(end_date.day)+'-'+str(end_date.month)+'-'+str(end_date.year)        
    

header1 ='nitgcode,gmwbroid,gldbroId,geographical_system,X,Y,verticaldatum,begin date,end date'
header2 = nitgcode+','+ gmwbroid +','+gldbroId +','+'epsg: '+epsg+','+Xstr+','+Ystr+','+verticaldatum+','+begin_date_string+','+end_date_string



abs_path = os.path.dirname(os.path.abspath(__file__))
abs_path_splitted = abs_path.split('\\')
path_to_parent_folder_elements = abs_path_splitted[:-1]
path_to_parent_folder = os.path.join(*path_to_parent_folder_elements)
splitted=path_to_parent_folder.split(':')#trick  to  repair  path
path_to_parent_folder = splitted[0]+':\\'+splitted[1]#trick  to  repair  path
path_to_resources_folder = os.path.join(path_to_parent_folder,'resources')



if path_to_parent_folder is not None:
    path_to_export_file =  os.path.join(path_to_parent_folder,gldbroId+'.txt')
else:
    curdir = os.getcwd()
    path_to_export_file = curdir+'//'+gldbroId+'.txt'

F=open(path_to_export_file,'w')
F.write(header1+'\n')
F.write(header2+'\n')
F.write('\n')
F.write('Date(dd-mm-yyyy HH:MM:SS),water level (m NAP)\n')
for i in range(0,len(tseries)):
    t = num2date(tseries[i,0])
    value = tseries[i,1]/100.
    record=str(t.day)+'-'+str(t.month)+'-'+str(t.year)+'  '+str(t.hour).zfill(2)+':'+str(t.minute).zfill(2)+':'+str(t.second).zfill(2)+','+str('%8.4f' %value)
    F.write(record+'\n')
F.close()
    





abs_path = os.path.dirname(os.path.abspath(__file__))
abs_path_splitted = abs_path.split('\\')
path_to_parent_folder_elements = abs_path_splitted[:-1]
path_to_parent_folder = os.path.join(*path_to_parent_folder_elements)
splitted=path_to_parent_folder.split(':')#trick  to  repair  path
path_to_parent_folder = splitted[0]+':\\'+splitted[1]#trick  to  repair  path
path_to_resources_folder = os.path.join(path_to_parent_folder,'resources')



if path_to_parent_folder is not None:
    path_to_export_file =  os.path.join(path_to_resources_folder,gldbroId+'.xml')
else:
    curdir = os.getcwd()
    path_to_export_file = curdir+'//'+gldbroId+'.xml'

F=open(path_to_export_file,'w')
F.write(gldsoup.prettify())
F.close()
