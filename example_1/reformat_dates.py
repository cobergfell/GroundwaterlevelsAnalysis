# -*- coding: utf-8 -*-
"""
Created on Sun Sep 17 20:48:06 2023

@author: christophe
"""
# -*- coding: utf-8 -*-
"""
Created on Mon May  1 21:06:32 2023

@author: christophe
"""
import numpy as np
from datetime import datetime
from matplotlib.dates import date2num, num2date
import matplotlib.pyplot as plt
import os
from matplotlib.pyplot import figure, show,savefig,clabel,contour,setp,gcf,getp,gca
from matplotlib.dates import DayLocator,YearLocator,DateFormatter #HourLocator,Yearlocator
from matplotlib.ticker import Formatter,FuncFormatter, NullLocator, FixedLocator
from matplotlib import rcParams

from logging import getLogger

rcParams['text.usetex'] = False#True
logger = getLogger(__name__)


input_file_name = 'prec_giersbergen.csv'
input_file_name = 'evap_eindhoven.csv'

initial_format = '%m/%d/%Y %H:%M'
transformed_format = '%d-%m-%Y %H:%M'
skip = 4
export_path = None
separator = ','


format_datetime_templates = ['%Y/%m/%d %H:%M:%S','%d/%m/%Y %H:%M:%S','%d-%m-%Y %H:%M:%S','%d-%m-%Y %H:%M',
                        '%Y/%m/%d','%d/%m/%Y','%d-%m-%Y','%d-%m-%Y', '%m-%d-%Y %H:%M', '%m/%d/%Y %H:%M']

if transformed_format not in format_datetime_templates:
    modulename = str(__name__)   
    message = (f'\nIn module {modulename}.py: input date format not recognized'
               f' Default format %d-%m-%Y %H:%M:%S will be used instead.\n ')                      
    logger.warning(message) 


F = open(input_file_name)
all_lines = F.readlines()
F.close()


if export_path is None:
    curdir = os.getcwd()
    export_path = curdir+'//'+'reformated_'+input_file_name

F=open(export_path,'w')
for i in range(0,skip):
    F.write(all_lines[i])
    

for i in range(skip,len(all_lines)):
    line_splitted = all_lines[i].rstrip().split(separator)
    datetime_str = line_splitted[0]
    value = float(line_splitted[1])
    dt = datetime.strptime(datetime_str,initial_format)
    timenum = date2num(dt)        
    t = num2date(timenum)
    if transformed_format == '%d-%m-%Y %H:%M:%S':
        record = str(t.day)+'-'+str(t.month)+'-'+str(t.year)+'  '+str(t.hour).zfill(2)+':'+str(t.minute).zfill(2)+':'+str(t.second).zfill(2)+','+str('%8.5f' %value)
    elif transformed_format == '%d-%m-%Y %H:%M':
        record = str(t.day)+'-'+str(t.month)+'-'+str(t.year)+'  '+str(t.hour).zfill(2)+':'+str(t.minute).zfill(2)+','+str('%8.5f' %value)
    elif transformed_format == '%d/%m/%Y %H:%M:%S':
        record = str(t.day)+'/'+str(t.month)+'/'+str(t.year)+'  '+str(t.hour).zfill(2)+':'+str(t.minute).zfill(2)+':'+str(t.second).zfill(2)+','+str('%8.5f' %value)
    elif transformed_format == '%d/%m/%Y %H:%M':
        record = str(t.day)+'/'+str(t.month)+'/'+str(t.year)+'  '+str(t.hour).zfill(2)+':'+str(t.minute).zfill(2)+','+str('%8.5f' %value)
    elif transformed_format == '%m/%d/%Y %H:%M:%S':
        record = str(t.month)+'/'+str(t.day)+'/'+str(t.year)+'  '+str(t.hour).zfill(2)+':'+str(t.minute).zfill(2)+':'+str(t.second).zfill(2)+','+str('%8.5f' %value)
    elif transformed_format == '%m/%d/%Y %H:%M':
        record = str(t.month)+'/'+str(t.day)+'/'+str(t.year)+'  '+str(t.hour).zfill(2)+':'+str(t.minute).zfill(2)+','+str('%8.5f' %value)
    else:
        record = str(t.day)+'-'+str(t.month)+'-'+str(t.year)+'  '+str(t.hour).zfill(2)+':'+str(t.minute).zfill(2)+':'+str(t.second).zfill(2)+','+str('%8.5f' %value)
        
    F.write(record+'\n')
    
F.close() 
