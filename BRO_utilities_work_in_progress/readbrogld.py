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






def read_from_BRO_GLD(path_to_file):
    
    """
    A function to read broservices GLD files
   

    Attributes
    ----------
    
    path_to_file: string
            string specifying the pth to the xml file to read
            
            
    Returns
    -------
    tseries: numpy array_like
        the time series to of groundwater levels in a numpy array 
        The time numbers in colomn 0 are time numbers as defined in matplotlib.dates 
        function date2num
            
            
            

    References
    ---------
    https://realpython.com/python-xml-parser/

    """    
 
    
    
    # with open('GLD000000014660_IMBRO.xml', 'r') as f:
    #     data = f.read()
        
    with open(path_to_file, 'r') as f:
        data = f.read()        
        
    
    
    gldsoup = BeautifulSoup(data, "xml")
    # print('75 gldsoup.prettify()', gldsoup.prettify())
    # input()
    
    monitoringPoint = gldsoup.find_all("ns4:monitoringPoint", limit=1)[0].text  
    tubeNumber = gldsoup.find_all("tubeNumber", limit=1)[0].text 
    
    
    gldbroId = gldsoup.find_all("broId", limit=1)[0].text
    
    gmwbroId = gldsoup.find_all("ns10:broId", limit=1)[0].text 
    
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
    Z = 0.5 * (screenbottomposition + screenbottomposition)        
           
    metadata = WellMetaData(name = nitgcode, X = X, Y = Y, Z = Z, surface_level = groundlevelposition, 
                 screen_top = screentopposition, screen_bottom = screenbottomposition , screen_number= tubeNumber)                  
               
    # first sort the measurements sequences chronologically (because they are not chronological in the BRO!)
    observations_sequences_list = gldsoup.find_all("observation")
    observations_sequences_indexes = np.empty((len(observations_sequences_list),2))
    observations_sequences_indexes[:,0]=np.arange(0,len(observations_sequences_list))
    j = -1
    for observations_sequence in observations_sequences_list:
        j = j+1
        # firstpoint = observations_sequence.find("point")
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

    return tseries,metadata



def plot(tseries,metadata, tminstr = '01-01-1900', tmaxstr = '31-12-2100', save_plot = False):
    
    """ 
    A function to plot the head time series (override the plot method of TimeSeriesBase)
    
    Parameters
    ----------    

            
            
    tminstr: string (optional)
            string specifying the minimum time on the plot (format dd-mm-yyyy)
            
    tmaxstr: string (optional)
            string specifying the maximum time on the plot (format dd-mm-yyyy)     
            
    save_plot: boolean (optional)
        Boolean variable specifying if the plot is to be saved on disk
        

    Returns
    -------
    ax: Matplotlib Axes object
            The Matplotlib Axes object

    """
    

    left_margin = 0.12
    right_margin = 0.06
    top_margin = 0.12
    bottom_margin = 0.44
    dy = 0.04
    dx = 0.05        

    fs = 12
 
    
    t = lambda t_str : date2num(datetime.strptime(t_str,'%d-%m-%Y')) # lamda function to transform time string into time num

    tmin = t(tminstr)
    tmax = t(tmaxstr)

    mask = (tseries[:,0] >= tmin) & (tseries[:,0] < tmax)
    tseries = tseries[mask]        


    yas_title = 'm above datum'

    plot_length = 1-left_margin-right_margin
    plot_height = 1-bottom_margin-top_margin
    xll = left_margin # x coordinate lower left corner
    yll = 1-top_margin-plot_height # y coordinate lower left corner
    
    
    fig=figure(num=None, figsize=(8,10),dpi=50,facecolor='w', edgecolor='k')# figsize=(width, height) in inches.

    ax = fig.add_axes([xll,yll,plot_length,plot_height],frameon=True, xscale=None, yscale=None)  
    
    leg_list = []
    ax.plot_date(tseries[:,0],tseries[:,1],'b-',linewidth = 1.0)     
    leg_list.append('observed')
            
    def format_date(dates, pos=None):
        return num2date(dates).strftime('%Y')
    
    # def format_date(dates, pos=None):
    #     return num2date(dates).strftime('%d-%m-%Y %H:%M')        

    ax.set_ylabel(yas_title,fontsize=fs)
  
    ax.set_xlabel('date (yyyy)',fontsize=fs)
    ax.xaxis.set_major_formatter(FuncFormatter(format_date))

    
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(fs)
        tick.label.set_rotation(0)
        # tick.label.set_rotation(45)
    
    for label in ax.xaxis.get_majorticklabels():
        #ha=label.get_horizontalalignment()
        label.set_horizontalalignment('center')
        
    name = metadata.name
    X = metadata.X
    Y = metadata.Y
    Z = metadata.Z
    surface_level = metadata.surface_level
    screen_top = metadata.screen_top
    screen_bottom = metadata.screen_bottom
    
    time = tseries[:,0]
    lower_expect, median, upper_expect = generate_gxg(tseries,metadata,tminstr = tminstr, tmaxstr = tmaxstr, alpha = 0.05)
    
    surface_level_array = surface_level*np.ones(len(tseries))
    
    upper_expect_array = upper_expect*np.ones(len(tseries))
    median_array = median*np.ones(len(tseries))
    lower_expect_array = lower_expect*np.ones(len(tseries))
    
    
    ax.plot_date(time,upper_expect_array,'c--',linewidth=2.0) 
    leg_list.append('95% upper quantile')
    ax.plot_date(time,median_array,'m-',linewidth=2.0)  
    leg_list.append('median')
    ax.plot_date(time,lower_expect_array,'g--',linewidth=2.0)         
    leg_list.append('95% lower quantile')
    
    
    leg = ax.legend((leg_list),ncol=3,bbox_to_anchor=(0.85,-0.24),frameon=False)
    
    Xleg = 0.05
    Yleg = -0.22
    
    ax.annotate('Legend', xy=(Xleg, Yleg), xycoords='axes fraction',horizontalalignment='left', verticalalignment='center',fontsize=fs)#this is for 1 column 2 rows
    leg = ax.legend((leg_list), ncol=2,shadow=False,bbox_to_anchor=[Xleg, Yleg-0.2], loc='lower left',frameon=False)
    ltext  = leg.get_texts() 
    plt.setp(ltext, fontsize=fs)    # the legend text fontsize

    title='Time series of groundwater heads for piezometer '+ name
    ax.set_title(title,fontsize=fs)       
    (x, y) = ax.title.get_position()
    #ax.title.set_y(0.95 * y)    
    ax.grid(True)               
    

    
    ax.set_ylabel('Meter above datum',fontsize=fs)
    ax.set_xlabel('Date (yyyy)',fontsize=fs)
    #ax.set_xlabel('Datum (yyyy)',fontsize=fs)#030919
    ax.xaxis.set_major_formatter(FuncFormatter(format_date))
    #fig.autofmt_xdate()
    #ax.xaxis.set_label_position('right')

    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(fs)
        tick.label.set_rotation(45)
    
    for label in ax.xaxis.get_majorticklabels():
        #ha=label.get_horizontalalignment()
        label.set_horizontalalignment('center')          
    
    ax.grid(True)
    
    #English version

    list_of_header_items=[
            'Name',
            'X',
            'Y',
            'Z (m above datum)',
            'Surface_level (m above datum)',
            '5% low quantile (m above datum)',
            'Median (m above datum)',
            '5% high quantile (m above datum)']


    # #Dutch version
    # list_of_header_items=[
    # 'Naam',
    # 'X-coordinaat (RD)',
    # 'Y-coordinaat (RD)',
    # 'Filterdiepte (m NAP)',
    # 'Maaiveldhoogte (m NAP)',
    # '5'+'%'+' laag kwantiel (m NAP)',
    # 'Mediaan (m NAP)',
    # '5'+'%'+' hoog kwantiel (m NAP)']            

    
    list_of_values=[
    name,
    str('%10.2f' %X),
    str('%10.2f' %Y),
    str('%10.2f' %Z),
    str('%10.2f' %surface_level),
    str('%10.2f' %lower_expect),
    str('%10.2f' %median),
    str('%10.2f' %upper_expect)]    
    
    
    
    # X=[0.30,0.60]
    # Y=[-0.37]
    X = [Xleg,Xleg+0.5]
    Y = [Yleg-0.3]
    
    
    delta_y=0.07
    for i in range(1,len(list_of_header_items)):    
        last_y=Y[-1]                      
        Y.append(last_y-dy)

    for i in range(0,len(list_of_header_items)):
        ax.annotate(list_of_header_items[i], xy=(X[0],Y[i]), xycoords='axes fraction',horizontalalignment='left', verticalalignment='center')
        ax.annotate(list_of_values[i], xy=(X[1],Y[i]), xycoords='axes fraction',horizontalalignment='left', verticalalignment='center')
    show()
    

        


    
    show()
        

    if save_plot == True:
        figname = "plot_of_{}".format(name)

        try: 
            abs_path=os.path.dirname(os.path.abspath(__file__))
            abs_path_splitted=abs_path.split('\\')
            path_to_parent_folder_elements=abs_path_splitted[:-1]
            path_to_parent_folder=os.path.join(*path_to_parent_folder_elements)
            splitted=path_to_parent_folder.split(':')#trick  to  repair  path
            path_to_parent_folder=splitted[0]+':\\'+splitted[1]#trick  to  repair  path
            path_to_file_folder=os.path.join(path_to_parent_folder,'results')
            savefig(path_to_file_folder+'\\'+figname, dpi=None, facecolor='w', edgecolor='w',orientation='portrait', format=None,transparent=False, bbox_inches=None, pad_inches=0.1)
          
        except:
            curdir=os.getcwd()
            savefig(curdir+'\\'+figname, dpi=None, facecolor='w', edgecolor='w',orientation='portrait', format=None,transparent=False, bbox_inches=None, pad_inches=0.1)
                   
    

        
    return ax       



def generate_gxg(tseries,metadata, tminstr = '01-01-1900', tmaxstr = '31-12-2100', alpha = 0.05):
    
    """
    To be used to estimate some standard quantiles of  groundwater levels time series
    
    (gxg refers to the Dutch habit to call these means GLG, GG, GHG which stand
     for gemiddelde laag grondwaterstand, gemiddelde grondwaterstand, 
     and gemiddelde hoog grondwaterstand)

    Parameters
    ----------
    
    tminstr: string (optional)
            string specifying the minimum time on the plot (format dd-mm-yyyy)
            
    tmaxstr: string (optional)
            string specifying the maximum time on the plot (format dd-mm-yyyy)    
        
    alpha: float
        quantile determining the low/high groundwater level
        


    Returns
    -------
    
    glg : float
        mean low groundwaterlevel ('gemiddelde lage grondwaterstand')
        
    gg: float
        mean  groundwaterlevel ('gemiddelde grondwaterstand')
        
    ghg : float
        mean high groundwaterlevel ('gemiddelde hoge grondwaterstand')
        
    """

    t = lambda t_str : date2num(datetime.strptime(t_str,'%d-%m-%Y')) # lamda function to transform time string into time num

    tmin = t(tminstr)
    tmax = t(tmaxstr)

    mask = (tseries[:,0] >= tmin) & (tseries[:,0] < tmax)
    tseries = tseries[mask]    
        
    
    mask = (tseries[:,0] >= tmin) & (tseries[:,0] < tmax)
    tseries = tseries[mask]
    
    glg = None
    gg = None
    ghg = None
    
   

    gg = np.median(tseries[:,1])

    #sort the time series in order to determine the quantiles
    b = np.argsort(tseries[:,1])#sort from low to high
    b = b[::-1]#sort from high to low
    sorted_serie = tseries[b,:]
    N = len(tseries)

    index = int((1-alpha)*N)# groundwater levels were sorted from high to low
    glg = sorted_serie[index,1]
    
    index=int(alpha*N)# groundwater levels were sorted from high to low
    ghg = sorted_serie[index,1]

    if (tmax-tmin) <365 or len(tseries) < 24:
        logger.warning("The selected time series is too short to estimates mean levels")        
        

    paddingright = 20
    justify = 16
    output_string = (f'\n\nMonitoring wells metadata\n'
            f'{"name =".ljust(justify)} {metadata.name:<{paddingright}}\n'     
            f'{"X =".ljust(justify)} {metadata.X:<{paddingright}.2f}\n'
            f'{"Y =".ljust(justify)} {metadata.Y:<{paddingright}.2f}\n'
            f'{"Z =".ljust(justify)} {metadata.Z:<{paddingright}.2f}\n'
            f'{"Surface level =".ljust(justify)} {metadata.surface_level:<{paddingright}.2f}\n'               
            )
    
    if ghg is not None:
        output_string += f'{"GHG =".ljust(justify)} {ghg:<{paddingright}.2f}\n'
    if gg is not None:
        output_string += f'{"GG =".ljust(justify)} {gg:<{paddingright}.2f}\n'
    if glg is not None:
        output_string += f'{"GLG =".ljust(justify)} {glg:<{paddingright}.2f}\n'            



    print(output_string)

    return glg,gg,ghg            
    




if __name__ == "__main__":
    
    abs_path = os.path.dirname(os.path.abspath(__file__))
    # abs_path_splitted = abs_path.split('\\')
    # path_to_parent_folder_elements = abs_path_splitted[:-1]
    # path_to_parent_folder = os.path.join(*path_to_parent_folder_elements)
    # splitted=path_to_parent_folder.split(':')#trick  to  repair  path
    # path_to_parent_folder = splitted[0]+':\\'+splitted[1]#trick  to  repair  path
    path_to_resources_folder = os.path.join(abs_path,'resources_BRO_GLD')


    files_names_list = os.listdir(path_to_resources_folder)
    for filename in files_names_list: 
        path_to_file = os.path.join(path_to_resources_folder,filename)
        tseries,metadata = read_from_BRO_GLD(path_to_file)    
        plot(tseries,metadata)
        glg,gg,ghg = generate_gxg(tseries,metadata)