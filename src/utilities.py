# -*- coding: utf-8 -*-
"""
Created on Wed Mar  8 15:33:59 2023

@author: christophe
"""

from datetime import datetime
from matplotlib.dates import date2num, num2date
from logging import getLogger
import numpy as np
import os

logger = getLogger(__name__)


"""  
This module contains a set of general utility functions to be used in different
modules

"""
    
def p_dict_copy(P_dict):
    """  
    A utility function to generate a deep copy of the p_dict dictionary
    
    """
    
    P_dict_copy = {}
    for key1 in P_dict:
        if isinstance(P_dict[key1],dict):
            P_dict_copy[key1] = {}
            for key2 in P_dict[key1]:
                if isinstance(P_dict[key1][key2],dict):
                    P_dict_copy[key1][key2] = {}
                    for key3 in P_dict[key1][key2]:
                        if isinstance(P_dict[key1][key2][key3],dict):
                            P_dict_copy[key1][key2][key3] = {}                        
                            for key4 in P_dict[key1][key2][key3]:
                                    P_dict_copy[key1][key2][key3][key4] = P_dict[key1][key2][key3][key4]                        
                        else:
                            P_dict_copy[key1][key2][key3] = P_dict[key1][key2][key3] 
                else:
                    P_dict_copy[key1][key2] = P_dict[key1][key2]      
        else:
            P_dict_copy[key1] = P_dict[key1]                           

  
    return P_dict_copy



def determine_format_datetime(datetimestring):
    """ 
    A utility function to try determining a date time format
    

    Parameters
    ----------
        
    datetimestring: string
        time given as string

    datetimeformat: string, optional
        Date format. Example:
        '%Y/%m/%d %H:%M:%S'
        '%d/%m/%Y %H:%M:%S'
        '%d-%m-%Y %H:%M:%S'
        '%d-%m-%Y %H:%M',
        '%Y/%m/%d',
        '%d/%m/%Y',
        '%d-%m-%Y',
        '%m-%d-%Y %H:%M',
        '%m/%d/%Y %H:%M',
        
   
    Return
    ------
    datetimestring: string
        the detected time format

        
    """  
    
    """
        Date stringes given as list (for possible use in a modified version of this function)
        format_datetime_templates = ['%Y/%m/%d %H:%M:%S','%d/%m/%Y %H:%M:%S','%d-%m-%Y %H:%M:%S','%d-%m-%Y %H:%M',
                        '%Y/%m/%d','%d/%m/%Y','%d-%m-%Y','%d-%m-%Y', '%m-%d-%Y %H:%M', '%m/%d/%Y %H:%M']
    """

    
    format_datetime = None

    if len(datetimestring.split(' ')) > 1:

        timestr = datetimestring.split(' ')[1]
        if len(timestr.split(':')) > 2:
            timeformat = ' %H:%M:%S'
        else:
            timeformat =' %H:%M'
    else:
        timeformat = ''

    if len(datetimestring.split('/')) > 1:
        sep = '/'
    else:
        sep = '-'
    
    first_date_element = datetimestring.split(sep)[0]

    if len(first_date_element) > 2:
        format_datetime = '%Y'+sep+'%m'+sep+'%d'+timeformat
    else:
        format_datetime = '%d'+sep+'%m'+sep+'%Y'+timeformat

    return format_datetime
            
            
            
            



def datestring2num(datetimestring, datetimeformat = None):

    """ 
    A utility function to transform a date string to a 
    matplotlib.dates time number
    

    Parameters
    ----------
        
    datetimestring: string
        time given as string

    datetimeformat: string, optional
        Date format. Example:
        '%Y/%m/%d %H:%M:%S'
        '%d/%m/%Y %H:%M:%S'
        '%d-%m-%Y %H:%M:%S'
        '%d-%m-%Y %H:%M',
        '%Y/%m/%d'
        '%d/%m/%Y'
        '%d-%m-%Y'
        '%d-%m-%Y'
        '%m-%d-%Y %H:%M'
        '%m/%d/%Y %H:%M'
        
   
    Return
    ------
    tmin : float
        time 
        (time number as defined in matplotlib.dates function date2num)

        
    """     
    
    
    #t = lambda t_str,format_date : date2num(datetime.strptime(t_str,format_date))
    
    datetimenum = None
    defaultformat = '%d-%m-%Y %H:%M:%S'
    
    if datetimeformat is not None:
        datetimeformat = datetimeformat

    try:
        dt = datetime.strptime(datetimestring,datetimeformat)
        datetimenum = date2num(dt)  
    except Exception:
        try:
            datetimeformat = determine_format_datetime(datetimestring)
            dt = datetime.strptime(datetimestring,datetimeformat) 
            datetimenum = date2num(dt) 
        except Exception:        
            message = (f'\nutilities.py module: Datetime format of {datetimestring} could not be parsed\n')
            logger.warning(message)

    return datetimenum






    
def unpack_time_boundaries(time_boundaries_dic):
    
    """ 
    A utility function to translate minimum and maximum dates and times 
    entered as input in different format into dates numbers

    Parameters
    ----------
    time_boundaries_dic: python object of data type 'dict'
            a dictionary of minimum times and dates
    
    time_boundaries_dic["tminnum"]: float, optional
        minimum time given as matplotlib.dates time number
        
    time_boundaries_dic["tminstr"]: string, optional
        minimum time given as string

    time_boundaries_dic["tmaxnum"]: float, optional
        maximum time given as matplotlib.dates time number
        
    time_boundaries_dic["tmaxstr"]: string, optional
        maximum time given as string
        
        
    Return
    ------
    tmin : float
        The minimum time 
        (time number as defined in matplotlib.dates function date2num)

    tmax : float
        The maximum time 
        (time number as defined in matplotlib.dates function date2num)     
        
    """    

    tminnum = time_boundaries_dic["tminnum"]
    tminstr = time_boundaries_dic["tminstr"]
    tmaxnum = time_boundaries_dic["tmaxnum"]
    tmaxstr = time_boundaries_dic["tmaxstr"]    

    
    t = lambda t_str : date2num(datetime.strptime(t_str,'%d-%m-%Y')) # lamda function to transform time string into time num

    if tminnum is not None and tminstr is not None:
        logger.warning("Minimum time is both given as string and numeric value."
                          "String value will be ignored")

    if tminnum is None and tminstr is None:
        tmin = t('01-01-1900')

    if tminstr is not None:
        tmin = t(tminstr)

    if tminnum is not None:   
         tmin = tminnum

    if tmaxnum is not None and tmaxstr is not None:
        logger.warning("Maximum time is both given as string and numeric value."
                          "String value will be ignored")
        
    if tmaxnum is None and tmaxstr is None:
        tmax = t('31-12-2100')
                    
    if tmaxstr is not None:
        tmax = t(tmaxstr)
     
    if tmaxnum is not None:   
         tmax = tmaxnum      
    
    return tmin,tmax



def update_p_dict_from_p(p_dict,p):
    
    """ 
    
    A utility function to update p_dict given a vector of optimized parameters
    

    Parameters
    ----------
    p_dict          python object of data type 'dict'
                    dictionary containing parameters names, initial values, optimized values,
                    whether parameetrs are logtansformed or native and whether parameters are varying of fixed
    
    p :             numpy array_like
                    vector of parameters optimized values 

        
    Returns
    -------
    p_dict          python object of data type 'dict'
                    updated P_dict
        
        
    """ 
    
    
    p_complete_updated = np.array(p_dict['p_init'])
    mask = np.array(p_dict['isvariable'])
    p_complete_updated[mask] = p
    p_dict['p'] = p_complete_updated.tolist()
    
    return p_dict


def get_path_to_parent_directory():
    
    """  
    A utility function to get the path of the parent directory of a file
    
    """
    
    abs_path = os.path.dirname(os.path.abspath(__file__))
    abs_path_splitted = abs_path.split('\\')
    path_to_parent_directory_elements = abs_path_splitted[:-1]
    path_to_parent_directory = os.path.join(*path_to_parent_directory_elements)
    splitted=path_to_parent_directory.split(':') # part 1 of a trick to repairpath
    path_to_parent_directory = splitted[0]+':\\'+splitted[1] # part  2 of a trick to repair path
    return path_to_parent_directory


def file_name_from_path(path_string):
    
    """  
    A utility function to extract the file name from its complete path
    
    """    
    path_splitted = path_string.split('\\')
    default_file_name = path_splitted[-1].split('.')[0]
    return default_file_name




def orderOfMagnitude(number):
    
    """  
    A utility function to determine the order of magnitude of a number
    
    """     
    import math
    return math.floor(math.log(number, 10))




def generate_ticks(vmin,vmax):
    
    """  
    A utility function to generate convenients y=as tick to be used for a plot
    
    """  
    
    
    import math
    valuesrange = vmax - vmin
    valuesrangemagnitude = orderOfMagnitude(valuesrange)
    #scaledvaluesrange = valuesrange/np.power(10,valuesrangemagnitude)
   
    if vmin > 0:
        magnitudevmin = orderOfMagnitude(vmin)
        tickmin = math.floor(float(vmin)/np.power(float(10),magnitudevmin) ) * np.power(float(10),magnitudevmin)
    elif vmin == 0:
        tickmin = 0
    else:
        magnitudevmin = orderOfMagnitude(-vmin)
        tickmin = - math.floor(float(-vmin)/np.power(float(10),magnitudevmin) ) * np.power(float(10),magnitudevmin)
    
    
    if vmax > 0:
        magnitudevmax = orderOfMagnitude(vmax)
        tickmax = math.floor(float(vmax)/np.power(float(10),magnitudevmax) ) * np.power(float(10),magnitudevmax)
    elif vmax == 0:
        tickmax = 0
    else:
        magnitudevmax = orderOfMagnitude(-vmax)
        tickmax = - math.floor(float(-vmax)/np.power(float(10),magnitudevmax) ) * np.power(float(10),magnitudevmax)    

    
    if valuesrangemagnitude >= 1:
        exponent = valuesrangemagnitude-1
        increment = 0.5 * np.power(float(10),exponent)
    else:
        increment = 0.5

    vmax = vmax + 2 *  increment    
    ticks = np.arange(tickmin,vmax,increment)

    
    # if scaledvaluesrange <=5 :
    #     ticks = np.arange(tickmin,tickmin + 5 * np.power(10,valuesrangemagnitude),0.5* np.power(10,valuesrangemagnitude))
    # else :
    #     ticks = np.arange(tickmin,tickmin + 10* np.power(10,valuesrangemagnitude),np.power(10,valuesrangemagnitude))

    return ticks


        

def stdev_back_from_log(mu,sigma):
    from numpy.random import normal,multivariate_normal,uniform,triangular#import more than necessary, just in case we want to try something else

    """ 
    
    A utility function to transform the log-transformed expected values of a parameter
    and associated estimated standard back to natural form, assuming a log normal distribution 
    of the expected value


    Parameters
    ----------
    mu :        float
                expected value of the log-transformed variable
    
    sigma :     float
                standard deviation of the log-transformed variable

        
    Returns
    -------
    mu_b :      float
                expected value of the back transformed variable
    
    sigma_b :   float
                standard deviation of the back transformed variable
                
    Based on:
    # https://math.stackexchange.com/questions/176196/calculate-the-expected-value-of-y-ex-where-x-sim-n-mu-sigma2
    # https://towardsdatascience.com/log-normal-distribution-a-simple-explanation-7605864fb67c                    
        
    """ 
    
    mu_b = np.exp(mu+0.5*sigma**2)
    variance = (np.exp(sigma**2)-1)*np.exp(2*mu+sigma**2)
    sigma_b = np.sqrt(variance)

    return mu_b,sigma_b


def stdev_back_from_log_stochastic(mu,sigma):
    from numpy.random import normal,multivariate_normal,uniform,triangular#import more than necessary, just in case we want to try something else

    """ 
    A utility function to transform the log-transformed expected values of a parameter
    and associated estimated standard back to natural form, assuming a log normal distribution 
    of the expected value, based on a stochastic generation of values


    Parameters
    ----------
    mu :        float
                expected value of the logtransformed variable
    
    sigma :     float
                standard deviation of the logtransformed variable

        
    Returns
    -------
    mu_b :      float
                expected value of the back transformed variable
    
    sigma_b :   float
                standard deviation of the back transformed variable
    
    """ 

    N = 1000
    location = np.atleast_1d(mu) # mean
    scale = np.atleast_1d(sigma) # stdev
    if len(location) >1:
        log_realizations = multivariate_normal(location, pow(scale,2), size=N)
    else:
        log_realizations = normal(loc=location, scale=scale, size=N)
        
    back_transformed = np.empty(N)
    for i in range(0,N):
        back_transformed[i] = np.exp(log_realizations[i])
    
    mu_b = np.mean(back_transformed)
    sigma_b = np.sqrt(np.var(back_transformed))
    
    return mu_b,sigma_b



        


def basename(filepath):#make a dictionary of basenames from fseries files

    """  
    A utility function to extract the basename of a stress time series 
    (path name stripped from extension and stripped from all element from the extended path)
    
    """ 

    a = filepath.split('\\')#from the whole file name, we extract the base name of the well or other stress
    a = a[len(a)-1].split('.')
    basename = a[0]
    return basename


def reformat_date(input_file_name = None, initial_format = None, transformed_format = None, 
                  skip = None, export_path = None, separator = None):


    # input_file_name = 'prec_giersbergen.csv'
    # initial_format = '%m/%d/%Y %H:%M'
    # transformed_format = '%d-%m-%Y %H:%M'
    # skip = 4
    # export_path = None
    # separator = ','
    
    
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
            record = str(t.day)+'-'+str(t.month)+'-'+str(t.year)+'  '+str(t.hour).zfill(2)+':'+str(t.minute).zfill(2)+':'+str(t.second).zfill(2)+','+str('%8.2f' %value)
        elif transformed_format == '%d-%m-%Y %H:%M':
            record = str(t.day)+'-'+str(t.month)+'-'+str(t.year)+'  '+str(t.hour).zfill(2)+':'+str(t.minute).zfill(2)+','+str('%8.2f' %value)
        elif transformed_format == '%d/%m/%Y %H:%M:%S':
            record = str(t.day)+'/'+str(t.month)+'/'+str(t.year)+'  '+str(t.hour).zfill(2)+':'+str(t.minute).zfill(2)+':'+str(t.second).zfill(2)+','+str('%8.2f' %value)
        elif transformed_format == '%d/%m/%Y %H:%M':
            record = str(t.day)+'/'+str(t.month)+'/'+str(t.year)+'  '+str(t.hour).zfill(2)+':'+str(t.minute).zfill(2)+','+str('%8.2f' %value)
        elif transformed_format == '%m/%d/%Y %H:%M:%S':
            record = str(t.month)+'/'+str(t.day)+'/'+str(t.year)+'  '+str(t.hour).zfill(2)+':'+str(t.minute).zfill(2)+':'+str(t.second).zfill(2)+','+str('%8.2f' %value)
        elif transformed_format == '%m/%d/%Y %H:%M':
            record = str(t.month)+'/'+str(t.day)+'/'+str(t.year)+'  '+str(t.hour).zfill(2)+':'+str(t.minute).zfill(2)+','+str('%8.2f' %value)
        else:
            record = str(t.day)+'-'+str(t.month)+'-'+str(t.year)+'  '+str(t.hour).zfill(2)+':'+str(t.minute).zfill(2)+':'+str(t.second).zfill(2)+','+str('%8.2f' %value)
            
        F.write(record+'\n')
        
    F.close() 

        