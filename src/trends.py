# -*- coding: utf-8 -*-
"""
Created on Sat Jan 21 13:34:46 2023

@author: christophe
see also https://towardsdatascience.com/solid-coding-in-python-1281392a6a94
"""
import numpy as np
from abc import ABC, abstractmethod
from collections import OrderedDict
import scipy
from utilities import datestring2num

class TrendsFunctions(ABC):
    """
    An abstract class of response functions - not necessary yet but could serve as blueprint later
    to apply a bit more some of the 'SOLID' coding recomendations such as
    abstract interfaces.

    see  https://towardsdatascience.com/solid-coding-in-python-1281392a6a94
    
    """
    
    @abstractmethod
    def get_initial_parameters(self, suffix = ""):
        """
        See description of the method in UnitResponseFunctionsParent

        """
    pass


    
class TrendsFunctionsParent(TrendsFunctions):
    
    """
        A parent class that provides instances of trends (drainage base included)
    
        ...
    

    Methods
    -------
    
    get_initial_parameters()
        Method that returns the initial parameters

    trend()
        Method that returns a trend
    """            
    
    
    def __init__(self,p_dict = None, trend_type = "db", direction = None, tminstr = None, tmaxstr = None, suffix = ""):
        
        self.p_dict = p_dict
        self.trend_type = trend_type
        
        
    
    def get_initial_parameters(self, suffix = ""):

        """
        Method that returns the initial parameters 
        
        Parameters
        ----------
        
        suffix : string (optional)
            An optional string to disambiguate a parameter name (if necessary) 
            
        regime: string (optional)
            Specify the fluctuation regime, when applicable (default is an empty string)


        Returns
        -------
        parameters: python object of data type 'dict'
            A dictionary of initial parameters definitions        
    
        """          
        
    pass




    def trend(self):
        """Return the trend function .

        Parameters
        ----------
        p_dict: dictionary
            Dictionary with the values as floats representing the
            model parameters.


        Returns
        -------
        trend: array-like
            time series of trend
                    

            
        """
        pass    
    
        
class DrainageBase(TrendsFunctionsParent):
    
    """
    A method to simulate a constant drainage base

    Parameters
    ----------
    db: scalar
        value of the drainage base


    Note
    -----
    The drainage base is usually not physically a drainage base, it is a constant
    that 
    
    """
    def __init__(self, p_dict = None, trend_type = "db", suffix = ""):
        
        TrendsFunctionsParent.__init__(self,p_dict = p_dict, trend_type = trend_type,
                                       suffix = suffix)
        self.p_dict = p_dict
        self.trend_type = trend_type
        self.unit_response_name = 'drainage base'
        
    def get_initial_parameters(self, suffix = ""):
        trend_type = self.trend_type
        if suffix != "":
            suffix = '_' + suffix 

            
        parameters = OrderedDict()
        
        if trend_type == 'constant': # alternatively, some kind of periodic function could be introduced
            pname = "db" + suffix
            parameters[pname] = {
                "isvariable": True,
                "logtransform": True,
                "pname": pname,
                "minvalue": 0.0001,
                "maxvalue": 2000,
                "initvalue": 11.43, 
                }  

        return parameters      

            
class LinearTrend(TrendsFunctionsParent):
    
    """
    A method to simulate a linear trend

    Parameters
    ----------

    
    slope: scalar
        value of the slope


    Attributes
    ----------
        
    time_array: numpy array_like
        A one-dimensional numpy array of time numbers (time numbers as defined in matplotlib.dates 
        function date2num)
    
    slope: scalar (optional because this parameters can also be entered through the parameters dictionary, see next line)
        value of the slope    
    
    p_dict : python object of data type 'dict'
        the dictionary defining the model parameters
        and associated specifications such as associated indexes, names, logtransformed or not, 
        variable or fixed etc. See the the parameterslogistic.py module for additional details.
   
    
    tminstr : string (optional)
        The begin time of the trend as string (dd-mm-yyyy HH:MM:SS)

    tmaxstr : string (optional)
        The  end time of the trend as string (dd-mm-yyyy HH:MM:SS)  

    
    direction : string
        A string specifying the vertical direction of the trend 'upwards' or 'downwards'
    
    """
    
    def __init__(self, time_array = None, slope = None, p_dict = None, tminstr = None, tmaxstr = None, trend_type = "linear", direction = 'downwards', suffix = ""):
    
        TrendsFunctionsParent.__init__(self,p_dict = p_dict, trend_type = trend_type, 
                                       suffix = suffix)
        self.time = time_array
        self.slope = slope
        self.p_dict = p_dict
        self.tminstr = tminstr
        self.tmaxstr = tmaxstr        
        self.trend_type = trend_type
        self.direction = direction
        self.name = 'lintrend'

        if tminstr is None:
            tminstr = "01-01-1900"
            
        if tmaxstr is None:
            tmaxstr = "31-12-2100"            

        self.tmin = datestring2num(tminstr)
        self.tmax = datestring2num(tmaxstr)
        
        
        
    def get_initial_parameters(self, suffix = ""):
        trend_type = self.trend_type
        if suffix != "":
            suffix = '_' + suffix 

            
        parameters = OrderedDict()

        pname = "s" + suffix
        parameters[pname] = {
            "isvariable": True,
            "logtransform": True,
            "pname": pname,
            "minvalue": 1.e-10,
            "maxvalue": 1.e-1,
            "initvalue": 5.e-3, #cm/day
            }  

        return parameters      

                                        
    def evaluate(self): 
    
        """ 
        
        A method that evaluates the linear trend
        
        Return
        -------
        output: numpy array_like  
                output is the trend
        
         """

        time = self.time
        slope = self.slope
        p_dict = self.p_dict
        tmin = self.tmin
        tmax = self.tmax
        

        if p_dict is not None:
            #local copy of mutable input dictionaries
            _p_dict = {}
            for key in p_dict:
                _p_dict[key] = p_dict[key]
        
            p_indexes = _p_dict['p_indexes']
            logtransform = _p_dict['logtransform']
    
            i = p_indexes['lintrend']['regime_1']['basicparam'][0]
            slope = p_dict['p'][i]
            
            if p_dict['logtransform'][i] == True:
                slope  = np.exp(slope)
        
        tmin = max(tmin,time[0]) 
        tmax = min(tmax,time[-1]) 

        
        mask = (time[:] >= tmin) & (time[:] <= tmax)
        
        indexes = np.arange(0,len(time))
        selected_indexes = indexes[mask]
        ibegin = selected_indexes[0]
        iend = selected_indexes[-1]+1
        time_diff = time[1:]-time[:-1]
        time_diff = time_diff[ibegin:iend]
        time_diff_cumsum = np.cumsum(time_diff)
        trend = np.zeros(np.shape(time))
        trend[ibegin+1:iend] = time_diff_cumsum[:] * slope
        if self.direction == 'downwards':
            trend = -trend
        

        return trend      
                                  
if __name__ == "__main__":
    #curdir=os.getcwd()
    drainageBase = DrainageBase()

    print('298 drainageBase', drainageBase)
    