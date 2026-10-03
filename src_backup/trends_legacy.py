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
    
    Attributes
    ----------
        
    trend_type : string
        A string specifying the type of trend ('db' stands for drainage base)   

    
    Methods
    -------
    
    get_initial_parameters()
        Method that returns the initial parameters

    trend()
        Method that returns a trend
    """            
    
    
    def __init__(self,p_dict = None, trend_type = "db", suffix = ""):
        
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
        
        if trend_type == 'constant':
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

            
 
    

                                  
if __name__ == "__main__":
    #curdir=os.getcwd()
    drainageBase = DrainageBase()

    print('298 drainageBase', drainageBase)
    