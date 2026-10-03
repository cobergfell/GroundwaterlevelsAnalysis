# -*- coding: utf-8 -*-
"""
Created on Tue Jul 11 15:40:25 2023

@author: christophe
"""
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 12 11:03:11 2023

@author: christophe
"""
import scipy
import os
import numpy as np
from collections import OrderedDict
from matplotlib import pyplot as plt
from matplotlib.pyplot import figure, show,savefig,clabel,contour,setp,gcf,getp,gca,close
from datetime import datetime
from matplotlib.dates import date2num, num2date
from matplotlib.ticker import Formatter,FuncFormatter, NullLocator, FixedLocator
from utilities import orderOfMagnitude,generate_ticks,datestring2num,stdev_back_from_log
from logging import getLogger



    
class WellMetaData():
    
        
    logger = getLogger(__name__)

    
    def __init__(self, name = None, X = None, Y = None, Z = None, L = None, surface_level = None,
                 screen_top = None, screen_bottom = None, screen_number= None):
        

        """ 
        
        Parent class for well metadata

        Attributes
        ----------      
        X: float 
            X coordinate (easting)           
            
        Y: float 
            Y coordinate (northing)       
            
        Z: float (optional)
            Z coordinate (m above datum)       
            
        L: float (optional)
            Layer number (to be used if the well is modeled in a 
                          groundwater model)
            
            
        surface_level: float (optional)
            Surface level altitude         
            
        screen_top: float (optional)
            Screen top elevation      
            
            
        screen_bottom: float (optional)
            Screen bottom elevation  

        screen_number: float (optional)
            Screen number  
             

        """
        
        self._name = name
        self._X = X
        self._Y = Y
        self._Z = Z
        self._L = L
        self._surface_level = surface_level
        self._screen_top = screen_top 
        self._screen_bottom = screen_bottom  
        self._screen_number = screen_number
        

    def __repr__(self):
        
        """String representation of the monitoring well metadata."""

        clsname = str(self.__class__.__name__)
        string = f'{clsname}'
        string += f'(name = {self.name},'
        string += f'X = {self.X},'
        string += f'Y = {self.Y},'
        string += f'Z = {self.Z},'
        string += f'surface_level = {self.surface_level},'
        string += f'screen_top = {self.screen_top},'
        string += f'screen_bottom = {self.screen_bottom}'
        string += f')' 
        
 
        return string


    @property
    def name(self):
        """getter for updated time series name."""
        return self._name

    @name.setter
    def name(self,name):
        """setter for updated time series name."""
        self._name = name 
        
        
    @property
    def X(self):
        """getter for updated X coordinate."""
        return self._X

    @X.setter
    def X(self,X):
        """setter for updated X coordinate."""
        self._X = X     
        
        

    @property
    def Y(self):
        """getter for updated Y coordinate."""
        return self._Y

    @Y.setter
    def Y(self,Y):
        """setter for updated Y coordinate."""
        self._Y = Y 
        
        
        
    @property
    def Z(self):
        """getter for updated Z coordinate."""
        return self._Z

    @Z.setter
    def Z(self,Z):
        """setter for updated Z coordinate."""
        self._Z = Z         
        
        
    @property
    def L(self):
        """getter for updated layer number (if applicable)."""
        return self._Z

    @L.setter
    def L(self,L):
        """setter for updated layer number (if applicable)."""
        self._L = L         
    
    
    @property
    def surface_level(self):
        """getter for updated surface_level."""
        return self._surface_level

    @surface_level.setter
    def surface_level(self,surface_level):
        """setter for updated surface_level."""
        self._surface_level = surface_level         

    @property
    def screen_top(self):
        """getter for updated screen top."""
        return self._screen_top

    @screen_top.setter
    def screen_top(self,screen_top):
        """setter for updated screen_top."""
        self._screen_top = screen_top         
    


    @property
    def screen_bottom(self):
        """getter for updated screen bottom."""
        return self._screen_bottom

    @screen_bottom.setter
    def screen_bottom(self,screen_bottom):
        """setter for updated screen_bottom."""
        self._screen_bottom = screen_bottom
        
        
    @property
    def screen_number(self):
        """getter for updated screen number."""
        return self._screen_bottom

    @screen_number.setter
    def screen_number(self,screen_number):
        """setter for updated screen number."""
        self._screen_number = screen_number        

    
            
    def wellmetadataTostring(self):
        """ 
        A method to export the well metadata as string
        """         

        paddingright = 20
        justify = 16
        output_string = (f'\n\nMonitoring wells metadata\n'
            
                f'{"name =".ljust(justify)} {self.name:<{paddingright}}\n'     
                f'{"X =".ljust(justify)} {self._X:<{paddingright}.2f}\n'
                f'{"Y =".ljust(justify)} {self._Y:<{paddingright}.2f}\n'
                f'{"Z =".ljust(justify)} {self._Z:<{paddingright}.2f}\n'
                )
        
        if self._surface_level is not None:
            output_string += f'{"Surface level =".ljust(justify)} {self._surface_level:<{paddingright}.2f}\n'
        
        if self._screen_top is not None:
            output_string += f'{"Screen top =".ljust(justify)} {self._screen_top:<{paddingright}.2f}\n'                             
            
        if self._screen_bottom is not None:
            output_string += f'{"Screen bottom =".ljust(justify)} {self._screen_bottom:<{paddingright}.2f}\n'            
            
            
        output_string+= (f'\n\n')                
                
                
        return output_string
    
    
 

if __name__ == "__main__":
    
    metadata = WellMetaData(name = '28AP0093_1', X = 224550.34, Y = 488100.53, 
                            Z = 10., surface_level = 11.6,L = None )
    
    
    print('121 metadata: ',metadata.wellmetadataTostring())
    