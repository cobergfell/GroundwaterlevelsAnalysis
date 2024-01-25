# -*- coding: utf-8 -*-
"""
Created on Fri Jan 27 12:29:58 2023

@author: christophe
"""


from logging import getLogger
import numpy as np
import os
from matplotlib.pyplot import figure, show,savefig,clabel,contour,setp,gcf,getp,gca,close
from matplotlib.dates import date2num, num2date,DayLocator,MonthLocator,YearLocator,DateFormatter #HourLocator,Yearlocator
from datetime import datetime
from matplotlib.ticker import Formatter,FuncFormatter, NullLocator, FixedLocator 

from timeseriesbase import TimeSeriesBase
from plotfunctions import simpleplot
from rcParams import rcParams
from utilities import orderOfMagnitude,generate_ticks,datestring2num,file_name_from_path, get_path_to_parent_directory

    
logger = getLogger(__name__)

class Stress(TimeSeriesBase):
    
    """
    
    A subclass of TimeSeriesBase that provides instances of time series of stress observations
    such as precipitation, evaporation, pumping or river stages.

    Attributes in supplement to the inherited ones 
    ----------------------------------------------
    
    X : float (also accessible via metadata dictionary)
        The X geographical coordinate of the measurement point (easting)
    
    Y : float (also accessible via metadata dictionary)
        The Y geographical coordinate of the measurement point (northing)  
        
    Z : float (only for pumping wells, also accessible via metadata dictionary)
        Indicative depth of the pumping well screen in height units above datum
        (typically the average between top and bottom of the pumping well screen)
    
    L : float (only for pumping wells if response function is a multilayer 
               groundwater flow model - not implemented yet)
        When the response function is a multi-layer groundwater flow model, L indicates the 
        layer number, with 0 corresponding to the top layer
        
    Methods
    -------
    read_from_csv(filepath = None,tminstr = None, tmaxstr = None, settings = None)
        A class method to be used to generate instances of the Stress class from simple csv ascii files
        

    """


             
    def __init__(self, name = None, X = None, Y = None, Z = None, L = None, observed = None,
                 tminstr = None, tmaxstr = None, tminnum = None, tmaxnum = None, cumulative = True, 
                 stress_type = None, settings = None, metadata = None):


            TimeSeriesBase.__init__(self, observed = observed, name = name, settings = settings, 
                                    metadata = metadata, cumulative = cumulative,
                                    tminstr = tminstr, tmaxstr = tmaxstr, tseries_type = stress_type)         

            self.name = name
            #self.path_to_parent_directory = get_path_to_parent_directory()
            self.stress_type = stress_type
            self.metadata = metadata
            self.cumulative = cumulative
            self._cumulated = None # cumulated values over time
            self.X = X
            self.Y = Y
            self.Z = Z #only for wells
            self.L = L #only for wells
            self._metadata = metadata    
            self._observed = observed     
            self._preprocessed = None
            self.settings = settings
            self.tmin = tminnum
            self.tmax = tmaxnum
            self.tminstr = tminstr
            self.tmaxstr = tmaxstr            
            


    def __repr__(self):
        
        """String representation of the stress."""
        string = f'class: {self.__class__.__name__}\n'
        string += f'stress_type: {self.stress_type}\n'
        string += f'name: {self.name}\n'
        string += f'attributes: {list(self.__dict__.keys())}'
        
        return string
          
           
    def format_pump_rate(x):
        return '%1.0f' % (x*1e-3)
                

    @property
    def cumulated(self):
        """getter for updated cumulated observations time series."""
        return self._cumulated

    @cumulated.setter
    def cumulated(self, cumulated):
        """Setter for updated cumulated observations time series."""
        self._cumulated = cumulated          
    
           
    @classmethod      
    def read_from_csv(cls ,filepath, tminstr = None, tmaxstr = None, settings = None, stress_type = None, cumulative = None):
        
        
        
        """
        a simple method to read prec series
        input: serie_file is an ascii time series formated according to following template:
  
        STRESS,Time serie name,X COORD,Y COORD,
        PREC,260_De_Bilt,0,0
        
        DATE,VALUE (meter per day)
        1-1-1906  23:45,  0.000000
        2-1-1906  23:45,  0.000000
        3-1-1906  23:45,  0.000000
        4-1-1906  23:45,  0.003600
            
            
        Parameters
        ----------
        filepath: csv file containing a list of observation time with associated oberved value
        name: str, optional
            String with the name of the time series which by default is the file name
        settings: dict, optional, with model detting other than default values
    
    
        Returns
        -------
        series: array_like
                array_like time series .

        """

        
        if cumulative is None:
            if stress_type in ['prec','evap']:
                cumulative = True
            else:
                cumulative = False 
        else: 
            cumulative = cumulative
  
        
        settings = settings
        tmin = -1e9
        tmax = +1e9
        
   
        if tminstr is not None:
            tmin = datestring2num(tminstr)
                        
        if tmaxstr is not None:
            tmax = datestring2num(tmaxstr)  
            
            
        F = open(filepath)
        all_lines = F.readlines()
        F.close()
        
        #read the header line and collect the metadata in a dictionary (a dictionnary is a Python data type)
        header_line = all_lines[1]
        header_line_splited = header_line.rstrip().split(',')

        metadata = {}
        metadata['type'] = stress_type
        name = header_line_splited[1]
        metadata['name'] = name
        metadata['X'] = header_line_splited[2]
        metadata['Y'] = header_line_splited[3]
        metadata['Z'] = None
        Z = None
        
        if stress_type =='pump':
            try:
                top_str = header_line_splited[4]
                bot_str = header_line_splited[5]
                metadata['screen top'] = top_str
                metadata['screen bottom'] = bot_str
                Z = 0.5 * (float(top_str) + float(bot_str))
                if Z is not None:           
                    metadata['Z'] = str(Z)
            except:
                message = (f'\nCould not find pump depth data for pump {name}.\n'
                f'Tip: you can enter it line 2 of your ascii file as metadata as follows:\n'
                f'stress_type,time_serie_name,Xcoord,Ycoord,well_screen_top_(meter_above_datum),'
                f'well_screen_botttom_(meter_above_datum),comments\n')
                logger.warning(message)
                

    
        #read the time series in an array of two columns: 1st column is the time (in hours since an arbitrary begin) 2nd column is the measured groundwater level (by default in meter NAP)
        N_lines=len(all_lines)# This is the total number of lines of the ascii file (standard format)
        N_lines_to_skip=4#number of lines to skip prior to time series records
        N_data=N_lines-N_lines_to_skip#this is the number of measurements, skipping the 4 upper lines
        tseries=np.empty((N_data,2))#declare an array of the desired size
        
        # if stress_type in ['prec','evap']:
        #     formatt = '%m/%d/%Y %H:%M'
        # else:
        #     formatt = '%d-%m-%Y %H:%M'                                                                                                                                                                                                                       
        for i in range(0,N_data):#read the records time series records from ascii file
            j=i+N_lines_to_skip
            data_line=all_lines[j].rstrip().split(',')
            date_time_string=data_line[0]
            value_str=data_line[1] 
            dt=datetime.strptime(date_time_string,'%d-%m-%Y %H:%M')
            # dt=datetime.strptime(date_time_string,'%m/%d/%Y %H:%M')
            time_num=date2num(dt)    
            # dt=datetime.strptime(date_time_string,formatt)
            time_num=date2num(dt)            
            
            # time_num = datestring2num(date_time_string) # this is slow! beter to uncheck the two lines above
            tseries[i,0]=time_num#time numerical value
            tseries[i,1]=float(value_str)#string to float for measured groundwater level
            
        mask = (tseries[:,0] > tmin) & (tseries[:,0] < tmax)
        tseries = tseries[mask]            

        #return cls(tseries,metadata=metadata)           
        # return cls(path_to_filepath, observed = tseries,stress_type = stress_type, settings = settings, metadata=metadata) 

        return cls(name = metadata['name'], X = float(metadata['X']), Y = float(metadata['Y']), Z = Z, 
                   metadata = metadata, observed = tseries, tminstr = tminstr, tmaxstr = tmaxstr, 
                   tminnum = tmin, tmaxnum = tmax, settings = settings, stress_type = stress_type, cumulative = cumulative)        
    
     


                           
if __name__ == "__main__":
    #curdir=os.getcwd()
    abs_path=os.path.dirname(os.path.abspath(__file__))
    abs_path_splitted=abs_path.split('\\')
    path_to_parent_folder_elements=abs_path_splitted[:-1]
    path_to_parent_folder=os.path.join(*path_to_parent_folder_elements)
    splitted=path_to_parent_folder.split(':')#trick  to  repair  path
    path_to_parent_folder=splitted[0]+':\\'+splitted[1]#trick  to  repair  path
    path_to_file_folder=os.path.join(path_to_parent_folder,'resources')
    path_to_file=os.path.join(path_to_file_folder,'260_De_Bilt_PREC_19010101_20200910.csv')
    path_to_file = os.path.join(path_to_file_folder,'672_PREC_Hellendoorn_19510101_20160731_corrected_for_interception2mm.csv')
    # path_to_file = os.path.join(path_to_file_folder,'prec_giersbergen.csv')
    
    
        
    tminstr = '01-01-1982 00:00:00'
    tmaxstr = '31-12-2005 00:00:00'

    ts = Stress.read_from_csv(path_to_file, stress_type = 'prec', tminstr = tminstr, tmaxstr = tmaxstr)
    ts.plot()
    harmonic = ts.harmonic_observed


    ax = ts.plot_harmonics()
    Ampl_prec = ts.harmonic_observed.amplitude[0]
    phi_prec = ts.harmonic_observed.phase_shift[0]  
    
    print('272 Ampl_prec',Ampl_prec)
    print('273 phi_prec',phi_prec)    
    
    print('276',repr(ts))