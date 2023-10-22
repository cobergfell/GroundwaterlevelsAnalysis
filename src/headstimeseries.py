# -*- coding: utf-8 -*-
"""
Created on Fri Jan 27 12:29:58 2023

@author: christophe
"""
from logging import getLogger
import numpy as np
import os
from timeseriesbase import TimeSeriesBase
from datetime import datetime
from matplotlib.dates import date2num, num2date
from matplotlib.ticker import Formatter,FuncFormatter, NullLocator, FixedLocator
from abc import ABC, abstractmethod
from collections import OrderedDict
from matplotlib import pyplot as plt
from matplotlib.pyplot import figure, show,savefig,clabel,contour,setp,gcf,getp,gca,close
from utilities import orderOfMagnitude,generate_ticks,datestring2num,stdev_back_from_log
from plotfunctions import simpleplot


logger = getLogger(__name__)

class Heads(TimeSeriesBase):
    
    """
    
    A subclass of TimeSeriesBase that provides instances of time series of head (= groundwater levels) observations

    Attributes in supplement to the inherited ones 
    ----------------------------------------------
    
    X : float (also accessible via metadata dictionary)
        The X geographical coordinate of the piezometer (easting)
    
    Y : float (also accessible via metadata dictionary)
        The Y geographical coordinate of the piezometer (northing)  
        
    Z : float (also accessible via metadata dictionary)
        Indicative depth of the well screen in height units above datum
        (typically the average between top and bottom of the well screen)
    
    L : float (only of use when the response function is a groundwater flow model)
        When the response function is a groundwater flow model, L indicates the 
        layer number, with 0 corresponding to the top layer
        
    Methods
    -------
    read_from_csv(csv_file = None,tminstr = None, tmaxstr = None, settings = None)
        A class method to be used to generate instances of the Heads class from simple csv ascii files
        
    read_from_DINO(file_name, tminstr = None, tmaxstr = None, separator=',', settings = None):        
        A class method to be used to generate instances of the Heads class from the DINOloket format
        (from the Dutch national geological survey institute)
        
    generate_gxg(tminstr= None, tmaxstr = None, alpha = 0.05)
        A method to calculate the mean groundwater levels (low, mean, high)

    generate_groundwaterlevels_stats_string()
        A method to output basic groundwater level statistics

    generate_output(model_definition = None)
        A method to output piezometer characteristics and modeling results (when available) 


    """

    
    def __init__(self, name = None, X = None, Y = None, Z = None, L = None, metadata = None,
                 observed = None,tminstr = None, tmaxstr = None,tminnum = None, tmaxnum = None, 
                 settings = None, p_dict = None, stresses_dict = None):
        


        TimeSeriesBase.__init__(self, observed = observed, name = name, settings = settings, 
                                    metadata = metadata,tminstr = tminstr, tmaxstr = tmaxstr, tseries_type = 'heads')      


        #local copy of mutable input dictionaries
        if p_dict is not None:
            _p_dict = {}
            for key in p_dict:
                _p_dict[key] = p_dict[key]
        else: 
            _p_dict = None
            
        if stresses_dict is not None:            
            _stresses_dict = {}
            for key in stresses_dict:
                _stresses_dict[key] = stresses_dict[key]      
        else:
            _stresses_dict = None

        self.name = name
        self.monitoringwell = None
        self.X = X
        self.Y = Y
        self.Z = Z
        self.L = L
        self._metadata = metadata    
        self._observed = observed
        self._simulated = None
        self._simulated_with_residuals = None
        self._residuals = None
        self._noise = None
        
        self._fitgoodness = None
        self._fitgoodness_with_residuals = None        
        self._preprocessed = None
        self._p_dict = _p_dict
        self._stresses_dict = _stresses_dict
        self.targets_selector = None
        self.targets = None
        self.weigths = None
        self.settings = settings
        self.tmin = tminnum
        self.tmax = tmaxnum
        self.tminstr = tminstr
        self.tmaxstr = tmaxstr
        
        
        fitgoodness = {}
        fitgoodness['expvar'] = None
        fitgoodness['pcov'] = None
        fitgoodness['pcor'] = None
        fitgoodness['pstdev'] = None
        self._fitgoodness = fitgoodness
        
        
        fitgoodness_with_residuals = {}
        fitgoodness_with_residuals['expvar'] = None
        fitgoodness_with_residuals['pcov'] = None
        fitgoodness_with_residuals['pcor'] = None
        fitgoodness_with_residuals['stdev'] = None
        self._fitgoodness_with_residuals = fitgoodness_with_residuals

        
        if metadata.name is not None:         
            self.name = metadata.name
        else:
            self.name = name
            
        if metadata.surface_level is not None:         
            self.surface_level = float(metadata.surface_level)
        else:
            self.surface_level = None
             
    
        tmin = -1e9
        tmax = +1e9
        
   
        if tminstr is not None:
            tmin = datestring2num(tminstr)
                        
        if tmaxstr is not None:
            tmax = datestring2num(tmaxstr)               
        
        if self.tmin is None:
            self.tmin = tmin
            
        if self.tmax is None:
            self.tmax = tmax          
  
       
        

    def __repr__(self):
        
        """String representation of the monitoring well."""

        clsname = str(self.__class__.__name__)
        string = f'{clsname}'
        string += f'(name = {self.name}, X = {self.X}, Y = {self.Y})\n'
        return string


    @property
    def observed(self):
        """getter for updated observed time series."""
        return self._observed

    @observed.setter
    def observed(self,observed):
        """setter for updated observed time series."""
        self._observed = observed 


    @property
    def simulated(self):
        """getter for updated simulated time series."""
        return self._simulated

    @simulated.setter
    def simulated(self,simulated):
        """setter for updated simulated time series."""
        self._simulated = simulated    
        
    @property
    def simulated_with_residuals(self):
        """getter for updated simulated time series with residuals."""
        return self._simulated

    @simulated_with_residuals.setter
    def simulated_with_residuals(self,simulated_with_residuals):
        """setter for updated simulated time series with residuals."""
        self._simulated_with_residuals = simulated_with_residuals    
        
        
    @property
    def residuals(self):
        """getter for updated residuals time series."""
        return self._residuals

    @residuals.setter
    def residuals(self,residuals):
        """setter for updated residuals."""
        self._residuals = residuals            
        
        
    @property
    def noise(self):
        """getter for updated noise time series."""
        return self._noise

    @noise.setter
    def noise(self,noise):
        """setter for updated noise."""
        self._noise = noise          
        
        
    @property
    def fitgoodness(self):
        """getter for goodness of fit."""
        return self._fitgoodness

    @fitgoodness.setter
    def fitgoodness(self,fitgoodness):
        """setter for goodness of fit."""

        self._fitgoodness = fitgoodness      

        
    @property
    def fitgoodness_with_residuals(self):
        """getter for goodness of fit with simulated residuals."""
        return self._fitgoodness_with_residuals

    @fitgoodness_with_residuals.setter
    def fitgoodness_with_residuals(self, fitgoodness_with_residuals):
        """setter for goodness of fit with simulated residuals."""
        self._fitgoodness_with_residuals = fitgoodness_with_residuals             
        
        
        
    @property
    def preprocessed(self):
        """getter for updated preprocessed time series."""
        return self._preprocessed

    @preprocessed.setter
    def preprocessed(self,preprocessed):
        """setter for updated preprocessed time series."""
        self._preprocessed = preprocessed      
        
        
    @property
    def preprocessed(self):
        """getter for updated preprocessed time series."""
        return self._preprocessed

    @preprocessed.setter
    def preprocessed(self,preprocessed):
        """setter for updated preprocessed time series."""
        self._preprocessed = preprocessed           
        
                
    @property
    def p_dict(self):
        """getter for updated p_dict."""
        return self._p_dict

    @p_dict.setter
    def p_dict(self,p_dict):
        """setter for updated p_dict."""
        self._p_dict = p_dict  
        
        
    @property
    def stresses_dict(self):
        """getter for updated stresses_dict."""
        return self._stresses_dict

    @stresses_dict.setter
    def stresses_dict(self,stresses_dict):
        """setter for updated stresses_dict."""
        self._stresses_dict = stresses_dict          
        
               
    @classmethod      
    def read_from_csv(cls, filepath,tminstr = None, tmaxstr = None, settings = None):
        
        
        """
        A class method to read groundwater head time series as ascci (csv) 
        files and generate a Heads object

        
        Parameters
        ----------
        filepath: string
            path to the csv file

        tminstr: string (optional)
                string specifying the minimum time on the plot (format dd-mm-yyyy)
                
        tmaxstr: string (optional)
                string specifying the maximum time on the plot (format dd-mm-yyyy)    
                
        settings : python object of data type 'dict' (optional)
            An optional dictionary of settings             
        
        Example
        --------
        Example of the csv file:
  
        screen name,screen number,X,Y,surface level (m NAP),screen top (m NAP), screen bottom (m NAP),begin date,end date
        28AP0093_1,1,224550.34,488100.53,11.6,2.6,1.6,,
    
        Datum(dd-mm-yyyy HH:MM:SS),groundwater level (m above datum)
        14-08-1978 12:00,8.44
        28-08-1978 12:00,8.4
        14-09-1978 12:00,8.33       
    
        Returns
        -------
        cls: Heads object
            An instance of Heads

        """
                        
        from wellmetadata import WellMetaData
        from matplotlib.pyplot import figure, show,savefig
        from datetime import datetime
        from matplotlib.dates import date2num, num2date
        from matplotlib.ticker import Formatter,FuncFormatter, NullLocator, FixedLocator
        import numpy as np
      
        settings = settings
        tmin = -1e9
        tmax = +1e9
        
   
        if tminstr is not None:
            tmin = datestring2num(tminstr)
                        
        if tmaxstr is not None:
            tmax = datestring2num(tmaxstr)  
        
        #open ascii file, read all lines as a list of lines and close the ascii file
        F = open(filepath)
        all_lines = F.readlines()
        F.close()
        
        #read the header line and collect the metadata in a dictionary (a dictionnary is a Python data type)
        header_line = all_lines[1]
        header_line_splited = header_line.rstrip().split(',')
        if len(header_line_splited) < 2:
            header_line_splited = header_line.rstrip().split(';')
        
        try:
            name = header_line_splited[0]
        except:
            message = (f'\nCould not read name from header line {header_line}\n')
            logger.error(message)
            
        try:
            screen_number = int(header_line_splited[1])
        except:
            message = (f'\nCould not read screen number from header line {header_line}\n'
                       f'\nScreen number is set to 1.\n' 
                       )
            logger.warning(message)
            screen_number = 1
                      
        try:
            X = float(header_line_splited[2])
        except:
            message = (f'\nCould not read X coordinate from header line {header_line}\n')
            logger.error(message)        
            
        try:
            Y = float(header_line_splited[3])
        except:
            message = (f'\nCould not read Y coordinate from header line {header_line}\n')
            logger.error(message)               
            
        try:
            surface_level = float(header_line_splited[4])
        except:
            message = (f'\nCould not read surface level from header line {header_line}\n')
            logger.error(message) 

        try:
            screen_top = float(header_line_splited[5])
        except:
            message = (f'\nCould not read screen top from header line {header_line}\n')
            logger.error(message) 
            
        try:
            screen_bottom = float(header_line_splited[6])
        except:
            message = (f'\nCould not read screen bottom from header line {header_line}\n')
            f'\nScreen bottom is set equal to screen top.\n' 
            logger.warning(message)  
            screen_bottom = screen_top
            
        Z = 0.5 * (screen_top + screen_bottom)
     
        metadata = WellMetaData(name = name, X = X, Y = Y, Z = Z, surface_level = surface_level )
    
        #read the time series in an array of two columns: 1st column is the time (in hours since an arbitrary begin) 2nd column is the measured groundwater level (by default in meter NAP)
        N_lines = len(all_lines) # This is the total number of lines of the ascii file (standard format)
        N_lines_to_skip = 4 # number of lines to skip prior to water levels measurements records
        N_data = N_lines-N_lines_to_skip # this is the number of measurements, skipping the 4 upper lines
        heads = np.empty((N_data,2)) # declare an array of the desired size
        
                                                                                                                                                                                                                                                    
        for i in range(0,N_data):#read the records time series records from ascii file
            j = i+N_lines_to_skip
            data_line = all_lines[j].rstrip().split(',')
            if len(data_line) < 2:
                data_line = all_lines[j].rstrip().split(';')
            try:
                # dt=datetime.strptime(data_line[0],'%d-%m-%Y %H:%M:%S')
                # dt=datetime.strptime(data_line[0],'%d-%m-%Y %H:%M')
                # time_num=date2num(dt)  
                time_num = datestring2num(data_line[0], datetimeformat = '%d-%m-%Y %H:%M:%S')
            except :
                logger.warning("Date-time format not recognized")            
            
            value_str = data_line[1]     
            heads[i,0] = time_num#time numerical value
            heads[i,1] = float(value_str)#string to float for measured groundwater level
    
   
        mask = heads[:,1] < 1e9
        heads = heads[mask]
        
        mask = (heads[:,0] > tmin) & (heads[:,0] < tmax)
        heads = heads[mask]
        
        
        if settings is not None:
            if settings['sampling_interval'] is not None:
                Nskip = settings['sampling_interval']
                indices = np.arange(0,len(heads),Nskip)
                heads = np.take(heads, indices, axis=0, out=None, mode='raise')
            
        if heads == []:
            logger.warning("The observation time series is empty, "
                                "check if time series and time boundaries are overlpping")
        
        return cls(name = name, X = X, Y = Y, Z = Z, metadata = metadata, 
                   observed = heads, tminstr = tminstr, tmaxstr = tmaxstr, tminnum = tmin, tmaxnum = tmax,settings = settings)           
    
       
    @classmethod      
    def read_from_DINO(cls,filepath, tminstr = None, tmaxstr = None, separator = ',', settings = None):
        
        """
        A class method to read groundwater head time series from am ascii file 
        in the DINOloket format (from the Dutch national geological survey institute)
        and generate a Heads object.

        
        Parameters
        ----------
        filepath: string
            path to the csv file

        tminstr: string (optional)
                string specifying the minimum time on the plot (format dd-mm-yyyy)
                
        tmaxstr: string (optional)
                string specifying the maximum time on the plot (format dd-mm-yyyy)    
                
        settings : python object of data type 'dict' (optional)
            An optional dictionary of settings             
        
         Example
        ----------   
        
        Example of the DINO format:
                
        Titel:,,,,,,,,,,,
        Gebruikersnaam:,,,,,,,,,,,
        Periode aangevraagd:,01-01-1800,tot:,15-05-2017,,,,,,,,
        Gegevens beschikbaar:,29-03-1971,tot:,27-05-2011,,,,,,,,
        Datum: ,15-05-2017,,,,,,,,,,
        Referentie:,NAP,,,,,,,,,,
        
        NAP:,Normaal Amsterdams Peil,,,,,,,,,,
        MV:,Maaiveld,,,,,,,,,,
        MP:,Meetpunt,,,,,,,,,,
        
        Locatie,Filternummer,Externe aanduiding,X-coordinaat,Y-coordinaat,Maaiveld (cm t.o.v. NAP),Datum maaiveld gemeten,Startdatum,Einddatum,Meetpunt (cm t.o.v. NAP),Meetpunt (cm t.o.v. MV),Bovenkant filter (cm t.o.v. NAP),Onderkant filter (cm t.o.v. NAP)
        B28C0314,001,28CP0022,227531,486117,1375,29-03-1971,29-03-1971,09-02-1987,1412,37,775,675
        B28C0314,001,28CP0022,227531,486117,1375,29-03-1971,09-02-1987,14-07-1988,1398,23,775,675
        B28C0314,001,28CP0022,227531,486117,1375,29-03-1971,14-07-1988,12-02-1997,1394,19,775,675
        B28C0314,001,28CP0022,227531,486117,1375,29-03-1971,12-02-1997,10-09-2003,1432,57,775,675
        B28C0314,001,28CP0022,227531,486117,1377,10-09-2003,10-09-2003,27-05-2011,1445,68,775,675
        
        
        Locatie,Filternummer,Peildatum,Stand (cm t.o.v. MP),Stand (cm t.o.v. MV),Stand (cm t.o.v. NAP),Bijzonderheid,Opmerking,,,
        B28C0314,001,29-03-1971,502,465,910,,,,,,
        B28C0314,001,14-04-1971,506,469,906,,,,,,
        B28C0314,001,28-04-1971,512,475,900,,,,,,
        B28C0314,001,14-05-1971,516,479,896,,,,,,    
    
    
        Returns
        -------
        cls: Heads object
            An instance of Heads

        """
        
        from wellmetadata import WellMetaData
        from datetime import datetime
        from matplotlib.dates import date2num, num2date
        from datetime import datetime
        from numpy import median
        
        
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
        
            
        #first determine if it is NAP or 'maaiveld'format
        stop = 0
        i = -1
        while stop == 0:#
            i = i+1
            if all_lines[i][:3] == 'Ref':
                stop = 1
                if 'NAP' in all_lines[i]:
                    reference = 'NAP'
                else:
                    reference = 'maaiveld'


        indices = []
        for i in range(0,30):#assume metadata bloc is within first 30 lines
            if all_lines[i][:3] in ["Loc","LOC"]:
                indices.append(i)
     
        metadata_bloc_begin = indices[0]+1
        
        condition = True
        i = metadata_bloc_begin
        while condition == True:
            i+=1
            if len(all_lines[i][0]) < 2 :
                condition = False
            if i > metadata_bloc_begin +30:
                condition = False
        metadata_bloc_end = i  
        
        # metadata_bloc_end = indices[1] - 1
        # data_bloc_begin = metadata_bloc_end + 3
        data_bloc_begin = indices[-1]
        
        
        data_bloc_end = len(all_lines)
        number_of_metadata_lines = metadata_bloc_end - metadata_bloc_begin
        Nrecords = len(all_lines)-data_bloc_begin

        
        if reference == 'NAP':
            
            list_of_metadata = []#create a list of dictionaries of metadat
            for i in range(metadata_bloc_begin,metadata_bloc_end):
                metadataline = all_lines[i]
                metadata_list = metadataline.rstrip().split(',')

            try:
                name = metadata_list[0]+'-'+metadata_list[1][1:]
            except:
                message = (f'\nCould not read name from line {metadataline}\n')
                logger.error(message)
                
            try:
                screen_number = int(metadata_list[1])
            except:
                message = (f'\nCould not read screen number from line {metadataline}\n'
                           f'\nScreen number is set to 1.\n' 
                           )
                logger.warning(message)
                screen_number = 1
                          
            try:
                X = float(metadata_list[3])
            except:
                message = (f'\nCould not read X coordinate from line {metadataline}\n')
                logger.error(message)        
                
            try:
                Y = float(metadata_list[4])
            except:
                message = (f'\nCould not read Y coordinate from line {metadataline}\n')
                logger.error(message)                
                
            try:
                surface_level = float(metadata_list[5])/100
            except:
                message = (f'\nCould not read surface level from line {metadataline}\n')
                logger.error(message) 
    
            try:
                screen_top = float(metadata_list[9]) / 100. 
            except:
                message = (f'\nCould not read screen top from line {metadataline}\n')
                logger.error(message) 
                
            try:
                screen_bottom = float(metadata_list[10]) / 100. 
            except:
                message = (f'\nCould not read screen bottom from line {metadataline}\n')
                f'\nScreen bottom is set equal to screen top.\n' 
                logger.warning(message)  
                screen_bottom = screen_top
                
            try:
                Z = 0.5 * (screen_top + screen_bottom)
            except:    
                message = (f'\nScreen mean altitude Z couls not be evaluated\n')
                logger.error(message)             
            
            
            metadata = WellMetaData(name = name, X = X, Y = Y, Z = Z, surface_level = surface_level )                   

            list_of_metadata.append(metadata)
            

                
        elif reference=='maaiveld':
            #to do check if metadata are different if reference is maaiveld
            pass
        
        
        metadata = list_of_metadata[-1]#the last metadata dictionaryis used later as reference when meta data are needed
    

        heads = np.empty((Nrecords,2))
        c = -1
        for i in range(data_bloc_begin,data_bloc_end):
            c = c+1
            line = all_lines[i]
            d = line.split(',')

            try:
                time_num = datestring2num(d[2])

            except :
                message = (f'\nDatetime format of {d[2]} could not be parsed\n')
                logger.warning(message)
                
            try:
                value = float(d[5])/100 #divide by hundred because values in DINO are given in cm

            except :
                value = 1e9
            heads[c,0] = time_num
            heads[c,1] = value

        mask = heads[:,1] < 1e9
        heads = heads[mask]
        
        mask = (heads[:,0] > tmin) & (heads[:,0] < tmax)
        heads = heads[mask]
        
        if settings is not None:
            if settings['sampling_interval'] is not None:
                Nskip = settings['sampling_interval']
                indices = np.arange(0,len(heads),Nskip)
                heads = np.take(heads, indices, axis=0, out=None, mode='raise')
        

        return cls(name = name, X = X, Y = Y, Z = Z, metadata = metadata, 
                   observed = heads, tminstr = tminstr, tmaxstr = tmaxstr, tminnum = tmin, tmaxnum = tmax, settings = settings)        
    
    
    def generate_gxg(self, tminstr= None, tmaxstr = None, alpha = 0.05):
        
        """
        A method to estimate the mean groundwater levels of a time series
        with levels expressed as height above some datum
        
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

        tseries = self._observed
        if tminstr is not None:
            tmin = datestring2num(tminstr)
        else:
            tmin = tseries[0,0]
            
            
        if tmaxstr is not None:
            tmax = datestring2num(tmaxstr)
        else:
            tmax = tseries[-1,0]       
  
            
        
        mask = (tseries[:,0] >= tmin) & (tseries[:,0] < tmax)
        tseries = tseries[mask]
        
        glg = None
        gg = None
        ghg = None
        
        if (tmax-tmin) >2*365 and len(tseries) > 24:
            
            

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
          
            
            
        else:
            self.logger.warning("The selected time series is too short to estimates mean levels")        
            
    
        return glg,gg,ghg            
    
    
    def generate_groundwaterlevels_stats_string(self):

        """
        A method to generate the groundwater level statistics as strings 
    
        """  

        paddingright = 20
        justify = 16
        output_string = (f'\n\nMonitoring wells metadata\n'
                f'{"name =".ljust(justify)} {self.name:<{paddingright}}\n'     
                f'{"X =".ljust(justify)} {self.X:<{paddingright}.2f}\n'
                f'{"Y =".ljust(justify)} {self.Y:<{paddingright}.2f}\n'
                f'{"Z =".ljust(justify)} {self.Z:<{paddingright}.2f}\n'
                f'{"Surface level =".ljust(justify)} {self.ground_surface_level:<{paddingright}.2f}\n'               
                )
        
        if self.ghg is not None:
            output_string += f'{"GHG =".ljust(justify)} {self.ghg:<{paddingright}.2f}\n'
        if self.gg is not None:
            output_string += f'{"GG =".ljust(justify)} {self.gg:<{paddingright}.2f}\n'
        if self.glg is not None:
            output_string += f'{"GLG =".ljust(justify)} {self.glg:<{paddingright}.2f}\n'        
        
        return output_string
            
        
        
        
        
        
            
            
    def generate_output(self, model_definition = None):

        """

        A method to generate the output string
        
        Parameters
        ----------
        model_definition : Modeldefinition object (optional)
            The model definition can be needed to select which items
            need to be outputed
    
        """        
        
        

        _stresses_dict = self._stresses_dict

        paddingright = 20
        justify = 16
        output_string = (f'\n\nMonitoring wells metadata\n'
            
                f'{"name =".ljust(justify)} {self.name:<{paddingright}}\n'     
                f'{"X =".ljust(justify)} {self.X:<{paddingright}.2f}\n'
                f'{"Y =".ljust(justify)} {self.Y:<{paddingright}.2f}\n'
                )
        
        if self.surface_level is not None:
            output_string += f'{"Z =".ljust(justify)} {self.Z:<{paddingright}.2f}\n'
           
        if self.surface_level is not None:
            output_string += f'{"Surface level =".ljust(justify)} {self.surface_level:<{paddingright}.2f}\n'
            
        if self._p_dict is not None:
            pnames = self._p_dict['pnames']
            popt = self._p_dict['p']
            isvariable = self._p_dict['isvariable']
            logtransform = self._p_dict['logtransform']

            if self._fitgoodness is not None:
                pcor = self._fitgoodness['pcor']
                pcov = self._fitgoodness['pcov']
                pstdev = self.fitgoodness['pstdev']
                expvar = self.fitgoodness['expvar']
                
                
                emptystr = ''
                paddingright = 12
                signspace = ' '
                
                
                output_string+=(f'\nExplained variance: {expvar:.2f}\n')
                output_string+=(f'\nOptimized parameters\n\n')

            
                j = -1
                for i in range(0,len(popt)):
                    output_string+= (f'{pnames[i]:<{paddingright}}:  ')
                    mu = popt[i]

                    if isvariable[i] == True:
                        j += 1
                        sigma = pstdev[j]
                        if logtransform[i] == True:
                            mu,sigma= stdev_back_from_log(mu,sigma)
                        
                        output_string+= (f'{mu:<{paddingright}.2e}')              
                        output_string+= (f'stdev = {sigma:<{paddingright}.2e}\n')    
                    else: 
                        if logtransform[i] == True:
                            mu = np.exp(mu)
                        output_string+= (f'{mu:<{paddingright}.2e}')
                        
                        mystring = 'fixed'
                        paddingright_modified = paddingright - len(mystring)
                        output_string+= (f'{mystring:<{paddingright_modified}}\n')       
                        
                output_string+= (f'\n\nCovariance matrix\n')  
                emptystr = ''
                paddingright = 12
                signspace = ' '
                
                output_string+= (f'{emptystr:<{paddingright}}')
                
                for i in range(0,len(popt)):
                    if isvariable[i] == True:
                        output_string+= (f'{pnames[i]:<{paddingright}}')
                k = -1        
                for i in range(0,len(popt)):
                    if isvariable[i] == True:
                        k += 1
                        l = -1
                        output_string+= (f'\n{pnames[i]:<{paddingright}}')
                        for j in range(0,len(popt)):
                            if isvariable[j] == True:
                                l += 1
                                cov_elem = pcov[k,l]
                                if k>=l:
                                    if cov_elem < 0:
                                        output_string+= (f'{cov_elem:<{paddingright}.2e}')
                                    else:
                                        paddingright_modified = paddingright - len(signspace)
                                        output_string+= (f'{signspace}{cov_elem:<{paddingright_modified}.2e}')
                                        
                output_string+= (f'\n\n')
                
                output_string+= (f'Correlation matrix\n')  
                paddingright = 8
                
                
                output_string+= (f'{emptystr:<{paddingright}}')
                
                for i in range(0,len(popt)):
                    if isvariable[i] == True:
                        output_string+= (f'{pnames[i]:<{paddingright}}')       
                k = -1        
                for i in range(0,len(popt)):
                    if isvariable[i] == True:
                        k += 1
                        l = -1
                        output_string+= (f'\n{pnames[i]:<{paddingright}}')
                        for j in range(0,len(popt)):
                            if isvariable[j] == True:
                                l += 1
                                cor_elem = pcor[k,l]
                                if k>=l:
                                    if cor_elem < 0:
                                        output_string+= (f'{cor_elem:<{paddingright}.2f}')
                                    else:
                                        paddingright_modified = paddingright - len(signspace)
                                        output_string+= (f'{signspace}{cor_elem:<{paddingright_modified}.2f}')

                                        
                output_string+= (f'\n\n')       
                

            
        if model_definition is not None:
            if model_definition['constrain_with_harmonics'] != []:
                
                emptystr = ''
                paddingright = 12
                signspace = ' '
                
                output_string+=(f'\nSeasonal harmonics\n\n')       
                output_string+=(f'\nIn observed heads\n')
                A = self.harmonic_observed.amplitude[0]
                sigma = self.harmonic_observed.amplitude[1]
                
                output_string+= (f'Amplitude {A:<{paddingright}.2e}')              
                output_string+= (f'stdev = {sigma:<{paddingright}.2e}\n')                
                
                phi = self.harmonic_observed.phase_shift[0]
                sigma = self.harmonic_observed.phase_shift[1]
                    
                output_string+= (f'Phase {phi:<{paddingright}.2e}')              
                output_string+= (f'stdev = {sigma:<{paddingright}.2e}\n')                 
                
                
                output_string+=(f'\n\nIn interpolated heads\n')
                A = self.harmonic_interpolated.amplitude[0]
                sigma = self.harmonic_interpolated.amplitude[1]
                
                output_string+= (f'Amplitude {A:<{paddingright}.2e}')              
                output_string+= (f'stdev = {sigma:<{paddingright}.2e}\n')                
                
                phi = self.harmonic_interpolated.phase_shift[0]
                sigma = self.harmonic_interpolated.phase_shift[1]
                    
                output_string+= (f'Phase shift {phi:<{paddingright}.2e}')              
                output_string+= (f'stdev = {sigma:<{paddingright}.2e}\n')            
                
                
                try :
                    output_string+=(f'\n\nIn modeled heads\n')
                    A = self.harmonic_modeled.amplitude[0]
                    sigma = self.harmonic_modeled.amplitude[1]
                    output_string+= (f'Amplitude {A:<{paddingright}.2e}')              
                    output_string+= (f'stdev = {sigma:<{paddingright}.2e}\n')   
                    phi = self.harmonic_modeled.phase_shift[0]
                    sigma = self.harmonic_modeled.phase_shift[1]
                    output_string+= (f'Phase shift {phi:<{paddingright}.2e}')              
                    output_string+= (f'stdev = {sigma:<{paddingright}.2e}\n')                     
                    
                except:
                    output_string+= (f'\nNo harmonic component could be fitted in modeled heads\n')              



                for stress_type in model_definition['constrain_with_harmonics']:
                    
                    if _stresses_dict is not None:
                        for e in _stresses_dict[stress_type]:
                            ts = _stresses_dict[stress_type][e]                    
    
                            output_string+=(f'\nIn observed {stress_type} {e}\n')
                            A = ts.harmonic_observed.amplitude[0]
                            sigma = ts.harmonic_observed.amplitude[1]
                            
                            output_string+= (f'Amplitude {A:<{paddingright}.2e}')              
                            output_string+= (f'stdev = {sigma:<{paddingright}.2e}\n')                
                            
                            phi = ts.harmonic_observed.phase_shift[0]
                            sigma = ts.harmonic_observed.phase_shift[1]
                                
                            output_string+= (f'Phase {phi:<{paddingright}.2e}')              
                            output_string+= (f'stdev = {sigma:<{paddingright}.2e}\n')                 
                            
                            
                            output_string+=(f'\n\nIn interpolated {stress_type} {e}\n')
                            A = ts.harmonic_interpolated.amplitude[0]
                            sigma = ts.harmonic_interpolated.amplitude[1]
                            
                            output_string+= (f'Amplitude {A:<{paddingright}.2e}')              
                            output_string+= (f'stdev = {sigma:<{paddingright}.2e}\n')                
                            
                            phi = ts.harmonic_interpolated.phase_shift[0]
                            sigma = ts.harmonic_interpolated.phase_shift[1]
                                
                            output_string+= (f'Phase shift {phi:<{paddingright}.2e}')              
                            output_string+= (f'stdev = {sigma:<{paddingright}.2e}\n')                       
                            
                            output_string+=(f'\n\nIn modeled response to {stress_type} {e}\n')
                            A = ts.harmonic_modeled.amplitude[0]
                            sigma = ts.harmonic_modeled.amplitude[1]
                            
                            output_string+= (f'Amplitude {A:<{paddingright}.2e}')              
                            output_string+= (f'stdev = {sigma:<{paddingright}.2e}\n')                
                            
                            phi = ts.harmonic_modeled.phase_shift[0]
                            sigma = ts.harmonic_modeled.phase_shift[1]
                                
                            output_string+= (f'Phase shift {phi:<{paddingright}.2e}')              
                            output_string+= (f'stdev = {sigma:<{paddingright}.2e}\n')                        
                            
                
        return output_string
    
    
    
                
                                  
if __name__ == "__main__":

        
    tminstr = '16-04-2023 08:00:00'
    tmaxstr = '26-04-2023 08:30:00'
    
    #curdir=os.getcwd()
    abs_path = os.path.dirname(os.path.abspath(__file__))
    abs_path_splitted = abs_path.split('\\')
    path_to_parent_folder_elements = abs_path_splitted[:-1]
    path_to_parent_folder = os.path.join(*path_to_parent_folder_elements)
    splitted=path_to_parent_folder.split(':')#trick  to  repair  path
    path_to_parent_folder = splitted[0]+':\\'+splitted[1]#trick  to  repair  path
    path_to_file_folder = os.path.join(path_to_parent_folder,'resources')
    
    # # DINO
    # path_to_file=os.path.join(path_to_file_folder,'PB101-1_20230508165127.csv')
    # heads =  MonitoringWell.read_from_DINO(path_to_file, tminstr = tminstr, tmaxstr = tmaxstr, settings = settings)
    # ax = heads.plot()    
    
    # # own program
    # path_to_file=os.path.join(path_to_file_folder,'PB101-1_20230508165127.csv')
    # heads =  MonitoringWell.read_from_DINO(path_to_file, tminstr = tminstr, tmaxstr = tmaxstr, settings = settings)
    # ax = heads.plot() 
    
    tminstr = '01-01-1900 08:00:00'
    tmaxstr = '31-12-2100 08:30:00'
    
    path_to_file = os.path.join(path_to_file_folder,'28AP0093_1.txt')
    heads = Heads.read_from_csv(path_to_file, tminstr = tminstr, tmaxstr = tmaxstr)
    ax = heads.plot()
    harmonic = heads.harmonic_observed

    ax = heads.plot_harmonics()
    
    
    
    output_string = heads.generate_output(model_definition = None)
    print('1210 output_string: ',output_string)
    
    print('1213 representer: ',repr(heads))
    