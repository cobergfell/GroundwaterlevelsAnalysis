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
# from abc import ABC, abstractmethod
# from collections import OrderedDict
from matplotlib import pyplot as plt
from matplotlib.pyplot import figure, show,savefig,clabel,contour,setp,gcf,getp,gca,close
from utilities import orderOfMagnitude,generate_ticks,datestring2num,stdev_back_from_log
from plotfunctions import simpleplot
from wellmetadata import WellMetaData
from bs4 import BeautifulSoup   
import requests
        
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
        
        metadata = WellMetaData(name = name, X = X, Y = Y, Z = Z, surface_level = surface_level, 
                        screen_top = screen_top, screen_bottom = screen_bottom )  

    
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
    def read_from_DINO(cls,filepath, tminstr = None, tmaxstr = None, separator = ',', settings = None, datetimeformat = '%d-%m-%Y', artdiverformat = False):
        
        """
        A class method to read groundwater head time series from am ascii file 
        in the DINOloket format (from the Dutch national geological survey institute)
        and generate a Heads object.

        
        Parameters
        ----------
        filepath: string
            Path to the csv file

        tminstr: string (optional)
                String specifying the minimum time on the plot (default format dd-mm-yyyy)
                
        tmaxstr: string (optional)
                String specifying the maximum time on the plot (default format dd-mm-yyyy)    
                
        settings : python object of data type 'dict' (optional)
            An optional dictionary of settings       
            
            
        datetimeformat: string (optional)
                String specifying the date-time format (default format dd-mm-yyyy)       
                
                
        artdiverformat: bool
            Boolean variable indicating if the DINO format is the Artdiver DINO format or not         

        
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
            An instance of the Heads class

        """
        
        # from wellmetadata import WellMetaData
        from datetime import datetime
        from matplotlib.dates import date2num, num2date
        from datetime import datetime
        from numpy import median
        
        
        settings = settings
        tmin = -1e9
        tmax = +1e9
        
        valueindex = 5
        if artdiverformat == True:
            datetimeformat = '%Y/%m/%d %H:%M:%S'
            valueindex = 3
   
        if tminstr is not None:
            tmin = datestring2num(tminstr, datetimeformat = datetimeformat)
                        
        if tmaxstr is not None:
            tmax = datestring2num(tmaxstr, datetimeformat = datetimeformat)               

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
        data_bloc_begin = indices[-1]+1
        
        
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
            
            
            metadata = WellMetaData(name = name, X = X, Y = Y, Z = Z, surface_level = surface_level, 
                                    screen_top = screen_top, screen_bottom = screen_bottom )                   

            list_of_metadata.append(metadata)
            

                
        elif reference=='maaiveld':
            #to do check if metadata are different if reference is surface level (in Dutch: 'maaiveld')instead of absolute datum 
            pass
        
        
        metadata = list_of_metadata[-1]#the last metadata dictionaryis used later as reference when meta data are needed

        heads = np.empty((Nrecords,2))
        c = -1
        for i in range(data_bloc_begin,data_bloc_end):
            c = c+1

            line = all_lines[i]
            d = line.split(',')
            
            try:
                time_num = datestring2num(d[2], datetimeformat = datetimeformat)

            except :
                message = (f'\nDatetime format of {d[2]} could not be parsed\n')
                logger.warning(message)
                
            try:
                value = float(d[valueindex])/100 #divide by hundred because values in DINO are given in cm

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
    
    


    @classmethod      
    def read_from_BRO_GLD(cls,filepath, tminstr = None, tmaxstr = None, separator = ',', settings = None, 
                          datetimeformat = '%d-%m-%Y'):
        
        """
        A class method to read 'GLD' type of groundwater head time series in xml format such as 
        issued by the by broservices  https://publiek.broservices.nl 
        (from the Dutch national geological survey institute)
        and generate a Heads object.
        
        Note that in these files, the measurements data are chronologcal within
        <measurement></measurement> elements, but the succession of these measurements
        element is not chronological, so sorting the measurements elements in
        ascending order is necessary

        
        Parameters
        ----------
        filepath: string
            Path to the csv file

        tminstr: string (optional)
                String specifying the minimum time on the plot (default format dd-mm-yyyy)
                
        tmaxstr: string (optional)
                String specifying the maximum time on the plot (default format dd-mm-yyyy)    
                
        settings : python object of data type 'dict' (optional)
            An optional dictionary of settings       
            
            
        datetimeformat: string (optional)
                String specifying the date-time format (default format dd-mm-yyyy)       
                 

        
    
        Returns
        -------
        cls: Heads object
            An instance of the Heads class
            
            
        References
        ---------
        
        Parsing using beautifullsoup is well explained in 
        https://realpython.com/python-xml-parser/            

        """


        settings = settings
        tmin = -1e9
        tmax = +1e9
        

        if tminstr is not None:
            tmin = datestring2num(tminstr, datetimeformat = datetimeformat)
                        
        if tmaxstr is not None:
            tmax = datestring2num(tmaxstr, datetimeformat = datetimeformat)   
              


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
        # print(gmwsoup.prettify())
    
        
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
            
        heads = np.array(records_list)
        
        mask = heads[:,1] < 1e9
        heads = heads[mask]
        
        mask = (heads[:,0] > tmin) & (heads[:,0] < tmax)
        heads = heads[mask]
        
        
    
            
        #export to my own format for time series analysis
        
        Xstr = str('%8.2f' %X)
        Ystr = str('%8.2f' %Y)
        
        
        begin_date = num2date(heads[0,0])
        begin_date_string= str(begin_date.day)+'-'+str(begin_date.month)+'-'+str(begin_date.year)
        end_date=num2date(heads[-1,0])
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
        for i in range(0,len(heads)):
            t = num2date(heads[i,0])
            value = heads[i,1]
            record=str(t.day)+'-'+str(t.month)+'-'+str(t.year)+'  '+str(t.hour).zfill(2)+':'+str(t.minute).zfill(2)+':'+str(t.second).zfill(2)+','+str('%8.4f' %value)
            F.write(record+'\n')
        F.close()

        return cls(name = metadata.name, X = X, Y = Y, Z = Z, metadata = metadata, 
                   observed = heads, tminstr = tminstr, tmaxstr = tmaxstr, tminnum = tmin, tmaxnum = tmax, settings = settings) 


    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    def generate_gxg(self,tminstr= None, tmaxstr = None, alpha = 0.05):
        
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
    
    
    def plot(self, tminstr = None, tmaxstr = None, save_plot = False):
        
        """ 
        A method to plot the head time series (override the plot method of TimeSeriesBase)
        
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
        
        clsname = str(self.__class__.__name__)
        modulename = str(__name__)   
        
        observed = self._observed
        left_margin = 0.12
        right_margin = 0.06
        top_margin = 0.12
        bottom_margin = 0.44
        dy = 0.04
        dx = 0.05        

        fs = 12
     
        
        if tminstr is not None:
            tmin = datestring2num(tminstr)
        else:
            tmin = observed[0,0]
                        
        if tmaxstr is not None:
            tmax = datestring2num(tmaxstr)
        else:
            tmax = observed[-1,0]

        mask = (observed[:,0] >= tmin) & (observed[:,0] < tmax)
        observed = observed[mask]        


        yas_title = 'm above datum'

        plot_length = 1-left_margin-right_margin
        plot_height = 1-bottom_margin-top_margin
        xll = left_margin # x coordinate lower left corner
        yll = 1-top_margin-plot_height # y coordinate lower left corner
        
        
        fig=figure(num=None, figsize=(8,10),dpi=50,facecolor='w', edgecolor='k')# figsize=(width, height) in inches.

        ax = fig.add_axes([xll,yll,plot_length,plot_height],frameon=True, xscale=None, yscale=None)  
        
        leg_list = []
        ax.plot_date(observed[:,0],observed[:,1],'b-',linewidth = 1.0)     
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
            
        
        metadata = self.metadata
        
        X = metadata.X
        Y = metadata.Y
        Z = metadata.Z
        surface_level = metadata.surface_level
        screen_top = metadata.screen_top
        screen_bottom = metadata.screen_bottom
        
        time = observed[:,0]
        lower_expect, median, upper_expect = self.generate_gxg(tminstr = tminstr, tmaxstr = tmaxstr, alpha = 0.05)
        
        surface_level_array = surface_level*np.ones(len(observed))
        
        upper_expect_array = upper_expect*np.ones(len(observed))
        median_array = median*np.ones(len(observed))
        lower_expect_array = lower_expect*np.ones(len(observed))
        
        
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
    
        title='Time series of groundwater heads for piezometer '+ self.name
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
        self.name,
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
            if self.name is not None:
                figname = "plot_of_{}".format(self.name)
            else:
                figname = "time_series_plot"
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
                
    def plot_multiple_series(self, tseries_list = None, legend_list = None, share_axes = True, plot_title = None, yas_title = None, 
             tminstr = None, tmaxstr = None, save_plot = False):
        
        """ 
        A method to plot the time series.
        
        Parameters
        ----------    
        tseries_list: list (optional)
                list of time series to plot. The time series are given in the 
                form of a 2 dimensional numpy array of float, with time numbers
                given in the first array axis (time numbers as defined in 
                matplotlib.dates function date2num), and observed values 
                given in the second array axis
                
        tseries_list_names: list (optional)
                list of time series names given as strings

        
        yas_title: string (optional)
                string specifying Y-axis title   
                
        plot_title: string (optional)
                string specifying plot title 
                
                
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
        
        clsname = str(self.__class__.__name__)
        modulename = str(__name__)   
        
        observed = self._observed
        interpolated = self._interpolated
        tseries_type = self.tseries_type

        left_margin=0.12
        right_margin=0.06
        top_margin=0.12
        bottom_margin=0.14
        dy=0.14
        dx=0.05        

        fs = 14
        colors = ['b-','m-','y-','c-','g-', 'tab:green-','tab:olive-', 'tab:purple-','tab:pink-','tab:orange-','tab:brown-']
            
        
        if tminstr is not None:
            tmin = datestring2num(tminstr)
        else:
            tmin = observed[0,0]
                        
        if tmaxstr is not None:
            tmax = datestring2num(tmaxstr)
        else:
            tmax = observed[-1,0]

        mask = (observed[:,0] >= tmin) & (observed[:,0] < tmax)
        observed = observed[mask]        

        if tseries_list is None:
            tseries_list = [observed]
            legend_item = self.name+' '+'observed'
            legend_list = [legend_item]
            
            
        if legend_list is None:
            default_list = []
            for i in range(0,len(tseries_list )):
                default_list.append('series_'+str(i))
            
            legend_list = default_list


        nl = len(tseries_list)
        if nl > 4:
            message = (f'\nIn class {clsname} of module {modulename}.py: '
                       f'The list of time series to plot should not exceed 4.\n')                          
            logger.warning(message) 
            
    
        if plot_title is None:  
            # title = (f'time series type: {(self.tseries_type} name: {self.name}') 
            plot_title = f'{self.name}'

              
        if yas_title is None:  
            yas_title = 'm above datum'


        if share_axes == False:
            nc = 1.0
            plot_length = (1-left_margin-right_margin-dx)/nc
            plot_height = (1-bottom_margin-top_margin-2*dy)/nl
            x1 = left_margin
            x2 = x1+plot_length+dx
            X = [x1]#list of x coordinates of bottom left corners
            y1 = 1-top_margin-plot_height
            
            Y=[y1]#list of x coordinates of bottom left corners  
            if nl>1:
                for i in range(1,nl):
                    Y.append(Y[i-1]-plot_height-dy)
            
            fig=figure(num=None, figsize=(8,8*nl),dpi=50,facecolor='w', edgecolor='k')# figsize=(width, height) in inches.

            for i in range(0,nl): 
                ax = fig.add_axes([X[0],Y[i],plot_length,plot_height],frameon=True, xscale=None, yscale=None)  
                ax.plot_date(tseries_list[i][:,0],tseries_list[i][:,1],colors[i],linewidth = 1.0)     
           
            ax.set_title(plot_title,fontsize=fs)
            (x, y) = ax.title.get_position()
            #ax.title.set_y(0.95 * y)    
            ax.grid(True)         
          
                
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
                
                
        
            ax.annotate('Legend', xy=(0.05, -0.14), xycoords='axes fraction',horizontalalignment='left', verticalalignment='center',fontsize=fs)#this is for 1 column 2 rows
            leg = ax.legend((legend_list), ncol=2,shadow=False,bbox_to_anchor=[0.045, -0.27], loc='lower left',frameon=False)
            ltext  = leg.get_texts() 
            plt.setp(ltext, fontsize=fs)    # the legend text fontsize
            
            show()
            
        else:
            plot_length = (1-left_margin-right_margin-dx)
            plot_height = (1-bottom_margin-top_margin-2*dy)
            x_lower_left = left_margin
            y_lower_left = 1-top_margin-plot_height

            fig=figure(num=None, figsize=(8,8),dpi=50,facecolor='w', edgecolor='k')# figsize=(width, height) in inches.
            ax = fig.add_axes([x_lower_left,y_lower_left,plot_length,plot_height],frameon=True, xscale=None, yscale=None)  

            for i in range(0,nl): 
                ax.plot_date(tseries_list[i][:,0],tseries_list[i][:,1],colors[i],linewidth = 1.0)     
    
            ax.set_title(plot_title,fontsize=fs)
            (x, y) = ax.title.get_position()
            #ax.title.set_y(0.95 * y)    
            ax.grid(True)         
        
        
    
                
            def format_date(dates, pos=None):
                return num2date(dates).strftime('%Y')
            
            def format_date(dates, pos=None):
                return num2date(dates).strftime('%d-%m-%Y %H:%M')    
            
            def format_date(dates, pos=None):
                return num2date(dates).strftime('%d-%m-%Y')             
            
        
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
                
                
        
            ax.annotate('Legend', xy=(0.05, -0.14), xycoords='axes fraction',horizontalalignment='left', verticalalignment='center',fontsize=fs)#this is for 1 column 2 rows
            leg = ax.legend((legend_list), ncol=nl,shadow=False,bbox_to_anchor=[0.045, -0.27], loc='lower left',frameon=False)
            ltext  = leg.get_texts() 
            plt.setp(ltext, fontsize=fs)    # the legend text fontsize
            
            show()            
            
            
            
            
            
            
            
            
        
        if save_plot == True:
            if plot_title is not None:
                figname = "plot_of_{}".format(plot_title)
            else:
                figname = "time_series_plot"
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
    path_to_resources_folder = os.path.join(path_to_parent_folder,'resources')
    
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
    
    path_to_file = os.path.join(path_to_resources_folder,'28AP0093_1.txt')
    heads = Heads.read_from_csv(path_to_file, tminstr = tminstr, tmaxstr = tmaxstr)
    ax = heads.plot()
    harmonic = heads.harmonic_observed
    
    output_string = heads.generate_output(model_definition = None)
    print('1852 output_string: ',output_string)    
    
    path_to_file = os.path.join(path_to_resources_folder,'HBpb001-1.csv')
    path_to_file = os.path.join(path_to_resources_folder,'HBpb003-1.csv')
    path_to_file = os.path.join(path_to_resources_folder,'HBpb004.csv')
    heads = Heads.read_from_DINO(path_to_file, tminstr = tminstr, tmaxstr = tmaxstr, artdiverformat = True )

    ax = heads.plot()
    ax = heads.plot_harmonics()

    
    output_string = heads.generate_output(model_definition = None)
    print('1868 output_string: ',output_string)
    
    print('1870 representer: ',repr(heads))
    
    

    path_to_file = os.path.join(path_to_resources_folder,'GLD000000014660_IMBRO.xml')
    heads = Heads.read_from_BRO_GLD(path_to_file, tminstr = tminstr, tmaxstr = tmaxstr)

    ax = heads.plot()
    ax = heads.plot_harmonics()

    
    output_string = heads.generate_output(model_definition = None)
    print('1868 output_string: ',output_string)    
    
    
    
    