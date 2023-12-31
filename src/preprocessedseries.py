# -*- coding: utf-8 -*-
"""
Created on Sun Jul 23 15:00:20 2023

@author: christophe
"""


import numpy as np
import os
from scipy import interpolate
from datetime import datetime
from matplotlib.dates import date2num, num2date
from logging import getLogger
from stresstimeseries import Stress
from headstimeseries import Heads
from utilities import datestring2num,basename
from plotfunctions import simpleplot
from harmonics import Harmonics
from modeldefinition import ModelDefinition

logger = getLogger(__name__)


class Preprocessed:
    
    """
    A    class that provides instances of pre-processed time equidistant series
    needed to apply Fast Fourier Transform convolutions or neural networks modeling

    Attributes
    ----------
    

    tmin : float
        The minimum time of the preprocessed time series 
        (time number as defined in matplotlib.dates function date2num)

    tmax : float
        The maximum time of the preprocessed time series 
        (time number as defined in matplotlib.dates function date2num)    

    tminstr : string
        The minimum time of the preprocessed time series as string (dd-mm-yyyy HH:MM:SS)

    tmaxstr : string
        The maximum time of the preprocessed time series as string (dd-mm-yyyy HH:MM:SS)   

    time : numpy array_like
        Arrays of time numbers generated by matplotlib.dates function date2num
        when prepocessing the input time series using the preprocessedseries.py module
                
    stresses_dict : python object of data type 'dict'
        A dictionary containing all the entered
        Stress objects used to explain the observed groundwater levels variations.
        
    heads : Heads object or list of Heads objects
        Single Heads object or list of Heads objects containing the observed groundwater levels,
        associated metadata and preprocessed versions of the groundwater level time series 

    settings : python object of data type 'dict' (optional)
        An optional dictionary of settings - not used in present version

    Nint_dict : python object of data type 'dict'
        A dictionary specifying per stress the memory of each stress, expressed 
        in number of time intervals (a time equals a time step)
        
    memory_dict : python object of data type 'dict'
        A dictionary specifying per stress the memory of each stress, expressed 
        in time units.      
                
    time_step : float
        The specified time step used for all heads and stresses time series     
          
    time_step_targets : float
        The specified time step used for the targets  
              
    time_boundaries_sim : python object of data type 'dict'
        A dictionary specifying the time boundaries of the simulation          
        
    model_definition : Modeldefinition object
        A Modeldefinition object specifying how the time series analysis is implemented 
        (which is the first step in the present time series analysis methodology) 
        
    all_data_for_neural_networks : python object of data type 'dict'
        A dictionary containing the data needed as input for a neural network model
        
    all_targets_for_neural_networks : python object of data type 'dict'
        A dictionary containing the data needed as input for a neural network model
        
    forcast : boolean
        A boolean flag to indicate if preprocessing applies to forcast groundwater heads 
        (in which case the simulation period extends the observation period)
        
        

    Methods
    -------
    cumulate(series, tmin = None, tmax = None)
        A method used to cumulate the values of a time series
        
    generate_time_gaps_list(obs_time_stamps,threshold):
        A method used to detect time gaps larger than threshold
        
   generate_targets_selector(time,time_gaps):
       A method that generates a boolean array indicating which time stamp corresponds to a target
   
    preprocess()
        The method that generates the preprocessed time series
        
        

    """
    

    

    def __init__(self, heads = None, stresses_dict = None,time_step = None,time_step_targets = None, model_definition = None,
                 memory_dict = None, tminstr = None, tmaxstr = None, tminnum = None, tmaxnum = None, 
                 Nint = None, settings = None, forcast = False):
        
        
        clsname = str(self.__class__.__name__)
        modulename = str(__name__)        

        if tminnum is not None and tminstr is not None:
            logger.warning("Minimum time is both given as string and numeric value."
                              "String value will be ignored")
            
        if tminnum is None and tminstr is None:
            tminstr_default = '01-01-1900'
            dtmin = datetime.strptime(tminstr_default,'%d-%m-%Y')
            tmin = date2num(dtmin)
            
        if tminstr is not None:
            tmin = datestring2num(tminstr)

            
        if tminnum is not None:   # entering a time number overrides entering a time string
             tmin = tminnum
                         
             
        if tmaxnum is not None and tmaxstr is not None:
            message = (f'\nIn class {clsname} of module {modulename}.py: maximum time is both given '
                       f'as string and numeric value. String value will be ignored.\n')                          
            logger.warning(message)                      

        if tmaxnum is None and tmaxstr is None:
            tmaxstr_default = '31-12-2100'
            dtmax = datetime.strptime(tmaxstr_default,'%d-%m-%Y')
            tmax = date2num(dtmax)
             
        if tmaxstr is not None:
            tmax = datestring2num(tmaxstr)
            
        if tmaxnum is not None:   # entering a time number overrides entering a time string
             tmax = tmaxnum    
                           
        self.tmin = tmin
        self.tmax = tmax
        self.tminstr = tminstr
        self.tmaxstr = tmaxstr        
        self.settings = settings
        self.heads = heads
        self.stresses_dict = stresses_dict
        
        if time_step is None:
            time_step=1.0
        self.time_step = time_step   
        self.model_definition = model_definition  

        if time_step_targets is None:
            time_step_targets = 14.0
        self.time_step_targets = time_step_targets
        self.memory_dict = memory_dict

        Nint_dict = {} 
        
        """Nint is a dictionary of number of time intervals needed for ech stress time series to take into account the system memory"""
        for key in memory_dict:
            memory = memory_dict[key]
            Nint = self.memory_to_time_intervals(memory, time_step)
            Nint_dict[key] = Nint
        self.Nint_dict = Nint_dict

        self.time = None
        self.all_data_for_neural_networks = None
        self.all_targets_for_neural_networks = None        
        self.forcast = forcast

    def __repr__(self):
        """Prints a simple string representation of the time series."""
        return f"{self.__class__.__name__}" \
               f"(name={self.name}, " \
               f"tmin={self.settings['tmin']}, " \
               f"tmax={self.settings['tmax']})"
        
        
    default_settings = {}

    @staticmethod    
    def memory_to_time_intervals(memory, time_step):
        # return 1 + int(memory / time_step)
        return int(memory / time_step)
               
    @staticmethod    
    def interpolate(series, time = None, time_step = None, tmin = None, tmax = None):
        
        """
        Function to be used to interpolate observation series
        
        Parameters
        ----------
        series: numpy array
                array of dimension 2 x Number of time series records
                
        time: numpy array
                array of dimension 1 x Number of interpolation time stamps
                
    
        Returns
        -------
        interpolated: numpy array
                array of dimension 2 x Number of interpolation time stamps
        
        """
        if time is None:
            tmin_ = series[0,0]
            tmax_ = series[-1,0]            
            time = np.arange(tmin_, tmax_, time_step)

        interpolated = np.empty((len(time),2))
        interpolated[:,0] = time[:]
        
        #selected_window=(serie[:,0]>=tmin)&(serie[:,0]<tmax)
        #serie=serie[selected_window]
    
        interp_func=interpolate.interp1d(series[:,0],series[:,1],kind='linear',bounds_error=False)
        interpolated[:,1]=interp_func(time)

        return interpolated

    def cumulate(self,series, tmin = None, tmax = None):
        """
        Method to be used to generate a time series of cumulated observations values
        from time series of quantities that can be cumulated such as precipitation
        or evaporation
        
        Parameters
        ----------
        series: numpy array_like
                array of dimension 2 x Number of time series records

        tmin : float
            The minimum time of the preprocessed time series 
            (time number as defined in matplotlib.dates function date2num)
    
        tmax : float
            The maximum time of the preprocessed time series 
            (time number as defined in matplotlib.dates function date2num)  
        
        
        Returns
        -------
        cumulated: numpy array_like
                array of dimension 2 x Number of time series records - 1
        
        """
        cumul=np.empty(np.shape(series),dtype=float)
        cumul[:,:]=series[:,:]
        for i in range(1,len(cumul)):
            cumul[i,1]=cumul[i-1,1]+cumul[i,1]
        
        if tmin is not None: 
            mask = cumul[:,0] <= tmin
            cumul = cumul[mask,:]
            
        if tmax is not None: 
            mask = cumul[:,0] >= tmax
            cumul = cumul[mask,:]            
            
        return cumul

  


    def generate_time_gaps_list(self,obs_time_stamps,threshold):
    
        """
        A method to be used to identify gaps in observation time series
        
        
        Parameters
        ----------
        obs_time_stamps: numpy array_like
                array of observations time stamps
                
        threshold : float
            probability density at quantile value. For example: if threshold 
            is 0.75, then a gap is detected if time interval > 75% of 
            all time intervals
                
        Returns
        -------
        time_gaps: list
                list of time gaps

        """       
        time_intervals_tuples=zip(obs_time_stamps[:-1],obs_time_stamps[1:])
        time_intervals_values = []
        for time_interval in time_intervals_tuples:
            time_intervals_values.append(time_interval[1] - time_interval[0])

        
        a = np.asarray(time_intervals_values)
        q = np.quantile(a, threshold)
        time_gaps=[]
        
        for time_interval in time_intervals_tuples:
            begin = time_interval[0]
            end = time_interval[1]
            if begin - end > q:
                time_gaps.append(time_interval)
                


        return time_gaps   


    
    def generate_targets_selector(self,time,time_gaps):
    
        """
        Method to be used to select calibration targets
        
        Parameters
        ----------
        time : numpy array_like
        Arrays of time numbers generated by matplotlib.dates function date2num
        when prepocessing the input time series using the preprocessedseries.py module

        time_gaps: list
                list of time gaps
                
        Returns
        -------
        selector: numpy array_like
                Array of boolean values of dimension 2 x number of heads time 
                series records

        """
        
        selector = np.zeros(len(time),dtype=bool)
        delt_obs  =self.time_step
        delt_targets = self.time_step_targets
        
        skip = int(delt_targets/delt_obs)
        target_indexes = np.arange(0,len(time))
        target_indexes = target_indexes[::skip]
        selector[target_indexes] = True 
     
        for gap in time_gaps:
            mask=((time > gap[0]) & (time < gap[1]))
            selector[mask] = False
   

        return selector    
    
    
               
    def preprocess(self): 
              
        """
        A quite fundamental method in the present time series modeling code that is 
        used to generate equidistant time series allowing to apply 
        Fast Fourier Transform convolutions or neural networks modeling
        
        
        Parameters
        ----------
        see class definition
                
        Returns
        -------
        self: Preprocessed object
        Returns an instance of the Preprocessed class

        """        
        
        
        clsname = str(self.__class__.__name__)
        modulename = str(__name__)
        
        stresses_dict = self.stresses_dict
        time_step = self.time_step 
        time_step_targets = self.time_step_targets
        settings = self.settings
        memory_dict = self.memory_dict
        Nint_dict = self.Nint_dict
        model_definition = self.model_definition
        
        
        t = lambda t_str : date2num(datetime.strptime(t_str,'%Y-%m-%d'))

        list_tmin = [self.tmin] # initialize list of minimum times 

        list_tmax = [self.tmax] # initialize list of minimum times 
        

        if not isinstance(self.heads, list):
            _heads = [self.heads]
        else:
            _heads = self.heads
            
            
        if self.forcast == False: 
            for heads in _heads:
                list_tmin.append(heads.observed[0,0])
                list_tmax.append(heads.observed[-1,0])   
                

        _observed_stresses_dict = {}
        """ make a local copy of stresses dictionary
        """
        for key in self.stresses_dict:
            _observed_stresses_dict[key] = {}
            for e in self.stresses_dict[key]:
                _observed_stresses_dict[key][e] = self.stresses_dict[key][e].observed
        

        if model_definition['root_zone']['apply_root_zone'] == True:
            """ this condition make sure that prec and evap have similar
            begin and end time when applying the root zone module
            
            """
            sublist_tmin = []
            sublist_tmax = []
            for key in ['prec','evap']:
                if key in stresses_dict:
                    for e in self.stresses_dict[key]:
                        sublist_tmin.append(_observed_stresses_dict[key][e][0,0])
                        sublist_tmax.append(_observed_stresses_dict[key][e][-1,0])
            subtmin = int(max(list_tmin)+1)
            subtmax = int(min(list_tmax)-1)    
            
            for key in ['prec','evap']:
                if key in  _observed_stresses_dict:
                    for e in _observed_stresses_dict[key]:
                        mask = (_observed_stresses_dict[key][e][:,0] > subtmin) & (_observed_stresses_dict[key][e][:,0] < subtmax)
                        _observed_stresses_dict[key][e] = _observed_stresses_dict[key][e][mask,:]

        for key in  _observed_stresses_dict:
            for e in  _observed_stresses_dict[key]:
                observed_stress =  _observed_stresses_dict[key][e]
                Nint = self.Nint_dict[key]
                
                warmup_period = Nint*time_step
                list_tmin.append( observed_stress[0,0]+warmup_period)
                list_tmax.append( observed_stress[-1,0])

        time = None
        # tmin = int(max(list_tmin)+1)
        # tmax = int(min(list_tmax)-1)
        tmin = max(list_tmin)
        tmax = min(list_tmax) 

        
        self.tmin = tmin
        self.tmax = tmax

        if tmax > tmin:

            time = np.arange(tmin,tmax,time_step)
            self.time = time
            
            for key in self.stresses_dict:
                for e in self.stresses_dict[key]:
                    stress = self.stresses_dict[key][e]
                    Nint = Nint_dict[key]
                    _time = np.arange(tmin-(Nint-1)*time_step,tmax,time_step) 
                    interpolated = np.empty((len(_time),2))
                    interpolated[:,0] = _time
                    _tmin = interpolated[0,0]
                    _tmax = interpolated[-1,0]
                    

                    condition = False
                    if key in ['prec','evap']: 
                        condition = True
                        
                    if key == 'pump':
                        if model_definition['pump']['given_as_pumping_rate'] == False:       
                            condition = True
  
                        else:
                            tseries = stress.observed
                            cumulated_volume = np.empty((len(tseries),2),dtype=float)
                            
                            
                            cumulated_volume[:,0] = tseries[:,0]
                            cumulated_volume[0,0] = 0

                            for i in range(1,len(tseries)):
                                cumulated_volume[i,1] = cumulated_volume[i-1,1] + tseries[i-1,1] * (tseries[i,0] - tseries[i-1,0]) # assume pumping rate time unit is the same as x axix time unit (days)
        
                            # if tmin is not None: 
                            #     mask = cumulated_volume[:,0] <= tmin
                            #     cumulated_volume = cumulated_volume[mask,:]
                                
                            # if tmax is not None: 
                            #     mask = cumulated_volume[:,0] >= tmax
                            #     cumulated_volume = cumulated_volume[mask,:]      
                        
                            interp_func = interpolate.interp1d(cumulated_volume[:,0],cumulated_volume[:,1],kind='linear',bounds_error=False)
                            interp_values = interp_func(_time) 
                            differentiated = interp_values[1:] - interp_values[:-1]
                            
                            interpolated[1:,1] = differentiated[:] / self.time_step   
                            interp_values = interp_func(_time) 
                            
                            
                            interp_func_nearest_neighbour = interpolate.interp1d(tseries[:,0],tseries[:,1],kind='nearest',bounds_error=False)
                            interp_values_nearest_neighbour = interp_func_nearest_neighbour(_time) 

                            
                            interpolated[0,1] = interp_values_nearest_neighbour[0] 
      

                    if condition == True:    
                        cumul = self.cumulate(stress.observed)
                        mask = cumul[:,0] >= _tmax                            
                        try:
                            after = cumul[mask,:]
                            lastvalue = after[0,1] - after[1,1]
                            delt = after[0,0] - after[1,0]
                        except: 
                            lastvalue = np.mean(stress.observed[:,1])
                            message = (f'\nIn class {clsname} of module {modulename}.py: '
                                       f'last interpolation value is set to mean observed.\n')                          
                            logger.warning(message)                            
                                                

                        interp_func = interpolate.interp1d(cumul[:,0],cumul[:,1],kind='linear',bounds_error=False)
                        interp_values = interp_func(_time) 
                        differentiated = interp_values[1:] - interp_values[:-1]
                        interpolated[-1,1] = lastvalue
                        interpolated[:-1,1] = differentiated[:]


                    else: 
                        interp_func = interpolate.interp1d(stress.observed[:,0],stress.observed[:,1],kind='linear',bounds_error=False)
                        interpolated[:,1] = interp_func(_time) 

                    self.stresses_dict[key][e].interpolated = interpolated

                    interpolated_normalized = np.empty(np.shape(interpolated))
                    interpolated_normalized[:,:] = interpolated[:,:]
                    mean = np.mean(interpolated_normalized[:,1])
                    
                    interpolated_normalized[:,1] -= mean                    
                    self.stresses_dict[key][e].interpolated_normalized = interpolated_normalized

                    harmonic_input = np.empty(np.shape(interpolated_normalized))
                    harmonic_input[:,:] = interpolated_normalized[:,:]
                        
                    message = (f'\nIn class {clsname} of module {modulename}.py: '
                               f'harmonic is multiplied by -1 for evap and well.\n')                          
                    logger.warning(message)                    
      
                    
                    if key in ['evap','pump']:
                        harmonic_input[:,1] *= -1

                    self.stresses_dict[key][e].harmonic_interpolated = Harmonics(harmonic_input, tmin = tmin, tmax = tmax, tseries_type = key).fit_harmonic()  
                      
                    if key in ['prec','evap']:

                        message = (f'\nIn class {clsname} of module {modulename}.py: '
                                   f'harmonic of interpolated is divided by time_step for prec, evap and well.\n')                          
                        logger.warning(message)
                        
                        self.stresses_dict[key][e].harmonic_interpolated.harmonic_component[:,1] /= time_step
                        self.stresses_dict[key][e].harmonic_interpolated.amplitude[0] /= time_step
                        self.stresses_dict[key][e].harmonic_interpolated.amplitude[1] /= time_step      
                        
                    if key in ['pump']:
                        if model_definition['pump']['given_as_pumping_rate'] == False: #
                            message = (f'\nIn class {clsname} of module {modulename}.py: '
                                       f'harmonic of interpolated is divided by time_step for prec, evap and well.\n')                          
                            logger.warning(message)
                            
                            self.stresses_dict[key][e].harmonic_interpolated.harmonic_component[:,1] /= time_step
                            self.stresses_dict[key][e].harmonic_interpolated.amplitude[0] /= time_step
                            self.stresses_dict[key][e].harmonic_interpolated.amplitude[1] /= time_step                           
        
                        
            preprocessed_heads = []
            
            for heads in _heads:
                
                heads.stresses_dict = self.stresses_dict
            
                interpolated = self.interpolate(heads.observed,time)
                heads.interpolated = interpolated

                interpolated_normalized = np.empty(np.shape(interpolated))
                interpolated_normalized[:,:] = interpolated[:,:]
                mean = np.mean(interpolated[:,1])
                interpolated_normalized[:,1] -= mean                    
                heads.interpolated_normalized = interpolated_normalized
                heads.harmonic_interpolated = Harmonics(interpolated_normalized, tmin = tmin, tmax = tmax).fit_harmonic() 
            
                obs_time_stamps = heads.observed[:,0]
                threshold = 0.9 # probability density at quantile value (so here a gap is detected if time interval > 90% of time intervals )
                time_gaps = self.generate_time_gaps_list(obs_time_stamps,threshold)
                heads.targets_selector = self.generate_targets_selector(time,time_gaps)
                heads.targets = heads.interpolated[heads.targets_selector]            

                # the weights vector below is given as default, more sophisticated should come later
                heads.weights = np.ones(len(heads.targets))
   
                preprocessed_heads.append(heads)
        

        else:
            message = (f'\nIn class {clsname} of module {modulename}.py: No valid time array was entered,'
                       f'check time overlap over all time series.\n')     
            logger.error(message)     
        
        # prepare inputs for neural networks training
        
        if not isinstance (self.heads,list):
            dim = 0
            for key in self.stresses_dict:
                for e in self.stresses_dict[key]:
                    dim += Nint_dict[key]
                        
            all_data_for_neural_networks = np.empty((len(time),dim))
            all_targets_for_neural_networks = np.empty(len(time))                    
                        
            
            
            # for i in range(0,len(self.heads.targets)):
            #     all_targets[i] = self.heads.targets[i,1]
            #     startindex = 0
            #     stopindex=0        
            
                # self.heads.targets_selector = self.generate_targets_selector(time,time_gaps)
                # self.heads.targets = self.heads.interpolated[self.heads.targets_selector]
            
            
            for i in range(0,len(time)):
                all_targets_for_neural_networks[i] = self.heads.interpolated[i,1]   
                startindex = 0
                stopindex=0
    
    
                for key in self.stresses_dict:
                    for e in self.stresses_dict[key]:
                        stress_data = self.stresses_dict[key][e].interpolated[:,1]
                        mean=np.mean(stress_data)
                        normalized_stress_data=stress_data-mean
                        Nint = Nint_dict[key]
                        stopindex = startindex+Nint
                        j = i+Nint
                        all_data_for_neural_networks[i,startindex:stopindex] = normalized_stress_data[j-Nint:j]
                        # all_data_for_neural_networks[i,startindex:stopindex] = stress_data[j-Nint:j]
                        startindex = stopindex

            #normalize targets
            mean = np.mean(all_targets_for_neural_networks)
            all_targets_for_neural_networks -= mean
            
        
        self.time = time
        self.stresses_dict = stresses_dict
        if isinstance(self.heads,list):
            self.heads = preprocessed_heads
            self.all_data_for_neural_networks = None
            self.all_targets_for_neural_networks = None
        else:
            self.heads = preprocessed_heads[0] # if dealing with one head series only
            self.all_data_for_neural_networks = all_data_for_neural_networks
            self.all_targets_for_neural_networks = all_targets_for_neural_networks  
            
        
            

        return self #time, heads, stresses_dict, Nint_dict, all_data_for_neural_networks, all_targets_for_neural_networks
    
    
if __name__ == "__main__":

    #curdir=os.getcwd()
    model_definition = ModelDefinition().model_definition
    abs_path = os.path.dirname(os.path.abspath(__file__))
    abs_path_splitted = abs_path.split('\\')
    path_to_parent_folder_elements = abs_path_splitted[:-1]
    path_to_parent_folder = os.path.join(*path_to_parent_folder_elements)
    splitted = path_to_parent_folder.split(':')#trick  to  repair  path
    path_to_parent_folder = splitted[0]+':\\'+splitted[1]#trick  to  repair  path
    path_to_file_folder = os.path.join(path_to_parent_folder,'resources')  

    stresses_dict = {}
    
    path_to_file = os.path.join(path_to_file_folder,'260_De_Bilt_PREC_19010101_20200910.csv')
    key = basename(path_to_file)
    stresses_dict['prec'] = {}
    stresses_dict['prec'][key] = Stress.read_from_csv(path_to_file, stress_type = 'prec')    
    

    path_to_file=os.path.join(path_to_file_folder,'260_De_Bilt_EVAP_19010101_20200910.csv')
    key = basename(path_to_file)
    stresses_dict['evap'] = {}
    stresses_dict['evap'][key] = Stress.read_from_csv(path_to_file, stress_type = 'evap')    


    tminstr = '01-01-1900 08:00:00'
    tmaxstr = '31-12-2100 08:30:00'
    
    path_to_file_folder = os.path.join(path_to_parent_folder,'resources')  
    path_to_file = os.path.join(path_to_file_folder,'28AP0093_1.txt')
    heads = Heads.read_from_csv(path_to_file, tminstr = tminstr, tmaxstr = tmaxstr)


    arguments_dic = {}

    for key in stresses_dict:
        for e in stresses_dict[key]:
            stresses_dict[key][e].plot()

    
        
    time_step = 1.#/(24.*2)
    
    memory_dict = {'prec':365*5,
              'evap':365*5,
              'riv':7,
              'pump':365*5,
              'noise':30}      
    
    preprocessed = Preprocessed(heads = heads,stresses_dict = stresses_dict, time_step = time_step, 
                                model_definition = model_definition, memory_dict = memory_dict,
                                tminstr = tminstr, tmaxstr = tmaxstr).preprocess()
    
    time = preprocessed.time
    heads = preprocessed.heads
    print('639 heads.harmonic_interpolated.amplitude',heads.harmonic_interpolated.amplitude)

    stresses_dict = preprocessed.stresses_dict 

    Nint_dict = preprocessed.Nint_dict
    all_data_for_neural_networks = preprocessed.all_data_for_neural_networks
    all_targets_for_neural_networks = preprocessed.all_targets_for_neural_networks

    
    heads.plot()
    ax = heads.plot_harmonics()
    ampl_h_obs = heads.harmonic_observed.amplitude[0]
    phi_h_obs = heads.harmonic_observed.phase_shift[0]
    print('492 Ampl_h_obs',ampl_h_obs)
    print('493 phi_h_obs',phi_h_obs)    

    
    ampl_h_interp = heads.harmonic_interpolated.amplitude[0]
    phi_h_interp = heads.harmonic_interpolated.phase_shift[0]
    print('497 Ampl_h_inter',ampl_h_interp)
    
    
    for stress_type in stresses_dict:
        for series_name in stresses_dict[stress_type]:
            ts = stresses_dict[stress_type][series_name]
            ts.plot()
            ax = ts.plot_harmonics()
            
            print('507 stress_type',stress_type)
            interpolated = ts.interpolated
            print('670 interpolated',interpolated)
            
            ampl_obs = ts.harmonic_observed.amplitude[0]
            phi_obs = ts.harmonic_observed.phase_shift[0]
            print('674 Ampl_obs',ampl_obs)
            print('675 phi_obs',phi_obs)    
            
            ampl_h_interp = ts.harmonic_interpolated.amplitude[0]
            phi_h_interp = ts.harmonic_interpolated.phase_shift[0]
            print('678 Ampl_inter',ampl_h_interp)
            print('679 phi_inter',phi_h_interp)
