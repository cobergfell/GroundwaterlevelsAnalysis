
"""
Created on Fri Jun 16 22:30:13 2023

@author: christophe
"""

import numpy as np
from collections import OrderedDict
import os
import scipy
import numpy as np
from scipy.optimize import least_squares, leastsq
from datetime import datetime
from matplotlib.dates import date2num, num2date
from logging import getLogger
from modelnoise import ResidualsDecayExponentially
from utilities import *
from matplotlib.pyplot import figure, show,savefig,clabel,contour,setp,gcf,getp,gca,close



    
class RootzoneFunctionsParent:


    """
    A parent class for classes that provide instances of root zone functions,
    aimed at simulating the percolation time of precipitation through the 
    unsaturated zone.
    
    Note:
    At the time of writing this doc (6-9-2023), this class is still not ready,
    and documentation might need some improvements/updates
    
    ...
    
    Attributes
    ----------
        
    time : numpy array_like
        An arrays of time numbers generated by matplotlib.dates function date2num
        when prepocessing the input time series using the preprocessedseries.py module   
        
    _p_dict : python object of data type 'dict'
        _p_dict is a local copy of the p_dict dictionary defining the model parameters
        and associated specifications such as associated indexes, names, logtransformed or not, 
        variable or fixed etc. See the the parameterslogistic.py module for additional details.
        
    _stresses_dict : python object of data type 'dict'
        _stresses_dict is a local copy of the stresses_dict dictionary containing all the entered
        Stress objects used to explain the observed groundwater levels variations.
        
    time_step : float (optional)
        The specified time step used for all heads and stresses time series          
        
    tmin : float (optional)
        The minimum time of the preprocessed time series 
        (time number as defined in matplotlib.dates function date2num)

    tmax : float
        The maximum time of the preprocessed time series 
        (time number as defined in matplotlib.dates function date2num)  

    tminstr : string (optional)
        The minimum time of the preprocessed time series as string (dd-mm-yyyy HH:MM:SS)

    tmaxstr : string (optional)
        The maximum time of the preprocessed time series as string (dd-mm-yyyy HH:MM:SS)          
        
    regime: string (optional)
        Specify the fluctuation regime, when applicable (default is an empty string)

    root_zone_model: string (optional)
        Specify the root zone model (default is Brooks_and_Core)        
        
    
    Methods
    -------
    
    get_initial_parameters()
        Method that returns the initial parameters
        
        
    estimate_percolation()
        Method that returns a time series of percolated precipitations quantities.    
    
    """        
    
    def __init__(self, time = None, p_dict = None, stress_dict = None, time_step = None, 
                 tminstr = None, tmaxstr = None, tminnum = None, tmaxnum = None,
                 regime = 'regime_1', root_zone_model = 'Brooks_and_Corey'):
    
        
        self.regime = regime
        self.name = root_zone_model
        self.root_zone_model = root_zone_model
        self._p_dict = p_dict
        self._stress_dict = stress_dict
        self.time = time
        self.time_step = time_step
        self.water_saturation = None # time series of degree of water saturation
        self.percolation_per_time_step = None # time series of height of percolated water per time step
        self.percolation = None # time series of percolation rate dimension L.T-1
        self.actual_evaporation = None # time series of actual evaporation rate dimension L.T-1
        self.cumulated_actual_evaporation = None # time series of cumulated actual evaporation
        self.cumulated_potential_evaporation = None # time series of cumulated potential evaporation
        self.cumulated_precipitation = None # time series of cumulated precipitation
        self.cumulated_runoff = None # time series of cumulated run off
        self.cumulated_percolation = None  # time series of cumulated percolation      
        
        
        t = lambda t_str : date2num(datetime.strptime(t_str,'%d-%m-%Y')) 

        tmin = None
        tmax = None
        
        if time is not None:
            tmin = time[0]
            tmax = time[-1]
        
        if tminnum is not None and tminstr is not None:
            self.logger.warning("Minimum time is both given as string and numeric value."
                              "String value will be ignored")
            
        if tmin is None and tminnum is None and tminstr is None:
            tmin = t('01-01-1900')

   
        if tminstr is not None:
            tmin = datestring2num(tminstr)
            
        if tminnum is not None:   
             tmin = tminnum
                         
             
        if tmaxnum is not None and tmaxstr is not None:
            self.logger.warning("Maximum time is both given as string and numeric value."
                              "String value will be ignored")
            
        if tmax is None and tmaxnum is None and tmaxstr is None:
            tmax = t('31-12-2100')
                        
        if tmaxstr is not None:
            tmax = datestring2num(tmaxstr)
         
        if tmaxnum is not None:   
             tmax = tmaxnum                  
            
        
        self.tmin = tmin
        self.tmax = tmax
        
        if time is None:
            time = np.arange(tmin,tmax,time_step)        

    
    def get_initial_parameters(self, suffix = ""):
        
        """
        Method that returns the initial parameters 
        
        Parameters
        ----------
        
        suffix : string (optional)
            An optional string to disambiguate a parameter name (if necessary) 
            

        Returns
        -------
        parameters: python object of data type 'dict'
            A dictionary of initial parameters definitions        
    
        """ 
    pass


    def estimate_percolation(self):
        """
        Method that returns a time series of percolated precipitations quantities
        
        
        Returns
        -------
        output: dictionary of root_zone variables
            dictionary of root_zone variables
                    

            
        """
        pass    
    
        
class Vangenuchten(RootzoneFunctionsParent):
    
    """ 
    
    A subclass of RootzoneFunctionsParent 
    following van Genuchten or Brooks and Corey model

    
    """
    
    
    def __init__(self, time = None, p_dict = None, stress_dict = None, time_step = None, 
                 stress_type = None,  tminstr = None, tmaxstr = None, tminnum = None, tmaxnum = None,
                 regime = 'regime_1', root_zone_model = 'Brooks_and_Corey',test=None):
        
        RootzoneFunctionsParent.__init__(self,time = time, p_dict = p_dict, stress_dict = stress_dict,
                                          tminstr = tminstr, tmaxstr = tmaxstr, tminnum = tminnum,
                                          tmaxnum = tmaxnum,regime = regime, root_zone_model = root_zone_model)

        self._simulated = None


    def __repr__(self):
        
        """String representation of the root zone."""

        string = 'root zone model: {root_zone_model}\n'
        # string += f'tmin: {tminstr}.\n'
        # string += f'tmax: {tmaxstr}.\n'
        return string
        

        
    def get_initial_parameters(self, suffix = ""):
        
        """
        See description of the method in the parent class
        
        """           
        
        if suffix != "":
            suffix = '_' + suffix 

            
        parameters = OrderedDict()
        
        
        pname = "D" + suffix
        parameters[pname] = {
            "isvariable": True,
            "logtransform": True,
            "pname": pname,
            "minvalue": 0.1,
            "maxvalue": 0.9,
            "initvalue":0.5, 
            }  
        
        
        pname = "Ks" + suffix
        parameters[pname] = {
            "isvariable": True,
            "logtransform": True,
            "pname": pname,
            "minvalue": 1.,
            "maxvalue": 0.001,
            "initvalue": 0.5, 
            }    


        pname = "lamda" + suffix
        parameters[pname] = {
            "isvariable": True,
            "logtransform": True,
            "pname": pname,
            "minvalue": 0.4,
            "maxvalue": 0.6,
            "initvalue": 0.5, 
            }    
         
        
        pname = "m" + suffix
        parameters[pname] = {
            "isvariable": True,
            "logtransform": True,
            "pname": pname,
            "minvalue": 0.3,
            "maxvalue": 0.8,
            "initvalue": 0.7, 
            }          
        
        
        pname = "f" + suffix
        parameters[pname] = {
            "isvariable": True,
            "logtransform": True,
            "pname": pname,
            "minvalue": 0.6,
            "maxvalue": 1.2,
            "initvalue": 1., 
            }  

        return parameters      

            
    def estimate_percolation(self): 
        
        """
        See description of the method in the parent class
        
        """      
        
        p_dict = self._p_dict
        stress_dict = self._stress_dict
        time = self.time
        time_step = self.time_step
        
        if self.time is None:
            key = list(stress_dict['prec'].keys())[0]
            time = stress_dict['prec'][key].interpolated[:,0]
            self.time = time
            
        if self.time_step is None:
            time_step = 1.0
        self.time_step = time_step
            

        p_indexes = p_dict['p_indexes']
        regime = self.regime
        indexes = p_indexes['root_zone'][regime]['basicparam']
        root_zone_model = self.root_zone_model
        
        _p = []
    
        for ind in indexes:
            p = p_dict['p'][ind]
            if p_dict['logtransform'][ind] == True:
                p = np.exp(p)
            _p.append(p)


        teta_min = 0.2
        teta_max = 0.7

        if root_zone_model == 'vangenuchten':    
            D = _p[0] #depth root zone
            ks = _p[1] #k sat
            lamda = _p[2] # lamda of Van Genuchten see Berendrecht p74
            m = _p[3] # m of Van Genuchten see Berendrecht p74
            f = _p[4] # interception coefficient
    
            def percolation_rate(teta,ks,lamda,m):
                
                term1 = 1-pow(teta,1.0/m)
                term2 = pow(term1,m)
                term3 = 1-term2
                term4 = pow(term3,2)
                output = ks*pow(teta,lamda)*term4
                return output
        
        
        elif root_zone_model == 'Brooks_and_Corey':    
            D = _p[0] #depth root zone
            ks = _p[1] #k sat
            lamda = _p[2] # lamda of Van Genuchten see Berendrecht p74
            m = 0 # m of Van Genuchten set to 0
            f = _p[4] # interception coefficient
            teta_init = 0.3 #initial degree of water saturation, not a parameter because too unsensitive
            
  
    
            def percolation_rate(teta,ks,lamda,m):
                
                # b = 0
                # b = 2.5
                # b = 3.0
                # output = ks*pow(teta,b+2/lamda)#attention, lamda is a different lamda than in case 1 based on Brooks and Correy, zie paper Van Genuchten
                output = ks*pow(teta,lamda)
    
                return output        
        
    
        else:    # for the moment, just make Brooks & Corey the default root zone model, so repeat the option root_zone_model == 'Brooks_and_Corey'
            D = _p[0] #depth root zone
            ks = _p[1] #k sat
            lamda = _p[2] # lamda of Van Genuchten see Berendrecht p74
            m = 0 # m of Van Genuchten set to 0
            f = _p[4] # interception coefficient
            teta_init = 0.3 # initial dimensionless water saturation of root zone (fraction)
    
            def percolation_rate(teta,ks,lamda,m):
    
                return ks*pow(teta,lamda)          

        cumulated_precipitation = np.empty((len(time),2))
        cumulated_precipitation[:,0] = time   
        
        cumulated_runoff = np.empty((len(time),2))
        cumulated_runoff[:,0] = time        

        teta_time_series = np.empty((len(time),2)) # dimensionless water saturation of root zone (fraction)
        teta_time_series[:,0] = time

        percolation_increments = np.empty((len(time),2))
        percolation_increments[:,0] = time      
        
                 
        cumulated_percolation = np.empty((len(time),2))
        cumulated_percolation[:,0] = time
        
        
        actual_evaporation = np.empty((len(time),2))
        actual_evaporation[:,0] = time
        
        
        cumulated_potential_evaporation = np.empty((len(time),2))
        cumulated_potential_evaporation[:,0] = time
        
        cumulated_actual_evaporation = np.empty((len(time),2))
        cumulated_actual_evaporation[:,0] = time        
        
        

        
        ft = np.empty((len(time),2)) # time series of fractions of reference evaporation
        ft[:,0] = time
        
        
        interception_model = 'constant'
        
        if interception_model == 'constant':
            ft[:,1] = f # for the moment f is constant
        
        elif interception_model == '2-seasons':
            seasonsnp.empty(len(time))
            for j in range(0,len(time)):
                time_num = time[j]
                t = num2date(time_num)
                month_num = t.month
                seasons[j] = month_num
    
            for i in range(0,len(time)):
                if seasons[i] in [11,12,1,2,3]:#winter
                    ft[i,1] = f1
    
                else:#summer
                    ft[i,1] = f2
        
        elif interception_model == 'sigmoid': # to do: implement the option of f as sigmoid to simulate land use change
            ft[i,1] = f1+(f2-f1)*1.0/(1+exp(-shape_fact*(float(i)/len(time)-switch)))

        #test
        tstr1 = '01-07-1980'
        dt=datetime.strptime(tstr1,'%d-%m-%Y')
        time_num1=date2num(dt)    

        for i in range(0,len(time)):
            if i == 0:
                teta_previous_time = teta_init
            else:
                teta_previous_time = teta_time_series[i-1,1]
    


            key = list(stress_dict['prec'].keys())[0]
            preprocessed_prec = stress_dict['prec'][key].interpolated 
            key = list(stress_dict['evap'].keys())[0]
            preprocessed_evap = stress_dict['evap'][key].interpolated 
        
            precipitation_increment = preprocessed_prec[i,1]
            potential_evaporation_increment =  ft[i,1]*preprocessed_evap[i,1]
            
            percolrate = percolation_rate(teta_previous_time,ks,lamda,m)
            percolation_increment = min(percolrate*time_step,teta_previous_time*D,precipitation_increment)
            
            # print('453 teta_previous_time*D',teta_previous_time*D)
            # print('453 percolrate*time_step',percolrate*time_step)
            # print('454 precipitation_increment',precipitation_increment)
                        
            

            teta = teta_previous_time + (precipitation_increment
                                         - percolation_increment 
                                         - potential_evaporation_increment) / D
            
            # if i<10:
            #     print('469 teta',teta)
            #     print('469 teta_previous_time',teta_previous_time)
            #     print('470 precipitation_increment',precipitation_increment)
            #     print('471 percolation_increment',percolation_increment)
            #     print('472 potential_evaporation_increment',potential_evaporation_increment)
            #     print('473 D',D)
            #     input()
            
            # if time[i]>=time_num1:
            #     print('477 date',num2date(time[i]))
            #     print('478 percol_per_time_step_root_zone',percolation_increment)
            #     print('456 evap[i,1]',preprocessed_evap[i,1])
            #     print('457 prec[i,:]',preprocessed_prec[i,:])
            #     print('458 teta_previous_time',teta_previous_time)
            #     print('458 teta',teta)
            #     input()
            
            
                        
            
            
            
            reduced_evaporation_increment = 0
            actual_evaporation_increment = potential_evaporation_increment
            runoff_increment = 0
    
            if teta < teta_min:
                teta = teta_min
                reduced_evaporation_increment = potential_evaporation_increment
                actual_evaporation_increment = 0
                percolation_increment = 0
                
            if teta > teta_max:
                runoff_increment = D*(teta - teta_max)                
                teta = teta_max
            
             #test to check if beside rootzone everything works fine  
            # percolation_increment = precipitation_increment- potential_evaporation_increment
            
            teta_time_series[i,1] = teta
            percolation_increments[i,1] = percolation_increment
            actual_evaporation[i,1] = actual_evaporation_increment

            if i == 0:
                cumulated_actual_evaporation[i,1] = actual_evaporation_increment
                cumulated_potential_evaporation[i,1] = potential_evaporation_increment
                cumulated_precipitation[i,1] = precipitation_increment
                cumulated_runoff[i,1] = runoff_increment
                cumulated_percolation[i,1] = percolation_increment
          
 
            else:
                cumulated_actual_evaporation[i,1] = cumulated_actual_evaporation[i-1,1] + actual_evaporation_increment
                cumulated_potential_evaporation[i,1] = cumulated_potential_evaporation[i-1,1] + potential_evaporation_increment
                cumulated_precipitation[i,1] = cumulated_precipitation[i-1,1] + precipitation_increment
                cumulated_runoff[i,1] = cumulated_runoff[i-1,1] = runoff_increment
                cumulated_percolation[i,1] = cumulated_percolation[i-1,1] + percolation_increment
   

    

        self.water_saturation = teta_time_series
        self.percolation_per_time_step = percolation_increments
        self.percolation = self.percolation_per_time_step / time_step
        self.actual_evaporation = actual_evaporation
        self.cumulated_actual_evaporation = cumulated_actual_evaporation
        self.cumulated_potential_evaporation = cumulated_potential_evaporation
        self.cumulated_precipitation = cumulated_precipitation
        self.cumulated_runoff = cumulated_runoff
        self.cumulated_percolation = cumulated_percolation
        
        
        # fig=figure()
        # legend_list = []
        # ax = fig.add_subplot(111)
        # ax.set_title('519 inspect water_saturation ')
        # ax.plot_date(teta_time_series[:,0],teta_time_series[:,1],'b')
        # legend_list.append('water_saturation')
        # leg = ax.legend((legend_list),loc='upper left',shadow=False)
        # ax.set_xlabel('time')
        # ax.grid(False)
        # show()


        # fig=figure()
        # legend_list = []
        # ax = fig.add_subplot(111)
        # ax.set_title('519 inspect percolation_per_time_step ')

        # ax.plot_date(preprocessed_prec[:,0],preprocessed_prec[:,1],'g')
        # legend_list.append('preprocessed_prec')            
        # ax.plot_date(percolation_increments[:,0],percolation_increments[:,1],'b')
        # legend_list.append('percolation_per_time_step')
        
        # leg = ax.legend((legend_list),loc='upper left',shadow=False)
        # ax.set_xlabel('time')
        # ax.grid(False)
        # show()
        
        
        
        # fig=figure()
        # legend_list = []
        # ax = fig.add_subplot(111)
        # ax.set_title('519 cumulated_potential_evaporation')
        # ax.plot_date(cumulated_potential_evaporation[:,0],cumulated_potential_evaporation[:,1],'b')
        # legend_list.append('cumulated_potential_evaporation')
        # ax.plot_date(cumulated_actual_evaporation[:,0],cumulated_actual_evaporation[:,1],'g')
        # legend_list.append('cumulated_actual_evaporation')            
        
        # leg = ax.legend((legend_list),loc='upper left',shadow=False)
        # ax.set_xlabel('time')
        # ax.grid(False)
        # show()
        
        
        

        return percolation_increments
    

                                  
if __name__ == "__main__":
    from stresstimeseries import Stress
    from headstimeseries import Heads
    from monitoringwell import MonitoringWell
    from preprocessedseries import Preprocessed
    from stressconvolutionV2 import StressConvolution
    from parameterslogistic import ParametersLogistic
    from modeldefinition import ModelDefinition
    from simulationV2 import Simulation

    abs_path = os.path.dirname(os.path.abspath(__file__))
    abs_path_splitted = abs_path.split('\\')
    path_to_parent_folder_elements = abs_path_splitted[:-1]
    path_to_parent_folder = os.path.join(*path_to_parent_folder_elements)
    splitted = path_to_parent_folder.split(':') #trick  to  repair  path
    path_to_parent_folder = splitted[0]+':\\'+splitted[1] #trick to repair path
    path_to_file_folder = os.path.join(path_to_parent_folder,'resources')  

    stresses_dict = {}
    
    # path_to_file = os.path.join(path_to_file_folder,'260_De_Bilt_PREC_19010101_20200910.csv')
    path_to_file = os.path.join(path_to_file_folder,'672_PREC_Hellendoorn_19510101_20160731_corrected_for_interception2mm.csv')
    # path_to_file = os.path.join(path_to_file_folder,'prec_giersbergen.csv')

    key = basename(path_to_file)
    stresses_dict['prec'] = {}
    stresses_dict['prec'][key] = Stress.read_from_csv(csv_file = path_to_file, stress_type = 'prec',cumulative = True)    
    

    path_to_file=os.path.join(path_to_file_folder,'260_De_Bilt_EVAP_19010101_20200910.csv')

    key = basename(path_to_file)
    stresses_dict['evap'] = {}
    stresses_dict['evap'][key] = Stress.read_from_csv(csv_file = path_to_file, stress_type = 'evap',cumulative = True)    

    tminstr = '01-01-1900 08:00:00'
    tmaxstr = '31-12-2100 08:30:00'
    
    path_to_file_folder = os.path.join(path_to_parent_folder,'resources')  
    path_to_file = os.path.join(path_to_file_folder,'28AP0093_1.txt')
    path_to_file = os.path.join(path_to_file_folder,'28CP0197_2.txt')
    
    heads = Heads.read_from_csv(path_to_file, tminstr = tminstr, tmaxstr = tmaxstr)

    arguments_dic = {}

    for key in stresses_dict:
        for e in stresses_dict[key]:
            stresses_dict[key][e].plot()

    md = ModelDefinition().model_definition
    memory_dict = ModelDefinition().memory_dict
    parametersLogistic = ParametersLogistic(md,stress_dict = stresses_dict, heads = heads)
    
    time_step = 1. #1./(24.*2)
    preprocessed = Preprocessed(heads = heads,stresses_dict = stresses_dict, time_step = time_step, 
                                model_definition = md, memory_dict = memory_dict,tminstr = tminstr, tmaxstr = tmaxstr)
    time, heads, stresses_dict, Nint_dic = preprocessed.preprocess()
    

    parametersLogistic = ParametersLogistic(md)
    p_dict = parametersLogistic.assemble_p()


    # arguments_dic = {}
    
    # p_dict = parametersLogistic.assemble_p()
    
    # if 'db' in p_dict['p_indexes']:
    #     if md['db']['funct_type'] == 'constant':
    #         if md['db']['initial values from observed median'] == True:
    #             index = p_dict['p_indexes']['db']['regime_1']['basicparam'][0]
    #             db = np.median(heads.observed[:,1])# - 1
    #             db_min = db - 2.
    #             db_max = db + 1.
    #             if p_dict['logtransform'][index] == True:
    #                 if db_min <= 0:
    #                     p_dict['logtransform'][index] = False
                
    #             if p_dict['logtransform'][index] == True:
    #                 db = np.log(db)
    #                 db_min = np.log(db_min)
    #                 db_max = np.log(db_max)
    
    #             p_dict['p'][index] = db
    #             p_dict['p_min'][index] = db_min
    #             p_dict['p_max'][index] = db_max

    # simulation = Simulation(time = time, p_dict = p_dict, stress_dict = stresses_dict, Nint_dic = Nint_dic,
    #                         time_step = time_step, heads = heads, simulate_residuals = False)
    # sim = simulation.simulate()
    # plotax = simulation.plot()
    if md['root_zone']['apply_root_zone'] == True:
        vg = Vangenuchten(p_dict = p_dict, stress_dict = stresses_dict, time_step = 1.0, 
                                  root_zone_model = 'Brooks_and_Corey')
        output = vg.estimate_percolation()

        print('583 percolation',output)

    