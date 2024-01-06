# -*- coding: utf-8 -*-
"""
Created on Sun Jul 16 13:38:49 2023

@author: christophe
"""
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  8 21:10:14 2023

@author: christophe
"""


import os
import numpy as np
from scipy.optimize import least_squares, leastsq
from datetime import datetime
from matplotlib.dates import date2num, num2date
from logging import getLogger



#curdir=os.getcwd()
abs_path = os.path.dirname(os.path.abspath(__file__))
abs_path_splitted = abs_path.split('\\')
path_to_parent_folder_elements = abs_path_splitted[:-1]
path_to_parent_folder = os.path.join(*path_to_parent_folder_elements)
splitted = path_to_parent_folder.split(':')#trick  to  repair  path
path_to_parent_folder = splitted[0]+':\\'+splitted[1]#trick  to  repair  path
path_to_resources_folder = os.path.join(path_to_parent_folder,'resources') 
path_to_src_folder = os.path.join(path_to_parent_folder,'src') 
# path_to_src_folder = 'D:\\dev\\_my_own_program\\src'


import sys
import os.path
# sys.path.append(
#     os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir)))

sys.path.append(path_to_src_folder)
# print('line 40 check sys.path', sys.path)
# input()

                                  
from stresstimeseries import Stress                                  
from headstimeseries import Heads
from parameterslogistic import ParametersLogistic
from modeldefinition import ModelDefinition
from modelnoise import ResidualsDecayExponentially
from utilities import p_dict_copy
from harmonics import Harmonics
from harmonics_toolkit import *
from matplotlib.pyplot import figure, show,savefig,clabel,contour,setp,gcf,getp,gca,close
from preprocessedseries import Preprocessed
from jacobian import Jacobian
from simulation import Simulation                      
from poptimizer_test import Poptimizer                               



stresses_dict = {}

# path_to_file = os.path.join(path_to_resources_folder,'260_De_Bilt_PREC_19010101_20200910.csv')
path_to_file = os.path.join(path_to_resources_folder,'reformated_prec_giersbergen.csv')

key = basename(path_to_file)
stresses_dict['prec'] = {}
stresses_dict['prec'][key] = Stress.read_from_csv(path_to_file, stress_type = 'prec',cumulative = True)    


path_to_file=os.path.join(path_to_resources_folder,'260_De_Bilt_EVAP_19010101_20200910.csv')
path_to_file=os.path.join(path_to_resources_folder,'reformated_evap_eindhoven.csv')

key = basename(path_to_file)
stresses_dict['evap'] = {}
stresses_dict['evap'][key] = Stress.read_from_csv(path_to_file, stress_type = 'evap',cumulative = True)    


    
stresses_dict['pump'] = {}
path_to_files_folder = os.path.join(path_to_resources_folder,'selected_wells_paper1_reformatted') 
# path_to_files_folder = os.path.join(path_to_resources_folder,'selected_wells_test')
list_of_files_names = os.listdir(path_to_files_folder)
for filename in list_of_files_names:
    path_to_file = os.path.join(path_to_files_folder,filename)
    key = basename(path_to_file)
    stresses_dict['pump'][key] = Stress.read_from_csv(path_to_file, stress_type = 'pump',cumulative = False)    


tminstr = '01-01-1900 08:00:00'
tmaxstr = '31-12-2100 08:30:00'


path_to_piezometers_directory = os.path.join(path_to_resources_folder,'piezometers_paper1') 

heads_list = []
list_of_piezometers_files_names = os.listdir(path_to_piezometers_directory)
for filename in list_of_piezometers_files_names: 
    path_to_file = os.path.join(path_to_piezometers_directory,filename)
    heads = Heads.read_from_csv(path_to_file, tminstr = tminstr, tmaxstr = tmaxstr) 
    
    heads_list.append(heads)
    #heads.plot()


for key in stresses_dict:
    for e in stresses_dict[key]:
        stresses_dict[key][e].plot()

# md = ModelDefinition().model_definition
# memory_dict = ModelDefinition().memory_dict
# parametersLogistic = ParametersLogistic(md)


memory_dict = {'prec':600,
      'evap':600,
      'riv':10,
      'pump':10,
      'noise':50}         


time_step = 1./(24.*2)
time_step_targets = 1.5

if isinstance(heads,list):
    _heads = heads[0]
else:
    _heads = heads
    
    
# read model definition
#curdir=os.getcwd()
abs_path = os.path.dirname(os.path.abspath(__file__))
abs_path_splitted = abs_path.split('\\')
working_directory_elements = abs_path_splitted[:-1]
working_directory = os.path.join(*working_directory_elements)
splitted = working_directory.split(':') # trick to repair path
working_directory = splitted[0]+':\\'+splitted[1] # trick to repair path

# curdir = os.getcwd()

mydir = os.path.join(path_to_parent_folder,'example_1')         
if not os.path.isdir(mydir):
    os.mkdir(mydir)  
    
filepath = os.path.join(mydir,'model_definition.json')

md = ModelDefinition().from_jason(filepath=filepath) # case 3: model definition specified in json json file in filepath


parametersLogistic = ParametersLogistic(md,stresses_dict = stresses_dict, heads = _heads)
p_dict = parametersLogistic.assemble_p()


preprocessed = Preprocessed(heads = heads_list,stresses_dict = stresses_dict, time_step = time_step,
                            time_step_targets = time_step_targets, memory_dict = memory_dict,tminstr = tminstr, 
                            tmaxstr = tmaxstr, model_definition = md).preprocess() 
    
time = preprocessed.time
heads = preprocessed.heads
time = preprocessed.time
time_step = preprocessed.time_step
time_step_targets = preprocessed.time_step_targets


stresses_dict = preprocessed.stresses_dict 
Nint_dict = preprocessed.Nint_dict

all_data_for_neural_networks = preprocessed.all_data_for_neural_networks
all_targets_for_neural_networks = preprocessed.all_targets_for_neural_networks        
    

settings = {}
settings['all_piezometers_share_same_model'] = True



potimizer = Poptimizer(heads = heads, time = time, p_dict = p_dict, time_step = time_step, time_step_targets = time_step_targets,
                        stresses_dict = stresses_dict, model_residuals = True, model_definition = md, maxiter = 100,
                        Nint_dict = Nint_dict, settings = settings, analytical_jacobian = True)    


popt, pcov, pcor, pstdev, p_dict, expvar, expvarnoise  = potimizer.poptimizer_home_brew()    




for _heads in heads:

    to_plot = [_heads.interpolated, _heads.modeled]
    legend_list = ['interpolated observations', 'modeled']
    
    plotax = _heads.plot_multiple_series(tseries_list = to_plot, legend_list = legend_list, share_axes = True )
    output_string = _heads.generate_output(model_definition = md)
    print('192 Check output_string: ',output_string)   


    ax = potimizer.residuals_check_diagrams(heads = _heads)

# sys.path.remove('D:\\dev\\_my_own_program\\src')
sys.path.remove(path_to_src_folder)

#check if sys,path is cleaned
# print('40 sys.path', sys.path)
# input()





