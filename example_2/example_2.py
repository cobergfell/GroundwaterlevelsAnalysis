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
print('40 sys.path', sys.path)
input()

                                  
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

path_to_file = os.path.join(path_to_resources_folder,'375_PREC_19510101_20130603.csv')
# path_to_file = os.path.join(path_to_resources_folder,'radar_RAAI2.csv')
key = basename(path_to_file)
stresses_dict['prec'] = {}
stresses_dict['prec'][key] = Stress.read_from_csv(path_to_file, stress_type = 'prec',cumulative = True)    

path_to_file = os.path.join(path_to_resources_folder,'375_EVAP_19510101_20130603.csv')
# path_to_file = os.path.join(path_to_file_folder,'evap_volkel_raaien_1_2_3.csv')
key = basename(path_to_file)
stresses_dict['evap'] = {}
stresses_dict['evap'][key] = Stress.read_from_csv(path_to_file, stress_type = 'evap',cumulative = True)    

path_to_file = os.path.join(path_to_resources_folder,'raai2_251A_20082012_repared.csv')
path_to_file = os.path.join(path_to_resources_folder,'raai2_251A_20082012_local_coordinates.csv')
key = basename(path_to_file)
stresses_dict['riv'] = {}
stresses_dict['riv'][key] = Stress.read_from_csv(path_to_file, stress_type = 'riv',cumulative = False)    


tminstr = '01-01-1900 08:00:00'
tmaxstr = '31-12-2100 08:30:00'


path_to_piezometers_directory = os.path.join(path_to_resources_folder,'piezometers_paper2') 
path_to_piezometers_directory = os.path.join(path_to_resources_folder,'piezometers_paper2_local_coordinates')

heads_list = []
list_of_piezometers_files_names = os.listdir(path_to_piezometers_directory)
for filename in list_of_piezometers_files_names: 
    path_to_file = os.path.join(path_to_piezometers_directory,filename)
    heads = Heads.read_from_csv(path_to_file, tminstr = tminstr, tmaxstr = tmaxstr) 

    if heads.name in ['_P12']:
        heads.L = 0
    else:
        heads.L = 1
    
    
    heads_list.append(heads)
    #heads.plot()


for key in stresses_dict:
    for e in stresses_dict[key]:
        stresses_dict[key][e].plot()

# md = ModelDefinition().model_definition
# memory_dict = ModelDefinition().memory_dict
# parametersLogistic = ParametersLogistic(md)


memory_dict = {'prec':1000,
      'evap':1000,
      'riv':1000,
      'pump':50,
      'noise':50}         


time_step = 1. #1.0/(2*24)
time_step_targets = 1.
    


###################   fit all piezometers together

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

mydir = os.path.join(path_to_parent_folder,'example_2')         
if not os.path.isdir(mydir):
    os.mkdir(mydir)  
    
filepath = os.path.join(mydir,'model_definition.json')

md = ModelDefinition().from_jason(filepath=filepath) # case 3: model definition specified in json json file in filepath




preprocessed = Preprocessed(heads = heads_list,stresses_dict = stresses_dict, time_step = time_step,
                            time_step_targets = time_step_targets, memory_dict = memory_dict,tminstr = tminstr, 
                            tmaxstr = tmaxstr, model_definition = md).preprocess() 
    
time = preprocessed.time
heads = preprocessed.heads
time_step = preprocessed.time_step
time_step_targets = preprocessed.time_step_targets


for _heads in heads:

    stresses_dict = preprocessed.stresses_dict 
    Nint_dict = preprocessed.Nint_dict
    all_data_for_neural_networks = preprocessed.all_data_for_neural_networks
    all_targets_for_neural_networks = preprocessed.all_targets_for_neural_networks        
        

settings = {}
settings['all_piezometers_share_same_model'] = True


parametersLogistic = ParametersLogistic(md,stresses_dict = stresses_dict, heads = _heads)
p_dict = parametersLogistic.assemble_p()



potimizer = Poptimizer(heads = heads, time = time, p_dict = p_dict, time_step = time_step, time_step_targets = time_step_targets,
                        stresses_dict = stresses_dict, model_residuals = True, model_definition = md,
                        Nint_dict = Nint_dict, settings = settings,
                        analytical_jacobian = True, maxiter = 5)    


popt, pcov, pcor, pstdev, p_dict, expvar, expvarnoise  = potimizer.poptimizer_home_brew(delta = 0.1)  

# to reproduce approximately the results of paper 2, set  maxiter=5 in Poptimizer and delta=0.1  in poptimizer_home_brew()
# but: this leads to a sub optimum, a fact I was not aware of at the time of the paper, however, it does not change the message of the paper and
# the fact that L is underestimated does not change the conclusion of the paper (in factt, we had already commented on th efact that the value of L
# was not important for the rest of the paper)
# To find what appears as the global optimum set  maxiter1=19 (or more) in Poptimizer and set poptimizer_home_brew(delta=2.0)


for _heads in heads:

    to_plot = [_heads.interpolated, _heads.modeled]
    legend_list = ['interpolated observations', 'modeled']
    
    plotax = _heads.plot_multiple_series(tseries_list = to_plot, legend_list = legend_list, share_axes = True )
    output_string = _heads.generate_output(model_definition = md)
    print('194 Check output_string: ',output_string)   
    ax = potimizer.residuals_check_diagrams(heads = _heads)

# sys.path.remove('D:\\dev\\_my_own_program\\src')
sys.path.remove(path_to_src_folder)

#check if sys,path is cleaned
# print('40 sys.path', sys.path)
# input()
