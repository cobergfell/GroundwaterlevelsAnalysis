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
path_to_src_folder = os.path.join(path_to_parent_folder,'src') 
curdir=os.getcwd()
path_to_resources_folder = os.path.join(curdir,'resources') 




import sys
import os.path
# sys.path.append(
#     os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir)))

sys.path.append(path_to_src_folder)
print('44 path_to_src_folder',path_to_src_folder)
print('45 sys.path', sys.path)
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
import export




stresses_dict = {}


######################## Test using a piezometer on de Sallands sand ridge in the Netherlands

path_to_file = os.path.join(path_to_resources_folder,'672_PREC_Hellendoorn_19510101_20160731_corrected_for_interception2mm.csv')
path_to_file = os.path.join(path_to_resources_folder,'672_PREC_Hellendoorn_19510101_20160731.csv')
path_to_file = os.path.join(path_to_resources_folder,'562_TIEL_PREC_19110602_20231010.csv')

key = basename(path_to_file)
stresses_dict['prec'] = {}
stresses_dict['prec'][key] = Stress.read_from_csv(path_to_file, stress_type = 'prec',cumulative = True)    


path_to_file=os.path.join(path_to_resources_folder,'260_De_Bilt_EVAP_19010101_20200910.csv')
path_to_file=os.path.join(path_to_resources_folder,'260_De_Bilt_EVAP_19010101_20131015.csv')
path_to_file=os.path.join(path_to_resources_folder,'260_De_Bilt_EVAP_19010103_20231024.csv')

   
key = basename(path_to_file)
stresses_dict['evap'] = {}
stresses_dict['evap'][key] = Stress.read_from_csv(path_to_file, stress_type = 'evap',cumulative = True)    



path_to_file = os.path.join(path_to_resources_folder,'Dodewaard-select.txt')
path_to_file = os.path.join(path_to_resources_folder,'Dodewaard.txt')

key = basename(path_to_file)
stresses_dict['riv'] = {}
stresses_dict['riv'][key] = Stress.read_from_csv(path_to_file, stress_type = 'riv',cumulative = True)    



tminstr = '01-01-1900 00:00:00'
tmaxstr = '31-12-2023 00:00:00'


tminstr = '01-01-1997 00:00:00'
tmaxstr = '31-12-2000 00:00:00'


# tminstr = '31-12-1981 00:00:00'
# tmaxstr = '31-12-1995 00:00:00'    


path_to_file = os.path.join(path_to_resources_folder,'B39G0027001_1.csv')
# path_to_file = os.path.join(path_to_resources_folder,'B39H0048001_1.csv')

# heads = Heads.read_from_csv(path_to_file, tminstr = tminstr, tmaxstr = tmaxstr)
heads = Heads.read_from_DINO(path_to_file, tminstr = tminstr, tmaxstr = tmaxstr)


path_to_file = os.path.join(path_to_resources_folder,'HBpb001-1.csv')
# path_to_file = os.path.join(path_to_resources_folder,'HBpb003-1.csv')
# path_to_file = os.path.join(path_to_resources_folder,'HBpb004.csv')
tminstr = '01-01-1900 00:00:00'
tmaxstr = '31-12-2100 00:00:00'


heads = Heads.read_from_DINO(path_to_file, tminstr = tminstr, tmaxstr = tmaxstr, artdiverformat = True )

heads.plot()

glg,gg,ghg = heads.generate_gxg()

print('125 glg,gg,ghg',glg,gg,ghg)


arguments_dict = {}

for key in stresses_dict:
    for e in stresses_dict[key]:
        stresses_dict[key][e].plot()
        

#md = ModelDefinition().model_definition # case 1: model definition specified in modeldefinition.py
#md = ModelDefinition(read_from_resources = True) # case 2: model definition specified in json file in resources map

#curdir=os.getcwd()
abs_path = os.path.dirname(os.path.abspath(__file__))
abs_path_splitted = abs_path.split('\\')
working_directory_elements = abs_path_splitted[:-1]
working_directory = os.path.join(*working_directory_elements)
splitted = working_directory.split(':') # trick to repair path
working_directory = splitted[0]+':\\'+splitted[1] # trick to repair path

# curdir = os.getcwd()

mydir = os.path.join(path_to_parent_folder,'example_5_Druten')         
if not os.path.isdir(mydir):
    os.mkdir(mydir)  
    
filepath = os.path.join(mydir,'model_definition.json')

md = ModelDefinition().from_jason(filepath=filepath) # case 3: model definition specified in json json file in filepath

memory_dict = ModelDefinition().memory_dict
# parametersLogistic = ParametersLogistic(md)

# time_step = 1./(24.*2)
time_step = 1.
time_step_targets = 14.


preprocessed = Preprocessed(heads = heads,stresses_dict = stresses_dict, time_step = time_step,
                            time_step_targets = time_step_targets, memory_dict = memory_dict,tminstr = tminstr, 
                            tmaxstr = tmaxstr, model_definition = md).preprocess() 
    
time = preprocessed.time
time_step = preprocessed.time_step
heads = preprocessed.heads
time_step_targets = preprocessed.time_step_targets

stresses_dict = preprocessed.stresses_dict 
Nint_dict = preprocessed.Nint_dict
all_data_for_neural_networks = preprocessed.all_data_for_neural_networks
all_targets_for_neural_networks = preprocessed.all_targets_for_neural_networks        
    
if isinstance(heads,list):
    _heads = heads[0]
else:
    _heads = heads    
 
parametersLogistic = ParametersLogistic(md,stresses_dict = stresses_dict, heads = _heads)
p_dict = parametersLogistic.assemble_p()


settings = {}
settings['all_piezometers_share_same_model'] = True



potimizer = Poptimizer(heads = heads, time = time, p_dict = p_dict, time_step = time_step,time_step_targets = time_step_targets,
                        stresses_dict = stresses_dict, model_residuals = False, model_definition = md,
                        Nint_dict = Nint_dict, settings = settings,
                        analytical_jacobian = False, maxiter = 100)    


popt, pcov, pcor, pstdev, p_dict, expvar, expvarnoise  = potimizer.poptimizer_home_brew()    


riv = stresses_dict[key][e]


to_plot = [heads.interpolated, heads.modeled, riv.interpolated]
legend_list = ['interpolated observations', 'modeled','riv']

plotax = heads.plot_multiple_series(tseries_list = to_plot, legend_list = legend_list, share_axes = True )
output_string = heads.generate_output(model_definition = md)
print('168 output_string: ',output_string)   
ax = potimizer.residuals_check_diagrams(heads = _heads)    


tminstr = '01-01-1997 00:00:00'
tmaxstr = '28-10-2023 00:00:00'


preprocessed = Preprocessed(heads = heads,stresses_dict = stresses_dict, time_step = time_step,
                            time_step_targets = time_step_targets, memory_dict = memory_dict,tminstr = tminstr, 
                            tmaxstr = tmaxstr, model_definition = md, forcast = True).preprocess() 
    
time = preprocessed.time
time_step = preprocessed.time_step
heads = preprocessed.heads
time_step_targets = preprocessed.time_step_targets

simulation = Simulation(time = time, p_dict = p_dict, stresses_dict = stresses_dict, 
heads = heads, time_step = time_step, Nint_dict = Nint_dict, 
tminstr = tminstr, tmaxstr = tmaxstr, settings = settings, 
model_residuals = False, model_definition = md)
          

sim = simulation.simulate()


to_plot = [heads.interpolated, heads.modeled, riv.interpolated,sim]
legend_list = ['interpolated observations', 'modeled','riv','sim']

plotax = heads.plot_multiple_series(tseries_list = to_plot, legend_list = legend_list, share_axes = True )

    
export_path = export.export_heads_nparray_to_ascii(sim, metadata = heads.metadata)



# heads = Heads.read_from_csv(path_to_file, tminstr = tminstr, tmaxstr = tmaxstr)
heads_sim = Heads.read_from_csv(export_path)
heads_sim.plot()

glg,gg,ghg = heads_sim.generate_gxg()

print('251 glg,gg,ghg',glg,gg,ghg)




# sys.path.remove('D:\\dev\\_my_own_program\\src')
sys.path.remove(path_to_src_folder)

# check if sys,path is cleaned
print('182 sys.path', sys.path)
input()
