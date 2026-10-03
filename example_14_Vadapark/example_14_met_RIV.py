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
# print('44 path_to_src_folder',path_to_src_folder)
# print('45 sys.path', sys.path)
# input()

                                  
from stresstimeseries import Stress                                  
from headstimeseries import Heads
from parameterslogistic import ParametersLogistic
from modeldefinition import ModelDefinition
from modelnoise import ResidualsDecayExponentially
from utilities import p_dict_copy
from harmonics import Harmonics
from harmonics_toolskit import *
from matplotlib.pyplot import figure, show,savefig,clabel,contour,setp,gcf,getp,gca,close
from preprocessedseries import Preprocessed
from jacobian import Jacobian
from simulation import Simulation                      
from poptimizer_test import Poptimizer                               
import export
import unitresp


stresses_dict = {}

path_to_file = os.path.join(path_to_resources_folder,'260_DeBilt_PREC_19800101_20241231.csv')
path_to_file = os.path.join(path_to_resources_folder,'280_Eelde_PREC_19700101_21001231.csv')
path_to_file = os.path.join(path_to_resources_folder,'344_Rotterdam_PREC_19700101_21001231.csv')
path_to_file = os.path.join(path_to_resources_folder,'286_NieuwBeerta_PREC_19700101_21001231.csv')
path_to_file = os.path.join(path_to_resources_folder,'275_Deelen_PREC_19700101_21001231.csv')
# path_to_file = os.path.join(path_to_resources_folder,'260_DeBilt_PREC_19700101_21001231.csv')


key = basename(path_to_file)
stresses_dict['prec'] = {}
stresses_dict['prec'][key] = Stress.read_from_csv(path_to_file, stress_type = 'prec',cumulative = True)    

path_to_file=os.path.join(path_to_resources_folder,'260_DeBilt_EVAP_19800101_20241231.csv')
path_to_file=os.path.join(path_to_resources_folder,'280_Eelde_EVAP_19700101_21001231.csv')
path_to_file=os.path.join(path_to_resources_folder,'344_Rotterdam_EVAP_19700101_21001231.csv')
path_to_file=os.path.join(path_to_resources_folder,'286_NieuwBeerta_EVAP_19700101_21001231.csv')
path_to_file = os.path.join(path_to_resources_folder,'275_Deelen_EVAP_19700101_21001231.csv')
# path_to_file = os.path.join(path_to_resources_folder,'260_DeBilt_EVAP_19700101_21001231.csv')


key = basename(path_to_file)
stresses_dict['evap'] = {}
stresses_dict['evap'][key] = Stress.read_from_csv(path_to_file, stress_type = 'evap',cumulative = True)    


# path_to_file = os.path.join(path_to_resources_folder,'Dodewaard-select.txt')
# path_to_file = os.path.join(path_to_resources_folder,'Dodewaard.txt')
path_to_file = os.path.join(path_to_resources_folder,'20222180_Groenzoom_v8_Bovenstrooms_Berkel.txt')
path_to_file = os.path.join(path_to_resources_folder,'20222180_Groenzoom_v8_Bovenstrooms_Pijnacker.txt')
# path_to_file = os.path.join(path_to_resources_folder,'20222180_Groenzoom_v8_Berkselsevaart_Pijnacker.txt')
path_to_file = os.path.join(path_to_resources_folder,'Rhenen_Grebbeberg_01-01-2012-01-10-2026.txt')

key = basename(path_to_file)
stresses_dict['riv'] = {}
stresses_dict['riv'][key] = Stress.read_from_csv(path_to_file, stress_type = 'riv',cumulative = False,multiply_values_by = 0.01)


path_to_file = os.path.join(path_to_resources_folder,'B13D0166001.csv')
path_to_file = os.path.join(path_to_resources_folder,'B13D1012001.csv')
path_to_file = os.path.join(path_to_resources_folder,'B13D0206001.csv')
path_to_file = os.path.join(path_to_resources_folder,'B13D1075001.csv')
path_to_file = os.path.join(path_to_resources_folder,'B13D1076001.csv')
path_to_file = os.path.join(path_to_resources_folder,'B39F2906_001_full.csv')
path_to_file = os.path.join(path_to_resources_folder,'B39F2906001.csv')
# path_to_file = os.path.join(path_to_resources_folder,'B39F2905001.csv')






# path_to_file = os.path.join(path_to_resources_folder,'B13D1014001.csv')
# path_to_file = os.path.join(path_to_resources_folder,'B13D1018001.csv')
# path_to_file = os.path.join(path_to_resources_folder,'B13D1019001.csv')
# path_to_file = os.path.join(path_to_resources_folder,'B13D1016001.csv')
# path_to_file = os.path.join(path_to_resources_folder,'B13D1021001.csv')
# path_to_file = os.path.join(path_to_resources_folder,'B13D0033001.csv')
# path_to_file = os.path.join(path_to_resources_folder,'S-BR1.01_20260917165935.csv')
# path_to_file = os.path.join(path_to_resources_folder,'S-BR1.02_20260917165935.csv')
# path_to_file = os.path.join(path_to_resources_folder,'S-BR2.03_20260917165935.csv')
# path_to_file = os.path.join(path_to_resources_folder,'S-BR6.04_20260917165937.csv')


# path_to_file = os.path.join(path_to_resources_folder,'B13D1020001.csv')

# path_to_file = os.path.join(path_to_resources_folder,'20222180_Groenzoom_v8_BE02.txt')

# path_to_file = os.path.join(path_to_resources_folder,'B39A0357001.csv')
# path_to_file = os.path.join(path_to_resources_folder,'B39A2335001.csv')

# heads = Heads.read_from_csv(path_to_file, tminstr = tminstr, tmaxstr = tmaxstr)

# tminstr = '01-01-1900 00:00:00'
# tmaxstr = '31-12-2100 00:00:00'
# heads = Heads.read_from_DINO(path_to_file, tminstr = tminstr, tmaxstr = tmaxstr, artdiverformat = True )



tminstr = '01-01-2008 00:00:00'
tmaxstr = '31-12-2011 00:00:00'

tminstr = '01-01-1980 00:00:00'
tmaxstr = '31-12-2026 00:00:00'
# heads = Heads.read_from_DINO(path_to_file, tminstr = tminstr, tmaxstr = tmaxstr)
heads = Heads.read_from_DINO(path_to_file, tminstr = tminstr, tmaxstr = tmaxstr, artdiverformat = False )

# heads = Heads.read_from_csv(path_to_file, tminstr = tminstr, tmaxstr = tmaxstr)

export_path = export.export_heads_nparray_to_ascii(heads.observed, metadata = heads.metadata)


glg,gg,ghg = heads.generate_gxg()

statistics_observed_time_series = {}
statistics_observed_time_series['GLG'] = glg
statistics_observed_time_series['GG'] = gg
statistics_observed_time_series['GHG'] = ghg
heads.statistics_observed_time_series = statistics_observed_time_series


print('144 glg,gg,ghg',glg,gg,ghg)

heads.plot()

arguments_dict = {}

for key in stresses_dict:
    for e in stresses_dict[key]:
        stresses_dict[key][e].plot()
        

#curdir=os.getcwd()
abs_path = os.path.dirname(os.path.abspath(__file__))
abs_path_splitted = abs_path.split('\\')
working_directory_elements = abs_path_splitted[:-1]
working_directory = os.path.join(*working_directory_elements)
splitted = working_directory.split(':') # trick to repair path
working_directory = splitted[0]+':\\'+splitted[1] # trick to repair path

# curdir = os.getcwd()

mydir = os.path.join(path_to_parent_folder,'example_12_Sellingen')         
mydir = os.path.join(path_to_parent_folder,'example_14_Vadapark')


if not os.path.isdir(mydir):
    os.mkdir(mydir)  
    
filepath = os.path.join(mydir,'model_options.json')
md = ModelDefinition().from_jason(filepath=filepath) 

memory_dict = ModelDefinition().memory_dict

# parametersLogistic = ParametersLogistic(md)

# time_step = 1./(24.*2)
time_step = 1.
time_step_targets = 1.


tminstr = '01-01-1980 00:00:00'
tmaxstr = '31-12-2026 00:00:00'

tminstr = '01-01-1996 00:00:00'
tmaxstr = '31-12-2026 00:00:00'

tminstr_targets = '01-01-2001 00:00:00'
tmaxstr_targets = '31-12-2003 00:00:00'

tminstr_targets = '14-08-2006 00:00:00'
tmaxstr_targets = '30-08-2010 00:00:00'


tminstr_targets = '21-12-2022 00:00:00'
# tminstr_targets = '04-03-2024 00:00:00'
# tminstr_targets = '21-02-2026 00:00:00'
tmaxstr_targets = '01-09-2026 00:00:00'


tminstr_targets = '01-01-1980 00:00:00'
tmaxstr_targets = '31-12-2026 00:00:00'

tminstr_targets = '01-01-2013 00:00:00'
tmaxstr_targets = '31-12-2020 00:00:00'







preprocessed = Preprocessed(heads = heads,stresses_dict = stresses_dict, time_step = time_step,
                            time_step_targets = time_step_targets, memory_dict = memory_dict,tminstr = tminstr, 
                            tmaxstr = tmaxstr, tminstr_targets = tminstr_targets, tmaxstr_targets = tmaxstr_targets, model_definition = md).preprocess() 
    
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

# tminstr = '31-12-2001 00:00:00'
# tmaxstr = '31-12-2003 00:00:00'    


potimizer = Poptimizer(heads = heads, time = time, p_dict = p_dict, time_step = time_step,time_step_targets = time_step_targets,
                        stresses_dict = stresses_dict, model_residuals = False, model_definition = md,
                        Nint_dict = Nint_dict, settings = settings,
                        analytical_jacobian = False, maxiter = 100)    


popt, pcov, pcor, pstdev, p_dict, expvar, expvarnoise  = potimizer.poptimizer_home_brew()    


# riv = stresses_dict[key][e]


# to_plot = [heads.interpolated, heads.modeled, riv.interpolated]
# legend_list = ['interpolated observations', 'modeled','riv']


# to_plot = [heads.observed,heads.modeled,heads.targets]
# colors = ['b-','r-','go']
# legend_list = ['observations', 'modeled','targets']
# plotax = heads.plot_multiple_series(tseries_list = to_plot, colors = colors, legend_list = legend_list, share_axes = True )

output_string = heads.generate_output(model_definition = md)
print('268 output_string: ',output_string)   


print('274 residual inspection line 275 is skipped')
input()
# ax = potimizer.residuals_check_diagrams(heads = _heads)   

export_path = export.export_heads_nparray_to_ascii(heads.observed, metadata = heads.metadata,suffix = '_obs')


tminstr = '01-01-2024 00:00:00'
tmaxstr = '31-08-2026 00:00:00'

tminstr = '01-01-1980 00:00:00'
tmaxstr = '31-12-2026 00:00:00'


tminstr = tminstr
tmaxstr = tmaxstr

if 'lintrend' in md:
    md['lintrend']['tminstr'] = tminstr_targets
    md['lintrend']['tmaxstr'] = tmaxstr_targets


preprocessed = Preprocessed(heads = heads,stresses_dict = stresses_dict, time_step = time_step,
                            time_step_targets = time_step_targets, memory_dict = memory_dict,tminstr = tminstr, 
                            tmaxstr = tmaxstr, model_definition = md, forcast = True).preprocess() 
    
time = preprocessed.time
time_step = preprocessed.time_step
# heads = preprocessed.heads
time_step_targets = preprocessed.time_step_targets


simulation = Simulation(time = time, p_dict = p_dict, stresses_dict = stresses_dict, 
heads = heads, time_step = time_step, Nint_dict = Nint_dict, 
tminstr = tminstr, tmaxstr = tmaxstr, settings = settings, 
model_residuals = False, model_definition = md, is_a_control_run = False)
          

sim = simulation.simulate()

heads.simulated = sim


to_plot = [heads.observed,sim ]
legend_list = ['observations','simulation']
colors = ['b-','c-']
plotax = heads.plot_multiple_series(tseries_list = to_plot, colors = colors, legend_list = legend_list,tminstr = tminstr, tmaxstr = tmaxstr, share_axes = True )


path_to_extrafile_1 = os.path.join(path_to_resources_folder,'pb_atkb.txt')
# path_to_extrafile_1 = os.path.join(path_to_resources_folder,'S-BR1.01_20260917165935.csv')
# path_to_extrafile_2 = os.path.join(path_to_resources_folder,'S-BR2.03_20260917165935.csv')
# path_to_extrafile_3 = os.path.join(path_to_resources_folder,'S-BR6.04_20260917165937.csv')

# extra_time_serie = Heads.read_from_DINO(path_to_extrafile, tminstr = tminstr, tmaxstr = tmaxstr)
# extra_time_serie_1= Heads.read_from_DINO(path_to_extrafile_1, tminstr = tminstr, tmaxstr = tmaxstr, artdiverformat = True )
extra_time_serie_1= Heads.read_from_csv(path_to_extrafile_1, tminstr = tminstr, tmaxstr = tmaxstr, datetimeformat = '%Y-%m-%d %H:%M:%S' )

# extra_time_serie_2= Heads.read_from_DINO(path_to_extrafile_2, tminstr = tminstr, tmaxstr = tmaxstr, artdiverformat = True )
# extra_time_serie_3= Heads.read_from_DINO(path_to_extrafile_3, tminstr = tminstr, tmaxstr = tmaxstr, artdiverformat = True )



to_plot = [heads.observed,sim,extra_time_serie_1.observed ]
# to_plot = [heads.observed,sim,extra_time_serie_1.observed ,extra_time_serie_2.observed,extra_time_serie_3.observed]
legend_list = ['observations','simulation','pb_ATKB']
# legend_list = ['observations','simulation','S-BR1-01','S-BR2.03','S-BR6.04']
colors = ['b-','c-','m-','y-','g-']
plotax = heads.plot_multiple_series(tseries_list = to_plot, colors = colors, legend_list = legend_list,tminstr = tminstr, tmaxstr = tmaxstr, share_axes = True )

    
export_path = export.export_heads_nparray_to_ascii(sim, metadata = heads.metadata,suffix = '_sim')









# #################################################################################
# # path_to_extra_file = os.path.join(path_to_resources_folder,'S-BR1-01_sim.txt')
# basename = 'S-BR1-01'
# # extra_time_series_name = basename+'_sim.txt'
# # path_to_extrafile = os.path.join(curdir,extra_time_series_name) 
# # extra_time_serie = Heads.read_from_csv(path_to_extrafile, tminstr = tminstr, tmaxstr = tmaxstr)

# # extra_time_series_name = basename+'.txt'
# # path_to_extrafile = os.path.join(curdir,extra_time_series_name) 

# path_to_extrafile = os.path.join(path_to_resources_folder,'S-BR1.01_20260917165935.csv')
# path_to_file = os.path.join(path_to_resources_folder,'S-BR2.03_20260917165935.csv')
# # path_to_file = os.path.join(path_to_resources_folder,'S-BR6.04_20260917165937.csv')

# # extra_time_serie = Heads.read_from_DINO(path_to_extrafile, tminstr = tminstr, tmaxstr = tmaxstr)
# extra_time_serie = Heads.read_from_DINO(path_to_extrafile, tminstr = tminstr, tmaxstr = tmaxstr, artdiverformat = True )

# to_plot = [heads.observed,sim,extra_time_serie.observed ]
# legend_list = ['observations','simulation','S-BR1-01']
# legend_list = ['observations','simulation','S-BR6.04']
# colors = ['b-','c-','m-']
# plotax = heads.plot_multiple_series(tseries_list = to_plot, colors = colors, legend_list = legend_list,tminstr = tminstr, tmaxstr = tmaxstr, share_axes = True )

    
# export_path = export.export_heads_nparray_to_ascii(sim, metadata = heads.metadata,suffix = '_sim')

# #####################################################################################

















############################################################################### 

# tminstr = '01-01-1980 00:00:00'
# tmaxstr = '31-12-2026 00:00:00'

# tminstr = tminstr
# tmaxstr = tmaxstr

# # preprocessed = Preprocessed(heads = heads,stresses_dict = stresses_dict, time_step = time_step,
# #                             time_step_targets = time_step_targets, memory_dict = memory_dict,tminstr = tminstr, 
# #                             tmaxstr = tmaxstr, model_definition = md, forcast = True).preprocess() 
    
# # time = preprocessed.time
# # time_step = preprocessed.time_step
# # heads = preprocessed.heads
# # time_step_targets = preprocessed.time_step_targets

# simulation = Simulation(time = time, p_dict = p_dict, stresses_dict = stresses_dict, 
# heads = heads, time_step = time_step, Nint_dict = Nint_dict, 
# tminstr = tminstr, tmaxstr = tmaxstr, settings = settings, 
# model_residuals = False, model_definition = md)
          

# sim = simulation.simulate()


# # to_plot = [heads.interpolated, heads.modeled, riv.interpolated,sim]
# # legend_list = ['interpolated observations', 'modeled','riv','sim']


# to_plot = [heads.interpolated, heads.modeled,sim]
# legend_list = ['interpolated observations', 'modeled','sim']

# plotax = heads.plot_multiple_series(tseries_list = to_plot, legend_list = legend_list, share_axes = True )

    
# export_path = export.export_heads_nparray_to_ascii(sim, metadata = heads.metadata)

# p_prec_indexes = p_dict['p_indexes']['prec']['regime_1']['basicparam']
# logtrans = np.array(p_dict['logtransform'])[p_prec_indexes]
# p_prec = popt[p_prec_indexes]
# ig = unitresp.IncompleteGamma(stress_type='prec',Nint=100, time_step=1)
# ig.plotblockfunction(p_prec,logtrans)


# # heads = Heads.read_from_csv(path_to_file, tminstr = tminstr, tmaxstr = tmaxstr)
# heads_sim = Heads.read_from_csv(export_path)
# heads_sim.plot()

# glg,gg,ghg = heads_sim.generate_gxg()

# print('251 glg,gg,ghg',glg,gg,ghg)



# # sys.path.remove('D:\\dev\\_my_own_program\\src')
# sys.path.remove(path_to_src_folder)

# # check if sys,path is cleaned
# # print('182 sys.path', sys.path)
# # input()
