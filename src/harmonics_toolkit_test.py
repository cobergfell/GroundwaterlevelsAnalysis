# -*- coding: utf-8 -*-
"""
Created on Fri Jul  7 15:02:49 2023

@author: christophe
"""
import os
from logging import getLogger
import numpy as np
from stresstimeseries import Stress
from headstimeseries import Heads
from parameterslogistic import ParametersLogistic
from modeldefinition import ModelDefinition
from preprocessedseries import Preprocessed
from utilities import basename
from scipy.optimize import least_squares, fsolve

logger = getLogger(__name__)

"""
This module bundles functions that are needed to constrain the optimisation of the 
parameters by requiring the consistency of observed and modeled seasonal harmonics

Reference
---------

Obergfell, C., M. Bakker, and K. Maas (2019), Estimation of Average Diffuse 
Aquifer Recharge Using Time Series Modeling of Groundwater Heads, 
Water Resources Research, 55(3), 2194-2210. https://doi.org/10.1029/2018WR024235
   
    
"""







def get_index(string,Pnames):

    """
    A commodity function to return the index of a parameter in the parameters vector given
    the index name
        
    """
    if string in Pnames:
        for i in range(0,len(Pnames)):
            if Pnames[i] == string:
                ind=i 
    return ind
                                
def func_Kees1(w,a,n):
    
    """
    A function to calculate the amplitude term
    in equation A9 in aboven mentionned paper
    
    Parameters
    ----------
    w: float
        angular frequency
        
    a: float
        the first parameter of the incomplete gamma function following the 
        notation in the module unitresp.py
        
    n: float
        the second parameter of the incomplete gamma function following the 
        notation in the module unitresp.py      
        
    """    
    term1 = 1+pow(w/a,2)
    
    return pow(term1,n/2) 
    
def func_Kees2(w,a,n):     

    """
    A function to calculate the phase term
    in equation A9 in above mentionned paper
    
    Parameters
    ----------
    see func_Kees1
        
    """      
    
    return n*np.arctan(w/a)/w       
    

def non_linear_equations1(trial_solution,w,A,a,p_dict, stresses_dict, heads):

    """
    A function to constitute the non linear equations
    needed to infer parameters n and f
    in equation A9 in above mentionned paper
    
    Parameters
    ----------
    trial_solution: numpy array_like
        the initial guess when solving the set of lon-linear equations

    w: float
        angular frequency
        
    A: float
        the scale parameter of the incomplete gamma function following the 
        notation in the module unitresp.py        
        
    a: float
        the first parameter of the incomplete gamma function following the 
        notation in the module unitresp.py
        
    n: float
        the second parameter of the incomplete gamma function following the 
        notation in the module unitresp.py         
    
    p_dict : python object of data type 'dict'
        Dictionary defining the model parameters
        and associated specifications such as associated indexes, names, logtransformed or not, 
        variable or fixed etc. See the the parameterslogistic.py module for additional details.
        
    stresses_dict : python object of data type 'dict'
        Dictionary containing all the entered Stress objects used to explain 
        the observed groundwater levels variations.
        
    heads : Heads object or list of Heads objects
        Single Heads object or list of Heads objects containing the observed groundwater levels,
        associated metadata and preprocessed versions of the groundwater level time series 
        
    """ 
 
    #local copy of mutable input dictionaries
    _p_dict = {}
    for key in p_dict:
        _p_dict[key] = p_dict[key]

    _stresses_dict = {}
    for key in stresses_dict:
        _stresses_dict[key] = stresses_dict[key]
    

    evap_tseries_name = list(_stresses_dict['evap'].keys())[0] # for the moment, in case of multiple time series per stress, pick the first one 
    prec_tseries_name = list(_stresses_dict['prec'].keys())[0] # for the moment, in case of multiple time series per stress, pick the first one 


    # try:
    p = _p_dict['p']
    p_indexes = _p_dict['p_indexes']
    logtransform = _p_dict['logtransform']
    ind_n = p_indexes['evap']['regime_1']['basicparam'][2]
    ind_f = p_indexes['evap']['regime_1']['basicparam'][3]
        
   
    # n = p[ind_n]
    # f = p[ind_f]
    
    n = trial_solution[0]
    f = trial_solution[1]
    
    if logtransform[ind_n] == True:
         n = np.exp(n)
    if logtransform[ind_f] == True:
         f = np.exp(f)    
    
    Ampl_h = heads.harmonic_observed.amplitude[0]
    Ampl_evap = stresses_dict['evap'][evap_tseries_name].harmonic_observed.amplitude[0]
    Ampl_prec = stresses_dict['prec'][prec_tseries_name].harmonic_observed.amplitude[0]
    phi_h = heads.harmonic_observed.phase_shift[0]
    phi_evap = stresses_dict['evap'][evap_tseries_name].harmonic_observed.phase_shift[0]
    phi_prec = stresses_dict['prec'][prec_tseries_name].harmonic_observed.phase_shift[0]    
     
    damping_term_Kees = func_Kees1(w,a,n)
    timelag_Kees = func_Kees2(w,a,n)
    
    #define equation 1
    # term11=sin(w*(phi_h-phi_prec-timelag_Kees))
    # term12=f*A*Ampl_evap*sin(w*(phi_evap-phi_prec))/(Ampl_h*damping_term_Kees)
    term11 = np.tan(w*(phi_h-phi_prec-timelag_Kees))
    term12 = f*Ampl_evap*np.sin(w*(phi_evap-phi_prec))/(Ampl_prec+f*Ampl_evap*np.cos(w*(phi_evap-phi_prec)))
    
    eq1 = term11-term12

    #define equation 2
    term21 = pow(Ampl_h*damping_term_Kees/A,2)
    term22 = pow(f*Ampl_evap,2)+pow(Ampl_prec,2)+2*f*Ampl_evap*Ampl_prec*np.cos(w*(phi_evap-phi_prec))
    eq2 = term21-term22 
    return eq1,eq2
    
    
    

def non_linear_equations2(trial_solution,p_dict, stresses_dict, heads):
    """
    Another function to constitute the non linear equations
    needed to infer parameters n and f
    in equation A9 in above mentionned paper. This version includes the
    effect of pumping.
    
    Parameters
    ----------
    trial_solution: numpy array_like
        the initial guess when solving the set of lon-linear equations

    w: float
        angular frequency
        
    A: float
        the scale parameter of the incomplete gamma function following the 
        notation in the module unitresp.py        
        
    a: float
        the first parameter of the incomplete gamma function following the 
        notation in the module unitresp.py
        
    n: float
        the second parameter of the incomplete gamma function following the 
        notation in the module unitresp.py         
    
    p_dict : python object of data type 'dict'
        Dictionary defining the model parameters
        and associated specifications such as associated indexes, names, logtransformed or not, 
        variable or fixed etc. See the the parameterslogistic.py module for additional details.
        
    stresses_dict : python object of data type 'dict'
        Dictionary containing all the entered Stress objects used to explain 
        the observed groundwater levels variations.
        
    heads : Heads object or list of Heads objects
        Single Heads object or list of Heads objects containing the observed groundwater levels,
        associated metadata and preprocessed versions of the groundwater level time series 
        
    """ 
 
    #local copy of mutable input dictionaries
    _p_dict = {}
    for key in p_dict:
        _p_dict[key] = p_dict[key]

    _stresses_dict = {}
    for key in stresses_dict:
        _stresses_dict[key] = stresses_dict[key]
    
        
    w = 2*np.pi/365.25
    evap_tseries_name = list(_stresses_dict['evap'].keys())[0] # for the moment, in case of multiple time series per stress, pick the first one 
    prec_tseries_name = list(_stresses_dict['prec'].keys())[0] # for the moment, in case of multiple time series per stress, pick the first one 

    
    # try:
    p = _p_dict['p']
    p_indexes = _p_dict['p_indexes']
    logtransform = _p_dict['logtransform']
    
    ind_A = p_indexes['prec']['regime_1']['basicparam'][0]
    ind_a = p_indexes['prec']['regime_1']['basicparam'][1]

    A = p[ind_A]
    a = p[ind_a]
    
    if logtransform[ind_A] == True:
          A = np.exp(A)
    if logtransform[ind_a] == True:
          a = np.exp(a)  

    ind_n = p_indexes['evap']['regime_1']['basicparam'][2]
    ind_f = p_indexes['evap']['regime_1']['basicparam'][3]
        
    n = trial_solution[0]
    f = trial_solution[1]
    

    if logtransform[ind_n] == True:
          n = np.exp(n)
    if logtransform[ind_f] == True:
          f = np.exp(f)   

          
    # Ampl_h = heads.harmonic_observed.amplitude[0]
    # Ampl_evap = stresses_dict['evap'][evap_tseries_name].harmonic_observed.amplitude[0]
    # Ampl_prec = stresses_dict['prec'][prec_tseries_name].harmonic_observed.amplitude[0]
    # phi_h = heads.harmonic_observed.phase_shift[0]
    # phi_evap = stresses_dict['evap'][evap_tseries_name].harmonic_observed.phase_shift[0]
    # phi_prec = stresses_dict['prec'][prec_tseries_name].harmonic_observed.phase_shift[0]  


    Ampl_h = heads.harmonic_interpolated.amplitude[0]
    Ampl_evap = stresses_dict['evap'][evap_tseries_name].harmonic_interpolated.amplitude[0]
    Ampl_prec = stresses_dict['prec'][prec_tseries_name].harmonic_interpolated.amplitude[0]
    phi_h = heads.harmonic_interpolated.phase_shift[0]
    phi_evap = stresses_dict['evap'][evap_tseries_name].harmonic_interpolated.phase_shift[0]
    phi_prec = stresses_dict['prec'][prec_tseries_name].harmonic_interpolated.phase_shift[0]  
         
    D_prec = func_Kees1(w,a,n) #damping term Kees Maas PREC and EVAP
    
    TL_prec = func_Kees2(w,a,n)  #time lag Kees Maas PREC and EVAP
    
    
    from datetime import datetime
    from matplotlib.dates import date2num, num2date
    
    
    hc = heads.harmonic_interpolated
    timenum_begin = hc.harmonic_component[0,0]
    timenum_begin = hc.harmonic_component[0,0]
    timenum_end = hc.harmonic_component[-1,0]
    date_begin= num2date(timenum_begin)
    date_end= num2date(timenum_end)
    print('345 heads.harmonic_interpolated date_begin',date_begin)
    print('346 heads._harmonic_interpolated date_end',date_end)    
    print('347 Ampl_h',Ampl_h) 
    print('348 phi_h',phi_h) 
    
    
    hc = stresses_dict['prec'][prec_tseries_name].harmonic_interpolated
    timenum_begin = hc.harmonic_component[0,0]
    timenum_begin = hc.harmonic_component[0,0]
    timenum_end = hc.harmonic_component[-1,0]
    date_begin= num2date(timenum_begin)
    date_end= num2date(timenum_end)
    print('357 prec._harmonic_interpolated date_begin',date_begin)
    print('358 prec._harmonic_interpolated date_end',date_end)
    print('359 Ampl_prec',Ampl_prec) 
    print('360 phi_prec',phi_prec)     
    
    
    
    
    hc = stresses_dict['evap'][evap_tseries_name].harmonic_interpolated
    timenum_begin = hc.harmonic_component[0,0]
    timenum_begin = hc.harmonic_component[0,0]
    timenum_end = hc.harmonic_component[-1,0]
    date_begin= num2date(timenum_begin)
    date_end= num2date(timenum_end)
    print('371 evap._harmonic_interpolated date_begin',date_begin)
    print('372 evap._harmonic_interpolated date_end',date_end)
    print('373 Ampl_evap',Ampl_evap) 
    print('374 phi_evap',phi_evap)      
    input()    
    
          

    if 'pump' in stresses_dict:
        pump_tseries_name = list(stresses_dict['pump'].keys())[0] # for the moment, in case of multiple time series per stress, pick the first one 

        ind_A_pump = p_indexes['pump']['regime_1']['basicparam'][0]
        ind_a_pump = p_indexes['pump']['regime_1']['basicparam'][1]
        ind_n_pump = p_indexes['pump']['regime_1']['basicparam'][2]
        
        A_pump = p[ind_A_pump]
        a_pump = p[ind_a_pump]
        n_pump = p[ind_n_pump]
        
        
        if logtransform[ind_A_pump] == True:
            A_pump = np.exp(A_pump)
        if logtransform[ind_a_pump] == True:
            a_pump = np.exp(a_pump)  
        if logtransform[ind_n_pump] == True:
            n_pump = np.exp(n_pump)     
        
        
    
        D_pump=func_Kees1(w,a_pump,n_pump)
        TL_pump=func_Kees2(w,a_pump,n_pump)      
        

        Ampl_pump = stresses_dict['pump'][pump_tseries_name].harmonic.amplitude[0]
        phi_pump = stresses_dict['pump'][pump_tseries_name].harmonic.phase_shift[0] 
        

    
    equations = []    
    
    G1 = A*Ampl_prec/D_prec
    G2 = f*A*Ampl_evap/D_prec
    

    delta_PE = w*(phi_evap-phi_prec)
    Ape = np.sqrt(pow(G1,2)+pow(G2,2)+2*G1*G2*np.cos(delta_PE))


    argument1 = G2*np.sin(delta_PE)/(G1+G2*np.cos(delta_PE))

    #eq2=tan(w*(phi_h+TL_prec+phi_prec))-argument1
    
    Tpe = np.arctan(argument1)/w   #0001

    if ('pump' not in stresses_dict) == True:
        eq1 = Ampl_h-Ape #equality of amplitudes 
        equations.append(eq1)
        eq2 = phi_h-(Tpe+TL_prec+phi_prec)#equality of phase
        equations.append(eq2)
        
    else:   
        G3 = A_pump*Ampl_pump/D_pump
        delta_PEW = w*(Tpe+TL_prec+phi_prec-TL_pump-phi_pump)
        Apew = np.sqrt(pow(G3,2)+pow(Ape,2)+2*G3*Ape*np.cos(delta_PEW))      
        eq1 = Ampl_h-Apew #equality of amplitudes        
        equations.append(eq1)
        argument2 = Ape*np.sin(delta_PEW)/(G3+Ape*np.cos(delta_PEW))
        Tpew = np.arctan(argument2)/w
        eq2 = phi_h-(Tpew+TL_pump+phi_pump)#equality of phase
        equations.append(eq2)
        #eq3=TL_prec-TL_pump-(phi_pump-phi_prec) THIS IS NOT CORRECT#relative time lag wel/prec or pump evap
        #equations.append(eq3)  
        
        
    return equations   



def update_P_given_harmonic_components(p_dict, stresses_dict, heads):

    """
    A function to update P given the solution to the system of non-linear 
    equations defined in non_linear_equations1 or non_linear_equations2.
    
    Parameters
    ----------

    
    p_dict : python object of data type 'dict'
        Dictionary defining the model parameters
        and associated specifications such as associated indexes, names, logtransformed or not, 
        variable or fixed etc. See the the parameterslogistic.py module for additional details.
        
    stresses_dict : python object of data type 'dict'
        Dictionary containing all the entered Stress objects used to explain 
        the observed groundwater levels variations.
        
    heads : Heads object or list of Heads objects
        Single Heads object or list of Heads objects containing the observed groundwater levels,
        associated metadata and preprocessed versions of the groundwater level time series 
        
    """ 
    

    dic_P_optimized = {}
    pnames = p_dict['pnames']
    p = p_dict['p']
    

    logtransform = p_dict['logtransform']
    parameter_names_list = ['A','a']
    

    for parameter_name in parameter_names_list: 
    
        ind=get_index(parameter_name,pnames)  
        parameter_value = p[ind] 
        if logtransform[ind] == True:
            parameter_value = np.exp(parameter_value)           
        dic_P_optimized[parameter_name] = parameter_value   

        
    trial_solution = []   
    # parameter_names_list = ['n_PIII_PREC','f_factor_EVAP']
    parameter_names_list = ['n','f']


    for parameter_name in parameter_names_list: 
        
        ind = get_index(parameter_name,pnames)  
        parameter_value = p[ind]             
        trial_solution.append(parameter_value)

    # solution = fsolve(non_linear_equations2,trial_solution,args=(p_dict, stresses_dict, heads))  #,full_output=1


    lb = np.array(p_dict['p_min'])
    ub = np.array(p_dict['p_max'])
    mask = np.zeros(len(lb),dtype = bool)
    
    p = p_dict['p']
    p_indexes = p_dict['p_indexes']
    
    ind_n = p_indexes['evap']['regime_1']['basicparam'][2]
    ind_f = p_indexes['evap']['regime_1']['basicparam'][3]

    mask[ind_n] = True
    mask[ind_f] = True
    lb = lb[mask]
    ub = ub[mask]
    bounds = (lb, ub)    


    outcome = least_squares(non_linear_equations2, x0 = trial_solution, bounds=bounds,
    args=(p_dict, stresses_dict, heads), kwargs={},method='trf', ftol=1e-08, xtol=1e-08, gtol=1e-08, 
    x_scale=1.0, loss='linear', f_scale=1.0, diff_step=None, tr_solver=None, tr_options={}, 
    jac_sparsity=None, max_nfev=None, verbose=0)

    solution = outcome['x']
    success = outcome['success']
    cost = outcome['cost']
    res = outcome['fun']
    jac = outcome['jac']
    message = outcome['message']
    
    
    for i in range(0,len(parameter_names_list)):
        parameter_name = parameter_names_list[i]
        ind = get_index(parameter_name,pnames)  
        parameter_value = solution[i] 

        # if logtransform[ind] == True:
        #     parameter_value = np.log(parameter_value)   
        p_dict['p'][ind] = parameter_value                    
      
        
    return p_dict



    

    
if __name__ == "__main__":
    #curdir=os.getcwd()
    abs_path = os.path.dirname(os.path.abspath(__file__))
    abs_path_splitted = abs_path.split('\\')
    path_to_parent_folder_elements = abs_path_splitted[:-1]
    path_to_parent_folder = os.path.join(*path_to_parent_folder_elements)
    splitted = path_to_parent_folder.split(':')#trick  to  repair  path
    path_to_parent_folder = splitted[0]+':\\'+splitted[1]#trick  to  repair  path
    path_to_file_folder = os.path.join(path_to_parent_folder,'resources')  

    stresses_dict = {}

    # path_to_file = os.path.join(path_to_file_folder,'260_De_Bilt_PREC_19010101_20200910.csv')
    path_to_file = os.path.join(path_to_file_folder,'672_PREC_Hellendoorn_19510101_20160731_corrected_for_interception2mm.csv')
    # path_to_file = os.path.join(path_to_file_folder,'prec_giersbergen.csv')

    key = basename(path_to_file)
    stresses_dict['prec'] = {}
    stresses_dict['prec'][key] = Stress.read_from_csv(csv_file = path_to_file, stress_type = 'prec',cumulative = True)    
    

    path_to_file=os.path.join(path_to_file_folder,'260_De_Bilt_EVAP_19010101_20200910.csv')
    # path_to_file=os.path.join(path_to_file_folder,'evap_eindhoven.csv')
    
    key = basename(path_to_file)
    stresses_dict['evap'] = {}
    stresses_dict['evap'][key] = Stress.read_from_csv(csv_file = path_to_file, stress_type = 'evap',cumulative = True)    
    
    

    
    tminstr = '01-01-1982 00:00:00'
    tmaxstr = '31-12-2005 00:00:00'
    
    path_to_file_folder = os.path.join(path_to_parent_folder,'resources')  
    path_to_file = os.path.join(path_to_file_folder,'28AP0093_1.txt')
    # heads = Heads.read_from_csv(path_to_file)
    heads = Heads.read_from_csv(path_to_file, tminstr = tminstr, tmaxstr = tmaxstr)
    
    # path_to_file_folder = os.path.join(path_to_parent_folder,'resources','selected_filters_paper1')  
    # list_of_piezometers_files_names = os.listdir(path_to_file_folder)
    # for filename in list_of_piezometers_files_names: 
    #     path_to_file = os.path.join(path_to_file_folder,filename)
    #     heads = Heads.read_from_csv(path_to_file, tminstr = tminstr, tmaxstr = tmaxstr) 
    #     #heads.plot()


    arguments_dict = {}

    for key in stresses_dict:
        for e in stresses_dict[key]:
            stresses_dict[key][e].plot()

    md = ModelDefinition().model_definition
    memory_dict = ModelDefinition().memory_dict
    parametersLogistic = ParametersLogistic(md)
    
    # time_step = 1./(24.*2)
    # time_step_targets = 1.5
    
    time_step = 1.
    time_step_targets = 1.
    
    preprocessed = Preprocessed(heads = heads,stresses_dict = stresses_dict, time_step = time_step,
                                time_step_targets = time_step_targets, memory_dict = memory_dict,tminstr = tminstr, 
                                tmaxstr = tmaxstr, model_definition = md).preprocess() 
        
    time = preprocessed.time
    time_step = preprocessed.time_step
    heads = preprocessed.heads

    stresses_dict = preprocessed.stresses_dict 
    Nint_dict = preprocessed.Nint_dict
    all_data_for_neural_networks = preprocessed.all_data_for_neural_networks
    all_targets_for_neural_networks = preprocessed.all_targets_for_neural_networks   
    

    parametersLogistic = ParametersLogistic(md)
    p_dict = parametersLogistic.assemble_p()

    p_dict = update_P_given_harmonic_components(p_dict, stresses_dict, heads)
    