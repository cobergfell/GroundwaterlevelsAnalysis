# -*- coding: utf-8 -*-
"""
Created on Sun Jan 22 22:44:22 2023

@author: christophe
"""
import numpy as np
import os
import json
from collections import OrderedDict
from  unitresp import IncompleteGamma, Hantush, FloodWaveModel2L
from  trends import DrainageBase
from  modelnoise import ResidualsDecayExponentially
from modeldefinition import ModelDefinition
import logging
from stresstimeseries import Stress
from weightingfunctions import Sigmoid
from rootzone import Vangenuchten
#from abc import ABC, abstractmethod




class ParametersLogistic:
    
    """
    Class that provides objects defining the model parameters, their indexes and choices
    such as variable or fixed, logtransformed or not, such as determined by the user in the 
    modeldefinition class.
    
    Attributes
    ----------

    model_definition: python object of data type 'dict'
        model_definition defines the modeling options and as such is the first 
        object considered in the modeling process

    
    stresses_dict : python object of data type 'dict'
        _p_stresses_dict is a local copy of the stresses_dict dictionary containing all the entered
        Stress objects used to explain the observed groundwater levels variations.
        
    heads : Heads object or list of Heads objects
        Single Heads object or list of Heads objects containing the observed groundwater levels,
        associated metadata and preprocessed versions of the groundwater level time series 
    
    
    settings : python object of data type 'dict' (optional)
        An optional dictionnary of settings - not used in present version
    
    
    Methods
    --------
    assemble_p()
        Method applied to construct the model_definition object
        
    initialize_db(p_dict, stresses_dict,heads)    
        Method applied to find an appropriate initial value for the drainage base
        depending on the available observations
    
    read_json_parameters_from_text_files(filepath=None)
        Method to read p_dict from a jason file
    
    write_json_parameters_to_text_files(dirpath=None)
        Method to write p_dict from a jason file

    """
    
    
    
        

    
    
    def __init__(self, model_definition, stresses_dict = None, heads = None, settings = None):
        
        self.model_definition = model_definition
        self.stresses_dict = stresses_dict
        self.heads = heads
        self.settings = settings
        self.p_dict = self.assemble_p()
    
    def assemble_p(self):
        
        """ 
        This method assembles the parameters vector based on the given modeling settings
        
        Parameters
        ----------
        model_definition: python object of data type 'dict'
            dictionary of modeling options

            
        Returns
        -------
        p_dict = {}     python object of data type 'dict'
                        dictionary containing parameters names, initial values, optimized values,
                        whether parameetrs are logtansformed or native and whether parameters are vatying of fixed
        
        p :             list
                        list of parameters optimized values  
        
        pnames :        list
                        list of parameters names
        
        
        p_init :        list
                        list of parameters initial values   
        
        logtransform:   list
                        list of boolean values indicating if a parameter is logtransformed
                        (logtransform = True) of not (logtransform = False)        
        
        
        isvariable :    list
                        list of boolean values indicating if a parameter is optimized 
                        (isvariable = True) of fixed (isvariable = Fale)
        
        p_indexes :     python object of data type 'dict'
                        dictionary with stress names as entries, associating stress 
                        names with list of parameters indexes which are
                        the indexes of the parameters in the time series model 
                        parameters vector
                        
                        
        componentsnames : python object of data type 'dict'
                        dictionary with model components names (contstants, trends and responses to stresses) 
        
        
        
        """

        contribution_types = {}
        contribution_types['trends'] = ['db']
        contribution_types['root_zone'] = ['root_zone']
        contribution_types['stresses'] = ['prec','evap','pump','riv']
        contribution_types['residuals_model'] = ['noise']
        
        
        
        model_definition = self.model_definition
                        
        p_dict = {}
        pnames = []
        pnames_roots = []
        p_init = []
        p_min = []
        p_max = []
        logtransform = []
        isvariable = []
        p_indexes = {}
        componentsnames = {}
        weightingfunctionsnames = {}
        number_of_regimes = {}
        
        # Assemble model parameters excluding weighting functions (if applicable)
        for contribution_type in contribution_types:
            for component in contribution_types[contribution_type]:
                if component in model_definition:
                    if model_definition[component] != {}:
                    

                        p_indexes[component] = {}
                        if 'number_of_regimes' in model_definition[component]:
                            number_of_regimes[component] = model_definition[component]['number_of_regimes']

                        else:
                            number_of_regimes[component] = 1

                        for i in range(0,number_of_regimes[component]):
                           
                            keyname = 'regime_'+str(i+1)
                            p_indexes[component][keyname] = {}
                            p_indexes[component][keyname]['basicparam'] = [] # impulse response function parameters
                            p_indexes[component][keyname]['weightingparam'] = []
                            # default suffix
                            if  component in ['prec','evap'] and model_definition['evap']['funct_type'] == 'incomplete_gamma_f':
                                #p_indexes[component] = p_indexes['prec'] #not so elegant,trick to share parameters between prec and evap
                                suffix = ''
                            elif  component in ['pump'] and model_definition[component]['funct_type'] == 'Hantush':
                                suffix = ''
                                                      
                            else:
                                suffix = component
                                
                            if i > 0:
                                suffix += '_'+str(i+1)
                                
    
                            if component == 'db':
                                
                                if model_definition[component]['funct_type'] == 'constant':
                                    # parameters_definition = DrainageBase(trend_type = 'constant').get_initial_parameters() 
                                    parameters_definition = {}
                                    parameters_definition['db'] = {}
                                    parameters_definition['db']['pname'] = 'db'
                                    parameters_definition['db']['minvalue'] = 1.e-4
                                    parameters_definition['db']['maxvalue'] = 1.e4
                                    parameters_definition['db']['initvalue'] = 0.5 #1.e-3
                                    parameters_definition['db']['logtransform'] = False
                                    parameters_definition['db']['isvariable'] = True
                                    
                                    if model_definition['db']['fixed'] == True:
                                        parameters_definition['db']['isvariable'] = False
                                        parameters_definition['db']['logtransform'] = False
                                        parameters_definition['db']['initvalue'] = model_definition['db']['user_entered_value']
                                        
                                    if model_definition['db']['mean_observed_heads_minus_mean_response'] == True:
                                        parameters_definition['db']['isvariable'] = False   
                                        parameters_definition['db']['logtransform'] = False
                                    
                                else:
                                    parameters_definition = {}

                            elif component == 'root_zone':
                                if model_definition['root_zone']['apply_root_zone'] == True:
                                    if model_definition[component]['funct_type'] == 'vangenuchten':
                                        parameters_definition = Vangenuchten(root_zone_model = component).get_initial_parameters() 
                                else:
                                    parameters_definition = {}
    
    
    
                            elif component == 'prec':           
                                
                                if model_definition[component]['funct_type'] == 'incomplete_gamma':
                                    componentsnames[component] = 'IncompleteGamma'
                                    ig = IncompleteGamma(stress_type = component)
                                    parameters_definition = ig.get_initial_parameters(suffix = suffix, regime = keyname)     
                                    shares_parameters_with = {}
                                    
                                    if 'prec' in model_definition['constrain_with_harmonics']:
                                        parameters_definition['n']['isvariable'] = False
                                        
                                elif model_definition[component]['funct_type'] == 'floodwavemodel2L_typeI':
                                    componentsnames[component] = 'FloodWaveModel2L'
                                    floodwavemodel2L = FloodWaveModel2L(stress_type = component)
                                    parameters_definition = floodwavemodel2L.get_initial_parameters(suffix = suffix, regime = keyname)     
                                    shares_parameters_with = {}
                                else: 
                                    parameters_definition = {}  
                                                             
    
                            elif component == 'evap':   
                                if model_definition['root_zone']['apply_root_zone'] == False:
                                    if model_definition[component]['funct_type'] == 'incomplete_gamma_f':
                                        componentsnames[component] = 'IncompleteGamma'
                                        ig = IncompleteGamma(stress_type = component)
                                        parameters_definition = ig.get_initial_parameters(suffix = suffix, regime = keyname)  
                                        if 'prec' in model_definition['constrain_with_harmonics']:
                                            parameters_definition['f']['isvariable'] = False

                                        if keyname == 'regime_1':
                                            p_indexes[component][keyname]['basicparam'] = p_indexes['prec'][keyname]['basicparam'][:3]
                                        
                                        elif keyname == 'regime_2':
                                            p_indexes[component][keyname]['basicparam'] = p_indexes['prec'][keyname]['basicparam'][:3]
                                            p_indexes[component][keyname]['basicparam'].append(p_indexes['evap']['regime_1']['basicparam'][-1])
                                        else:
                                            pass # ignore more than 2 regimes for the moment
       
                                    elif model_definition[component]['funct_type'] == 'floodwavemodel2L_typeI':
                                        componentsnames[component] = 'FloodWaveModel2L'
                                        floodwavemodel2L = FloodWaveModel2L(stress_type = component)
                                        parameters_definition = floodwavemodel2L.get_initial_parameters(suffix = suffix, regime = keyname)     
                                        shares_parameters_with = {}  
                                        
                                        if model_definition['prec']['funct_type'] == 'floodwavemodel2L_typeI':
                                            shares_parameters_with = {"prec":[0,1,2,3,4]} # indices here refer to the indices of indices stored in p_indexes["prec"]
                                            for stress in shares_parameters_with:
                                                for j in range(0, len(shares_parameters_with[stress])):
                                                    keyname2 = 'regime_'+str(i+1)
                                                    p_indexes[component][keyname]['basicparam'].append(p_indexes[stress][keyname2]['basicparam'][j])
                               
                                else: 
                                    parameters_definition = {}
                                
                                
            
                            elif component == 'pump':           
                                if model_definition[component]['funct_type'] == 'incomplete_gamma':
                                    componentsnames[component] = 'IncompleteGamma'
                                    ig = IncompleteGamma(stress_type = component)
                                    parameters_definition = ig.get_initial_parameters(suffix = suffix)    
                                    shares_parameters_with = {}    
                                elif model_definition[component]['funct_type'] == 'Hantush':
                                    componentsnames[component] = 'Hantush'
                                    hantush = Hantush(stress_type = component)
                                    parameters_definition = hantush.get_initial_parameters(suffix = suffix)  
                                    shares_parameters_with = {}    
                                    
                                else: 
                                    parameters_definition = {}
                                                                    
                                    
    
                            elif component == 'riv':           
                                if model_definition[component]['funct_type'] == 'incomplete_gamma':
                                    componentsnames[component] = 'IncompleteGamma'
                                    ig = IncompleteGamma(stress_type = component)
                                    parameters_definition = ig.get_initial_parameters(suffix = suffix)    
                                    shares_parameters_with = {}    
                                    
                                elif model_definition[component]['funct_type'] == 'floodwavemodel2L_typeI':
                                    componentsnames[component] = 'FloodWaveModel2L'
                                    floodwavemodel2L = FloodWaveModel2L(stress_type = component)
                                    parameters_definition = floodwavemodel2L.get_initial_parameters(suffix = suffix)  
                                    
                                    if model_definition['prec']['funct_type'] == 'floodwavemodel2L_typeI':
                                        shares_parameters_with = {"prec":[0,1,2,3,4]} # indices here refer to the indices of indices stored in p_indexes["prec"]
                                        for stress in shares_parameters_with:
                                            for j in range(0, len(shares_parameters_with[stress])):
                                                keyname2 = 'regime_'+str(i+1)
                                                p_indexes[component][keyname]['basicparam'].append(p_indexes[stress][keyname2]['basicparam'][j])
                                                                           
                                else: 
                                    parameters_definition = {}                            
                                
                          
                                
                            elif component == 'noise':   
                                if model_definition[component]['funct_type'] == 'residuals_decay_exponentially':
                                    componentsnames[component] = 'ResidualsDecayExponentially'
                                    parameters_definition = ResidualsDecayExponentially().get_initial_parameters() 
                                else: 
                                    parameters_definition = {} 
                                    
                            else:
                                parameters_definition = {}                                
                                
                                
                                                                       

                            for pname in parameters_definition:

                                try:
                                    splitted = pname.split('_')
                                    pname_root = splitted[0]
                                    suffix = splitted[1]

                                except:
                                    suffix = None
       
                                condition_1 = True
                                
                                if component == 'evap':
                                    if model_definition[component]['funct_type'] == 'incomplete_gamma_f':
   
                                        if pname_root in pnames_roots:
                                            
                                            condition_1 = False
                                            
                                    elif model_definition[component]['funct_type'] == 'floodwavemodel2L_typeI':
                                        if pname_root in pnames_roots:
                                            condition_1 = False        
                                            
                                if component == 'riv':
                                            
                                    if model_definition[component]['funct_type'] == 'floodwavemodel2L_typeI':
                                        if pname_root in pnames_roots:
                                            condition_1 = False                                         

                                if condition_1 == True:
                                    parameters_sub_definition = parameters_definition[pname]
                                    pnames_roots.append(pname_root)
                                    pnames.append(pname)
                                    logtrans = parameters_sub_definition["logtransform"]
                                    p_value = parameters_sub_definition["initvalue"]
                                    p_min_value = parameters_sub_definition["minvalue"]
                                    p_max_value = parameters_sub_definition["maxvalue"]
                                    if logtrans == True:
                                        if p_value == 0:
                                            p_value = 1e-10
                                        if p_min_value == 0:
                                            p_min_value = 1e-10    
                                        if p_max_value == 0:
                                            p_max_value = 1e-10                                 
                                        p_value = np.log(p_value)
                                        p_min_value = np.log(p_min_value)
                                        p_max_value = np.log(p_max_value)
                                    p_init.append(p_value)
                                    p_min.append(p_min_value)
                                    p_max.append(p_max_value)
                                    logtransform.append(logtrans)
                                    isvariable.append(parameters_sub_definition["isvariable"])
                                    
                                    p_indexes[component][keyname]['basicparam'].append(len(p_init)-1)

                            if i > 0: # if stress responses are superposed  

                                condition_2 = True
    
                                weightingfunctionname = model_definition[component]['weighting']
                                if weightingfunctionname in ['prec','evap','evap']: # in which case the weights are shared with one of these stresses 
                                    condition_2 == False

                                weightingfunctionsnames[component] = weightingfunctionname
                                if weightingfunctionname == 'sigmoid': 
                                    parameters_definition = Sigmoid().get_initial_parameters(suffix = suffix) 
            
    
                                if (condition_1 == True) and (condition_2 == True):
                                    for pname in parameters_definition:
                                        parameters_sub_definition = parameters_definition[pname]
                                        pnames_roots.append(pname_root)
                                        pnames.append(pname)
                                        logtrans = parameters_sub_definition["logtransform"]
                                        p_value = parameters_sub_definition["initvalue"]
                                        p_min_value = parameters_sub_definition["minvalue"]
                                        p_max_value = parameters_sub_definition["maxvalue"]
                                        if logtrans == True:
                                            if p_value == 0:
                                                p_value = 1e-10
                                            if p_min_value == 0:
                                                p_min_value = 1e-10    
                                            if p_max_value == 0:
                                                p_max_value = 1e-10                                 
                                            p_value = np.log(p_value)
                                            p_min_value = np.log(p_min_value)
                                            p_max_value = np.log(p_max_value)
                                        p_init.append(p_value)
                                        p_min.append(p_min_value)
                                        p_max.append(p_max_value)
                                        logtransform.append(logtrans)
                                        isvariable.append(parameters_sub_definition["isvariable"])
                                        p_indexes[component][keyname]['weightingparam'].append(len(p_init)-1)   
                               

        p_dict["p_indexes"] = p_indexes 
        p_dict["p"] = p_init
        p_dict["p_min"] = p_min
        p_dict["p_max"] = p_max
        p_dict["pnames"] = pnames
        p_dict["p_init"] = p_init
        p_dict["logtransform"] = logtransform
        p_dict["isvariable"] = isvariable
        p_dict["componentsnames"] = componentsnames
        p_dict["weightingfunctionsnames"] = weightingfunctionsnames
        p_dict["number_of_regimes"] = number_of_regimes
        
        if model_definition['db']['mean_observed_heads_minus_mean_response'] == False:
            if model_definition['db']['fixed'] == False:
                if (self.stresses_dict is not None) & (self.heads is not None):
                    p_dict = self.initialize_db(p_dict, self.stresses_dict, self.heads)
                self.p_dict = p_dict
                
        # print('447 check p_dict',p_dict)
        # input()
                
        
        return p_dict
    


    def initialize_db(self,p_dict, stresses_dict,heads):
        

        model_definition = self.model_definition

        if 'db' in p_dict['p_indexes']:
            if model_definition['db']['funct_type'] == 'constant':
                index_db = p_dict['p_indexes']['db']['regime_1']['basicparam'][0]
                if model_definition['db']['initial_value_from'] == 'observed_heads_median':  
                    db = np.median(heads.observed[:,1])
                    
                elif 'user_entered_value' in model_definition['db']:
                    db = model_definition['db']['user_entered_value']
                else:
                    db = 5. #default value
                    
                for stress_type in ['prec','evap','pump','riv']:
                    if stress_type in model_definition:
                        if 'use_normalized_time_series' in model_definition[stress_type] :
                            if model_definition[stress_type]['use_normalized_time_series'] == True :
                                if stress_type in stresses_dict:
                                    for key in stresses_dict[stress_type]:
                                        stress = stresses_dict[stress_type][key].interpolated
                                        
                                        if stress_type == 'prec':
                                            scale_factor = 1000.
                                        elif stress_type == 'evap':
                                            scale_factor = -1000.
                                        elif stress_type == 'pump':
                                            scale_factor = -1e-3    
                                        elif stress_type == 'riv':
                                            scale_factor = 1.0
            
                                        db += np.mean(stress[:,1])*scale_factor

                
                db_min = db - 2.
                db_max = db + 1.
                if p_dict['logtransform'][index_db] == True:
                    if db_min <= 0:
                        p_dict['logtransform'][index_db] = False
                
                if p_dict['logtransform'][index_db] == True:
                    db = np.log(db)
                    db_min = np.log(db_min)
                    db_max = np.log(db_max)
    
                p_dict['p'][index_db] = db
                p_dict['p_min'][index_db] = db_min
                p_dict['p_max'][index_db] = db_max
                if model_definition['db']['fixed'] == True:
                    p_dict['isvariable'][index_db] = False
                else:
                    p_dict['isvariable'][index_db] = True
            return p_dict
                

        

    def read_json_parameters_from_text_files(self,filepath=None):
        """ TO DO: implement this function to read parameters from json text files"""
        pass

    def write_json_parameters_to_text_files(self,dirpath=None):
        """ TO DO: implement this function to write parameters to json text files """
        pass
    
    

                                  
if __name__ == "__main__":
    md = ModelDefinition()
    parametersLogistic = ParametersLogistic(md.model_definition)
    p_dict = parametersLogistic.assemble_p()
    print('426 p_dict', p_dict)    
    