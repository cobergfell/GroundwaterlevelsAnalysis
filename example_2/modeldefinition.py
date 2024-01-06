# -*- coding: utf-8 -*-
"""
Created on Sun Jan 22 22:44:22 2023

@author: christophe
"""
import numpy as np
import os
import json
from collections import OrderedDict
from logging import getLogger


logger = getLogger(__name__)




class ModelDefinition:
    
    """ 
    A class that produces an object defining the time series model. That is,
    it specifies which stresses are considered to explain the observed groundwater level 
    fluctuations with the corresponding reponse functions. It also defines the trends and constants
    to be added, whether the residuals are modeled as well, and how many groundwater regimes are simulated.
    
    Attributes
    ----------
    working_directory: string (optional)
        string defining the path to a working directory
    
    model_definition_from_json : boolean
        Boolean value indicating if the model definition is to be read from 
        a jason file

    
    Examples
    --------
    
    """    
    
    
    
    

    
    def __init__(self, working_directory  = None, read_from_resources = False):

        self._model_definition = None

        model_definition = {}
        model_definition['db'] = {}
        model_definition['db']['funct_type'] = 'constant'
        model_definition['db']['fixed'] = True # False
        model_definition['db']['equal_to_mean_observed_heads_minus_mean_response'] = False #True 
        model_definition['db']['default_initial_value'] = 0 # 11.43
        model_definition['db']['equal_to_observed_heads_median'] = False
        model_definition['db']['equal_to_mean_river_stage'] = False
        
        model_definition['constrain_with_harmonics'] = [] #['prec','evap'] # ['prec','evap']       
        model_definition['root_zone'] = {}
        model_definition['root_zone']['funct_type'] = 'vangenuchten'
        model_definition['root_zone']['apply_root_zone'] = False #True

        # model_definition['prec'] = {}
        # model_definition['prec']['funct_type'] = 'incomplete_gamma'
        # model_definition['prec']['number_of_regimes'] = 1 # if > 1, response will be weighted and superposed
        # model_definition['prec']['weighting'] = 'sigmoid'

        model_definition['prec'] = {}
        # model_definition['prec']['funct_type'] = 'incomplete_gamma' #'floodwavemodel2L_typeI' #
        model_definition['prec']['funct_type'] = 'floodwavemodel2L_typeI' #'incomplete_gamma' #
        model_definition['prec']['number_of_regimes'] = 1 # if > 1, response will be weighted and superposed
        model_definition['prec']['weighting'] = 'sigmoid'        
        model_definition['prec']['use_normalized_time_series'] = False #True #  if true then add mean stress to the drainage base

        
        if model_definition['root_zone']['apply_root_zone'] == True: # resolve incompatibility of options
            model_definition['constrain_with_harmonics'] = []
        
        model_definition['evap'] = {} 
        model_definition['evap']['funct_type'] = 'floodwavemodel2L_typeI' # 'incomplete_gamma_f' #
        model_definition['evap']['number_of_regimes'] = model_definition['prec']['number_of_regimes'] # for the moment prec and evap share the same number of regimes
        model_definition['evap']['weighting'] = 'prec' # this indicates that the weighing is similar to prec
        model_definition['evap']['use_normalized_time_series'] = False#True #  if true then add mean stress to the drainage base
 
        # model_definition['pump'] = {} 
        # # model_definition['pump']['funct_type'] = 'incomplete_gamma'
        # model_definition['pump']['funct_type'] = 'Hantush'
                    
        model_definition['riv'] = {} 
        # # model_definition['riv']['funct_type'] = 'incomplete_gamma'
        model_definition['riv']['funct_type'] = 'floodwavemodel2L_typeI'
        model_definition['riv']['use_normalized_time_series'] = True # if true then add mean stress to the drainage base
        
        model_definition['noise'] = {} 
        model_definition['noise']['funct_type'] = 'residuals_decay_exponentially'            
        model_definition['noise']['model_residuals'] = True #    #noise can also be switched off in simulations options

                
        
        if working_directory == None:
            #curdir=os.getcwd()
            abs_path = os.path.dirname(os.path.abspath(__file__))
            abs_path_splitted = abs_path.split('\\')
            working_directory_elements = abs_path_splitted[:-1]
            working_directory = os.path.join(*working_directory_elements)
            splitted = working_directory.split(':') # trick to repair path
            working_directory = splitted[0]+':\\'+splitted[1] # trick to repair path
    
        if read_from_resources == True: # overrules above defined model definition
            try:
                mydir = os.path.join(working_directory,'resources')
                filepath = os.path.join(mydir,'model_definition.json')
                
                fp = open(filepath,'r')
                model_definition = json.load(fp)
                fp.close()              
            except FileNotFoundError:
                message = (f'\n\nFilepath {filepath} does not exist\n')
                logger.warning(message)               

        memory_dict = {'prec':365,
                      'evap':365,
                      'riv':365,
                      'pump':100,
                      'noise':100}        
  
        self.working_directory = working_directory 
        self._model_definition = model_definition     
        self._memory_dict = memory_dict
            
          
                
    @property
    def model_definition(self):
        """getter for model_definition."""
        return self._model_definition

    @model_definition.setter
    def model_definition(self,model_definition):
        """setter for updated observed time series."""
        self._model_definition = model_definition       
        
        
        
    @property
    def memory_dict(self):
        """getter for memory_dict."""
        return self._memory_dict

    @memory_dict.setter
    def memory_dict(self,memory_dict):
        """setter for updated observed time series."""
        self._memory_dict = memory_dict           
        
        
        
        
    def to_json(self,dirpath=None):
        
        """ write parameters to json text files """

        # https://stackoverflow.com/questions/17043860/how-to-dump-a-dict-to-a-json-file
        import json
        
        if dirpath is None:
            with open('model_definition.json', 'w') as fp:
                dumped = json.dump(self._model_definition, fp, indent=4)
                
        else:

            if isinstance(dirpath,str):
                try:
                    filepath = os.path.join(dirpath,'model_definition.json')
                    with open(filepath,"w") as fp:
                        json.dump(self._model_definition, fp, indent=4)    
                except FileNotFoundError:
                            
                    clsname = str(self.__class__.__name__)
                    modulename = str(__name__)
        
                    message = (f'\nIn class {clsname} of module {modulename}.py: '
                               f'Filepath {filepath} does not exist.\n')                          
                    logger.warning(message)                            
                                
                    with open('model_definition.json', 'w') as fp:
                        json.dump(self._model_definition, fp, indent=4)    
                        
            else:
                #curdir=os.getcwd()
                abs_path = os.path.dirname(os.path.abspath(__file__))
                abs_path_splitted = abs_path.split('\\')
                working_directory_elements = abs_path_splitted[:-1]
                working_directory = os.path.join(*working_directory_elements)
                splitted = working_directory.split(':') # trick to repair path
                working_directory = splitted[0]+':\\'+splitted[1] # trick to repair path
                mydir = os.path.join(working_directory,'resources')
                # mydir = os.path.join(path_to_parent_folder,'model_definition')
                
                
                if not os.path.isdir(mydir):
                    os.mkdir(mydir)
                    
                # curdir = os.getcwd()
                # filepath = os.path.join(curdir,'resources','model_definition.json')
                filepath = os.path.join(mydir,'model_definition.json')
    
                with open(filepath, 'w') as fp:
                    json.dump(self._model_definition, fp, indent=4)  
                    
            

    
    
    def from_jason(self,filepath=None):
        """ read model settings  from json text files"""

        import json
        if isinstance(filepath,str):
            try:
                # filepath = os.path.join(dirpath,'model_definition.json')
                fp = open(filepath,'r')
                model_definition = json.load(fp)
                fp.close() 
                
            except FileNotFoundError:
                        
                clsname = str(self.__class__.__name__)
                modulename = str(__name__)
    
                message = (f'\nIn class {clsname} of module {modulename}.py: '
                           f'Filepath {filepath} does not exist.\n')                          
                logger.warning(message)   
                model_definition = None                         

        return model_definition
        
        
                                  
if __name__ == "__main__":
 
    md = ModelDefinition()
    print('290 model_definition', md.model_definition)
    
    md.to_json()
    
    curdir = os.getcwd()
    filepath = os.path.join(curdir,'model_definition.json')
    md2 = ModelDefinition().from_jason(filepath=filepath) # case 3: model definition specified in json json file in filepath
    # md2 = ModelDefinition(read_from_resources = True)
    print('303 model_definition', md2)
    