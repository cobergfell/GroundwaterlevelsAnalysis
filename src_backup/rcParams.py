# -*- coding: utf-8 -*-
"""
Created on Wed Mar 22 14:28:30 2023

@author: christophe
"""
from matplotlib import rcParams

params = {'legend.fontsize': 12}
rcParams = {'legend.fontsize': 12}
rcParams['font.size'] = 12#This is to change default font size
rcParams['text.usetex'] = False#True
rcParams.update(params)    
