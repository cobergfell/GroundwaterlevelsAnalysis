# -*- coding: utf-8 -*-
"""
Created on Mon Jul 24 13:17:45 2023

@author: christophe
"""


import scipy
import numpy as np
from preprocessedseries import Preprocessed
from matplotlib.pyplot import figure, show,savefig,clabel,contour,setp,gcf,getp,gca,close
from logging import getLogger


""" 

Description
----------
This module contains functions needed to apply a special type of response functions
that can be used when groundwater fluctuations are measured in the direct vicinity 
of a river. It has been refered to in the literature as the 'flood wave method' because
the river impulse is comparable to a floodwave.

The solutions to the differential equations describing the flow in the aquifer system 
are given for the step responses in Laplace domain. It is applicable for stresses 
such as precipitation or river stage fluctuations.
The modeled aquifer system at the moment is an aquifer system drained by a river,
representaed in cross-section. The aquifer system consists of a highly permeable 
layer underlying a low permeable layer.
The flow in the highly permeable layer is assumed to be only horizontal while the
flow in the poorly permeable layer is assumed to be only vertical.
The time dependent groundwaterhead is calculated in the aquifer and in the
semi-pervious layer as a function of the distance to the river


Parameters
----------

p: model parameter vector where

p[0] = T represnts the transmitivity of the aquifer (dimension: L2.T-1)
p[1] = c represents the resistance to vertical flow of the semi-pervious layer (dimension: T)
p[2] = S represnets the storage coefficient of the aquifer (dimensionless)
p[3] = w represents the specific river bottom resistance to flow (dimension: T.L-1)
p[4] = L represents the aquifer length (dimension: L)

x  : distance between the piezometer and the closest point of the river bench

lp : Laplace parameter


Returns
-------
phi1_laplace, phi2_Laplace : float
head fluctuation in Laplace domain resp. in the semi-pervious layer
and in the underlying aquifer
       

Further reading
--------------
Paper Obergfell, C., M. Bakker, and K. Maas (2016), 
A time-series analysis framework for the flood-wave method to estimate 
groundwater model parameters, Hydrogeology Journal, 1-13. 
https://doi.org/10.1007/s10040-016-1436-5


Thesis: https://doi.org/10.4233/uuid:40454512-e67c-41c5-963b-5862a1b94ac3




"""    

from scipy.special import factorial

def Stehfest_weights(N):
    """ 
    A utility function to calculate Stehfest weights for N points when 
    numerically back transforming the solution in Laplace domain
    to the time domain.
    
    Parameters
    ----------
    
    N: integer
        The number of weights (as far as I understand,
                               the more weights the more accurate 
                               the solution but also the slower 
                               the computation time)
          
    """     


    weights = np.zeros(N)
    for i in range(1,N+1):
        summ = 0
        for k in range(int((i+1)/2),int(min(i,N/2))+1):
            summ += pow(k,N/2.)*factorial(2*k)/(factorial(N/2.-k)*factorial(k)*factorial(k-1)*factorial(i-k)*factorial(2*k-i))

        weights[i-1] = pow(-1,i+N/2.)*summ
    return weights



    

def step_prec_floodwave_model(p,x,lp):

    """ 
    A function to calculate response to a unit step of precipitation
        
    Parameters
    ----------
    
    p: model parameter vector where
    
    p[0] = T represnts the transmitivity of the aquifer (dimension: L2.T-1)
    p[1] = c represents the resistance to vertical flow of the semi-pervious layer (dimension: T)
    p[2] = S represnets the storage coefficient of the aquifer (dimensionless)
    p[3] = w represents the specific river bottom resistance to flow (dimension: T.L-1)
    p[4] = L represents the aquifer length (dimension: L)
    
    x  : distance between the piezometer and the closest point of the river bank
    lp : Laplace parameter        
        
        
    """      
    


    T = p[0]
    c = p[1]
    S = p[2]
    w = p[3]
    L = p[4]

    gamm = np.sqrt(S*lp/(T*(c*S*lp+1)))
    term1 = gamm*(1-np.exp(2*gamm*L))
    term2 = 1.*(1+np.exp(2*gamm*L))/(T*w)
    term3 = 1./(T*w*S*lp**2)
    beta = term3/(term1-term2)
    constant = 1./(S*lp**2)

    phi2_Laplace = beta*(np.exp(gamm*(2*L-x))+np.exp(gamm*x))+constant
    phi1_Laplace = (lp*phi2_Laplace+c)/(lp*(c*S*lp+1))

    return [phi1_Laplace,phi2_Laplace]
    
    

def step_riv_floodwave_model_bc_type_1(p,x,lp):
    
    
    """
    
    A function to calculate the response to a unit step of river stage 
    variation, with L the distance between the river and a closed boundary
    (water divide)
    
    Parameters
    ----------
    
    p: model parameter vector where
    
    p[0] = T represnts the transmitivity of the aquifer (dimension: L2.T-1)
    p[1] = c represents the resistance to vertical flow of the semi-pervious layer (dimension: T)
    p[2] = S represnets the storage coefficient of the aquifer (dimensionless)
    p[3] = w represents the specific river bottom resistance to flow (dimension: T.L-1)
    p[4] = L represents the aquifer length (dimension: L)
    
    x  : distance between the piezometer and the closest point of the river bank
    lp : Laplace parameter    

    
    (see derivation page 97-100 of my notebook)
    
    """    
    T = p[0]
    c = p[1]
    S = p[2]
    w = p[3]
    L = p[4]

    gamm = np.sqrt(S*lp/(T*(c*S*lp+1)))
    term1 = gamm*(-1+np.exp(2*gamm*L))
    term2 = 1.0*(1+np.exp(2*gamm*L))/(T*w)
    term3 = 1.0/(T*w*lp)
    beta = term3/(term1+term2)
    

    phi2_Laplace = beta*(np.exp(gamm*(2*L-x))+np.exp(gamm*x))
    phi1_Laplace = phi2_Laplace/(c*S*lp+1)
    
   
    return [phi1_Laplace,phi2_Laplace]





def step_riv_floodwave_model_bc_type_2(p,x,lp):
    
    """
    An alternative utility function (not used now) to calculate the response 
    to a unit step of river stage variation, with L the distance between the river and 
    a boundary ideally represented as constant head.
    
    Parameters
    ----------
    
    p: model parameter vector where
    
    p[0] = T represnts the transmitivity of the aquifer (dimension: L2.T-1)
    p[1] = c represents the resistance to vertical flow of the semi-pervious layer (dimension: T)
    p[2] = S represnets the storage coefficient of the aquifer (dimensionless)
    p[3] = w represents the specific river bottom resistance to flow (dimension: T.L-1)
    p[4] = L represents the aquifer length (dimension: L)
    
    x  : distance between the piezometer and the closest point of the river bank
    lp : Laplace parameter    

    (see my derivation page 101-102 of my notebook)
    
    """   
    
    T = p[0]
    c = p[1]
    S = p[2]
    w = p[3]
    L = p[4]

    gamma = np.sqrt(S*lp/(T*(c*S*lp+1)))
    denominator = lp*(T*w*gamma*np.sinh(2*gamma*L)+np.cosh(2*gamma*L))
    phi2_Laplace = np.sinh(gamma*(2*L-x))/denominator
    phi1_Laplace = phi2_Laplace/(c*S*lp+1)        

    return [phi1_Laplace,phi2_Laplace]



if __name__ == "__main__":
   
    T = 100.
    c = 110.
    S = 0.1
    w = 0.05
    L = 700.
    p = np.array([T,c,S,w,L])
    x = 50.
    lp = 0.69314718055994529
    
    [phi1_Laplace,phi2_Laplace] = step_prec_floodwave_model(p,x,lp)
    print('176 [phi1_Laplace,phi2_Laplace] as step_prec_floodwave_model(p,x): ',[phi1_Laplace,phi2_Laplace])
    
    [phi1_Laplace,phi2_Laplace] = step_riv_floodwave_model_bc_type_1(p,x,lp)
    print('179 [phi1_Laplace,phi2_Laplace] as step_riv_floodwave_model_bc_type_1(p,x): ',[phi1_Laplace,phi2_Laplace])
    
    [phi1_Laplace,phi2_Laplace] = step_riv_floodwave_model_bc_type_2(p,x,lp)
    print('182 [phi1_Laplace,phi2_Laplace] as step_riv_floodwave_model_bc_type_2(p,x): ',[phi1_Laplace,phi2_Laplace])
    
    