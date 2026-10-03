# -*- coding: utf-8 -*-
"""
Created on Thu Jul  6 20:38:10 2023

@author: christophe
"""

from logging import getLogger
import numpy as np
from utilities import orderOfMagnitude,generate_ticks,datestring2num,file_name_from_path
from matplotlib.pyplot import figure, show,savefig,clabel,contour,setp,gcf,getp,gca,close
from matplotlib.dates import date2num, num2date
from matplotlib.ticker import Formatter,FuncFormatter, NullLocator, FixedLocator
from matplotlib import pyplot as plt
from utilities import stdev_back_from_log, stdev_back_from_log_stochastic

logger = getLogger(__name__)


class Harmonics:


    
    """
    A class that provides instances of seasonal harmonic of time series model
    
        ...
    
    Attributes
    ----------
        
    tseries : numpy array_like
        The time series for which the harmonic component is to be evaluated
        The time numbers in colomn 0 are time numbers as defined in matplotlib.dates 
        function date2num

    tminstr : string (optional)
        The minimum time of the preprocessed time series as string (dd-mm-yyyy HH:MM:SS)

    tmaxstr : string (optional)
        The maximum time of the preprocessed time series as string (dd-mm-yyyy HH:MM:SS)  

    tmin : float
        The minimum time of the preprocessed time series 
        (time number as defined in matplotlib.dates function date2num)

    tmax : float
        The maximum time of the preprocessed time series 
        (time number as defined in matplotlib.dates function date2num)     

    harmonic_component : numpy array_like
        harmonic_component the evaluated harmonic componenet
        
    tseries_type : string
        A string specifying the type of time series ('heads', 'prec', 'evap', 'pump', 'riv')             
        
    amplitude : float
        The amplitude of the harmonic component
        
    phase_shift : float
        The phase_shift of the harmonic component        
        
    popt : numpy array_like
        The optimum parameters vector of the harmonic component  

    residuals : numpy array_like
        The residuals of the fit of the harmonic component            

    pcov : numpy array_like
        The covariance matrix of the fit of the harmonic component  
        
    pcor : numpy array_like
        The correlation matrix of the fit of the harmonic component     
        
    mae : float
        The mean square error of the fit of the harmonic component         
        
    explained_variance : float
        The explained_variance of the fit of the harmonic component         

    
    
    Methods
    -------
    
    harmonic(A,phi,t,freq=1./365.25)
        A method to evalute a harmonic function
        
    harmonic_residuals(p,t,yobs,freq=1./365.25)
        A method to evalute the residuals of a harmonic function
    
    fit_harmonic()
        A method to fit of a harmonic function    
        
    
    """        
    
    
    def __init__(self,tseries, tminstr = None, tmaxstr = None,tmin = None, tmax= None, tseries_type = None):
        
        
        self.tseries = tseries
        self.tminstr = tminstr 
        self.tmaxstr = tmaxstr
        self.tmin = tmin 
        self.tmax = tmax       
        self.harmonic_component = None
        self.amplitude = None
        self.phase_shift = None
        self.popt = None
        self.residuals = None
        self.explained_variance = None
        self.mae = None
        self.pcov = None
        self.pcor = None 
        self.tseries_type = tseries_type
        


    def harmonic(self,A,phi,t,freq=1./365.25):
        
        """ 
        Method to evaluate the harmonic (= sinus) function of time
        given an amplitude A, phase shift phi and frequency f
        
        Parameters
        ----------
        
        A: float
        The amplitude of the harmonic component
        
        phi: float
        The phase shift of the harmonic component        
        
        t: numpy array_like
        The time at which the harmonic component is evaluated
        
        freq: float
            the frequency of the harmonic component
        
        """
    
        return A*np.sin(2*np.pi* freq*(t-phi)) 
    

    def harmonic_residuals(self,p,t,yobs,freq=1./365.25):#func is the function to optimize, p parameters vector, x independant variables, yobs the observations
        
        """ 
        Method to estimate the difference between observed and calculated 
        seasonal harmonic

        
        Parameters
        ----------
        
        p: numpy array_like
        The vector of parameters
        
        t: numpy array_like
        The time at which the harmonic component is evaluated    
        
        yobs: numpy array_like
        The values of the time series for which the harmonic component is evaluated
        
        freq: float
            the frequency of the harmonic component
        
        
        """
        A = np.exp(p[0])
        phi = np.exp(p[1])
        
        y = self.harmonic(A,phi,t, freq = freq)
        residuals = yobs-y
        
        return residuals

        

    def fit_harmonic(self): 

        
        """ 
        Method to fit the harmonic component

        
        Parameters
        ----------
        
        see class attributes
        
        
        """        
        
        tseries_type = self.tseries_type
        tseries = self.tseries
        tminstr = self.tminstr
        tmaxstr = self.tmaxstr
        
        tmin = self.tmin
        tmax = self.tmax      

        from scipy.optimize import least_squares, leastsq

        if tmin is None:
            if tminstr is not None:
                tmin = datestring2num(tminstr)
            else:
                tmin = tseries[0,0]
            
        if tmax is None:    
            if tmaxstr is not None:
                tmax = datestring2num(tmaxstr)
            else:
                tmax = tseries[-1,0]        
            
        mask = (tseries[:,0] >= tmin) & (tseries[:,0] <= tmax)
        tseries = tseries[mask,:]        
        
        yobs = tseries[:,1]
        t = tseries[:,0] 
        mean_level = np.mean(yobs)
        yobs = yobs - mean_level
        
          
        A0 = abs(max(yobs))
        phi0 = 180.
        
        Amin = 1e-9
        Amax = 1e1
        
        phimin = 1e-9
        phimax = 360.
        
        p0 = np.array([np.log(A0),np.log(phi0)])
        lb = np.array([np.log(Amin),np.log(phimin)])
        ub = np.array([np.log(Amax),np.log(phimax)])
          
        bounds = (lb, ub)
 
        # # if using least_squares
        # outcome = least_squares(self.harmonic_residuals, x0 = p0, bounds=bounds,args=(t,yobs), kwargs={},method='trf', ftol=1e-08, xtol=1e-08, gtol=1e-08, x_scale=1.0, loss='linear', f_scale=1.0, diff_step=None, tr_solver=None, tr_options={}, jac_sparsity=None, max_nfev=None, verbose=0)
        # popt = outcome['x']
        # success = outcome['success']
        # cost = outcome['cost']
        # res = outcome['fun']
        # jac = outcome['jac']
        # message = outcome['message']
        # mae = np.mean(abs(res))        
        
        # try:  
        #     A = np.exp(popt[0])
        #     phi = np.exp(popt[1])
        #     JTJ = np.dot(np.transpose(jac),jac)
        #     pcov = np.linalg.inv(JTJ)*np.var(res)
        #     pvar = np.diag(pcov) #parameters variance vector
        #     stdev = np.sqrt(pvar) #parameters standard deviation vector
        #     outer_product = np.outer(stdev,stdev) #outer product of standard error vector
        #     A, stdev_A = stdev_back_from_log_stochastic(popt[0],stdev[0])
        #     phi, stdev_phi = stdev_back_from_log_stochastic(popt[1],stdev[1])
            
        #     harmonic_component = np.empty(np.shape(tseries))
        #     harmonic_component[:,0] = t 
        #     harmonic_component[:,1] = self.harmonic(A,phi,t)
        
        #     observed_deviations_from_mean = tseries[:,1]-np.mean(tseries[:,1])  
        #     SSreg = sum(pow(res,2))
        #     SStot = sum(pow(observed_deviations_from_mean,2))
        #     explained_variance = 100*(1-SSreg/SStot)        
            
            
        #     self.harmonic_component = harmonic_component
        #     self.amplitude = [A,stdev_A]
        #     self.phase_shift = [phi,stdev_phi]
        #     self.popt = popt
        #     self.residuals = res
        #     self.explained_variance = explained_variance
        #     self.mae = mae
        #     self.pcov = pcov
        #     self.pcor = pcor            
            
        #if using leastsq
        try:  
            outcome = leastsq(self.harmonic_residuals,p0,args=(t,yobs),full_output=1,Dfun=None, col_deriv=0, ftol=1.49012e-10, xtol=1.49012e-12, gtol=0.0, maxfev=0, epsfcn=0.0, factor=0.1, diag=None)
            popt = outcome[0]
            res = outcome[2]['fvec']
            pcov = outcome[1]*np.var(res)
            mae = np.mean(abs(res))                
            pvar = np.diag(pcov) #parameters variance vector
            stdev = np.sqrt(pvar) #parameters standard deviation vector
            outer_product = np.outer(stdev,stdev) #outer product of standard error vector              
            pcor = pcov / outer_product  
            A, stdev_A = stdev_back_from_log_stochastic(popt[0],stdev[0])
    
            phi, stdev_phi = stdev_back_from_log_stochastic(popt[1],stdev[1])
            
            phi = phi%365.25
    
            harmonic_component = np.empty(np.shape(tseries))
            harmonic_component[:,0] = t 
            harmonic_component[:,1] = self.harmonic(A,phi,t)
    
    
            observed_deviations_from_mean = tseries[:,1]-np.mean(tseries[:,1])  
            SSreg = sum(pow(res,2))
            SStot = sum(pow(observed_deviations_from_mean,2))
            explained_variance = 100*(1-SSreg/SStot)
                            
            
            self.harmonic_component = harmonic_component
            self.amplitude = [A,stdev_A]
            self.phase_shift = [phi,stdev_phi]
            self.popt = popt
            self.residuals = res
            self.explained_variance = explained_variance
            self.mae = mae
            self.pcov = pcov
            self.pcor = pcor


        except :
            selfname = str(self.__class__.__name__)
            modulename = str(__name__)   

            message = (f'\nIn class {selfname} of module {modulename}.py: '
                        f'Failed fitting an harmonic component to time series.\n'
                        f'outcome: {outcome}\n')                          
            logger.warning(message)       
            
            harmonic_component = np.empty(np.shape(tseries))
            harmonic_component[:,0] = t 
            harmonic_component[:,1] = 0
            
            
            self.harmonic_component = harmonic_component
            self.amplitude = [0,1e9]
            self.phase_shift = [360.,360.]
            self.popt = [1e9,1e9]
            self.residuals = harmonic_component
            self.explained_variance = -1e9
            self.mae = 1.e9
            self.pcov = 1.e9*np.ones((2,2))
            self.pcor = 1.e9*np.ones((2,2))          
            
            
            
            
            

        return self