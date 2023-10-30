# -*- coding: utf-8 -*-
"""
Created on Sat Jan 21 13:34:46 2023

@author: christophe

"""
import scipy
import numpy as np
from abc import ABC, abstractmethod
from collections import OrderedDict
from preprocessedseries import Preprocessed
from matplotlib.pyplot import figure, show,savefig,clabel,contour,setp,gcf,getp,gca,close
from scipy import integrate
from logging import getLogger
from floodwave import (Stehfest_weights,
                       step_prec_floodwave_model,
                       step_riv_floodwave_model_bc_type_1,
                       step_riv_floodwave_model_bc_type_2)


logger = getLogger(__name__)


class UnitResponseFunctionsAbstract(ABC):

    """
    An abstract class of response functions - not necessary yet but could serve as blueprint later
    to apply a bit more some of the 'SOLID' coding recomendations such as
    abstract interfaces.

    see  https://towardsdatascience.com/solid-coding-in-python-1281392a6a94
    
    """
    
    @abstractmethod
    def get_initial_parameters():
        """
        See description of the method in UnitResponseFunctionsParent
        
        """         
    pass



    
class UnitResponseFunctionsParent(UnitResponseFunctionsAbstract):
    
    
    """
        A parent class that provides instances of response functions
    
        ...
    
    Attributes
    ----------
        
    t : float
        The relative time ( that is. the time spent after the stress impulse which corresponds to t = 0 )
        at which the response function is evaluated
   
    stresses_dict : python object of data type 'dict'
        _p_stresses_dict is a local copy of the stresses_dict dictionary containing all the entered
        Stress objects used to explain the observed groundwater levels variations.

    Nint: int
        The memory of the stress considered expressed 
        in number of time intervals (a time interval equals a time step)
        
    time_step : float
        The specified time step used for all heads and stresses time series     
          
    arguments_dict : python object of data type 'dict'
        A dictionary specifying additional arguments to pass to the response functions.
        For example: X,Y coordinates of pumping wells           
        
    stress_type : string
        A string specifying the type of stress ('prec', 'evap', 'pump', 'riv')   
        
    pnumber: int
        pnumber is the rank of a parameter in the parameter vector       
        

    
    Methods
    -------
    
    get_initial_parameters()
        Method that returns the initial parameters

    stepfunction()
         Method that returns the step function of a given response function

    blockfunction()
        Method that returns the blokfunction of a given response function
        
    functionmoments()
        Method that returns the function moments (not implemented yet)
    
    """        
    
    
    
    
    
    def __init__(self, stress_type = None, Nint = None,time_step = None, 
                 pnumber = None, arguments_dict = None, t = None):

        self.stress_type = stress_type     
        self.Nint = Nint
        self.time_step = time_step
        self.arguments_dict = arguments_dict
        self.stress_type = stress_type
        self.pnumber = pnumber
        self.t = t

        if (Nint is not None) and (time_step is not None): # overrules argument t
            #t = np.arange(self.time_step,(Nint+1)*self.time_step,self.time_step) # block function time
            t = np.arange(self.time_step,(Nint+2)*self.time_step,self.time_step) # step function time
            self.t = t


    def get_initial_parameters(self, suffix = "", regime = ""):
        

        """
        Method that returns the initial parameters 
        
        Parameters
        ----------
        
        suffix : string (optional)
            An optional string to disambiguate a parameter name (if necessary) 
            
        regime: string (optional)
            Specify the fluctuation regime, when applicable (default is an empty string)


        Returns
        -------
        parameters: python object of data type 'dict'
            A dictionary of initial parameters definitions        
    
        """          
        
        
        
    pass


    def stepfunction(self, p, logtransform):
        
        """
        
        Method that returns the step function

        Parameters
        ----------
        p:  numpy array_like of floats
            Vector of response function parameters.
            
        logtransform:  numpy array-like of booleans
            vector of boolean values indicating if parameters are logtransformed  


        Return
        -------
        output: numpy array_like  
                output is either the step function of its partial derivatives

        """
        pass    
    
    
    def blockfunction(self,p ,logtransform): 

        """ 
        Method that evaluates the time dependant step function 
     
        Parameters
        ----------
        p: array of floats
            the parameter vector
            
        logtransform: array of booleans
            vector of booleans values indicating if a parameter is logtransformed   


        Return
        ----------
        output: numpy array_like  
            output is either the block function of its partial derivatives

         """
        time_step = self.time_step
        step = self.stepfunction(p ,logtransform)
        

        block = np.empty(len(step)-1)
        block[:] = step[0:-1]
        block[1:] = block[1:] - block[:-1]
        block[:] = block[:] / time_step #  to keep block magnitude to one
        
                  
        return block 
 

    def functionmoments(self, p ,logtransform): 
        """
        Method that returns the function moments (not implemented yet)
        
        Parameters
        ----------
        p: array of floats
            the parameter vector
            
        logtransform: array of booleans
            vector of booleans values indicating if a parameter is logtransformed   
        
        
        Return
        ----------
        output: numpy array_like  
            output is either the block function of its partial derivatives
        
        """
        pass    
    
            
     
        
class IncompleteGamma(UnitResponseFunctionsParent):
    
    """
    A subclass of UnitResponseFunctionsParent to generate instances of 
    response functions of the type Incomplete Gamma function with 3 parameters 
    A, a, and n.

    """
    
    def __init__(self, stress_type = None, Nint = None, time_step = None, 
                 pnumber = None, arguments_dict = None, t= None):
        
        UnitResponseFunctionsParent.__init__(self, stress_type = stress_type, Nint = Nint, time_step = time_step, 
                                             pnumber = pnumber, arguments_dict = arguments_dict, t = t)
        
        self.funcname = 'Incomplete_Gamma_Function'
        
        

    def get_initial_parameters(self, suffix = "", regime = ""):
        
        """
        See description of the method in the parent class
        
        """        
        

        stress_type = self.stress_type
        
        if suffix != "":
            suffix = suffix 

        parameters = OrderedDict()

        if stress_type in ["prec","evap"] :
            
            if regime == 'regime_2':
                pname = "A" + suffix
                parameters[pname] = {
                    "isvariable": True,
                    "logtransform": True,
                    "pname": pname,
                    "minvalue": 50.,
                    "maxvalue": 4000.,
                    "initvalue": 1300., 
                    }                    
                
                pname = "a" + suffix
                parameters[pname] = {
                    "isvariable": True,
                    "logtransform": True,
                    "pname": pname,
                    "minvalue": 1e-5,#1e-4,
                    "maxvalue": 1e-2,
                    "initvalue": 0.002
                    }                   
                    
                    
                pname = "n"+suffix
                parameters[pname] = {   
                    "isvariable": True,
                    "logtransform": True,
                    "pname": pname,
                    "minvalue": 0.7,
                    "maxvalue": 3.0,
                    "initvalue": 1.0}           
                

            else:
                pname = "A" + suffix
                parameters[pname] = {
                    "isvariable": True,
                    "logtransform": True,
                    "pname": pname,
                    "minvalue": 100.,
                    "maxvalue": 4000.,
                    "initvalue": 500., 
                    }                    
                
                pname = "a" + suffix
                parameters[pname] = {
                    "isvariable": True,
                    "logtransform": True,
                    "pname": pname,
                    "minvalue": 1e-5,#1e-4,
                    "maxvalue": 1e-2,
                    "initvalue": 0.002
                    
                    }                   
                    
                    
                pname = "n"+suffix
                parameters[pname] = {   
                    "isvariable": True,
                    "logtransform": True,
                    "pname": pname,
                    "minvalue": 0.3,
                    "maxvalue": 4.0,
                    "initvalue": 0.99}           
                
            if stress_type == 'evap':
                pname= "f"+suffix
                parameters[pname] = { 
                    "isvariable": True,
                    "logtransform": True,
                    "pname": pname,
                    "minvalue": 0.5,
                    "maxvalue": 1.5,
                    "initvalue": 1.0, 
                    }  
            
     

        elif stress_type == "pump" :
            pname = "A" + suffix
            parameters[pname] = {
                "isvariable": True,
                "logtransform": True,
                "pname": pname,
                "minvalue": 1.00E-3,
                "maxvalue": 2.,
                "initvalue": 0.05, 
                }                    
            
            pname = "a" + suffix
            parameters[pname] = {
                "isvariable": True,
                "logtransform": True,
                "pname": pname,
                "minvalue": 5e-5,#1e-4,
                "maxvalue": 1e-2,
                "initvalue": 0.0031
                }                   
                
                
            pname = "n"+suffix
            parameters[pname] = {   
                "isvariable": True,
                "logtransform": True,
                "pname": pname,
                "minvalue": 0.7,
                "maxvalue": 3.0,
                "initvalue": 2.3}              
            
        elif stress_type == "riv" :
            pname = "A" + suffix
            parameters[pname] = {
                "isvariable": True,
                "logtransform": True,
                "pname": pname,
                "minvalue": 1.00E-3,
                "maxvalue": 2.,
                "initvalue": 0.85, 
                }                    
            
            pname = "a" + suffix
            parameters[pname] = {
                "isvariable": True,
                "logtransform": True,
                "pname": pname,
                "minvalue": 5e-5,#1e-4,
                "maxvalue": 1e-1,
                "initvalue": 1e-2
                }                   
                
                
            pname = "n"+suffix
            parameters[pname] = {   
                "isvariable": True,
                "logtransform": True,
                "pname": pname,
                "minvalue": 0.3,
                "maxvalue": 3.0,
                "initvalue": 0.5}             
             
        return parameters
   

      
    def integral(self,p,b):
        
        """
         A method that evaluates an integral needed to
        evaluate the partial derivative of the incomplete gamma function
        with regard to parameter n
 
         Parameters
         ----------
             
         p:  numpy array_like
             Vector of parameters of the response function                
         
         t: float
             time since gebin of the impulse

         
         Return
         ----------
         output: float
             Evaluated intergral
  
        """

        n = p[2]
        f = lambda x, n: np.log(x)*pow(x,n-1)*np.exp(-x)


        lb = 0 # lower bound
        ub = b # upper bound
        
        output = scipy.integrate.quad(f, lb, ub, args=(n))[0]

        return output
    
    
    def integrate_derivative_to_n(self,p, ubounds):
        
        """
         A method that performs the numerical intergrations needed
        to evaluate the partial derivative of the incomplete gamma function
        with regard to parameter n
 
         Parameters
         ----------
          p:  array of float
             vector of parameters of the response function
             
         ubounds: numpy array of float values representing the upper bound
         of the numerical integral
         
         
         Return
         ----------
         output: numpy array_like
         Returns the partial derivative with regrads to n. evaluated for each time stamp
            
    
        """
        

        t = self.t
        
        output = np.empty(len(ubounds))
        
        for i in range(0,len(ubounds)):
            output[i] = self.integral(p, p[1]*t[i])  - self.integral(p, np.Inf)*scipy.special.gammainc(p[2],p[1]*t[i])
                 
        
        return output     
                  
                
    def stepfunction(self, p ,logtransform): 
    
        """ 
        
        A method that evaluates the step function 
         
        Parameters
        ----------
        p:  numpy array_like of floats
            Vector of response function parameters.
            
        logtransform:  numpy array-like of booleans
            vector of boolean values indicating if parameters are logtransformed  
        
        The method evaluates step(t) = p[0] *scipy.special.gammainc(p[2],p[1]*t)
        
        or using the commonly used symbols
        
        p[0] = A
        p[1] = a
        p[2] = n
        
        step(t) = A *scipy.special.gammainc(n,a*t)
        
        where:
        
        A is scale factor for the value of the step as time goes to infinity large
        a is a time scale factor that modify the recession curve of the raised groundwater head after recharge
        n is a scale factor for the value of the step as time goes to infinity large
        
        When the response function is applied to model recharge, a fourth parameter p[3] is introduced, 
        
        or using the commonly used symbols
        
        p[3] = f
        
        where f is a scale factor to estimate in a simplified way recharge from gross precipitation
        
        In that case:
        
        step(t) = f * A *scipy.special.gammainc(n,a*t)
             

        
        Return
        -------
        output: numpy array_like  
                output is either the step function of its partial derivatives
        
         """

        t = self.t
        pnumber = self.pnumber
        
        P = np.empty(len(p)) #P capital letter stands for retransformed p (so not log transformed anymore)
        indexes = np.flatnonzero(logtransform)
        P[indexes] = np.exp(p[indexes])
        

        if pnumber == 0: # partial derivative with regards to parameter A
            output = scipy.special.gammainc(P[2],P[1]*t)
                            
            if logtransform[pnumber] == True:
                output  = P[pnumber] * output
 
        elif pnumber == 1: # partial derivative with regards to parameter a
            output  = P[0]*((P[2]*scipy.special.gammainc(P[2],t*P[1])-scipy.special.gammainc(P[2]+1,t*P[1]))/P[1])
            # below: numerical derivative for checking
            # d = 0.01
            # output = P[0]*(scipy.special.gammainc(P[2],t*P[1]*(1+d))-scipy.special.gammainc(P[2],t*P[1]*(1-d)))/(2*P[1]*d)

            if logtransform[pnumber] == True:
                output  = P[pnumber] * output
                
                

        elif pnumber == 2: # partial derivative with regards to parameter n (we do it here numerically)
            # output  = P[0] * self.integrate_derivative_to_n(P, P[1]*t)  # this numerical intergration is too slow and not more accurate
            d = 0.01
            output  = derBR = P[0]*(scipy.special.gammainc(P[2]*(1+d),t*P[1])-scipy.special.gammainc(P[2]*(1-d),t*P[1]))/(2*P[2]*d)
            if logtransform[pnumber] == True:
                output  = P[pnumber] * output            
 
        elif pnumber == 3:
            output  = P[0] * scipy.special.gammainc(P[2],P[1]*t)  
            
            if logtransform[pnumber] == True:
                output  = P[pnumber] * output
            
            
        else: # normal flow (no partial derivative)
            output = P[0] *scipy.special.gammainc(P[2],P[1]*t)    
            
            if len(p) == 4 :
                output  = output * P[3]  
                  
        return output       



class Hantush(UnitResponseFunctionsParent):
    
    """
    A subclass of UnitResponseFunctionsParent to generate instances of 
    response functions of the type Hantush well function

    """
    def __init__(self, stress_type = None, Nint = None, time_step = None, 
                 arguments_dict = None, pnumber = None, t = None):
        
        UnitResponseFunctionsParent.__init__(self, stress_type = stress_type,Nint = Nint, time_step = time_step,
                                             pnumber = pnumber, arguments_dict = arguments_dict, t=t)
        
        self.funcname = 'Hantush_well_function'

        self.Xwell = None
        self.Ywell = None
        self.Xpiezo = None
        self.Ypiezo = None
        self.pnumber = pnumber
        
        
    def get_initial_parameters(self, suffix = "", regime = ""):
        
        """
        See description of the method in the parent class
        
        """         

        stress_type = self.stress_type
        
        if suffix != "":
            suffix = '_' + suffix 

        parameters = OrderedDict()

        pname = "KD" 
        parameters[pname] = {
            "isvariable": True,
            "logtransform": True,
            "pname": pname,
            "minvalue": 10,
            "maxvalue": 4500, 
            "initvalue": 1500.
            
            }     

        pname = "c"
        parameters[pname] = {
            "isvariable": True,
            "logtransform":  True,
            "pname": pname,
            "minvalue": 1e-2,
            "maxvalue": 1e4,
            "initvalue": 500.
            }                   
                    
        pname = "S"
        parameters[pname] = {
            "isvariable": True,
            "logtransform":  True,
            "pname": pname,
            "minvalue": 1.e-9,
            "maxvalue": 4.e-1,
            "initvalue": 0.001
            }                   
                
        
        return parameters
         

                  
                
    def stepfunction(self, p ,logtransform): 
        
              
            
        """ evaluate the step function 
         
             Parameters
             ----------
             p: array of floats
                 the parameter vector
                 
             logtransform: array of booleans
                 vector of booleans values indicating if a parameter is logtransformed   
                 
                 
            arguments_dict: dictionary
                optional dictionary of extra function arguments
   
             
             evaluate the well function of Hantush where 
             
             
             where:
             
             p[0] = KD represents the transmissivity (dimension: L2.T-1)
             p[1] = c represents the leaky top layer resistance to vertical flow (dimension: T)
             p[2] = S represents the stotage coefficient (dimensionless)
             
             

             Return
             ----------
             output: numpy array_like  
                     step function
        
         """

        t = self.t
        arguments_dict = self.arguments_dict
        P = np.empty(len(p)) #P capital letter stands for retransformed p (so not log transformed anymore)
        indexes = np.flatnonzero(logtransform)
        P[indexes] = np.exp(p[indexes])
        pnumber = self.pnumber
        

        # try: 
        Xwell = arguments_dict['Xwell']
        Ywell = arguments_dict['Ywell']
        Xpiezo = arguments_dict['Xpiezo']
        Ypiezo = arguments_dict['Ypiezo']
        
        self.Xwell = Xwell
        self.Ywell = Ywell
        self.Xpiezo = Xpiezo
        self.Ypiezo = Ypiezo
      
        # except:
        #     clsname = str(self.__class__.__name__)
        #     modulename = str(__name__)
        #     message = (f'\nIn class {clsname} line 525 of module {modulename}.py:\n'
        #                f'Xwell and Ywell should be entered in arguments_dict.\n')     
        #     logger.error(message)                
            
            
        r = np.sqrt((Xpiezo - Xwell)**2 +(Ypiezo - Ywell)**2)

        P = np.empty(len(p)) #P capital letter stands for retransformed p (so not log transformed anymore)
        P[:] = p[:]

        indexes = np.flatnonzero(logtransform)
        P[indexes] = np.exp(p[indexes])
        


        from numpy import sign
        from scipy.special import exp1,expi,kn
        KD = P[0]
        c = P[1]    
        S = P[2]

        output = np.zeros(len(t))
        # try:
        rho = r/np.sqrt(KD*c)
        tau = np.log(2./rho*t/(c*S))
        h_inf = kn(0.,rho)
        if h_inf > 0:
            w = (-expi(-rho)-h_inf)/(-expi(-rho)+expi(-rho/2.)) #expi is taken from -inf to x, wheras in exp1 from x to +inf
            I = h_inf+w*expi(-rho/2.*np.exp(abs(tau)))-(w-1.)*expi(-rho*np.cosh(tau))
            F = h_inf+sign(tau)*I
            
            # w = (scipy.special.exp1(rho)-h_inf)/(scipy.special.exp1(rho)-scipy.special.exp1(rho/2))
            # I = h_inf-w*scipy.special.exp1(rho/2*np.exp(abs(tau)))+(w-1)*scipy.special.exp1(rho*np.cosh(tau))
            # F = h_inf+sign(tau)*I  

        # except:
        #     message = (f'\nIn unitresp.py line 667: Failed evaluating Hantush formula\n')
        #     logger.warning(message)
        
        d = 0.05
        #partial derivatives of rho
        # analytical
        lamda = np.sqrt(KD*c)
        DrhoDKD = -r/pow(lamda,2)*1./np.sqrt(c/(4*KD))   # analytical partial derivative of rho with respect to KD appears incorect, check why. 
        DrhoDc = -r/pow(lamda,2)*1./np.sqrt(KD/(4*c))   #partial derivative of tau with respect to c

        # replace analytical DrhoDKD and DrhoDc by their numericl counterpart
        DrhoDKD = (r/np.sqrt(KD*(1+d)*c)-r/np.sqrt(KD*(1-d)*c))/(2*d*KD)
        DrhoDc = (r/np.sqrt(KD*c*(1+d))-r/np.sqrt(KD*c*(1-d)))/(2*d*c)


        #partial derivatives of tau
        # analytical
        DtauDKD = 1./(2*KD)
        DtauDc = -1./(2*c)
        DtauDS = -1./S

           
        # # numerical
        # DtauDKD = (np.log(2*t*np.sqrt(KD*(1+d)*c)/(r*c*S))-np.log(2*t*np.sqrt(KD*(1-d)*c)/(r*c*S)))/(2*d*KD)
        # DtauDc = (np.log(2*t*np.sqrt(KD*c*(1+d))/(r*c*(1+d)*S))-np.log(2*t*np.sqrt(KD*c*(1-d))/(r*c*(1-d)*S)))/(2*d*c)
        # DtauDS = (np.log(2*t*np.sqrt(KD*c)/(r*c*S*(1+d)))-np.log(2*t*np.sqrt(KD*c)/(r*c*S*(1-d))))/(2*d*S)

        
        #partial derivative of w with respect to rho
        term1 = max((scipy.special.kn(1,rho)-np.exp(-rho)/rho)*(scipy.special.exp1(rho)-scipy.special.exp1(rho/2)),1e-50)#to avoid error messages if <1e-200
        term2 = max((np.exp(-rho/2)-np.exp(-rho))/rho*(scipy.special.exp1(rho)-scipy.special.kn(0,rho)),1e-50)
        term3 = max(pow(scipy.special.exp1(rho)-scipy.special.exp1(rho/2),2),1e-50)
        
        DwDrho = (term1-term2)/term3

        #partial derivative of I with respect to rho
        DIDrho = np.empty(len(tau))
        term1taupositive = DwDrho*scipy.special.exp1(rho*np.exp(tau)/2)-w/rho*np.exp(-rho/2*np.exp(tau))
        term1taunegative = DwDrho*scipy.special.exp1(rho*np.exp(-tau)/2)-w/rho*np.exp(-rho/2*np.exp(-tau))
        term2 = DwDrho*scipy.special.exp1(rho*np.cosh(tau))-(w-1)/rho*np.exp(-rho*np.cosh(tau))
        iftaupositive = -scipy.special.kn(1,rho)-term1taupositive+term2
        DIDrho[np.flatnonzero(tau>=0)] = iftaupositive[np.flatnonzero(tau>=0)]
        iftaunegative = -kn(1,rho)-term1taunegative+term2
        DIDrho[np.flatnonzero(tau<0)] = iftaunegative[np.flatnonzero(tau<0)]


        #partial derivative of I with respect to tau
        DIDtau = np.empty(len(tau))
        term1taupositive = -w*np.exp(-rho/2*np.exp(tau))
        term1taunegative = w*np.exp(-rho/2*np.exp(-tau))
        term2 = -(w-1)*np.exp(-rho*np.cosh(tau))*np.sinh(tau)/np.cosh(tau)
        iftaupositive = -term1taupositive+term2
        DIDtau[np.flatnonzero(tau>=0)] = iftaupositive[np.flatnonzero(tau>=0)]
        iftaupositive = -term1taupositive+term2
        DIDtau[np.flatnonzero(tau<0)] = iftaunegative[np.flatnonzero(tau<0)]
        

        #partial derivative of h with respect to rho
        DFDrho = -scipy.special.kn(1,rho)+sign(tau)*DIDrho

        #partial derivative of h with respect to tau
        DFDtau = scipy.special.kn(0,rho)+sign(tau)*DIDtau        

        if pnumber == 0:
            output = 1.0/(4.0*np.pi*pow(KD,2))*F-1.0/(4.0*np.pi*KD)*(DFDrho*DrhoDKD+DFDtau*DtauDKD)   #analytical

            # below is the numerical derivative for comtrol purpose
            # output = -(1.0/(4.0*np.pi*KD*(1+d))-1.0/(4.0*np.pi*KD*(1-d)))/(2*d*KD)*F-1.0/(4.0*np.pi*KD)*(DFDrho*DrhoDKD+DFDtau*DtauDKD) #numerical

            if logtransform[pnumber] == True:
                output  = P[pnumber] * output  
                
            output = -output # change the sign because the negative sifn will be attributed automatically in the module stressconvolution.py
       
        elif pnumber == 1:      
            output = -1.0/(4.0*np.pi*KD)*(DFDrho*DrhoDc+DFDtau*DtauDc)
            if logtransform[pnumber] == True:
                output  = P[pnumber] * output
            output = -output # change the sign because the negative sifn will be attributed automatically in the module stressconvolution.py
   
        elif pnumber == 2: 
            output = -1.0/(4.0*np.pi*KD)*(DFDtau*DtauDS)#analytical
            if logtransform[pnumber] == True:
                output  = P[pnumber] * output
            
            output = -output # change the sign because the negative sifn will be attributed automatically in the module stressconvolution.py
   
        else: # normal flow
            
            output = 1.0/(4.0*np.pi*KD)*F
            
            # mask = np.isnan(output)
            # if len(mask)>0:
            #     output[mask] = 0          
  

        return output       
  

    

    # overwite initial blockfunction because we do not divide it by time step
    def blockfunction(self,p ,logtransform): 

        """ 
        Method that evaluates the time dependant step function 
     
        Parameters
        ----------
        p: array of floats
            the parameter vector
            
        logtransform: array of booleans
            vector of booleans values indicating if a parameter is logtransformed   


        Return
        ----------
        output: numpy array_like  
            output is either the block function of its partial derivatives


         """
        time_step = self.time_step
        step = self.stepfunction(p ,logtransform)
        block = np.empty(len(step)-1)
        block[:] = step[0:-1]
        block[1:] = block[1:] - block[:-1]
               
        return block 
 
    
     
    def funtionmoments(self, p ,logtransform): 
    # To do: evaluate function moments

        pass



    

class FloodWaveModel2L(UnitResponseFunctionsParent):

    
    """ 

    A subclass of UnitResponseFunctionsParent to generate instances of 
    response functions describing in cross section the interaction between a river 
    and the underlying aquifer.

    The aquifer system consists of an aquifer underlying an semi-pervious layer
    The time dependent groundwaterhead is calculated in the aquifer and in the
    semi-pervious layer as a function of the distance to the river. See  the 
    described in cited literature below.

    Further reading
    --------------
    Obergfell, C., M. Bakker, and K. Maas (2016), 
    A time-series analysis framework for the flood-wave method to estimate 
    groundwater model parameters, Hydrogeology Journal, 1-13. 
    https://doi.org/10.1007/s10040-016-1436-5
    
    
    Thesis: https://doi.org/10.4233/uuid:40454512-e67c-41c5-963b-5862a1b94ac3
    
 
    
    """    
    def __init__(self, stress_type = None, Nint = None, time_step = None, 
                 arguments_dict = None, pnumber = None, t = None):

        
        UnitResponseFunctionsParent.__init__(self, stress_type = stress_type, Nint = Nint, time_step = time_step, 
                                             pnumber = pnumber, arguments_dict = arguments_dict, t=t)
        

        self.funcname = 'Floodwave_function_type1'       

        self.Xriv = None
        self.Yriv = None
        self.Xpiezo = None
        self.Ypiezo = None
        self.arguments_dict = arguments_dict
        self.pnumber = pnumber

        
        
    def get_initial_parameters(self, suffix = "", regime = ""):
        
        """
        See description of the method in the parent class
        
        """ 
        stress_type = self.stress_type
        
        if suffix != "":
            suffix = '_' + suffix 

        parameters = OrderedDict()
  
        

        pname = "T" 
        parameters[pname] = {
            "isvariable": True,
            "logtransform": True,
            "pname": pname,
            "minvalue": 10,
            "maxvalue": 4500, 
            "initvalue": 1.00e2
            }     

        pname = "c"
        parameters[pname] = {
            "isvariable": True,
            "logtransform":  True,
            "pname": pname,
            "minvalue": 1e1,
            "maxvalue": 2e3,
            "initvalue": 1.1e2
            }                   
                    
        pname = "S"
        parameters[pname] = {
            "isvariable": True,
            "logtransform":  True,
            "pname": pname,
            "minvalue": 1.e-2,
            "maxvalue": 2.e-1,
            "initvalue": 1.0e-1
            }   

        pname = "w"
        parameters[pname] = {
            "isvariable": True,
            "logtransform":  True,
            "pname": pname,
            "minvalue": 1e-2,
            "maxvalue": 1e-1,
            "initvalue": 5.0e-2
            }                   
                    
        pname = "L"
        parameters[pname] = {
            "isvariable": True,
            "logtransform":  True,
            "pname": pname,
            "minvalue": 3e2,
            "maxvalue": 3e3,
            "initvalue": 7.00e2
            }                  
                
        return parameters
         

                  
                
    def stepfunction(self, p ,logtransform): 
            
        """ 
        Method that evaluates the step function in the floodwave method
        
        Parameters
        ----------
        T = p[0]: transmitivity of the aquifer (L2.T-1)
        c = p[1]: resistance to vertical flow of the semi-pervios layer (T)
        S = p[2]: storage coefficient of the aquifer (-)
        w = p[3]: specific river bottom resistance to flow (T.L-1)
        L = p[4]: aquifer length (L)

         Return
         ----------
         output: numpy array_like  
                 step function
        
         """
         
        t = self.t
        pnumber = self.pnumber
        P = np.empty(len(p)) #P capital letter stands for back-transformed p (so not log transformed anymore)
        indexes = np.flatnonzero(logtransform)
        P[indexes] = np.exp(p[indexes])
  
    
        arguments_dict = self.arguments_dict
    
            
        try:
            self.Xriv = arguments_dict['Xriv']
            self.Yriv = arguments_dict['Yriv']
            self.Xpiezo = arguments_dict['Xpiezo']
            self.Ypiezo = arguments_dict['Ypiezo']
            self.Zpiezo = arguments_dict['Zpiezo'] # depth piezometer in m above datum
            self.Lpiezo = arguments_dict['Lpiezo'] # groundwater model layer number of piezometer   
            
        except:
            clsname = str(self.__class__.__name__)
            modulename = str(__name__)
            message = (f'\nIn class {clsname} line 809 of module {modulename}.py:\n'
                        f'Xriv and Yriv should be entered in arguments_dict.\n')     
            logger.error(message)    
 
        
        stress_type = self.stress_type
        
        try:
            x = np.sqrt((self.Xpiezo - self.Xriv)**2 +(self.Ypiezo - self.Yriv)**2)

        except:
            clsname = str(self.__class__.__name__)
            modulename = str(__name__)
            message = (f'\nIn class {clsname} line 922 of module {modulename}.py:'
                       f'Enter Xriv and Yriv in arguments_dict.\n')     
                        
            logger.warning(message)   
            
 
        P = np.empty(len(p)) #P capital letter stands for retransformed p (so not log transformed anymore)
        P[:] = p[:]

        indexes = np.flatnonzero(logtransform)
        P[indexes] = np.exp(p[indexes])
        
        T = P[0]
        c = P[1]
        S = P[2]
        w = P[3]
        L = P[4]
        

        output = np.zeros(len(t))

        # try:

        N = 10
        weight = Stehfest_weights(N)

        d = 0.01

        indices = [0,1,2,3,4]

        if pnumber in indices:
            Pplus = np.array(P)
            Pplus[pnumber] *= (1+d)
            Pmin = np.array(P)
            Pmin[pnumber] /= (1+d)
            
        
            sr1plus = np.zeros(len(t))
            sr1min = np.zeros(len(t))     
            sr2plus = np.zeros(len(t))
            sr2min = np.zeros(len(t)) 
            
            
            for j in range(0,len(t)):
                sum1plus = 0
                sum1min = 0
                sum2plus = 0
                sum2min = 0
                
                for i in range (1,N+1):
    
                    lp = i*np.log(2)/t[j] # lp is the Laplace variable
                    
                    if stress_type in ['prec','evap']:
                        sr1plus_Laplace,sr2plus_Laplace = step_prec_floodwave_model(Pplus,x,lp)
                        sr1min_Laplace,sr2min_Laplace = step_prec_floodwave_model(Pmin,x,lp)
                        
                    elif stress_type =='riv':
                        sr1plus_Laplace,sr2plus_Laplace = step_riv_floodwave_model_bc_type_1(Pplus,x,lp)
                        sr1min_Laplace,sr2min_Laplace = step_riv_floodwave_model_bc_type_1(Pmin,x,lp)                    
                        
         
                    sum1plus += weight[i-1]*sr1plus_Laplace*lp/i
                    sum2plus += weight[i-1]*sr2plus_Laplace*lp/i
                    sum1min += weight[i-1]*sr1min_Laplace*lp/i
                    sum2min += weight[i-1]*sr2min_Laplace*lp/i                
                
                sr1plus[j] = sum1plus
                sr1min[j] = sum1min
                sr2plus[j] = sum2plus
                sr2min[j] = sum2min        
            
            
            derivative_sr1 = (sr1plus - sr1min)/(2*d*P[pnumber])
            derivative_sr2 = (sr2plus - sr2min)/(2*d*P[pnumber])         
    
                 
        
            if logtransform[pnumber] == True:
                derivative_sr1 *= P[pnumber]
                derivative_sr2 *= P[pnumber]
            
            if self.Lpiezo == 0:
                output = derivative_sr1
            else:
                output = derivative_sr2
                
                

        else: # normal flow   
        
            sr1 = np.zeros(len(t))
            sr2 = np.zeros(len(t))
    
            for j in range(0,len(t)):
                sum1 = 0
                sum2 = 0
                
                for i in range (1,N+1):
    
                    lp = i*np.log(2)/t[j] # lp is the Laplace variable
                    
                    if stress_type in ['prec','evap']:
                        sr1_Laplace,sr2_Laplace = step_prec_floodwave_model(P,x,lp)
    
                    elif stress_type =='riv':
    
                        sr1_Laplace,sr2_Laplace = step_riv_floodwave_model_bc_type_1(P,x,lp)
         
                    sum1 += weight[i-1]*sr1_Laplace*lp/i
                    sum2 += weight[i-1]*sr2_Laplace*lp/i
                
                sr1[j] = sum1
                sr2[j] = sum2
            

            if self.Lpiezo == 0:
                output = sr1
            else:
                output = sr2
            

        # except:
        #     clsname = str(self.__class__.__name__)
        #     modulename = str(__name__)
        
        #     message = (f'\nIn class {clsname} of module {modulename}.py:'
        #                 f'\nFailed evaluating floodwave model.\n')
        #     logger.warning(message)   
            
              
        
        return output       
  

    
    # overwite initial blockfunction because we do not divide it by time step
    def blockfunction(self,p ,logtransform): 


        """ 
        Method that evaluates the time dependant step function 
     
        Parameters
        ----------
        p: array of floats
            the parameter vector
            
        logtransform: array of booleans
            vector of booleans values indicating if a parameter is logtransformed   


        Return
        ----------
        output: numpy array_like  
            output is either the block function of its partial derivatives

         """
         
        time_step = self.time_step
        step = self.stepfunction(p ,logtransform)

        block = np.empty(len(step)-1)
        block[:] = step[0:-1]
        block[1:] = block[1:] - block[:-1]
             
        return block 


   

                                  
if __name__ == "__main__":
    #from matplotlib.pyplot import figure, show,savefig,clabel,contour,setp,gcf,getp,gca,close
    arguments_dict = {}
    arguments_dict['Xwell'] = 0
    arguments_dict['Ywell'] = 0
    arguments_dict['Xpiezo'] = 211.
    arguments_dict['Ypiezo'] = 0

    time_step = 1./(24.*2)
    
    p = np.array([1500,500,0.001])
    logtransform = np.zeros(len(p),dtype=bool)

    hantush = Hantush(arguments_dict = arguments_dict, Nint = 480, time_step=time_step)
    t = hantush.t
    step = hantush.stepfunction( p ,logtransform)
    block = hantush.blockfunction( p ,logtransform)
    
    print('1272 step',step[:5])
    print('1273 block',block[:5])
    input()
        
    fig=figure()
    ax = fig.add_subplot(111)
    ax.set_title('Step function hantush plotted for control')
    ax.set_xlabel('days')
    ax.set_ylabel('m')
    ax.plot(t,step,'b')
    ax.set_xlim(t[0],t[-1])
    show()
    
    
    hantush = Hantush(arguments_dict = arguments_dict, Nint = 480, time_step=time_step, pnumber = 0)
    step = hantush.stepfunction( p ,logtransform)
    block = hantush.blockfunction( p ,logtransform)
    
    print('1280 step pnumber = 0',step[:5])
    print('1281 block pnumber = 0',block[:5])
    input()
    
    
    hantush = Hantush(arguments_dict = arguments_dict, Nint = 480, time_step=time_step, pnumber = 1)
    step = hantush.stepfunction( p ,logtransform)
    block = hantush.blockfunction( p ,logtransform)
    
    print('1299 step pnumber = 1',step[:5])
    print('1300 block pnumber = 1',block[:5])
    input()    
    
    hantush = Hantush(arguments_dict = arguments_dict, Nint = 480, time_step=time_step, pnumber = 2)
    step = hantush.stepfunction( p ,logtransform)
    block = hantush.blockfunction( p ,logtransform)
    
    print('1307 step pnumber = 2',step[:5])
    print('1308 block pnumber = 2',block[:5])
    input()     
    
    
    
    
    