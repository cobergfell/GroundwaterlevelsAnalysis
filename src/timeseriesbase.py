# -*- coding: utf-8 -*-
"""
Created on Fri Jan 27 12:29:58 2023

@author: christophe
"""
from logging import getLogger
import numpy as np
import os
from abc import ABC, abstractmethod
from matplotlib.pyplot import figure, show,savefig,clabel,contour,setp,gcf,getp,gca,close
from matplotlib.dates import date2num, num2date
from matplotlib.ticker import Formatter,FuncFormatter, NullLocator, FixedLocator
from matplotlib import pyplot as plt
from harmonics import Harmonics
from utilities import orderOfMagnitude,generate_ticks,datestring2num,file_name_from_path


logger = getLogger(__name__)

class TimeSeriesBase:
    
    """
    A class provides the basic blueprint shared by all time series

    Attributes
    ----------

    _observed: numpy array_like (optional)
        A time series can be passed directly as classs initialisation input
        in the form of a 2 dimensional numpy array of float, with time numbers
        given in the first array axis (time numbers as defined in matplotlib.dates 
        function date2num), and observed values given in the second array axis
        
    _metadata : python object of data type 'dict' (optional)
        An dictionary of metadata which can be entered together with the
        time series array        
        
    name: string (optional)
        String with the name of the time series, if None, a default name is provided
        
    tseries_type : string
        A string specifying the type of time series ('heads', 'prec', 'evap', 'pump', 'riv')           
        
    tminstr : string (optional)
        The minimum time of the preprocessed time series as string (dd-mm-yyyy HH:MM:SS)

    tmaxstr : string (optional)
        The maximum time of the preprocessed time series as string (dd-mm-yyyy HH:MM:SS)  

    settings : python object of data type 'dict' (optional)
        An optional dictionary of settings - not used in present version
        
    cumulative : boolean
        A boolean variable specifying if the time series corresponds to a quantity 
        that can be cumulated such as precipitation or evaporation (unlike heads or river stages)     

    _interpolated: numpy array_like (optional)
        The interpolated version of the time series      
        
    _interpolated_normalized: numpy array_like (optional)
        The normalized version of the interpolated time series        
        
    _modeled: numpy array_like (optional)
        The modeled version of the time series 
        
   _harmonic_observed: numpy array_like (optional)
        The seasonal harmonic component present in the observedtime series   
    
   _harmonic_interpolated: numpy array_like (optional)
        The seasonal harmonic component present in the interpolated time series  
        
   _harmonic_modeled: numpy array_like (optional)
        The seasonal harmonic component present in the modeled time series      
        
        

    Methods (other than getters and setters)
    -------
    plot()
        A method to plot the time series

    plot_harmonics()
        A method to plot the seasonal harmonic of the time series
        
    """


    def __init__(self, observed = None, name = None, settings = None, metadata = None, tminstr = None, 
                 tmaxstr = None, cumulative = None, tseries_type = None):


            self._observed = observed 
            self._metadata = metadata
            self.name = name 
            self.tseries_type = tseries_type
            self.tminstr = tminstr 
            self.tmaxstr = tmaxstr 
            self.settings = settings    
            self.cumulative = cumulative    
            self._interpolated = None
            self._interpolated_normalized = None
            self._modeled = None
            self._modeled_normalized = None
            self._harmonic_observed = Harmonics(observed).fit_harmonic()
            self._harmonic_interpolated = None
            self._harmonic_modeled = None


    def __repr__(self):
        """Prints a representation of the class instance."""
        return f"{self.__class__.__name__}" \
               f"(name={self.name}, " \
               f"tmin={self.settings['tmin']}, " \
               f"tmax={self.settings['tmax']})"
            

    @property
    def observed(self):
        """getter for updated observed time series."""
        return self._observed

    @observed.setter
    def observed(self, observed):
        """Setter for updated observed time series."""
        self._observed = observed
        
    @property
    def metadata(self):
        """getter for updated metadata of observed time series."""
        return self._metadata

    @metadata.setter
    def metadata(self, metadata):
        """Setter for updated metadata of observed time series."""
        self._metadata = metadata

    @property
    def interpolated(self):
        """getter for updated interpolated time series."""
        return self._interpolated

    @interpolated.setter
    def interpolated(self,interpolated):
        """setter for updated interpolated time series."""
        self._interpolated = interpolated       
        
        
    @property
    def interpolated_normalized(self):
        """getter for updated interpolated and normalized time series."""
        return self._interpolated_normalized

    @interpolated_normalized.setter
    def interpolated_normalized(self,interpolated_normalized):
        """setter for updated interpolated and normalized time series."""
        self._interpolated_normalized = interpolated_normalized            
        
        
               
    @property
    def modeled(self):
        """getter for updated modeled time series."""
        return self._modeled

    @modeled.setter
    def modeled(self,modeled):
        """setter for updated modeled time series."""
        self._modeled = modeled
        
        
    @property
    def modeled_normalized(self):
        """getter for updated modeled and normalized time series."""
        return self._modeled_normalized

    @modeled_normalized.setter
    def modeled_normalized(self,modeled_normalized):
        """setter for updated modeled and normalized time series."""
        self._modeled_normalized = modeled_normalized      
   
     
    @property
    def harmonic_observed(self):
        """getter for seasonal harmonic in observed series."""
        return self._harmonic_observed
     
    
    @harmonic_observed.setter
    def harmonic_observed(self,harmonic_observed):
        """setter for updated observed time series harmonic."""
        self._harmonic_observed = harmonic_observed    


    @property
    def harmonic_interpolated(self):
        """getter for interpolated seasonal harmonic."""
        return self._harmonic_interpolated
    
    
    @harmonic_interpolated.setter
    def harmonic_interpolated(self,harmonic_interpolated):
        """setter for updated interpolated time series harmonic."""
        self._harmonic_interpolated = harmonic_interpolated        
    

    @property
    def harmonic_modeled(self):
        """getter for seasonal harmonic of modeled sime series."""
        return self._harmonic_modeled
    
    
    @harmonic_modeled.setter
    def harmonic_modeled(self,harmonic_modeled):
        """setter for updated harmonic of modeled time series."""
        self._harmonic_modeled = harmonic_modeled        
        

         
    def plot(self, tseries_list = None, legend_list = None, share_axes = True, plot_title = None, yas_title = None, 
             tminstr = None, tmaxstr = None, save_plot = False):
        
        """ 
        A method to plot the time series.
        
        Parameters
        ----------    
        tseries_list: list (optional)
                list of time series to plot. The time series are given in the 
                form of a 2 dimensional numpy array of float, with time numbers
                given in the first array axis (time numbers as defined in 
                matplotlib.dates function date2num), and observed values 
                given in the second array axis
                
        tseries_list_names: list (optional)
                list of time series names given as strings

        
        yas_title: string (optional)
                string specifying Y-axis title   
                
        plot_title: string (optional)
                string specifying plot title 
                
                
        tminstr: string (optional)
                string specifying the minimum time on the plot (format dd-mm-yyyy)
                
        tmaxstr: string (optional)
                string specifying the maximum time on the plot (format dd-mm-yyyy)     
                
        save_plot: boolean (optional)
            Boolean variable specifying if the plot is to be saved on disk
            
    
        Returns
        -------
        ax: Matplotlib Axes object
                The Matplotlib Axes object

        """
        
        clsname = str(self.__class__.__name__)
        modulename = str(__name__)   
        
        observed = self._observed
        interpolated = self._interpolated
        tseries_type = self.tseries_type

        left_margin=0.12
        right_margin=0.06
        top_margin=0.12
        bottom_margin=0.14
        dy=0.14
        dx=0.05        

        fs = 12
        colors = ['b-','m-','y-','c-','g-', 'tab:green-','tab:olive-', 'tab:purple-','tab:pink-','tab:orange-','tab:brown-']
            
        
        if tminstr is not None:
            tmin = datestring2num(tminstr)
        else:
            tmin = observed[0,0]
                        
        if tmaxstr is not None:
            tmax = datestring2num(tmaxstr)
        else:
            tmax = observed[-1,0]

        mask = (observed[:,0] >= tmin) & (observed[:,0] < tmax)
        observed = observed[mask]        

        if tseries_list is None:
            tseries_list = [observed]
            legend_item = self.name+' '+'observed'
            legend_list = [legend_item]
            
            
        if legend_list is None:
            default_list = []
            for i in range(0,len(tseries_list )):
                default_list.append('series_'+str(i))
            
            legend_list = default_list


        nl = len(tseries_list)
        if nl > 4:
            message = (f'\nIn class {clsname} of module {modulename}.py: '
                       f'The list of time series to plot should not exceed 4.\n')                          
            logger.warning(message) 
            
    
        if plot_title is None:  
            # title = (f'time series type: {(self.tseries_type} name: {self.name}') 
            plot_title = f'{self.name}'

              
        if yas_title is None:  
            yas_title = 'm above datum'


        if share_axes == False:
            nc = 1.0
            plot_length = (1-left_margin-right_margin-dx)/nc
            plot_height = (1-bottom_margin-top_margin-2*dy)/nl
            x1 = left_margin
            x2 = x1+plot_length+dx
            X = [x1]#list of x coordinates of bottom left corners
            y1 = 1-top_margin-plot_height
            
            Y=[y1]#list of x coordinates of bottom left corners  
            if nl>1:
                for i in range(1,nl):
                    Y.append(Y[i-1]-plot_height-dy)
            
            fig=figure(num=None, figsize=(8,8*nl),dpi=50,facecolor='w', edgecolor='k')# figsize=(width, height) in inches.

            for i in range(0,nl): 
                ax = fig.add_axes([X[0],Y[i],plot_length,plot_height],frameon=True, xscale=None, yscale=None)  
                ax.plot_date(tseries_list[i][:,0],tseries_list[i][:,1],colors[i],linewidth = 1.0)     
           
            ax.set_title(plot_title,fontsize=fs)
            (x, y) = ax.title.get_position()
            #ax.title.set_y(0.95 * y)    
            ax.grid(True)         
          
                
            def format_date(dates, pos=None):
                return num2date(dates).strftime('%Y')
            
            # def format_date(dates, pos=None):
            #     return num2date(dates).strftime('%d-%m-%Y %H:%M')        
        
            ax.set_ylabel(yas_title,fontsize=fs)
          
            ax.set_xlabel('date (yyyy)',fontsize=fs)
            ax.xaxis.set_major_formatter(FuncFormatter(format_date))
    
            
            for tick in ax.xaxis.get_major_ticks():
                tick.label.set_fontsize(fs)
                tick.label.set_rotation(0)
                # tick.label.set_rotation(45)
            
            for label in ax.xaxis.get_majorticklabels():
                #ha=label.get_horizontalalignment()
                label.set_horizontalalignment('center')
                
                
        
            ax.annotate('Legend', xy=(0.05, -0.14), xycoords='axes fraction',horizontalalignment='left', verticalalignment='center',fontsize=fs)#this is for 1 column 2 rows
            leg = ax.legend((legend_list), ncol=2,shadow=False,bbox_to_anchor=[0.045, -0.27], loc='lower left',frameon=False)
            ltext  = leg.get_texts() 
            plt.setp(ltext, fontsize=fs)    # the legend text fontsize
            
            show()
            
        else:
            plot_length = (1-left_margin-right_margin-dx)
            plot_height = (1-bottom_margin-top_margin-2*dy)
            x_lower_left = left_margin
            y_lower_left = 1-top_margin-plot_height

            fig=figure(num=None, figsize=(8,8),dpi=50,facecolor='w', edgecolor='k')# figsize=(width, height) in inches.
            ax = fig.add_axes([x_lower_left,y_lower_left,plot_length,plot_height],frameon=True, xscale=None, yscale=None)  

            for i in range(0,nl): 
                ax.plot_date(tseries_list[i][:,0],tseries_list[i][:,1],colors[i],linewidth = 1.0)     
    
            ax.set_title(plot_title,fontsize=fs)
            (x, y) = ax.title.get_position()
            #ax.title.set_y(0.95 * y)    
            ax.grid(True)         
        
        
    
                
            def format_date(dates, pos=None):
                return num2date(dates).strftime('%Y')
            
            def format_date(dates, pos=None):
                return num2date(dates).strftime('%d-%m-%Y %H:%M')    
            
            def format_date(dates, pos=None):
                return num2date(dates).strftime('%d-%m-%Y')             
            
        
            ax.set_ylabel(yas_title,fontsize=fs)
          
            ax.set_xlabel('date (yyyy)',fontsize=fs)
            ax.xaxis.set_major_formatter(FuncFormatter(format_date))
    
            
            for tick in ax.xaxis.get_major_ticks():
                tick.label.set_fontsize(fs)
                tick.label.set_rotation(0)
                # tick.label.set_rotation(45)
            
            for label in ax.xaxis.get_majorticklabels():
                #ha=label.get_horizontalalignment()
                label.set_horizontalalignment('center')
                
                
        
            ax.annotate('Legend', xy=(0.05, -0.14), xycoords='axes fraction',horizontalalignment='left', verticalalignment='center',fontsize=fs)#this is for 1 column 2 rows
            leg = ax.legend((legend_list), ncol=nl,shadow=False,bbox_to_anchor=[0.045, -0.27], loc='lower left',frameon=False)
            ltext  = leg.get_texts() 
            plt.setp(ltext, fontsize=fs)    # the legend text fontsize
            
            show()            
            
            
            
            
            
            
            
            
        
        if save_plot == True:
            if self.name is not None:
                figname = "plot_of_{}".format(self.name)
            else:
                figname = "time_series_plot"
            try: 
                abs_path=os.path.dirname(os.path.abspath(__file__))
                abs_path_splitted=abs_path.split('\\')
                path_to_parent_folder_elements=abs_path_splitted[:-1]
                path_to_parent_folder=os.path.join(*path_to_parent_folder_elements)
                splitted=path_to_parent_folder.split(':')#trick  to  repair  path
                path_to_parent_folder=splitted[0]+':\\'+splitted[1]#trick  to  repair  path
                path_to_file_folder=os.path.join(path_to_parent_folder,'results')
                savefig(path_to_file_folder+'\\'+figname, dpi=None, facecolor='w', edgecolor='w',orientation='portrait', format=None,transparent=False, bbox_inches=None, pad_inches=0.1)
              
            except:
                curdir=os.getcwd()
                savefig(curdir+'\\'+figname, dpi=None, facecolor='w', edgecolor='w',orientation='portrait', format=None,transparent=False, bbox_inches=None, pad_inches=0.1)
                       
        

            
        return ax   



    def plot_harmonics(self, yas_title = None,plot_title = None, tminstr = None, tmaxstr = None,
                       plot_observed = True, plot_interpolated = True, plot_modeled = True):
        
        """ 
        A method to plot the seasonal harmonic component of a time series.
        
        Parameters
        ----------    

        yas_title: string (optional)
                string specifying Y-axis title   
                           
        plot_title: string (optional)
                string specifying plot title 
                                            
        plot_observed: boolean (optional)
            Boolean variable specifying if the seasonal harmonic in the observed heads
            is to be ploted
                
        plot_interpolated: boolean (optional)
            Boolean variable specifying if the seasonal harmonic in the interpolated heads
            is to be ploted                

        plot_modeled: boolean (optional)
            Boolean variable specifying if the seasonal harmonic in the modeled heads
            is to be ploted    
            
     
        tminstr: string (optional)
                string specifying the minimum time on the plot (format dd-mm-yyyy)
                
        tmaxstr: string (optional)
                string specifying the maximum time on the plot (format dd-mm-yyyy)     
                
           
    
        Returns
        -------
        ax: Matplotlib Axes object
                The Matplotlib Axes object

        """

        observed = self._observed
        interpolated = self._interpolated
        modeled = self._modeled
        name = self.name
        tseries_type = self.tseries_type

        tseries_list = []
        harmonics_list = []
        legend_list_of_list = []


        def tseries_cut_off(tseries,tminstr=None,tmaxstr=None ):
            """
            A function to select the part of a time series 
            between a minimum time tminstr and maximum time tmaxstr
            
            """
            if tminstr is not None:
                tmin = datestring2num(tminstr)
            else: 
                tmin = tseries[0,0] 
                    
            if tmaxstr is not None:
                tmax = datestring2num(tmaxstr)
            else:
                tmax = tseries[-1,0]
                
            mask = (tseries[:,0] >= tmin) & (tseries[:,0] < tmax)
            
            return tseries[mask]

  
        if self._observed is not None:
            if plot_observed is True:
                normalized_series = np.empty(np.shape(self._observed))
                normalized_series[:,0] = self._observed[:,0]
                normalized_series[:,1] = self._observed[:,1] - np.mean(self._observed[:,1])
                normalized_series = tseries_cut_off(normalized_series,tminstr=tminstr,tmaxstr=tmaxstr)
                tseries_list.append(normalized_series)
                harmonics_list.append(self.harmonic_observed.harmonic_component)
                legend_list_of_list.append(['Observed fluctuations','Harmonic of observed series'])
                
                
                
        if self._interpolated is not None:
            if plot_interpolated is True:
                normalized_series = np.empty(np.shape(self._interpolated))
                normalized_series[:,0] = self._interpolated[:,0]
                normalized_series[:,1] = self._interpolated[:,1] - np.mean(self._interpolated[:,1])
                normalized_series = tseries_cut_off(normalized_series,tminstr=tminstr,tmaxstr=tmaxstr)
                tseries_list.append(normalized_series)
                harmonics_list.append(self.harmonic_interpolated.harmonic_component)
                legend_list_of_list.append(['Interpolated fluctuations','Harmonic of interpolated series'])                
                
          
        if self._modeled is not None:
            if plot_modeled is True:
                normalized_series = np.empty(np.shape(self._modeled))
                normalized_series[:,0] = self._modeled[:,0]
                normalized_series[:,1] = self._modeled[:,1] - np.mean(self._modeled[:,1])
                normalized_series = tseries_cut_off(normalized_series,tminstr=tminstr,tmaxstr=tmaxstr)
                tseries_list.append(normalized_series)
                harmonics_list.append(self.harmonic_modeled.harmonic_component)
                legend_list_of_list.append(['Modeled fluctuations','Harmonic of modeled series'])    

                

        # try:
        left_margin = 0.1
        right_margin = 0.06
        top_margin = 0.06
        bottom_margin = 0.14
        dy = 0.14
        dx = 0.05        
        
        nl = len(tseries_list)
        nc = 1
        plot_length = (1-left_margin-right_margin-dx)/nc
        plot_height = (1-bottom_margin-top_margin-2*dy)/nl
        x1 = left_margin
        x2 = x1+plot_length+dx
        #X=[x1,x2,x1,x2,x1,x2,x1,x2]#list of x coordinates of bottom left corners
        X = [x1,x1,x1]#list of x coordinates of bottom left corners
        y1 = 1-top_margin-plot_height
        y2 = y1-plot_height-dy
        y3 = y2-plot_height-dy
        y4 = y3-plot_height-dy
        Y = [y1,y2,y3]             
    
        fs = 12
    
        fig = figure(num=None, figsize=(12,nl*8),dpi=50,facecolor='w', edgecolor='k')#figures with all plots in a colomn
    
        axnum = -1
        for i in range(0,nl):

            # axnum = axnum+1
            # ax = 'ax'+str(axnum)
            
            tseries = tseries_list[i]
            harmonics = harmonics_list[i]
            leg_list = legend_list_of_list[i]
            ax = fig.add_axes([X[i],Y[i],plot_length,plot_height],frameon=True, xscale=None, yscale=None)
            ax.plot_date(tseries[:,0],tseries[:,1],'b-',linewidth = 1.0)
            ax.plot_date(harmonics[:,0],harmonics[:,1],'r-',linewidth = 1.0)
    
            if plot_title is None:
                title = (f'Harmonic component of {tseries_type} time series {name}')
            else:
                title = plot_title
    
            ax.set_title(title,fontsize=fs)
                
            (x, y) = ax.title.get_position()
            #ax.title.set_y(0.95 * y)    
            ax.grid(True)            

    
            def format_date(dates, pos=None):
                return num2date(dates).strftime('%Y')
            
            # def format_date(dates, pos=None):
            #     return num2date(dates).strftime('%d-%m-%Y %H:%M')        
        
        
            # ticks = generate_ticks(vmin,vmax)
            # ax.set_ylim(ticks[0],ticks[-1]
            # majorLocator   = FixedLocator(ticks) 
            # ax.yaxis.set_major_locator(majorLocator)
            ax.set_ylabel(yas_title,fontsize=fs)
          
            # ax.set_xlabel('date (dd-mm-yyyy)',fontsize=fs, rotation=20)
            ax.set_xlabel('date (yyyy)',fontsize=fs)
            ax.xaxis.set_major_formatter(FuncFormatter(format_date))
            # ax.set_ylabel('values')
            # ax.fmt_xdata = DateFormatter('%Y-%m-%d')   
            # fig.autofmt_xdate()
            #ax.xaxis.set_label_position('right')
            
            for tick in ax.xaxis.get_major_ticks():
                tick.label.set_fontsize(fs)
                tick.label.set_rotation(0)
                # tick.label.set_rotation(45)
            
            for label in ax.xaxis.get_majorticklabels():
                #ha=label.get_horizontalalignment()
                label.set_horizontalalignment('center')
                
                
        
            ax.annotate('Legend', xy=(0.05, -0.12), xycoords='axes fraction',horizontalalignment='left', verticalalignment='center',fontsize=fs)#this is for 1 column 2 rows
            leg = ax.legend((leg_list), ncol=2,shadow=False,bbox_to_anchor=[0.045, -0.22], loc='lower left',frameon=False)
            ltext  = leg.get_texts() 
            plt.setp(ltext, fontsize=fs)    # the legend text fontsize
                             
        show()  

                    
        return ax
                                                         