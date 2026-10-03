# -*- coding: utf-8 -*-
"""
Created on Sun Jun 11 14:08:30 2023

@author: christophe
"""


from matplotlib.pyplot import figure, show,savefig,clabel,contour,setp,gcf,getp,gca
from matplotlib.dates import date2num, num2date,DayLocator,MonthLocator,YearLocator,DateFormatter #HourLocator,Yearlocator
from datetime import datetime
from matplotlib.ticker import Formatter,FuncFormatter, NullLocator, FixedLocator
from matplotlib import rcParams
from utilities import *
from matplotlib import pyplot as plt



rcParams['text.usetex'] = False#True



""" 
Description
----------
This module contains plot utilities

It is still in construction
    
""" 


def simpleplot(tseries, tminstr = '01-01-1900', tmaxstr = '31-12-2100', Yticks = None, 
         tseriesname = None, title = None):
    

    """ 
    Description
    ----------
    a simple function to plot a time serie
    
    Parameters
    ----------    

    tseries: numpy array_like
            the time series to be plotted
            The time numbers in colomn 0 are time numbers as defined in matplotlib.dates 
            function date2num

    tminstr: string (optional)
            string specifying the minimum time on the plot (format dd-mm-yyyy)
            
    tmaxstr: string (optional)
            string specifying the maximum time on the plot (format dd-mm-yyyy)     
            
    Yticks: numpy array_like (optional)
            Y axis tick labels                   

    Returns
    -------
    ax: Matplotlib Axes object
            The Matplotlib Axes object    
    
    """ 
            
   

    fs = 12
    
    if tseriesname == None:
         tseriesname = 'tseries'
         
    if title == None:
         title = 'tseries'    

    t = lambda t_str : date2num(datetime.strptime(t_str,'%d-%m-%Y'))
    tmin = max(t(tminstr),tseries[0,0])  
    tmax = min(t(tmaxstr),tseries[-1,0])  
 
    left_margin=0.12
    right_margin=0.06
    top_margin=0.05
    bottom_margin=0.14
    dy=0.05
    dx=0.05

    nl = 1 #number of rows
    nc = 1 #number of columns
    plot_length = (1-left_margin-right_margin-dx)/nc
    plot_height = (1-bottom_margin-top_margin-2*dy)/nl
    x1 = left_margin
    x2 = x1+plot_length+dx # in case of double columns 
    X = [x1,x2,x1,x2,x1,x2,x1,x2]#list of x coordinates of bottom left corners
    y1 = 1-top_margin-plot_height
    y2 = y1-plot_height-dy # in case of 2 rows
    y3 = y2-plot_height-dy # in case of 3 rows
    y4 = y3-plot_height-dy # in case of 4 rows
    Y = [y1,y1,y2,y2,y3,y3,y4,y4]#list of x coordinates of bottom left corners    
    
    
    fig=figure(num=None, figsize=(8,8),dpi=50,facecolor='w', edgecolor='k')#figures with all plots in a colomn
    #ax = fig.add_subplot(111)
    ax=fig.add_axes([X[0],Y[0],plot_length,plot_height],frameon=True, xscale=None, yscale=None)    #where rect=[left, bottom, width, height] in normalized (0,1) units. axisbg is the background color for the axis, default white


    mask = (tseries[:,0] >= tmin) & (tseries[:,0] < tmax)
    tseries = tseries[mask]
    
    
    legend_list = []
    ax.plot_date(tseries[:,0],tseries[:,1],'b')
    legend_list.append(tseriesname)

 
    
    def format_pump_rate(x):
        return '%1.0f' % (x*1e-3)

    def format_date(dates, pos=None):
        return num2date(dates).strftime('%d-%m-%Y')
    
    def format_date(dates, pos=None):
        return num2date(dates).strftime('%Y')
    

    ax.set_ylabel('head (m above datum)',fontsize=fs)
    ax.set_ylabel('tseries y values',fontsize=fs)
  
    # ax.set_xlabel('date (dd-mm-yyyy)',fontsize=fs, rotation=20)
    ax.set_xlabel('date (yyyy)',fontsize=fs)
    ax.xaxis.set_major_formatter(FuncFormatter(format_date))
    #fig.autofmt_xdate()
    #ax.xaxis.set_label_position('right')
    
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(fs)
        tick.label.set_rotation(0)
        # tick.label.set_rotation(45)
    
    for label in ax.xaxis.get_majorticklabels():
        #ha=label.get_horizontalalignment()
        label.set_horizontalalignment('center')
            

    ax.set_title(title,fontsize=fs)
    
    
    (x, y) = ax.title.get_position()
    ax.title.set_y(0.87 * y)    

    ax.grid(True)
    
    # ax.fmt_xdata = DateFormatter('%Y-%m-%d')                        # %H:%M:%S
    # fig.autofmt_xdate()
    
                    

 
    ax.annotate('Legend', xy=(0.05, -0.15), xycoords='axes fraction',horizontalalignment='left', verticalalignment='center',fontsize=fs)
    leg = ax.legend((legend_list), ncol=3,shadow=False,bbox_to_anchor=[0.045, -0.25], loc='lower left',frameon=False)
    ltext  = leg.get_texts() 
    plt.setp(ltext, fontsize=fs)    # the legend text fontsize

    show()
    curdir=os.getcwd()
    #savefig(curdir+'\\simulations', dpi=None, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None,transparent=False, bbox_inches=None, pad_inches=0.1)


    return ax    



def multipleplot(tseries_list, tminstr = '01-01-1900', tmaxstr = '31-12-2100', Yticks = None, 
         tseriesname = None, title = None):
    

    """ 
    Description
    ----------
    a simple function to plot a time serie
    
    Parameters
    ----------    

    tseries: numpy array_like
            the time series to be plotted
            The time numbers in colomn 0 are time numbers as defined in matplotlib.dates 
            function date2num

    tminstr: string (optional)
            string specifying the minimum time on the plot (format dd-mm-yyyy)
            
    tmaxstr: string (optional)
            string specifying the maximum time on the plot (format dd-mm-yyyy)     
            
    Yticks: numpy array_like (optional)
            Y axis tick labels                   

    Returns
    -------
    ax: Matplotlib Axes object
            The Matplotlib Axes object    
    
    """ 
            
   

    fs = 12
    
    if tseriesname == None:
         tseriesname = 'tseries'
         
    if title == None:
         title = 'tseries'    

    t = lambda t_str : date2num(datetime.strptime(t_str,'%d-%m-%Y'))
    tmin = max(t(tminstr),tseries[0,0])  
    tmax = min(t(tmaxstr),tseries[-1,0])  
 
    left_margin=0.12
    right_margin=0.06
    top_margin=0.05
    bottom_margin=0.14
    dy=0.05
    dx=0.05

    nl = 1 #number of rows
    nc = 1 #number of columns
    plot_length = (1-left_margin-right_margin-dx)/nc
    plot_height = (1-bottom_margin-top_margin-2*dy)/nl
    x1 = left_margin
    x2 = x1+plot_length+dx # in case of double columns 
    X = [x1,x2,x1,x2,x1,x2,x1,x2]#list of x coordinates of bottom left corners
    y1 = 1-top_margin-plot_height
    y2 = y1-plot_height-dy # in case of 2 rows
    y3 = y2-plot_height-dy # in case of 3 rows
    y4 = y3-plot_height-dy # in case of 4 rows
    Y = [y1,y1,y2,y2,y3,y3,y4,y4]#list of x coordinates of bottom left corners    
    
    
    fig=figure(num=None, figsize=(8,8),dpi=50,facecolor='w', edgecolor='k')#figures with all plots in a colomn
    #ax = fig.add_subplot(111)
    ax=fig.add_axes([X[0],Y[0],plot_length,plot_height],frameon=True, xscale=None, yscale=None)    #where rect=[left, bottom, width, height] in normalized (0,1) units. axisbg is the background color for the axis, default white


    mask = (tseries[:,0] >= tmin) & (tseries[:,0] < tmax)
    tseries = tseries[mask]
    
    
    legend_list = []
    ax.plot_date(tseries[:,0],tseries[:,1],'b')
    legend_list.append(tseriesname)

 
    
    def format_pump_rate(x):
        return '%1.0f' % (x*1e-3)

    def format_date(dates, pos=None):
        return num2date(dates).strftime('%d-%m-%Y')
    
    def format_date(dates, pos=None):
        return num2date(dates).strftime('%Y')
    

    ax.set_ylabel('head (m above datum)',fontsize=fs)
    ax.set_ylabel('tseries y values',fontsize=fs)
  
    # ax.set_xlabel('date (dd-mm-yyyy)',fontsize=fs, rotation=20)
    ax.set_xlabel('date (yyyy)',fontsize=fs)
    ax.xaxis.set_major_formatter(FuncFormatter(format_date))
    #fig.autofmt_xdate()
    #ax.xaxis.set_label_position('right')
    
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(fs)
        tick.label.set_rotation(0)
        # tick.label.set_rotation(45)
    
    for label in ax.xaxis.get_majorticklabels():
        #ha=label.get_horizontalalignment()
        label.set_horizontalalignment('center')
            

    ax.set_title(title,fontsize=fs)
    
    
    (x, y) = ax.title.get_position()
    ax.title.set_y(0.87 * y)    

    ax.grid(True)
    
    # ax.fmt_xdata = DateFormatter('%Y-%m-%d')                        # %H:%M:%S
    # fig.autofmt_xdate()
    
                    

 
    ax.annotate('Legend', xy=(0.05, -0.15), xycoords='axes fraction',horizontalalignment='left', verticalalignment='center',fontsize=fs)
    leg = ax.legend((legend_list), ncol=3,shadow=False,bbox_to_anchor=[0.045, -0.25], loc='lower left',frameon=False)
    ltext  = leg.get_texts() 
    plt.setp(ltext, fontsize=fs)    # the legend text fontsize

    show()
    curdir=os.getcwd()
    #savefig(curdir+'\\simulations', dpi=None, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None,transparent=False, bbox_inches=None, pad_inches=0.1)


    return ax    



    def plot_multiple_series(tseries_list = None, legend_list = None, share_axes = True, plot_title = None, yas_title = None, 
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
        

        left_margin=0.12
        right_margin=0.06
        top_margin=0.12
        bottom_margin=0.14
        dy=0.14
        dx=0.05        

        fs = 14
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

