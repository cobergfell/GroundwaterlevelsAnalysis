# -*- coding: utf-8 -*-
"""
Created on Mon Oct 30 10:41:10 2023

@author: christophe
"""
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 28 20:09:40 2023

@author: christophe
"""
import os
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure, show,savefig
from datetime import datetime
from matplotlib.dates import date2num, num2date
from matplotlib.ticker import Formatter,FuncFormatter, NullLocator, FixedLocator
from logging import getLogger    
import numpy as np
import pandas as pd

logger = getLogger(__name__)


def export_heads_nparray_to_ascii(tseries, metadata = None, path = None, plot = True, tminstr = None, tmaxstr = None, suffix = '_'):#read head serie from own format
    

    
    """
    Export a numpy head time series to simple ascii file
    
    Parameters
    ----------
    tseries: numpy array_like
        array of time series of  values
        
    tminstr: string (optional)
            string specifying the minimum time on the plot (format dd-mm-yyyy)
            
    tmaxstr: string (optional)
            string specifying the maximum time on the plot (format dd-mm-yyyy)   
            
    path :  string
            Path to export directory
        
    plot : boolean
        A boolean flag to indicate that the times series is to be ploted (if plot == True)  
        
        
    Return
    ----------        
        
    export_path :  string
        Path to export file
    
    """    
        
        
    from matplotlib.pyplot import figure, show,savefig
    from datetime import datetime
    from matplotlib.dates import date2num, num2date
    from matplotlib.ticker import Formatter,FuncFormatter, NullLocator, FixedLocator
    from logging import getLogger    
    import numpy as np
    import pandas as pd
    import os


    name = metadata.name+suffix
    X = metadata.X
    Y = metadata.Y
    Z = metadata.Z
    surface_level = metadata.surface_level
    screen_top = metadata.screen_top
    screen_bottom = metadata.screen_bottom
    screen_number = metadata.screen_number

          
    N_data = len(tseries)
    if X is not None:
        Xstr = str('%8.2f' %X)
    else:
        Xstr = '0'
        
    if Y is not None:
        Ystr = str('%8.2f' %Y)
    else:
        Ystr = '0'        
        
    if screen_number is not None:
        screen_number_str = str('%8.2f' %screen_number)
    else:
        screen_number_str = '1'
        
    if surface_level is not None:
        surface_level_str = str('%8.2f' %surface_level)
    else:
        surface_level_str = '-999'          
        
        
    if screen_top is not None:
        screen_top_str = str('%8.2f' %screen_top)
    else:
        screen_top_str = '-999'        
        
    if screen_bottom is not None:
        screen_bottom_str = str('%8.2f' %screen_bottom)
    else:
        screen_bottom_str = '-999'            
    

    begin_date = num2date(tseries[0,0])
    begin_date_str = str(begin_date.day)+'-'+str(begin_date.month)+'-'+str(begin_date.year)
    end_date = num2date(tseries[-1,0])
    end_date_str = str(end_date.day)+'-'+str(end_date.month)+'-'+str(end_date.year)        
    header1 = 'name, screen_number, X, Y, surface_level_m_above_datum, screen_top_m_above_datum, screen_bottom_m_above_datum,begin_date,end_date'
    header2 = name+','+screen_number_str+','+Xstr+','+Ystr+','+surface_level_str+','+screen_top_str+','+screen_bottom_str+','+begin_date_str+','+end_date_str


    if path is not None:
        export_path = path+'//'+name+'.txt'
    else:
        curdir = os.getcwd()
        export_path = curdir+'//'+name+'.txt'
    
    F=open(export_path,'w')
    F.write(header1+'\n')
    F.write(header2+'\n')
    F.write('\n')
    F.write('Date(dd-mm-yyyy HH:MM:SS),head(m above datum)\n')
    for i in range(0,len(tseries)):
        t = num2date(tseries[i,0])
        value = tseries[i,1]
        record = str(t.day).zfill(2)+'-'+str(t.month).zfill(2)+'-'+str(t.year)+'  '+str(t.hour).zfill(2)+':'+str(t.minute).zfill(2)+':'+str(t.second).zfill(2)+','+str('%8.4f' %value)
        F.write(record+'\n')
    F.close()

    if plot == True:
        
        left_margin = 0.12
        right_margin = 0.06
        top_margin = 0.12
        bottom_margin = 0.44
        dy = 0.04
        dx = 0.05        
        
        fs = 12
         
        t = lambda t_str : date2num(datetime.strptime(t_str,'%d-%m-%Y')) # lamda function to transform time string into time num
        
        if tminstr is not None:
            tmin = t(tminstr)
        else:
            tmin = tseries[0,0]
                        
        if tmaxstr is not None:
            tmax = t(tmaxstr)
        else:
            tmax = tseries[-1,0]
        
        mask = (tseries[:,0] >= tmin) & (tseries[:,0] < tmax)
        tseries = tseries[mask]        
        
        
        yas_title = 'm above datum'
        
        plot_length = 1-left_margin-right_margin
        plot_height = 1-bottom_margin-top_margin
        xll = left_margin # x coordinate lower left corner
        yll = 1-top_margin-plot_height # y coordinate lower left corner
        
        
        fig=figure(num=None, figsize=(8,10),dpi=50,facecolor='w', edgecolor='k')# figsize=(width, height) in inches.
        
        ax = fig.add_axes([xll,yll,plot_length,plot_height],frameon=True, xscale=None, yscale=None)  
        
        leg_list = []
        ax.plot_date(tseries[:,0],tseries[:,1],'b-',linewidth = 1.0)     
        leg_list.append('observed')
                
        def format_date(dates, pos=None):
            return num2date(dates).strftime('%Y')
        
        # def format_date(dates, pos=None):
        #     return num2date(dates).strftime('%d-%m-%Y %H:%M')        
        
        ax.set_ylabel(yas_title,fontsize=fs)
          
        ax.set_xlabel('date (yyyy)',fontsize=fs,ha='center')
        ax.xaxis.set_major_formatter(FuncFormatter(format_date))
        ax.xaxis.set_tick_params(rotation=0, labelsize=fs)
         
        
        # for tick in ax.xaxis.get_major_ticks():
        #     tick.label.set_fontsize(fs)
        #     tick.label.set_rotation(0)
        #     # tick.label.set_rotation(45)
              
        
        for label in ax.xaxis.get_majorticklabels():
            #ha=label.get_horizontalalignment()
            label.set_horizontalalignment('center')
            
        time = tseries[:,0]
        
        
        lower_expect, median, upper_expect = generate_gxg(tseries, tminstr = tminstr, tmaxstr = tmaxstr, alpha = 0.05)
        
        
        upper_expect_array = upper_expect*np.ones(len(tseries))
        median_array = median*np.ones(len(tseries))
        lower_expect_array = lower_expect*np.ones(len(tseries))
        
        
        ax.plot_date(time,upper_expect_array,'c--',linewidth=2.0) 
        leg_list.append('95% upper quantile')
        ax.plot_date(time,median_array,'m-',linewidth=2.0)  
        leg_list.append('median')
        ax.plot_date(time,lower_expect_array,'g--',linewidth=2.0)         
        leg_list.append('95% lower quantile')
        
        
        leg = ax.legend((leg_list),ncol=3,bbox_to_anchor=(0.85,-0.24),frameon=False)
        
        Xleg = 0.05
        Yleg = -0.22
        
        ax.annotate('Legend', xy=(Xleg, Yleg), xycoords='axes fraction',horizontalalignment='left', verticalalignment='center',fontsize=fs)#this is for 1 column 2 rows
        leg = ax.legend((leg_list), ncol=2,shadow=False,bbox_to_anchor=[Xleg, Yleg-0.2], loc='lower left',frameon=False)
        ltext  = leg.get_texts() 
        plt.setp(ltext, fontsize=fs)    # the legend text fontsize
        
        title='Time series of water levels at  '+ name
        ax.set_title(title,fontsize=fs)       
        (x, y) = ax.title.get_position()
        #ax.title.set_y(0.95 * y)    
        ax.grid(True)               
        
        
        
        ax.set_ylabel('Meter above datum',fontsize=fs)
        ax.set_xlabel('Date (yyyy)',fontsize=fs,ha="center")
        #ax.set_xlabel('Datum (yyyy)',fontsize=fs)#030919
        ax.xaxis.set_major_formatter(FuncFormatter(format_date))
        #fig.autofmt_xdate()
        #ax.xaxis.set_label_position('right')
        
        ax.xaxis.set_tick_params(rotation=45, labelsize=fs)
  
        
        # for tick in ax.xaxis.get_major_ticks():
        #     tick.label.set_fontsize(fs)
        #     tick.label.set_rotation(45)
        
        for label in ax.xaxis.get_majorticklabels():
            #ha=label.get_horizontalalignment()
            label.set_horizontalalignment('center')          
        
        ax.grid(True)
    
        list_of_header_items=[
                'Location name',
                'X',
                'Y',
                '5% low quantile (m above datum)',
                'Median (m above datum)',
                '5% high quantile (m above datum)']
        
        
        
        list_of_values=[
        name,
        str('%10.2f' %X),
        str('%10.2f' %Y),
        str('%10.2f' %lower_expect),
        str('%10.2f' %median),
        str('%10.2f' %upper_expect)]    
        
    
        X = [Xleg,Xleg+0.5]
        Y = [Yleg-0.3]
        
        
        delta_y=0.07
        for i in range(1,len(list_of_header_items)):    
            last_y=Y[-1]                      
            Y.append(last_y-dy)
        
        for i in range(0,len(list_of_header_items)):
            ax.annotate(list_of_header_items[i], xy=(X[0],Y[i]), xycoords='axes fraction',horizontalalignment='left', verticalalignment='center')
            ax.annotate(list_of_values[i], xy=(X[1],Y[i]), xycoords='axes fraction',horizontalalignment='left', verticalalignment='center')
        show()
        

        if name is not None:
            figname = "plot_of_{}".format(name)
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
                    

    return export_path
    





def generate_gxg(tseries, tminstr= None, tmaxstr = None, alpha = 0.05):
    
    """
    This is the static method version of the generate_gxg method already defined in head timeseries.
    (for future, consider keep only the static method version)
    
    
    To be used to estimate some standard quantiles of  groundwater levels time series

    (gxg refers to the Dutch habit to call these means GLG, GG, GHG which stand
     for gemiddelde laag grondwaterstand, gemiddelde grondwaterstand, 
     and gemiddelde hoog grondwaterstand)

    Parameters
    ----------
    tseries: numpy array_like
        array of time series of  values
        
    tminstr: string (optional)
            string specifying the minimum time on the plot (format dd-mm-yyyy)
            
    tmaxstr: string (optional)
            string specifying the maximum time on the plot (format dd-mm-yyyy)    
        
    alpha: float
        quantile determining the low/high groundwater level
        


    Returns
    -------
    
    glg : float
        mean low groundwaterlevel ('gemiddelde lage grondwaterstand')
        
    gg: float
        mean  groundwaterlevel ('gemiddelde grondwaterstand')
        
    ghg : float
        mean high groundwaterlevel ('gemiddelde hoge grondwaterstand')
        
    """

    t = lambda t_str : date2num(datetime.strptime(t_str,'%d-%m-%Y')) # lamda function to transform time string into time num


    if tminstr is not None:
        tmin = t(tminstr)
    else:
        tmin = tseries[0,0]
        
        
    if tmaxstr is not None:
        tmax = t(tmaxstr)
    else:
        tmax = tseries[-1,0]       
  
        
    
    mask = (tseries[:,0] >= tmin) & (tseries[:,0] < tmax)
    tseries = tseries[mask]
    
    glg = None
    gg = None
    ghg = None


    gg = np.median(tseries[:,1])

    #sort the time series in order to determine the quantiles
    b = np.argsort(tseries[:,1])#sort from low to high
    b = b[::-1]#sort from high to low
    sorted_serie = tseries[b,:]
    N = len(tseries)

    index = int((1-alpha)*N)# groundwater levels were sorted from high to low
    glg = sorted_serie[index,1]
    
    index=int(alpha*N)# groundwater levels were sorted from high to low
    ghg = sorted_serie[index,1]
  
        
        
    if len(tseries) < 24:
        logger.warning("The selected time series is too short to estimates mean levels")        
        

    return glg,gg,ghg 



def generate_output_from_heads_list(Heads_list, export_model_goodness_of_fit = False, exportpath = None):

    """
    A method to generate an output string for a list of Head objects

    """  
    
    output_string = f'\n\nSummary of observed groundwater levels\n'
    
    for head in Heads_list:
    
        paddingright = 20
        justify = 16
          
            
    # create structured array
    if export_model_goodness_of_fit is True:
        dtype=[
            ('name', 'S20'),
            ('X', 'f4'),
            ('Y', 'f4'),
            ('Z', 'f4'),
            ('Surface level', 'f4'),
            ('GLG', 'f4'),
            ('GG', 'f4'),
            ('GHG', 'f4'),
            ('Expl.var.', 'f4')]  
    else:
        dtype=[
            ('name', 'S20'),
            ('X', 'f4'),
            ('Y', 'f4'),
            ('Z', 'f4'),
            ('Surface level', 'f4'),
            ('GLG', 'f4'),
            ('GG', 'f4'),
            ('GHG', 'f4')]              
      
    structured_array= np.empty(len(Heads_list), dtype=dtype)             

    for i in range(0,len(Heads_list)):
        head = Heads_list[i]
        structured_array[i][0] = head.name
        structured_array[i][1] = head.X
        structured_array[i][2] = head.Y
        structured_array[i][3] = head.Z
        structured_array[i][4] = head.surface_level
        structured_array[i][5] = head.glg
        structured_array[i][6] = head.gg
        structured_array[i][7] = head.ghg
        if export_model_goodness_of_fit is True:
            try: 
                structured_array[i][8] = head.fitgoodness['expvar']
            except:
                structured_array[i][8] = None     

            
   
     # create pandas dataframe as alternative to structured array, we keep both options for the moment

    columns=[
    'Name', 
    'X', 
    'Y', 
    'Z', 
    'Surface level',
    'GLG',
    'GG', 
    'GHG']
    if export_model_goodness_of_fit is True:
        columns.append('Expvar')
 
    df = pd.DataFrame(columns=columns) 
 
    data = [] 
    for i in range(0,len(Heads_list)):
        head = Heads_list[i]
        new_record = []
        new_record.append(head.name)
        new_record.append(head.X)
        new_record.append(head.Y)
        new_record.append(head.Z)
        new_record.append(head.surface_level)
        new_record.append(head.glg)
        new_record.append(head.gg)
        new_record.append(head.ghg)         
        if export_model_goodness_of_fit is True:
            new_record.append(head.fitgoodness['expvar'])

        
    data.append(new_record)
    
    df = pd.DataFrame(data,columns=columns) 
   
            
    output_string_csv = f'name,X,Y,Z,Surface_level,GLG,GG,GHG'  
    if export_model_goodness_of_fit is True:
        output_string_csv += ',Expvar\n'  
    else: 
        output_string_csv += '\n'  

    for i in range(0,len(Heads_list)):
        head = Heads_list[i]
        output_string_csv += f'{head.name},'
        output_string_csv += f'{head.X:.2f},'
        output_string_csv += f'{head.Y:.2f},'
        output_string_csv += f'{head.Z:.2f},'
        output_string_csv += f'{head.surface_level:.2f},'
        output_string_csv += f'{head.glg:.2f},'
        output_string_csv += f'{head.gg:.2f},'
        output_string_csv += f'{head.ghg:.2f}'    
        if export_model_goodness_of_fit is True:
            expvar=head.fitgoodness['expvar']
            output_string_csv += f',{expvar:.2f}\n'
        else:
            output_string_csv += f'\n'    

    if exportpath is None:
        abs_path = os.path.dirname(os.path.abspath(__file__))
        abs_path_splitted = abs_path.split('\\')
        path_to_parent_folder_elements = abs_path_splitted[:-1]
        path_to_parent_folder = os.path.join(*path_to_parent_folder_elements)
        splitted=path_to_parent_folder.split(':')#trick  to  repair  path
        path_to_parent_folder = splitted[0]+':\\'+splitted[1]#trick  to  repair  path
        dirpath = os.path.join(path_to_parent_folder,'export')
        
        # curdir = os.getcwd()
        # curdir = os.getcwd()
        # exportpath = curdir+'//'+'summary_heads_time_series.txt'
        
        dirpath = dirpath
        
        if not os.path.exists(dirpath):
           os.mkdir(dirpath)

    exportpath =  os.path.join(dirpath,'summary_heads_time_series.txt')
    
    F=open(exportpath,'w')
    F.write(output_string_csv)
    F.close()   
    
    print(output_string_csv)

    exportpath =  os.path.join(dirpath,'summary_heads_time_series.xlsx')
    df.to_excel(exportpath) 
        
        
    return output_string




def multipleplots(Heads_list , share_axes = True, plot_title = None, yas_title = None, 
         tminstr = '01-01-1900', tmaxstr = '31-12-2100', save_plot = False):
    
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
    
    left_margin=0.12
    right_margin=0.06
    top_margin=0.12
    bottom_margin=0.14
    dy=0.14
    dx=0.05        

    fs = 12
    colors = ['b-','m-','y-','c-','g-', 'b--','m--','y--','c--','g--','b*','m*','y*','c*','g*']
        
    
    t = lambda t_str : date2num(datetime.strptime(t_str,'%d-%m-%Y')) # lamda function to transform time string into time num

    tmin = t(tminstr)
    tmax = t(tmaxstr)



    nl = len(Heads_list)
    nseriesmax = 10
    if nl > nseriesmax:
        message = (f'The list of time series to plot should not exceed {nseriesmax}.\n')                        
        logger.warning(message) 
        

    if plot_title is None:  
        plot_title = ''

          
    if yas_title is None:  
        yas_title = 'm above datum'



        
    plot_length = (1-left_margin-right_margin-dx)
    plot_height = (1-bottom_margin-top_margin-2*dy)
    x_lower_left = left_margin
    y_lower_left = 1-top_margin-plot_height

    fig=figure(num=None, figsize=(8,8),dpi=50,facecolor='w', edgecolor='k')# figsize=(width, height) in inches.
    ax = fig.add_axes([x_lower_left,y_lower_left,plot_length,plot_height],frameon=True, xscale=None, yscale=None)  

    leglist = []
    for i in range(0,nl): 
        heads = Heads_list[i]
        tseries = heads.observed
        mask = (tseries[:,0] >= tmin) & (tseries[:,0] < tmax)
        tseries = tseries[mask]          
        ax.plot_date(tseries[:,0],tseries[:,1],colors[i],linewidth = 1.0)   
        leglist.append(heads.name)

    ax.set_title(plot_title,fontsize=fs)
    (x, y) = ax.title.get_position()
    #ax.title.set_y(0.95 * y)    
    ax.grid(True) 


    longtime = 365 * 3
    mediumtime = 365 * 1
    shorttime = 32
    deltime = tmax-tmin
    
    if deltime > longtime:    
        def format_date(dates, pos=None):
            return num2date(dates).strftime('%Y')
        
    elif (deltime>shorttime) & (deltime<mediumtime):
        def format_date(dates, pos=None):
            return num2date(dates).strftime('%d-%m-%Y')     
    else:
        def format_date(dates, pos=None):
            return num2date(dates).strftime('%d-%m-%Y %H:%M')  
        
    

    ax.set_ylabel(yas_title,fontsize=fs)
  
    ax.set_xlabel('date (yyyy)',fontsize=fs,ha='center')
    ax.xaxis.set_major_formatter(FuncFormatter(format_date))

    ax.xaxis.set_tick_params(rotation=0, labelsize=fs)
  
        
    # for tick in ax.xaxis.get_major_ticks():
    #     tick.label.set_fontsize(fs)
    #     tick.label.set_rotation(0)
    #     # tick.label.set_rotation(45)
    
    # for label in ax.xaxis.get_majorticklabels():
    #     #ha=label.get_horizontalalignment()
    #     label.set_horizontalalignment('center')
        
    Xleg = 0.05
    Yleg = -0.22
    
    xul = Xleg # x upper left of bbox
    yul = Yleg-0.05 # y upper left of bbox      
    
    nitems = len(Heads_list)
    
    ncol = 3
    
    nrows = nitems // ncol # // is the floor division
    if nitems % ncol > 0: # if a row is not full
        nrows +=1
        
    xll = xul    
    yll = yul-nrows*0.08*14.0/fs
         
  
    # https://stackoverflow.com/questions/39803385/what-does-a-4-element-tuple-argument-for-bbox-to-anchor-mean-in-matplotlib
    # In bbox_to_anchor=(xll, yll, width, height), xll,yll are the lower left corner coordinates of the bounding box.

    ax.annotate('Legend', xy=(Xleg, Yleg), xycoords='axes fraction',horizontalalignment='left', verticalalignment='center',fontsize=fs)#this is for 1 column 2 rows
    leg = ax.legend((leglist), ncol=ncol,shadow=False,bbox_to_anchor=[xll,yll], loc='lower left',frameon=False,fontsize=fs)
    ltext  = leg.get_texts() 
    plt.setp(ltext, fontsize=fs)    # the legend text fontsize
    
    show()            
    


    if save_plot == True:
        if plot_title is not None:
            figname = "plot_of_{}".format(plot_title)
            f'{plot_title}'
            
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

