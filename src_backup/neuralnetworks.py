# -*- coding: utf-8 -*-
"""
Created on Sun Jul 23 15:33:23 2023

@author: christophe
"""


from keras import models
from keras import layers
from preprocessedseries import Preprocessed
from stresstimeseries import Stress
from headstimeseries import Heads
from logging import getLogger
import os
import numpy as np
from matplotlib.pyplot import figure, show,savefig,clabel,contour,setp,gcf,getp,gca,close
from utilities import basename
from modeldefinition import ModelDefinition

"""  
This module contains the functions to apply a neural network model based of the
keras library

It is still in development, so documentation is not ready yet

"""


# def fit_goodness(P,M,observed,prec):

#     e=residuals(P,M,observed,prec)
#     observed_deviations_from_mean=observed[:,1]-np.mean(observed[:,1])  
#     SSreg=sum(pow(e,2))
#     SStot=sum(pow(observed_deviations_from_mean,2))
#     explained_variance=100*(1-SSreg/SStot)
    
#     print('48 SSreg',SSreg)
#     print('49 SStot',SStot)
    
    
#     return e,observed,explained_variance



def explained_variance(observed,simulated):

    observed=observed.reshape(len(observed),)
    simulated=simulated.reshape(len(simulated),)
    
    e=observed-simulated
    observed_deviations_from_mean=observed-np.mean(observed)  
    SSreg=sum(pow(e,2))
    SStot=sum(pow(observed_deviations_from_mean,2))
    
    print('66 SSreg',SSreg)
    print('67 SStot',SStot)
    explained_variance=100*(1-SSreg/SStot)

    return explained_variance

def build_model():
    #Multi Layers Perceptron (2 hidden layers)
    model = models.Sequential()
    model.add(layers.Dense(64, activation='relu',input_shape=(train_data.shape[1],)))
    model.add(layers.Dense(64, activation='relu'))
    model.add(layers.Dense(1))
    model.compile(optimizer='rmsprop', loss='mse', metrics=['mae'])
    return model


def build_model2():
    #Multi Layers Perceptron (3 hidden layers)
    model = models.Sequential()
    model.add(layers.Dense(64, activation='relu',input_shape=(train_data.shape[1],)))
    model.add(layers.Dense(64, activation='relu'))
    model.add(layers.Dense(64, activation='relu'))
    model.add(layers.Dense(64, activation='relu'))
    model.add(layers.Dense(1))
    model.compile(optimizer='rmsprop', loss='mse', metrics=['mae'])
    return model


def build_model3():
    #Mrecurrent Neural Network
    from keras.layers import SimpleRNN
    model = models.Sequential()
    #model.add(SimpleRNN(32,input_shape=(train_data.shape[1],),return_sequences=True))
    n_steps=train_data.shape[1]
    n_features=1
    input_shape=(n_steps, n_features)
    model.add(layers.GRU(64, input_shape=input_shape))         
    #model.add(layers.LSTM(32, input_shape=(None, train_data.shape[-1])))
    model.add(layers.Dense(1))
    model.compile(optimizer='rmsprop', loss='mse', metrics=['mae'])
    return model




if __name__ == "__main__":

    #curdir=os.getcwd()
    model_definition = ModelDefinition().model_definition
    abs_path = os.path.dirname(os.path.abspath(__file__))
    abs_path_splitted = abs_path.split('\\')
    path_to_parent_folder_elements = abs_path_splitted[:-1]
    path_to_parent_folder = os.path.join(*path_to_parent_folder_elements)
    splitted = path_to_parent_folder.split(':')#trick  to  repair  path
    path_to_parent_folder = splitted[0]+':\\'+splitted[1]#trick  to  repair  path
    path_to_file_folder = os.path.join(path_to_parent_folder,'resources')  

    stresses_dict = {}
    
    path_to_file = os.path.join(path_to_file_folder,'260_De_Bilt_PREC_19010101_20200910.csv')
    key = basename(path_to_file)
    stresses_dict['prec'] = {}
    stresses_dict['prec'][key] = Stress.read_from_csv(path_to_file, stress_type = 'prec')    
    

    path_to_file=os.path.join(path_to_file_folder,'260_De_Bilt_EVAP_19010101_20200910.csv')
    key = basename(path_to_file)
    stresses_dict['evap'] = {}
    stresses_dict['evap'][key] = Stress.read_from_csv(path_to_file, stress_type = 'evap')    


    tminstr = '01-01-1900 08:00:00'
    tmaxstr = '31-12-2100 08:30:00'
    
    path_to_file_folder = os.path.join(path_to_parent_folder,'resources')  
    path_to_file = os.path.join(path_to_file_folder,'28AP0093_1.txt')
    heads = Heads.read_from_csv(path_to_file, tminstr = tminstr, tmaxstr = tmaxstr)

    arguments_dic = {}

    for key in stresses_dict:
        for e in stresses_dict[key]:
            stresses_dict[key][e].plot()


    time_step = 1.#/(24.*2)
    
    memory_dict = {'prec':365*5,
              'evap':365*5,
              'riv':7,
              'pump':365*5,
              'noise':30}      
    
    
    time_step = 1. #1./(24.*2)
    preprocessed = Preprocessed(heads = heads,stresses_dict = stresses_dict, 
                                time_step = time_step, memory_dict = memory_dict,tminstr = tminstr, tmaxstr = tmaxstr,model_definition = md).preprocess()
    time = preprocessed.time
    heads = preprocessed.heads
    stresses_dict = preprocessed.stresses_dict 
    Nint_dict = preprocessed.Nint_dict
    all_data_for_neural_networks = preprocessed.all_data_for_neural_networks
    all_targets_for_neural_networks = preprocessed.all_targets_for_neural_networks

    


    #normalize data
    # =============================================================================
    # mean = train_data.mean(axis=0)
    # train_data -= mean
    # std = train_data.std(axis=0)
    # train_data /= std
    # test_data -= mean
    # test_data /= std
    # =============================================================================
    
    #building network - model definition    
    

    test_split=0.2
    # all_data,all_targets,time_deep_learn,obs_deep_learn,Nstresses=ml.generate_deep_learning_input()
    
    x = all_data_for_neural_networks
    y = all_targets_for_neural_networks
    
    # method is copied from from boston_housing.py to be found in C:\Users\Gebruiker\anaconda3\Lib\site-packages\keras\datasets
    # =============================================================================
    # rng = np.random.RandomState(seed=None)
    # indices = np.arange(len(x))
    # rng.shuffle(indices)
    # x = x[indices]
    # y = y[indices]
    # =============================================================================
    
    x_train = np.array(x[:int(len(x) * (1 - test_split))])
    y_train = np.array(y[:int(len(x) * (1 - test_split))])
    x_test = np.array(x[int(len(x) * (1 - test_split)):])
    y_test = np.array(y[int(len(x) * (1 - test_split)):])
    
    
    
    # x_train = np.array(x[int(len(x) * (1 - test_split)):])
    # y_train = np.array(y[int(len(x) * (1 - test_split)):])
    # x_test = np.array(x[:int(len(x) * (1 - test_split))])
    # y_test = np.array(y[:int(len(x) * (1 - test_split))])
    
    
    
    train_data = x_train
    train_targets = y_train
    test_data = x_train
    test_targets = y_train

    
    
    

    model = build_model()
    
    
    
    
    num_epochs = 10
    model.fit(train_data, train_targets, validation_data=(test_data, test_targets),epochs=num_epochs, batch_size=1, verbose=0)
    val_mse, val_mae = model.evaluate(train_data, train_targets, verbose=0)
    print('1912 Mean Absolute Error Deep Learning',val_mae)
    y_predicted = np.empty(np.shape(all_targets_for_neural_networks))
    
    #model= build_model()
    y_predicted = model.predict(all_data_for_neural_networks)
    
    
    
    EV_Deep_Learning = explained_variance(all_targets_for_neural_networks,y_predicted)
    print('1921 EV_Deep_Learning',EV_Deep_Learning)
    
    
    
    obs_deep_learning=np.empty((len(time),2))
    obs_deep_learning[:,0] = time
    obs_deep_learning[:,1] = heads.interpolated_normalized[:,1] 
    
    model_deep_learning_MLP=np.empty((len(time),2))
    model_deep_learning_MLP[:,0] = time 
    model_deep_learning_MLP[:,1] = y_predicted.reshape(len(y_predicted),)                 
    
    # plot_time_series_models(model_deep_learning_MLP,model_name='Deep Learning Multi Layers Perceptron',obs=obs_deep_learning,targets=None,stress=None,stress_name=None,begin_date=None,end_date=None,response_func=None)
    
    fig=figure()
    ax = fig.add_subplot(111)
    ax.set_title('261 inspect deep learning results ')
    ax.plot_date(obs_deep_learning[:,0],obs_deep_learning[:,1],'b')
    ax.plot_date(model_deep_learning_MLP[:,0],model_deep_learning_MLP[:,1],'r')
    leg = ax.legend((["observed","MLP"]),loc='upper left',shadow=False)
    ax.set_xlabel('time')
    ax.grid(False)
    show()
    

    
