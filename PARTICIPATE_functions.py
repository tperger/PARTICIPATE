# -*- coding: utf-8 -*-
"""
Created on Thu Feb 10 14:23:50 2022

@author: perger
"""
import pandas as pd
import numpy as np
import math
from pathlib import Path
import statistics
import pyam
from tslearn.clustering import TimeSeriesKMeans

def custom_reshape(timeseries, hours):
    """ Auxilary function for k-means"""
    datapoints = int(len(timeseries)/hours)
    return np.reshape(timeseries[:(datapoints*hours)],
                      (datapoints, hours))

def custom_norm(timeseries):
    """ Auxilary function for k-means"""
    norm = np.amax(timeseries.to_numpy())
    if norm > 0:
        norm_array = timeseries.to_numpy()/norm
    else:
        norm_array = timeseries.to_numpy()
    return norm_array, norm

def cluster_input(prosumer, grid_data, load, PV, k, hours):
    """ Clustering of hourly input data into representative time periods using 
    k-means clustering algorithm
    
    Keyword arguments:
        prosumer ... list of prosumer names
        emissions ... marginal emissions (DataFrame(index=time_steps,
                                                    columns='Emissions'))
        p_G_out ... DA prices (DataFrame(index=time_steps,
                                         columns='DA')) 
        load ... electricity demand (DataFrame(index=time_steps,
                                               columns=prosumer))
        PV ... PV generation (DataFrame(index=time_steps,
                                        columns=prosumer))
        k ... number of clusters in the k means algorithm / number of 
              representative time periods (e.g. days, weeks, ...)
        hours ... number of hours in one representative time period
    """
    
    # Normalize and create input matrix (cluster_matrix) for k-means algorithm
    cluster_matrix = []
    norm = []
    
    # emissions
    _input, _norm = custom_norm(grid_data['Emissions'])
    cluster_matrix = custom_reshape(_input, hours)
    norm.append(_norm)
    
    # p_G_out
    _input, _norm = custom_norm(grid_data['DA'])
    cluster_matrix = np.append(cluster_matrix,
                               custom_reshape(_input, hours), 
                               axis=1)
    norm.append(_norm)
    
    # p_G_in
    _input, _norm = custom_norm(grid_data['Residential'])
    cluster_matrix = np.append(cluster_matrix,
                               custom_reshape(_input, hours), 
                               axis=1)
    norm.append(_norm)
    
    # load
    for i in prosumer:
        _input, _norm = custom_norm(load[i])
        cluster_matrix = np.append(cluster_matrix,
                                   custom_reshape(_input, hours), 
                                   axis=1)
        norm.append(_norm)
        
    # PV
    for i in prosumer:
        _input, _norm = custom_norm(PV[i])
        cluster_matrix = np.append(cluster_matrix,
                                   custom_reshape(_input, hours), 
                                   axis=1)
        norm.append(_norm)
    
    # Run k-means algorithm
    kmeans = TimeSeriesKMeans(n_clusters=k).fit(cluster_matrix)
    
    # Results of k-means
    result_kmeans = np.array(kmeans.cluster_centers_)[:,:,0]
    labels = np.array(kmeans.labels_)
    (unique, counts) = np.unique(labels, return_counts=True)
    
    # Scale k-means results according to weights 
    #for i in unique:
    #    result_kmeans[i,:] *= counts[i] / (sum(counts[j] for j in unique))
    
    # Re-scale and create data frames with representative times
    j = 0
    grid_data_cluster = pd.DataFrame()
    load_cluster = pd.DataFrame()
    PV_cluster = pd.DataFrame()
    
    grid_data_cluster['Emissions'] = np.reshape(
        result_kmeans[:,(j*hours):((j+1)*hours)],
        k*hours) * norm[j]    
    j += 1
    
    grid_data_cluster['DA'] = np.reshape(
        result_kmeans[:,(j*hours):((j+1)*hours)],
        k*hours) * norm[j]                                   
    j += 1
    
    grid_data_cluster['Residential'] = np.reshape(
        result_kmeans[:,(j*hours):((j+1)*hours)],
        k*hours) * norm[j]                                   
    j += 1
        
    for i in prosumer:
        load_cluster[i] = np.reshape(
            result_kmeans[:,(j*hours):((j+1)*hours)],
            k*hours) * norm[j]
        j += 1
            
    for i in prosumer:
        PV_cluster[i] = np.reshape(
            result_kmeans[:,(j*hours):((j+1)*hours)],
            k*hours) * norm[j]
        j += 1
    
    time_clustered = grid_data_cluster.index
        
    return (grid_data_cluster, load_cluster, PV_cluster, 
            time_clustered, counts)

def define_community(settlement_pattern=None,
                     buildings_per_SP=None,
                     model_name=None,   
                     scenario_name=None,
                     region_name=None,
                     year=None,
                     clustering=True):
    
    
    """
    Settlement pattern options: city, town, suburban, rural
    Region options: Austria, Greece, Norway, Spain, UK
        
    """
    
    if region_name not in ['Austria', 'Greece', 'Norway', 'Spain', 'UK']:
        raise Exception('Selected country not in list of available countries')
                
    # Read Input Data (from the IAMC Format)

    # input data of prosumer
    
    PATH_FILES = Path(__file__).parent / 'Community data'
            
    n = buildings_per_SP[settlement_pattern]
    prosumer = (['Prosumer SH '+str(i+1) for i in range(n['SH'])]
                + ['Prosumer SAB '+str(i+1) for i in range(n['SAB'])]
                + ['Prosumer LAB '+str(i+1) for i in range(n['LAB'])]
                )    
                   
    # IAMC variable names: Electricity demand, PV generation, other prosumer data
    # load_var = 'Final Energy|Residential and Commercial|Electricity'
    # PV_var = 'Secondary Energy|Electricity|Solar|PV'
    SoC_max = 'Maximum Storage|Electricity|Energy Storage System'
    SoC_min = 'Minimum Storage|Electricity|Energy Storage System'
    q_bat_max = 'Maximum Charge|Electricity|Energy Storage System'
    q_bat_min = 'Maximum Discharge|Electricity|Energy Storage System'
    PV_capacity = 'Maximum Active power|Electricity|Solar'
    w = 'Price|Carbon'
    prosumer_var = [w, SoC_max, SoC_min, q_bat_max, q_bat_min, PV_capacity]
    
    load = pd.DataFrame()
    PV = pd.DataFrame()
    prosumer_data = pd.DataFrame()
    grid_data = pd.DataFrame()
    
    # Prosumer data
    for i in prosumer:
        _filename = i+'.csv'
        _df = pyam.IamDataFrame(PATH_FILES / _filename)
        _data = (_df
            .filter(region=region_name)
            .filter(model=model_name)
            .filter(scenario=scenario_name)
            .filter(year=year))
        load[i] = (_data
            .filter(
                variable='Final Energy|Residential and Commercial|Electricity')
            .as_pandas().set_index('time').value)
        PV[i] = (_data
            .filter(
                variable='Secondary Energy|Electricity|Solar|PV')
            .as_pandas().set_index('time').value)
        # prosumer data DataFrame
        prosumer_data[i] = (_data
            .filter(variable=prosumer_var)
            .as_pandas().set_index('variable').value)
        
    # Grid data
    _df = pyam.IamDataFrame(data='Grid_data.csv', sep=';')
    _data = (_df
            .filter(region=region_name)
            .filter(model=model_name)
            .filter(scenario=scenario_name)
            .filter(year=year))
    grid_data['DA'] = (_data
                 .filter(variable='Price|Secondary Energy|Electricity')
                 .as_pandas().set_index('time').value/1000) # price EUR/kWh
    grid_data['Residential'] = (_data
              .filter(
                  variable='Price|Final Energy|Residential|Electricity')['value']
              .values[0]/1000) # price EUR/kWh
    grid_data['Emissions'] = (_data
                 .filter(variable='Emissions|CO2')
                 .as_pandas().set_index('time').value)
    grid_data['Emissions'] = grid_data['Emissions'] + 0.1
    
    time_steps = load.index.tolist()
    
    if clustering:
        k = 3 # number of representative days
        hours = 24 # time steps of representative days       
        grid_data, load, PV, time_steps, counts = cluster_input(
                                                        prosumer, 
                                                        grid_data,
                                                        load, 
                                                        PV, 
                                                        k, 
                                                        hours)
        _data = np.repeat(counts, k*[hours])
        weight = pd.DataFrame(_data, index=time_steps, columns=['weight'])
    else:
        _data = [1]*8760
        weight = pd.DataFrame(_data, index=time_steps, columns=['weight'])
        
    # Other values
    # _file_name = 'Distances_'+settlement_pattern+'.csv'
    # distances = pd.read_csv(PATH_FILES / _file_name,
    #                       sep=';',
    #                       header=0, 
    #                       index_col='Prosumer')
    
    M = len(prosumer)
    
    b = np.zeros((M,1), int)
    U = np.random.uniform(low=0, high=1.0, size=(M, M))
    S = np.tril(U) + np.tril(U, -1).T
    np.fill_diagonal(S, b)
    
    distances = pd.DataFrame(index=prosumer, 
                             columns=prosumer, 
                             data=S)
    distances=distances.round(2)
        
    return load, PV, prosumer_data, grid_data, weight, distances