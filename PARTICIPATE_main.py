# -*- coding: utf-8 -*-
"""
Created on Thu May 27 16:59:44 2021

@author: perger
"""

import numpy as np
from numpy import matlib
import pandas as pd
from pyomo.environ import *
from pathlib import Path
import time

import PARTICIPATE_functions as pf
import PARTICIPATE_KKT_stochastic
import PARTICIPATE_KKT_deterministic  
import PARTICIPATE_plots
import FRESH_LP

solver_name = 'gurobi'
battery = True

country = 'Austria'
scenario_name = 'Default scenario'
model_name = 'FRESH:COM v2.0'

building_types = ['SH', 'SAB', 'LAB']

settlement_pattern = 'rural'

buildings_per_SP ={'city': {'SH': 0, 'SAB': 0, 'LAB': 10}, #city 
                   'town': {'SH': 0, 'SAB': 10, 'LAB': 0}, #town
                   'suburban': {'SH': 10, 'SAB': 0, 'LAB': 2}, #suburban
                   'rural': {'SH': 10, 'SAB': 8, 'LAB': 0, 'SME': 2} #rural
                   }  

load, PV, prosumer_data, grid_data, weight, distances = pf.define_community(
    settlement_pattern,
    buildings_per_SP,
    model_name,   
    scenario_name,
    country,
    year=2019,
    clustering=True)

PV_peak = [5,5,0,3,5,5,0,3,5,0,
           8,8,0,5,8,8,0,5,8,0]        
load *= [3335.59,4537.53,5252.55,5823.80,6336.76,
          6833.25,7346.21,7917.46,8632.48,9834.42,
          6257.90,8512.85,9854.29,10926.02,11888.38,
          12819.84,13782.20,14853.93,16195.37,18450.32]
storage = [1,0,0,1,0,1,0,1,0,0,
            1,0,0,1,0,1,0,1,0,0] 
prosumer_data.loc['Maximum Storage|Electricity|Energy Storage System'] *= storage
prosumer_data.loc['Maximum Charge|Electricity|Energy Storage System'] *= storage
prosumer_data.loc['Maximum Discharge|Electricity|Energy Storage System'] *= storage
prosumer_data.loc['Maximum Storage|Electricity|Energy Storage System'] *= 3
prosumer_data.loc['Maximum Charge|Electricity|Energy Storage System'] *= 2
prosumer_data.loc['Maximum Discharge|Electricity|Energy Storage System'] *= 2 

load /= 1000
PV *= PV_peak

N = 5
years = np.arange(1,N+1).tolist()

# Generate scenarios
num_scen = 4
scenarios = ['Scenario '+str(i) for i in range(1,num_scen+1)]

prosumer = load.columns.tolist() 

# x_0 = [1,2,3,5,2,1,0,0,2,2,
#        1,2,3,5,2,1,0,0,2,2]
x_0 = [3,3,3,0,0,2,2,0,0,0,
       0,0,2,2,2,0,3,0,3,0]
x_0 = dict(zip(prosumer, x_0))
prosumer = load.columns.tolist() 
values = [1,1,1,0,1,1,1,0,0,0, # SH 5 will dabei sein
          0,0,0,1,1,0,0,1,1,1] # for s_1
s_1 = dict(zip(prosumer, values))

s = {}
for w in scenarios:
    s[w] = {}
    for n in years[1:]:
        if n == years[1]:
            _values = values 
        else:
            if w == scenarios[0]:
                _values = [1,1,1,0,1,1,1,0,0,0,
                          1,1,1,1,1,1,1,1,1,1]
            elif w == scenarios[1]:
                _values = [1,1,1,0,0,1,1,0,0,0,
                          0,0,0,0,0,0,0,0,0,0]
            elif w == scenarios[2]:
                _values = [1,1,1,1,1,1,1,1,1,1,
                          0,0,0,1,1,0,0,1,1,1]
            elif w == scenarios[3]:
                _values = [0,0,0,0,0,0,0,0,0,0,
                          0,0,0,1,1,0,0,1,1,1]
        s[w][n] = dict(zip(prosumer, _values))
           
      

df = prosumer_data.copy(deep=True)

results, q_share_total, community_welfare = FRESH_LP.run_LP(
    load, PV, df, grid_data, weight, 
    distances, battery, solver_name, sharing=False)

emissions_old = results['emissions']
costs_old = results['costs']

# run model without forcast


scen = 'Scenario 1'
print("Solving deterministic problem ...")
t_start = time.time()
 
b_det, u_1_det, x_1_det, emissions_det = PARTICIPATE_KKT_deterministic.run_KKT(
          load, PV, df, grid_data, weight, 
          distances, emissions_old, battery, solver_name, [1], s_1, x_0, [scen])

print("... took ", str(time.time() - t_start), " sec")

b_det_scen = {}
u_1_det_scen = {}
x_1_det_scen = {}
emissions_det_scen = {}

for scen in scenarios:

    b_det = {}
    u_1_det = {}
    x_1_det = {}
    emissions_det = {}
    
    for n in years:   
        if n == 1:
            _x_0 = x_0
            _s_1 = s_1
        else:
            _x_0 = _x_1
            _s_1 = s[scen][n]
        _b, _u_1, _x_1, _emissions = PARTICIPATE_KKT_deterministic.run_KKT(
            load, PV, df, grid_data, weight, 
            distances, emissions_old, battery, solver_name, [n], _s_1, _x_0, [scen])
        b_det[n] = _b[scen][n]
        u_1_det[n] = _u_1
        x_1_det[n] = _x_1
        emissions_det[n] = _emissions
        
    b_det_scen[scen] = b_det
    u_1_det_scen[scen]  = u_1_det
    x_1_det_scen[scen]  = x_1_det
    emissions_det_scen[scen]  = emissions_det

print("Solving stochastic problem ...")
t_start = time.time()
 
b, u_1, x_1, x, emissions_new, emissions_scen = PARTICIPATE_KKT_stochastic.run_KKT(
    load, PV, df, grid_data, weight, distances, emissions_old, 
    battery, solver_name, years, s, s_1, x_0)

print("... took ", str(time.time() - t_start), " sec")

# run model with forecast

b_sto = {}
u_1_sto= {}
x_1_sto = {}
emissions_sto = {}

for n in years:   
    if n == 1:
        _x_0 = x_0
        _s_1 = s_1
        print('n is 1')
    else:
        _x_0 = x_1
        _s_1 = s['Scenario 1'][n]
        print('n is',n)
    if len(years[n-1:]) == 1:
        b, u_1, x_1, emissions = PARTICIPATE_KKT_deterministic.run_KKT(
            load, PV, df, grid_data, weight, 
            distances, emissions_old, battery, solver_name, [n], _s_1, _x_0)
        print('use KKT deterministic')
    else: 
        b, u_1, x_1, emissions = PARTICIPATE_KKT_stochastic.run_KKT(
            load, PV, df, grid_data, weight, 
            distances, emissions_old, battery, solver_name, 
            years[n-1:], s, _s_1, _x_0)
    b_sto[n] = b['Scenario 1'][n]
    u_1_sto[n] = u_1
    x_1_sto[n] = x_1
    emissions_sto[n] = emissions

plot_heatmap = False
if plot_heatmap:
    prosumer_old = []
    for i in prosumer:
        if x_0[i] > 0:
            prosumer_old.append(i)
    _results, q_share_old, _community_welfare = FRESH_LP.run_LP(
       load[prosumer_old], PV, df, grid_data, weight, 
       distances, battery, solver_name, sharing=True) 
    PARTICIPATE_plots.heatmap(q_share_old, filename='q_share_old.pdf')