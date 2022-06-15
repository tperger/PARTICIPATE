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

import PARTICIPATE_functions as pf
import PARTICIPATE_KKT_neu
import PARTICIPATE_KKT_deterministic  
import FRESH_LP

solver_name = 'gurobi'
battery = False

country = 'Austria'
scenario_name = 'Default scenario'
model_name = 'FRESH:COM v2.0'

building_types = ['SH', 'SAB', 'LAB']

settlement_pattern = 'rural'

buildings_per_SP ={'city': {'SH': 0, 'SAB': 0, 'LAB': 10}, #city 
                   'town': {'SH': 0, 'SAB': 10, 'LAB': 0}, #town
                   'suburban': {'SH': 10, 'SAB': 0, 'LAB': 2}, #suburban
                   'rural': {'SH': 10, 'SAB': 10, 'LAB': 0} #rural
                   }  

load, PV, prosumer_data, grid_data, weight, distances = pf.define_community(
    settlement_pattern,
    buildings_per_SP,
    model_name,   
    scenario_name,
    country,
    year=2019,
    clustering=True)

N = 3
years = np.arange(1,N+1).tolist()

# Generate scenarios
num_scen = 4
scenarios = ['Scenario '+str(i) for i in range(1,num_scen+1)]

prosumer = load.columns.tolist() 

# x_0 = [1,2,3,5,2,1,0,0,2,2,
#        1,2,3,5,2,1,0,0,2,2]
x_0 = [3,3,3,0,0,2,2,0,0,0,
       0,0,0,2,2,0,0,3,3,3]
x_0 = dict(zip(prosumer, x_0))
prosumer = load.columns.tolist() 
values = [1,1,1,0,0,1,1,0,0,0,
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
                _values = [1,1,1,0,0,1,1,0,0,0,
                          1,1,1,1,1,1,1,1,1,1]
            if w == scenarios[1]:
                _values = [1,1,1,0,0,1,1,0,0,0,
                          0,0,0,0,0,0,0,0,0,0]
            if w == scenarios[2]:
                _values = [1,1,1,1,1,1,1,1,1,1,
                          0,0,0,1,1,0,0,1,1,1]
            if w == scenarios[3]:
                _values = [0,0,0,0,0,0,0,0,0,0,
                          0,0,0,1,1,0,0,1,1,1]
        s[w][n] = dict(zip(prosumer, _values))
           
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

df = prosumer_data.copy(deep=True)

results, q_share_total, community_welfare = FRESH_LP.run_LP(
    load, PV, df, grid_data, weight, 
    distances, battery, solver_name, sharing=False)

emissions_old = results['emissions']
costs_old = results['costs']

# run model without forcast

b_det = {}
u_1_det = {}
x_1_det = {}
emissions_det = {}

for n in years:   
    if n == 1:
        _x_0 = x_0
        _s_1 = s_1
    else:
        _x_0 = x_1
        _s_1 = s['Scenario 1'][n]
    b, u_1, x_1, emissions = PARTICIPATE_KKT_deterministic.run_KKT(
        load, PV, df, grid_data, weight, 
        distances, emissions_old, battery, solver_name, [n], _s_1, _x_0)
    b_det[n] = b['Scenario 1'][n]
    u_1_det[n] = u_1
    x_1_det[n] = x_1
    emissions_det[n] = emissions

b, u, u1, q_share = PARTICIPATE_KKT_neu.run_KKT(
    load, PV, df, grid_data, weight, 
    distances, costs_old, battery, solver_name, years, s, s_1, x_0)
