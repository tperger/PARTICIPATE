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
import FRESH_KKT
import PARTICIPATE_KKT_neu 
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
                   'rural': {'SH': 10, 'SAB': 0, 'LAB': 0} #rural
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
num_scen = 5
scenarios = ['Scenario '+str(i) for i in range(1,num_scen+1)]

# x_0 = [1,2,3,5,2,1,0,0,2,2,
#        1,2,3,5,2,1,0,0,2,2]
x_0 = [1,2,3,5,2,1,0,0,2,2]
prosumer = load.columns.tolist() 
s_1 = {}
for i in prosumer:
    s_1[i] = np.random.choice([0,1])
s = {}
for w in scenarios:
    s[w] = {}
    for n in years[1:]:
       s[w][n] = {}
       for i in prosumer:
           s[w][n][i] = np.random.choice([0,1])
           
# PV_peak = [5,5,0,3,5,5,0,3,5,0,
#            8,8,0,5,8,8,0,5,8,0]  
# load *= [3335.594403,4537.52957,5252.549817,5823.804207,6336.761219,
#         6833.249609,7346.206622,7917.461011,8632.481259,9834.416426,
#         6257.898815,8512.845835,9854.293211,10926.02189,11888.37902,
#         12819.83943,13782.19656,14853.92524,16195.37261,18450.31963
#         ]
PV_peak = [5,5,0,3,5,5,0,3,5,0]  
load *= [3335.594403,4537.52957,5252.549817,5823.804207,6336.761219,
         6833.249609,7346.206622,7917.461011,8632.481259,9834.416426] 
prosumer_data.loc['Maximum Storage|Electricity|Energy Storage System'] *= 3
prosumer_data.loc['Maximum Charge|Electricity|Energy Storage System'] *= 2
prosumer_data.loc['Maximum Discharge|Electricity|Energy Storage System'] *= 2 

load /= 1000
PV *= PV_peak      

df = prosumer_data.copy(deep=True)

results, q_share_total, community_welfare = FRESH_LP.run_LP(
    load, PV, df, grid_data, weight, 
    distances, battery, solver_name, sharing=False)

emissions = results['emissions']

# costs, b, u, q_share = PARTICIPATE_KKT.run_KKT(load, PV, df, grid_data, weight, 
#               distances, emissions, battery, solver_name, years, d, x_0)
