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
x_0 = [1,2,3,5,2,1,0,0,2,2]
d = {}

costs, b, q_share = FRESH_KKT.run_KKT(load, PV, prosumer_data, grid_data, weight, 
              distances, battery, solver_name, years, d, x_0)
