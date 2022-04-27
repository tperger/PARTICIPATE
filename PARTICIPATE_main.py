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

import FRESH_define_community as cm
import FRESH_KKT
import FRESH_plots

solver_name = 'gurobi'
battery = True

prosumer = (['Prosumer 1',
             'Prosumer 2',
             'Prosumer 3',
             'Prosumer 4', 
             'Prosumer 5', 
             'Prosumer 6'])


if cm.clustering == True:
    hours = len(cm.time_steps) / len(cm.counts)
    weight = np.matlib.repmat(cm.counts,int(hours),1).transpose().reshape(len(cm.time_steps),1)[:,0]

N = 5
years = np.arange(1,N+1)
x_0 = [0,0,3,5,2,1]

costs = FRESH_KKT.run_KKT(prosumer,
                          cm, 
                          weight, 
                          battery, 
                          solver_name, 
                          years, 
                          d, 
                          x_0)
    

