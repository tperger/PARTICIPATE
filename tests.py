# -*- coding: utf-8 -*-
"""
Created on Mon Apr 25 09:11:45 2022

@author: perger
"""

import numpy as np
from numpy import matlib
import pandas as pd
from pyomo.environ import *

solver_name = 'gurobi'

N = 5
years = np.arange(1,N+1)

prosumer = (['Prosumer 1',
             'Prosumer 2',
             'Prosumer 3',
             'Prosumer 4', 
             'Prosumer 5', 
             'Prosumer 6'])

M = 10 # must be > n

x_0 = [0,0,3,5,2,1]
x_0 = dict(zip(prosumer, x_0))
e_0 = {}
for i in prosumer:
    if x_0[i] > 0:
      e_0[i] = 0
    else:
      e_0[i] = 1


# Generate scenarios

num_scen = 5
scenarios = ['Scenario '+str(i) for i in range(1,num_scen+1)]

d = {}
for w in scenarios:
    d[w] = {}
    for n in years:
        d[w][n] = {}
        for i in prosumer:
            d[w][n][i] = 1

# Define model as concrete model
model = ConcreteModel()

model.x = Var(years,
              prosumer,
              within=NonNegativeIntegers) # state variable

model.u = Var(years,
              prosumer,
              within=NonNegativeIntegers) # control variable

model.b = Var(years,
              prosumer,
              within=Binary) # binary auxilary variable

model.e = Var(years,
              prosumer,
              within=Binary) # binary auxilary variable

# objective function



model.obj = Objective(expr = sum(model.x[n,i] for n in years for i in prosumer), 
                      sense = minimize)

# transition function

def transition_fct_rule(model, i, n):
    if n == 1:
        return (model.x[n,i] == x_0[i] - 1 + model.u[n,i] +e_0[i])
    else:
        return (model.x[n,i] == model.x[n-1,i] - 1 + model.u[n,i] + model.e[n,i])
model.transition_fct_con = Constraint(prosumer, 
                                      years,
                                      rule=transition_fct_rule)

# constraints

# maximum length of contracts = length of horizon 
def max_x_con_rule(model, i, n):
    return (model.x[n,i] <= N)
model.max_x_con = Constraint(prosumer,
                             years,
                             rule=max_x_con_rule)

# auxilary variable constraints (big-M style)
def bin1_con_rule(model, i, n):
    return (model.x[n,i] >= 1 - M * (1 - model.b[n,i]))
model.bin1_con = Constraint(prosumer, 
                            years, 
                            rule = bin1_con_rule)
def bin2_con_rule(model, i, n):
    return (model.x[n,i] <= M * model.b[n,i])
model.bin2_con = Constraint(prosumer, 
                            years, 
                            rule = bin2_con_rule)
def bin3_con_rule(model, i, n):
    return model.b[n,i] + model.e[n,i] == 1
model.bin3_con = Constraint(prosumer,
                            years,
                            rule=bin3_con_rule)

opt = SolverFactory(solver_name)
opt_success = opt.solve(model)