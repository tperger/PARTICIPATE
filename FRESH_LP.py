# -*- coding: utf-8 -*-
"""
Created on Wed Mar 10 11:51:30 2021

@author: perger
"""

import numpy as np
import pandas as pd
from pyomo.environ import *

def run_LP(load, PV, prosumer_data, grid_data, weight, 
            distances, battery, solver_name, sharing=False):

    time_steps = load.index.tolist()
    index_time = list(range(len(time_steps)))
    prosumer = load.columns.tolist() 
    weight=weight['weight'].tolist()
    
    # Define some parameters and variables
    SoC_max = 'Maximum Storage|Electricity|Energy Storage System'
    SoC_min = 'Minimum Storage|Electricity|Energy Storage System'
    q_bat_max = 'Maximum Charge|Electricity|Energy Storage System'
    q_bat_min = 'Maximum Discharge|Electricity|Energy Storage System'
    w = 'Price|Carbon'


    # deactivate BESS
    # if battery == False:
    #     for i in prosumer:
    #         prosumer_data.loc[SoC_max,i] = 0
    #         prosumer_data.loc[q_bat_max,i] = 0            
    
    eta_battery = 0.9
    
    #index_time = list(range(len(time_steps)))
    
    wtp={}
    for i in prosumer:
        wtp[i] = {}
        for j in prosumer:
            wtp[i][j] = (
                grid_data['Residential']
                + prosumer_data.loc[w,j] 
                * (1 - distances.loc[i,j])
                * grid_data['Emissions'] 
                / 1000000)
    
    # Define model as concrete model
    model = ConcreteModel()
    
    #Define optimization variables 
    model.q_G_in = Var(time_steps, 
                       prosumer, 
                       within = NonNegativeReals)
    model.q_G_out = Var(time_steps, 
                        prosumer, 
                        within = NonNegativeReals)
    model.q_share = Var(time_steps, 
                        prosumer, 
                        prosumer, 
                        within = NonNegativeReals)
    model.q_B_in = Var(time_steps, 
                       prosumer, 
                       within = NonNegativeReals)
    model.q_B_out = Var(time_steps, 
                        prosumer, 
                        within = NonNegativeReals)
    model.SoC = Var(time_steps, 
                    prosumer, 
                    within = NonNegativeReals)
    model.community_welfare = Var(time_steps,
                                  prosumer,
                                  within = Reals)
    
    # Define constraints
     
    def no_sharing_rule(model, i, j, t):
         if i is not j:
             return(model.q_share[t,i,j] == 0)
         else:
             return (model.q_share[t,i,j] >= 0)
     
    if sharing == False:
         model.no_sharing_con = Constraint(prosumer,
                                           prosumer,
                                           time_steps,
                                           rule = no_sharing_rule)
         
    def load_constraint_rule(model, i, t):    
         return (model.q_G_in[t,i] 
                 + model.q_B_out[t,i] 
                 + sum(model.q_share[t,j,i] for j in prosumer)
                 - load.loc[t,i] == 0)
    model.load_con = Constraint(prosumer, 
                                 time_steps, 
                                 rule = load_constraint_rule)
     
    def PV_constraint_rule(model, i, t):    
         return (model.q_G_out[t,i] 
                 + model.q_B_in[t,i] 
                 + sum(model.q_share[t,i,j] for j in prosumer) 
                 - PV.loc[t,i] == 0)
    model.PV_con = Constraint(prosumer, 
                               time_steps, 
                               rule = PV_constraint_rule)
     
    def SoC_min_constraint_rule(model, i, t):
         return (model.SoC[t,i] >= prosumer_data.loc[SoC_min][i])
    model.SoC_min_con = Constraint(prosumer, 
                                    time_steps, 
                                    rule = SoC_min_constraint_rule)
     
    def SoC_max_constraint_rule(model, i, t):
         return (model.SoC[t,i] <= prosumer_data.loc[SoC_max][i])
    model.SoC_max_con = Constraint(prosumer, 
                                    time_steps, 
                                    rule = SoC_max_constraint_rule)
     
    def q_B_in_constraint_rule(model, i, t):
         return (model.q_B_in[t,i] <= prosumer_data.loc[q_bat_max][i])
    model.q_B_in_con = Constraint(prosumer, 
                                  time_steps, 
                                  rule = q_B_in_constraint_rule)
     
    def q_B_out_constraint_rule(model, i, t):
         return (model.q_B_out[t,i] <= prosumer_data.loc[q_bat_min][i])
    model.q_B_out_con = Constraint(prosumer, 
                                   time_steps, 
                                   rule = q_B_out_constraint_rule)
     
    def SoC_constraint_rule(model, i, t):
         if t == 0:
             return (model.SoC[time_steps[-1],i] 
                     + model.q_B_in[time_steps[t],i] * eta_battery 
                     - model.q_B_out[time_steps[t],i] / eta_battery
                     - model.SoC[time_steps[t],i] == 0)
         elif t > 0:
             return (model.SoC[time_steps[t-1],i] 
                     + model.q_B_in[time_steps[t],i] * eta_battery 
                     - model.q_B_out[time_steps[t],i] / eta_battery
                     - model.SoC[time_steps[t],i] == 0)
    model.SoC_con = Constraint(prosumer, 
                               index_time, 
                               rule = SoC_constraint_rule)
     
    def CW_constraint_rule(model, i, t):
         return model.community_welfare[t,i] == (
             - grid_data.loc[t,'Residential'] * model.q_G_in[t,i] * weight[t]
             + grid_data.loc[t, 'DA'] * model.q_G_out[t,i] * weight[t]
             + sum(wtp[j][i][t]
                   * model.q_share[t,j,i]
                   * weight[t]
                   for j in prosumer))
    model.CW_con = Constraint(prosumer, 
                              time_steps, 
                              rule = CW_constraint_rule)
    
    # Objective function
    community_welfare = {new_list: [] for new_list in prosumer}
    prosumer_welfare = {new_list: [] for new_list in prosumer}
    prosumer_welfare2 = {new_list: [] for new_list in prosumer}
       
    
    for i in prosumer:
        community_welfare[i] = sum(- grid_data.loc[t,'Residential'] * model.q_G_in[t,i]*weight[t]
                                    + grid_data.loc[t, 'DA'] * model.q_G_out[t,i]*weight[t] 
                                    for t in time_steps)
        prosumer_welfare[i] = sum(wtp[i][j][t]
                                  * model.q_share[t,i,j]
                                  * weight[t]
                                  for j in prosumer 
                                  for t in time_steps)
        prosumer_welfare2[i] = sum(wtp[j][i][t]
                                   * model.q_share[t,j,i]
                                   * weight[t]
                                   for j in prosumer 
                                   for t in time_steps)
    
    #     # 1. prosumer i sells to prosumer j
    #     # 2. prosumer i buys from prosumer j
    
    # model.obj = Objective(
    #     expr = sum(model.community_welfare[t,i] 
    #                 for t in time_steps 
    #                 for i in prosumer), 
    #     sense = maximize)
    
    model.obj = Objective(
        expr = sum(community_welfare[i] + prosumer_welfare2[i]
                   for i in prosumer), 
        sense = maximize)
    
    opt = SolverFactory(solver_name)
    opt_success = opt.solve(model)
    
    # Evaluate the results
    community_welfare = value(model.obj)
    
    q_share_total = pd.DataFrame(index=prosumer)
    for j in prosumer:
        a = []
        for i in prosumer:
            a.append(value(sum(model.q_share[t,i,j]*weight[t] 
                               for t in time_steps)))
        q_share_total[j] = a
    
    results= pd.DataFrame(index=prosumer)
    for i in prosumer:
        results.loc[i,'buying grid'] = value(sum(model.q_G_in[t,i]
                                                 * weight[t] 
                                                 for t in time_steps))
        results.loc[i,'selling grid'] = value(sum(model.q_G_out[t,i]
                                                  * weight[t] 
                                                  for t in time_steps))
        results.loc[i,'battery charging'] = value(sum(model.q_B_in[t,i]
                                                      * weight[t] 
                                                      for t in time_steps))
        results.loc[i,'battery discharging'] = value(sum(model.q_B_out[t,i]
                                                         * weight[t] 
                                                         for t in time_steps))
        results.loc[i,'self-consumption'] = q_share_total.loc[i,i]
        results.loc[i,'buying community'] = (sum(q_share_total.loc[j,i] 
                                                 for j in prosumer) 
                                             - q_share_total.loc[i,i])
        results.loc[i,'selling community'] = (sum(q_share_total.loc[i,j] 
                                                  for j in prosumer) 
                                              - q_share_total.loc[i,i])
        results.loc[i,'emissions'] = value(sum(model.q_G_in[t,i]
                                               * weight[t]
                                               * grid_data.loc[t,'Emissions']
                                               / 1000000 
                                                for t in time_steps))
        results.loc[i,'costs'] = value(sum(
            grid_data.loc[t,'Residential'] * model.q_G_in[t,i] * weight[t] 
            + sum(wtp[j][i][t] * model.q_share[t,j,i] * weight[t] for j in prosumer)
            - grid_data.loc[t,'DA'] * model.q_G_out[t,i] * weight[t] 
            - sum(wtp[i][j][t] * model.q_share[t,i,j] * weight[t] for j in prosumer)
            for t in time_steps))
        
    
    return results, q_share_total, community_welfare