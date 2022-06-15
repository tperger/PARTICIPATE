# -*- coding: utf-8 -*-
"""
Created on Wed May 18 10:47:13 2022

@author: perger
"""

import numpy as np
from numpy import matlib
import pandas as pd
from pyomo.environ import *
  
def run_KKT(load, PV, prosumer_data, grid_data, weight, 
            distances, old, battery, solver_name, 
            years, s, s_1, x_0):
        
    time_steps = load.index.tolist()
    index_time = list(range(len(time_steps)))
    prosumer = load.columns.tolist() 
    weight=weight['weight'].tolist()
    scenarios = list(s.keys()) 
    # scenarios = ['Scenario 1']
    p = 1/len(scenarios)
    
    # Define some parameters and variables
    SoC_max = 'Maximum Storage|Electricity|Energy Storage System'
    SoC_min = 'Minimum Storage|Electricity|Energy Storage System'
    q_bat_max = 'Maximum Charge|Electricity|Energy Storage System'
    q_bat_min = 'Maximum Discharge|Electricity|Energy Storage System'
    w_carb = 'Price|Carbon'
    
    # deactivate BESS
    if battery == False:
        for i in prosumer:
            prosumer_data.loc[SoC_max,i] = 0
            prosumer_data.loc[q_bat_max,i] = 0   
    eta_battery = 0.9
    SoC_init = 0
       
    N = len(years)     
    b_0 = {}
    for i in prosumer:
        if x_0[i] > 0:
          b_0[i] = 1
        else:
          b_0[i] = 0
          
    wtp={}
    for i in prosumer:
        wtp[i] = {}
        for j in prosumer:
            wtp[i][j] = (
                grid_data['Residential']
                + prosumer_data.loc[w_carb,j] 
                * (1 - distances.loc[i,j])
                * grid_data['Emissions'] 
                / 1000000)
        
    # Define model as concrete model
    model = ConcreteModel()
    
    # Define decision variables: upper level
    model.x_1 = Var(prosumer,
                    within=NonNegativeReals) # state variable

    model.u_1 = Var(prosumer,
                    within=NonNegativeReals) # control variable
    
    model.x = Var(years[1:],
                  prosumer,
                  scenarios,
                  within=NonNegativeReals) # state variable

    model.u = Var(years[1:],
                  prosumer,
                  scenarios,
                  within=NonNegativeReals) # control variable

    model.b = Var(years,
                  prosumer,
                  scenarios,
                  within=Binary) # binary auxilary variable    
    
    # Define decision variables: lower level
    model.q_G_in = Var(time_steps, 
                       prosumer, 
                       years,
                       scenarios,
                       within = NonNegativeReals)
    model.q_G_out = Var(time_steps, 
                        prosumer, 
                        years,
                        scenarios,
                        within = NonNegativeReals)
    model.q_share = Var(time_steps, 
                        prosumer, 
                        prosumer, 
                        years,
                        scenarios,
                        within = NonNegativeReals)
    model.q_B_in = Var(time_steps, 
                       prosumer, 
                       years,
                       scenarios,
                       within = NonNegativeReals)
    model.q_B_out = Var(time_steps, 
                        prosumer, 
                        years,
                        scenarios,
                        within = NonNegativeReals)
    model.SoC = Var(time_steps, 
                    prosumer, 
                    years,
                    scenarios,
                    within = NonNegativeReals)
    
    # Define dual variables
    
    model.lambda_load = Var(time_steps, prosumer, years, scenarios)
    model.lambda_PV = Var(time_steps, prosumer, years, scenarios)
    model.lambda_SoC = Var(time_steps, prosumer, years, scenarios)
    model.lambda_SoC_init = Var(prosumer, years, scenarios)
    
    model.mu_SoC = Var(time_steps, prosumer, years, scenarios,
                       within = NonNegativeReals)
    model.mu_B_in = Var(time_steps, prosumer, years, scenarios,
                        within = NonNegativeReals)
    model.mu_B_out = Var(time_steps, prosumer, years, scenarios,
                         within = NonNegativeReals)
    
    # Define auxilary variables for big-M
    
    model.u_G_in = Var(time_steps, prosumer, years, scenarios,
                       within=Binary)
    model.u_G_out = Var(time_steps, prosumer, years, scenarios,
                        within=Binary)
    model.u_share = Var(time_steps, prosumer, prosumer, years, scenarios, 
                        within=Binary)
    model.u_B_in = Var(time_steps, prosumer, years, scenarios,
                       within=Binary)
    model.u_B_out = Var(time_steps, prosumer, years, scenarios,
                        within=Binary)
    model.u_SoC = Var(time_steps, prosumer, years, scenarios,
                      within=Binary)
    model.u_SoC_max = Var(time_steps, prosumer, years, scenarios,
                          within=Binary)
    model.u_B_max_in = Var(time_steps, prosumer, years, scenarios,
                           within=Binary)
    model.u_B_max_out = Var(time_steps, prosumer, years, scenarios,
                            within=Binary)
    
    
    model.M1 = 5000000
    model.M2 = 2000000
    model.M3 = 200
    
    print('variable declaration done')
    
    # Upper level constraints

    # transition function
    
    def transition_fct_rule(model, i, n, w):
        if n == 1:
            if s_1[i] == 1:
                return (x_0[i] 
                        + s_1[i] * model.u_1[i] 
                        - b_0[i] == model.x_1[i])
            if s_1[i] == 0:
                return (model.x_1[i] == 0)
        if n ==2: 
            if s[w][n][i] == 1:
                return (model.x_1[i] 
                        + s[w][n][i] * model.u[n,i,w] 
                        - model.b[n-1,i,w] == model.x[n,i,w])
            if s[w][n][i] == 0:
                return (model.x[n,i,w] == 0)
        else:
            if s[w][n][i] == 1:
                return (model.x[n-1,i,w] 
                        + s[w][n][i] * model.u[n,i,w] 
                        - model.b[n-1,i,w] == model.x[n,i,w])
            if s[w][n][i] == 0:
                return (model.x[n,i,w] == 0)
    model.transition_fct_con = Constraint(prosumer, 
                                          years,
                                          scenarios,
                                          rule=transition_fct_rule)

    # maximum length of contracts = length of horizon 
    def max_x_con_rule(model, i, n, w):
        if n == 1:
            return (model.u_1[i] <= N)
        else:
            return (model.u[n,i,w] <= N)
    model.max_x_con = Constraint(prosumer,
                                 years,
                                 scenarios,
                                 rule=max_x_con_rule)

    # auxilary variable constraints (big-M style)
    def bin1_con_rule(model, i, n, w):
        if n == 1:
            return (model.x_1[i] >= 1 - model.M3 * (1 - model.b[n,i,w]))
        else:
            return (model.x[n,i,w] >= 1 - model.M3 * (1 - model.b[n,i,w]))
    model.bin1_con = Constraint(prosumer, 
                                years,
                                scenarios, 
                                rule = bin1_con_rule)
    def bin2_con_rule(model, i, n, w):
        if n == 1:
            return (model.x_1[i] <= model.M3 * model.b[n,i,w])
        else:
            return (model.x[n,i,w] <= model.M3 * model.b[n,i,w])
    model.bin2_con = Constraint(prosumer, 
                                years,
                                scenarios, 
                                rule = bin2_con_rule)

    # Lower level constraints
    
    # From KKT
    def q_G_in_complementarity_rule_1(model, i, t, n, w):
        return (grid_data.loc[t,'Residential'] * weight[t] * p
                + model.lambda_load[t,i,n,w] >= 0)
    model.q_G_in_compl_1 = Constraint(prosumer, time_steps, years, scenarios,
                                      rule = q_G_in_complementarity_rule_1)
    
    def q_G_in_complementarity_rule_2(model, i, t, n, w):
        return (grid_data.loc[t,'Residential'] * weight[t] * p
                + model.lambda_load[t,i,n,w] <= (1-model.u_G_in[t,i,n,w])*model.M1)
    model.q_G_in_compl_2 = Constraint(prosumer, time_steps, years, scenarios,
                                      rule = q_G_in_complementarity_rule_2)
    
    def q_G_in_complementarity_rule_3(model, i, t, n, w):
        return (model.q_G_in[t,i,n,w] <= model.u_G_in[t,i,n,w]*model.M2)
    model.q_G_in_compl_3 = Constraint(prosumer, time_steps, years, scenarios, 
                                      rule = q_G_in_complementarity_rule_3)
    
    def q_G_out_complementarity_rule_1(model, i, t, n, w):
        return (-grid_data.loc[t,'DA'] * weight[t] * p 
                + model.lambda_PV[t,i,n,w] >= 0)
    model.q_G_out_compl_1 = Constraint(prosumer, time_steps, years, scenarios, 
                                       rule = q_G_out_complementarity_rule_1)
    
    def q_G_out_complementarity_rule_2(model, i, t, n, w):
        return (-grid_data.loc[t,'DA'] * weight[t] * p
                + model.lambda_PV[t,i,n,w] <= (1-model.u_G_out[t,i,n,w])*model.M1)
    model.q_G_out_compl_2 = Constraint(prosumer, time_steps, years, scenarios, 
                                       rule = q_G_out_complementarity_rule_2)
    
    def q_G_out_complementarity_rule_3(model, i, t, n, w):
        return (model.q_G_out[t,i,n,w] <= model.u_G_out[t,i,n,w]*model.M2)
    model.q_G_out_compl_3 = Constraint(prosumer, time_steps, years, scenarios, 
                                       rule = q_G_out_complementarity_rule_3)
    
    def q_share_complementarity_rule_1(model, i, j, t, n, w):
        return (-wtp[i][j][t] * weight[t]  * p
                + model.lambda_load[t,j,n,w] 
                + model.lambda_PV[t,i,n,w] >= 0)
    model.q_share_compl_1 = Constraint(prosumer, prosumer, 
                                       time_steps, years, scenarios, 
                                       rule = q_share_complementarity_rule_1)
    
    def q_share_complementarity_rule_2(model, i, j, t, n, w):
        return (-wtp[i][j][t] * weight[t] * p
                + model.lambda_load[t,j,n,w] 
                + model.lambda_PV[t,i,n,w] <= (1-model.q_share[t,i,j,n,w])*model.M1)
    model.q_share_compl_2 = Constraint(prosumer, prosumer, 
                                       time_steps, years, scenarios,
                                       rule = q_share_complementarity_rule_2)
    
    def q_share_complementarity_rule_3(model, i, j, t, n, w):
        return (model.q_share[t,i,j,n,w] <= model.q_share[t,i,j,n,w] * model.M2)
    model.q_share_compl_3 = Constraint(prosumer, prosumer, 
                                       time_steps, years, scenarios, 
                                       rule = q_share_complementarity_rule_3)
    
    def q_B_in_complementarity_rule_1(model, i, t, n, w):
        return (model.lambda_PV[t,i,n,w] 
                + model.lambda_SoC[t,i,n,w] * eta_battery 
                + model.mu_B_in[t,i,n,w] >= 0)
    model.q_B_in_compl_1 = Constraint(prosumer, time_steps, years, scenarios,
                                      rule = q_B_in_complementarity_rule_1)
    
    def q_B_in_complementarity_rule_2(model, i, t, n, w):
        return (model.lambda_PV[t,i,n,w] 
                + model.lambda_SoC[t,i,n,w] * eta_battery 
                + model.mu_B_in[t,i,n,w] <= (1-model.u_B_in[t,i,n,w]) * model.M1)
    model.q_B_in_compl_2 = Constraint(prosumer, time_steps, years, scenarios, 
                                      rule = q_B_in_complementarity_rule_2)
    
    def q_B_in_complementarity_rule_3(model, i, t, n, w):
        return (model.q_B_in[t,i,n,w] <= model.u_B_in[t,i,n,w] * model.M2)
    model.q_B_in_compl_3 = Constraint(prosumer, time_steps, years, scenarios, 
                                      rule = q_B_in_complementarity_rule_3)
    
    def q_B_out_complementarity_rule_1(model, i, t, n, w):
        return (model.lambda_load[t,i,n,w] 
                - model.lambda_SoC[t,i,n,w] / eta_battery 
                + model.mu_B_out[t,i,n,w] >= 0)
    model.q_B_out_compl_1 = Constraint(prosumer, time_steps, years, scenarios, 
                                       rule = q_B_out_complementarity_rule_1)
    
    def q_B_out_complementarity_rule_2(model, i, t, n, w):
        return (model.lambda_load[t,i,n,w] 
                - model.lambda_SoC[t,i,n,w] / eta_battery 
                + model.mu_B_out[t,i,n,w] <= (1-model.u_B_out[t,i,n,w]) * model.M1)
    model.q_B_out_compl_2 = Constraint(prosumer, time_steps, years, scenarios, 
                                       rule = q_B_out_complementarity_rule_2)
    
    def q_B_out_complementarity_rule_3(model, i, t, n, w):
        return (model.q_B_out[t,i,n,w] <= model.u_B_out[t,i,n,w] * model.M2)
    model.q_B_out_compl_3 = Constraint(prosumer, time_steps, years, scenarios, 
                                       rule = q_B_out_complementarity_rule_3)
    
    def SoC_complementarity_rule_1(model, i, t, n, w):
        if t < index_time[-1]:
            return (-model.lambda_SoC[time_steps[t],i,n,w] 
                    + model.lambda_SoC[time_steps[t+1],i,n,w] 
                    + model.mu_SoC[time_steps[t],i,n,w] >= 0)
        elif t == index_time[-1]:
            return (-model.lambda_SoC[time_steps[t],i,n,w] 
                    + model.lambda_SoC[time_steps[0],i,n,w] 
                    + model.mu_SoC[time_steps[t],i,n,w]
                    + model.lambda_SoC_init[i,n,w]  >= 0)
    model.SoC_compl_1 = Constraint(prosumer, index_time, years, scenarios, 
                                   rule = SoC_complementarity_rule_1)
    
    def SoC_complementarity_rule_2(model, i, t, n, w):
        if t < index_time[-1]:
            return (-model.lambda_SoC[time_steps[t],i,n,w] 
                    + model.lambda_SoC[time_steps[t+1],i,n,w] 
                    + model.mu_SoC[time_steps[t],i,n,w] 
                    <= (1-model.u_SoC[time_steps[t],i,n,w]) * model.M1)
        elif t == index_time[-1]:
            return (-model.lambda_SoC[time_steps[t],i,n,w] 
                    + model.lambda_SoC[time_steps[0],i,n,w] 
                    + model.mu_SoC[time_steps[t],i,n,w] 
                    + model.lambda_SoC_init[i,n,w]
                    <= (1-model.u_SoC[time_steps[t],i,n,w]) * model.M1)
    model.SoC_compl_2 = Constraint(prosumer, index_time, years, scenarios, 
                                   rule = SoC_complementarity_rule_2)
    
    def SoC_complementarity_rule_3(model, i, t, n, w):
            return (model.SoC[t,i,n,w] <= model.u_SoC[t,i,n,w] * model.M2)
    model.SoC_compl_3 = Constraint(prosumer, time_steps, years, scenarios, 
                                   rule = SoC_complementarity_rule_3)
    
    # Equality constraints: lambda
    
    def load_constraint_rule(model, i, t, n, w):    
        return (model.q_G_in[t,i,n,w] 
                + model.q_B_out[t,i,n,w] 
                + sum(model.q_share[t,j,i,n,w] for j in prosumer) 
                - model.b[n,i,w]*load.loc[t,i] == 0)
    model.load_con = Constraint(prosumer, time_steps, years, scenarios,
                                rule = load_constraint_rule)
    
    
    def PV_constraint_rule(model, i, t, n, w):    
        return (model.q_G_out[t,i,n,w] 
                + model.q_B_in[t,i,n,w] 
                + sum(model.q_share[t,i,j,n,w] for j in prosumer) 
                - model.b[n,i,w]*PV.loc[t,i] == 0)
    model.PV_con = Constraint(prosumer, time_steps, years, scenarios, 
                              rule = PV_constraint_rule)
    
    def SoC_constraint_rule(model, i, t, n, w):
        if t == 0:
            return (model.SoC[time_steps[-1],i,n,w] 
                    + model.q_B_in[time_steps[t],i,n,w] * eta_battery 
                    - model.q_B_out[time_steps[t],i,n,w] / eta_battery
                    - model.SoC[time_steps[t],i,n,w] == 0)
        elif t > 0:
            return (model.SoC[time_steps[t-1],i,n,w] 
                    + model.q_B_in[time_steps[t],i,n,w] * eta_battery 
                    - model.q_B_out[time_steps[t],i,n,w] / eta_battery
                    - model.SoC[time_steps[t],i,n,w] == 0)
    model.SoC_con = Constraint(prosumer, index_time, years, scenarios, 
                               rule = SoC_constraint_rule)
   
    def SoC_init_constraint_rule(model, i, n, w):
            return (model.SoC[time_steps[-1],i,n,w] == SoC_init)
    model.SoC_init_con = Constraint(prosumer, years, scenarios, 
                                    rule = SoC_init_constraint_rule)
    
    
    # Inequality constraints: mu
    
    def mu_SoC_max_complementarity_rule_1(model, i, t, n, w):
        return (model.b[n,i,w] * prosumer_data.loc[SoC_max][i] 
                - model.SoC[t,i,n,w] >= 0)
    model.mu_SoC_max_compl_1 = Constraint(prosumer, time_steps, 
                                          years, scenarios,
                                          rule = mu_SoC_max_complementarity_rule_1)
    
    def mu_SoC_max_complementarity_rule_2(model, i, t, n, w):
        return (model.b[n,i,w] * prosumer_data.loc[SoC_max][i] 
                - model.SoC[t,i,n,w] <= (1-model.u_SoC_max[t,i,n,w])*model.M1)
    model.mu_SoC_max_compl_2 = Constraint(prosumer, time_steps, 
                                          years, scenarios, 
                                          rule = mu_SoC_max_complementarity_rule_2)
    
    def mu_SoC_max_complementarity_rule_3(model, i, t, n, w):
        return (model.mu_SoC[t,i,n,w] <= model.u_SoC_max[t,i,n,w]*model.M2)
    model.mu_SoC_max_compl_3 = Constraint(prosumer, time_steps, 
                                          years, scenarios, 
                                          rule = mu_SoC_max_complementarity_rule_3)
    
    def mu_q_B_in_complementarity_rule_1(model, i, t, n, w):
        return (model.b[n,i,w] * prosumer_data.loc[q_bat_max,i] 
                - model.q_B_in[t,i,n,w] >= 0)
    model.mu_q_B_in_compl_1 = Constraint(prosumer, time_steps, 
                                         years, scenarios, 
                                         rule = mu_q_B_in_complementarity_rule_1)
    
    def mu_q_B_in_complementarity_rule_2(model, i, t, n, w):
        return (model.b[n,i,w] * prosumer_data.loc[q_bat_max,i] 
                - model.q_B_in[t,i,n,w] <= (1-model.u_B_max_in[t,i,n,w])*model.M1)
    model.mu_q_B_in_compl_2 = Constraint(prosumer, time_steps, 
                                         years, scenarios, 
                                         rule = mu_q_B_in_complementarity_rule_2)
    
    def mu_q_B_in_complementarity_rule_3(model, i, t, n, w):
        return (model.mu_B_in[t,i,n,w] <= model.u_B_max_in[t,i,n,w]*model.M2)
    model.mu_q_B_in_compl_3 = Constraint(prosumer, time_steps, 
                                         years, scenarios, 
                                         rule = mu_q_B_in_complementarity_rule_3)
    
    def mu_q_B_out_complementarity_rule_1(model, i, t, n, w):
        return (model.b[n,i,w] * prosumer_data.loc[q_bat_min,i] 
                - model.q_B_out[t,i,n,w] >= 0)
    model.mu_q_B_out_compl_1 = Constraint(prosumer, time_steps, 
                                          years, scenarios, 
                                          rule = mu_q_B_out_complementarity_rule_1)
    
    def mu_q_B_out_complementarity_rule_2(model, i, t, n, w):
        return (model.b[n,i,w] * prosumer_data.loc[q_bat_min,i] 
                - model.q_B_out[t,i,n,w] <= (1-model.u_B_max_out[t,i,n,w])*model.M1)
    model.mu_q_B_out_compl_2 = Constraint(prosumer, time_steps, 
                                          years, scenarios, 
                                          rule = mu_q_B_out_complementarity_rule_2)
    
    def mu_q_B_out_complementarity_rule_3(model, i, t, n, w):
        return (model.mu_B_out[t,i,n,w] <= model.u_B_max_out[t,i,n,w]*model.M2)
    model.mu_q_B_out_compl_3 = Constraint(prosumer, time_steps, 
                                          years, scenarios, 
                                          rule = mu_q_B_out_complementarity_rule_3)
        
    
    
    
    emissions_1 = {}
    for i in prosumer:
        emissions_1[i] = sum(model.q_G_in[t,i,years[0],scenarios[0]]
                              * weight[t]
                              * grid_data.loc[t,'Emissions']
                              / 1000000 
                              for t in time_steps) 
        
    emissions = {}
    for w in scenarios:
        emissions[w] = {}
        for n in years[1:]:
            emissions[w][n] = {}
            for i in prosumer:
                emissions[w][n][i] = sum(model.q_G_in[t,i,n,w]
                                          * weight[t]
                                          * grid_data.loc[t,'Emissions']
                                          / 1000000 
                                          for t in time_steps)
                
    # costs_1 = {}
    # for i in prosumer:
    #     costs_1[i] = sum((
    #         model.q_G_in[t,i,years[0],scenarios[0]] * grid_data.loc[t,'Residential']
    #         - model.q_G_out[t,i,years[0],scenarios[0]] * grid_data.loc[t,'DA']
    #         + sum(model.q_share[t,j,i,years[0],scenarios[0]] * wtp[j][i][t]
    #         - model.q_share[t,i,j,years[0],scenarios[0]] * wtp[i][j][t] for j in prosumer))
    #         * weight[t]
    #         for t in time_steps)    
    
    # costs = {}
    # for w in scenarios:
    #     costs[w] = {}
    #     for n in years[1:]:
    #         costs[w][n] = {}
    #         for i in prosumer:
    #             costs[w][n][i] = sum((
    #                 model.q_G_in[t,i,n,w] * grid_data.loc[t,'Residential']
    #                 - model.q_G_out[t,i,n,w] * grid_data.loc[t,'DA']
    #                 + sum(model.q_share[t,j,i,n,w] * wtp[j][i][t]
    #                 - model.q_share[t,i,j,n,w] * wtp[i][j][t] for j in prosumer))
    #                 * weight[t]
    #                 for t in time_steps)
    
   
    # objective functions
 

    # F1 = (sum((costs_1[i] - old[i]) * s_1[i] * b_0[i] 
    #           for i in prosumer) 
    #       + p * sum((costs[w][n][i] - model.b[n,i,w] * old[i]) * s[w][n][i] * b_0[i]
    #                 for i in prosumer
    #                 for n in years[1:]
    #                 for w in scenarios))
    
    F2 = (sum((emissions_1[i] - old[i]) * s_1[i] * b_0[i] 
              for i in prosumer) 
          + p * sum((emissions[w][n][i] - model.b[n,i,w] * old[i]) * s[w][n][i] * b_0[i]
                    for i in prosumer
                    for n in years[1:]
                    for w in scenarios))
    

    
    # choose one of the objective functions F1, F2, ... defined above
    
    print('Solving...')
    model.obj = Objective(expr = F2, 
                          sense = minimize)
    
    opt = SolverFactory(solver_name)
    opt_success = opt.solve(model)
    
    # Evaluate the results
    # social_welfare = value(sum(community_welfare[i] 
    #                            + prosumer_welfare[i] for i in prosumer))
    
    q_share = pd.DataFrame(index=prosumer)
    n = years[-1]
    for j in prosumer:
        a = []
        for i in prosumer:
            a.append(value(sum(model.q_share[t,i,j,n,w]*weight[t]  
                               for t in time_steps
                               for w in scenarios)))
        q_share[j] = a
    
    # results = pd.DataFrame(index=prosumer)
    # for i in prosumer:
    #     results.loc[i,'buying grid'] = value(sum(model.q_G_in[t,i]*weight[t]  
    #                                              for t in time_steps))
    #     results.loc[i,'selling grid'] = value(sum(model.q_G_out[t,i]*weight[t]  
    #                                               for t in time_steps))
    #     results.loc[i,'battery charging'] = value(sum(model.q_B_in[t,i]*weight[t]  
    #                                                   for t in time_steps))
    #     results.loc[i,'battery discharging'] = value(sum(model.q_B_out[t,i]*weight[t]  
    #                                                      for t in time_steps))
    #     results.loc[i,'self-consumption'] = q_share.loc[i,i]
    #     results.loc[i,'buying community'] = (sum(q_share.loc[j,i] 
    #                                              for j in prosumer) 
    #                                          - q_share.loc[i,i])
    #     results.loc[i,'selling community'] = (sum(q_share.loc[i,j] 
    #                                               for j in prosumer) 
    #                                           - q_share.loc[i,i])
    #     results.loc[i,'emissions'] = (value(sum(model.q_G_in[t,i]*weight[t] 
    #                                             *grid_data.Emissions.loc[t]
    #                                             / 1000000
    #                                             for t in time_steps)))
    #     results.loc[i,'costs'] = (value(-community_welfare[i]) 
    #                               - value(prosumer_welfare[i]) 
    #                               + value(prosumer_welfare2[i])) 
    
    # parameter = pd.DataFrame(index=prosumer_new)
    # for i in prosumer_new:
    #     parameter.loc[i, 'PV'] = value(model.PV_new[i])
    #     parameter.loc[i, 'load'] = value(model.load_new[i])
        
    b = {}
    for w in scenarios:
        b[w] = {}
        for n in years:
            b[w][n] = {}
            for i in prosumer:
                b[w][n][i] = (value(model.b[n,i,w]))
                
    u = {}
    for w in scenarios:
        u[w] = {}
        for n in years[1:]:
            u[w][n] = {}
            for i in prosumer:
                u[w][n][i] = (value(model.u[n,i,w]))
    u_1 = {}
    for i in prosumer:
        u_1[i] = (value(model.u_1[i]))


    return b, u, u_1, q_share