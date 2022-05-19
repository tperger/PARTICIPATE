# -*- coding: utf-8 -*-
"""
Created on Wed Mar  3 14:54:50 2021

@author: perger
"""

import numpy as np
from numpy import matlib
import pandas as pd
from pyomo.environ import *
  
def run_KKT(load, PV, prosumer_data, grid_data, weight, 
            distances, emissions_wo_comm, battery, solver_name, years, d, x_0):
    
    
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
    if battery == False:
        for i in prosumer:
            prosumer_data.loc[SoC_max,i] = 0
            prosumer_data.loc[q_bat_max,i] = 0   
    eta_battery = 0.9
    SoC_init = 0
       
    N = len(years)     
    x_0 = dict(zip(prosumer, x_0))
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
                + prosumer_data.loc[w,j] 
                * (1 - distances.loc[i,j])
                * grid_data['Emissions'] 
                / 1000000)
        
    # Define model as concrete model
    model = ConcreteModel()
    
    # Define decision variables: upper level
    model.x = Var(years,
                  prosumer,
                  within=NonNegativeReals) # state variable

    model.u = Var(years,
                  prosumer,
                  within=NonNegativeReals) # control variable

    model.b = Var(years,
                  prosumer,
                  within=Binary) # binary auxilary variable    
    
    # Define decision variables: lower level
    model.q_G_in = Var(time_steps, 
                       prosumer, 
                       years,
                       within = NonNegativeReals)
    model.q_G_out = Var(time_steps, 
                        prosumer, 
                        years,
                        within = NonNegativeReals)
    model.q_share = Var(time_steps, 
                        prosumer, 
                        prosumer, 
                        years,
                        within = NonNegativeReals)
    model.q_B_in = Var(time_steps, 
                       prosumer, 
                       years,
                       within = NonNegativeReals)
    model.q_B_out = Var(time_steps, 
                        prosumer, 
                        years,
                        within = NonNegativeReals)
    model.SoC = Var(time_steps, 
                    prosumer, 
                    years,
                    within = NonNegativeReals)
    
    # Define dual variables
    
    model.lambda_load = Var(time_steps, prosumer, years)
    model.lambda_PV = Var(time_steps, prosumer, years)
    model.lambda_SoC = Var(time_steps, prosumer, years)
    model.lambda_SoC_init = Var(prosumer, years)
    
    model.mu_SoC = Var(time_steps, prosumer, years,
                       within = NonNegativeReals)
    model.mu_B_in = Var(time_steps, prosumer, years, 
                        within = NonNegativeReals)
    model.mu_B_out = Var(time_steps, prosumer, years, 
                         within = NonNegativeReals)
    
    # Define auxilary variables for big-M
    
    model.u_G_in = Var(time_steps, prosumer, years, 
                       within=Binary)
    model.u_G_out = Var(time_steps, prosumer, years, 
                        within=Binary)
    model.u_share = Var(time_steps, prosumer, prosumer, years, 
                        within=Binary)
    model.u_B_in = Var(time_steps, prosumer, years, 
                       within=Binary)
    model.u_B_out = Var(time_steps, prosumer, years, 
                        within=Binary)
    model.u_SoC = Var(time_steps, prosumer, years, 
                      within=Binary)
    model.u_SoC_max = Var(time_steps, prosumer, years, 
                          within=Binary)
    model.u_B_max_in = Var(time_steps, prosumer, years, 
                           within=Binary)
    model.u_B_max_out = Var(time_steps, prosumer, years, 
                            within=Binary)
    
    
    model.M1 = 5000000
    model.M2 = 2000000
    model.M3 = 200
    
    print('variable declaration done')
    
    # Upper level constraints

    # transition function
    
    def transition_fct_rule(model, i, n):
        if d[n][i] == 1:
            if n == 1:
                return (model.x[n,i] == x_0[i] + d[n][i] * model.u[n,i] - b_0[i])
            else:
                return (model.x[n,i] == model.x[n-1,i] + model.u[n,i] - model.b[n,i])
        if d[n][i] == 0:
            return (model.x[n,i] == 0)
    model.transition_fct_con = Constraint(prosumer, 
                                          years,
                                          rule=transition_fct_rule)

    # maximum length of contracts = length of horizon 
    def max_x_con_rule(model, i, n):
        return (model.u[n,i] <= N)
    model.max_x_con = Constraint(prosumer,
                                  years,
                                  rule=max_x_con_rule)

    # auxilary variable constraints (big-M style)
    def bin1_con_rule(model, i, n):
        return (model.x[n,i] >= 1 - model.M3 * (1 - model.b[n,i]))
    model.bin1_con = Constraint(prosumer, 
                                years, 
                                rule = bin1_con_rule)
    def bin2_con_rule(model, i, n):
        return (model.x[n,i] <= model.M3 * model.b[n,i])
    model.bin2_con = Constraint(prosumer, 
                                years, 
                                rule = bin2_con_rule)

    # Lower level constraints
    
    # From KKT
    def q_G_in_complementarity_rule_1(model, i, t, n):
        return (grid_data.loc[t,'Residential'] * weight[t] 
                + model.lambda_load[t,i,n] >= 0)
    model.q_G_in_compl_1 = Constraint(prosumer, time_steps, years,
                                      rule = q_G_in_complementarity_rule_1)
    
    def q_G_in_complementarity_rule_2(model, i, t, n):
        return (grid_data.loc[t,'Residential'] * weight[t]
                + model.lambda_load[t,i,n] <= (1-model.u_G_in[t,i,n])*model.M1)
    model.q_G_in_compl_2 = Constraint(prosumer, time_steps, years,
                                      rule = q_G_in_complementarity_rule_2)
    
    def q_G_in_complementarity_rule_3(model, i, t, n):
        return (model.q_G_in[t,i,n] <= model.u_G_in[t,i,n]*model.M2)
    model.q_G_in_compl_3 = Constraint(prosumer, time_steps, years, 
                                      rule = q_G_in_complementarity_rule_3)
    
    def q_G_out_complementarity_rule_1(model, i, t, n):
        return (-grid_data.loc[t,'DA'] * weight[t] 
                + model.lambda_PV[t,i,n] >= 0)
    model.q_G_out_compl_1 = Constraint(prosumer, time_steps, years, 
                                       rule = q_G_out_complementarity_rule_1)
    
    def q_G_out_complementarity_rule_2(model, i, t, n):
        return (-grid_data.loc[t,'DA'] * weight[t]
                + model.lambda_PV[t,i,n] <= (1-model.u_G_out[t,i,n])*model.M1)
    model.q_G_out_compl_2 = Constraint(prosumer, time_steps, years, 
                                       rule = q_G_out_complementarity_rule_2)
    
    def q_G_out_complementarity_rule_3(model, i, t, n):
        return (model.q_G_out[t,i,n] <= model.u_G_out[t,i,n]*model.M2)
    model.q_G_out_compl_3 = Constraint(prosumer, time_steps, years, 
                                       rule = q_G_out_complementarity_rule_3)
    
    def q_share_complementarity_rule_1(model, i, j, t, n):
        return (-wtp[i][j][t] * weight[t] 
                + model.lambda_load[t,j,n] 
                + model.lambda_PV[t,i,n] >= 0)
    model.q_share_compl_1 = Constraint(prosumer, prosumer, time_steps, years, 
                                       rule = q_share_complementarity_rule_1)
    
    def q_share_complementarity_rule_2(model, i, j, t, n):
        return (-wtp[i][j][t] * weight[t]
                + model.lambda_load[t,j,n] 
                + model.lambda_PV[t,i,n] <= (1-model.u_share[t,i,j,n])*model.M1)
    model.q_share_compl_2 = Constraint(prosumer, prosumer, time_steps, years, 
                                       rule = q_share_complementarity_rule_2)
    
    def q_share_complementarity_rule_3(model, i, j, t, n):
        return (model.q_share[t,i,j,n] <= model.u_share[t,i,j,n] * model.M2)
    model.q_share_compl_3 = Constraint(prosumer, prosumer, time_steps, years, 
                                       rule = q_share_complementarity_rule_3)
    
    def q_B_in_complementarity_rule_1(model, i, t, n):
        return (model.lambda_PV[t,i,n] 
                + model.lambda_SoC[t,i,n] * eta_battery 
                + model.mu_B_in[t,i,n] >= 0)
    model.q_B_in_compl_1 = Constraint(prosumer, time_steps, years,
                                      rule = q_B_in_complementarity_rule_1)
    
    def q_B_in_complementarity_rule_2(model, i, t, n):
        return (model.lambda_PV[t,i,n] 
                + model.lambda_SoC[t,i,n] * eta_battery 
                + model.mu_B_in[t,i,n] <= (1-model.u_B_in[t,i,n]) * model.M1)
    model.q_B_in_compl_2 = Constraint(prosumer, time_steps, years, 
                                      rule = q_B_in_complementarity_rule_2)
    
    def q_B_in_complementarity_rule_3(model, i, t, n):
        return (model.q_B_in[t,i,n] <= model.u_B_in[t,i,n] * model.M2)
    model.q_B_in_compl_3 = Constraint(prosumer, time_steps, years, 
                                      rule = q_B_in_complementarity_rule_3)
    
    def q_B_out_complementarity_rule_1(model, i, t, n):
        return (model.lambda_load[t,i,n] 
                - model.lambda_SoC[t,i,n] / eta_battery 
                + model.mu_B_out[t,i,n] >= 0)
    model.q_B_out_compl_1 = Constraint(prosumer, time_steps, years, 
                                       rule = q_B_out_complementarity_rule_1)
    
    def q_B_out_complementarity_rule_2(model, i, t, n):
        return (model.lambda_load[t,i,n] 
                - model.lambda_SoC[t,i,n] / eta_battery 
                + model.mu_B_out[t,i,n] <= (1-model.u_B_out[t,i,n]) * model.M1)
    model.q_B_out_compl_2 = Constraint(prosumer, time_steps, years, 
                                       rule = q_B_out_complementarity_rule_2)
    
    def q_B_out_complementarity_rule_3(model, i, t, n):
        return (model.q_B_out[t,i,n] <= model.u_B_out[t,i,n] * model.M2)
    model.q_B_out_compl_3 = Constraint(prosumer, time_steps, years, 
                                       rule = q_B_out_complementarity_rule_3)
    
    def SoC_complementarity_rule_1(model, i, t, n):
        if t < index_time[-1]:
            return (-model.lambda_SoC[time_steps[t],i,n] 
                    + model.lambda_SoC[time_steps[t+1],i,n] 
                    + model.mu_SoC[time_steps[t],i,n] >= 0)
        elif t == index_time[-1]:
            return (-model.lambda_SoC[time_steps[t],i,n] 
                    + model.lambda_SoC[time_steps[0],i,n] 
                    + model.mu_SoC[time_steps[t],i,n]
                    + model.lambda_SoC_init[i,n]  >= 0)
    model.SoC_compl_1 = Constraint(prosumer, index_time, years, 
                                   rule = SoC_complementarity_rule_1)
    
    def SoC_complementarity_rule_2(model, i, t, n):
        if t < index_time[-1]:
            return (-model.lambda_SoC[time_steps[t],i,n] 
                    + model.lambda_SoC[time_steps[t+1],i,n] 
                    + model.mu_SoC[time_steps[t],i,n] 
                    <= (1-model.u_SoC[time_steps[t],i,n]) * model.M1)
        elif t == index_time[-1]:
            return (-model.lambda_SoC[time_steps[t],i,n] 
                    + model.lambda_SoC[time_steps[0],i,n] 
                    + model.mu_SoC[time_steps[t],i,n] 
                    + model.lambda_SoC_init[i,n]
                    <= (1-model.u_SoC[time_steps[t],i,n]) * model.M1)
    model.SoC_compl_2 = Constraint(prosumer, index_time, years, 
                                   rule = SoC_complementarity_rule_2)
    
    def SoC_complementarity_rule_3(model, i, t, n):
            return (model.SoC[t,i,n] <= model.u_SoC[t,i,n] * model.M2)
    model.SoC_compl_3 = Constraint(prosumer, time_steps, years, 
                                   rule = SoC_complementarity_rule_3)
    
    # Equality constraints: lambda
    
    def load_constraint_rule(model, i, t, n):    
        return (model.q_G_in[t,i,n] 
                + model.q_B_out[t,i,n] 
                + sum(model.q_share[t,j,i,n] for j in prosumer) 
                - model.b[n,i]*load.loc[t,i] == 0)
    model.load_con = Constraint(prosumer, time_steps, years,
                                rule = load_constraint_rule)
    
    
    def PV_constraint_rule(model, i, t, n):    
        return (model.q_G_out[t,i,n] 
                + model.q_B_in[t,i,n] 
                + sum(model.q_share[t,i,j,n] for j in prosumer) 
                - model.b[n,i]*PV.loc[t,i] == 0)
    model.PV_con = Constraint(prosumer, time_steps, years, 
                              rule = PV_constraint_rule)
    
    def SoC_constraint_rule(model, i, t, n):
        if t == 0:
            return (model.SoC[time_steps[-1],i,n] 
                    + model.q_B_in[time_steps[t],i,n] * eta_battery 
                    - model.q_B_out[time_steps[t],i,n] / eta_battery
                    - model.SoC[time_steps[t],i,n] == 0)
        elif t > 0:
            return (model.SoC[time_steps[t-1],i,n] 
                    + model.q_B_in[time_steps[t],i,n] * eta_battery 
                    - model.q_B_out[time_steps[t],i,n] / eta_battery
                    - model.SoC[time_steps[t],i,n] == 0)
    model.SoC_con = Constraint(prosumer, index_time, years, 
                               rule = SoC_constraint_rule)
   
    def SoC_init_constraint_rule(model, i, n):
            return (model.SoC[time_steps[-1],i,n] == SoC_init)
    model.SoC_init_con = Constraint(prosumer, years, 
                                    rule = SoC_init_constraint_rule)
    
    
    # Inequality constraints: mu
    
    def mu_SoC_max_complementarity_rule_1(model, i, t, n):
        return (model.b[n,i] * prosumer_data.loc[SoC_max][i] 
                - model.SoC[t,i,n] >= 0)
    model.mu_SoC_max_compl_1 = Constraint(prosumer, time_steps, years,
                                          rule = mu_SoC_max_complementarity_rule_1)
    
    def mu_SoC_max_complementarity_rule_2(model, i, t, n):
        return (model.b[n,i] * prosumer_data.loc[SoC_max][i] 
                - model.SoC[t,i,n] <= (1-model.u_SoC_max[t,i,n])*model.M1)
    model.mu_SoC_max_compl_2 = Constraint(prosumer, time_steps, years, 
                                          rule = mu_SoC_max_complementarity_rule_2)
    
    def mu_SoC_max_complementarity_rule_3(model, i, t, n):
        return (model.mu_SoC[t,i,n] <= model.u_SoC_max[t,i,n]*model.M2)
    model.mu_SoC_max_compl_3 = Constraint(prosumer, time_steps, years, 
                                          rule = mu_SoC_max_complementarity_rule_3)
    
    def mu_q_B_in_complementarity_rule_1(model, i, t, n):
        return (model.b[n,i] * prosumer_data.loc[q_bat_max,i] 
                - model.q_B_in[t,i,n] >= 0)
    model.mu_q_B_in_compl_1 = Constraint(prosumer, time_steps, years, 
                                          rule = mu_q_B_in_complementarity_rule_1)
    
    def mu_q_B_in_complementarity_rule_2(model, i, t, n):
        return (model.b[n,i] * prosumer_data.loc[q_bat_max,i] 
                - model.q_B_in[t,i,n] <= (1-model.u_B_max_in[t,i,n])*model.M1)
    model.mu_q_B_in_compl_2 = Constraint(prosumer, time_steps, years, 
                                          rule = mu_q_B_in_complementarity_rule_2)
    
    def mu_q_B_in_complementarity_rule_3(model, i, t, n):
        return (model.mu_B_in[t,i,n] <= model.u_B_max_in[t,i,n]*model.M2)
    model.mu_q_B_in_compl_3 = Constraint(prosumer, time_steps, years, 
                                          rule = mu_q_B_in_complementarity_rule_3)
    
    def mu_q_B_out_complementarity_rule_1(model, i, t, n):
        return (model.b[n,i] * prosumer_data.loc[q_bat_min,i] 
                - model.q_B_out[t,i,n] >= 0)
    model.mu_q_B_out_compl_1 = Constraint(prosumer, time_steps, years, 
                                          rule = mu_q_B_out_complementarity_rule_1)
    
    def mu_q_B_out_complementarity_rule_2(model, i, t, n):
        return (model.b[n,i] * prosumer_data.loc[q_bat_min,i] 
                - model.q_B_out[t,i,n] <= (1-model.u_B_max_out[t,i,n])*model.M1)
    model.mu_q_B_out_compl_2 = Constraint(prosumer, time_steps, years, 
                                          rule = mu_q_B_out_complementarity_rule_2)
    
    def mu_q_B_out_complementarity_rule_3(model, i, t, n):
        return (model.mu_B_out[t,i,n] <= model.u_B_max_out[t,i,n]*model.M2)
    model.mu_q_B_out_compl_3 = Constraint(prosumer, time_steps, years, 
                                          rule = mu_q_B_out_complementarity_rule_3)
    
    
    # define welfare (with different parts  for simplification in the code)
    # community_welfare = {new_list: [] for new_list in prosumer}
    # prosumer_welfare = {new_list: [] for new_list in prosumer}
    # prosumer_welfare2 = {new_list: [] for new_list in prosumer}
    
    
    # costs = sum(
    #     grid_data.loc[t,'Residential'] * model.q_G_in[t,i,n]
    #     + sum(wtp[j][i][t] * model.q_share[t,j,i,n] for j in prosumer)
    #     - grid_data.loc[t,'DA'] * model.q_G_out[t,i,n]
    #     - sum(wtp[i][j][t] * model.q_share[t,i,j,n] for j in prosumer)
    #     for i in prosumer
    #     for t in time_steps
    #     for n in years)
    
    emissions = {}
    for n in years:
        emissions[n] = {}
        for i in prosumer:
            emissions[n][i] = sum(model.q_G_in[t,i,n]
                                  * weight[t]
                                  * grid_data.loc[t,'Emissions']
                                  / 1000000 
                                  for t in time_steps)
    
    # for i in prosumer:
    #     community_welfare[i] = sum(- grid_data.loc[t,'Residential']*model.q_G_in[t,i]*weight[t]
    #                                + grid_data.loc[t,'DA']*model.q_G_out[t,i]*weight[t] 
    #                                for t in time_steps)
    #     prosumer_welfare[i] = sum((grid_data.loc[t,'Residential'] 
    #                                + (prosumer_data.loc[w,j]
    #                                   * (1 - distances.loc[i,j]))
    #                                * grid_data.Emissions.loc[t] / 1000000)
    #                               * model.q_share[t,i,j]*weight[t] 
    #                               for j in prosumer 
    #                               for t in time_steps)
    #     prosumer_welfare2[i] = sum((grid_data.loc[t,'Residential'] 
    #                                 + (prosumer_data.loc[w,i]
    #                                    * (1 - distances.loc[j,i]))
    #                                 * grid_data.Emissions.loc[t] / 1000000)
    #                                * model.q_share[t,j,i]*weight[t] 
    #                                for j in prosumer 
    #                                for t in time_steps)
    
    #     # prosumer_welfare[i]: prosumer i sells to prosumer j
    #     # prosumer_welfare2[i]: prosumer i buys from prosumer j
    
    # costs = {new_list: [] for new_list in prosumer}    
    # for i in prosumer:
    #     costs[i] = (-community_welfare[i] - prosumer_welfare[i] 
    #                 + prosumer_welfare2[i])
        
    # Delta_costs = {new_list: [] for new_list in prosumer_old}    
    # for i in prosumer_old:
    #     Delta_costs[i] = ((costs[i]-results_old.loc[i,'costs']))
    #                       #/ abs(results_old.loc[i,'costs']))
    
    # emissions = {new_list: [] for new_list in prosumer}
    # for i in prosumer:
    #     emissions[i] = (sum(model.q_G_in[t,i]*weight[t]
    #                         *grid_data.Emissions.loc[t] / 1000000
    #                         for t in time_steps))
        
    # Delta_emissions = {new_list: [] for new_list in prosumer_old}    
    # for i in prosumer_old:
    #     Delta_emissions[i] = ((emissions[i]-results_old.loc[i,'emissions'])*1000)
    #                           #/ abs(results_old.loc[i,'emissions']))
        
    
    # objective functions
    
    # F1 ... overall emissions
    # F2 ... individual emissions
    # F3 ... individual costs
    # F4 ... individual weights on individual emissions and costs
    
    
    F1 = sum(emissions[n][i] - model.b[n,i] * emissions_wo_comm[i]
             for i in prosumer
             for n in years)
    
    F2 = sum(sum((emissions[n][i] - emissions_wo_comm[i]) * d[n][i] 
                 for n in years)
             * b_0[i]
             for i in prosumer)
    
    # F2 = sum(Delta_emissions[i] for i in prosumer_old)
    
    # F3 = sum(Delta_costs[i] for i in prosumer_old)
    
    # F4 = sum((alpha.loc[i,'alpha'] * Delta_costs[i])
    #          +((1-alpha.loc[i,'alpha']) * Delta_emissions[i])
    #          for i in prosumer_old)
    
    # F5 = sum(costs[i] for i in prosumer_old)
    
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
            a.append(value(sum(model.q_share[t,i,j,n]*weight[t]  
                               for t in time_steps)))
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
        
    costs_value = value(model.obj)
    b = []
    for n in years:
        for i in prosumer:
            b.append(value(model.b[n,i]))
    u = []
    for n in years:
        for i in prosumer:
            u.append(value(model.u[n,i]))

    return costs_value, b, u, q_share