# -*- coding: utf-8 -*-
"""
Created on Wed Jun  1 21:29:13 2022

@author: perger
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import seaborn as sns
sns.set_theme(palette='muted')

def plot_drop_in_out(scenarios, prosumer, years, b_0, b):

    fig, axs = plt.subplots(4,
                            figsize=(16, 12), 
                            tight_layout=True)
    
    fig.suptitle('Scenarios', fontsize = 18)
    j = 0
    
    # fig, ax = plt.subplots()
    
    for w in scenarios:
    
        dat = b[w]
        
        neu = pd.DataFrame(data=dat)
        accepted = [0,0,0,0,0]
        existing = [0,0,0,0,0]
        no_member = [0,0,0,0,0]
        dropping = [0,0,0,0,0]
        
        for i in prosumer:
            for n in years:
                if n == 1:
                    a = 2*neu.loc[i,n] - b_0[i]
                else:
                    a = 2*neu.loc[i,n] - neu.loc[i,n-1]
                if a == 2:
                    accepted[n-1] += 1
                elif a == 1:
                    existing[n-1] += 1
                elif a == 0:
                    no_member[n-1] -= 1
                elif a == -1:
                    dropping[n-1] -= 1
                    
        
        bar1 = list(np.array(accepted) + np.array(existing))
        bar2 = existing
        bar3 = list(np.array(dropping) + np.array(no_member))
        bar4 = no_member
        axs[j].bar(years, bar1, color = 'g')
        axs[j].bar(years, bar2, color = 'b')
        axs[j].bar(years, bar3, color = 'r')
        axs[j].bar(years, bar4, color = 'y')
        
        if j == 0:
            axs[j].set(ylabel='$\omega_1$')
        if j == 1:
            axs[j].set(ylabel='$\omega_2$')
        if j == 2:
            axs[j].set(ylabel='$\omega_3$')
        if j == 3:
            axs[j].set(ylabel='$\omega_4$')
            
        j += 1
    
        
    plt.xlabel('years', fontsize = 18)
        
    plt.legend(['joining', 'existing', 'leaving', 'no member'],
                    loc = 'lower left', 
                    bbox_to_anchor=(1, 0),
                    fontsize = 18)
        
    matplotlib.rcParams.update({'font.size': 22})
        
    plt.savefig('plot2.svg',  
                    bbox_inches='tight',  
                    dpi = 1200)
    
    return None

def plot_year_one_comparison(scenarios, prosumer, years, b_0, b_stoch, b_det, x_0):
    
    fig, axs = plt.subplots(2,
                            figsize=(16, 8), 
                            tight_layout=True)
    
    x_0 = dict(zip(prosumer, x_0))
    b_0 = {}
    for i in prosumer:
        if x_0[i] > 0:
          b_0[i] = 1
        else:
          b_0[i] = 0
    
    # without stochastic (DETERMINISTIC)
    
    b_1 = b_det[scenarios[0]][years[0]]
    
    accepted = len(prosumer)*[0]
    existing = len(prosumer)*[0]
    no_member = len(prosumer)*[0]
    dropping = len(prosumer)*[0]
    
    k = 0
    
    for i in prosumer:
        a = 2 * b_1[i] - b_0[i]
        
        if a == 2:
            accepted[k] = 1
        elif a == 1:
            existing[k] = 1
        elif a == 0:
            no_member[k] = 1
        elif a == -1:
            dropping[k] = 1
        
        k +=1
  
    axs[0].bar(prosumer, accepted, color = 'g')
    axs[0].bar(prosumer, existing, color = 'b')
    axs[0].bar(prosumer, dropping , color = 'r')
    axs[0].bar(prosumer, no_member, color = 'y')
    
    axs[0].axes.xaxis.set_ticks([])
    axs[0].set(ylabel='without scenarios')
    axs[0].axes.yaxis.set_ticks([])
    
    # with stochastic 
    
    b_1 = b_stoch[scenarios[0]][years[0]]
    
    accepted = len(prosumer)*[0]
    existing = len(prosumer)*[0]
    no_member = len(prosumer)*[0]
    dropping = len(prosumer)*[0]
    
    k = 0
    
    for i in prosumer:
        a = 2 * b_1[i] - b_0[i]
        
        if a == 2:
            accepted[k] = 1
        elif a == 1:
            existing[k] = 1
        elif a == 0:
            no_member[k] = 1
        elif a == -1:
            dropping[k] = 1
        
        k +=1
        
    axs[1].bar(prosumer, accepted, color = 'g')
    axs[1].bar(prosumer, existing, color = 'b')
    axs[1].bar(prosumer, dropping , color = 'r')
    axs[1].bar(prosumer, no_member, color = 'y')
    
    axs[1].set(ylabel='with scenarios')
    axs[1].axes.yaxis.set_ticks([])
    
    plt.xticks(rotation=90)
        
    plt.legend(['joining', 'existing', 'leaving', 'no member'],
               loc = 'lower left', 
               bbox_to_anchor=(1, 0),
               fontsize = 18)
    
    return None 