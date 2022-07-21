# -*- coding: utf-8 -*-
"""
Created on Wed Jun  1 21:29:13 2022

@author: perger
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
from pathlib import Path
import seaborn as sns
sns.set_theme(palette='muted')

def heatmap(q_share, filename=None):
    
    sns.set_theme(palette='muted')   
    
    # labels_old = q_share.columns.tolist() 
    # labels = []
    # for i in labels_old:
    #     labels.append(i.removeprefix('Prosumer '))
        
    fig, ax = plt.subplots(figsize=(6,6))
    white = [(0.9764705882352941,0.9764705882352941,0.9764705882352941)]
    cmap = sns.color_palette("RdPu", 15)
    ax = sns.heatmap(q_share.round(decimals=0), 
                     annot=True, 
                     cmap=cmap,  
                     fmt='g', square=True, cbar=False) # xticklabels
    sns.heatmap(q_share.round(decimals=0), 
                cmap=white, 
                linewidths=0.5, 
                vmin=0, 
                vmax=2, 
                mask=q_share.round(decimals=0) > 0, 
                cbar=False, 
                ax=ax, square=True)
    #heatmap.set_xticklabels(heatmap.get_xticklabels(), rotation=90)
           
    if filename:
            plt.savefig(Path('graphics') / filename,  
                        bbox_inches='tight', 
                        format = "pdf", 
                        dpi = 1200)
    else:
            plt.show()
            
def costs_emissions(results_new, results_old, results,
                    prosumer_old, prosumer_new, prosumer,
                    filename=None):
    
    colors = sns.color_palette("muted", 10).as_hex()
    
    prs = []
    for i in prosumer:
        if i in prosumer_old: 
            prs.append(i)
        elif i in prosumer_new:
            prs.append(i)
            
    delta_costs = []
    delta_emissions = []
    for i in prs:
        if i in prosumer_old and i in prosumer_new: 
            delta_costs.append((results_new.loc[i,'costs']
                            -results_old.loc[i,'costs']))
            delta_emissions.append((results_new.loc[i,'emissions']
                                -results_old.loc[i,'emissions']))
        elif i in prosumer_old:
            delta_costs.append((results.loc[i,'costs']
                            -results_old.loc[i,'costs']))
            delta_emissions.append((results.loc[i,'emissions']
                                -results_old.loc[i,'emissions']))
        elif i in prosumer_new:
            delta_costs.append((results_new.loc[i,'costs']
                                -results.loc[i,'costs']))
            delta_emissions.append((results_new.loc[i,'emissions']
                                    -results.loc[i,'emissions']))
    # Prepare data for bars
    barWidth = 0.4
    x = np.arange(1, len(delta_costs)+1)
    x1 = [_x - 0.5*barWidth for _x in x]
    x2 = [_x + 0.5*barWidth for _x in x]
    
    y1max = max(abs(i) for i in delta_costs)*1.1
    y2max = max(abs(i) for i in delta_emissions)*1.1

    fig, ax1 = plt.subplots()
    xlabels = []
    for i in prs:
        xlabels.append(i.removeprefix('Prosumer '))
    
    
    # Set x-axis 
    plt.xticks(x, xlabels, rotation=90)
    # plt.xlabel('Prosumer')
     
    # Plot left axis
    ax1.bar(x1, delta_costs, width = barWidth, color = colors[3], label='Costs')
    
    # Plot right axis
    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
    ax2.bar(x2, delta_emissions, width = barWidth, color = colors[2], label='Emissions')
        
    # Adding labels
    ax1.set_ylabel('Increase/decrease in EUR', color = colors[3])
    ax2.set_ylabel('Increase/decrease in tCO2', color = colors[2])
    
    # Color y-axes
    ax1.tick_params(axis='y', labelcolor=colors[3])
    ax2.tick_params(axis='y', labelcolor=colors[2])
    
    # Setting axes limits
    ax1.set_ylim(-y1max,y1max)
    ax2.set_ylim(-y2max,y2max)

    # Setting grid
    ax1.xaxis.grid(True)
    ax2.yaxis.grid(False)
    
    # Adding legend
    # ax1.legend(loc='lower left') # costs
    # ax2.legend(loc='upper left') # emissions
    
    lines_1, labels_1 = ax1.get_legend_handles_labels()
    lines_2, labels_2 = ax2.get_legend_handles_labels()

    lines = lines_1 + lines_2
    labels = labels_1 + labels_2
    
    plt.legend(lines, labels, loc='lower center', bbox_to_anchor=(.5, 1), ncol=5, prop={'size': 10})
    # plt.legend(loc='lower center', bbox_to_anchor=(.5, 1), ncol=5, prop={'size': 10})
    
    #ax1.set_ylim(-300,500)
    #ax2.set_ylim(-0.54,0.9)

    if filename:
            plt.savefig(Path('graphics') / filename,  
                        bbox_inches='tight', 
                        format = "pdf", 
                        dpi = 1200)
    else:
            plt.show()
            
            
def plot_emissions(emissions1, emissions, emissions_old, 
                   years, x_0, scenarios,
                   filename=None):
    
    sns.set_theme(palette='muted')
    colors = (sns.color_palette('cool', len(scenarios)))
    
    prosumer = list(x_0.keys())
    em_old = 0
    for i in prosumer:
        if x_0[i] > 0:
            em_old += emissions_old[i]

    fig, ax = plt.subplots(figsize=(8,5))
    ax.plot(years, [em_old]*len(years), linestyle='--', color='k', label='w/o community')
    
    j = 0
    for w in scenarios:
                   
        em = {}
        for n in years:
            if n == 1:
                em[n] = 0
                for i in prosumer:
                    if x_0[i] > 0:
                        # em[n] += emissions1[i]
                        em[n] += emissions[w][n][i]
            else:
                em[n] = 0
                for i in prosumer:
                    if x_0[i] > 0:
                        em[n] += emissions[w][n][i]
    
        ax.plot(years, list(em.values()), linestyle='-', marker = 'o', 
                markerfacecolor = colors[j], color=colors[j],
                label='$\omega$'+str(j+1))
        
        j += 1
    
    plt.legend(loc='lower center', bbox_to_anchor=(.5, 1), ncol=5, prop={'size': 10})
    plt.xlabel('Years', fontsize=12)
    plt.ylabel('Annual emissions in tCO2', fontsize=12)
    ax.set_xticks(years)
    
    
    
    if filename:
            plt.savefig(Path('graphics') / filename,  
                        bbox_inches='tight', 
                        format = "pdf", 
                        dpi = 1200)
    else:
            plt.show()
    

def plot_drop_in_out(scenarios, prosumer, years, x_0, b, filename=None):

    b_0 = {}
    for i in prosumer:
        if x_0[i] > 0:
          b_0[i] = 1
        else:
          b_0[i] = 0

    fig, axs = plt.subplots(4, figsize=(8,6),
                            tight_layout=True)

    j = 0    
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
        
        axs[j].set_ylabel('$\omega$'+str(j+1), fontsize=12)
           
        j += 1
            
    plt.xlabel('Years', fontsize=12)
    
    axs[0].legend(['joining', 'member', 'leaving', 'no member (portfolio)'], 
               loc='lower center', 
               bbox_to_anchor=(0.5, 1),
               ncol=4,
               prop={'size': 12}) 
                
    if filename:
            plt.savefig(Path('graphics') / filename,  
                        bbox_inches='tight', 
                        format = "pdf", 
                        dpi = 1200)
    else:
            plt.show()
    
    return None

def plot_one_year_comparison(scenarios, prosumer, years, 
                             b_prev1, b_next1, 
                             b_prev2, b_next2, 
                             filename=None):
    
    fig, axs = plt.subplots(2,
                            figsize=(8,6), 
                            tight_layout=True)
    
    
    # first graph
    
    accepted = len(prosumer)*[0]
    existing = len(prosumer)*[0]
    no_member = len(prosumer)*[0]
    dropping = len(prosumer)*[0]
    
    k = 0
    
    for i in prosumer:
        a = 2 * b_next1[i] - b_prev1[i]
        
        if a == 2:
            accepted[k] = 1
        elif a == 1:
            existing[k] = 1
        elif a == 0:
            no_member[k] = -1
        elif a == -1:
            dropping[k] = -1
        
        k +=1
  
    axs[0].bar(prosumer, accepted, color = 'g')
    axs[0].bar(prosumer, existing, color = 'b')
    axs[0].bar(prosumer, dropping , color = 'r')
    axs[0].bar(prosumer, no_member, color = 'y')
    
    axs[0].axes.xaxis.set_ticks([])
    axs[0].set(ylabel='stochastic')
    axs[0].axes.yaxis.set_ticks([])
    plt.grid()
    
    axs[0].legend(['joining', 'member', 'leaving', 'no member (portfolio)'], 
               loc='lower center', 
               bbox_to_anchor=(0.5, 1),
               ncol=4,
               prop={'size': 12}) 
    
    # second graph 
       
    accepted = len(prosumer)*[0]
    existing = len(prosumer)*[0]
    no_member = len(prosumer)*[0]
    dropping = len(prosumer)*[0]
    
    k = 0
    
    for i in prosumer:
        a = 2 * b_next2[i] - b_prev2[i]
        
        if a == 2:
            accepted[k] = 1
        elif a == 1:
            existing[k] = 1
        elif a == 0:
            no_member[k] = -1
        elif a == -1:
            dropping[k] = -1
        
        k +=1
        
    axs[1].bar(prosumer, accepted, color = 'g')
    axs[1].bar(prosumer, existing, color = 'b')
    axs[1].bar(prosumer, dropping , color = 'r')
    axs[1].bar(prosumer, no_member, color = 'y')
    
    axs[1].set(ylabel='deterministic')
    axs[1].axes.yaxis.set_ticks([])
    
    xticks = range(0,len(prosumer))
    xlabels = []
    for i in prosumer:
        xlabels.append(i.removeprefix('Prosumer '))
    
    plt.xticks(xticks, xlabels, rotation=90)
    
    axs[0].set_ylim(-1,1)
    axs[1].set_ylim(-1,1)
    axs[0].set_xlim(-0.5,len(prosumer)-0.5)
    axs[1].set_xlim(-0.5,len(prosumer)-0.5)
    
    if filename:
            plt.savefig(Path('graphics') / filename,  
                        bbox_inches='tight', 
                        format = "pdf", 
                        dpi = 1200)
    else:
            plt.show()
    
    return None

def plot_drop_in_out_years(scenarios, prosumer, years, x_0, b, filename=None):

    b_0 = {}
    for i in prosumer:
        if x_0[i] > 0:
          b_0[i] = 1
        else:
          b_0[i] = 0

    fig, axs = plt.subplots(4, figsize=(8,6),
                            tight_layout=True)
    
    # fig.suptitle('Scenarios', fontsize = 18)

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
        
        # if j == 0:
        axs[j].set_ylabel('$\omega$'+str(j+1), fontsize=12)
        # if j == 1:
        #     axs[j].set_ylabel('$\omega_2$', fontsize=12)
        # if j == 2:
        #     axs[j].set_ylabel('$\omega_3$', fontsize=12)
        # if j == 3:
        #     axs[j].set_ylabel('$\omega_4$', fontsize=12)
            
        # axs[j].set(ylabel='$\omega_4$')
            
        j += 1
    
        
    plt.xlabel('Years', fontsize=12)
        
    # plt.legend(['joining', 'existing', 'leaving', 'no member'],
    #                 loc = 'lower left', 
    #                 bbox_to_anchor=(1, 0),
    #                 fontsize = 18)
    
    axs[0].legend(['joining', 'member', 'leaving', 'no member (portfolio)'], 
               loc='lower center', 
               bbox_to_anchor=(0.5, 1),
               ncol=4,
               facecolor="white",
               prop={'size': 12}) 
        
        
    if filename:
            plt.savefig(Path('graphics') / filename,  
                        bbox_inches='tight', 
                        format = "pdf", 
                        dpi = 1200)
    else:
            plt.show()
    
    return None
