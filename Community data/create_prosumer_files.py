# -*- coding: utf-8 -*-
"""
Created on Tue Feb  8 14:27:25 2022

@author: perger

This programm creates prosumer files for different building types 
(and settlement patterns) using the IAMC format.
"""

import pandas as pd
from datetime import timedelta, datetime

IAMC_columns = ['model',
                'scenario',
                'region',
                'variable',
                'unit',
                'time',
                'value']

model = 'FRESH:COM v2.0'
scenario = 'Default scenario'
countries = ['Austria', 'Greece', 'Spain', 'Norway', 'UK']
year = 2019
building_type = 'LAB'
wtp = [100, 0.1, 90, 30, 50, 60, 40, 80, 20, 100]
num_prosumer = 10

if building_type in ['SAB', 'LAB']:
    prosumer_East = [3,5,7]
    prosumer_West = [2,6,10]
else:
    prosumer_East = []
    prosumer_West = []

load_norm = pd.read_excel('Settlement_10_SH.xlsx', sheet_name='SH_norm')
PV_South = pd.read_csv('PV_data_South.csv', sep=';') #  1kWp data
PV_East = pd.read_csv('PV_data_East.csv', sep=';') #  1kWp data
PV_West = pd.read_csv('PV_data_West.csv', sep=';') #  1kWp data

for i in range(1, num_prosumer+1):

    df = pd.DataFrame(columns = IAMC_columns)

    for j in countries:
        
        if j in ['Austria', 'Spain', 'Norway']:
            time_zone = '+01:00' # deviation from UTC (+01:00 is CET)
        if j == 'Greece':
            time_zone = '+02:00' # deviation from UTC (+02:00 is EET)
        if j == 'UK':
            time_zone = '+00:00' # deviation from UTC (+00:00 is GMT)
        start_date = str(year)+'-01-01 00:00' # YYYY-MM-DD HH:MM
        number_days = 365
        delta = timedelta(hours=1) # resolution ... hourly
    
        time_steps = []
        for t in range(24*number_days):
            time_steps.append((datetime
                               .fromisoformat(start_date+time_zone)+t*delta))
            
        # prosumer data
        _df_data = pd.DataFrame(columns = IAMC_columns)
        _df_data.value = [1,0,1,1,1,wtp[i-1]]
        _df_data.model = model
        _df_data.scenario = scenario
        _df_data.region = j
        _df_data.variable = [
            'Maximum Storage|Electricity|Energy Storage System',
            'Minimum Storage|Electricity|Energy Storage System',
            'Maximum Charge|Electricity|Energy Storage System',
            'Maximum Discharge|Electricity|Energy Storage System',
            'Maximum Active power|Electricity|Solar',
            'Price|Carbon'
            ]
        _df_data.unit = [
            'kWh',
            'kWh',
            'kW',
            'kW',
            'kW',
            'EUR/tCO2'
            ]
        _df_data.time = year
        
        df = df.append(_df_data, ignore_index = True)
    
        # electricity demand
        _df_load = pd.DataFrame(columns = IAMC_columns)
        _df_load.value = load_norm['SH '+str(i)]
        _df_load.model = model
        _df_load.scenario = scenario
        _df_load.region = j
        _df_load.variable = 'Final Energy|Residential and Commercial|Electricity'
        _df_load.unit = 'kWh'
        _df_load.time = time_steps
        
        df = df.append(_df_load, ignore_index = True)
    
        # PV generation
        _df_PV = pd.DataFrame(columns = IAMC_columns)
        if i in prosumer_East:
            _df_PV.value = PV_East[j]
        elif i in prosumer_West:
            _df_PV.value = PV_West[j]
        else:
            _df_PV.value = PV_South[j]
        _df_PV.model = model
        _df_PV.scenario = scenario
        _df_PV.region = j
        _df_PV.variable = 'Secondary Energy|Electricity|Solar|PV'
        _df_PV.unit = 'kWh'
        _df_PV.time = time_steps
        
        df = df.append(_df_PV, ignore_index = True)
        df.to_csv('Prosumer '+building_type+' '+str(i)+'.csv', index=False)
