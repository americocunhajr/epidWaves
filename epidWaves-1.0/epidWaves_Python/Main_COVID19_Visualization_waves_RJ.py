# -*- coding: utf-8  -*-
from myfunctions import *
import numpy as np
from datetime import datetime, timedelta
import locale
locale.setlocale(locale.LC_ALL, 'en_EN.utf8')

# simulation information
case_name = 'COVID19_waves_RJ'

#load surveillance data
df = load_data('COVID19_Data_RJ_Jan_01_2020_to_Dec_31_2021.csv')

# raw data start/end
RawDataStart = '2020-01-01'
RawDataEnd = '2021-12-31'

# number of new events per day (incidence)
Data_Cases = df.loc[RawDataStart:RawDataEnd]['data_inicio_sintomas'].astype(np.float64)
Data_Deaths = df.loc[RawDataStart:RawDataEnd]['data_obito'].astype(np.float64)
time_Cases = df.set_index('time')['data_inicio_sintomas']
time_Deaths = df.set_index('time')['data_obito']
time = df.loc[RawDataStart:RawDataEnd]['time'].astype(np.int32)
date = df.loc[RawDataStart:RawDataEnd]['date']

#cumulative number of events (prevalence)
Data_Cases_cum = np.cumsum(Data_Cases)
Data_Deaths_cum = np.cumsum(Data_Deaths)

# number of new events per week (incidence)
Data_Cases_w  = pd.DataFrame(data={'cases':np.zeros(2*52)})
Data_Deaths_w = pd.DataFrame(data={'deaths':np.zeros(2*52)})

for id in range(0, 2*52):
    Data_Cases_w['cases'][id]  = Data_Cases[id*7:id*7+7].sum()
    Data_Deaths_w['deaths'][id] = Data_Deaths[id*7:id*7+7].sum()

# cumulative number of events (prevalence)
Data_Cases_cum_w = np.cumsum(Data_Cases_w['cases'])
Data_Deaths_cum_w = np.cumsum(Data_Deaths_w['deaths'])

# moving average (7 days) to remove fluctuations
Data_Cases_MA = df.loc[RawDataStart:RawDataEnd]['data_inicio_sintomas'].rolling(7).mean()
Data_Deaths_MA = df.loc[RawDataStart:RawDataEnd]['data_obito'].rolling(7).mean()

Data_Cases_cum_MA  = np.cumsum(Data_Cases_MA)
Data_Deaths_cum_MA = np.cumsum(Data_Deaths_MA)

Data_Cases_MA_w  = Data_Cases_w['cases'].rolling(7).mean()
Data_Deaths_MA_w = Data_Deaths_w['deaths'].rolling(7).mean()

Data_Cases_cum_MA_w  = np.cumsum(Data_Cases_MA_w)
Data_Deaths_cum_MA_w = np.cumsum(Data_Deaths_MA_w)

graphobj = {}
# Figure 1 - incidence cases
graphobj['gname']  = str(case_name)+'__cases'
graphobj['gtitle'] = ''
graphobj['leg1']   = ' surveillance data'
graphobj['leg2']   = ' 7d moving average'
graphobj['ymin']   = 0
graphobj['ymax']   = 2500
graphobj['xlab']   = []
graphobj['ylab']   = 'new reported cases per day'
graph_I_raw(date, Data_Cases, Data_Cases_MA, graphobj)

#Figure 2 - incidence deaths
graphobj['gname']  = str(case_name)+'__deaths'
graphobj['gtitle'] = ''
graphobj['leg1']   = ' surveillance data'
graphobj['leg2']   = ' 7d moving average'
graphobj['ymin']   = 0
graphobj['ymax']   = 160
graphobj['xlab']   = []
graphobj['ylab']   = 'new reported deaths per day'
graph_I_raw(date, Data_Deaths, Data_Deaths_MA, graphobj)

# Figure 3 - incidence vs prevalence cases
graphobj['gname']  = str(case_name)+'__total_vs_new_cases'
graphobj['gtitle'] = ''
graphobj['leg1']   = ' surveillance data'
graphobj['leg2']   = ' 4w moving average'
graphobj['xmin']   = 1e3
graphobj['xmax']   = 1e6
graphobj['ymin']   = 1e0
graphobj['ymax']   = 1e4
graphobj['xlab']   = 'total reported cases'
graphobj['ylab']   = 'new reported cases per week'
graph_I_vs_C_raw(Data_Cases_cum_w, Data_Cases_w, Data_Cases_cum_MA_w, Data_Cases_MA_w, graphobj)

# Figure 4 - incidence vs prevalence deaths
graphobj['gname']  = str(case_name)+'__total_vs_new_deaths'
graphobj['gtitle'] = ''
graphobj['leg1']   = ' surveillance data'
graphobj['leg2']   = ' 4w moving average'
graphobj['xmin']   = 1e2
graphobj['xmax']   = 1e5
graphobj['ymin']   = 1e0
graphobj['ymax']   = 1e3
graphobj['xlab']   = 'total reported deaths'
graphobj['ylab']   = 'new reported deaths per week'
graph_I_vs_C_raw(Data_Deaths_cum_w, Data_Deaths_w, Data_Deaths_cum_MA_w, Data_Deaths_MA_w, graphobj)

