# -*- coding: utf-8  -*-
from myfunctions import *
import numpy as np
import matplotlib.pyplot as plt
from lmfit import Model, printfuncs, conf_interval
from scipy.optimize import curve_fit
import locale
locale.setlocale(locale.LC_ALL, 'en_EN.utf8')

# simulation information
case_name = 'COVID19_multi_waves_RJ'

df = load_data('COVID19_Data_RJ_Jan_01_2020_to_Dec_31_2021.csv')

# raw data start/end
RawDataStart = '2020-01-01' #1
RawDataEnd = '2021-12-31' #731

#  training data start/end
TrainDataStart = '2020-04-01' #122
TrainDataEnd   = '2021-12-01' #701

# new deaths per day (incidence)
Data_I_raw = df.loc[RawDataStart:RawDataEnd]['data_obito'].astype(np.float64)

#total deaths (prevalence)
Data_C_raw = np.cumsum(Data_I_raw)

# training data (incidence)
Data_I_train = df.loc[TrainDataStart:TrainDataEnd]['data_obito'].astype(np.float64)

# maximum value for the prevalence curve
Cmax = max(Data_C_raw)

# raw dataset size
N_data = Data_I_raw.size

# training dataset size
N_train = Data_I_train.size

# time vector
time = df.loc[RawDataStart:RawDataEnd]['time'].astype(np.int32)
date = df.loc[RawDataStart:RawDataEnd]['date']

# time vector for training
time_train = time[TrainDataStart:TrainDataEnd]

# moving average (7 days) to remove fluctuations
Data_I_MA = df.loc[RawDataStart:RawDataEnd]['data_obito'].rolling(7).mean()
Data_C_MA = np.cumsum(Data_I_MA)

# -----------------------------------------------------------
# statistical model estimation
# -----------------------------------------------------------
# admissible values interval for the initial guess
# Remarks:
# 1 - pay attention to the choice, it makes a lot of difference!
# 2 - choose values that make epidemiological sense!
# 3 - this choice is often trial and error game!
# 4 - need to be changed for every new dataset!

K0   = [ 9700,  4480,  6540,  6360,  8110,  4700]
r0   = [.053, .027, .052, .065, .032, .051]
t0 = [  121,   268,   355,   465,   505,   608]

# starting dates of the epidemic waves
#tau_ast = [60, 188, 300, 416, 430, 545]
tau_ast = ['2020-02-28', '2020-07-06', '2020-10-26', '2021-02-19', '2021-05-03', '2021-06-28']

HyperParam = {}
# range of admissible values for model parameters
HyperParam['p'] = ['K1', 'K2', 'K3', 'K4', 'K5', 'K6', 'r1', 'r2', 'r3', 'r4', 'r5', 'r6', 'tau1', 'tau2', 'tau3', 'tau4', 'tau5', 'tau6']
HyperParam['lb'] = np.array([[0.5*K0[0]], [0.5*K0[1]], [0.5*K0[2]], [0.5*K0[3]], [0.5*K0[4]], [0.5*K0[5]],
                             [0.10*r0[0]],[0.10*r0[1]],[0.10*r0[2]],[0.10*r0[3]],[0.10*r0[4]],[0.10*r0[5]],
                             [0.9*t0[0]], [0.9*t0[1]], [0.9*t0[2]], [0.9*t0[3]], [0.9*t0[4]], [0.9*t0[5]]
                             ], dtype=np.float64)
HyperParam['ub'] = np.array([[1.5*K0[0]], [1.5*K0[1]], [1.5*K0[2]], [1.5*K0[3]], [1.5*K0[4]], [1.5*K0[5]],
                             [10.0*r0[0]],[10.0*r0[1]],[10.0*r0[2]],[10.0*r0[3]],[10.0*r0[4]],[10.0*r0[5]],
                             [1.1*t0[0]], [1.1*t0[1]], [1.1*t0[2]], [1.1*t0[3]], [1.1*t0[4]], [1.1*t0[5]]
                             ], dtype=np.float64)

# number of initial guesses to fit the model
HyperParam['Ns'] = 30

# logistic model for new notifications per day
MyModel_I = Model(LogisticPDF6w)

# incidence curve fitting via Monte Carlo simulation
Result_I, ErrorObj_I = RegressionMC(time_train, Data_I_train, MyModel_I, HyperParam)

# parameters
print("-----------------------")
print('RMSE= ' + str(np.round(ErrorObj_I['rmse'], 1)))
print('RSquare= ' + str(np.round(ErrorObj_I['rsquare'], 2)))
print("-----------------------")

# evaluate the model
p = Result_I.best_values
MyFit_I_pred = LogisticPDF6w(time,
                             p['K1'], p['K2'], p['K3'], p['K4'], p['K5'], p['K6'],
                             p['r1'], p['r2'], p['r3'], p['r4'], p['r5'], p['r6'],
                             p['tau1'], p['tau2'], p['tau3'], p['tau4'], p['tau5'], p['tau6'])

# confidence envelope
p = np.array([p['K1'], p['K2'], p['K3'], p['K4'], p['K5'], p['K6'],
              p['r1'], p['r2'], p['r3'], p['r4'], p['r5'], p['r6'],
              p['tau1'], p['tau2'], p['tau3'], p['tau4'], p['tau5'], p['tau6']])

MyFit_I_env_lower, MyFit_I_env_upper = predband(time, time_train, Data_I_train, p, LogisticPDF6w, conf=0.95)

# legend labels
graphobj = {}
graphobj['leg1'] = ' surveillance data'
graphobj['leg2'] = ' 7d moving average'
graphobj['leg3'] = ' statistical model'
graphobj['leg4'] = ' 95% confidence   '

#  Figure 1 - prevalence
graphobj['gname'] = str(case_name) + '__I_vs_time_6w'
graphobj['gtitle'] = ''
graphobj['ymin']   = 0
graphobj['ymax']   = 200
graphobj['tauy']   = 0.5
graphobj['xlab']   = []
graphobj['ylab']   = 'new reported deaths per day'

graph_I_6w(date, Data_I_raw, Data_I_MA, MyFit_I_env_upper, MyFit_I_env_lower, MyFit_I_pred, tau_ast, graphobj)



