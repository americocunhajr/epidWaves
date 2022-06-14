# -*- coding: utf-8  -*-
from myfunctions import *
import numpy as np
from lmfit import Model
import locale
locale.setlocale(locale.LC_ALL, 'en_EN.utf8')

# simulation information
case_name = 'COVID19_1st_wave_RJ'

df = load_data('COVID19_Data_RJ_Jan_01_2020_to_Dec_31_2021.csv')

# raw data start/end
RawDataStart = '2020-01-01' #1
RawDataEnd = '2020-07-01' #183

#  training data start/end
TrainDataStart = '2020-05-01' #122;
TrainDataEnd   = '2020-07-01' #183;

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

# maximum value for the incidence curve
tImax = time[Data_I_raw.idxmax()]+1

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

K0   = Cmax
ExpDataStart = '2020-04-01' #92
ExpDataEnd = '2020-05-01' #122
r0, B = np.polyfit(time[ExpDataStart:ExpDataEnd], np.log(Data_I_raw[ExpDataStart:ExpDataEnd]), 1)
tau0 = tImax

HyperParam = {}
# fixed model hyperparameter
tau = 119
HyperParam['tau'] = tau

# range of admissible values for model parameters
HyperParam['ub'] = np.array([[1.5*K0],
                             [10.0*r0]], dtype=np.float64)
HyperParam['lb'] = np.array([[0.5*K0],
                             [0.10*r0]], dtype=np.float64)

# number of initial guesses to fit the model
HyperParam['Ns'] = 30

# logistic model for total notifications
MyModel_C = Model(LogisticCDF)

# logistic model for new notifications per day
MyModel_I = Model(LogisticPDF)

# incidence curve fitting via Monte Carlo simulation
Result_I, ErrorObj_I = RegressionMC(time_train, Data_I_train, MyModel_I, HyperParam)
K_best = Result_I.best_values['K']
r_best = Result_I.best_values['r']
tau_best = Result_I.best_values['tau']

# model PDF
ModelPDF = LogisticPDF_model(tau)

# initial guess for MLE estimator
p0 = np.array([K_best, r_best])

# Akaike and Bayesian information criteria
[AIC, BIC] = AkaikeBIC(ModelPDF, time_train, p0)

# parameters
print("-----------------------")
print('K_best= ' + str(np.round(K_best, 0)))
print('r_best= ' + str(np.round(r_best, 3)))
print('tau_best= ' + str(tau_best))
print('RMSE= ' + str(np.round(ErrorObj_I['rmse'], 1)))
print('RSquare= ' + str(np.round(ErrorObj_I['rsquare'], 2)))
print('AIC= ' + str(np.round(AIC, 2)))
print('BIC= ' + str(np.round(BIC, 2)))
print("-----------------------")

# evaluate the model
MyFit_I_pred = LogisticPDF(time, K_best, r_best, tau)
MyFit_C_pred = LogisticCDF(time, K_best, r_best, tau)

# confidence envelope
p = np.array([K_best, r_best, tau_best])
MyFit_I_env_lower, MyFit_I_env_upper = predband(time, time_train, Data_I_train, p, LogisticPDF, conf=0.95)
MyFit_C_env_lower, MyFit_C_env_upper = np.cumsum(MyFit_I_env_lower), np.cumsum(MyFit_I_env_upper)

# upper bound estimation for starting date
tau_ast = np.where(MyFit_I_env_lower > 0.1)

# legend labels
graphobj = {}
graphobj['leg1'] = ' surveillance data'
graphobj['leg2'] = ' 7d moving average'
graphobj['leg3'] = ' statistical model'
graphobj['leg4'] = ' 95% confidence   '

#  Figure 1 - prevalence
graphobj['gname'] = str(case_name) + '__C_vs_time_tau_' + str(tau)
graphobj['gtitle'] = ''
graphobj['ymin']   = 0
graphobj['ymax']   = 15000
graphobj['xlab']   = []
graphobj['ylab']   = 'total reported deaths'

graph_C_1w(date, Data_C_raw, Data_C_MA, MyFit_C_env_upper, MyFit_C_env_lower, MyFit_C_pred, graphobj)

#  Figure 2 - prevalence
graphobj['gname'] = str(case_name) + '__I_vs_time_tau_' + str(tau)
graphobj['gtitle'] = ''
graphobj['ymin']   = 0
graphobj['ymax']   = 160
graphobj['taux']   = date[tau_ast[0][0]]
graphobj['tauy']   = 0.55
graphobj['tauh']   = 'right'
graphobj['xlab']   = []
graphobj['ylab']   = 'new reported deaths per day'

graph_I_1w(date, Data_I_raw, Data_I_MA, MyFit_I_env_upper, MyFit_I_env_lower, MyFit_I_pred, graphobj)



