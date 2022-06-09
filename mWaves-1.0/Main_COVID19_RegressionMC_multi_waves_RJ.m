% -----------------------------------------------------------------
%  Main_COVID19_RegressionMC_multi_waves_RJ.m
% -----------------------------------------------------------------
%  This program constructs a data-driven algebraic statistical 
%  model to describe the evolution of COVID-19 notifications on 
%  a location which survailance data is available.
%  
%  Reference:
%  PRL Gianfelice, RS Oyarzabal, A Cunha Jr,
%  JMV Grzybowski, FC Batista, EEN Macau
%  The starting dates of COVID-19 multiple waves
%  Preprint, 2022
% -----------------------------------------------------------------
%  programmer: Americo Cunha Jr (UERJ)
%               
%  last update: Jan 17, 2021
% -----------------------------------------------------------------


clc
clear
close all

% program execution start time
% -----------------------------------------------------------
timeStart = tic();
% -----------------------------------------------------------


% program header
% -----------------------------------------------------------
disp(' ---------------------------------------------------------- ')
disp(' COVID-19 Model Estimation via Nonlinear Regression         ')
disp('                                                            ')
disp(' by                                                         ')
disp(' P. R. L. Gianfelice                                        ')
disp(' R. S. Oyarzabal                                            ')
disp(' A. Cunha Jr                                                ')
disp(' J. M. V. Grzybowski                                        ')
disp(' F. C. Batista                                              ')
disp(' E. E. N. Macau                                             ')
disp(' ---------------------------------------------------------- ')
% -----------------------------------------------------------



% simulation information
% -----------------------------------------------------------
case_name = 'COVID19_multi_waves_RJ';

disp(' '); 
disp([' Case Name: ',num2str(case_name)]);
disp(' ');
% -----------------------------------------------------------


% load surveillance data
% -----------------------------------------------------------
tic
disp(' '); 
disp(' --- loading surveillance data --- ');
disp(' ');
disp('    ... ');
disp(' ');

load('COVID19_Data_RJ_Jan_01_2020_to_Dec_31_2021.mat')

% range of dates
DateStart = datenum('01-01-2020');
DateEnd   = datenum('01-01-2022');

% indices to access the dates
% Jan 1, 2020 -   1   |   Jan 1, 2021 - 367
% Fev 1, 2020 -  32   |   Fev 1, 2021 - 398
% Mar 1, 2020 -  61   |   Mar 1, 2021 - 426
% Apr 1, 2020 -  92   |   Apr 1, 2021 - 457
% May 1, 2020 - 122   |   May 1, 2021 - 487
% Jun 1, 2020 - 153   |   Jun 1, 2021 - 518
% Jul 1, 2020 - 183   |   Jul 1, 2021 - 548
% Ago 1, 2020 - 214   |   Ago 1, 2021 - 579
% Sep 1, 2020 - 245   |   Sep 1, 2021 - 610
% Oct 1, 2020 - 275   |   Oct 1, 2021 - 640
% Nov 1, 2020 - 306   |   Nov 1, 2021 - 671
% Dec 1, 2020 - 336   |   Dec 1, 2021 - 701

% raw data start/end
RawDataStart =   1;  % Jan 01, 2020
RawDataEnd   = 731;  % Dec 31, 2021

% training data start/end
TrainDataStart = 122;  % Apr 01, 2020
TrainDataEnd   = 701;  % Dec 01, 2021

% new deaths per day (incidence)
Data_I_raw = data_deaths(RawDataStart:RawDataEnd);

% total deaths (prevalence)
Data_C_raw = cumsum(Data_I_raw);

% training data (incidence)
Data_I_train = data_deaths(TrainDataStart:TrainDataEnd);

% raw dataset size
N_data = length(Data_I_raw);

% training dataset size
N_train = length(Data_I_train);

% time vector
time = (1:N_data)';

% time vector for training
t_train_Start = TrainDataStart - RawDataStart;
t_train_End   = TrainDataEnd   - RawDataStart;
time_train    = (t_train_Start:t_train_End)';

toc
% -----------------------------------------------------------



% process raw data
% -----------------------------------------------------------
tic
disp(' '); 
disp(' --- processing raw-data --- ');
disp(' ');
disp('    ... ');
disp(' ');

% moving average (7 days) to remove fluctuations
Data_I_MA = movmean(Data_I_raw,[6 0]);
Data_C_MA = cumsum(Data_I_MA);

toc
% -----------------------------------------------------------



% statistical model estimation
% -----------------------------------------------------------
tic
disp(' '); 
disp(' --- statistical model estimation --- ');
disp(' ');
disp('    ... ');
disp(' ');

% admissible values interval for the initial guess
% Remarks:
% 1 - pay attention to the choice, it makes a lot of difference!
% 2 - choose values that make epidemiological sense!
% 3 - this choice is often trial and error game!
% 4 - need to be changed for every new dataset!

K0   = [ 9.70;  4.48;  6.54;  6.36;  8.11;  4.70]*1e3;
r0   = [53.00; 27.00; 52.00; 65.00; 32.00; 51.00]*1e-3;
tau0 = [  121;   268;   355;   465;   505;   608];

% starting dates of the epidemic waves
tau_ast = [60; 188; 300; 416; 430; 545];

% range of admissible values for model parameters
HyperParam.lb = [0.5*K0; 0.10*r0; 0.9*tau0];
HyperParam.ub = [1.5*K0; 10.0*r0; 1.1*tau0];

% number of initial guesses to fit the model
HyperParam.Ns = 30;

% logistic model for total notifications
MyModel_C = ...
    fittype(@(K1,K2,K3,K4,K5,K6,r1,r2,r3,r4,r5,r6,tau1,tau2,tau3,tau4,tau5,tau6,x) ...
            LogisticCDF6w(x,K1,K2,K3,K4,K5,K6,r1,r2,r3,r4,r5,r6,tau1,tau2,tau3,tau4,tau5,tau6),...
            'independent', 'x');

% logistic model for new notifications per day
MyModel_I = ...
    fittype(@(K1,K2,K3,K4,K5,K6,r1,r2,r3,r4,r5,r6,tau1,tau2,tau3,tau4,tau5,tau6,x) ...
            LogisticPDF6w(x,K1,K2,K3,K4,K5,K6,r1,r2,r3,r4,r5,r6,tau1,tau2,tau3,tau4,tau5,tau6),...
            'independent', 'x');

% incidence curve fitting via Monte Carlo simulation
[MyFit_I,ErrorObj_I] = RegressionMC(time_train,Data_I_train,...
                                     MyModel_I,HyperParam)
                                 
MyFit_I.K1

% prevalence curve fitting 
MyFit_C = cfit(MyModel_C,...
               MyFit_I.K1,MyFit_I.K2,MyFit_I.K3,MyFit_I.K4,MyFit_I.K5,MyFit_I.K6,...
               MyFit_I.r1,MyFit_I.r2,MyFit_I.r3,MyFit_I.r4,MyFit_I.r5,MyFit_I.r6,...
               MyFit_I.tau1,MyFit_I.tau2,MyFit_I.tau3,MyFit_I.tau4,MyFit_I.tau5,MyFit_I.tau6);
           
% incidence curve fitting 
% MyFit_I = cfit(MyModel_I,...
%                K0(1),K0(2),K0(3),K0(4),K0(5),...
%                r0(1),r0(2),r0(3),r0(4),r0(5));

% % model PDF
% ModelPDF = @(x,p) LogisticPDF(x,p(1:Nw),p(Nw+1:2*Nw),tau0,Nw);
% 
% % initial guess for MLE estimator
% p0 = [MyFit_I.K; MyFit_I.r];
% 
% % Akaike and Bayesian information criteria
% [AIC,BIC] = AkaikeBIC(ModelPDF,time_train,p0)

toc
% -----------------------------------------------------------



% statistical model evaluation
% -----------------------------------------------------------
tic
disp(' '); 
disp(' --- statistical model evaluation --- ');
disp(' ');
disp('    ... ');
disp(' ');


% evaluate the model
MyFit_I_pred = MyFit_I(time);
MyFit_C_pred = MyFit_C(time);

% confidence envelope
MyFit_I_env = predint(MyFit_I,time,0.95,'observation','on'); %'functional' 'observation'
MyFit_C_env = cumsum(MyFit_I_env);

% upper bound estimation for starting date
%tau_ast = DateStart + find(MyFit_I_env(:,1)>0,1);

toc
% -----------------------------------------------------------



% post-processing
% -----------------------------------------------------------
tic
disp(' ')
disp(' --- post-processing --- ');
disp(' ');
disp('    ... ');
disp(' ');


% adjust time vector for date format
% ..........................................................
time    = linspace(DateStart,DateEnd,length(time))';
tau_ast = DateStart + tau_ast;
% ..........................................................


% Figure 1 - incidence of deaths
% ..........................................................
graphobj.gname  = [num2str(case_name),'__I_vs_time_6w'];
graphobj.gtitle = '';
graphobj.leg1   = ' surveillance data';
graphobj.leg2   = ' 7d moving average';
graphobj.leg3   = ' statistical model';
graphobj.leg4   = ' 95% confidence   ';
graphobj.xmin   = DateStart;
graphobj.xmax   = DateEnd;
graphobj.ymin   = 0;
graphobj.ymax   = 250;
graphobj.tauy   = 0.7;
graphobj.xlab   = [];
graphobj.ylab   = 'total reported deaths';
graphobj.flag   = 'eps';

fig1 = graph_I_6w(time,Data_I_raw,...
                       Data_I_MA,...
                       MyFit_I_env(:,2),...
                       MyFit_I_env(:,1),...
                       MyFit_I_pred,...
                       tau_ast,...
                       graphobj);
% ..........................................................

toc
% -----------------------------------------------------------


% program execution time
% -----------------------------------------------------------
disp(' ');
disp(' -----------------------------');
disp('            THE END!          ');
disp(' -----------------------------');
disp('  Total execution time:       ');
disp(['  ',num2str(toc(timeStart)),' seconds']);
disp(' -----------------------------');
% -----------------------------------------------------------
