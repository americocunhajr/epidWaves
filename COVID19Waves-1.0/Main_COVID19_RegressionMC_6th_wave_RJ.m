% -----------------------------------------------------------------
%  Main_COVID19_RegressionMC_5th_wave_RJ.m
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
case_name = 'COVID19_6th_wave_RJ';

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
DateStart = datenum('07-01-2021');
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
RawDataStart = 548;  % Jul 01, 2021
RawDataEnd   = 731;  % Dec 31, 2021

% training data start/end
TrainDataStart = 579;  % Ago 01, 2021
TrainDataEnd   = 701;  % Dec 01, 2021

% new deaths per day (incidence)
Data_I_raw = data_deaths(RawDataStart:RawDataEnd);

% total deaths (prevalence)
Data_C_raw = cumsum(Data_I_raw);

% training data (incidence)
Data_I_train = data_deaths(TrainDataStart:TrainDataEnd);

% maximum value for the incidence curve
[Imax,tImax] = max(Data_I_raw);

% maximum value for the prevalence curve
Cmax = max(Data_C_raw);

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

K0   = Cmax;
r0   = RegressionExp(time(1:50),Data_I_raw(1:50));
tau0 = tImax;

% fixed model hyperparameter
tau = 58;

% range of admissible values for model parameters
HyperParam.lb = [0.5*K0; 0.10*r0];
HyperParam.ub = [1.5*K0; 10.0*r0];

% number of initial guesses to fit the model
HyperParam.Ns = 30;

% logistic model for total notifications
MyModel_C = fittype( @(K,r,x) LogisticCDF(x,K,r,tau),'independent', 'x');

% logistic model for new notifications per day
MyModel_I = fittype( @(K,r,x) LogisticPDF(x,K,r,tau),'independent', 'x');

% incidence curve fitting via Monte Carlo simulation
[MyFit_I,ErrorObj_I] = RegressionMC(time_train,Data_I_train,...
                                    MyModel_I,HyperParam)
                           
% prevalence curve fitting 
MyFit_C = cfit(MyModel_C,MyFit_I.K,MyFit_I.r);

% model PDF
ModelPDF = @(x,p) LogisticPDF(x,p(1),p(2),tau);

% initial guess for MLE estimator
p0 = [MyFit_I.K; MyFit_I.r];

% Akaike and Bayesian information criteria
[AIC,BIC] = AkaikeBIC(ModelPDF,time_train,p0)

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
tau_ast = DateStart + find(MyFit_I_env(:,1)>0,1);

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
time = linspace(DateStart,DateEnd,length(time))';
% ..........................................................

% custom colors
% ..........................................................
MyGray = [0.8 0.8 0.8];
% ..........................................................

% legend labels
% ..........................................................
label_data  = ' surveillance data';
label_MA    = ' 7d moving average';
label_model = ' statistical model';
label_env   = ' 95% confidence   ';
% ..........................................................



% ..........................................................
figure(1)
fig1a = fill([time' fliplr(time')],...
             [MyFit_C_env(:,2)' fliplr(MyFit_C_env(:,1)')],MyGray);
hold on
fig1b = plot(time,Data_C_raw  ,'o-m','LineWidth',1);
fig1c = plot(time,Data_C_MA   ,'.-g','LineWidth',3);
fig1d = plot(time,MyFit_C_pred,' -b','LineWidth',3);
hold off
set(fig1a,'DisplayName',label_env  );
set(fig1b,'DisplayName',label_data );
set(fig1c,'DisplayName',label_MA   );
set(fig1d,'DisplayName',label_model);
leg = [fig1b; fig1c; fig1d; fig1a];
leg = legend(leg,'Location','Best');
ylabel('total reported deaths');
datetick('x',28,'keeplimits');
xlim([DateStart DateEnd]); ylim([0 max(Cmax,max(MyFit_C_env(:,2)))]);
set(leg,'FontSize',18);
set(0,'DefaultAxesFontSize',18);
gname = [num2str(case_name),'__C_vs_time_tau_',num2str(tau0)];
saveas(gcf,gname,'epsc2');
% ..........................................................


% ..........................................................
figure(2)
%ymax = max(Imax,max(MyFit_I_env(:,2)));
ymax = 160;
fig2a = fill([time' fliplr(time')],...
             [MyFit_I_env(:,2)' fliplr(MyFit_I_env(:,1)')],MyGray);
hold on
fig2b = plot(time,Data_I_raw  ,'o-m','LineWidth',1);
fig2c = plot(time,Data_I_MA   ,'.-g','LineWidth',3);
fig2d = plot(time,MyFit_I_pred,' -b','LineWidth',3);
fig2e = plot([tau_ast tau_ast],[0 0.85*ymax],'--k','LineWidth',2);
hold off
set(fig2a,'DisplayName',label_env  );
set(fig2b,'DisplayName',label_data );
set(fig2c,'DisplayName',label_MA   );
set(fig2d,'DisplayName',label_model);
leg = [fig2b; fig2c; fig2d; fig2a];
leg = legend(leg,'Location','Best');
ylabel('new reported deaths per day');
datetick('x',28,'keeplimits');
xlim([DateStart DateEnd]); ylim([0 ymax]);
xpos = find(MyFit_I_env(:,1)>0,1)/length(time);
ypos = 0.87;
text(xpos,ypos, datestr(tau_ast), ...
     'Units', 'normalized', ...
     'HorizontalAlignment','left', ...
     'VerticalAlignment'  ,'bottom',...
     'FontSize', 16);
set(leg,'FontSize',18);
set(0,'DefaultAxesFontSize',18);
gname = [num2str(case_name),'__I_vs_time_tau_',num2str(tau0)];
saveas(gcf,gname,'epsc2');
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
