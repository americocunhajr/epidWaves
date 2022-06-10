% -----------------------------------------------------------------
%  Main_COVID19_Visualization_waves_RJ.m
% -----------------------------------------------------------------
%  This program plots the COVID-19 cases and deaths curves from 
%  raw date for Rio de Janeiro city.
%  
%  Reference:
%  PRL Gianfelice, RS Oyarzabal, A Cunha Jr,
%  JMV Grzybowski, FC Batista, EEN Macau
%  The starting dates of COVID-19 multiple waves
%  Preprint, 2021
% -----------------------------------------------------------------
%  programmer: Americo Cunha Jr (UERJ)
%               
%  last update: Jan 17, 2022
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
disp(' COVID-19 raw data visualization                            ')
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
case_name = 'COVID19_waves_RJ';

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

% date range
DateStart = datenum('01-01-2020');
DateEnd   = datenum('01-01-2022');

% number of new events per day (incidence)
Data_Cases  = data_cases_by_symptoms; % cases organized by first symptoms
%Data_Cases  = data_cases_by_notifications; % cases organized by notification
Data_Deaths = data_deaths; % deaths by event ocurrency

% cumulative number of events (prevalence)
Data_Cases_cum  = cumsum(Data_Cases);
Data_Deaths_cum = cumsum(Data_Deaths);

% number of new events per week (incidence)
Data_Cases_w  = zeros(2*52,1);
Data_Deaths_w = zeros(2*52,1);
for i=1:2*52
    idx = (7*i-6):(7*i);
    Data_Cases_w(i)  = sum(Data_Cases(idx));
    Data_Deaths_w(i) = sum(Data_Deaths(idx));
end

% cumulative number of events (prevalence)
Data_Cases_cum_w  = cumsum(Data_Cases_w);
Data_Deaths_cum_w = cumsum(Data_Deaths_w);

% dataset size
N_data = length(Data_Cases);

% time vector
time = (1:N_data)';
% -----------------------------------------------------------



% process raw-data
% -----------------------------------------------------------
tic
disp(' '); 
disp(' --- processing raw-data --- ');
disp(' ');
disp('    ... ');
disp(' ');

% moving average (7 days) to remove fluctuations
Data_Cases_MA  = movmean(Data_Cases ,[6 0]);
Data_Deaths_MA = movmean(Data_Deaths,[6 0]);

Data_Cases_cum_MA  = cumsum(Data_Cases_MA);
Data_Deaths_cum_MA = cumsum(Data_Deaths_MA);

Data_Cases_MA_w  = movmean(Data_Cases_w ,[4 0]);
Data_Deaths_MA_w = movmean(Data_Deaths_w,[4 0]);

Data_Cases_cum_MA_w  = cumsum(Data_Cases_MA_w);
Data_Deaths_cum_MA_w = cumsum(Data_Deaths_MA_w);

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


% Figure 1 - incidence cases
% ..........................................................
graphobj.gname  = [num2str(case_name),'__cases'];
graphobj.gtitle = '';
graphobj.leg1   = ' surveillance data';
graphobj.leg2   = ' 7d moving average';
graphobj.xmin   = DateStart;
graphobj.xmax   = DateEnd;
graphobj.ymin   = 0;
graphobj.ymax   = 2500;
graphobj.xlab   = [];
graphobj.ylab   = 'new reported cases per day';
graphobj.flag   = 'eps';

fig1 = graph_I_raw(time,Data_Cases,...
                        Data_Cases_MA,...
                        graphobj);
% ..........................................................


% Figure 2 - incidence deaths
% ..........................................................
graphobj.gname  = [num2str(case_name),'__deaths'];
graphobj.gtitle = '';
graphobj.leg1   = ' surveillance data';
graphobj.leg2   = ' 7d moving average';
graphobj.xmin   = DateStart;
graphobj.xmax   = DateEnd;
graphobj.ymin   = 0;
graphobj.ymax   = 160;
graphobj.xlab   = [];
graphobj.ylab   = 'new reported deaths per day';
graphobj.flag   = 'eps';

fig2 = graph_I_raw(time,Data_Deaths,...
                        Data_Deaths_MA,...
                        graphobj);
% ..........................................................


% Figure 3 - incidence vs prevalence cases
% ..........................................................
graphobj.gname  = [num2str(case_name),'__total_vs_new_cases'];
graphobj.gtitle = '';
graphobj.leg1   = ' surveillance data';
graphobj.leg2   = ' 4w moving average';
graphobj.xmin   = 1e3;
graphobj.xmax   = 1e6;
graphobj.ymin   = 1e0;
graphobj.ymax   = 1e4;
graphobj.xlab   = 'total reported cases';
graphobj.ylab   = 'new reported cases per week';
graphobj.flag   = 'eps';

fig3 = graph_I_vs_C_raw(Data_Cases_cum_w   ,Data_Cases_w,...
                        Data_Cases_cum_MA_w,Data_Cases_MA_w,...
                        graphobj);
% ..........................................................


% Figure 4 - incidence vs prevalence deaths
% ..........................................................
graphobj.gname  = [num2str(case_name),'__total_vs_new_deaths'];
graphobj.gtitle = '';
graphobj.leg1   = ' surveillance data';
graphobj.leg2   = ' 4w moving average';
graphobj.xmin   = 1e2;
graphobj.xmax   = 1e5;
graphobj.ymin   = 1e0;
graphobj.ymax   = 1e3;
graphobj.xlab   = 'total reported deaths';
graphobj.ylab   = 'new reported deaths per week';
graphobj.flag   = 'eps';

fig4 = graph_I_vs_C_raw(Data_Deaths_cum_w   ,Data_Deaths_w,...
                        Data_Deaths_cum_MA_w,Data_Deaths_MA_w,...
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