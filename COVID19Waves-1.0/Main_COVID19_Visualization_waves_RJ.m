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
%Data_Cases_cum  = cumsum(Data_Cases);
Data_Deaths_cum = cumsum(Data_Deaths);

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
label_end   = ' delayed data     ';
% ..........................................................

% ..........................................................
figure(1)
fig1a = plot(time,Data_Cases   ,'om','LineWidth',1);
hold on
fig1b = plot(time,Data_Cases_MA,'.-g','LineWidth',3);
hold off
set(fig1a,'DisplayName',label_data );
set(fig1b,'DisplayName',label_MA   );
leg = [fig1a; fig1b];
leg = legend(leg,'Location','Best');
ylabel('new reported cases per day');
datetick('x',28,'keeplimits');
xlim([DateStart DateEnd]); %ylim([0 600]);
set(leg,'FontSize',18);
set(0,'DefaultAxesFontSize',18);
gname = [num2str(case_name),'__cases'];
saveas(gcf,gname,'epsc2');
% ..........................................................


% ..........................................................
figure(2)
fig2a = plot(time,Data_Deaths  ,'om','LineWidth',1);
hold on
fig2b = plot(time,Data_Deaths_MA,'.-g','LineWidth',3);
hold off
set(fig2a,'DisplayName',label_data );
set(fig2b,'DisplayName',label_MA   );
leg = [fig2a; fig2b];
leg = legend(leg,'Location','Best');
ylabel('new reported deaths per day');
datetick('x',28,'keeplimits');
xlim([DateStart DateEnd]); %ylim([0 300]);
set(leg,'FontSize',18);
set(0,'DefaultAxesFontSize',18);
gname = [num2str(case_name),'__deaths'];
saveas(gcf,gname,'epsc2');
% ..........................................................


% ..........................................................
figure(3)
fig3a = loglog(Data_Cases_cum(90:end),Data_Cases(90:end)   ,'om','LineWidth',1);
hold on
fig3b = loglog(Data_Cases_cum_MA(90:end),Data_Cases_MA(90:end),'.-g','LineWidth',3);
hold off
set(fig3a,'DisplayName',label_data );
set(fig3b,'DisplayName',label_MA   );
leg = [fig3a; fig3b];
leg = legend(leg,'Location','Best');
xlabel('total reported cases');
ylabel('new reported cases per day');
xlim([1e3 1e6]); ylim([1 1e4]);
set(leg,'FontSize',18);
set(0,'DefaultAxesFontSize',18);
gname = [num2str(case_name),'__total_vs_new_cases'];
saveas(gcf,gname,'epsc2');
% ..........................................................



% ..........................................................
figure(4)
fig4a = loglog(Data_Deaths_cum(90:end),Data_Deaths(90:end)   ,'om','LineWidth',1);
hold on
fig4b = loglog(Data_Deaths_cum_MA(90:end),Data_Deaths_MA(90:end),'.-g','LineWidth',3);
hold off
set(fig4a,'DisplayName',label_data );
set(fig4b,'DisplayName',label_MA   );
leg = [fig4a; fig4b];
leg = legend(leg,'Location','Best');
xlabel('total reported deaths');
ylabel('new reported deaths per day');
xlim([1e2 1e5]); ylim([1 1e3]);
set(leg,'FontSize',18);
set(0,'DefaultAxesFontSize',18);
gname = [num2str(case_name),'__total_vs_new_deaths'];
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