% -----------------------------------------------------------------
%  Main_COVID19_HypothesisTest_1st_wave_RJ.m
% -----------------------------------------------------------------
%  This program run a hypothesis test.
%  
%  Reference:
%  PRL Gianfelice, RS Oyarzabal, A Cunha Jr,
%  JMV Grzybowski, FC Batista, EEN Macau
%  The starting dates of COVID-19 multiple waves
%  Preprint, 2022
% -----------------------------------------------------------------
%  programmer: Americo Cunha Jr (UERJ)
%               
%  last update: Jan 22, 2021
% -----------------------------------------------------------------


clc
clear
close all

% program header
% -----------------------------------------------------------
disp(' ---------------------------------------------------------- ')
disp(' Estimation of 1sdt wave starting date                      ')
disp(' (Hypothesis Test)                                          ')
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



% load surveillance data
% -----------------------------------------------------------
tic
disp(' '); 
disp(' --- Hypothesis Test --- ');
disp(' ');
disp('    ... ');
disp(' ');

% Jan 1, 2020 -   1
% Fev 1, 2020 -  32
% Mar 1, 2020 -  61

% null-hypothesis (starting date <= March 1, 2020)

ref_date   = 61;
%start_date = [55 56 60 63 68 76];
start_date = [61 60 60 60];

% hypothesis Test 
[h,p,ci,stats] = ttest(start_date,ref_date,'Tail','left')
%[h,p,ci,stats] = ttest(start_date,ref_date,'Tail','right')


