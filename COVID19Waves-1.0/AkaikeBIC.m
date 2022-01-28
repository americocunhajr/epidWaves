% -----------------------------------------------------------------
%  AkaikeBIC.m
% -----------------------------------------------------------------
%  This function computes Akaike and Bayesian information 
%  criteria scores given the dataset and a statistical model.
%  
%  Input:
%  ModelPDF - statistical model PDF
%  data     - dataset
%  p0       - initial guess for the MLE estimator
%  
%  Output:
%  AIC - Akaike score
%  BIC - BIC score
%  
%  Reference:
%  PRL Gianfelice, RS Oyarzabal, A Cunha Jr, JMV Grzybowski,
%  FC Batista, EEN Macau
%  The starting dates of COVID-19 multiple waves,
%  Preprint, 2022
% -----------------------------------------------------------------
%  programmer: Americo Cunha Jr (UERJ)
%               
%  last update: Jan 14, 2022
% -----------------------------------------------------------------

% -----------------------------------------------------------------
function [AIC,BIC] = AkaikeBIC(ModelPDF,data,p0)

    % check number of arguments
    if nargin < 3
        error('Too few inputs.')
    elseif nargin > 3
        error('Too many inputs.')
    end

	% number of data points
    n = length(data);
    
    % number of parameters
    k = length(p0);
    
    % log-likelihood function
    LogLikelihood = @(p) sum(log(abs(ModelPDF(data,p))));
    
    % solver optionals
    opt = optimoptions(@fsolve,...
                       'Algorithm','levenberg-marquardt',...
                       'Display','off');

    % maximum likelihood estimator
    p_mle = fsolve(LogLikelihood,p0,opt);
    
    %  maximum log-likelihood value
    LL = LogLikelihood(p_mle);
    
    % Akaike information criterion
    AIC = 2*k - 2*LL;
    
	% Bayesian information criterion
    BIC = k*log(n) - 2*LL;
end
% -----------------------------------------------------------------