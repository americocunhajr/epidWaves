% -----------------------------------------------------------------
%  LogisticCDF.m
% -----------------------------------------------------------------
%  This routine defines the logistic function.
%  
%  Input:
%  x   - independent variable
%  K   - curve's maximum value
%  r   - growth rate
%  tau - point of inflection
%  
%  Output:
%  C - logistic function value
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
function C = LogisticCDF(x,K,r,tau)

    % check number of arguments
    if nargin < 4
        error('Too few inputs.')
    elseif nargin > 4
        error('Too many inputs.')
    end
    
    % logistic function
	C = K./(1+exp(-r.*(x-tau)));
end
% -----------------------------------------------------------------