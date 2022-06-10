% -----------------------------------------------------------------
%  LogisticCDF6w.m
% -----------------------------------------------------------------
%  This routine defines the logistic function.
%  
%  Input:
%  x   - independent variable
%  K   - curve's maximum value
%  r   - growth rate
%  tau - point of inflection
%  N   - number of waves
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
function C = LogisticCDF6w(x,K1,K2,K3,K4,K5,K6,...
                             r1,r2,r3,r4,r5,r6,...
                             tau1,tau2,tau3,tau4,tau5,tau6)

    % check number of arguments
    if nargin < 19
        error('Too few inputs.')
    elseif nargin > 19
        error('Too many inputs.')
    end
    
    % multiples waves logistic function
	C = K1./(1+exp(-r1.*(x-tau1))) + ...
        K2./(1+exp(-r2.*(x-tau2))) + ...
        K3./(1+exp(-r3.*(x-tau3))) + ...
        K4./(1+exp(-r4.*(x-tau4))) + ...
        K5./(1+exp(-r5.*(x-tau5))) + ...
        K6./(1+exp(-r6.*(x-tau6)));
end
% -----------------------------------------------------------------