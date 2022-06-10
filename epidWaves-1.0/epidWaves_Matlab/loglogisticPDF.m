% -----------------------------------------------------------------
%  loglogisticPDF.m
% -----------------------------------------------------------------
%  This routine defines the logistic function derivative.
%  
%  Input:
%  x   - independent variable
%  K   - curve's maximum value
%  r   - growth rate
%  tau - point of inflection
%  N   - number of waves
%  
%  Output:
%  dCdx - logistic function derivative value
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
function dCdx = loglogisticPDF(x,K,r,tau,N)

    % check number of arguments
    if nargin < 5
        error('Too few inputs.')
    elseif nargin > 5
        error('Too many inputs.')
    end
    
    % initialize dCdx value
	dCdx = 0.0;
    
    % multiples waves logistic function derivative
    for n=1:N
        
        % loglogistic function derivative
        dCdx = dCdx + log(r(n)) + log(K(n)) ...
            - r(n)*(x-tau(n)) -2*log(1+exp(-r(n)*(x-tau(n))));
    end
end
% -----------------------------------------------------------------