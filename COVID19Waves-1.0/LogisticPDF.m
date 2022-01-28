% -----------------------------------------------------------------
%  LogisticPDF.m
% -----------------------------------------------------------------
%  This routine defines the logistic function derivative.
%  
%  Input:
%  x   - independent variable
%  K   - curve's maximum value
%  r   - growth rate
%  tau - point of inflection
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
function dCdx = LogisticPDF(x,K,r,tau)

    % check number of arguments
    if nargin < 4
        error('Too few inputs.')
    elseif nargin > 4
        error('Too many inputs.')
    end
        
	% logistic function derivative
	dCdx = r.*K.*exp(-r.*(x-tau))./(1+exp(-r.*(x-tau))).^2;
end
% -----------------------------------------------------------------