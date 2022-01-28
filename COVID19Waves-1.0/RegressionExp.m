% -----------------------------------------------------------------
%  RegressionExp.m
% -----------------------------------------------------------------
%  This routine computes an exponential regression y = B*exp(r*t).
%  
%  Input:
%  xdata      - independent parameter data
%  ydata      - dependent   parameter data
%  
%  Output:
%  r - exponential rate
%  A - exponential amplitude
%  
%  Reference:
%  PRL Gianfelice, RS Oyarzabal, A Cunha Jr,
%  JMV Grzybowski, FC Batista, EEN Macau
%  The starting dates of COVID-19 multiple waves,
%  Preprint, 2022
% -----------------------------------------------------------------
%  programmer: Americo Cunha Jr (UERJ)
%               
%  last update: Jan 17, 2022
% -----------------------------------------------------------------

% -----------------------------------------------------------------
function [r,B] = RegressionExp(xdata,ydata)

    % check number of arguments
    if nargin < 2
        error('Too few inputs.')
    elseif nargin > 2
        error('Too many inputs.')
    end
    
    % check for consistency in xdata and ydata
    if size(xdata) ~= size(ydata)
        error('xdata and ydata must have the same dimensions')
    end
    
    % convert xdata to a column vector (if necessary)
    [s1,s2] = size(xdata);
    if isrow(xdata)
        xdata = xdata(:);
    elseif s2 ~= 1
        error('xdata must be a column vector')
    end
    
    % convert ydata to a column vector (if necessary)
    [s1,s2] = size(ydata);
    if isrow(ydata)
        ydata = ydata(:);
    elseif s2 ~= 1
        error('xdata must be a column vector')
    end
    
    % number of data points
    N = length(xdata);
    
    % sample matrix
    A = [xdata ones(N,1)];
    
    % right hand side vector
    b = log(ydata);
    
    % regression coefficients
    x = A\b;
    
    % model parameters
    r = x(1);
    B = exp(x(2));
end
% -----------------------------------------------------------------