% -----------------------------------------------------------------
%  RegressionMC.m
% -----------------------------------------------------------------
%  This routine combines Monte Carlo simulation and a nonlinear 
%  regression algorithm to estimate an algebraic statistical model
%  that fits a given dataset. Monte Carlo is employed to reduced
%  the result sensibility to the choice of an inital guess in the
%  regression process.
%  
%  Input:
%  xdata      - independent parameter data
%  ydata      - dependent   parameter data
%  MyModel    - algebraic model structure
%  ModelParam - algebraic model parameters
%  
%  Output:
%  MyFit    - fitting model object
%  ErrorObj - fitting error object
%  
%  Reference:
%  PRL Gianfelice, RS Oyarzabal, A Cunha Jr,
%  JMV Grzybowski, FC Batista, EEN Macau
%  The starting dates of COVID-19 multiple waves,
%  Preprint, 2022
% -----------------------------------------------------------------
%  programmer: Americo Cunha Jr (UERJ)
%               
%  last update: Jan 14, 2022
% -----------------------------------------------------------------

% -----------------------------------------------------------------
function [MyFit,ErrorObj] = RegressionMC(xdata,ydata,MyModel,HyperParam)

    % check number of arguments
    if nargin < 4
        error('Too few inputs.')
    elseif nargin > 4
        error('Too many inputs.')
    end
    
    % check for consistency in xdata and ydata
    if size(xdata) ~= size(ydata)
        error('xdata and ydata must have the same dimensions')
    end
    
    % range of admissible values for model parameters
    lb = HyperParam.lb;
	ub = HyperParam.ub;
    
    % check for consistency in lb and ub
    if size(lb) ~= size(ub)
        error('lb and ub must have the same dimensions')
    end
    if isempty(lb) || isempty(ub)
        error('lb and ub must be non-empty')
    end
    if any(isnan(lb)) || any(~isfinite(lb))
        error('lb must have finite values')
    end
    if any(isnan(ub)) || any(~isfinite(ub))
        error('ub must have finite values')
    end
    if any(lb > ub)
        error('lb < ub for all components')
    end
    
    % convert lb to a column vector (if necessary)
    [s1,s2] = size(lb);
    if isrow(lb)
        lb = lb(:);
    elseif s2 ~= 1
        error('lb must be a column vector')
    end

    % convert ub to a column vector (if necessary)
    [s1,s2] = size(ub);
    if isrow(ub)
        ub = ub(:);
    elseif s2 ~= 1 
        error('ub must be a column vector')
    end
    
    % number of samples for Monte Carlo simulation
    Ns = HyperParam.Ns;
    if Ns <= 0
        error('Ns must be a positive integer')
    end

    % estimation for squared sum of errors
    ErrorObj.rmse = inf;

    % ensemble of random initial guesses
    x0 =  lb + (ub-lb)*rand(1,Ns);
    
    % find the best curve fit via Monte Carlo
    for n = 1:Ns

        % n-th curve fit
        [curvefit,gof] = fit(xdata,ydata,MyModel,...
                               'Lower',lb,...
                               'Upper',ub,...
                               'start',x0(:,n));

        % update the model with small error
        if gof.rmse < ErrorObj.rmse
            MyFit            = curvefit;
            ErrorObj.sse     = gof.sse;
            ErrorObj.rsquare = gof.rsquare;
            ErrorObj.rmse    = gof.rmse;
        end
    end
end
% -----------------------------------------------------------------