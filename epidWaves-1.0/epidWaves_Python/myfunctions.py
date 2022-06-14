# -*- coding: utf-8  -*-
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime
import matplotlib.dates as mdates
from lmfit import Parameters
import sklearn.metrics as metrics
import sympy as sym
from sympy import symbols, lambdify, hessian, Matrix, ordered
from scipy import optimize
from pandas.plotting import deregister_matplotlib_converters


def load_data(file_name):
    custom_date_parser = lambda x: datetime.strptime(x, "%m/%d/%y")
    df = pd.read_csv(file_name,
                     parse_dates=['data'],
                     date_parser=custom_date_parser)
    df['time'] = np.arange(1, len(df['data']) + 1, 1)
    df['data'] = pd.to_datetime(df['data'], format='%Y-%m-%d')
    df['date'] = df['data']
    df = df.set_index('data')
    return df


def LogisticPDF(x, K, r, tau):
    #  This routine defines the logistic function derivative.
    #  Input:
    #  x   - independent variable
    #  K   - curve's maximum value
    #  r   - growth rate
    #  tau - point of inflection
    #  Output:
    #  logistic function derivative value
    return r * K * np.exp(-r * (x - tau)) / (1 + np.exp(-r * (x - tau))) ** 2


def LogisticPDF6w(x, K1, K2, K3, K4, K5, K6, r1, r2, r3, r4, r5, r6, tau1, tau2, tau3, tau4, tau5, tau6):
    #  This routine defines the logistic function derivative.
    #  Input:
    #  x   - independent variable
    #  K   - curve's maximum value
    #  r   - growth rate
    #  tau - point of inflection
    #  Output:
    #  logistic function derivative value
    dCdx  = r1 * K1 * np.exp(-r1 * (x - tau1)) / (1 + np.exp(-r1 * (x - tau1))) ** 2
    dCdx += r2 * K2 * np.exp(-r2 * (x - tau2)) / (1 + np.exp(-r2 * (x - tau2))) ** 2
    dCdx += r3 * K3 * np.exp(-r3 * (x - tau3)) / (1 + np.exp(-r3 * (x - tau3))) ** 2
    dCdx += r4 * K4 * np.exp(-r4 * (x - tau4)) / (1 + np.exp(-r4 * (x - tau4))) ** 2
    dCdx += r5 * K5 * np.exp(-r5 * (x - tau5)) / (1 + np.exp(-r5 * (x - tau5))) ** 2
    dCdx += r6 * K6 * np.exp(-r6 * (x - tau6)) / (1 + np.exp(-r6 * (x - tau6))) ** 2
    return dCdx


def LogisticCDF(x, K, r, tau):
    #  This routine defines the logistic function.
    #  Input:
    #  x   - independent variable
    #  K   - curve's maximum value
    #  r   - growth rate
    #  tau - point of inflection
    #  Output:
    #  logistic function value
    return K / (1 + np.exp(-r * (x - tau)))


def RegressionMC(xdata, ydata, MyModel, HyperParam):
    #  This routine combines Monte Carlo simulation and a nonlinear
    #  regression algorithm to estimate an algebraic statistical model
    #  that fits a given dataset. Monte Carlo is employed to reduced
    #  the result sensibility to the choice of an inital guess in the
    #  regression process.
    #  Input:
    #  xdata      - independent parameter data
    #  ydata      - dependent   parameter data
    #  MyModel    - algebraic model structure
    #  HyperParam - algebraic model parameters
    #  Output:
    #  Result_last - fitting model object
    #  ErrorObj - fitting error object

    # range of admissible values for model parameters

    ub = HyperParam['ub']
    lb = HyperParam['lb']

    # number of samples for Monte Carlo simulation
    Ns = HyperParam['Ns']

    ErrorObj = {}
    # estimation for squared sum of errors
    ErrorObj['rmse'] = 10 ** 10

    # ensemble of random initial guesses
    x0 = lb + (ub - lb) * np.random.rand(1, Ns)

    # find the best curve fit via Monte Carlo
    for n in range(0, Ns):
        # n-th curve fit
        print('Monte Carlo: ' + str(n+1))
        params = Parameters()
        #add with tuples:(NAME VALUE        VARY     MIN        MAX    EXPR  BRUTE_STEP)
        if len(ub) == 2:
            tau = HyperParam['tau']
            params.add_many(('tau',      tau,  False,    None,      None),
                            ('K'  , x0[0, n],   True, lb[0, 0], ub[0, 0]),
                            ('r'  , x0[1, n],   True, lb[1, 0], ub[1, 0]))
        else:
            pp = HyperParam['p']
            for i in range(0, len(pp)):
                params.add_many((pp[i], x0[i, n], True, lb[i, 0], ub[i, 0]))

        Result = MyModel.fit(ydata, params, x=xdata)
        yhat = Result.best_fit
        mse = metrics.mean_squared_error(ydata, yhat)
        rmse = np.sqrt(mse)
        rsquare = metrics.r2_score(ydata, yhat)

        # update the model with small error
        if rmse < ErrorObj['rmse']:
            ErrorObj['mse'] = mse
            ErrorObj['rmse'] = rmse
            ErrorObj['rsquare'] = rsquare
            Result_last = Result

    return Result_last, ErrorObj


def predband(x, xd, yd, p, func, conf=0.95):
    from scipy import stats
    x, xd, yd = np.array(x), np.array(xd), np.array(yd)
    alpha = 1.0 - conf  # significance
    N = xd.size  # data sample size
    var_n = len(p)  # number of parameters
    # Quantile of Student's t distribution for p=(1-alpha/2)
    q = stats.t.ppf(1.0 - alpha / 2.0, N - var_n)  # Student's t
    #q = stats.norm.ppf(1.0 - alpha / 2.0, loc=0, scale=1)
    se = np.sqrt(1. / (N - var_n) * np.sum((yd - func(xd, *p)) ** 2))
    # print(se)
    # print(yd.std(ddof=1))
    # print(stats.sem(yd))
    sx = (x - xd.mean()) ** 2
    sxd = np.sum((xd - xd.mean()) ** 2)
    yp = func(x, *p)
    dy = q * se * np.sqrt(1.0 + (1.0 / N) + (sx / sxd))
    # dy = q * se * np.sqrt(1.0 + (1.0/N) )
    lpb, upb = yp - dy, yp + dy
    return lpb, upb


# -------------------------------------
def LogisticPDF_model(tau):
    K = symbols('K', real=True)
    r = symbols('r', real=True)
    t = symbols('t', real=True)
    return r * K * sym.exp(-r * (t - tau)) / (1 + sym.exp(-r * (t - tau))) ** 2


def AkaikeBIC(ModelPDF, time_train, p0):
    t = symbols('t', real=True)
    np.seterr(divide = 'ignore')

    def LogLikelihood_Model_sympy(Model, t, time):
        LogLikelihood_sympy = sym.log(Model.subs({t: time[0]}))
        for id in range(1, len(time)):
            LogLikelihood_sympy = LogLikelihood_sympy + sym.log(Model.subs({t: time[id]}))
        return -LogLikelihood_sympy

    def LogLikelihood_Model_numpy(Model, t, time):
        LogLikelihood_sympy = LogLikelihood_Model_sympy(Model, t, time)
        v = list(ordered(LogLikelihood_sympy.free_symbols))
        return lambdify(v, LogLikelihood_sympy, 'numpy')

    def LogLikelihood_numpy(x):
        return LogLikelihood_numpy2(x[0], x[1]).astype(object)

    def gradient_hessian_sympy(f_sympy):
        v = list(ordered(f_sympy.free_symbols))
        gradient = lambda f_sympy, v: Matrix([f_sympy]).jacobian(v)
        return gradient(f_sympy, v), hessian(f_sympy, v)

    def gradient_hessian_numpy(f_sympy):
        v = list(ordered(f_sympy.free_symbols))
        [grad_sympy, hess_sympy] = gradient_hessian_sympy(f_sympy)
        return lambdify(v, grad_sympy, 'numpy'), lambdify(v, hess_sympy, 'numpy')

    def grad_numpy(x):
        return grad_numpy2(x[0], x[1]).astype(object)

    LogLikelihood_sympy = LogLikelihood_Model_sympy(ModelPDF, t, time_train)
    LogLikelihood_numpy2 = LogLikelihood_Model_numpy(ModelPDF, t, time_train)
    [grad_numpy2, hess_numpy2] = gradient_hessian_numpy(LogLikelihood_sympy)

    # ------ method='SLSQP' ------------------------------
    minimum = optimize.minimize(LogLikelihood_numpy, p0,
                                bounds=[[0, 1e7], [0, 1]],
                                method='SLSQP',
                                jac=grad_numpy,
                                #options={'maxiter':400, 'disp':True},
                                options={'maxiter': 400},
                                tol=1e-4)
    # maximum log-likelihood value
    LL = LogLikelihood_numpy(minimum.x)

    # Akaike information criterion
    AIC = 2 * len(p0) - 2 * LL

    # Bayesian information criterion
    BIC = len(p0) * np.log(len(time_train)) - 2 * LL

    return AIC, BIC


def set_xmargin(ax, left=0.0, right=0.3):
    ax.set_xmargin(0)
    ax.autoscale_view()
    lim = ax.get_xlim()
    delta = np.diff(lim)
    left = lim[0] - delta * left
    right = lim[1] + delta * right
    ax.set_xlim(left, right)


def graph_C_1w(date, yraw, yMA, yupp, ylow, ypred, graphobj):
    deregister_matplotlib_converters()
    fig, host = plt.subplots(figsize=(6, 6), dpi=100, facecolor='w', edgecolor='k')
    plt.subplots_adjust(wspace=0, hspace=0, left=0.10, right=0.99, top=0.96, bottom=0.09)
    host.grid(True, which="both", ls="-")
    host.scatter(date, yraw, color="m", s=12, label=graphobj['leg1'])
    host.plot(date, yraw, color="m", linewidth=1)
    host.plot(date, yMA, color="g", linewidth=2, label=graphobj['leg2'])
    host.plot(date, ypred, color="b", linewidth=2, label=graphobj['leg3'])
    host.fill_between(date, yupp, ylow, alpha=0.3, color='gray', label=graphobj['leg4'])
    host.set_ylim((graphobj['ymin'], graphobj['ymax']))
    host.set_yticks(np.arange(graphobj['ymin'], graphobj['ymax'] + .1, graphobj['ymax'] / 10))
    plt.ylabel(graphobj['ylab'], fontsize=14)
    host.xaxis.set_major_formatter(mdates.DateFormatter('%b%Y'))
    host.xaxis.set_major_locator(mdates.MonthLocator(interval=1))
    set_xmargin(host, left=0, right=0)
    legend = host.legend(fontsize='large', loc='upper left',
                         borderpad=0.5, labelspacing=.5,
                         framealpha=1, facecolor='white',
                         title='',
                         frameon=True)
    legend.get_title().set_fontsize('10')
    legend._legend_box.align = "left"
    plt.tight_layout()
    fig.savefig(graphobj['gname'] + "_py.png")
    plt.show()
    return


def graph_I_1w(date, yraw, yMA, yupp, ylow, ypred, graphobj):
    deregister_matplotlib_converters()
    fig, host = plt.subplots(figsize=(6, 6), dpi=100, facecolor='w', edgecolor='k')
    plt.subplots_adjust(wspace=0, hspace=0, left=0.10, right=0.99, top=0.96, bottom=0.09)
    host.grid(True, which="both", ls="-")
    host.scatter(date, yraw, color="m", s=12, label=graphobj['leg1'])
    host.plot(date, yraw, color="m", linewidth=1)
    host.plot(date, yMA, color="g", linewidth=2, label=graphobj['leg2'])
    host.plot(date, ypred, color="b", linewidth=2, label=graphobj['leg3'])
    host.fill_between(date, yupp, ylow, alpha=0.3, color='gray', label=graphobj['leg4'])
    host.plot([date[graphobj['taux']], date[graphobj['taux']]],
              [0, .85 * graphobj['tauy'] * graphobj['ymax']], color="k", ls="--", linewidth=2)
    text_date = datetime.strptime(str(graphobj['taux']), "%Y-%m-%d %H:%M:%S")
    host.text(date[graphobj['taux']], .85 * graphobj['tauy'] * graphobj['ymax'],
              str(datetime.strftime(text_date, "%d-%b-%Y")), horizontalalignment=graphobj['tauh'], fontsize=12,
              color='k',
              bbox=dict(facecolor='white', edgecolor='none', pad=2))
    host.set_ylim((graphobj['ymin'], graphobj['ymax']))
    plt.ylabel(graphobj['ylab'], fontsize=14)
    host.xaxis.set_major_formatter(mdates.DateFormatter('%b%Y'))
    host.xaxis.set_major_locator(mdates.MonthLocator(interval=1))
    set_xmargin(host, left=0, right=0)
    legend = host.legend(fontsize='large', loc='upper left',
                         borderpad=0.5, labelspacing=.5,
                         framealpha=1, facecolor='white',
                         title='',
                         frameon=True)
    legend.get_title().set_fontsize('10')
    legend._legend_box.align = "left"
    plt.tight_layout()
    fig.savefig(graphobj['gname'] + "_py.png")
    plt.show()
    return


def graph_I_6w(date, yraw, yMA, yupp, ylow, ypred, tau_ast, graphobj):
    deregister_matplotlib_converters()
    fig, host = plt.subplots(figsize=(6, 6), dpi=100, facecolor='w', edgecolor='k')
    plt.subplots_adjust(wspace=0, hspace=0, left=0.10, right=0.99, top=0.96, bottom=0.09)
    host.grid(True, which="both", ls="-")
    host.scatter(date, yraw, color="m", s=12, label=graphobj['leg1'])
    host.plot(date, yraw, color="m", linewidth=1)
    host.plot(date, yMA, color="g", linewidth=2, label=graphobj['leg2'])
    host.plot(date, ypred, color="b", linewidth=2, label=graphobj['leg3'])
    host.fill_between(date, yupp, ylow, alpha=0.3, color='gray', label=graphobj['leg4'])
    for i in range(0, len(tau_ast)):
        host.plot([date[tau_ast[i]], date[tau_ast[i]]],
              [0, .85 * graphobj['tauy'] * graphobj['ymax']], color="k", ls="--", linewidth=1.5)
    host.set_ylim((graphobj['ymin'], graphobj['ymax']))
    plt.ylabel(graphobj['ylab'], fontsize=14)
    host.xaxis.set_major_formatter(mdates.DateFormatter('%b%Y'))
    host.xaxis.set_major_locator(mdates.MonthLocator(interval=4))
    set_xmargin(host, left=0, right=0)
    legend = host.legend(fontsize='large', loc='upper right',
                         borderpad=0.5, labelspacing=.5,
                         framealpha=1, facecolor='white',
                         title='',
                         frameon=True)
    legend.get_title().set_fontsize('10')
    legend._legend_box.align = "left"
    plt.tight_layout()
    fig.savefig(graphobj['gname'] + "_py.png")
    plt.show()
    return

def graph_I_raw(date, yraw, yMA, graphobj):
    deregister_matplotlib_converters()
    fig, host = plt.subplots(figsize=(6, 6), dpi=100, facecolor='w', edgecolor='k')
    plt.subplots_adjust(wspace=0, hspace=0, left=0.10, right=0.99, top=0.96, bottom=0.09)
    host.grid(True, which="both", ls="-")
    host.scatter(date, yraw, color="w", edgecolor="m", s=14, label=graphobj['leg1'])
    #host.plot(date, yraw, color="m", linewidth=1)
    host.plot(date, yMA, color="g", linewidth=2, label=graphobj['leg2'])
    host.set_ylim((graphobj['ymin'], graphobj['ymax']))
    plt.ylabel(graphobj['ylab'], fontsize=14)
    host.xaxis.set_major_formatter(mdates.DateFormatter('%b%Y'))
    #host.xaxis.set_major_locator(mdates.DayLocator(bymonthday=range(1, 2), interval=1))
    host.xaxis.set_major_locator(mdates.MonthLocator(interval=4))
    set_xmargin(host, left=0, right=0)
    legend = host.legend(fontsize='large', loc='upper left',
                         borderpad=0.5, labelspacing=.5,
                         framealpha=1, facecolor='white',
                         title='',
                         frameon=True)
    legend.get_title().set_fontsize('10')
    legend._legend_box.align = "left"
    plt.tight_layout()
    fig.savefig(graphobj['gname'] + "_py.png")
    plt.show()
    return

def graph_I_vs_C_raw(x1, y1, x2, y2, graphobj):
    deregister_matplotlib_converters()
    fig, host = plt.subplots(figsize=(6, 6), dpi=100, facecolor='w', edgecolor='k')
    plt.subplots_adjust(wspace=0, hspace=0, left=0.10, right=0.99, top=0.96, bottom=0.09)
    host.grid(True, which="both", ls="-")
    host.scatter(x1, y1, color="w", edgecolor="m", s=14, label=graphobj['leg1'])
    host.plot(x2, y2, color="g", linewidth=2, label=graphobj['leg2'])
    host.set_xlim((graphobj['xmin'], graphobj['xmax']))
    host.set_ylim((graphobj['ymin'], graphobj['ymax']))
    set_xmargin(host, left=0, right=0)
    plt.xlabel(graphobj['xlab'], fontsize=14)
    plt.ylabel(graphobj['ylab'], fontsize=14)
    legend = host.legend(fontsize='large', loc='lower center',
                         borderpad=0.5, labelspacing=.5,
                         framealpha=1, facecolor='white',
                         title='',
                         frameon=True)
    legend.get_title().set_fontsize('10')
    legend._legend_box.align = "left"
    host.set_xscale('log')
    host.set_yscale('log')
    plt.tight_layout()
    fig.savefig(graphobj['gname'] + "_py.png")
    plt.show()
    return
