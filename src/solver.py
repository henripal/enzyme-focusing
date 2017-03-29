import numpy as np
import XixiLib as xl
from scipy.integrate import odeint
from imp import reload
reload(xl)


def solve_system(Nx, xMax, t, tSyringe, enzyme_conc, substrate_conc,
                 konarray, koffarray):

    # CONSTANTS #
    # D = {"E": 70 "S": 620"ES": 70 "P": 70
    Darray = np.array([70, 620, 70, 70])

    #    GRID     #
    # x is in micrometers, t in s
    x = np.linspace(0, xMax, Nx)
    dx = x[1]-x[0]
    Nt = t.size

    # INITIAL CONDITIONS #
    # so 0.2 = 0.2 uM = 200 nM
    # 50e3 = 50e3 uM = 50 mM
    C0_E = enzyme_conc*np.ones(Nx)
    C0_S = substrate_conc * (np.tanh((x-dx*np.floor(Nx/3))/20) -
                   np.tanh((x-dx*np.floor(2*Nx/3))/20))/2
    # C0_S = 50e3 * np.ones(Nx)
    C0_ES = np.zeros(Nx)
    C0_P = np.zeros(Nx)

    # concatenating all this in a 1D vector
    y0 = np.concatenate((C0_E, C0_S, C0_ES, C0_P), 0)

    # CALLING THE SOLVER #
    y = odeint(xl.ct_OneStep_ATP, y0, t,
               (Nx, dx, Darray, konarray, koffarray, tSyringe),
               atol=0.001, rtol=.0001)

    return np.reshape(y, (Nt, 4, Nx))


def tochannel(y):
    """ integrates the result of solver to three channels
    Returns:
    the aggregated 3 dimensional tensor corresponding to
    the total concentration of product per channel.
    dim1: time dim2: product dim3: channel
    """
    Nx = y.shape[2]
    Nt = y.shape[1]
    Np = y.shape[0]
    lim_inf = np.floor(Nx/3)
    lim_sup = np.floor(2*Nx/3)
    yout = np.zeros((Np, Nt, 3))

    for i, time in enumerate(y):
        for j, product in enumerate(time):
            yout[i, j, :] = [product[0:lim_inf].mean(),
                             product[lim_inf+1:lim_sup].mean(),
                             product[lim_sup+1:Nx].mean()]

    return yout


def tototal(y):
    """ integrates the result of solver to total concentrations
    Returns:
    the aggregated 2 dimensional matrix corresponding to
    the total concentration of product.
    dim1: time dim2: product
    """
    Nt = y.shape[1]
    Np = y.shape[0]
    yout = np.zeros((Np, Nt))

    for i, time in enumerate(y):
        for j, product in enumerate(time):
            yout[i, j] = product.mean()

    return yout
