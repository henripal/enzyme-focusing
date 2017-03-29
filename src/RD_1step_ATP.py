import numpy as np
import XixiLib as xl
from scipy.integrate import odeint


# get rid of this in final version
from imp import reload
reload(xl)
# mpl.style.use('classic')


def solve_system():

    Darray = np.array([70, 620, 70, 70])

    # kon and koff are dictionaries with kinetic constants
    # shorthand for a reaction is the list of RH terms in order
    # for example ESS is the reaction where EA + S -> ESA
    # EP_D is the reaction where EPD -> EP + D
    # units are normalized to microM

    # kon = {"ES": 2,
    #        "ESA": 1,
    #        "EDP": 800,
    #        "EA_S": 10,
    #        "EP_D": 400,
    #        "E_P": 6000,
    #        "ED": 40}

    # konarray = np.array([2, 1, 8000, 10, 400, 6000, 40])
    # trying it with lower glucose dissociation rates
    konarray = np.array([2,  8000*0])

    # koff = {"ES": 60,
    #         "ESA": 200,
    #         "ES_A": 0,
    #         "EDP": 800,
    #         "EA_S": 0,
    #         "EP_D": 0.2,
    #         "E_P": 2,
    #         "ED": 0}

    koffarray = np.array([60, 8000*0])

    ###############
    #    GRID     #
    ###############

    # defining the grid in x and t
    # x is in micrometers, t in s

    Nx = 200
    xMax = 400
    x = np.linspace(0, xMax, Nx)
    dx = x[1]-x[0]

    tSyringe = 0*60

    # t = np.arange(0, tSyringe, .001)

    t = np.array([0, tSyringe+10,
                  tSyringe+20, tSyringe+30, tSyringe+40,
                  tSyringe+50, tSyringe+60, tSyringe+70])

    # t = np.arange(tSyringe, tSyringe + 60, .1)

    ######################
    # INITIAL CONDITIONS #
    ######################
    # remember all units in microM
    # so 0.2 = 0.2 uM = 200 nM
    # 50e3 = 50e3 uM = 50 mM

    # EXPERIMENT 1
    # C0_E = 0.2*np.ones(Nx)*(np.tanh((x-dx*np.floor(Nx/3))/10) -
    #                         np.tanh((x-dx*np.floor(2*Nx/3))/10))/2
    # C0_S = 1*50e3 * np.ones(Nx)*(1-(np.tanh((x-dx*np.floor(Nx/3))/10) -
    #                                 np.tanh((x-dx*np.floor(2*Nx/3))/10))/2)
    # C0_A = 2e3 * np.ones(Nx)

    # # EXPERIMENT 2
    C0_E = 0.2*np.ones(Nx)
    C0_S = 50e3 * (np.tanh((x-dx*np.floor(Nx/3))/20) -
                   np.tanh((x-dx*np.floor(2*Nx/3))/20))/2

    C0_ES = np.zeros(Nx)
    C0_P = np.zeros(Nx)

    # concatenating all this in a 1D vector

    y0 = np.concatenate((C0_E, C0_S, C0_ES, C0_P), 0)

    ######################
    # CALLING THE SOLVER #
    ######################

    y = odeint(xl.ct_OneStep_ATP, y0, t,
               (Nx, dx, Darray, konarray, koffarray, tSyringe),
               atol=0.001, rtol=.0001)

    return y
