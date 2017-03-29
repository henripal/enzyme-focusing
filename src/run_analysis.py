import numpy as np
import solver as slv
from imp import reload
reload(slv)

tSyringe = 0
Nx = 100
xMax = 400
# t = np.arange(0, 1, 0.1)
# t = np.append(t1, np.arange(.1, 50, .1))  
t = np.arange(0, 50, 1)

Nt = t.size

konarray = np.array([2,  8000*1])

koffarray = np.array([60, 8000*0])
enzyme_conc = 0.2
substrate_conc = 50e3

y = slv.solve_system(Nx, xMax, t, tSyringe,
                     enzyme_conc, substrate_conc, konarray, koffarray)


enzyme_conc = enzyme_conc*1.01 

y2 = slv.solve_system(Nx, xMax, t, tSyringe,
                      enzyme_conc, substrate_conc, konarray, koffarray)

# this will return a 3 dimensional tensor with
# i is the time
# j is the species
# k is the space
