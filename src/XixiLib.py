from numba import jit
import numpy as np


# jit to compile with numba
@jit
def ct_OneStep_ATP(y, t, Nx, dx, D_c, kon, koff, tSyringe):
    # this function is the RHS of the RD equation
    # y: start concentration 1D vector. will be split

    # This has XD of ES toward A and E toward S
    # Ke1 is the equilibirum connstant for ES and A
    Ke1 = kon[1]/koff[1]

    # Ke is the equilibirum constant for ES and S
    Ke0 = kon[0]/koff[0]

    if t < tSyringe:
        D_c_temp = D_c*0
    else:
        D_c_temp = D_c

    # splitting y
    [E, S, ES, P] = np.reshape(y, (4, Nx))

    return np.concatenate(subloop(E, S, ES, P,
                                  kon, koff, dx, Ke0, Ke1, Nx, D_c_temp), 0)


@jit(nopython=True, cache=True)
def subloop(E, S, ES, P,
            kon, koff, dx, Ke0, Ke1, Nx, D_c):
    # subloop is the nopython-friendly part of the RHS
    # we start by initializing the arrays:
    Et = np.zeros(Nx)
    St = np.zeros(Nx)
    ESt = np.zeros(Nx)
    Pt = np.zeros(Nx)

    # then loop over the x dimension to discretize
    for i in np.arange(0, Nx):
        if i == 0:
            # boundary condition at x=0
            # the derivatives are 0 so we simplfy expressions of the
            # second derivatives
            Exx = (E[i+1]-2*E[i]+E[i+1])/dx**2
            Sxx = (S[i+1]-2*S[i]+S[i+1])/dx**2
            ESxx = (ES[i+1]-2*ES[i]+ES[i+1])/dx**2
            Pxx = (P[i+1]-2*P[i]+P[i+1])/dx**2
            # factorES = 0
            ESx = 0
            Ex = 0
            Sx = 0
        elif i == Nx-1:
            # boundary condition at x=Nx
            # same simplification as x=0
            Exx = (E[i-1]-2*E[i]+E[i-1])/dx**2
            Sxx = (S[i-1]-2*S[i]+S[i-1])/dx**2
            ESxx = (ES[i-1]-2*ES[i]+ES[i-1])/dx**2
            Pxx = (P[i-1]-2*P[i]+P[i-1])/dx**2
            # factorES = 0
            ESx = 0
            Ex = 0
            Sx = 0
        else:
            # traditional discretization of second derivatives
            Exx = (E[i+1]-2*E[i]+E[i-1])/dx**2
            Sxx = (S[i+1]-2*S[i]+S[i-1])/dx**2
            ESxx = (ES[i+1]-2*ES[i]+ES[i-1])/dx**2
            Pxx = (P[i+1]-2*P[i]+P[i-1])/dx**2
        # factorES = 1/(2*dx) * (ES[i+1]/(1+Ke*A[i+1])-ES[i-1]/(1+Ke*A[i-1]))
            ESx = (ES[i+1]-ES[i-1])/dx/2
            Ex = (E[i+1]-E[i-1])/dx/2
            Sx = (S[i+1]-S[i-1])/dx/2

        # once all derivatives have been completed we can explicit the
        # time derivatives of each species

        Et[i] = D_c[0]*Exx-kon[0]*E[i]*S[i]+koff[0]*ES[i]+ \
            - koff[1]*E[i]*P[i] \
            + kon[1]*ES[i] \
            - 1*D_c[0]*Ke0 \
            * (Ex*Sx-Ke0*E[i]*(Sx**2)/(1+Ke0*S[i])+E[i]*Sxx)/(1+Ke0*S[i])
        # cross diffusion term above ^^^^^^
        # cross diffusion of E toward S

        St[i] = D_c[1]*Sxx-kon[0]*E[i]*S[i]+koff[0]*ES[i]

        ESt[i] = D_c[2]*ESxx + kon[0]*E[i]*S[i]-koff[0]*ES[i] \
            - kon[1]*ES[i] + koff[1]*E[i]*P[i]

        Pt[i] = D_c[3]*Pxx+kon[1]*ES[i]-koff[1]*E[i]*P[i]

    return (Et, St, ESt, Pt)
