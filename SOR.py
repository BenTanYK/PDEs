"""
Programme to calculate the electric and potential field for a
positive point charge.
"""

import matplotlib
matplotlib.use('TKAgg')

import sys
import math
import random
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import time
import pandas as pd

lx=int(sys.argv[1])
convergence=float(sys.argv[2])

from poisson import point_charge, error, E_field, E_pot, jacobi, gauss_seidel, sor

def main():

    omega_vals=0.01*np.arange(100, 200)
    n_conv=[]
    rho=point_charge(lx)

    for m in range(len(omega_vals)):

        phi=np.zeros((lx,lx,lx),dtype=float) 

        err=10
        n=0

        while err>convergence:

            n+=1

            omega=omega_vals[m]
            phi_new=sor(np.copy(phi), rho, omega)

            err=error(phi, phi_new)
            phi=phi_new

        n_conv.append(n)
        print('Omega=' + str(omega))

    df=pd.DataFrame()
    df['omega']=omega_vals
    df['Updates to convergence']=np.array(n_conv)
    df.to_excel('SOR.xlsx')

main()



