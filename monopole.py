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

cmap = plt.colormaps['plasma']

lx=int(sys.argv[1])
convergence=float(sys.argv[2])
algorithm=str(sys.argv[3])

from poisson import point_charge, error, E_field, E_pot, jacobi, gauss_seidel

def main():

    rho=point_charge(lx) #initialise charge distribution (point charge)

    phi=np.zeros((lx,lx,lx),dtype=float) 

    err=10
    n=0

    while err>convergence:

        n+=1
        if n%10==0: 
            print(n)
            print(err)

        if algorithm=='jacobi' or algorithm=='Jacobi':
            phi_new=jacobi(phi,rho)

        elif algorithm=='gauss_seidel' or algorithm=='Gauss_seidel':
            phi_new=gauss_seidel(np.copy(phi),rho)

        else:
            print('select jacobi or gauss_seidel algorithm')

        err=error(phi, phi_new)
        phi=phi_new

    plt.cla()
    im=plt.imshow(phi[int(lx/2)], cmap=cmap, animated=True)
    plt.colorbar()
    plt.draw()
    plt.pause(5.0)
    # plt.savefig('pot_cut.png')
    plt.clf()

    # data_cut=phi[int(lx/2)]
    # df_cut=pd.DataFrame(data_cut)
    # df_cut.to_excel('potential_field.xlsx')

    # data=E_pot(phi)
    # columns=np.array(['x', 'y', 'z', 'r', 'pot', 'Ex', 'Ey', 'Ez', '|E|'])
    # df=pd.DataFrame(np.transpose(data), index=None, columns=columns)
    # df.to_excel('poisson_data.xlsx')

start_time = time.time() 
main()
print('Programme run time is:' + "%s seconds" % (time.time() - start_time))


