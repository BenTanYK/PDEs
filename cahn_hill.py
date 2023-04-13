"""
Functions for use in PDE simulations
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

cmap = plt.colormaps['coolwarm']

#numerical constants
a=0.1 
b=0.1
kappa=0.1
M=0.1

dx=1
dt=1

lx=int(sys.argv[1])
phi_o=float(sys.argv[2])
nstep=100000

def init_phi(lx, phi_o):
    """
    Returns randomly generated lx*ly array based on the specified initial condition, phi_o.
    
    :param lx: dimensionality of desired output lattice.
    :param phi_o: initial condition (0.0 or 0.5)
    :return: lx*ly array of randomly generated sites 
    """    
    ly=lx

    #initialise empty NxN array 
    phi=np.zeros((lx,ly),dtype=float) 

    #initialise values of phi randomly
    for i in range(lx):
        for j in range(ly):
            r=random.random()*0.2-0.1
            phi[i,j]=phi_o+r    

    return phi

def grad(lattice): 
    """
    Uses dell operator on NxN lattice to estimate the gradient field 
    
    :param lattice: NxN lattice
    :return: gradient vector field of NxN lattice 
    """    

    grad_x, grad_y = np.gradient(lattice)

    return np.array([grad_x, grad_y])

def grad_sq(lattice, dx):
    """
    Estimates laplacian of a NxN lattice 
    
    :param lattice: NxN lattice
    :param dx: spatial discretisation step 
    :return: lattice of laplacian values
    """  
    laplacian = (1/dx**2)*np.roll(lattice,1,axis=1) + np.roll(lattice,-1,axis=1) + np.roll(lattice,1,axis=0) + np.roll(lattice,-1,axis=0) - 4*lattice
    return laplacian

def chem_pot(phi, a, b, kappa, dx):
    """
    Calculates a lattice of chemical potentials, given an initial array of phi values
    
    :param phi: NxN array of phi values
    :params a, b, kappa: numerical constants
    :param dx: spatial timestep
    :return: lattice of chemical potentials
    """

    return -a*phi+b*np.power(phi, 3)-kappa*grad_sq(phi, dx)

def phi_new(phi, M, a, b, kappa, dx, dt):
    """
    Updates potential lattice using Euler algorithm.
    
    :param phi: NxN array of phi values
    :params a, b, kappa: numerical constants
    :param dx: spatial timestep
    :return: lattice of chemical potentials
    """
    mu=chem_pot(phi, a, b, kappa, dx)

    phi_updated=np.copy(phi)+M*dt*grad_sq(mu, dx)

    return phi_updated

def free_energy(phi, a, kappa, dx):
    """
    Calculates the free energy of a phi given phi lattice
    
    :param phi: NxN array of phi values
    :params a, kappa: numerical constants
    :param dx: spatial timestep
    :return: lattice containing free energy values
    """
    lx=len(phi[0])

    f=-0.5*a*np.square(phi)+0.25*a*np.power(phi, 4)+0.5*kappa*np.square(np.gradient(phi, axis=(0,1)))

    return np.sum(f)

def main():
    phi=init_phi(lx, phi_o) #initialise initial condition

    plt.cla() #show initial condition
    im=plt.imshow(phi, cmap=cmap, animated=True, vmin=-1, vmax=1)
    plt.draw()
    plt.pause(0.0001)

    F=[]

    for n in range(nstep):

        phi_updated=phi_new(np.copy(phi), M, a, b, kappa, dx, dt)  #update phi    
        phi=phi_updated

        #calculate free energy every 10 updates
        if n%10==0: 
            F.append(free_energy(phi, a, kappa, dx))

        #show animation every 250 sweeps
        if n%250==0:
            plt.cla()
            im=plt.imshow(phi, cmap=cmap, animated=True, vmin=-1, vmax=1)
            plt.colorbar()
            plt.draw()
            plt.pause(0.0001)
            plt.clf()

            print('Number of runs completed: ' + str(n))

    df=pd.DataFrame()
    df['Number of updates']=10*np.arange(nstep*0.1)
    df['Free energy']=np.array(F)

    print(df)
    df.to_excel('data.xlsx') #save datafile

main()