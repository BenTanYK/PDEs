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

def point_charge(lx):
    """
    Initialises charge density array for a point charge situation at the centre of the lx*lx*lx cubic lattice.
    
    :param lx: dimensionality of desired output lattice.
    :return: lx*lx*lx array of zeros, with a single point of value 1.0 at the centre. 
    """    
    #initialise empty NxNxN array 
    lattice=np.zeros((lx,lx,lx),dtype=float) 

    c=int(lx/2)

    lattice[c,c,c]=1.0

    return lattice

def error(phi, phi_new):
    """
    Calculates the lattice error for a single update of the potential lattice, phi
    
    :param phi: lx*lx*lx potential lattice
    :param phi_new: updated potential lattice
    :return: float giving the error in the lattice update
    """    
    diff=np.abs(phi_new-phi)
    return np.sum(diff)

def E_field(phi):
    """
    Calculates a lattice of electric field vectors calculated from a lattice of electric potentials
    
    :param phi: lx*lx*lx potential lattice
    :return: 3xlx*lx*lx lattice containing electric field vectors at each lattice point
    """  

    Ex,Ey,Ez= np.gradient(phi)

    return np.array([-1*Ex, -1*Ey, -1*Ez])

def E_pot(phi):
    """
    Returns dataframe containing lattice coordinates, and the corresponding potential and electric field
    at each point
    
    :param phi: lx*lx*lx potential lattice
    :return: pandas dataframe containing electric field and potential data
    """  

    elec_field=E_field(phi) #generate electric field array

    lx=len(phi[0])
    c=int(lx/2)

    X=[] #cartesian coordinates, origin at lattice centre
    Y=[]
    Z=[]
    R=[] #radial distance to centre
    pot=[] #potential
    Ex=[] #electric field
    Ey=[]
    Ez=[]
    E_abs=[]

    for i in range(lx):
        for j in range(lx):
            for k in range(lx):
                x=c-i #set centre of lattice as origin
                y=c-j
                z=c-k

                r=math.sqrt(x**2+y**2+z**2)

                X.append(x) #coordinates
                Y.append(y)
                Z.append(z)
                R.append(r) #radial distance to origin

                pot.append(phi[i,j,k]) #potential

                E_x=-1*elec_field[0][i,j,k] #electric field components
                E_y=-1*elec_field[1][i,j,k]
                E_z=-1*elec_field[2][i,j,k]
                E=np.sqrt(E_x**2+E_y**2+E_z**2)

                if E!=0:
                    E_x=E_x/E
                    E_y=E_y/E
                    E_z=E_z/E

                Ex.append(E_x)
                Ey.append(E_y)
                Ez.append(E_z)
                E_abs.append(E)
                
    data=np.vstack((X,Y,Z,R,pot,Ex,Ey,Ez,E_abs))

    return data

def jacobi(phi, rho):
    """
    Calculates an updated potential lattice based on the Jacobi rule
    
    :param phi: lx*lx*lx potential lattice
    :param rho: lattice containing charge distribution
    :return: updated phi lattice
    """  

    phi_new=(1/6)*(np.roll(phi,1,axis=0)+np.roll(phi,-1,axis=0)+np.roll(phi,1,axis=1)+np.roll(phi,-1,axis=1)+np.roll(phi,1,axis=2)+np.roll(phi,-1,axis=2)+rho)

    phi_new[0,:,:] = 0.0
    phi_new[:,0,:] = 0.0
    phi_new[:,:,0] = 0.0
    phi_new[-1,:,:] = 0.0
    phi_new[:,-1,:] = 0.0
    phi_new[:,:,-1] = 0.0

    return phi_new

def gauss_seidel(phi, rho):
    """
    Calculates an updated potential lattice based on the Gauss-Seidel rule
    
    :param phi: lx*lx*lx potential lattice
    :param rho: lattice containing charge distribution
    :return: updated phi lattice
    """  
    lx=len(phi[0])

    for i in range(1,lx-1):
        for j in range(1,lx-1):
            for k in range(1,lx-1):
                phi[i,j,k] = 1/6 * (phi[i+1,j,k] + phi[i-1,j,k]+ phi[i,j+1,k] + phi[i,j-1,k] + phi[i,j,k+1] + phi[i,j,k-1] + rho[i,j,k])

    return phi

def sor(phi, rho, omega):
    """
    Calculates an updated potential lattice based on successive over-relaxation

    :param phi: lx*lx*lx potential lattice
    :param rho: lattice containing charge distribution
    :param omega: SOR parameter
    :return: updated phi lattice
    """  
    lx=len(phi[0])

    for i in range(1,lx-1):
        for j in range(1,lx-1):
            for k in range(1,lx-1):
                phi[i,j,k] = (1-omega)*phi[i,j,k]+omega/6 * (phi[i+1,j,k] + phi[i-1,j,k]+ phi[i,j+1,k] + phi[i,j-1,k] + phi[i,j,k+1] + phi[i,j,k-1] + rho[i,j,k])

    return phi
 
