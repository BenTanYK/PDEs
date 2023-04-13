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
algorithm=str(sys.argv[3])

from poisson import error, jacobi, gauss_seidel

def init_J(lx):
    """
    Initialises current density array for a straight wire parallel to the x-axis.
    
    :param lx: dimensionality of desired output lattice.
    :return: lx*lx*lx array of zeros, with a single array of length lx passing through the
    centre of the yz plane. 
    """    
    J=np.zeros((lx, lx, lx))
    c=int(lx/2)
    J[:,c,c]=np.ones(lx)

    return J

def B_field(A):
    """
    Calculates the magnetic field based on a lattice describing the magnetic potential field.
    
    :param A: lx*lx*lx magnetic potential field.
    :return: 2*lx*lx*lx array containing the magnetic field y and z components (the x component
    is zero at all points, by symmetry). 
    """    
    dx,dy,dz= np.gradient(A)

    B_y=dz
    B_z=-1*dy

    return np.array([B_y, B_z])

def B_pot(A):
    """
    Returns dataframe containing lattice coordinates, and the corresponding potential and magnetic field
    at each point
    
    :param phi: lx*lx*lx potential lattice
    :return: pandas dataframe containing magnetic field and potential data
    """  
    mag_field=B_field(A) #generate electric field array

    lx=len(A[0])
    c=int(lx/2)

    X=[] #cartesian coordinates, origin at lattice centre
    Y=[]
    Z=[]
    R=[] #radial distance to centre
    pot=[] #potential
    Bx=[] #magnetic field
    By=[]
    Bz=[]
    B_abs=[]

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

                pot.append(A[i,j,k]) #potential

                B_x=0                       #magnetic field components
                B_y=-1*mag_field[0][i,j,k] 
                B_z=-1*mag_field[1][i,j,k]
                B=np.sqrt(B_x**2+B_y**2+B_z**2)

                if B!=0:
                    B_x=B_x/B
                    B_y=B_y/B
                    B_z=B_z/B

                Bx.append(B_x)
                By.append(B_y)
                Bz.append(B_z)
                B_abs.append(B)
                
    data=np.vstack((X,Y,Z,R,pot,Bx,By,Bz,B_abs))

    return data

def main():

    A=np.zeros((lx,lx,lx),dtype=float) 
    J=init_J(lx)

    err=10
    n=0

    while err>convergence:

        n+=1
        if n%10==0: 
            print(n)
            print(err)

        if algorithm=='jacobi' or algorithm=='Jacobi':
            A_new=jacobi(A,J)

        elif algorithm=='gauss_seidel' or algorithm=='Gauss_seidel':
            A_new=gauss_seidel(np.copy(A),J)

        else:
            print('select jacobi or gauss_seidel algorithm')

        err=error(A, A_new)
        A=A_new

    # plt.cla()
    # im=plt.imshow(A[int(lx/2), :, :], cmap='viridis', animated=True)
    # plt.colorbar()
    # plt.draw()
    # plt.pause(5)
    # #plt.savefig('mag_cut.png')
    # plt.clf()

    # data_cut=A[int(lx/2), :, :]
    # df_cut=pd.DataFrame(data_cut)
    # df_cut.to_excel('mag_pot.xlsx')

    data=B_pot(A)
    columns=np.array(['x', 'y', 'z', 'r', 'pot', 'Bx', 'By', 'Bz', '|B|'])
    df=pd.DataFrame(np.transpose(data), index=None, columns=columns)
    df.to_excel('mag_data.xlsx')
    print(df)

start_time = time.time() 
main()
print('Programme run time is:' + "%s seconds" % (time.time() - start_time))



