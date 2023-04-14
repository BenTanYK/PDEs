# PDEs
This repository contains files with code and data for the simulation of the Cahn-Hilliard and Poisson equations.  
**Please see the notebook 'CP3_analysis.ipynb' for a summary of all figures and datafiles**

### Command line arguments
- cahn_hill.py lx phi_o 
- monopole.py lx convergence algorithm
- mag.py lx convergence algorithm
- SOR.py lx convergence

For all simulations, lx was taken as 50 or 100. A value of convergence=0.01 was used. 

### The Cahn-Hilliard Equation
All code for simulation of the Cahn-Hilliard equation using the explicit Euler algorthim is housed in the file cahn_hill.py. The 2D lattice size, lx, is taken as the first command-line argument. The initial state, phi_o, is specified as the second command-line argument. phi_o=0.0 leads to spinodal decomposition, while phi_o=0.5 results in Ostwald ripening. 10^6 updates were performed on the latttice.

-Spinodal decomposition: python3 cahn_hill.py 50 0.0

-Ostwald ripening: python3 cahn_hill.py 50 0.5

### Poisson's equation
All functions to be used in simulations of Poisson's equations are housed in the file poisson.py. This includes functions for the 
Jacobi and Gauss-Seidel update rules. 

The programme monopole.py initialises a 3D point charge distribution and calculates the potential and electric fields, based on 
either the Jacobi or Gauss-Seidel algorithms. The lattice is updated until the error is less than the specified convergence limit.

-Jacobi: python3 monopole.py 50 0.01 jacobi
-Gauss-Seidel: python3 monopole.py 50 0.01 gauss_seidel

The successive over-relaxation method is carried out in the programme SOR.py. The SOR parameter omega is varied between 1.0 and 2.0,
over 0.01 intervals. The programme returns a dataframe containing the tested values of omega, and the corresponding number of updates
required for convergence.

-Command line arguments: python3 SOR.py 50 0.01

Finally, the programme mag.py contains the simulation of the potential and magnetic field for a straight wire. The programme initialises
a straight wire parallel to the x-axis, and uses either the Jacobi or Gauss-Seidel update rule. A dataframe containing lattice coordinates, and the corresponding potential and magnetic fields is returned. 

Jacobi: python3 mag.py 50 0.01 jacobi
Gauss-Seidel: python3 mag.py 50 0.01 gauss_seidel

### Data analysis
- The free-energy data for the Cahn-Hilliard equations is kept in the excel files 'spinodal_50.xlsx' and 'ostwald_50.xlsx'.
- Data for the magnetic field and electric field simulations is held in the files 'magnetic.xlsx' and 'electric.xlsx'.
- Data for the optimisation parameter, omega, is found in the file 'SOR_data.xlsx'.
- Cuts of the magnetic and electric potentials fields are shown in the images 'mag_cut.png' and 'elec_cut.png'.
- The notebook 'CP3_analysis.ipynb' summarises all findings - PLEASE READ.
