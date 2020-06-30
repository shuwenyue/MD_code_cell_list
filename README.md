
## 3D Lennard Jones Particles Molecular Dynamics Simulation.

**APC 523 Term Project**

*Shuwen Yue*
*Sassan Hajirezaie*
*Andrew Rock*

## RUNNING

To run the program, follow the below instruction:
(Please follow recommended values for inputs)

1. make
2. ./output
3. Enter number of particles (only perfect cube numbers are allowed) (e.g. 512)
4. Enter temperature (e.g. 1)
5. Enter density (e.g. 0.5)
6. Enter cutoff distance (e.g. 3)
7. Enter run time (e.g. 5)
8. Enter time step (e.g. 0.005)
9. Enter 0 for no cell list, 1 for cell list

Screen output:
-------------
The program will run and output a series of "timestep ------ temperature ------- energy"
When the program ends, CPU time for the run will be printed.

Output files
------------
A logfile containing "timestep ------ temperature ------- energy" will be printed to "test.log", this data can be plotted - see graphic subdirectory.
A xyz file of trajectories at every timestep will be printed to "traj.xyz", this data can be used to view a movie of the simulation - see graphics subdirectory.

One can run the simulation using the same parameters both using the cell list and without the cell list to observe the CPU time advantage of adding a cell list

----------------------------------------------------------------------
DIRECTORY CONTENTS
----------------------------------------------------------------------

The directory contains the source code for the simulation.

Files
-----

main.cpp
	- This is the main program with input parameters and outline of how the simulation works

force.cpp
	- This class computes forces of particle interactions from particle positions
	- There are two methods to approach the particle interactions: no cell list and with cell list, both can be done in this class.
	
force.h

particle.cpp
	- This class initiate the particles on a lattice and holds the verlet algorithm integrator to move particles
	- Particle positions and particle velocities are held in this class.

particle.h

box.cpp
	- This class provides periodic boundary conditions for the particles in energy/force calculations and for placing particles back inside the box

box.h

VectorMath.h
	- This is a tool to allow for creation of vectors that hold 3 double precision variables

Makefile

traj.xyz
	- example trajectory xyz file, can be used in VMD

output
	- executable


Subdirectories:
---------------
analysis: this holds python scripts used for analysis of data produced by simulation

isotherms: this holds isotherms of energy outputs 

graphics: this holds graphs for isotherms, CPUtime, snapshots from simulation movie, and movie file in .mpg



