#include <fstream>
#include <math.h>
#include <iostream>
#include <stdio.h> 
#include <stdlib.h>
#include <iomanip>
#include <utility>
#include <vector>
#include <string>
#include <cmath>
#include <ctime>
#include "particle.h"
#include "VectorMath.h"
#include "force.h"
#include "box.h"

using namespace std;

int main () 
{
	// define parameters
	int NParticles; // number of particles in box
	double T; // temperature in reduced units
	double density; // density = Number of Particles / Volume of Box
	double cutoff; // maximum distance between any two particles where interaction energy is explicitly calculated
	double cellist; // algorithm that reduces calculation time from N2 to N
	double BoxLength;
	double tmax; // total run time
	double dt; // time step

	// set configuration parameters in reduced units 
	cout << "-------------------------------------------------------------------" << endl;
	cout << "Enter Parameters for 3D Lennard Jones Molecular Dynamics Simulation" << endl;
	cout << "(in reduced units)" << endl;
	cout << "-------------------------------------------------------------------" << endl;
	cout << "Number of Particles (only perfect cube numbers e.g. 512, 10648, 27000):";
	cin >> NParticles;
	cout << "Temperature (e.g. 1):";
	cin >> T;
	cout << "Density (e.g. 0.5):";
	cin >> density;
	cout << "Cutoff distance (e.g. 3):";
	cin >> cutoff;
	cout << "Run Time (e.g. 5 seconds):";
	cin >> tmax;
	cout << "Time Step (e.g. 0.005):";
	cin >> dt;
	cout << "Use cellist? (0 for NO, 1 for YES)";
	cin >> cellist;
	BoxLength = pow((NParticles/density), 1.0/3.0);
	cout << "BoxLength: " << BoxLength << endl;
	cout << "-------------------------------------------------------------------" << endl;

	// pre-assigned parameters
    /*T = 1.6; 
    NParticles = 512; 
	density = 0.7;
    cutoff = 3;
	BoxLength = pow((NParticles/density), 1.0/3.0);
	cellist = 1; // 0 for no cellist, 1 for cellist
	tmax = 1; // total run time
    dt = 0.005; // time step*/

	cout << "BoxLength: " << BoxLength << endl;
    
	double energy;
	double ETotal;
	double rtotal;

	// print to output file
    ofstream trajFile;
    trajFile.open("traj.xyz"); // trajectory file (position coordinates of every particle at every timestep, can view movie in VMD)
    ofstream logFile;
    logFile.open("test.log"); // output of temperature and energy at each timestep

	// clock (to calculate total run time of simulation)
	clock_t start, end;

	// define objects to use functions in other classes 
	particle particleObj;
	box boxObj;
	force forceObj;

	// initialize positions of particles on a lattice in a box and initialize velocities
    particleObj.initialize(BoxLength, NParticles, T, dt); 
	vector<double3> initial = particleObj.getpositionVector(); // save initial positions of particles

	// start clock
	start = clock();

	// print header for screen output
	cout << "-------------------------------------------------------------------" << endl;
	cout << endl;
	cout << "Timestep" << "          Temperature" << "         Energy" << endl;
	cout << endl;
	cout << "-------------------------------------------------------------------" << endl;

    for (double t = 0; t<tmax; t=t+dt)
    {

		// call getForce (no cellist) or getForcecellist (with cellist) to calculate forces on all particles (in matrix called forceVector) and return total energy from particle interactions
		if (cellist == 0)
		{
			energy = forceObj.getForce(NParticles, BoxLength, cutoff, density, particleObj.getpositionVector());
		}

		if (cellist == 1)
		{
			energy = forceObj.getForcecellist(t, NParticles, BoxLength, cutoff, density, particleObj.getpositionVector());
		}
    
        // access forceVector from force class and place in variable called forceV
		vector<double3> forceV = forceObj.getforceVector();

		// integrate forces using verlet algorithm and return total velocity^2 squared (to be used in Temperature and Energy calculation later)
		double totalV2 = particleObj.verlet(dt, NParticles, BoxLength, forceV);

		// velocity verlet algorithm
		//double totalV2 = particleObj.velocityverlet(NParticles, BoxLength, cutoff, density, dt, cellist, forceV);

		// temperature calculation using averaged velocity: velocity^2/NParticles
		T = totalV2 / (3*NParticles);

		// total energy calculating using total velocity^2 and energy using averaged velocity
		ETotal = (energy + 0.5*totalV2)/NParticles;

		// print to screen the (1) timestep, (2) temperature, (3) energy 
		cout << "Timestep: " << t << " ---  T: " << T << " ---  E: " << ETotal << endl;

		// print to output file the (1) timestep, (2) temperature, (3) energy 
		logFile << t << "  " << T << "  " << ETotal << endl;

		// print particle positions at every time step to trajectory xyz file 
		vector<double3> positionOutput = particleObj.getpositionVector();
		trajFile << NParticles << endl;
		trajFile << endl;
		for (int i = 0; i < NParticles; i++)
		{
			trajFile << "1 " << positionOutput[i].x << " " << positionOutput[i].y << " " << positionOutput[i].z << endl;
		}
    }

	end = clock();
	// print total CPU time from simulation
	cout << "CPU time: " << ((float)(end - start)/CLOCKS_PER_SEC) << endl;

	logFile.close();
	trajFile.close();

    return 0;
}
