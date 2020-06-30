#include <fstream>
#include <math.h>
#include <iostream>
#include <stdio.h> 
#include <stdlib.h>
#include <iomanip>
#include <utility>
#include <vector>
#include <string>
#include "particle.h"
#include "force.h"
#include "VectorMath.h"
#include "box.h"

using namespace std;

// initialize particles positions on a lattice box and initialize velocities
void particle::initialize(double BoxLength, int NParticles, double T, double dt)
{
	int particles = 0;

	// resize NParticles x 3 matrices to hold positions, previous positions, and velocities
	positionVector.resize(NParticles);
    positionPrevVector.resize(NParticles);
    velocityVector.resize(NParticles);

	// calculate length of lattice box, or alternatively the length a particle "takes up"
	double ParticleLength = BoxLength/pow(NParticles,1.0/3.0);

	// calculate number of lattice boxes on each side of big box, also the number of particles in 1 dimension of box
	double ParticlesPerSide = pow(NParticles, 1.0/3.0);	

	// initialize positions, looping over x, y, z directions of lattice to place particle at each point
	for (int x = 0; x < ParticlesPerSide; x++)
	{
		for (int y = 0; y < ParticlesPerSide; y++)
		{
			for (int z = 0; z < ParticlesPerSide; z++) 
			{
				// x,y,z coordinates for 1 single particle
				// position = cooridinate * lattice box length
				position.x = x * ParticleLength; 
				position.y = y * ParticleLength;
				position.z = z * ParticleLength;

				// put that position of 1 particle in positionVector corresponding to particle number
				positionVector[particles] = position;

				particles = particles + 1;

			}
		}
	}

	// initialize velocity variables to 0
    double totalVx = 0;
    double totalVy = 0;
    double totalVz = 0;
    double totalV2 = 0;

	// initialize velocities using random numbers
	for (int j = 0; j < NParticles; j++)
	{
		// draw random number from -0.5 and 0.5 for velocity of 1 particle
        velocityVector[j].x = (rand() / double(RAND_MAX)) - 0.5;
        velocityVector[j].y = (rand() / double(RAND_MAX)) - 0.5;
        velocityVector[j].z = (rand() / double(RAND_MAX)) - 0.5;
       
        // add random new velocity to total velocity 
        totalVx = totalVx + velocityVector[j].x;
        totalVy = totalVy + velocityVector[j].y;
        totalVz = totalVz + velocityVector[j].z;
       
		// total velocity^2 calculated from all 3 components 
        totalV2 = totalV2 + (velocityVector[j].x)*(velocityVector[j].x) + (velocityVector[j].y)*(velocityVector[j].y) + (velocityVector[j].z)*(velocityVector[j].z);

	}

	// average total velocities to each particle
	double averageVx;
	double averageVy;
	double averageVz;
	double averageV2;

	// average velocities over all particles
    averageVx = totalVx/NParticles;
    averageVy = totalVy/NParticles;
    averageVz = totalVz/NParticles;
    averageV2 = totalV2/NParticles;
    
    // scale factor for velocities
    double fScale = sqrt((3*T)/averageV2);
	box boxObj;
    
    for (int k = 0; k < NParticles; k++)
    {
		// subtracting by averageV gives center of mass = 0 so box does not shift
		// fscale assigns the correct velocity to produce the input temperature
        velocityVector[k].x = (velocityVector[k].x - averageVx) * fScale;
        velocityVector[k].y = (velocityVector[k].y - averageVy) * fScale;
        velocityVector[k].z = (velocityVector[k].z - averageVz) * fScale;
       
		// initialize matrix of 'previous' positions by moving 1 time step
        positionPrevVector[k].x = positionVector[k].x - (velocityVector[k].x)*dt;
        positionPrevVector[k].y = positionVector[k].y - (velocityVector[k].y)*dt;
        positionPrevVector[k].z = positionVector[k].z - (velocityVector[k].z)*dt;

		// periodic boundary conditions
		positionPrevVector[k].x = boxObj.pbc(positionPrevVector[k].x, BoxLength);
		positionPrevVector[k].y = boxObj.pbc(positionPrevVector[k].y, BoxLength);
		positionPrevVector[k].z = boxObj.pbc(positionPrevVector[k].z, BoxLength);
    }

}

// verlet algorithm to integrate forces and move particles
double particle::verlet(double dt, double NParticles, double BoxLength, vector<double3> forceV)
{
	// define variable: total velocities in each direction, total velocity^2, new positions, particle velocities
    double totalVx = 0;
    double totalVy = 0;
    double totalVz = 0;
    double totalV2 = 0;
    double newX = 0;
    double newY = 0;
    double newZ = 0;
    double velocityx = 0;
    double velocityy = 0;
    double velocityz = 0;
    
    // define objects to use functions in other classes 
    force forceObj;
	box boxObj;

	// loop through all particles to determine new positions
    for (int a = 0; a < NParticles; a++)
    {
		// calculate new position coordinates from verlet algorithm
        newX = 2 * positionVector[a].x - positionPrevVector[a].x + dt*dt * forceV[a].x;
        newY = 2 * positionVector[a].y - positionPrevVector[a].y + dt*dt * forceV[a].y;
        newZ = 2 * positionVector[a].z - positionPrevVector[a].z + dt*dt * forceV[a].z;

		// periodic boundary conditions
		newX = boxObj.pbc(newX,BoxLength);
		newY = boxObj.pbc(newY,BoxLength);
		newZ = boxObj.pbc(newZ,BoxLength);

        // calculate new velocity from verlet algorithm
        velocityx = (boxObj.boundaries(newX - positionPrevVector[a].x,BoxLength))/(2*dt);
        velocityy = (boxObj.boundaries(newY - positionPrevVector[a].y,BoxLength))/(2*dt);
        velocityz = (boxObj.boundaries(newZ - positionPrevVector[a].z,BoxLength))/(2*dt);

        // velocity center of mass
        //totalVx = totalVx + velocityx;
        //totalVy = totalVy + velocityy;
        //totalVz = totalVz + velocityz;
        
        // use velocity^2 to calculate total kinetic energy
        totalV2 = totalV2 + velocityx*velocityx + velocityy*velocityy + velocityz*velocityz;

		// define new positions as current positions, current positions as previous positions
        positionPrevVector[a].x = positionVector[a].x;
        positionPrevVector[a].y = positionVector[a].y;
        positionPrevVector[a].z = positionVector[a].z;
        positionVector[a].x = boxObj.pbc(newX,BoxLength);
        positionVector[a].y = boxObj.pbc(newY,BoxLength);
        positionVector[a].z = boxObj.pbc(newZ,BoxLength);
    }
    
	// return total velocity^2
    return totalV2;
}

double particle::velocityverlet(int NParticles, double BoxLength, double cutoff, double density, double dt, int cellist, vector<double3> forceV)
{
	// define variable: total velocities in each direction, total velocity^2, new positions, particle velocities
    double totalVx = 0;
    double totalVy = 0;
    double totalVz = 0;
    double totalV2 = 0;
    double velocityx = 0;
    double velocityy = 0;
    double velocityz = 0;
    
	// define objects to use functions in other classes 
    force forceObj;
	box boxObj;

	double energy;

    // calculate new position coordinates 
    for (int a = 0; a < NParticles; a++)
    {
        positionVector[a].x = positionVector[a].x + dt*velocityVector[a].x + dt*dt*forceV[a].x/2;
        positionVector[a].y = positionVector[a].y + dt*velocityVector[a].y + dt*dt*forceV[a].y/2;
        positionVector[a].z = positionVector[a].z + dt*velocityVector[a].z + dt*dt*forceV[a].z/2;

		// periodic boundary conditions
		positionVector[a].x = boxObj.pbc(positionVector[a].x, BoxLength);
		positionVector[a].y = boxObj.pbc(positionVector[a].y, BoxLength);
		positionVector[a].z = boxObj.pbc(positionVector[a].z, BoxLength);
	
		// first half velocity calculation using velocity verlet algorithm 
		velocityVector[a].x = velocityVector[a].x + forceV[a].x*dt/2;
        velocityVector[a].y = velocityVector[a].y + forceV[a].y*dt/2;
        velocityVector[a].z = velocityVector[a].z + forceV[a].z*dt/2;
    }

	// get new force from new positions, using either a cellist or no cellist
	if (cellist == 0)
	{
		energy = forceObj.getForce(NParticles,BoxLength,cutoff,density,positionVector);
	}
	else if (cellist == 1)
	{
		energy = forceObj.getForcecellist(dt, NParticles,BoxLength,cutoff,density,positionVector);
	}

	// obtain new forces in force vector from force class
	vector<double3> newForce = forceObj.getforceVector();

	// calculate new second half velocity calculation from new force
	for (int a = 0; a < NParticles; a++)
	{
		// calculate new velocity from verlet algorithm
		velocityVector[a].x = velocityVector[a].x + forceV[a].x*dt/2;
		velocityVector[a].y = velocityVector[a].y + forceV[a].y*dt/2;
		velocityVector[a].z = velocityVector[a].z + forceV[a].z*dt/2;
        
		// velocity center of mass
		//totalVx = totalVx + velocityVector[a].x;
		//totalVy = totalVy + velocityVector[a].y;
		//totalVz = totalVz + velocityVector[a].z;
        
		// use velocity^2 to calculate total kinetic energy
		totalV2 = totalV2 + velocityVector[a].x*velocityVector[a].x + velocityVector[a].y*velocityVector[a].y + velocityVector[a].z*velocityVector[a].z;
    }
    // return total velocity^2
    return totalV2;
}


vector<double3> particle::getpositionVector()
{
  return positionVector;
}



