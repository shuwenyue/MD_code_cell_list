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

// function to calculate total energy and forces on each particle using a radial cutoff distance
double force::getForce(int NParticles, double BoxLength, double cutoff, double density, vector<double3> coord)
{
    
	// define force vector initialized to 0
    forceVector.resize(NParticles);
    for (int f=0; f < NParticles; f++)
    {
        forceVector[f].x = 0;
        forceVector[f].y = 0;
        forceVector[f].z = 0;
    }

	// define initial energy to 0
	double energy = 0;

	// define objects to use functions in other classes 
	particle particleObj;
	box boxObj;

	// calculations for 'long range interaction' (mean field energy outside of cutoff), to be used in total Energy calculation
    double rcutoff6i = 1/(cutoff*cutoff*cutoff*cutoff*cutoff*cutoff); // 1/(cutoff^6) 
    double rcutoff12i = rcutoff6i*rcutoff6i; // 1/(cutoff^12)

	// define length variables
	double xlength, ylength, zlength, r2;

	// loop through every interaction of particle i and j, without double counting for each pair
	for (int i = 0; i < (NParticles-1); i++)
    {
        for (int j = i+1; j < NParticles; j++)
        {
            xlength = coord[i].x - coord[j].x; // interaction length between particle i and j in the x direction
            ylength = coord[i].y - coord[j].y; // interaction length between particle i and j in the y direction
            zlength = coord[i].z - coord[j].z; // interaction length between particle i and j in the z direction

			// apply periodic boundary conditions on interaction lengths
            xlength = boxObj.boundaries(xlength, BoxLength);
            ylength = boxObj.boundaries(ylength, BoxLength);
            zlength = boxObj.boundaries(zlength, BoxLength);

            // calculate radius^2
            r2 = (xlength*xlength) + (ylength*ylength) + (zlength*zlength);

			// define formulas for LJ force calculation 
			// LJ force formula is force = (48/cutoff^2)*(1/cutoff^12 - 1/cutoff^6) - we can break this down to multiple components
            double r6i = 1/(r2*r2*r2);
            double r12i = r6i*r6i;
            double lj = (48/r2)*(r12i-0.5*r6i);	// lj = force
           
			// if r < cutoff radius, then add explicitly calculated interaction force to cumulative force
			// otherwise, this interaction is outside of cutoff radius and not included in calculation
            if (r2 < (cutoff*cutoff))
            {
				// update force
                forceVector[i].x = forceVector[i].x + lj*xlength;
                forceVector[i].y = forceVector[i].y + lj*ylength;
                forceVector[i].z = forceVector[i].z + lj*zlength;
                
                forceVector[j].x = forceVector[j].x - lj*xlength;
                forceVector[j].y = forceVector[j].y - lj*ylength;
                forceVector[j].z = forceVector[j].z - lj*zlength;
               
				// interactions energy outside of cutoff distance represented by in 'long range correction' (8/3*pi*density...)
				// energy = energy + LJ energy of current interaction long range correction
				energy = energy + 4*(r12i-r6i) - (8/3)*3.14159265*density*((1/3)*(1/pow(cutoff,9.0))-(1/pow(cutoff,6.0)));
				
            }
        }
    }

	// returns total energy 
	return energy;
}

// function to calculate total energy and forces on each particle using cell list
double force::getForcecellist(double t, int NParticles, double BoxLength, double cutoff, double density, vector<double3> coord)
{
    // define force vector initialized to 0
    forceVector.resize(NParticles);
    for (int f=0; f < NParticles; f++)
    {
        forceVector[f].x = 0;
        forceVector[f].y = 0;
        forceVector[f].z = 0;
    }

	// define initial energy to 0
	double energy = 0;

	// define objects to use functions in other classes   
	particle particleObj;
	box boxObj;

	// calculations for 'long range interaction' (mean field energy outside of cutoff), to be used in total Energy calculation
    double rcutoff6i = 1/(cutoff*cutoff*cutoff*cutoff*cutoff*cutoff); // 1/(cutoff^6) 
    double rcutoff12i = rcutoff6i*rcutoff6i; // 1/(cutoff^12)

	// define length variables
	double xlength, ylength, zlength, r2;

	// build Cell List, determine number of cells long on each side of box (each cell length > cutoff)
	int CellsPerSide = round(BoxLength/cutoff);

	// length of each cell
	double cellLength = BoxLength/CellsPerSide;

	// total number of cells
	int Ncells = CellsPerSide*CellsPerSide*CellsPerSide;
	
	// ParticleCell vector contains the cell number location for each particle
	ParticleCell.resize(NParticles);	
	for (int a=0; a<NParticles; a++)
	{
		// convert position coordinates to cell coordinates
		double cx = floor(coord[a].x/cellLength); 
		double cy = floor(coord[a].y/cellLength);
		double cz = floor(coord[a].z/cellLength);
		ParticleCell[a] = CellsPerSide*CellsPerSide*cz+cy*CellsPerSide+cx; // algorithm to determine cell number location of a particle
	}
	
	// create head vector and linked list and initialize to -1 (NULL)
	HeadList.resize(Ncells);
	LinkedList.resize(NParticles);
	for (int m=0; m<Ncells; m++)
	{
		HeadList[m] = -1;
	}
	for (int m=0; m<NParticles; m++)
	{
		LinkedList[m] = -1;
	}
		
	// create cell coordinate matrix: NParticlesx3 matrix of cell coordinate for each particle
	CellCoord.resize(NParticles);
	int cellCounter = 0;
	for (int z = 0; z < CellsPerSide ; z++)
	{
		for (int y = 0; y < CellsPerSide ; y++)
		{
			for(int x = 0; x < CellsPerSide ; x++)
			{
				CellCoord[cellCounter].x = x;
				CellCoord[cellCounter].y = y;
				CellCoord[cellCounter].z = z;
				cellCounter = cellCounter + 1;
			}
		}
	}

	//loop over particle list and assign values to headlist and linked list
	for (int i = 0; i<NParticles; i++)
	{
		double cellnumber = ParticleCell[i];
		LinkedList[i]=HeadList[cellnumber];
		HeadList[cellnumber] = i;
	}


	// loop over all cells for interactions, where any 1 cell is designated 'center cell' for surrounding interactions
	for (int cell = 0; cell < Ncells; cell++)
	{
		// loop over all 26 neighbors of 'center cell' and 'center cell' itself (27 boxes total)
		for (int cx1 = CellCoord[cell].x-1; cx1 < CellCoord[cell].x+2; cx1++)
		{
			for (int cy1 = CellCoord[cell].y-1; cy1 < CellCoord[cell].y+2; cy1++)
			{
				for (int cz1 = CellCoord[cell].z-1; cz1 < CellCoord[cell].z+2; cz1++)
				{

					// periodic boundary conditions (correct cell locations of neighbors)
					int cx11 = cx1;
					int cy11 = cy1;
					int cz11 = cz1;

					if (cx1 < 0)
					{
						cx11 = CellsPerSide - 1;
					}
					else if (cx1 > (CellsPerSide-1))
					{
						cx11 = 0;
					}
					if (cy1 < 0)
					{
						cy11 = CellsPerSide - 1;
					}
					else if (cy1 > (CellsPerSide-1))
					{
						cy11 = 0;
					}
					if (cz1 < 0)
					{
						cz11 = CellsPerSide - 1;
					}
					else if(cz1 > (CellsPerSide-1))
					{
						cz11 = 0;
					}
					
					// use algorithm to determine cell number of neighbor cells
					int neighborCell = CellsPerSide*CellsPerSide*cz11+cy11*CellsPerSide+cx11;

					// define atom i as HeadList value of 'center cell'
					int i = HeadList[cell];
					
					// scan for atom i in center cell (cell)
					int counter = 0;
					while (i > -1) // while HeadList vector space is null
					{
						counter = counter + 1;

						// scan atom j in cell neighbor (neighborCell)
						int j = HeadList[neighborCell];
						
						while (j > -1) // while LinkedList vector space is null
						{
						
							if (i < j) // avoid double counting of each pair
							{
            					xlength = coord[i].x - coord[j].x; // interaction length between particle i and j in the x direction
            					ylength = coord[i].y - coord[j].y; // interaction length between particle i and j in the y direction
           						zlength = coord[i].z - coord[j].z; // interaction length between particle i and j in the z direction

								// apply periodic boundary conditions on interaction lengths
								xlength = boxObj.boundaries(xlength, BoxLength);
								ylength = boxObj.boundaries(ylength, BoxLength);
								zlength = boxObj.boundaries(zlength, BoxLength);

								// calculate radius^2
								r2 = (xlength*xlength) + (ylength*ylength) + (zlength*zlength);

								// define formulas for LJ force calculation 
								// LJ force formula is force = (48/cutoff^2)*(1/cutoff^12 - 1/cutoff^6) - we can break this down to multiple components
								double r6i = 1/(r2*r2*r2);
								double r12i = r6i*r6i;
								double lj = (48/r2)*(r12i-0.5*r6i);	// used for derivation of energy = force

								// if r < cutoff radius, then add explicitly calculated interaction force to cumulative force
								// otherwise, this interaction is outside of cutoff radius and not included in calculation
								if (r2 < (cutoff*cutoff))
								{
									// update force
									forceVector[i].x = forceVector[i].x + lj*xlength;
									forceVector[i].y = forceVector[i].y + lj*ylength;
                					forceVector[i].z = forceVector[i].z + lj*zlength;
                
               						forceVector[j].x = forceVector[j].x - lj*xlength;
               						forceVector[j].y = forceVector[j].y - lj*ylength;
                					forceVector[j].z = forceVector[j].z - lj*zlength;
							

									// interactions energy outside of cutoff distance represented by in 'long range correction' (8/3*pi*density...)
									// energy = energy + LJ energy of current interaction long range correction
									energy = energy + 4*(r12i-r6i) - (8/3)*3.14159265*density*((1/3)*(1/pow(cutoff,9.0))-(1/pow(cutoff,6.0)));
								}
							}
							j = LinkedList[j];
						}
						i = LinkedList[i];
					}
				}
			}
		}
	}
	// returns total energy
	return energy;
}



vector<double3> force::getforceVector()
{
	return forceVector;
}
