#ifndef FORCE_H
#define FORCE_H

#include <fstream>
#include <math.h>
#include <iostream>
#include <stdio.h> 
#include <stdlib.h>
#include <iomanip>
#include <utility>
#include <vector>
#include <string>
#include "VectorMath.h"
#include "particle.h"
#include "box.h"

using namespace std;

class force
{
	public:
        double getForce(int NParticles, double BoxLength, double cutoff, double density, vector<double3> coord);
	double getForcecellist(double t, int NParticles, double BoxLength, double cutoff, double density, vector<double3> coord);
	vector<double3> getforceVector();
	protected:
			
	private:
	vector<double3> forceVector;
        double3 forceCoord;
	vector<int> ParticleCell;
	vector<int> HeadList;
	vector<int> LinkedList;
	vector<double3> CellCoord;
};

#endif
