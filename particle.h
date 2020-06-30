#ifndef PARTICLE_H
#define PARTICLE_H

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
#include "force.h"
#include "box.h"

using namespace std;

class particle
{
	public:
		void initialize(double BoxLength, int NParticles, double T, double dt);
		double verlet(double dt, double NParticles, double BoxLength, vector<double3> forceV);
		double velocityverlet(int NParticles, double BoxLength, double cutoff, double density, double dt, int cellist, vector<double3> forceV);
        vector<double3> getpositionVector();

	protected:
			
	private:
		double3 position;
		vector<double3> positionVector;
        vector<double3> positionPrevVector;
		vector<double3> velocityVector;
};

#endif
