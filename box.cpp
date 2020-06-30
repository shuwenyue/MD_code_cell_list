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

// periodic boundary conditions when calculating energy and forces
double box::boundaries (double length, double BoxLength)
{

	while (length > (BoxLength/2))
	{
		length = length - BoxLength;
	}

	while (length < (-BoxLength/2))
	{
		length = length + BoxLength;
	}
			
	return length;
}

// periodic boundary conditions to physically put particle back in box
double box::pbc (double position, double BoxLength)
{
	if (position > BoxLength)
    {
        position = position - BoxLength;                                    
    }

	if (position < 0)
    {
        position = position + BoxLength;                                    
    }
	
	return position;
}
