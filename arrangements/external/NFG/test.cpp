#include "numerics.h"
#include <chrono>

typedef interval_number Interval;

int main(int argc, char *argv[])
{
	// The following is deprecated. No longer necessary
	// initFPU();				
	
	// Set the rounding mode to +infinity to use intervals
	setFPUModeToRoundUP();	

	// Define three very different numbers as intervals
	Interval a(0.0);		
	Interval b(1.0e+12);
	Interval c(-1.9e+6);

	// Calculates an expression on intervals
	Interval d = (a+b) * c + c;

	// Print result
	std::cout << "Resulting interval: " << d << "\n";
	std::cout << "Interval width (must be 256): " << d.sup()-d.inf() << "\n";

	// We no longer do computation on intervals.
	// Return the rounding mode to 'nearest' (default)
	setFPUModeToRoundNEAR();

	return 0;
}
