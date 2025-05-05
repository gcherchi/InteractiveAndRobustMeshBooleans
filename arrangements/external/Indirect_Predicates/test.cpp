#include "implicit_point.h"

int main(int argc, char *argv[])
{
	explicitPoint2D a(1, 1), b(3, 3), c(2, 1), d(1, 2);
	implicitPoint2D_SSI i(a, b, c, d);

	if (genericPoint::orient2D(a, i, b) == 0) std::cout << "Collinear - Test succeeded.\n";
	else std::cout << "Not collinear - Test failed.\n";

	return 0;
}
