#include "f2c.h"

#define log10e 0.43429448190325182765

#include "math.h"

double d_lg10(double *x)
{
return( log10e * log(*x) );
}
