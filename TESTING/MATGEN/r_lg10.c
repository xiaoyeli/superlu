#define log10e 0.43429448190325182765

#include "math.h"

double r_lg10(float *x)
{
return( log10e * log(*x) );
}
