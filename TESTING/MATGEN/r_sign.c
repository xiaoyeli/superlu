double r_sign(float *a, float *b)
{
double x;
x = (*a >= 0 ? *a : - *a);
return( *b >= 0 ? x : -x);
}
