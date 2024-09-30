/*! \file
 * \ingroup TestingMatgen
 * \brief Power with integer exponent
 */

/*!
 * Power with integer exponent
 *
 * Implemented as Exponentiation by Squaring for speed.
 *
 * \param[in] base base as double value
 * \param[in] exp  exponent as integer value
 * \return Computes the value of \p base raised to the power \p exp
 */
double powi(double base, int exp)
{
    double x = base;
    int n = exp;
    double pow = 1.0;

    if (n != 0) {
        if(n < 0) {
            n = -n;
            x = 1 / x;
        }
        for( ; ; ) {
            if(n & 01)
                pow *= x;
            if(n >>= 1)
                x *= x;
            else
                break;
        }
    }
    return pow;
}

/*!
 * Power with integer exponent
 *
 * Implemented as Exponentiation by Squaring for speed.
 *
 * \param[in] base base as float value
 * \param[in] exp  exponent as integer value
 * \return Computes the value of \p base raised to the power \p exp
 */
float powif(float base, int exp)
{
    float x = base;
    int n = exp;
    float pow = 1.0;

    if (n != 0)  {
        if (n < 0) {
            n = -n;
            x = 1 / x;
        }
        for ( ; ; ) {
            if (n & 01)
                pow *= x;
            if (n >>= 1)
                x *= x;
            else
                break;
        }
    }
    return pow;
}
