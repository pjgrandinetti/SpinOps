#include <stdio.h>

double fac(double x)
{
    if (x < 0) {
        fprintf(stderr, "illegal argument x = %g in factorial...\n",x);
        return 0;
    }
    int ix = (int) x;
    double sum = 1;
    for (; ix > 1; ix--) sum *= ix;
    return sum;
}

