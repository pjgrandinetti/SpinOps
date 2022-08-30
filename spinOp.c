#include <stdio.h>

double fac(double x)
{
    double sum = 1;
    int ix;
    
    if (x < 0) {
        fprintf(stderr, "illegal argument x = %g in factorial...\n",x);
        return 0;
    }
    ix = (int) x;
    for (; ix > 1; ix--) sum *= ix;
    return sum;
}

