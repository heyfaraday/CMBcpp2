#include "math.h"

#include "functionals.hpp"

double area(double** map, int npix, double level)
{
    double sum_without_norm = 0.0;
    double norm = 0.0;
    double theta = 0.0;

    for (int i = 0; i < npix; ++i) {
        for (int j = 0; j < npix/2; ++j) {

            if (map[i][j] > level)
            {
                sum_without_norm += sin(theta);
            }

            norm += sin(theta);

        }
    }

    return sum_without_norm / norm;
}
