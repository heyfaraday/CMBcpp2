#include "math.h"

#include "functionals.hpp"

#include "parameters.hpp"
#include "constants.hpp"

double area(long double** map, long double level) {

    double sum_without_norm = 0.0;
    double norm = 0.0;
    double theta;

    for (unsigned int j = 1; j < npix / 2; ++j) {

        theta = 2.0 * j * PI / npix;

        for (unsigned int i = 0; i < npix; ++i) {

            if (map[i][j] > level)
                sum_without_norm += sin(theta);

            if (map[i][j] != 0.0)
                norm = norm + sin(theta);
        }
    }

    return sum_without_norm / norm;
}

double length(long double** map, long double level) {

    return 0.0;
}
