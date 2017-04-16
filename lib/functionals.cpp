#include "math.h"

#include "functionals.hpp"

#include "parameters.hpp"
#include "constants.hpp"

long double area(long double** map, long double level) {

    long double sum_without_norm = 0.0L;
    long double norm = 0.0L;
    long double theta;

    for (unsigned int j = 1; j < npix / 2; ++j) {

        long double long_j = static_cast<long double>(j);

        theta = 2.0L * long_j * PI / long_npix;

        for (unsigned int i = 0; i < npix; ++i) {

            if (map[i][j] > level)
                sum_without_norm += sinl(theta);

            if (map[i][j] != 0.0L)
                norm = norm + sinl(theta);
        }
    }

    return sum_without_norm / norm;
}

long double length(long double** map, long double level) {

    long double l = 0.0L;
    long double h = 2.0L * PI / long_npix;

    for (unsigned int i = 0; i < npix; ++i) {
        for (unsigned int j = 1; j < npix / 2 - 1; ++j) {

            map[i][j] = map[i][j] - level;

            if (map[i][j] == 0.0L and (map[i][j+1] == 0.0L or map[i+1][j] == 0.0L)) {
                l += 0.0L;
            }


//            if f[i][j] * f[i][j + 1] < 0.0:
//
//            if f[i][j] * f[i + 1][j] < 0.0:
//
//            phi1 = x[i][j]
//            theta1 = y[i][j] + h_theta * fabs(f[i][j]) / (fabs(f[i][j]) + fabs(f[i][j + 1]))
//
//            phi2 = x[i][j] + h_phi * fabs(f[i][j]) / (fabs(f[i][j]) + fabs(f[i + 1][j]))
//            theta2 = y[i][j]
//
//            l += s2(phi1, phi2, theta1, theta2)
//
//            elif f[i + 1][j] * f[i + 1][j + 1] < 0.0:
//
//            phi1 = x[i][j]
//            theta1 = y[i][j] + h_theta * fabs(f[i][j]) / (fabs(f[i][j]) + fabs(f[i][j + 1]))
//
//            phi2 = x[i + 1][j]
//            theta2 = y[i + 1][j] + h_theta * fabs(f[i + 1][j]) / (fabs(f[i + 1][j]) + fabs(f[i + 1][j + 1]))
//
//            l += s2(phi1, phi2, theta1, theta2)
//
//            elif f[i][j + 1] * f[i + 1][j + 1] < 0.0:
//
//            phi1 = x[i][j]
//            theta1 = y[i][j] + h_theta * fabs(f[i][j]) / (fabs(f[i][j]) + fabs(f[i][j + 1]))
//
//            phi2 = x[i][j + 1] + h_phi * fabs(f[i][j + 1]) / (fabs(f[i][j + 1]) + fabs(f[i + 1][j + 1]))
//            theta2 = y[i][j + 1]
//
//            l += s2(phi1, phi2, theta1, theta2)
//
//            elif f[i][j] * f[i + 1][j] <= 0.0:
//
//            if f[i + 1][j] * f[i + 1][j + 1] < 0.0:
//
//            phi1 = x[i][j] + h_phi * fabs(f[i][j]) / (fabs(f[i][j]) + fabs(f[i + 1][j]))
//            theta1 = y[i][j]
//
//            phi2 = x[i + 1][j]
//            theta2 = y[i + 1][j] + h_theta * fabs(f[i + 1][j]) / (fabs(f[i + 1][j]) + fabs(f[i + 1][j + 1]))
//
//            l += s2(phi1, phi2, theta1, theta2)
//
//            elif f[i][j + 1] * f[i + 1][j + 1] < 0.0:
//
//            phi1 = x[i][j] + h_phi * fabs(f[i][j]) / (fabs(f[i][j]) + fabs(f[i + 1][j]))
//            theta1 = y[i][j]
//
//            phi2 = x[i][j + 1] + h_phi * fabs(f[i][j + 1]) / (fabs(f[i][j + 1]) + fabs(f[i + 1][j + 1]))
//            theta2 = y[i][j + 1]
//
//            l += s2(phi1, phi2, theta1, theta2)
//
//            elif f[i + 1][j] * f[i + 1][j + 1] < 0.0:
//
//            if f[i][j + 1] * f[i + 1][j + 1] < 0.0:
//            phi1 = x[i + 1][j]
//            theta1 = y[i + 1][j] + h_theta * fabs(f[i + 1][j]) / (fabs(f[i + 1][j]) + fabs(f[i + 1][j + 1]))
//
//            phi2 = x[i][j + 1] + h_phi * fabs(f[i][j + 1]) / (fabs(f[i][j + 1]) + fabs(f[i + 1][j + 1]))
//            theta2 = y[i][j + 1]
//
//            l += s2(phi1, phi2, theta1, theta2)

            map[i][j] = map[i][j] + level;
        }
    }

    return l / (4.0L * PI);
}
