#include "math.h"

#include "functionals.hpp"

#include "parameters.hpp"
#include "constants.hpp"
#include "distance.hpp"
#include "utils.hpp"

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

    long double phi;
    long double theta;

    long double phi1;
    long double phi2;
    long double theta1;
    long double theta2;

    for (unsigned int i = 0; i < npix; ++i) {
        for (unsigned int j = 1; j < npix / 2 - 1; ++j) {

            phi = 2.0L * PI * i / long_npix;
            theta = 2.0L * PI * j / long_npix;

            if ((map[i][j] - level) == 0.0L and ((map[i][j + 1] - level) == 0.0L or (map[i + 1][j] - level) == 0.0L)) {
                l += 0.0L;
            }

            if ((map[i][j] - level) * (map[i][j + 1] - level) < 0.0L) {

                if ((map[i][j] - level) * (map[i + 1][j] - level) < 0.0L) {

                    phi1 = phi;
                    theta1 = theta + h * fabsl((map[i][j] - level)) / (fabsl((map[i][j] - level)) + fabsl((map[i][j + 1] - level)));

                    phi2 = phi + h * fabsl((map[i][j] - level)) / (fabsl((map[i][j] - level)) + fabsl((map[i + 1][j] - level)));
                    theta2 = theta;

                    l += s2(phi1, phi2, theta1, theta2);

                } else if ((map[i + 1][j] - level) * (map[i + 1][j + 1] - level) < 0.0L) {

                    phi1 = phi;
                    theta1 = theta + h * fabsl((map[i][j] - level)) / (fabsl((map[i][j] - level)) + fabsl((map[i][j + 1] - level)));

                    phi2 = phi + h;
                    theta2 = theta + h * fabsl((map[i + 1][j] - level)) / (fabsl((map[i + 1][j] - level)) + fabsl((map[i + 1][j + 1] - level)));

                    l += s2(phi1, phi2, theta1, theta2);

                } else if ((map[i][j + 1] - level) * (map[i + 1][j + 1] - level) < 0.0L) {

                    phi1 = phi;
                    theta1 = theta + h * fabsl((map[i][j] - level)) / (fabsl((map[i][j] - level)) + fabsl((map[i][j + 1] - level)));

                    phi2 = phi + h * fabsl((map[i][j + 1] - level)) / (fabsl((map[i][j + 1] - level)) + fabsl((map[i + 1][j + 1] - level)));
                    theta2 = theta + h;

                    l += s2(phi1, phi2, theta1, theta2);
                }
            } else if (map[i][j] * map[i + 1][j] < 0.0L) {

                if (map[i + 1][j] * map[i + 1][j + 1] < 0.0L) {

                    phi1 = phi + h * fabsl((map[i][j] - level)) / (fabsl((map[i][j] - level)) + fabsl((map[i + 1][j] - level)));
                    theta1 = theta;

                    phi2 = phi + h;
                    theta2 = theta + h * fabsl((map[i + 1][j] - level)) / (fabsl((map[i + 1][j] - level)) + fabsl((map[i + 1][j + 1] - level)));

                    l += s2(phi1, phi2, theta1, theta2);

                } else if (map[i][j + 1] * map[i + 1][j + 1] < 0.0L) {

                    phi1 = phi + h * fabsl((map[i][j] - level)) / (fabsl((map[i][j] - level)) + fabsl((map[i + 1][j] - level)));
                    theta1 = theta;

                    phi2 = phi + h * fabsl((map[i][j + 1] - level)) / (fabsl((map[i][j + 1] - level)) + fabsl((map[i + 1][j + 1] - level)));
                    theta2 = theta + h;

                    l += s2(phi1, phi2, theta1, theta2);
                }
            } else if (map[i + 1][j] * map[i + 1][j + 1] < 0.0L) {

                if (map[i][j + 1] * map[i + 1][j + 1] < 0.0L) {

                    phi1 = phi + h;
                    theta1 = theta + h * fabsl((map[i + 1][j] - level)) / (fabsl((map[i + 1][j] - level)) + fabsl((map[i + 1][j + 1] - level)));

                    phi2 = phi + h * fabsl((map[i][j + 1] - level)) / (fabsl((map[i][j + 1] - level)) + fabsl((map[i + 1][j + 1] - level)));
                    theta2 = theta + h;

                    l += s2(phi1, phi2, theta1, theta2);
                }
            }
        }
    }

    return l / (PI);
}

int condition_1(long double xx, long double yy, long double xy) {
    if ((xx * yy - xy * xy >= 0.0L and xx >= 0.0L) or (xx * yy - xy * xy >= 0.0L and yy >= 0.0L)) {
        return 0;
    } else if ((xx * yy - xy * xy >= 0.0L and xx < 0.0L) or (xx * yy - xy * xy >= 0.0L and yy < 0.0L)) {
        return 2;
    } else {
        return 0;
    }
}

int condition_2(long double qx, long double qy, long double ux, long double uy) {

    long double root1;
    long double root2;

    d_solver(&root1, &root2, qx, qy, ux, uy);
    long double mean_root1 = mean_t_solver(root1, qx, qy, ux, uy);
    long double mean_root2 = mean_t_solver(root2, qx, qy, ux, uy);

    if (det(qx, qy, ux, uy) < 0) {
        return 1;
    } else if ((mean_root1 > 0.0L and mean_root2 < 0.0L) or (mean_root1 < 0.0L and  mean_root2 > 0.0L)) {
        return 2;
    } else {
        return 3;
    }
}

void points_classifier(long double** map, long double** map_x, long double** map_y,
                       long double** map_xx, long double** map_yy, long double** map_xy) {

}

void singular_points_classifier(long double** map, long double** map_x, long double** map_y,
                       long double** map_xx, long double** map_yy, long double** map_xy) {

}
