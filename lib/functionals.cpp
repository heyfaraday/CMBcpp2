#include "math.h"
#include <fstream>
#include <iostream>

#include "functionals.hpp"

#include "parameters.hpp"
#include "constants.hpp"
#include "distance.hpp"
#include "utils.hpp"
#include "fft.hpp"

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
            } else if ((map[i][j] - level) * (map[i + 1][j] - level) < 0.0L) {

                if ((map[i + 1][j] - level) * (map[i + 1][j + 1] - level) < 0.0L) {

                    phi1 = phi + h * fabsl((map[i][j] - level)) / (fabsl((map[i][j] - level)) + fabsl((map[i + 1][j] - level)));
                    theta1 = theta;

                    phi2 = phi + h;
                    theta2 = theta + h * fabsl((map[i + 1][j] - level)) / (fabsl((map[i + 1][j] - level)) + fabsl((map[i + 1][j + 1] - level)));

                    l += s2(phi1, phi2, theta1, theta2);

                } else if ((map[i][j + 1] - level) * (map[i + 1][j + 1] - level) < 0.0L) {

                    phi1 = phi + h * fabsl((map[i][j] - level)) / (fabsl((map[i][j] - level)) + fabsl((map[i + 1][j] - level)));
                    theta1 = theta;

                    phi2 = phi + h * fabsl((map[i][j + 1] - level)) / (fabsl((map[i][j + 1] - level)) + fabsl((map[i + 1][j + 1] - level)));
                    theta2 = theta + h;

                    l += s2(phi1, phi2, theta1, theta2);
                }
            } else if ((map[i + 1][j] - level) * (map[i + 1][j + 1] - level) < 0.0L) {

                if ((map[i][j + 1] - level) * (map[i + 1][j + 1] - level) < 0.0L) {

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
        return 1;
    }
}

int condition_2(long double qx, long double qy, long double ux, long double uy) {

    long double root1;
    long double root2;

    long double mean_root1;
    long double mean_root2;

    if (det(qx, qy, ux, uy) < 0) {
        return 1;
    } else {

        d_solver(root1, root2, qx, qy, ux, uy);
        mean_root1 = mean_t_solver(root1, qx, qy, ux, uy);
        mean_root2 = mean_t_solver(root2, qx, qy, ux, uy);

        if ((mean_root1 > 0.0L and mean_root2 < 0.0L) or (mean_root1 < 0.0L and  mean_root2 > 0.0L)) {
            return 2;
        } else {
            return 3;
        }
    }
}

void points_classifier(long double** map_x, long double** map_y, long double** cos_lm, long double** sin_lm,
                       std::string name) {

    std::ofstream out_file;
    out_file.open(name);
    typedef std::numeric_limits<long double> dbl;
    out_file.precision(dbl::max_digits10);

    long double h = 2.0L * PI / long_npix;

    long double phi;
    long double theta;

    int z_x = 0;
    int z_y = 0;
    int flag = 0;

    long double phi1a = 0.0L;
    long double phi2a = 0.0L;
    long double theta1a = 0.0L;
    long double theta2a = 0.0L;
    long double phi1b = 0.0L;
    long double phi2b = 0.0L;
    long double theta1b = 0.0L;
    long double theta2b = 0.0L;

    long double phi_answer;
    long double theta_answer;

    long double fxx;
    long double fyy;
    long double fxy;

    int condition_answer = 0;

    for (unsigned int i = 0; i < npix; ++i) {
        for (unsigned int j = 1; j < npix / 2 - 1; ++j) {

            z_x = 0;
            z_y = 0;

            phi = 2.0L * PI * i / long_npix;
            theta = 2.0L * PI * j / long_npix;

            if (map_x[i][j] * map_x[i][j + 1] < 0.0L) {

                if (map_x[i][j] * map_x[i + 1][j] < 0.0L) {

                    phi1a = phi;
                    theta1a = theta + h * fabsl(map_x[i][j]) / (fabsl(map_x[i][j]) + fabsl(map_x[i][j + 1]));

                    phi1b = phi + h * fabsl(map_x[i][j]) / (fabsl(map_x[i][j]) + fabsl(map_x[i + 1][j]));
                    theta1b = theta;

                    z_x = 1;

                } else if (map_x[i + 1][j] * map_x[i + 1][j + 1] < 0.0L) {

                    phi1a = phi;
                    theta1a = theta + h * fabsl(map_x[i][j]) / (fabsl(map_x[i][j]) + fabsl(map_x[i][j + 1]));

                    phi1b = phi + h;
                    theta1b = theta + h * fabsl(map_x[i + 1][j]) / (fabsl(map_x[i + 1][j]) + fabsl(map_x[i + 1][j + 1]));

                    z_x = 1;

                } else if (map_x[i][j + 1] * map_x[i + 1][j + 1] < 0.0L) {

                    phi1a = phi;
                    theta1a = theta + h * fabsl(map_x[i][j]) / (fabsl(map_x[i][j]) + fabsl(map_x[i][j + 1]));

                    phi1b = phi + h * fabsl(map_x[i][j + 1]) / (fabsl(map_x[i][j + 1]) + fabsl(map_x[i + 1][j + 1]));
                    theta1b = theta + h;

                    z_x = 1;

                }
            } else if (map_x[i][j] * map_x[i + 1][j] < 0.0L) {

                if (map_x[i + 1][j] * map_x[i + 1][j + 1] < 0.0L) {

                    phi1a = phi + h * fabsl(map_x[i][j]) / (fabsl(map_x[i][j]) + fabsl(map_x[i + 1][j]));
                    theta1a = theta;

                    phi1b = phi + h;
                    theta1b = theta + h * fabsl(map_x[i + 1][j]) / (fabsl(map_x[i + 1][j]) + fabsl(map_x[i + 1][j + 1]));

                    z_x = 1;

                } else if (map_x[i][j + 1] * map_x[i + 1][j + 1] < 0.0L) {

                    phi1a = phi + h * fabsl(map_x[i][j]) / (fabsl(map_x[i][j]) + fabsl(map_x[i + 1][j]));
                    theta1a = theta;

                    phi1b = phi + h * fabsl(map_x[i][j + 1]) / (fabsl(map_x[i][j + 1]) + fabsl(map_x[i + 1][j + 1]));
                    theta1b = theta + h;

                    z_x = 1;
                }
            } else if (map_x[i + 1][j] * map_x[i + 1][j + 1] < 0.0L) {

                if (map_x[i][j + 1] * map_x[i + 1][j + 1] < 0.0L) {

                    phi1a = phi + h;
                    theta1a = theta + h * fabsl(map_x[i + 1][j]) / (fabsl(map_x[i + 1][j]) + fabsl(map_x[i + 1][j + 1]));

                    phi1b = phi + h * fabsl(map_x[i][j + 1]) / (fabsl(map_x[i][j + 1]) + fabsl(map_x[i + 1][j + 1]));
                    theta1b = theta + h;

                    z_x = 1;
                }
            }

            if (map_y[i][j] * map_y[i][j + 1] < 0.0L) {

                if (map_y[i][j] * map_y[i + 1][j] < 0.0L) {

                    phi2a = phi;
                    theta2a = theta + h * fabsl(map_y[i][j]) / (fabsl(map_y[i][j]) + fabsl(map_y[i][j + 1]));

                    phi2b = phi + h * fabsl(map_y[i][j]) / (fabsl(map_y[i][j]) + fabsl(map_y[i + 1][j]));
                    theta2b = theta;

                    z_y = 1;

                } else if (map_y[i + 1][j] * map_y[i + 1][j + 1] < 0.0L) {

                    phi2a = phi;
                    theta2a = theta + h * fabsl(map_y[i][j]) / (fabsl(map_y[i][j]) + fabsl(map_y[i][j + 1]));

                    phi2b = phi + h;
                    theta2b = theta + h * fabsl(map_y[i + 1][j]) / (fabsl(map_y[i + 1][j]) + fabsl(map_y[i + 1][j + 1]));

                    z_y = 1;

                } else if (map_y[i][j + 1] * map_y[i + 1][j + 1] < 0.0L) {

                    phi2a = phi;
                    theta2a = theta + h * fabsl(map_y[i][j]) / (fabsl(map_y[i][j]) + fabsl(map_y[i][j + 1]));

                    phi2b = phi + h * fabsl(map_y[i][j + 1]) / (fabsl(map_y[i][j + 1]) + fabsl(map_y[i + 1][j + 1]));
                    theta2b = theta + h;

                    z_y = 1;

                }
            } else if (map_y[i][j] * map_y[i + 1][j] < 0.0L) {

                if (map_y[i + 1][j] * map_y[i + 1][j + 1] < 0.0L) {

                    phi2a = phi + h * fabsl(map_y[i][j]) / (fabsl(map_y[i][j]) + fabsl(map_y[i + 1][j]));
                    theta2a = theta;

                    phi2b = phi + h;
                    theta2b = theta + h * fabsl(map_y[i + 1][j]) / (fabsl(map_y[i + 1][j]) + fabsl(map_y[i + 1][j + 1]));

                    z_y = 1;

                } else if (map_y[i][j + 1] * map_y[i + 1][j + 1] < 0.0L) {

                    phi2a = phi + h * fabsl(map_y[i][j]) / (fabsl(map_y[i][j]) + fabsl(map_y[i + 1][j]));
                    theta2a = theta;

                    phi2b = phi + h * fabsl(map_y[i][j + 1]) / (fabsl(map_y[i][j + 1]) + fabsl(map_y[i + 1][j + 1]));
                    theta2b = theta + h;

                    z_y = 1;
                }
            } else if (map_y[i + 1][j] * map_y[i + 1][j + 1] < 0.0L) {

                if (map_y[i][j + 1] * map_y[i + 1][j + 1] < 0.0L) {

                    phi2a = phi + h;
                    theta2a = theta + h * fabsl(map_y[i + 1][j]) / (fabsl(map_y[i + 1][j]) + fabsl(map_y[i + 1][j + 1]));

                    phi2b = phi + h * fabsl(map_y[i][j + 1]) / (fabsl(map_y[i][j + 1]) + fabsl(map_y[i + 1][j + 1]));
                    theta2b = theta + h;

                    z_y = 1;
                }
            }

            if (z_x == 1 and z_y == 1) {

                flag = 0;
                phi_answer = 0.0L;
                theta_answer = 0.0L;

                cross(phi_answer, theta_answer, phi1a, theta1a, phi1b, theta1b, phi2a, theta2a, phi2b, theta2b);

                if ((phi <= phi_answer) and (phi_answer <= phi + h) and (theta <= theta_answer) and (theta_answer <= theta + h)) {
                    flag = 1;
                } else if ((phi <= phi_answer) and (phi_answer <= phi + h) and (theta <= - theta_answer + PI) and (- theta_answer + PI <= theta + h)) {
                    theta_answer = - theta_answer + PI;
                    flag = 1;
                } else if ((phi <= phi_answer + PI) and (phi_answer + PI <= phi + h) and (theta <= - theta_answer + PI) and (- theta_answer + PI <= theta + h)) {
                    phi_answer = phi_answer + PI;
                    theta_answer = - theta_answer + PI;
                    flag = 1;
                } else if ((phi <= phi_answer + PI) and (phi_answer + PI <= phi + h) and (theta <= theta_answer) and (theta_answer <= theta + h)) {
                    phi_answer = phi_answer + PI;
                    flag = 1;
                } else if ((phi <= phi_answer + 2 * PI) and (phi_answer + 2 * PI <= phi + h) and (theta <= - theta_answer + PI) and (- theta_answer + PI <= theta + h)) {
                    phi_answer = phi_answer + 2 * PI;
                    theta_answer = - theta_answer + PI;
                    flag = 1;
                } else if ((phi <= phi_answer + 2 * PI) and (phi_answer + 2 * PI <= phi + h) and (theta <= theta_answer) and (theta_answer <= theta + h)) {
                    phi_answer = phi_answer + 2 * PI;
                    flag = 1;
                }

                if (flag == 1) {

                    fxx = fft_point_xx_forward(phi_answer, theta_answer, cos_lm, sin_lm);
                    fyy = fft_point_yy_forward(phi_answer, theta_answer, cos_lm, sin_lm);
                    fxy = fft_point_xy_forward(phi_answer, theta_answer, cos_lm, sin_lm);

                    condition_answer = condition_1(fxx, fyy, fxy);

                    out_file << std::scientific << phi_answer << " "
                             << std::scientific << theta_answer << " "
                             << condition_answer << std::endl;
                }
            }
        }
    }

    out_file.close();
}

void singular_points_classifier(long double** q, long double** u, long double** q_cos_lm, long double** q_sin_lm,
                                long double** u_cos_lm, long double** u_sin_lm, std::string name) {
    std::ofstream out_file;
    out_file.open(name);
    typedef std::numeric_limits<long double> dbl;
    out_file.precision(dbl::max_digits10);

    long double h = 2.0L * PI / long_npix;

    long double phi;
    long double theta;

    int z_x = 0;
    int z_y = 0;
    int flag = 0;

    long double phi1a = 0.0L;
    long double phi2a = 0.0L;
    long double theta1a = 0.0L;
    long double theta2a = 0.0L;
    long double phi1b = 0.0L;
    long double phi2b = 0.0L;
    long double theta1b = 0.0L;
    long double theta2b = 0.0L;

    long double phi_answer;
    long double theta_answer;

    long double qx;
    long double qy;
    long double ux;
    long double uy;

    int condition_answer = 0;

    for (unsigned int i = 0; i < npix; ++i) {
        for (unsigned int j = 1; j < npix / 2 - 1; ++j) {

            z_x = 0;
            z_y = 0;

            phi = 2.0L * PI * i / long_npix;
            theta = 2.0L * PI * j / long_npix;

            if (q[i][j] * q[i][j + 1] < 0.0L) {

                if (q[i][j] * q[i + 1][j] < 0.0L) {

                    phi1a = phi;
                    theta1a = theta + h * fabsl(q[i][j]) / (fabsl(q[i][j]) + fabsl(q[i][j + 1]));

                    phi1b = phi + h * fabsl(q[i][j]) / (fabsl(q[i][j]) + fabsl(q[i + 1][j]));
                    theta1b = theta;

                    z_x = 1;

                } else if (q[i + 1][j] * q[i + 1][j + 1] < 0.0L) {

                    phi1a = phi;
                    theta1a = theta + h * fabsl(q[i][j]) / (fabsl(q[i][j]) + fabsl(q[i][j + 1]));

                    phi1b = phi + h;
                    theta1b = theta + h * fabsl(q[i + 1][j]) / (fabsl(q[i + 1][j]) + fabsl(q[i + 1][j + 1]));

                    z_x = 1;

                } else if (q[i][j + 1] * q[i + 1][j + 1] < 0.0L) {

                    phi1a = phi;
                    theta1a = theta + h * fabsl(q[i][j]) / (fabsl(q[i][j]) + fabsl(q[i][j + 1]));

                    phi1b = phi + h * fabsl(q[i][j + 1]) / (fabsl(q[i][j + 1]) + fabsl(q[i + 1][j + 1]));
                    theta1b = theta + h;

                    z_x = 1;

                }
            } else if (q[i][j] * q[i + 1][j] < 0.0L) {

                if (q[i + 1][j] * q[i + 1][j + 1] < 0.0L) {

                    phi1a = phi + h * fabsl(q[i][j]) / (fabsl(q[i][j]) + fabsl(q[i + 1][j]));
                    theta1a = theta;

                    phi1b = phi + h;
                    theta1b = theta + h * fabsl(q[i + 1][j]) / (fabsl(q[i + 1][j]) + fabsl(q[i + 1][j + 1]));

                    z_x = 1;

                } else if (q[i][j + 1] * q[i + 1][j + 1] < 0.0L) {

                    phi1a = phi + h * fabsl(q[i][j]) / (fabsl(q[i][j]) + fabsl(q[i + 1][j]));
                    theta1a = theta;

                    phi1b = phi + h * fabsl(q[i][j + 1]) / (fabsl(q[i][j + 1]) + fabsl(q[i + 1][j + 1]));
                    theta1b = theta + h;

                    z_x = 1;
                }
            } else if (q[i + 1][j] * q[i + 1][j + 1] < 0.0L) {

                if (q[i][j + 1] * q[i + 1][j + 1] < 0.0L) {

                    phi1a = phi + h;
                    theta1a = theta + h * fabsl(q[i + 1][j]) / (fabsl(q[i + 1][j]) + fabsl(q[i + 1][j + 1]));

                    phi1b = phi + h * fabsl(q[i][j + 1]) / (fabsl(q[i][j + 1]) + fabsl(q[i + 1][j + 1]));
                    theta1b = theta + h;

                    z_x = 1;
                }
            }

            if (u[i][j] * u[i][j + 1] < 0.0L) {

                if (u[i][j] * u[i + 1][j] < 0.0L) {

                    phi2a = phi;
                    theta2a = theta + h * fabsl(u[i][j]) / (fabsl(u[i][j]) + fabsl(u[i][j + 1]));

                    phi2b = phi + h * fabsl(u[i][j]) / (fabsl(u[i][j]) + fabsl(u[i + 1][j]));
                    theta2b = theta;

                    z_y = 1;

                } else if (u[i + 1][j] * u[i + 1][j + 1] < 0.0L) {

                    phi2a = phi;
                    theta2a = theta + h * fabsl(u[i][j]) / (fabsl(u[i][j]) + fabsl(u[i][j + 1]));

                    phi2b = phi + h;
                    theta2b = theta + h * fabsl(u[i + 1][j]) / (fabsl(u[i + 1][j]) + fabsl(u[i + 1][j + 1]));

                    z_y = 1;

                } else if (u[i][j + 1] * u[i + 1][j + 1] < 0.0L) {

                    phi2a = phi;
                    theta2a = theta + h * fabsl(u[i][j]) / (fabsl(u[i][j]) + fabsl(u[i][j + 1]));

                    phi2b = phi + h * fabsl(u[i][j + 1]) / (fabsl(u[i][j + 1]) + fabsl(u[i + 1][j + 1]));
                    theta2b = theta + h;

                    z_y = 1;

                }
            } else if (u[i][j] * u[i + 1][j] < 0.0L) {

                if (u[i + 1][j] * u[i + 1][j + 1] < 0.0L) {

                    phi2a = phi + h * fabsl(u[i][j]) / (fabsl(u[i][j]) + fabsl(u[i + 1][j]));
                    theta2a = theta;

                    phi2b = phi + h;
                    theta2b = theta + h * fabsl(u[i + 1][j]) / (fabsl(u[i + 1][j]) + fabsl(u[i + 1][j + 1]));

                    z_y = 1;

                } else if (u[i][j + 1] * u[i + 1][j + 1] < 0.0L) {

                    phi2a = phi + h * fabsl(u[i][j]) / (fabsl(u[i][j]) + fabsl(u[i + 1][j]));
                    theta2a = theta;

                    phi2b = phi + h * fabsl(u[i][j + 1]) / (fabsl(u[i][j + 1]) + fabsl(u[i + 1][j + 1]));
                    theta2b = theta + h;

                    z_y = 1;
                }
            } else if (u[i + 1][j] * u[i + 1][j + 1] < 0.0L) {

                if (u[i][j + 1] * u[i + 1][j + 1] < 0.0L) {

                    phi2a = phi + h;
                    theta2a = theta + h * fabsl(u[i + 1][j]) / (fabsl(u[i + 1][j]) + fabsl(u[i + 1][j + 1]));

                    phi2b = phi + h * fabsl(u[i][j + 1]) / (fabsl(u[i][j + 1]) + fabsl(u[i + 1][j + 1]));
                    theta2b = theta + h;

                    z_y = 1;
                }
            }

            if (z_x == 1 and z_y == 1) {

                flag = 0;
                phi_answer = 0.0L;
                theta_answer = 0.0L;

                cross(phi_answer, theta_answer, phi1a, theta1a, phi1b, theta1b, phi2a, theta2a, phi2b, theta2b);

                if ((phi <= phi_answer) and (phi_answer <= phi + h) and (theta <= theta_answer) and (theta_answer <= theta + h)) {
                    flag = 1;
                } else if ((phi <= phi_answer) and (phi_answer <= phi + h) and (theta <= - theta_answer + PI) and (- theta_answer + PI <= theta + h)) {
                    theta_answer = - theta_answer + PI;
                    flag = 1;
                } else if ((phi <= phi_answer + PI) and (phi_answer + PI <= phi + h) and (theta <= - theta_answer + PI) and (- theta_answer + PI <= theta + h)) {
                    phi_answer = phi_answer + PI;
                    theta_answer = - theta_answer + PI;
                    flag = 1;
                } else if ((phi <= phi_answer + PI) and (phi_answer + PI <= phi + h) and (theta <= theta_answer) and (theta_answer <= theta + h)) {
                    phi_answer = phi_answer + PI;
                    flag = 1;
                } else if ((phi <= phi_answer + 2 * PI) and (phi_answer + 2 * PI <= phi + h) and (theta <= - theta_answer + PI) and (- theta_answer + PI <= theta + h)) {
                    phi_answer = phi_answer + 2 * PI;
                    theta_answer = - theta_answer + PI;
                    flag = 1;
                } else if ((phi <= phi_answer + 2 * PI) and (phi_answer + 2 * PI <= phi + h) and (theta <= theta_answer) and (theta_answer <= theta + h)) {
                    phi_answer = phi_answer + 2 * PI;
                    flag = 1;
                }

                if (flag == 1) {

                    qx = fft_point_x_forward_n(phi_answer, theta_answer, q_cos_lm, q_sin_lm);
                    qy = fft_point_y_forward_n(phi_answer, theta_answer, q_cos_lm, q_sin_lm);
                    ux = fft_point_x_forward_n(phi_answer, theta_answer, u_cos_lm, u_sin_lm);
                    uy = fft_point_y_forward_n(phi_answer, theta_answer, u_cos_lm, u_sin_lm);

                    condition_answer = condition_2(qx, qy, ux, uy);

                    out_file << std::scientific << phi_answer << " "
                             << std::scientific << theta_answer << " "
                             << condition_answer << std::endl;
                }
            }
        }
    }

    out_file.close();
}
