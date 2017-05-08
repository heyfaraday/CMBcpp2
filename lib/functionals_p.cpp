#include "parameters.hpp"
#include <cmath>
#include <string>
#include <fstream>
#include "constants.hpp"
#include "distance.hpp"
#include "functionals.hpp"

#include "functionals_p.hpp"

void points_classifier_p(long double** map, long double** map_x, long double** map_y, long double** map_xx,
                         long double** map_yy, long double** map_xy, std::string name, int** whitelist,
                         long double sigma_p) {

    FILE* fp = fopen(name.c_str(), "w");

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

    long double p_answer;

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
                    theta1a = theta + map_parameter * fabsl(map_x[i][j])
                                      / (fabsl(map_x[i][j]) + fabsl(map_x[i][j + 1]));

                    phi1b = phi + map_parameter * fabsl(map_x[i][j])
                                  / (fabsl(map_x[i][j]) + fabsl(map_x[i + 1][j]));
                    theta1b = theta;

                    z_x = 1;

                } else if (map_x[i + 1][j] * map_x[i + 1][j + 1] < 0.0L) {

                    phi1a = phi;
                    theta1a = theta + map_parameter * fabsl(map_x[i][j])
                                      / (fabsl(map_x[i][j]) + fabsl(map_x[i][j + 1]));

                    phi1b = phi + map_parameter;
                    theta1b = theta + map_parameter * fabsl(map_x[i + 1][j])
                                      / (fabsl(map_x[i + 1][j]) + fabsl(map_x[i + 1][j + 1]));

                    z_x = 1;

                } else if (map_x[i][j + 1] * map_x[i + 1][j + 1] < 0.0L) {

                    phi1a = phi;
                    theta1a = theta + map_parameter * fabsl(map_x[i][j])
                                      / (fabsl(map_x[i][j]) + fabsl(map_x[i][j + 1]));

                    phi1b = phi + map_parameter * fabsl(map_x[i][j + 1])
                                  / (fabsl(map_x[i][j + 1]) + fabsl(map_x[i + 1][j + 1]));
                    theta1b = theta + map_parameter;

                    z_x = 1;

                }
            } else if (map_x[i][j] * map_x[i + 1][j] < 0.0L) {

                if (map_x[i + 1][j] * map_x[i + 1][j + 1] < 0.0L) {

                    phi1a = phi + map_parameter * fabsl(map_x[i][j])
                                  / (fabsl(map_x[i][j]) + fabsl(map_x[i + 1][j]));
                    theta1a = theta;

                    phi1b = phi + map_parameter;
                    theta1b = theta + map_parameter * fabsl(map_x[i + 1][j])
                                      / (fabsl(map_x[i + 1][j]) + fabsl(map_x[i + 1][j + 1]));

                    z_x = 1;

                } else if (map_x[i][j + 1] * map_x[i + 1][j + 1] < 0.0L) {

                    phi1a = phi + map_parameter * fabsl(map_x[i][j])
                                  / (fabsl(map_x[i][j]) + fabsl(map_x[i + 1][j]));
                    theta1a = theta;

                    phi1b = phi + map_parameter * fabsl(map_x[i][j + 1])
                                  / (fabsl(map_x[i][j + 1]) + fabsl(map_x[i + 1][j + 1]));
                    theta1b = theta + map_parameter;

                    z_x = 1;
                }
            } else if (map_x[i + 1][j] * map_x[i + 1][j + 1] < 0.0L) {

                if (map_x[i][j + 1] * map_x[i + 1][j + 1] < 0.0L) {

                    phi1a = phi + map_parameter;
                    theta1a = theta + map_parameter * fabsl(map_x[i + 1][j])
                                      / (fabsl(map_x[i + 1][j]) + fabsl(map_x[i + 1][j + 1]));

                    phi1b = phi + map_parameter * fabsl(map_x[i][j + 1])
                                  / (fabsl(map_x[i][j + 1]) + fabsl(map_x[i + 1][j + 1]));
                    theta1b = theta + map_parameter;

                    z_x = 1;
                }
            }

            if (map_y[i][j] * map_y[i][j + 1] < 0.0L) {

                if (map_y[i][j] * map_y[i + 1][j] < 0.0L) {

                    phi2a = phi;
                    theta2a = theta + map_parameter * fabsl(map_y[i][j])
                                      / (fabsl(map_y[i][j]) + fabsl(map_y[i][j + 1]));

                    phi2b = phi + map_parameter * fabsl(map_y[i][j])
                                  / (fabsl(map_y[i][j]) + fabsl(map_y[i + 1][j]));
                    theta2b = theta;

                    z_y = 1;

                } else if (map_y[i + 1][j] * map_y[i + 1][j + 1] < 0.0L) {

                    phi2a = phi;
                    theta2a = theta + map_parameter * fabsl(map_y[i][j])
                                      / (fabsl(map_y[i][j]) + fabsl(map_y[i][j + 1]));

                    phi2b = phi + map_parameter;
                    theta2b = theta + map_parameter * fabsl(map_y[i + 1][j])
                                      / (fabsl(map_y[i + 1][j]) + fabsl(map_y[i + 1][j + 1]));

                    z_y = 1;

                } else if (map_y[i][j + 1] * map_y[i + 1][j + 1] < 0.0L) {

                    phi2a = phi;
                    theta2a = theta + map_parameter * fabsl(map_y[i][j])
                                      / (fabsl(map_y[i][j]) + fabsl(map_y[i][j + 1]));

                    phi2b = phi + map_parameter * fabsl(map_y[i][j + 1])
                                  / (fabsl(map_y[i][j + 1]) + fabsl(map_y[i + 1][j + 1]));
                    theta2b = theta + map_parameter;

                    z_y = 1;

                }
            } else if (map_y[i][j] * map_y[i + 1][j] < 0.0L) {

                if (map_y[i + 1][j] * map_y[i + 1][j + 1] < 0.0L) {

                    phi2a = phi + map_parameter * fabsl(map_y[i][j])
                                  / (fabsl(map_y[i][j]) + fabsl(map_y[i + 1][j]));
                    theta2a = theta;

                    phi2b = phi + map_parameter;
                    theta2b = theta + map_parameter * fabsl(map_y[i + 1][j])
                                      / (fabsl(map_y[i + 1][j]) + fabsl(map_y[i + 1][j + 1]));

                    z_y = 1;

                } else if (map_y[i][j + 1] * map_y[i + 1][j + 1] < 0.0L) {

                    phi2a = phi + map_parameter * fabsl(map_y[i][j])
                                  / (fabsl(map_y[i][j]) + fabsl(map_y[i + 1][j]));
                    theta2a = theta;

                    phi2b = phi + map_parameter * fabsl(map_y[i][j + 1])
                                  / (fabsl(map_y[i][j + 1]) + fabsl(map_y[i + 1][j + 1]));
                    theta2b = theta + map_parameter;

                    z_y = 1;
                }
            } else if (map_y[i + 1][j] * map_y[i + 1][j + 1] < 0.0L) {

                if (map_y[i][j + 1] * map_y[i + 1][j + 1] < 0.0L) {

                    phi2a = phi + map_parameter;
                    theta2a = theta + map_parameter * fabsl(map_y[i + 1][j])
                                      / (fabsl(map_y[i + 1][j]) + fabsl(map_y[i + 1][j + 1]));

                    phi2b = phi + map_parameter * fabsl(map_y[i][j + 1])
                                  / (fabsl(map_y[i][j + 1]) + fabsl(map_y[i + 1][j + 1]));
                    theta2b = theta + map_parameter;

                    z_y = 1;
                }
            }

            if (z_x == 1 and z_y == 1 and whitelist[i][j] != 1) {

                flag = 0;
                phi_answer = 0.0L;
                theta_answer = 0.0L;

                cross(phi_answer, theta_answer, phi1a, theta1a, phi1b, theta1b, phi2a, theta2a, phi2b, theta2b);

                if ((phi <= phi_answer) and (phi_answer <= phi + map_parameter) and (theta <= theta_answer)
                    and (theta_answer <= theta + map_parameter)) {
                    flag = 1;
                } else if ((phi <= phi_answer) and (phi_answer <= phi + map_parameter)
                           and (theta <= - theta_answer + PI) and (- theta_answer + PI <= theta + map_parameter)) {
                    theta_answer = - theta_answer + PI;
                    flag = 1;
                } else if ((phi <= phi_answer + PI) and (phi_answer + PI <= phi + map_parameter)
                           and (theta <= - theta_answer + PI) and (- theta_answer + PI <= theta + map_parameter)) {
                    phi_answer = phi_answer + PI;
                    theta_answer = - theta_answer + PI;
                    flag = 1;
                } else if ((phi <= phi_answer + PI) and (phi_answer + PI <= phi + map_parameter)
                           and (theta <= theta_answer) and (theta_answer <= theta + map_parameter)) {
                    phi_answer = phi_answer + PI;
                    flag = 1;
                } else if ((phi <= phi_answer + 2 * PI) and (phi_answer + 2 * PI <= phi + map_parameter)
                           and (theta <= - theta_answer + PI) and (- theta_answer + PI <= theta + map_parameter)) {
                    phi_answer = phi_answer + 2 * PI;
                    theta_answer = - theta_answer + PI;
                    flag = 1;
                } else if ((phi <= phi_answer + 2 * PI) and (phi_answer + 2 * PI <= phi + map_parameter)
                           and (theta <= theta_answer) and (theta_answer <= theta + map_parameter)) {
                    phi_answer = phi_answer + 2 * PI;
                    flag = 1;
                }

                if (flag == 1) {

                    fxx = (map_xx[i][j] * (phi + map_parameter - phi_answer) / (map_parameter) + map_xx[i + 1][j] * (phi_answer - phi) / (map_parameter))
                          * (theta + map_parameter - theta_answer) / (map_parameter)
                          + (map_xx[i][j + 1] * (phi + map_parameter - phi_answer) / (map_parameter) + map_xx[i + 1][j + 1] * (phi_answer - phi) / (map_parameter))
                            * (theta_answer - theta) / map_parameter;


                    fyy = (map_yy[i][j] * (phi + map_parameter - phi_answer) / (map_parameter) + map_yy[i + 1][j] * (phi_answer - phi) / (map_parameter))
                          * (theta + map_parameter - theta_answer) / (map_parameter)
                          + (map_yy[i][j + 1] * (phi + map_parameter - phi_answer) / (map_parameter) + map_yy[i + 1][j + 1] * (phi_answer - phi) / (map_parameter))
                            * (theta_answer - theta) / map_parameter;

                    fxy = (map_xy[i][j] * (phi + map_parameter - phi_answer) / (map_parameter) + map_xy[i + 1][j] * (phi_answer - phi) / (map_parameter))
                          * (theta + map_parameter - theta_answer) / (map_parameter)
                          + (map_xy[i][j + 1] * (phi + map_parameter - phi_answer) / (map_parameter) + map_xy[i + 1][j + 1] * (phi_answer - phi) / (map_parameter))
                            * (theta_answer - theta) / map_parameter;


                    p_answer = (map[i][j] * (phi + map_parameter - phi_answer) / (map_parameter) + map[i + 1][j] * (phi_answer - phi) / (map_parameter))
                               * (theta + map_parameter - theta_answer) / (map_parameter)
                               + (map[i][j + 1] * (phi + map_parameter - phi_answer) / (map_parameter) + map[i + 1][j + 1] * (phi_answer - phi) / (map_parameter))
                                 * (theta_answer - theta) / map_parameter;

                    condition_answer = condition_1(fxx, fyy, fxy);

                    fprintf(fp, "%.21Le %.21Le %.21Le %u \n", phi_answer, theta_answer, p_answer / sigma_p, condition_answer);

                }
            }
        }
    }

    fclose(fp);
}

void singular_points_classifier(long double** q, long double** u, long double** qx, long double** ux,
                                long double** qy, long double** uy, std::string name, int** whitelist) {

    FILE* fp = fopen(name.c_str(), "w");

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

    long double qx_answer;
    long double qy_answer;
    long double ux_answer;
    long double uy_answer;

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
                    theta1a = theta + map_parameter * fabsl(q[i][j])
                                      / (fabsl(q[i][j]) + fabsl(q[i][j + 1]));

                    phi1b = phi + map_parameter * fabsl(q[i][j])
                                  / (fabsl(q[i][j]) + fabsl(q[i + 1][j]));
                    theta1b = theta;

                    z_x = 1;

                } else if (q[i + 1][j] * q[i + 1][j + 1] < 0.0L) {

                    phi1a = phi;
                    theta1a = theta + map_parameter * fabsl(q[i][j])
                                      / (fabsl(q[i][j]) + fabsl(q[i][j + 1]));

                    phi1b = phi + map_parameter;
                    theta1b = theta + map_parameter * fabsl(q[i + 1][j])
                                      / (fabsl(q[i + 1][j]) + fabsl(q[i + 1][j + 1]));

                    z_x = 1;

                } else if (q[i][j + 1] * q[i + 1][j + 1] < 0.0L) {

                    phi1a = phi;
                    theta1a = theta + map_parameter * fabsl(q[i][j])
                                      / (fabsl(q[i][j]) + fabsl(q[i][j + 1]));

                    phi1b = phi + map_parameter * fabsl(q[i][j + 1])
                                  / (fabsl(q[i][j + 1]) + fabsl(q[i + 1][j + 1]));
                    theta1b = theta + map_parameter;

                    z_x = 1;

                }
            } else if (q[i][j] * q[i + 1][j] < 0.0L) {

                if (q[i + 1][j] * q[i + 1][j + 1] < 0.0L) {

                    phi1a = phi + map_parameter * fabsl(q[i][j])
                                  / (fabsl(q[i][j]) + fabsl(q[i + 1][j]));
                    theta1a = theta;

                    phi1b = phi + map_parameter;
                    theta1b = theta + map_parameter * fabsl(q[i + 1][j])
                                      / (fabsl(q[i + 1][j]) + fabsl(q[i + 1][j + 1]));

                    z_x = 1;

                } else if (q[i][j + 1] * q[i + 1][j + 1] < 0.0L) {

                    phi1a = phi + map_parameter * fabsl(q[i][j])
                                  / (fabsl(q[i][j]) + fabsl(q[i + 1][j]));
                    theta1a = theta;

                    phi1b = phi + map_parameter * fabsl(q[i][j + 1])
                                  / (fabsl(q[i][j + 1]) + fabsl(q[i + 1][j + 1]));
                    theta1b = theta + map_parameter;

                    z_x = 1;
                }
            } else if (q[i + 1][j] * q[i + 1][j + 1] < 0.0L) {

                if (q[i][j + 1] * q[i + 1][j + 1] < 0.0L) {

                    phi1a = phi + map_parameter;
                    theta1a = theta + map_parameter * fabsl(q[i + 1][j])
                                      / (fabsl(q[i + 1][j]) + fabsl(q[i + 1][j + 1]));

                    phi1b = phi + map_parameter * fabsl(q[i][j + 1])
                                  / (fabsl(q[i][j + 1]) + fabsl(q[i + 1][j + 1]));
                    theta1b = theta + map_parameter;

                    z_x = 1;
                }
            }

            if (u[i][j] * u[i][j + 1] < 0.0L) {

                if (u[i][j] * u[i + 1][j] < 0.0L) {

                    phi2a = phi;
                    theta2a = theta + map_parameter * fabsl(u[i][j])
                                      / (fabsl(u[i][j]) + fabsl(u[i][j + 1]));

                    phi2b = phi + map_parameter * fabsl(u[i][j])
                                  / (fabsl(u[i][j]) + fabsl(u[i + 1][j]));
                    theta2b = theta;

                    z_y = 1;

                } else if (u[i + 1][j] * u[i + 1][j + 1] < 0.0L) {

                    phi2a = phi;
                    theta2a = theta + map_parameter * fabsl(u[i][j])
                                      / (fabsl(u[i][j]) + fabsl(u[i][j + 1]));

                    phi2b = phi + map_parameter;
                    theta2b = theta + map_parameter * fabsl(u[i + 1][j])
                                      / (fabsl(u[i + 1][j]) + fabsl(u[i + 1][j + 1]));

                    z_y = 1;

                } else if (u[i][j + 1] * u[i + 1][j + 1] < 0.0L) {

                    phi2a = phi;
                    theta2a = theta + map_parameter * fabsl(u[i][j])
                                      / (fabsl(u[i][j]) + fabsl(u[i][j + 1]));

                    phi2b = phi + map_parameter * fabsl(u[i][j + 1])
                                  / (fabsl(u[i][j + 1]) + fabsl(u[i + 1][j + 1]));
                    theta2b = theta + map_parameter;

                    z_y = 1;

                }
            } else if (u[i][j] * u[i + 1][j] < 0.0L) {

                if (u[i + 1][j] * u[i + 1][j + 1] < 0.0L) {

                    phi2a = phi + map_parameter * fabsl(u[i][j])
                                  / (fabsl(u[i][j]) + fabsl(u[i + 1][j]));
                    theta2a = theta;

                    phi2b = phi + map_parameter;
                    theta2b = theta + map_parameter * fabsl(u[i + 1][j])
                                      / (fabsl(u[i + 1][j]) + fabsl(u[i + 1][j + 1]));

                    z_y = 1;

                } else if (u[i][j + 1] * u[i + 1][j + 1] < 0.0L) {

                    phi2a = phi + map_parameter * fabsl(u[i][j])
                                  / (fabsl(u[i][j]) + fabsl(u[i + 1][j]));
                    theta2a = theta;

                    phi2b = phi + map_parameter * fabsl(u[i][j + 1])
                                  / (fabsl(u[i][j + 1]) + fabsl(u[i + 1][j + 1]));
                    theta2b = theta + map_parameter;

                    z_y = 1;
                }
            } else if (u[i + 1][j] * u[i + 1][j + 1] < 0.0L) {

                if (u[i][j + 1] * u[i + 1][j + 1] < 0.0L) {

                    phi2a = phi + map_parameter;
                    theta2a = theta + map_parameter * fabsl(u[i + 1][j])
                                      / (fabsl(u[i + 1][j]) + fabsl(u[i + 1][j + 1]));

                    phi2b = phi + map_parameter * fabsl(u[i][j + 1])
                                  / (fabsl(u[i][j + 1]) + fabsl(u[i + 1][j + 1]));
                    theta2b = theta + map_parameter;

                    z_y = 1;
                }
            }

            if (z_x == 1 and z_y == 1) {

                flag = 0;
                phi_answer = 0.0L;
                theta_answer = 0.0L;

                cross(phi_answer, theta_answer, phi1a, theta1a, phi1b, theta1b, phi2a, theta2a, phi2b, theta2b);

                if ((phi <= phi_answer) and (phi_answer <= phi + map_parameter) and (theta <= theta_answer)
                    and (theta_answer <= theta + map_parameter)) {
                    flag = 1;
                } else if ((phi <= phi_answer) and (phi_answer <= phi + map_parameter)
                           and (theta <= - theta_answer + PI) and (- theta_answer + PI <= theta + map_parameter)) {
                    theta_answer = - theta_answer + PI;
                    flag = 1;
                } else if ((phi <= phi_answer + PI) and (phi_answer + PI <= phi + map_parameter)
                           and (theta <= - theta_answer + PI) and (- theta_answer + PI <= theta + map_parameter)) {
                    phi_answer = phi_answer + PI;
                    theta_answer = - theta_answer + PI;
                    flag = 1;
                } else if ((phi <= phi_answer + PI) and (phi_answer + PI <= phi + map_parameter)
                           and (theta <= theta_answer) and (theta_answer <= theta + map_parameter)) {
                    phi_answer = phi_answer + PI;
                    flag = 1;
                } else if ((phi <= phi_answer + 2 * PI) and (phi_answer + 2 * PI <= phi + map_parameter)
                           and (theta <= - theta_answer + PI) and (- theta_answer + PI <= theta + map_parameter)) {
                    phi_answer = phi_answer + 2 * PI;
                    theta_answer = - theta_answer + PI;
                    flag = 1;
                } else if ((phi <= phi_answer + 2 * PI) and (phi_answer + 2 * PI <= phi + map_parameter)
                           and (theta <= theta_answer) and (theta_answer <= theta + map_parameter)) {
                    phi_answer = phi_answer + 2 * PI;
                    flag = 1;
                }

                if (flag == 1) {

                    qx_answer = (qx[i][j] * (phi + map_parameter - phi_answer) / (map_parameter) + qx[i + 1][j] * (phi_answer - phi) / (map_parameter))
                                * (theta + map_parameter - theta_answer) / (map_parameter)
                                + (qx[i][j + 1] * (phi + map_parameter - phi_answer) / (map_parameter) + qx[i + 1][j + 1] * (phi_answer - phi) / (map_parameter))
                                  * (theta_answer - theta) / map_parameter;

                    ux_answer = (ux[i][j] * (phi + map_parameter - phi_answer) / (map_parameter) + ux[i + 1][j] * (phi_answer - phi) / (map_parameter))
                                * (theta + map_parameter - theta_answer) / (map_parameter)
                                + (ux[i][j + 1] * (phi + map_parameter - phi_answer) / (map_parameter) + ux[i + 1][j + 1] * (phi_answer - phi) / (map_parameter))
                                  * (theta_answer - theta) / map_parameter;

                    qy_answer = (qy[i][j] * (phi + map_parameter - phi_answer) / (map_parameter) + qy[i + 1][j] * (phi_answer - phi) / (map_parameter))
                                * (theta + map_parameter - theta_answer) / (map_parameter)
                                + (qy[i][j + 1] * (phi + map_parameter - phi_answer) / (map_parameter) + qy[i + 1][j + 1] * (phi_answer - phi) / (map_parameter))
                                  * (theta_answer - theta) / map_parameter;

                    uy_answer = (uy[i][j] * (phi + map_parameter - phi_answer) / (map_parameter) + uy[i + 1][j] * (phi_answer - phi) / (map_parameter))
                                * (theta + map_parameter - theta_answer) / (map_parameter)
                                + (uy[i][j + 1] * (phi + map_parameter - phi_answer) / (map_parameter) + uy[i + 1][j + 1] * (phi_answer - phi) / (map_parameter))
                                  * (theta_answer - theta) / map_parameter;

                    condition_answer = condition_2(qx_answer, qy_answer, ux_answer, uy_answer);

                    fprintf(fp, "%.21Le %.21Le %u \n", phi_answer, theta_answer, condition_answer);

                    whitelist[i][j] = 1;
                }
            }
        }
    }

    fclose(fp);
}

void new_points_classifier_p(long double** map, long double** map_x, long double** map_y, long double** map_xx,
                         long double** map_yy, long double** map_xy, std::string name, int** whitelist,
                         long double sigma_p) {

    FILE* fp = fopen(name.c_str(), "w");

    long double phi;
    long double theta;

    int z_x = 0;
    int z_y = 0;
    int flag = 0;

    long double phi_answer_1;
    long double theta_answer_1;
    long double phi_answer_2;
    long double theta_answer_2;
    long double p_answer;
    long double under_sqrt;

    long double a00;
    long double a01;
    long double a10;
    long double a11;

    long double b00;
    long double b01;
    long double b10;
    long double b11;

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

                    z_x = 1;

                } else if (map_x[i + 1][j] * map_x[i + 1][j + 1] < 0.0L) {

                    z_x = 1;

                } else if (map_x[i][j + 1] * map_x[i + 1][j + 1] < 0.0L) {

                    z_x = 1;

                }
            } else if (map_x[i][j] * map_x[i + 1][j] < 0.0L) {

                if (map_x[i + 1][j] * map_x[i + 1][j + 1] < 0.0L) {

                    z_x = 1;

                } else if (map_x[i][j + 1] * map_x[i + 1][j + 1] < 0.0L) {

                    z_x = 1;
                }
            } else if (map_x[i + 1][j] * map_x[i + 1][j + 1] < 0.0L) {

                if (map_x[i][j + 1] * map_x[i + 1][j + 1] < 0.0L) {

                    z_x = 1;
                }
            }

            if (map_y[i][j] * map_y[i][j + 1] < 0.0L) {

                if (map_y[i][j] * map_y[i + 1][j] < 0.0L) {

                    z_y = 1;

                } else if (map_y[i + 1][j] * map_y[i + 1][j + 1] < 0.0L) {

                    z_y = 1;

                } else if (map_y[i][j + 1] * map_y[i + 1][j + 1] < 0.0L) {

                    z_y = 1;

                }
            } else if (map_y[i][j] * map_y[i + 1][j] < 0.0L) {

                if (map_y[i + 1][j] * map_y[i + 1][j + 1] < 0.0L) {

                    z_y = 1;

                } else if (map_y[i][j + 1] * map_y[i + 1][j + 1] < 0.0L) {

                    z_y = 1;
                }
            } else if (map_y[i + 1][j] * map_y[i + 1][j + 1] < 0.0L) {

                if (map_y[i][j + 1] * map_y[i + 1][j + 1] < 0.0L) {

                    z_y = 1;
                }
            }

            if (z_x == 1 and z_y == 1 and whitelist[i][j] != 1) {

                flag = 0;
                phi_answer_1 = 0.0L;
                theta_answer_1 = 0.0L;
                phi_answer_2 = 0.0L;
                theta_answer_2 = 0.0L;

                a00 = map_x[i][j];
                a01 = map_x[i][j + 1];
                a10 = map_x[i + 1][j];
                a11 = map_x[i + 1][j + 1];

                b00 = map_y[i][j];
                b01 = map_y[i][j + 1];
                b10 = map_y[i + 1][j];
                b11 = map_y[i + 1][j + 1];

                under_sqrt = a11 * a11 * b00 * b00 + a10 * a10 * b01 * b01 + (a01 * b10 - a00 * b11) * (a01 * b10 - a00 * b11) -
                        2.0l * a11 * (a10 * b00 * b01 + a01 * b00 * b10 - 2.0l * a00 * b01 * b10 + a00 * b00 * b11) -
                        2.0l * a10 * (a01 * b01 * b10 - 2.0l * a01 * b00 * b11 + a00 * b01 * b11);

                if (under_sqrt > 0.0l) {
                    phi_answer_1 = -((-2.0l * a01 * b00 + a11 * b00 + 2.0l * a00 * b01 - a10 * b01 + a01 * b10 -
                                      a00 * b11 + sqrtl(under_sqrt))) / (2.0l * ((a01 - a11) * (b00 - b10) - (a00 - a10) * (b01 - b11)));
                    theta_answer_1 = (2.0l * a10 * b00 - a11 * b00 - a10 * b01 - 2 * a00 * b10 + a01 * b10 +
                                      a00 * b11 + sqrtl(under_sqrt))/(2.0l * ((a10 - a11) * (b00 - b01) - (a00 - a01) * (b10 - b11)));
                    phi_answer_2 = -((-2.0l * a01 * b00 + a11 * b00 + 2.0l * a00 * b01 - a10 * b01 + a01 * b10 -
                                      a00 * b11 - sqrtl(under_sqrt))) / (2.0l * ((a01 - a11) * (b00 - b10) - (a00 - a10) * (b01 - b11)));
                    theta_answer_2 = (2.0l * a10 * b00 - a11 * b00 - a10 * b01 - 2 * a00 * b10 + a01 * b10 +
                                      a00 * b11 - sqrtl(under_sqrt))/(2.0l * ((a10 - a11) * (b00 - b01) - (a00 - a01) * (b10 - b11)));


                    if ((0.0l < phi_answer_1) and (phi_answer_1 < 1.0l) and (0.0l < theta_answer_1)
                        and (theta_answer_1 < 1.0l)) {
                        flag = 1;
                    }

                    if ((0.0l < phi_answer_2) and (phi_answer_2 < 1.0l) and (0.0l < theta_answer_2)
                        and (theta_answer_2 < 1.0l)) {
                        flag = 2;
                    }
                }

                if (flag == 1) {

                    fxx = map_xx[i][j] + phi_answer_1 * (map_xx[i + 1][j] - map_xx[i][j]) + theta_answer_1 * (map_xx[i][j + 1] - map_xx[i][j]) +
                            phi_answer_1 * theta_answer_1 * (map_xx[i][j] - map_xx[i + 1][j] - map_xx[i][j + 1] + map_xx[i + 1][j + 1]);


                    fyy = map_yy[i][j] + phi_answer_1 * (map_yy[i + 1][j] - map_yy[i][j]) + theta_answer_1 * (map_yy[i][j + 1] - map_yy[i][j]) +
                          phi_answer_1 * theta_answer_1 * (map_yy[i][j] - map_yy[i + 1][j] - map_yy[i][j + 1] + map_yy[i + 1][j + 1]);

                    fxy = map_xy[i][j] + phi_answer_1 * (map_xy[i + 1][j] - map_xy[i][j]) + theta_answer_1 * (map_xy[i][j + 1] - map_xy[i][j]) +
                          phi_answer_1 * theta_answer_1 * (map_xy[i][j] - map_xy[i + 1][j] - map_xy[i][j + 1] + map_xy[i + 1][j + 1]);


                    p_answer = map[i][j] + phi_answer_1 * (map[i + 1][j] - map[i][j]) + theta_answer_1 * (map[i][j + 1] - map[i][j]) +
                               phi_answer_1 * theta_answer_1 * (map[i][j] - map[i + 1][j] - map[i][j + 1] + map[i + 1][j + 1]);

                    condition_answer = condition_1(fxx, fyy, fxy);

                    phi_answer_1 = phi + phi_answer_1 * map_parameter;

                    theta_answer_1 = theta + theta_answer_1 * map_parameter;

                    fprintf(fp, "%.21Le %.21Le %.21Le %u \n", phi_answer_1, theta_answer_1, p_answer / sigma_p, condition_answer);

                }

                if (flag == 2) {

                    fxx = map_xx[i][j] + phi_answer_2 * (map_xx[i + 1][j] - map_xx[i][j]) + theta_answer_2 * (map_xx[i][j + 1] - map_xx[i][j]) +
                          phi_answer_2 * theta_answer_2 * (map_xx[i][j] - map_xx[i + 1][j] - map_xx[i][j + 1] + map_xx[i + 1][j + 1]);


                    fyy = map_yy[i][j] + phi_answer_2 * (map_yy[i + 1][j] - map_yy[i][j]) + theta_answer_2 * (map_yy[i][j + 1] - map_yy[i][j]) +
                          phi_answer_2 * theta_answer_2 * (map_yy[i][j] - map_yy[i + 1][j] - map_yy[i][j + 1] + map_yy[i + 1][j + 1]);

                    fxy = map_xy[i][j] + phi_answer_2 * (map_xy[i + 1][j] - map_xy[i][j]) + theta_answer_2 * (map_xy[i][j + 1] - map_xy[i][j]) +
                          phi_answer_2 * theta_answer_2 * (map_xy[i][j] - map_xy[i + 1][j] - map_xy[i][j + 1] + map_xy[i + 1][j + 1]);


                    p_answer = map[i][j] + phi_answer_2 * (map[i + 1][j] - map[i][j]) + theta_answer_2 * (map[i][j + 1] - map[i][j]) +
                               phi_answer_2 * theta_answer_2 * (map[i][j] - map[i + 1][j] - map[i][j + 1] + map[i + 1][j + 1]);

                    condition_answer = condition_1(fxx, fyy, fxy);

                    phi_answer_2 = phi + phi_answer_2 * map_parameter;

                    theta_answer_2 = theta + theta_answer_2 * map_parameter;

                    fprintf(fp, "%.21Le %.21Le %.21Le %u \n", phi_answer_2, theta_answer_2, p_answer / sigma_p, condition_answer);

                }
            }
        }
    }

    fclose(fp);
}
