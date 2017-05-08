#include <cmath>
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>

#include "functionals.hpp"

#include "parameters.hpp"
#include "constants.hpp"
#include "distance.hpp"
#include "utils.hpp"
#include "fft.hpp"

long double area(long double** map, long double level) {

    long double sum_without_norm = 0.0L;
    long double norm = 0.0L;

    for (unsigned int i = 0; i < npix; ++i) {
        for (unsigned int j = 1; j < npix / 2 - 1; ++j) {

            if (!is_equal(map[i][j], level) and !is_equal(map[i + 1][j], level)
                and !is_equal(map[i][j + 1], level) and !is_equal(map[i + 1][j + 1], level)) {

                if (map[i][j] + map[i + 1][j] + map[i][j + 1] + map[i + 1][j + 1] > level * 4.0L)
                    sum_without_norm = sum_without_norm +
                                       sinl(map_parameter * static_cast<long double>(j) + map_parameter / 2.0L);
                norm = norm + sinl(map_parameter * static_cast<long double>(j) + map_parameter / 2.0L);
            }
        }
    }
    return sum_without_norm / norm;
}

long double area(long double** map, long double level, unsigned int l_1, unsigned int l_2) {

    long double sum_without_norm = 0.0L;
    long double norm = 0.0L;

    if (l_1 < 1 or l_1 > npix / 2 - 2 or l_2 < 1 or l_2 > npix / 2 - 2 or l_1 >= l_2) {
        std::cout << "Wrong range! l_1: " << l_1 << " l_2: " << l_2 << std::endl;
        exit (EXIT_FAILURE);
    }

    for (unsigned int i = 0; i < npix; ++i) {
        for (unsigned int j = l_1; j < l_2 + 1; ++j) {

            if (!is_equal(map[i][j], level) and !is_equal(map[i + 1][j], level)
                and !is_equal(map[i][j + 1], level) and !is_equal(map[i + 1][j + 1], level)) {

                if (map[i][j] + map[i + 1][j] + map[i][j + 1] + map[i + 1][j + 1] > level * 4.0L)
                    sum_without_norm = sum_without_norm +
                            sinl(map_parameter * static_cast<long double>(j) + map_parameter / 2.0L);
                norm = norm + sinl(map_parameter * static_cast<long double>(j) + map_parameter / 2.0L);
            }
        }
    }
    return sum_without_norm / norm;
}

long double area(long double** map, long double level, unsigned int l_1, unsigned int l_2,
                 unsigned int phi_1, unsigned int phi_2) {

    long double sum_without_norm = 0.0L;
    long double norm = 0.0L;

    if (l_1 < 1 or l_1 > npix / 2 - 2 or l_2 < 1 or l_2 > npix / 2 - 2 or l_1 >= l_2) {
        std::cout << "Wrong range! l_1: " << l_1 << " l_2: " << l_2 << std::endl;
        exit (EXIT_FAILURE);
    }

    if (phi_1 < 0 or phi_1 > npix - 1 or phi_2 < 0 or phi_2 > npix - 1 or phi_1 >= phi_2) {
        std::cout << "Wrong range! phi_1: " << phi_1 << " phi_2: " << phi_2 << std::endl;
        exit (EXIT_FAILURE);
    }

    for (unsigned int i = phi_1; i < phi_2; ++i) {
        for (unsigned int j = l_1; j < l_2 + 1; ++j) {

            if (!is_equal(map[i][j], level) and !is_equal(map[i + 1][j], level)
                and !is_equal(map[i][j + 1], level) and !is_equal(map[i + 1][j + 1], level)) {

                if (map[i][j] + map[i + 1][j] + map[i][j + 1] + map[i + 1][j + 1] > level * 4.0L)
                    sum_without_norm = sum_without_norm +
                                       sinl(map_parameter * static_cast<long double>(j) + map_parameter / 2.0L);
                norm = norm + sinl(map_parameter * static_cast<long double>(j) + map_parameter / 2.0L);
            }
        }
    }
    return sum_without_norm / norm;
}

long double length(long double** map, long double level) {

    long double l = 0.0L;

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
                    theta1 = theta + map_parameter * fabsl((map[i][j] - level))
                                     / (fabsl((map[i][j] - level)) + fabsl((map[i][j + 1] - level)));

                    phi2 = phi + map_parameter * fabsl((map[i][j] - level))
                                 / (fabsl((map[i][j] - level)) + fabsl((map[i + 1][j] - level)));
                    theta2 = theta;

                    l += s2(phi1, phi2, theta1, theta2);

                } else if ((map[i + 1][j] - level) * (map[i + 1][j + 1] - level) < 0.0L) {

                    phi1 = phi;
                    theta1 = theta + map_parameter * fabsl((map[i][j] - level))
                                     / (fabsl((map[i][j] - level)) + fabsl((map[i][j + 1] - level)));

                    phi2 = phi + map_parameter;
                    theta2 = theta + map_parameter * fabsl((map[i + 1][j] - level))
                                     / (fabsl((map[i + 1][j] - level)) + fabsl((map[i + 1][j + 1] - level)));

                    l += s2(phi1, phi2, theta1, theta2);

                } else if ((map[i][j + 1] - level) * (map[i + 1][j + 1] - level) < 0.0L) {

                    phi1 = phi;
                    theta1 = theta + map_parameter * fabsl((map[i][j] - level))
                                     / (fabsl((map[i][j] - level)) + fabsl((map[i][j + 1] - level)));

                    phi2 = phi + map_parameter * fabsl((map[i][j + 1] - level))
                                 / (fabsl((map[i][j + 1] - level)) + fabsl((map[i + 1][j + 1] - level)));
                    theta2 = theta + map_parameter;

                    l += s2(phi1, phi2, theta1, theta2);
                }
            } else if ((map[i][j] - level) * (map[i + 1][j] - level) < 0.0L) {

                if ((map[i + 1][j] - level) * (map[i + 1][j + 1] - level) < 0.0L) {

                    phi1 = phi + map_parameter * fabsl((map[i][j] - level))
                                 / (fabsl((map[i][j] - level)) + fabsl((map[i + 1][j] - level)));
                    theta1 = theta;

                    phi2 = phi + map_parameter;
                    theta2 = theta + map_parameter * fabsl((map[i + 1][j] - level))
                                     / (fabsl((map[i + 1][j] - level)) + fabsl((map[i + 1][j + 1] - level)));

                    l += s2(phi1, phi2, theta1, theta2);

                } else if ((map[i][j + 1] - level) * (map[i + 1][j + 1] - level) < 0.0L) {

                    phi1 = phi + map_parameter * fabsl((map[i][j] - level))
                                 / (fabsl((map[i][j] - level)) + fabsl((map[i + 1][j] - level)));
                    theta1 = theta;

                    phi2 = phi + map_parameter * fabsl((map[i][j + 1] - level))
                                 / (fabsl((map[i][j + 1] - level)) + fabsl((map[i + 1][j + 1] - level)));
                    theta2 = theta + map_parameter;

                    l += s2(phi1, phi2, theta1, theta2);
                }
            } else if ((map[i + 1][j] - level) * (map[i + 1][j + 1] - level) < 0.0L) {

                if ((map[i][j + 1] - level) * (map[i + 1][j + 1] - level) < 0.0L) {

                    phi1 = phi + map_parameter;
                    theta1 = theta + map_parameter * fabsl((map[i + 1][j] - level))
                                     / (fabsl((map[i + 1][j] - level)) + fabsl((map[i + 1][j + 1] - level)));

                    phi2 = phi + map_parameter * fabsl((map[i][j + 1] - level))
                                 / (fabsl((map[i][j + 1] - level)) + fabsl((map[i + 1][j + 1] - level)));
                    theta2 = theta + map_parameter;

                    l += s2(phi1, phi2, theta1, theta2);
                }
            }
        }
    }

    return l / PI;
}

long double length(long double** map, long double level, unsigned int l_1, unsigned int l_2) {

    long double l = 0.0L;

    long double phi;
    long double theta;

    long double phi1;
    long double phi2;
    long double theta1;
    long double theta2;

    if (l_1 < 1 or l_1 > npix / 2 - 2 or l_2 < 1 or l_2 > npix / 2 - 2 or l_1 >= l_2) {
        std::cout << "Wrong range! l_1: " << l_1 << " l_2: " << l_2 << std::endl;
        exit (EXIT_FAILURE);
    }

    for (unsigned int i = 0; i < npix; ++i) {
        for (unsigned int j = l_1; j < l_2 + 1; ++j) {

            phi = 2.0L * PI * i / long_npix;
            theta = 2.0L * PI * j / long_npix;

            if ((map[i][j] - level) == 0.0L and ((map[i][j + 1] - level) == 0.0L or (map[i + 1][j] - level) == 0.0L)) {
                l += 0.0L;
            }

            if ((map[i][j] - level) * (map[i][j + 1] - level) < 0.0L) {

                if ((map[i][j] - level) * (map[i + 1][j] - level) < 0.0L) {

                    phi1 = phi;
                    theta1 = theta + map_parameter * fabsl((map[i][j] - level))
                                     / (fabsl((map[i][j] - level)) + fabsl((map[i][j + 1] - level)));

                    phi2 = phi + map_parameter * fabsl((map[i][j] - level))
                                 / (fabsl((map[i][j] - level)) + fabsl((map[i + 1][j] - level)));
                    theta2 = theta;

                    l += s2(phi1, phi2, theta1, theta2);

                } else if ((map[i + 1][j] - level) * (map[i + 1][j + 1] - level) < 0.0L) {

                    phi1 = phi;
                    theta1 = theta + map_parameter * fabsl((map[i][j] - level))
                                     / (fabsl((map[i][j] - level)) + fabsl((map[i][j + 1] - level)));

                    phi2 = phi + map_parameter;
                    theta2 = theta + map_parameter * fabsl((map[i + 1][j] - level))
                                     / (fabsl((map[i + 1][j] - level)) + fabsl((map[i + 1][j + 1] - level)));

                    l += s2(phi1, phi2, theta1, theta2);

                } else if ((map[i][j + 1] - level) * (map[i + 1][j + 1] - level) < 0.0L) {

                    phi1 = phi;
                    theta1 = theta + map_parameter * fabsl((map[i][j] - level))
                                     / (fabsl((map[i][j] - level)) + fabsl((map[i][j + 1] - level)));

                    phi2 = phi + map_parameter * fabsl((map[i][j + 1] - level))
                                 / (fabsl((map[i][j + 1] - level)) + fabsl((map[i + 1][j + 1] - level)));
                    theta2 = theta + map_parameter;

                    l += s2(phi1, phi2, theta1, theta2);
                }
            } else if ((map[i][j] - level) * (map[i + 1][j] - level) < 0.0L) {

                if ((map[i + 1][j] - level) * (map[i + 1][j + 1] - level) < 0.0L) {

                    phi1 = phi + map_parameter * fabsl((map[i][j] - level))
                                 / (fabsl((map[i][j] - level)) + fabsl((map[i + 1][j] - level)));
                    theta1 = theta;

                    phi2 = phi + map_parameter;
                    theta2 = theta + map_parameter * fabsl((map[i + 1][j] - level))
                                     / (fabsl((map[i + 1][j] - level)) + fabsl((map[i + 1][j + 1] - level)));

                    l += s2(phi1, phi2, theta1, theta2);

                } else if ((map[i][j + 1] - level) * (map[i + 1][j + 1] - level) < 0.0L) {

                    phi1 = phi + map_parameter * fabsl((map[i][j] - level))
                                 / (fabsl((map[i][j] - level)) + fabsl((map[i + 1][j] - level)));
                    theta1 = theta;

                    phi2 = phi + map_parameter * fabsl((map[i][j + 1] - level))
                                 / (fabsl((map[i][j + 1] - level)) + fabsl((map[i + 1][j + 1] - level)));
                    theta2 = theta + map_parameter;

                    l += s2(phi1, phi2, theta1, theta2);
                }
            } else if ((map[i + 1][j] - level) * (map[i + 1][j + 1] - level) < 0.0L) {

                if ((map[i][j + 1] - level) * (map[i + 1][j + 1] - level) < 0.0L) {

                    phi1 = phi + map_parameter;
                    theta1 = theta + map_parameter * fabsl((map[i + 1][j] - level))
                                     / (fabsl((map[i + 1][j] - level)) + fabsl((map[i + 1][j + 1] - level)));

                    phi2 = phi + map_parameter * fabsl((map[i][j + 1] - level))
                                 / (fabsl((map[i][j + 1] - level)) + fabsl((map[i + 1][j + 1] - level)));
                    theta2 = theta + map_parameter;

                    l += s2(phi1, phi2, theta1, theta2);
                }
            }
        }
    }

    return l / PI;
}

long double length(long double** map, long double level, unsigned int l_1, unsigned int l_2,
                   unsigned int phi_1, unsigned int phi_2) {

    long double l = 0.0L;

    long double phi;
    long double theta;

    long double phi1;
    long double phi2;
    long double theta1;
    long double theta2;

    if (l_1 < 1 or l_1 > npix / 2 - 2 or l_2 < 1 or l_2 > npix / 2 - 2 or l_1 >= l_2) {
        std::cout << "Wrong range! l_1: " << l_1 << " l_2: " << l_2 << std::endl;
        exit (EXIT_FAILURE);
    }

    for (unsigned int i = phi_1; i < phi_2; ++i) {
        for (unsigned int j = l_1; j < l_2 + 1; ++j) {

            phi = 2.0L * PI * i / long_npix;
            theta = 2.0L * PI * j / long_npix;

            if ((map[i][j] - level) == 0.0L and ((map[i][j + 1] - level) == 0.0L or (map[i + 1][j] - level) == 0.0L)) {
                l += 0.0L;
            }

            if ((map[i][j] - level) * (map[i][j + 1] - level) < 0.0L) {

                if ((map[i][j] - level) * (map[i + 1][j] - level) < 0.0L) {

                    phi1 = phi;
                    theta1 = theta + map_parameter * fabsl((map[i][j] - level))
                                     / (fabsl((map[i][j] - level)) + fabsl((map[i][j + 1] - level)));

                    phi2 = phi + map_parameter * fabsl((map[i][j] - level))
                                 / (fabsl((map[i][j] - level)) + fabsl((map[i + 1][j] - level)));
                    theta2 = theta;

                    l += s2(phi1, phi2, theta1, theta2);

                } else if ((map[i + 1][j] - level) * (map[i + 1][j + 1] - level) < 0.0L) {

                    phi1 = phi;
                    theta1 = theta + map_parameter * fabsl((map[i][j] - level))
                                     / (fabsl((map[i][j] - level)) + fabsl((map[i][j + 1] - level)));

                    phi2 = phi + map_parameter;
                    theta2 = theta + map_parameter * fabsl((map[i + 1][j] - level))
                                     / (fabsl((map[i + 1][j] - level)) + fabsl((map[i + 1][j + 1] - level)));

                    l += s2(phi1, phi2, theta1, theta2);

                } else if ((map[i][j + 1] - level) * (map[i + 1][j + 1] - level) < 0.0L) {

                    phi1 = phi;
                    theta1 = theta + map_parameter * fabsl((map[i][j] - level))
                                     / (fabsl((map[i][j] - level)) + fabsl((map[i][j + 1] - level)));

                    phi2 = phi + map_parameter * fabsl((map[i][j + 1] - level))
                                 / (fabsl((map[i][j + 1] - level)) + fabsl((map[i + 1][j + 1] - level)));
                    theta2 = theta + map_parameter;

                    l += s2(phi1, phi2, theta1, theta2);
                }
            } else if ((map[i][j] - level) * (map[i + 1][j] - level) < 0.0L) {

                if ((map[i + 1][j] - level) * (map[i + 1][j + 1] - level) < 0.0L) {

                    phi1 = phi + map_parameter * fabsl((map[i][j] - level))
                                 / (fabsl((map[i][j] - level)) + fabsl((map[i + 1][j] - level)));
                    theta1 = theta;

                    phi2 = phi + map_parameter;
                    theta2 = theta + map_parameter * fabsl((map[i + 1][j] - level))
                                     / (fabsl((map[i + 1][j] - level)) + fabsl((map[i + 1][j + 1] - level)));

                    l += s2(phi1, phi2, theta1, theta2);

                } else if ((map[i][j + 1] - level) * (map[i + 1][j + 1] - level) < 0.0L) {

                    phi1 = phi + map_parameter * fabsl((map[i][j] - level))
                                 / (fabsl((map[i][j] - level)) + fabsl((map[i + 1][j] - level)));
                    theta1 = theta;

                    phi2 = phi + map_parameter * fabsl((map[i][j + 1] - level))
                                 / (fabsl((map[i][j + 1] - level)) + fabsl((map[i + 1][j + 1] - level)));
                    theta2 = theta + map_parameter;

                    l += s2(phi1, phi2, theta1, theta2);
                }
            } else if ((map[i + 1][j] - level) * (map[i + 1][j + 1] - level) < 0.0L) {

                if ((map[i][j + 1] - level) * (map[i + 1][j + 1] - level) < 0.0L) {

                    phi1 = phi + map_parameter;
                    theta1 = theta + map_parameter * fabsl((map[i + 1][j] - level))
                                     / (fabsl((map[i + 1][j] - level)) + fabsl((map[i + 1][j + 1] - level)));

                    phi2 = phi + map_parameter * fabsl((map[i][j + 1] - level))
                                 / (fabsl((map[i][j + 1] - level)) + fabsl((map[i + 1][j + 1] - level)));
                    theta2 = theta + map_parameter;

                    l += s2(phi1, phi2, theta1, theta2);
                }
            }
        }
    }

    return l / PI;
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
        mean_root1 = mean_t_solver(qx, qy, ux, uy, root1);
        mean_root2 = mean_t_solver(qx, qy, ux, uy, root2);

        if (is_equal(root1, 0.0l) and is_equal(root2, 0.0l)) {
            return 3;
        } else {
            if ((mean_root1 > 0.0L and mean_root2 < 0.0L) or (mean_root1 < 0.0L and  mean_root2 > 0.0L)) {
                return 2;
            } else {
                return 3;
            }
        }
    }
}

void points_classifier(long double** map_x, long double** map_y, long double** cos_ml, long double** sin_ml,
                       std::string name) {

    std::ofstream out_file;
    out_file.open(name);
    typedef std::numeric_limits<long double> dbl;
    out_file.precision(dbl::max_digits10);

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

                    fxx = fft_point_xx_forward(phi_answer, theta_answer, cos_ml, sin_ml);
                    fyy = fft_point_yy_forward(phi_answer, theta_answer, cos_ml, sin_ml);
                    fxy = fft_point_xy_forward(phi_answer, theta_answer, cos_ml, sin_ml);

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

void points_classifier(long double** map_x, long double** map_y, long double** cos_ml, long double** sin_ml,
                       std::string name, int** whitelist) {

    std::ofstream out_file;
    out_file.open(name);
    typedef std::numeric_limits<long double> dbl;
    out_file.precision(dbl::max_digits10);

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

            // whitelist[i][j] != 1 через is_equal

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

                    fxx = fft_point_xx_forward(phi_answer, theta_answer, cos_ml, sin_ml);
                    fyy = fft_point_yy_forward(phi_answer, theta_answer, cos_ml, sin_ml);
                    fxy = fft_point_xy_forward(phi_answer, theta_answer, cos_ml, sin_ml);

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
//void points_classifier_p(long double** map_x, long double** map_y, long double** q_cos_ml, long double** q_sin_ml,
//                         long double** u_cos_ml, long double** u_sin_ml, std::string name, long double** whitelist,
//                         unsigned int n, long double sigma_p)
void points_classifier_p(long double** map, long double** map_x, long double** map_y, long double** map_xx,
                         long double** map_yy, long double** map_xy, std::string name, long double** whitelist,
                         unsigned int n, long double sigma_p) {

    std::ofstream out_file;
    out_file.open(name);
    typedef std::numeric_limits<long double> dbl;
    out_file.precision(dbl::max_digits10);

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

//    long double q_answer;
//    long double u_answer;
//    long double q_x_answer;
//    long double q_y_answer;
//    long double u_x_answer;
//    long double u_y_answer;
//    long double q_xx_answer;
//    long double q_yy_answer;
//    long double q_xy_answer;
//    long double u_xx_answer;
//    long double u_yy_answer;
//    long double u_xy_answer;
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

//                    q_answer = fft_point_forward(phi_answer, theta_answer, q_cos_ml,  q_sin_ml, n);
//                    u_answer = fft_point_forward(phi_answer, theta_answer, u_cos_ml,  u_sin_ml, n);
//                    q_x_answer = fft_point_x_forward(phi_answer, theta_answer, q_cos_ml,  q_sin_ml, n);
//                    q_y_answer = fft_point_y_forward(phi_answer, theta_answer, q_cos_ml,  q_sin_ml, n);
//                    u_x_answer = fft_point_x_forward(phi_answer, theta_answer, u_cos_ml,  u_sin_ml, n);
//                    u_y_answer = fft_point_y_forward(phi_answer, theta_answer, u_cos_ml,  u_sin_ml, n);
//                    q_xx_answer = fft_point_xx_forward(phi_answer, theta_answer, q_cos_ml,  q_sin_ml, n);
//                    q_yy_answer = fft_point_yy_forward(phi_answer, theta_answer, q_cos_ml,  q_sin_ml, n);
//                    q_xy_answer = fft_point_xy_forward(phi_answer, theta_answer, q_cos_ml,  q_sin_ml, n);
//                    u_xx_answer = fft_point_xx_forward(phi_answer, theta_answer, u_cos_ml,  u_sin_ml, n);
//                    u_yy_answer = fft_point_yy_forward(phi_answer, theta_answer, u_cos_ml,  u_sin_ml, n);
//                    u_xy_answer = fft_point_xy_forward(phi_answer, theta_answer, u_cos_ml,  u_sin_ml, n);
//                    p_answer = sqrtl(q_answer * q_answer + u_answer * u_answer);

//                    fxx = (q_x_answer * q_x_answer + u_x_answer * u_x_answer + q_answer * q_xx_answer + u_answer * u_xx_answer);
//                    fyy = (q_y_answer * q_y_answer + u_y_answer * u_y_answer + q_answer * q_yy_answer + u_answer * u_yy_answer);
//                    fxy = (q_y_answer * q_x_answer + u_y_answer * u_x_answer + q_answer * q_xy_answer + u_answer * u_xy_answer);

                    fxx = map_xx[i][j] ;
//                          + ((map_xx[i + 1][j] - map_xx[i][j]) * ((theta_answer - theta) / map_parameter)
//                                            + (map_xx[i + 1][j + 1] - map_xx[i][j + 1]) * ((theta + map_parameter - theta_answer) / map_parameter)) *
//                                           s2(phi, phi_answer, theta_answer, theta_answer) / s2(phi, phi + map_parameter, theta_answer, theta_answer)
//                                + ((map_xx[i][j + 1] - map_xx[i][j]) * ((phi_answer - phi) / map_parameter)
//                                   + (map_xx[i + 1][j + 1] - map_xx[i + 1][j]) * (phi_answer + map_parameter - phi) / map_parameter)
//                                  * s2(phi_answer, phi_answer, theta, theta_answer) / s2(phi_answer, phi_answer, theta, theta + map_parameter);

                    fyy = map_yy[i][j] ;
//                          + ((map_yy[i + 1][j] - map_yy[i][j]) * ((theta_answer - theta) / map_parameter)
//                                            + (map_yy[i + 1][j + 1] - map_yy[i][j + 1]) * ((theta + map_parameter - theta_answer) / map_parameter)) *
//                                           s2(phi, phi_answer, theta_answer, theta_answer) / s2(phi, phi + map_parameter, theta_answer, theta_answer)
//                                + ((map_yy[i][j + 1] - map_yy[i][j]) * ((phi_answer - phi) / map_parameter)
//                                   + (map_yy[i + 1][j + 1] - map_yy[i + 1][j]) * (phi_answer + map_parameter - phi) / map_parameter)
//                                  * s2(phi_answer, phi_answer, theta, theta_answer) / s2(phi_answer, phi_answer, theta, theta + map_parameter);

                    fxy = map_xy[i][j];
//                          + ((map_xy[i + 1][j] - map_xy[i][j]) * ((theta_answer - theta) / map_parameter)
//                                            + (map_xy[i + 1][j + 1] - map_xy[i][j + 1]) * ((theta + map_parameter - theta_answer) / map_parameter)) *
//                                           s2(phi, phi_answer, theta_answer, theta_answer) / s2(phi, phi + map_parameter, theta_answer, theta_answer)
//                                + ((map_xy[i][j + 1] - map_xy[i][j]) * ((phi_answer - phi) / map_parameter)
//                                   + (map_xy[i + 1][j + 1] - map_xy[i + 1][j]) * (phi_answer + map_parameter - phi) / map_parameter)
//                                  * s2(phi_answer, phi_answer, theta, theta_answer) / s2(phi_answer, phi_answer, theta, theta + map_parameter);

                    p_answer = map[i][j];
//                               + ((map[i + 1][j] - map[i][j]) * ((theta_answer - theta) / map_parameter)
//                                      + (map[i + 1][j + 1] - map[i][j + 1]) * ((theta + map_parameter - theta_answer) / map_parameter)) *
//                                     s2(phi, phi_answer, theta_answer, theta_answer) / s2(phi, phi + map_parameter, theta_answer, theta_answer)
//                          + ((map[i][j + 1] - map[i][j]) * ((phi_answer - phi) / map_parameter)
//                             + (map[i + 1][j + 1] - map[i + 1][j]) * (phi_answer + map_parameter - phi) / map_parameter)
//                            * s2(phi_answer, phi_answer, theta, theta_answer) / s2(phi_answer, phi_answer, theta, theta + map_parameter);

                    condition_answer = condition_1(fxx, fyy, fxy);

                    out_file << std::scientific << phi_answer << " "
                             << std::scientific << theta_answer << " "
                             << std::scientific << p_answer / sigma_p << " "
                             << condition_answer << std::endl;
                }
            }
        }
    }

    out_file.close();
}

void level_points_classifier(long double** map_x, long double** map_y, long double** cos_ml, long double** sin_ml,
                       long double level, std::string name) {

    std::ofstream out_file;
    out_file.open(name, std::ios::app);
    typedef std::numeric_limits<long double> dbl;
    out_file.precision(dbl::max_digits10);

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

    long double f;
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

                    f = fft_point_forward(phi_answer, theta_answer, cos_ml, sin_ml);

                    fxx = fft_point_xx_forward(phi_answer, theta_answer, cos_ml, sin_ml);
                    fyy = fft_point_yy_forward(phi_answer, theta_answer, cos_ml, sin_ml);
                    fxy = fft_point_xy_forward(phi_answer, theta_answer, cos_ml, sin_ml);

                    condition_answer = condition_1(fxx, fyy, fxy);

                    if (f > level) {
                        out_file << std::scientific << phi_answer << " "
                                 << std::scientific << theta_answer << " "
                                 << condition_answer << " "
                                 << std::scientific << level << std::endl;
                    }
                }
            }
        }
    }

    out_file.close();
}

void points_classifier(long double** map_x, long double** map_y, long double** cos_ml, long double** sin_ml,
                       std::string name, unsigned int l_1, unsigned int l_2) {

    if (l_1 < 1 or l_1 > npix / 2 - 2 or l_2 < 1 or l_2 > npix / 2 - 2 or l_1 >= l_2) {
        std::cout << "Wrong range! l_1: " << l_1 << " l_2: " << l_2 << std::endl;
        exit (EXIT_FAILURE);
    }

    std::ofstream out_file;
    out_file.open(name, std::ios::app);
    typedef std::numeric_limits<long double> dbl;
    out_file.precision(dbl::max_digits10);

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
        for (unsigned int j = l_1; j < l_2 + 1; ++j) {

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

                    fxx = fft_point_xx_forward(phi_answer, theta_answer, cos_ml, sin_ml);
                    fyy = fft_point_yy_forward(phi_answer, theta_answer, cos_ml, sin_ml);
                    fxy = fft_point_xy_forward(phi_answer, theta_answer, cos_ml, sin_ml);

                    condition_answer = condition_1(fxx, fyy, fxy);

                    out_file << std::scientific << phi_answer << " "
                             << std::scientific << theta_answer << " "
                             << condition_answer << " "
                             << l_1 << " " << l_2
                             << std::endl;
                }
            }
        }
    }

    out_file.close();
}

void singular_points_classifier(long double** q, long double** u, long double** q_cos_ml, long double** q_sin_ml,
                                long double** u_cos_ml, long double** u_sin_ml, std::string name) {
    std::ofstream out_file;
    out_file.open(name);
    typedef std::numeric_limits<long double> dbl;
    out_file.precision(dbl::max_digits10);

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

                    qx = fft_point_x_forward(phi_answer, theta_answer, q_cos_ml, q_sin_ml, nback);
                    qy = fft_point_y_forward(phi_answer, theta_answer, q_cos_ml, q_sin_ml, nback);
                    ux = fft_point_x_forward(phi_answer, theta_answer, u_cos_ml, u_sin_ml, nback);
                    uy = fft_point_y_forward(phi_answer, theta_answer, u_cos_ml, u_sin_ml, nback);

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

//void singular_points_classifier(long double** q, long double** u, long double** q_cos_ml, long double** q_sin_ml,
//                                long double** u_cos_ml, long double** u_sin_ml, std::string name, long double** whitelist,
//                                unsigned int n) {
//
//}
void singular_points_classifier(long double** q, long double** u, long double** qx, long double** ux,
                                long double** qy, long double** uy, std::string name, long double** whitelist,
                                unsigned int n) {
//    std::ofstream out_file;
//    out_file.open(name);
//    typedef std::numeric_limits<long double> dbl;
//    out_file.precision(dbl::max_digits10);

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

//                    qx = fft_point_x_forward(phi_answer, theta_answer, q_cos_ml, q_sin_ml, n);
//                    qy = fft_point_y_forward(phi_answer, theta_answer, q_cos_ml, q_sin_ml, n);
//                    ux = fft_point_x_forward(phi_answer, theta_answer, u_cos_ml, u_sin_ml, n);
//                    uy = fft_point_y_forward(phi_answer, theta_answer, u_cos_ml, u_sin_ml, n);

                    qx_answer = qx[i][j] + ((qx[i + 1][j] - qx[i][j]) * ((theta_answer - theta) / map_parameter)
                         + (qx[i + 1][j + 1] - qx[i][j + 1]) * ((theta + map_parameter - theta_answer) / map_parameter)) *
                         s2(phi, phi_answer, theta_answer, theta_answer) / s2(phi, phi + map_parameter, theta_answer, theta_answer)
                         + ((qx[i][j + 1] - qx[i][j]) * ((phi_answer - phi) / map_parameter)
                         + (qx[i + 1][j + 1] - qx[i + 1][j]) * (phi_answer + map_parameter - phi) / map_parameter)
                         * s2(phi_answer, phi_answer, theta, theta_answer) / s2(phi_answer, phi_answer, theta, theta + map_parameter);

                    ux_answer = ux[i][j] + ((ux[i + 1][j] - ux[i][j]) * ((theta_answer - theta) / map_parameter)
                                            + (ux[i + 1][j + 1] - ux[i][j + 1]) * ((theta + map_parameter - theta_answer) / map_parameter)) *
                                           s2(phi, phi_answer, theta_answer, theta_answer) / s2(phi, phi + map_parameter, theta_answer, theta_answer)
                                + ((ux[i][j + 1] - ux[i][j]) * ((phi_answer - phi) / map_parameter)
                                   + (ux[i + 1][j + 1] - ux[i + 1][j]) * (phi_answer + map_parameter - phi) / map_parameter)
                                  * s2(phi_answer, phi_answer, theta, theta_answer) / s2(phi_answer, phi_answer, theta, theta + map_parameter);

                    qy_answer = qy[i][j] + ((qy[i + 1][j] - qy[i][j]) * ((theta_answer - theta) / map_parameter)
                                            + (qy[i + 1][j + 1] - qy[i][j + 1]) * ((theta + map_parameter - theta_answer) / map_parameter)) *
                                           s2(phi, phi_answer, theta_answer, theta_answer) / s2(phi, phi + map_parameter, theta_answer, theta_answer)
                                + ((qy[i][j + 1] - qy[i][j]) * ((phi_answer - phi) / map_parameter)
                                   + (qy[i + 1][j + 1] - qy[i + 1][j]) * (phi_answer + map_parameter - phi) / map_parameter)
                                  * s2(phi_answer, phi_answer, theta, theta_answer) / s2(phi_answer, phi_answer, theta, theta + map_parameter);

                    uy_answer = uy[i][j] + ((uy[i + 1][j] - uy[i][j]) * ((theta_answer - theta) / map_parameter)
                                            + (uy[i + 1][j + 1] - uy[i][j + 1]) * ((theta + map_parameter - theta_answer) / map_parameter)) *
                                           s2(phi, phi_answer, theta_answer, theta_answer) / s2(phi, phi + map_parameter, theta_answer, theta_answer)
                                + ((uy[i][j + 1] - uy[i][j]) * ((phi_answer - phi) / map_parameter)
                                   + (uy[i + 1][j + 1] - uy[i + 1][j]) * (phi_answer + map_parameter - phi) / map_parameter)
                                  * s2(phi_answer, phi_answer, theta, theta_answer) / s2(phi_answer, phi_answer, theta, theta + map_parameter);


                    condition_answer = condition_2(qx_answer, qy_answer, ux_answer, uy_answer);

                    fprintf(fp, "%.21Le %.21Le %u \n", phi_answer, theta_answer, condition_answer);

//                    out_file << std::scientific << phi_answer << " "
//                             << std::scientific << theta_answer << " "
//                             << condition_answer << std::endl;

                    whitelist[i][j] = 1;
                }
            }
        }
    }

    fclose(fp);
//    out_file.close();
}

void level_singular_points_classifier(long double** q, long double** u, long double** q_cos_ml, long double** q_sin_ml,
                                long double** u_cos_ml, long double** u_sin_ml, long double level, std::string name) {
    std::ofstream out_file;
    out_file.open(name, std::ios::app);
    typedef std::numeric_limits<long double> dbl;
    out_file.precision(dbl::max_digits10);

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

    long double q_pr;
    long double u_pr;
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

                    q_pr = fft_point_forward(phi_answer, theta_answer, q_cos_ml, q_sin_ml);
                    u_pr = fft_point_forward(phi_answer, theta_answer, u_cos_ml, u_sin_ml);

                    qx = fft_point_x_forward(phi_answer, theta_answer, q_cos_ml, q_sin_ml, nback);
                    qy = fft_point_y_forward(phi_answer, theta_answer, q_cos_ml, q_sin_ml, nback);
                    ux = fft_point_x_forward(phi_answer, theta_answer, u_cos_ml, u_sin_ml, nback);
                    uy = fft_point_y_forward(phi_answer, theta_answer, u_cos_ml, u_sin_ml, nback);

                    condition_answer = condition_2(qx, qy, ux, uy);

                    if (sqrtl(q_pr * q_pr + u_pr * u_pr) > level) {
                        out_file << std::scientific << phi_answer << " "
                                 << std::scientific << theta_answer << " "
                                 << condition_answer << " "
                                 << std::scientific << level << std::endl;
                    }
                }
            }
        }
    }

    out_file.close();
}

void singular_points_classifier(long double** q, long double** u, long double** q_cos_ml, long double** q_sin_ml,
                                long double** u_cos_ml, long double** u_sin_ml, std::string name,
                                unsigned int l_1, unsigned int l_2) {

    if (l_1 < 1 or l_1 > npix / 2 - 2 or l_2 < 1 or l_2 > npix / 2 - 2 or l_1 >= l_2) {
        std::cout << "Wrong range! l_1: " << l_1 << " l_2: " << l_2 << std::endl;
        exit (EXIT_FAILURE);
    }

    std::ofstream out_file;
    out_file.open(name, std::ios::app);
    typedef std::numeric_limits<long double> dbl;
    out_file.precision(dbl::max_digits10);

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
        for (unsigned int j = l_1; j < l_2 + 1; ++j) {

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

                    qx = fft_point_x_forward(phi_answer, theta_answer, q_cos_ml, q_sin_ml, nback);
                    qy = fft_point_y_forward(phi_answer, theta_answer, q_cos_ml, q_sin_ml, nback);
                    ux = fft_point_x_forward(phi_answer, theta_answer, u_cos_ml, u_sin_ml, nback);
                    uy = fft_point_y_forward(phi_answer, theta_answer, u_cos_ml, u_sin_ml, nback);

                    condition_answer = condition_2(qx, qy, ux, uy);

                    out_file << std::scientific << phi_answer << " "
                             << std::scientific << theta_answer << " "
                             << condition_answer << " "
                             << l_1 << " " << l_2
                             << std::endl;
                }
            }
        }
    }

    out_file.close();
}
