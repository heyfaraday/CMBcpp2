#include <cmath>

#include "pml.hpp"

#include "utils.hpp"
#include "constants.hpp"
#include "parameters.hpp"

long double coef1(unsigned int l, unsigned int m);

void pml_gen(unsigned int j, long double** array) {

    long double long_j = static_cast<long double>(j);
    long double theta = map_parameter * long_j;
    long double long_m;
    long double long_l;

    for (unsigned int m = 0; m < nmod; ++m) {
        for (unsigned int l = 0; l < nmod; ++l) {
            array[m][l] = 0.0L;
        }
    }

    array[0][0] = 1.0L / sqrtl(4.0L * PI);

    for (unsigned int m = 0; m < nmod - 1; ++m) {

        long_m = static_cast<long double>(m);

        array[m + 1][m + 1] = - array[m][m] * sinl(theta) * sqrtl(2.0L * long_m + 3.0L)
                              / sqrtl(2.0L * long_m + 2.0L);
    }

    for (unsigned int m = 0; m < nmod - 1; ++m) {

        long_m = static_cast<long double>(m);

        array[m][m + 1] = array[m][m] * cosl(theta) * sqrtl(2.0L * long_m + 3.0L);
    }

    for (unsigned int m = 0; m < nmod - 2; ++m) {
        for (unsigned int l = m + 2; l < nmod; ++l) {

            long_m = static_cast<long double>(m);
            long_l = static_cast<long double>(l);

            array[m][l] = ((2.0L * long_l - 1.0L) * sqrtl((long_l - long_m) * (2.0L * long_l + 1.0L))
                           / sqrtl((long_l + long_m) * (2.0L * l - 1.0L)) * array[m][l - 1] * cosl(theta)
                           - (long_l + long_m - 1.0L) * sqrtl((long_l - long_m) * (long_l - 1.0L - long_m)
                           * (2.0L * long_l + 1.0L)) / sqrtl((long_l + long_m) * (long_l - 1.0L + long_m)
                           * (2.0L * long_l - 3.0L)) * array[m][l - 2]) / (long_l - long_m);
        }
    }

    for (unsigned int m = 1; m < nmod ; ++m) {
        for (unsigned int l = m; l < nmod; ++l) {
            array[m][l] = array[m][l];// * sqrtl(2.0L);
        }
    }
}

void pml_gen(long double theta, long double** array) {

    long double long_m;
    long double long_l;

    for (unsigned int m = 0; m < nmod; ++m) {
        for (unsigned int l = 0; l < nmod; ++l) {
            array[m][l] = 0.0L;
        }
    }

    array[0][0] = 1.0L / sqrtl(4.0L * PI);

    for (unsigned int m = 0; m < nmod - 1; ++m) {

        long_m = static_cast<long double>(m);

        array[m + 1][m + 1] = - array[m][m] * sinl(theta) * sqrtl(2.0L * long_m + 3.0L)
                              / sqrtl(2.0L * long_m + 2.0L);
    }

    for (unsigned int m = 0; m < nmod - 1; ++m) {

        long_m = static_cast<long double>(m);

        array[m][m + 1] = array[m][m] * cosl(theta) * sqrtl(2.0L * long_m + 3.0L);
    }

    for (unsigned int m = 0; m < nmod - 2; ++m) {
        for (unsigned int l = m + 2; l < nmod; ++l) {

            long_m = static_cast<long double>(m);
            long_l = static_cast<long double>(l);

            array[m][l] = ((2.0L * long_l - 1.0L) * sqrtl((long_l - long_m) * (2.0L * long_l + 1.0L))
                           / sqrtl((long_l + long_m) * (2.0L * l - 1.0L)) * array[m][l - 1] * cosl(theta)
                           - (long_l + long_m - 1.0L) * sqrtl((long_l - long_m) * (long_l - 1.0L - long_m)
                           * (2.0L * long_l + 1.0L)) / sqrtl((long_l + long_m) * (long_l - 1.0L + long_m)
                           * (2.0L * long_l - 3.0L)) * array[m][l - 2]) / (long_l - long_m);
        }
    }

    for (unsigned int m = 1; m < nmod ; ++m) {
        for (unsigned int l = m; l < nmod; ++l) {
            array[m][l] = array[m][l] * sqrtl(2.0L);
        }
    }
}

long double coef1(unsigned int l, unsigned int m) {

    long double long_l = static_cast<long double>(l);
    long double long_m = static_cast<long double>(m);

    if (l != 0) {
        return sqrtl((long_l - long_m) * (2.0L * long_l + 1.0L)
                     / ((long_l + long_m) * (2.0L * long_l - 1.0L)));
    } else {
        return 0.0L;
    }
}

void pml_x_gen(unsigned int j, long double** array) {

    long double long_j = static_cast<long double>(j);
    long double theta = 2.0L * long_j * PI / long_npix;

    long double** pml = n_matrix_generator(nmod, nmod);
    pml_gen(j, pml);

    for (unsigned int l = 2; l < nmod ; ++l) {
        for (unsigned int m = 0; m <= l; ++m) {
            array[m][l] = pml[m][l] * static_cast<long double>(m) / sinl(theta);
        }
    }
    n_matrix_destroyer(pml, nmod);
}

void pml_x_gen(long double theta, long double** array) {

    long double** pml = n_matrix_generator(nmod, nmod);
    pml_gen(theta, pml);

    for (unsigned int l = 2; l < nmod ; ++l) {
        for (unsigned int m = 0; m <= l; ++m) {
            array[m][l] = pml[m][l] * static_cast<long double>(m) / sinl(theta);
        }
    }
    n_matrix_destroyer(pml, nmod);
}

void pml_y_gen(unsigned int j, long double** array) {

    long double long_j = static_cast<long double>(j);
    long double theta = 2.0L * long_j * PI / long_npix;
    long double long_m;
    long double long_l;

    long double** pml = n_matrix_generator(nmod, nmod);
    pml_gen(j, pml);

    for (unsigned int l = 2; l < nmod ; ++l) {
        long_l = static_cast<long double>(l);
        for (unsigned int m = 0; m <= l; ++m) {
            long_m = static_cast<long double>(m);
            array[m][l] = long_l * cosl(theta) / sinl(theta) * pml[m][l]
                          - 1.0L / sinl(theta) * (long_l + long_m)
                          * coef1(l, m) * pml[m][l - 1];
        }
    }
    n_matrix_destroyer(pml, nmod);
}

void pml_y_gen(long double theta, long double** array) {

    long double long_m;
    long double long_l;

    long double** pml = n_matrix_generator(nmod, nmod);
    pml_gen(theta, pml);

    for (unsigned int l = 2; l < nmod ; ++l) {
        long_l = static_cast<long double>(l);
        for (unsigned int m = 0; m <= l; ++m) {
            long_m = static_cast<long double>(m);
            array[m][l] = long_l * cosl(theta) / sinl(theta) * pml[m][l]
                          - 1.0L / sinl(theta) * (long_l + long_m)
                          * coef1(l, m) * pml[m][l - 1];
        }
    }
    n_matrix_destroyer(pml, nmod);
}

void pml_xx_gen(unsigned int j, long double** array) {

    long double long_j = static_cast<long double>(j);
    long double theta = 2.0L * long_j * PI / long_npix;
    long double long_m;
    long double long_l;

    long double** pml = n_matrix_generator(nmod, nmod);
    pml_gen(j, pml);

    for (unsigned int l = 2; l < nmod ; ++l) {
        long_l = static_cast<long double>(l);
        for (unsigned int m = 0; m <= l; ++m) {
            long_m = static_cast<long double>(m);
            array[m][l] = - long_m * long_m * pml[m][l] / (sinl(theta) * sinl(theta))
                          + (long_l * cosl(theta) / sinl(theta) * pml[m][l] - 1.0L / sinl(theta) * (long_l + long_m)
                          * coef1(l, m) * pml[m][l - 1]) * cosl(theta) / sinl(theta);
        }
    }
    n_matrix_destroyer(pml, nmod);
}

void pml_xx_gen(long double theta, long double** array) {

    long double long_m;
    long double long_l;

    long double** pml = n_matrix_generator(nmod, nmod);
    pml_gen(theta, pml);

    for (unsigned int l = 2; l < nmod ; ++l) {
        long_l = static_cast<long double>(l);
        for (unsigned int m = 0; m <= l; ++m) {
            long_m = static_cast<long double>(m);
            array[m][l] = - long_m * long_m * pml[m][l] / (sinl(theta) * sinl(theta))
                          + (long_l * cosl(theta) / sinl(theta) * pml[m][l] - 1.0L / sinl(theta) * (long_l + long_m)
                          * coef1(l, m) * pml[m][l - 1]) * cosl(theta) / sinl(theta);
        }
    }
    n_matrix_destroyer(pml, nmod);
}

void pml_yy_gen(unsigned int j, long double** array) {

    long double long_j = static_cast<long double>(j);
    long double theta = 2.0L * long_j * PI / long_npix;
    long double long_m;
    long double long_l;

    long double** pml = n_matrix_generator(nmod, nmod);
    pml_gen(j, pml);

    for (unsigned int l = 2; l < nmod ; ++l) {
        long_l = static_cast<long double>(l);
        for (unsigned int m = 0; m <= l; ++m) {
            long_m = static_cast<long double>(m);
            array[m][l] = 0.5L / sinl(theta) * ((1.0L / sinl(theta)) * (long_l * long_l * cosl(2.0L * theta)
                          - (long_l + 2.0L) * long_l + 2.0L * long_m * long_m) * pml[m][l]
                          + 2.0L * (long_l + long_m) * cosl(theta) / sinl(theta) * coef1(l, m) * pml[m][l - 1]);
        }
    }
    n_matrix_destroyer(pml, nmod);
}

void pml_yy_gen(long double theta, long double** array) {

    long double long_m;
    long double long_l;

    long double** pml = n_matrix_generator(nmod, nmod);
    pml_gen(theta, pml);

    for (unsigned int l = 2; l < nmod ; ++l) {
        long_l = static_cast<long double>(l);
        for (unsigned int m = 0; m <= l; ++m) {
            long_m = static_cast<long double>(m);
            array[m][l] = 0.5L / sinl(theta) * ((1.0L / sinl(theta)) * (long_l * long_l * cosl(2.0L * theta)
                          - (long_l + 2.0L) * long_l + 2.0L * long_m * long_m) * pml[m][l]
                          + 2.0L * (long_l + long_m) * cosl(theta) / sinl(theta) * coef1(l, m) * pml[m][l - 1]);
        }
    }
    n_matrix_destroyer(pml, nmod);
}

void pml_xy_gen(unsigned int j, long double** array) {

    long double long_j = static_cast<long double>(j);
    long double theta = 2.0L * long_j * PI / long_npix;
    long double long_m;
    long double long_l;

    long double** pml = n_matrix_generator(nmod, nmod);
    pml_gen(j, pml);

    for (unsigned int l = 2; l < nmod ; ++l) {
        long_l = static_cast<long double>(l);
        for (unsigned int m = 0; m <= l; ++m) {
            long_m = static_cast<long double>(m);
            array[m][l] = long_m / sinl(theta) * ((1.0L / sinl(theta)) * (l + m) * pml[m][l - 1] * coef1(l, m) -
                                            (long_l - 1.0L) * cosl(theta) / sinl(theta) * pml[m][l]);
        }
    }
    n_matrix_destroyer(pml, nmod);
}

void pml_xy_gen(long double theta, long double** array) {

    long double long_m;
    long double long_l;

    long double** pml = n_matrix_generator(nmod, nmod);
    pml_gen(theta, pml);

    for (unsigned int l = 2; l < nmod ; ++l) {
        long_l = static_cast<long double>(l);
        for (unsigned int m = 0; m <= l; ++m) {
            long_m = static_cast<long double>(m);
            array[m][l] = long_m / sinl(theta) * ((1.0L / sinl(theta)) * (l + m) * pml[m][l - 1] * coef1(l, m) -
                                                  (long_l - 1.0L) * cosl(theta) / sinl(theta) * pml[m][l]);
        }
    }
    n_matrix_destroyer(pml, nmod);
}

void x_1_ml_gen(unsigned int j, long double** array) {

    long double long_j = static_cast<long double>(j);
    long double theta = 2.0L * long_j * PI / long_npix;

    long double** pml_xx = n_matrix_generator(nmod, nmod);
    long double** pml_yy = n_matrix_generator(nmod, nmod);
    pml_xx_gen(j, pml_xx);
    pml_yy_gen(j, pml_yy);

    for (unsigned int l = 2; l < nmod ; ++l) {
        for (unsigned int m = 0; m <= l; ++m) {
            array[m][l] = pml_xx[m][l] - pml_yy[m][l];
        }
    }
    n_matrix_destroyer(pml_xx, nmod);
    n_matrix_destroyer(pml_yy, nmod);
}

void x_2_ml_gen(unsigned int j, long double** array) {

    long double long_j = static_cast<long double>(j);
    long double theta = 2.0L * long_j * PI / long_npix;
    long double long_m;
    long double long_l;

    long double** pml = n_matrix_generator(nmod, nmod);
    pml_gen(j, pml);

    for (unsigned int l = 2; l < nmod ; ++l) {
        long_l = static_cast<long double>(l);
        for (unsigned int m = 0; m <= l; ++m) {
            long_m = static_cast<long double>(m);
            array[m][l] = 2.0L * long_m / sinl(theta) * ((1.0L / sinl(theta)) * (l + m) * pml[m][l - 1] * coef1(l, m) -
                                                  (long_l - 1.0L) * cosl(theta) / sinl(theta) * pml[m][l]);
        }
    }
    n_matrix_destroyer(pml, nmod);
}
