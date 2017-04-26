#include <fftw3.h>
#include <cmath>

#include "fft.hpp"

#include "parameters.hpp"
#include "constants.hpp"
#include "pml.hpp"
#include "utils.hpp"

void fft_map_forward(long double** map, long double** cos_ml, long double** sin_ml) {

    fftwl_complex* in = static_cast<fftwl_complex*>(fftwl_malloc(sizeof(*in) * npix));
    fftwl_complex* out = static_cast<fftwl_complex*>(fftwl_malloc(sizeof(*out) * npix));
    fftwl_plan p = fftwl_plan_dft_1d(int(npix), in, out, FFTW_FORWARD, FFTW_ESTIMATE);

    for (unsigned int i = 0; i < npix; ++i) {
        in[i][0] = 0.0L;
        in[i][1] = 0.0L;
        out[i][0] = 0.0L;
        out[i][1] = 0.0L;
    }

    for (unsigned int i = 0; i <= npix; ++i) {
        for (unsigned int j = 0; j <= npix / 2; ++j) {
            map[i][j] = 0.0L;
        }
    }

    long double** pml = n_matrix_generator(nmod, nmod);

    for (unsigned int j = 1; j < npix / 2; ++j) {

        pml_gen(j, pml);

        in[0][0] = 0.0L;
        in[0][1] = 0.0L;

        for (unsigned int l = 0; l < nmod; ++l) {
            in[0][0] = in[0][0] + cos_ml[0][l] * pml[0][l];
        }

        for (unsigned int m = 1; m < nmod; ++m) {
            in[m][0] = 0.0L;
            in[m][1] = 0.0L;
            for (unsigned int l = m; l < nmod; ++l) {
                in[m][0] = in[m][0] + 2.0L * cos_ml[m][l] * pml[m][l];
            }
        }

        fftwl_execute(p);

        for (unsigned int i = 0; i < npix; ++i) {
            map[i][j] = out[i][0];
        }
        map[npix][j] = out[0][0];

        in[0][0] = 0.0L;
        in[0][1] = 0.0L;

        // бесполезно вообще то
        for (unsigned int l = 0; l < nmod; ++l) {
            in[0][0] = in[0][0] + sin_ml[0][l] * pml[0][l];
        }

        for (unsigned int m = 1; m < nmod; ++m) {
            in[m][0] = 0.0L;
            in[m][1] = 0.0L;
            for (unsigned int l = m; l < nmod; ++l) {
                in[m][0] = in[m][0] + 2.0L * sin_ml[m][l] * pml[m][l];
            }
        }

        fftwl_execute(p);

        for (unsigned int i = 0; i < npix; ++i) {
            map[i][j] = map[i][j] + out[i][1];
        }
        map[npix][j] = map[npix][j] + out[0][1];
    }

    n_matrix_destroyer(pml, nmod);

    fftwl_destroy_plan(p);
    fftwl_free(in);
    fftwl_free(out);
}

void fft_map_forward(long double** map, long double** cos_ml, long double** sin_ml, unsigned int n)
{
    fftwl_complex* in = static_cast<fftwl_complex*>(fftwl_malloc(sizeof(*in) * npix));
    fftwl_complex* out = static_cast<fftwl_complex*>(fftwl_malloc(sizeof(*out) * npix));
    fftwl_plan p = fftwl_plan_dft_1d(int(npix), in, out, FFTW_FORWARD, FFTW_ESTIMATE);

    for (unsigned int i = 0; i < npix; ++i) {
        in[i][0] = 0.0L;
        in[i][1] = 0.0L;
        out[i][0] = 0.0L;
        out[i][1] = 0.0L;
    }

    for (unsigned int i = 0; i <= npix; ++i) {
        for (unsigned int j = 0; j <= npix / 2; ++j) {
            map[i][j] = 0.0L;
        }
    }

    long double** pml = n_matrix_generator(n, n);

    for (unsigned int j = 1; j < npix / 2; ++j) {

        pml_gen(j, pml);

        in[0][0] = 0.0L;
        in[0][1] = 0.0L;

        for (unsigned int l = 0; l < n; ++l) {
            in[0][0] = in[0][0] + cos_ml[0][l] * pml[0][l];
        }

        for (unsigned int m = 1; m < n; ++m) {
            in[m][0] = 0.0L;
            in[m][1] = 0.0L;
            for (unsigned int l = m; l < n; ++l) {
                in[m][0] = in[m][0] + 2.0L * cos_ml[m][l] * pml[m][l];
            }
        }

        fftwl_execute(p);

        for (unsigned int i = 0; i < npix; ++i) {
            map[i][j] = out[i][0];
        }
        map[npix][j] = out[0][0];

        in[0][0] = 0.0L;
        in[0][1] = 0.0L;

        for (unsigned int l = 0; l < n; ++l) {
            in[0][0] = in[0][0] + sin_ml[0][l] * pml[0][l];
        }

        for (unsigned int m = 1; m < n; ++m) {
            in[m][0] = 0.0L;
            in[m][1] = 0.0L;
            for (unsigned int l = m; l < n; ++l) {
                in[m][0] = in[m][0] + 2.0L * sin_ml[m][l] * pml[m][l];
            }
        }

        fftwl_execute(p);

        for (unsigned int i = 0; i < npix; ++i) {
            map[i][j] = map[i][j] + out[i][1];
        }
        map[npix][j] = map[npix][j] + out[0][1];
    }

    n_matrix_destroyer(pml, n);

    fftwl_destroy_plan(p);
    fftwl_free(in);
    fftwl_free(out);
}


long double fft_point_forward(long double phi, long double theta, long double** cos_ml, long double** sin_ml) {

    long double f = 0.0L;
    long double long_m;

    long double** pml = n_matrix_generator(nmod, nmod);

    pml_gen(theta, pml);

    for (unsigned int m = 1; m < nmod; ++m) {

        long_m = static_cast<long double>(m);

        for (unsigned int l = m; l < nmod; ++l) {
            f = f + 2.0L * (cos_ml[m][l] * pml[m][l] * cosl(long_m * phi) -
                 sin_ml[m][l] * pml[m][l] * sinl(long_m * phi));
        }
    }

    for (unsigned int l = 0; l < nmod; ++l)
        f = f + cos_ml[0][l] * pml[0][l];

    n_matrix_destroyer(pml, nmod);

    return f;
}

long double fft_point_forward(long double phi, long double theta, long double** cos_ml, long double** sin_ml,
                                unsigned int n) {

    long double f = 0.0L;
    long double long_m;

    long double** pml = n_matrix_generator(n, n);

    pml_gen(theta, pml);

    for (unsigned int m = 1; m < n; ++m) {

        long_m = static_cast<long double>(m);

        for (unsigned int l = m; l < n; ++l) {
            f = f + 2.0L * (cos_ml[m][l] * pml[m][l] * cosl(long_m * phi) -
                         sin_ml[m][l] * pml[m][l] * sinl(long_m * phi));
        }
    }

    for (unsigned int l = 0; l < n; ++l)
        f = f + cos_ml[0][l] * pml[0][l];

    n_matrix_destroyer(pml, n);

    return f;
}



void fft_map_backward(long double** map, long double** cos_ml, long double** sin_ml) {

    fftwl_complex* in = static_cast<fftwl_complex*>(fftwl_malloc(sizeof(*in) * npix));
    fftwl_complex* out = static_cast<fftwl_complex*>(fftwl_malloc(sizeof(*out) * npix));
    fftwl_plan p = fftwl_plan_dft_1d(int(npix), in, out, FFTW_FORWARD, FFTW_ESTIMATE);

    long double theta;
    long double long_j;
    long double norm = 0.0L;

    for (unsigned int i = 0; i < npix; ++i) {
        in[i][0] = 0.0L;
        in[i][1] = 0.0L;
        out[i][0] = 0.0L;
        out[i][1] = 0.0L;
    }

    for (unsigned int m = 0; m < nmod; ++m) {
        for (unsigned int l = 0; l < nmod; ++l) {
            cos_ml[m][l] = 0.0L;
            sin_ml[m][l] = 0.0L;
        }
    }

    long double** pml = n_matrix_generator(nmod, nmod);

    for (unsigned int j = 1; j < npix / 2; ++j) {

        long_j = static_cast<long double>(j);

        theta =  map_parameter * long_j;
        norm = norm + sinl(theta);

        pml_gen(j, pml);

        for (unsigned int i = 0; i < npix; ++i) {
            in[i][0] = map[i][j];
        }

        fftwl_execute(p);

        for (unsigned int i = 1; i < npix / 2; ++i) {
            out[i][0] = out[i][0] * 2.0L;
        }

        for (unsigned int i = 1; i < npix / 2 - 1; ++i) {
            out[i][1] = 2.0L * out[i][1];
        }

        for (unsigned int m = 0; m < nmod; ++m) {
            for (unsigned int l = 0; l < nmod; ++l) {
                cos_ml[m][l] = cos_ml[m][l] + out[m][0] * pml[m][l] * sinl(theta) * 4.0L * PI;
                sin_ml[m][l] = sin_ml[m][l] + out[m][1] * pml[m][l] * sinl(theta) * 4.0L * PI;
            }
        }
    }

    for (unsigned int m = 1; m < nmod; ++m) {
        for (unsigned int l = 0; l < nmod; ++l) {
            cos_ml[m][l] = cos_ml[m][l] / (norm * long_npix * 2.0L);
            sin_ml[m][l] = sin_ml[m][l] / (norm * long_npix * 2.0L);
        }
    }

    for (unsigned int l = 0; l < nmod; ++l) {
        cos_ml[0][l] = cos_ml[0][l] / (norm * long_npix);
        sin_ml[0][l] = sin_ml[0][l] / (norm * long_npix);
    }

    n_matrix_destroyer(pml, nmod);

    fftwl_destroy_plan(p);
    fftwl_free(in);
    fftwl_free(out);
}

void fft_map_backward(long double** map, long double** cos_ml, long double** sin_ml, unsigned int n) {

    fftwl_complex* in = static_cast<fftwl_complex*>(fftwl_malloc(sizeof(*in) * npix));
    fftwl_complex* out = static_cast<fftwl_complex*>(fftwl_malloc(sizeof(*out) * npix));
    fftwl_plan p = fftwl_plan_dft_1d(int(npix), in, out, FFTW_FORWARD, FFTW_ESTIMATE);

    long double theta;
    long double long_j;
    long double norm = 0.0L;

    for (unsigned int i = 0; i < npix; ++i) {
        in[i][0] = 0.0L;
        in[i][1] = 0.0L;
        out[i][0] = 0.0L;
        out[i][1] = 0.0L;
    }

    for (unsigned int m = 0; m < n; ++m) {
        for (unsigned int l = 0; l < n; ++l) {
            cos_ml[m][l] = 0.0L;
            sin_ml[m][l] = 0.0L;
        }
    }

    long double** pml = n_matrix_generator(n, n);

    for (unsigned int j = 1; j < npix / 2; ++j) {

        long_j = static_cast<long double>(j);

        theta =  map_parameter * long_j;
        norm = norm + sinl(theta);

        pml_gen(j, pml);

        for (unsigned int i = 0; i < npix; ++i) {
            in[i][0] = map[i][j];
        }

        fftwl_execute(p);

        for (unsigned int i = 1; i < npix / 2; ++i) {
            out[i][0] = out[i][0] * 2.0L;
        }

        for (unsigned int i = 1; i < npix / 2 - 1; ++i) {
            out[i][1] = 2.0L * out[i][1];
        }

        for (unsigned int m = 0; m < n; ++m) {
            for (unsigned int l = 0; l < n; ++l) {
                cos_ml[m][l] = cos_ml[m][l] + out[m][0] * pml[m][l] * sinl(theta) * 4.0L * PI;
                sin_ml[m][l] = sin_ml[m][l] + out[m][1] * pml[m][l] * sinl(theta) * 4.0L * PI;
            }
        }
    }

    for (unsigned int m = 1; m < n; ++m) {
        for (unsigned int l = 0; l < n; ++l) {
            cos_ml[m][l] = cos_ml[m][l] / (norm * long_npix * 2.0L);
            sin_ml[m][l] = sin_ml[m][l] / (norm * long_npix * 2.0L);
        }
    }

    for (unsigned int l = 0; l < n; ++l) {
        cos_ml[0][l] = cos_ml[0][l] / (norm * long_npix);
        sin_ml[0][l] = sin_ml[0][l] / (norm * long_npix);
    }

    n_matrix_destroyer(pml, n);

    fftwl_destroy_plan(p);
    fftwl_free(in);
    fftwl_free(out);
}



void fft_map_x_forward(long double** map, long double** cos_ml, long double** sin_ml)
{
    fftwl_complex* in = static_cast<fftwl_complex*>(fftwl_malloc(sizeof(*in) * npix));
    fftwl_complex* out = static_cast<fftwl_complex*>(fftwl_malloc(sizeof(*out) * npix));
    fftwl_plan p = fftwl_plan_dft_1d(int(npix), in, out, FFTW_FORWARD, FFTW_ESTIMATE);

    for (unsigned int i = 0; i < npix; ++i) {
        in[i][0] = 0.0L;
        in[i][1] = 0.0L;
        out[i][0] = 0.0L;
        out[i][1] = 0.0L;
    }

    for (unsigned int i = 0; i <= npix; ++i) {
        for (unsigned int j = 0; j <= npix / 2; ++j) {
            map[i][j] = 0.0L;
        }
    }

    long double** pml = n_matrix_generator(nmod, nmod);

    for (unsigned int j = 1; j < npix / 2; ++j) {

        pml_x_gen(j, pml);

        in[0][0] = 0.0L;
        in[0][1] = 0.0L;

        for (unsigned int l = 0; l < nmod; ++l) {
            in[0][0] = in[0][0] + cos_ml[0][l] * pml[0][l];
        }

        for (unsigned int m = 1; m < nmod; ++m) {
            in[m][0] = 0.0L;
            in[m][1] = 0.0L;
            for (unsigned int l = m; l < nmod; ++l) {
                in[m][0] = in[m][0] + 2.0L * cos_ml[m][l] * pml[m][l];
            }
        }

        fftwl_execute(p);

        for (unsigned int i = 0; i < npix; ++i) {
            map[i][j] = out[i][1];
        }
        map[npix][j] = out[0][1];

        in[0][0] = 0.0L;
        in[0][1] = 0.0L;

        for (unsigned int l = 0; l < nmod; ++l) {
            in[0][0] = in[0][0] + sin_ml[0][l] * pml[0][l];
        }

        for (unsigned int m = 1; m < nmod; ++m) {
            in[m][0] = 0.0L;
            in[m][1] = 0.0L;
            for (unsigned int l = m; l < nmod; ++l) {
                in[m][0] = in[m][0] + 2.0L * sin_ml[m][l] * pml[m][l];
            }
        }

        fftwl_execute(p);

        for (unsigned int i = 0; i < npix; ++i) {
            map[i][j] = map[i][j] - out[i][0];
        }
        map[npix][j] = map[npix][j] - out[0][0];
    }

    n_matrix_destroyer(pml, nmod);

    fftwl_destroy_plan(p);
    fftwl_free(in);
    fftwl_free(out);
}

void fft_map_x_forward(long double** map, long double** cos_ml, long double** sin_ml, unsigned int n)
{
    fftwl_complex* in = static_cast<fftwl_complex*>(fftwl_malloc(sizeof(*in) * npix));
    fftwl_complex* out = static_cast<fftwl_complex*>(fftwl_malloc(sizeof(*out) * npix));
    fftwl_plan p = fftwl_plan_dft_1d(int(npix), in, out, FFTW_FORWARD, FFTW_ESTIMATE);

    for (unsigned int i = 0; i < npix; ++i) {
        in[i][0] = 0.0L;
        in[i][1] = 0.0L;
        out[i][0] = 0.0L;
        out[i][1] = 0.0L;
    }

    for (unsigned int i = 0; i <= npix; ++i) {
        for (unsigned int j = 0; j <= npix / 2; ++j) {
            map[i][j] = 0.0L;
        }
    }

    long double** pml = n_matrix_generator(n, n);

    for (unsigned int j = 1; j < npix / 2; ++j) {

        pml_x_gen(j, pml);

        in[0][0] = 0.0L;
        in[0][1] = 0.0L;

        for (unsigned int l = 0; l < n; ++l) {
            in[0][0] = in[0][0] + cos_ml[0][l] * pml[0][l];
        }

        for (unsigned int m = 1; m < n; ++m) {
            in[m][0] = 0.0L;
            in[m][1] = 0.0L;
            for (unsigned int l = m; l < n; ++l) {
                in[m][0] = in[m][0] + 2.0L * cos_ml[m][l] * pml[m][l];
            }
        }

        fftwl_execute(p);

        for (unsigned int i = 0; i < npix; ++i) {
            map[i][j] = out[i][1];
        }
        map[npix][j] = out[0][1];

        in[0][0] = 0.0L;
        in[0][1] = 0.0L;

        for (unsigned int l = 0; l < n; ++l) {
            in[0][0] = in[0][0] + sin_ml[0][l] * pml[0][l];
        }

        for (unsigned int m = 1; m < n; ++m) {
            in[m][0] = 0.0L;
            in[m][1] = 0.0L;
            for (unsigned int l = m; l < n; ++l) {
                in[m][0] = in[m][0] + 2.0L * sin_ml[m][l] * pml[m][l];
            }
        }

        fftwl_execute(p);

        for (unsigned int i = 0; i < npix; ++i) {
            map[i][j] = map[i][j] - out[i][0];
        }
        map[npix][j] = map[npix][j] - out[0][0];
    }

    n_matrix_destroyer(pml, n);

    fftwl_destroy_plan(p);
    fftwl_free(in);
    fftwl_free(out);
}

long double fft_point_x_forward(long double phi, long double theta, long double** cos_ml, long double** sin_ml) {

    long double f = 0.0L;
    long double long_m;

    long double** pml = n_matrix_generator(nmod, nmod);

    pml_x_gen(theta, pml);

    for (unsigned int m = 1; m < nmod; ++m) {

        long_m = static_cast<long double>(m);

        for (unsigned int l = m; l < nmod; ++l) {
            f = f + 2.0L * (- cos_ml[m][l] * pml[m][l] * sinl(long_m * phi) -
                            sin_ml[m][l] * pml[m][l] * cosl(long_m * phi));
        }
    }

    for (unsigned int l = 0; l < nmod; ++l)
        f = f + cos_ml[0][l] * pml[0][l];

    n_matrix_destroyer(pml, nmod);

    return f;
}

long double fft_point_x_forward(long double phi, long double theta, long double** cos_ml, long double** sin_ml,
                                unsigned int n) {

    long double f = 0.0L;
    long double long_m;

    long double** pml = n_matrix_generator(n, n);

    pml_x_gen(theta, pml);

    for (unsigned int m = 1; m < n; ++m) {

        long_m = static_cast<long double>(m);

        for (unsigned int l = m; l < n; ++l) {
            f = f + 2.0L * (- cos_ml[m][l] * pml[m][l] * sinl(long_m * phi) -
                            sin_ml[m][l] * pml[m][l] * cosl(long_m * phi));
        }
    }

    for (unsigned int l = 0; l < n; ++l)
        f = f + cos_ml[0][l] * pml[0][l];

    n_matrix_destroyer(pml, n);

    return f;
}



void fft_map_y_forward(long double** map, long double** cos_ml, long double** sin_ml)
{
    fftwl_complex* in = static_cast<fftwl_complex*>(fftwl_malloc(sizeof(*in) * npix));
    fftwl_complex* out = static_cast<fftwl_complex*>(fftwl_malloc(sizeof(*out) * npix));
    fftwl_plan p = fftwl_plan_dft_1d(int(npix), in, out, FFTW_FORWARD, FFTW_ESTIMATE);

    for (unsigned int i = 0; i < npix; ++i) {
        in[i][0] = 0.0L;
        in[i][1] = 0.0L;
        out[i][0] = 0.0L;
        out[i][1] = 0.0L;
    }

    for (unsigned int i = 0; i <= npix; ++i) {
        for (unsigned int j = 0; j <= npix / 2; ++j) {
            map[i][j] = 0.0L;
        }
    }

    long double** pml = n_matrix_generator(nmod, nmod);

    for (unsigned int j = 1; j < npix / 2; ++j) {

        pml_y_gen(j, pml);

        in[0][0] = 0.0L;
        in[0][1] = 0.0L;

        for (unsigned int l = 0; l < nmod; ++l) {
            in[0][0] = in[0][0] + cos_ml[0][l] * pml[0][l];
        }

        for (unsigned int m = 1; m < nmod; ++m) {
            in[m][0] = 0.0L;
            in[m][1] = 0.0L;
            for (unsigned int l = m; l < nmod; ++l) {
                in[m][0] = in[m][0] + 2.0L * cos_ml[m][l] * pml[m][l];
            }
        }

        fftwl_execute(p);

        for (unsigned int i = 0; i < npix; ++i) {
            map[i][j] = out[i][0];
        }
        map[npix][j] = out[0][0];

        in[0][0] = 0.0L;
        in[0][1] = 0.0L;

        for (unsigned int l = 0; l < nmod; ++l) {
            in[0][0] = in[0][0] + sin_ml[0][l] * pml[0][l];
        }

        for (unsigned int m = 1; m < nmod; ++m) {
            in[m][0] = 0.0L;
            in[m][1] = 0.0L;
            for (unsigned int l = m; l < nmod; ++l) {
                in[m][0] = in[m][0] + 2.0L * sin_ml[m][l] * pml[m][l];
            }
        }

        fftwl_execute(p);

        for (unsigned int i = 0; i < npix; ++i) {
            map[i][j] = map[i][j] + out[i][1];
        }
        map[npix][j] = map[npix][j] + out[0][1];
    }

    n_matrix_destroyer(pml, nmod);

    fftwl_destroy_plan(p);
    fftwl_free(in);
    fftwl_free(out);
}

void fft_map_y_forward(long double** map, long double** cos_ml, long double** sin_ml, unsigned int n)
{
    fftwl_complex* in = static_cast<fftwl_complex*>(fftwl_malloc(sizeof(*in) * npix));
    fftwl_complex* out = static_cast<fftwl_complex*>(fftwl_malloc(sizeof(*out) * npix));
    fftwl_plan p = fftwl_plan_dft_1d(int(npix), in, out, FFTW_FORWARD, FFTW_ESTIMATE);

    for (unsigned int i = 0; i < npix; ++i) {
        in[i][0] = 0.0L;
        in[i][1] = 0.0L;
        out[i][0] = 0.0L;
        out[i][1] = 0.0L;
    }

    for (unsigned int i = 0; i <= npix; ++i) {
        for (unsigned int j = 0; j <= npix / 2; ++j) {
            map[i][j] = 0.0L;
        }
    }

    long double** pml = n_matrix_generator(n, n);

    for (unsigned int j = 1; j < npix / 2; ++j) {

        pml_y_gen(j, pml);

        in[0][0] = 0.0L;
        in[0][1] = 0.0L;

        for (unsigned int l = 0; l < n; ++l) {
            in[0][0] = in[0][0] + cos_ml[0][l] * pml[0][l];
        }

        for (unsigned int m = 1; m < n; ++m) {
            in[m][0] = 0.0L;
            in[m][1] = 0.0L;
            for (unsigned int l = m; l < n; ++l) {
                in[m][0] = in[m][0] + 2.0L * cos_ml[m][l] * pml[m][l];
            }
        }

        fftwl_execute(p);

        for (unsigned int i = 0; i < npix; ++i) {
            map[i][j] = out[i][0];
        }
        map[npix][j] = out[0][0];

        in[0][0] = 0.0L;
        in[0][1] = 0.0L;

        for (unsigned int l = 0; l < n; ++l) {
            in[0][0] = in[0][0] + sin_ml[0][l] * pml[0][l];
        }

        for (unsigned int m = 1; m < n; ++m) {
            in[m][0] = 0.0L;
            in[m][1] = 0.0L;
            for (unsigned int l = m; l < n; ++l) {
                in[m][0] = in[m][0] + 2.0L * sin_ml[m][l] * pml[m][l];
            }
        }

        fftwl_execute(p);

        for (unsigned int i = 0; i < npix; ++i) {
            map[i][j] = map[i][j] + out[i][1];
        }
        map[npix][j] = map[npix][j] + out[0][1];
    }

    n_matrix_destroyer(pml, n);

    fftwl_destroy_plan(p);
    fftwl_free(in);
    fftwl_free(out);
}

long double fft_point_y_forward(long double phi, long double theta, long double** cos_ml, long double** sin_ml) {

    long double f = 0.0L;
    long double long_m;

    long double** pml = n_matrix_generator(nmod, nmod);

    pml_y_gen(theta, pml);

    for (unsigned int m = 1; m < nmod; ++m) {

        long_m = static_cast<long double>(m);

        for (unsigned int l = m; l < nmod; ++l) {
            f = f + 2.0L * (cos_ml[m][l] * pml[m][l] * cosl(long_m * phi) -
                            sin_ml[m][l] * pml[m][l] * sinl(long_m * phi));
        }
    }

    for (unsigned int l = 0; l < nmod; ++l)
        f = f + cos_ml[0][l] * pml[0][l];

    n_matrix_destroyer(pml, nmod);

    return f;
}

long double fft_point_y_forward(long double phi, long double theta, long double** cos_ml, long double** sin_ml,
                                unsigned int n) {

    long double f = 0.0L;
    long double long_m;

    long double** pml = n_matrix_generator(n, n);

    pml_y_gen(theta, pml);

    for (unsigned int m = 1; m < n; ++m) {

        long_m = static_cast<long double>(m);

        for (unsigned int l = m; l < n; ++l) {
            f = f + 2.0L * (cos_ml[m][l] * pml[m][l] * cosl(long_m * phi) -
                            sin_ml[m][l] * pml[m][l] * sinl(long_m * phi));
        }
    }

    for (unsigned int l = 0; l < n; ++l)
        f = f + cos_ml[0][l] * pml[0][l];

    n_matrix_destroyer(pml, n);

    return f;
}



void fft_map_xx_forward(long double** map, long double** cos_ml, long double** sin_ml)
{
    fftwl_complex* in = static_cast<fftwl_complex*>(fftwl_malloc(sizeof(*in) * npix));
    fftwl_complex* out = static_cast<fftwl_complex*>(fftwl_malloc(sizeof(*out) * npix));
    fftwl_plan p = fftwl_plan_dft_1d(int(npix), in, out, FFTW_FORWARD, FFTW_ESTIMATE);

    for (unsigned int i = 0; i < npix; ++i) {
        in[i][0] = 0.0L;
        in[i][1] = 0.0L;
        out[i][0] = 0.0L;
        out[i][1] = 0.0L;
    }

    for (unsigned int i = 0; i <= npix; ++i) {
        for (unsigned int j = 0; j <= npix / 2; ++j) {
            map[i][j] = 0.0L;
        }
    }

    long double** pml = n_matrix_generator(nmod, nmod);

    for (unsigned int j = 1; j < npix / 2; ++j) {

        pml_xx_gen(j, pml);

        in[0][0] = 0.0L;
        in[0][1] = 0.0L;

        for (unsigned int l = 0; l < nmod; ++l) {
            in[0][0] = in[0][0] + cos_ml[0][l] * pml[0][l];
        }

        for (unsigned int m = 1; m < nmod; ++m) {
            in[m][0] = 0.0L;
            in[m][1] = 0.0L;
            for (unsigned int l = m; l < nmod; ++l) {
                in[m][0] = in[m][0] + 2.0L * cos_ml[m][l] * pml[m][l];
            }
        }

        fftwl_execute(p);

        for (unsigned int i = 0; i < npix; ++i) {
            map[i][j] = out[i][0];
        }
        map[npix][j] = out[0][0];

        in[0][0] = 0.0L;
        in[0][1] = 0.0L;

        for (unsigned int l = 0; l < nmod; ++l) {
            in[0][0] = in[0][0] + sin_ml[0][l] * pml[0][l];
        }

        for (unsigned int m = 1; m < nmod; ++m) {
            in[m][0] = 0.0L;
            in[m][1] = 0.0L;
            for (unsigned int l = m; l < nmod; ++l) {
                in[m][0] = in[m][0] + 2.0L * sin_ml[m][l] * pml[m][l];
            }
        }

        fftwl_execute(p);

        for (unsigned int i = 0; i < npix; ++i) {
            map[i][j] = map[i][j] + out[i][1];
        }
        map[npix][j] = map[npix][j] + out[0][1];
    }

    n_matrix_destroyer(pml, nmod);

    fftwl_destroy_plan(p);
    fftwl_free(in);
    fftwl_free(out);
}

void fft_map_xx_forward(long double** map, long double** cos_ml, long double** sin_ml, unsigned int n)
{
    fftwl_complex* in = static_cast<fftwl_complex*>(fftwl_malloc(sizeof(*in) * npix));
    fftwl_complex* out = static_cast<fftwl_complex*>(fftwl_malloc(sizeof(*out) * npix));
    fftwl_plan p = fftwl_plan_dft_1d(int(npix), in, out, FFTW_FORWARD, FFTW_ESTIMATE);

    for (unsigned int i = 0; i < npix; ++i) {
        in[i][0] = 0.0L;
        in[i][1] = 0.0L;
        out[i][0] = 0.0L;
        out[i][1] = 0.0L;
    }

    for (unsigned int i = 0; i <= npix; ++i) {
        for (unsigned int j = 0; j <= npix / 2; ++j) {
            map[i][j] = 0.0L;
        }
    }

    long double** pml = n_matrix_generator(n, n);

    for (unsigned int j = 1; j < npix / 2; ++j) {

        pml_xx_gen(j, pml);

        in[0][0] = 0.0L;
        in[0][1] = 0.0L;

        for (unsigned int l = 0; l < n; ++l) {
            in[0][0] = in[0][0] + cos_ml[0][l] * pml[0][l];
        }

        for (unsigned int m = 1; m < n; ++m) {
            in[m][0] = 0.0L;
            in[m][1] = 0.0L;
            for (unsigned int l = m; l < n; ++l) {
                in[m][0] = in[m][0] + 2.0L * cos_ml[m][l] * pml[m][l];
            }
        }

        fftwl_execute(p);

        for (unsigned int i = 0; i < npix; ++i) {
            map[i][j] = out[i][0];
        }
        map[npix][j] = out[0][0];

        in[0][0] = 0.0L;
        in[0][1] = 0.0L;

        for (unsigned int l = 0; l < n; ++l) {
            in[0][0] = in[0][0] + sin_ml[0][l] * pml[0][l];
        }

        for (unsigned int m = 1; m < n; ++m) {
            in[m][0] = 0.0L;
            in[m][1] = 0.0L;
            for (unsigned int l = m; l < n; ++l) {
                in[m][0] = in[m][0] + 2.0L * sin_ml[m][l] * pml[m][l];
            }
        }

        fftwl_execute(p);

        for (unsigned int i = 0; i < npix; ++i) {
            map[i][j] = map[i][j] + out[i][1];
        }
        map[npix][j] = map[npix][j] + out[0][1];
    }

    n_matrix_destroyer(pml, n);

    fftwl_destroy_plan(p);
    fftwl_free(in);
    fftwl_free(out);
}


long double fft_point_xx_forward(long double phi, long double theta, long double** cos_ml, long double** sin_ml) {

    long double f = 0.0L;
    long double long_m;

    long double** pml = n_matrix_generator(nmod, nmod);

    pml_xx_gen(theta, pml);

    for (unsigned int m = 1; m < nmod; ++m) {

        long_m = static_cast<long double>(m);

        for (unsigned int l = m; l < nmod; ++l) {
            f = f + 2.0L * (cos_ml[m][l] * pml[m][l] * cosl(long_m * phi) -
                            sin_ml[m][l] * pml[m][l] * sinl(long_m * phi));
        }
    }

    for (unsigned int l = 0; l < nmod; ++l)
        f = f + cos_ml[0][l] * pml[0][l];

    n_matrix_destroyer(pml, nmod);

    return f;
}

long double fft_point_xx_forward(long double phi, long double theta, long double** cos_ml, long double** sin_ml,
                                 unsigned int n) {

    long double f = 0.0L;
    long double long_m;

    long double** pml = n_matrix_generator(n, n);

    pml_xx_gen(theta, pml);

    for (unsigned int m = 1; m < n; ++m) {

        long_m = static_cast<long double>(m);

        for (unsigned int l = m; l < n; ++l) {
            f = f + 2.0L * (cos_ml[m][l] * pml[m][l] * cosl(long_m * phi) -
                            sin_ml[m][l] * pml[m][l] * sinl(long_m * phi));
        }
    }

    for (unsigned int l = 0; l < n; ++l)
        f = f + cos_ml[0][l] * pml[0][l];

    n_matrix_destroyer(pml, n);

    return f;
}


void fft_map_yy_forward(long double** map, long double** cos_ml, long double** sin_ml)
{
    fftwl_complex* in = static_cast<fftwl_complex*>(fftwl_malloc(sizeof(*in) * npix));
    fftwl_complex* out = static_cast<fftwl_complex*>(fftwl_malloc(sizeof(*out) * npix));
    fftwl_plan p = fftwl_plan_dft_1d(int(npix), in, out, FFTW_FORWARD, FFTW_ESTIMATE);

    for (unsigned int i = 0; i < npix; ++i) {
        in[i][0] = 0.0L;
        in[i][1] = 0.0L;
        out[i][0] = 0.0L;
        out[i][1] = 0.0L;
    }

    for (unsigned int i = 0; i <= npix; ++i) {
        for (unsigned int j = 0; j <= npix / 2; ++j) {
            map[i][j] = 0.0L;
        }
    }

    long double** pml = n_matrix_generator(nmod, nmod);

    for (unsigned int j = 1; j < npix / 2; ++j) {

        pml_yy_gen(j, pml);

        in[0][0] = 0.0L;
        in[0][1] = 0.0L;

        for (unsigned int l = 0; l < nmod; ++l) {
            in[0][0] = in[0][0] + cos_ml[0][l] * pml[0][l];
        }

        for (unsigned int m = 1; m < nmod; ++m) {
            in[m][0] = 0.0L;
            in[m][1] = 0.0L;
            for (unsigned int l = m; l < nmod; ++l) {
                in[m][0] = in[m][0] + 2.0L * cos_ml[m][l] * pml[m][l];
            }
        }

        fftwl_execute(p);

        for (unsigned int i = 0; i < npix; ++i) {
            map[i][j] = out[i][0];
        }
        map[npix][j] = out[0][0];

        in[0][0] = 0.0L;
        in[0][1] = 0.0L;

        for (unsigned int l = 0; l < nmod; ++l) {
            in[0][0] = in[0][0] + sin_ml[0][l] * pml[0][l];
        }

        for (unsigned int m = 1; m < nmod; ++m) {
            in[m][0] = 0.0L;
            in[m][1] = 0.0L;
            for (unsigned int l = m; l < nmod; ++l) {
                in[m][0] = in[m][0] + 2.0L * sin_ml[m][l] * pml[m][l];
            }
        }

        fftwl_execute(p);

        for (unsigned int i = 0; i < npix; ++i) {
            map[i][j] = map[i][j] + out[i][1];
        }
        map[npix][j] = map[npix][j] + out[0][1];
    }

    n_matrix_destroyer(pml, nmod);

    fftwl_destroy_plan(p);
    fftwl_free(in);
    fftwl_free(out);
}

void fft_map_yy_forward(long double** map, long double** cos_ml, long double** sin_ml, unsigned int n)
{
    fftwl_complex* in = static_cast<fftwl_complex*>(fftwl_malloc(sizeof(*in) * npix));
    fftwl_complex* out = static_cast<fftwl_complex*>(fftwl_malloc(sizeof(*out) * npix));
    fftwl_plan p = fftwl_plan_dft_1d(int(npix), in, out, FFTW_FORWARD, FFTW_ESTIMATE);

    for (unsigned int i = 0; i < npix; ++i) {
        in[i][0] = 0.0L;
        in[i][1] = 0.0L;
        out[i][0] = 0.0L;
        out[i][1] = 0.0L;
    }

    for (unsigned int i = 0; i <= npix; ++i) {
        for (unsigned int j = 0; j <= npix / 2; ++j) {
            map[i][j] = 0.0L;
        }
    }

    long double** pml = n_matrix_generator(n, n);

    for (unsigned int j = 1; j < npix / 2; ++j) {

        pml_yy_gen(j, pml);

        in[0][0] = 0.0L;
        in[0][1] = 0.0L;

        for (unsigned int l = 0; l < n; ++l) {
            in[0][0] = in[0][0] + cos_ml[0][l] * pml[0][l];
        }

        for (unsigned int m = 1; m < n; ++m) {
            in[m][0] = 0.0L;
            in[m][1] = 0.0L;
            for (unsigned int l = m; l < n; ++l) {
                in[m][0] = in[m][0] + 2.0L * cos_ml[m][l] * pml[m][l];
            }
        }

        fftwl_execute(p);

        for (unsigned int i = 0; i < npix; ++i) {
            map[i][j] = out[i][0];
        }
        map[npix][j] = out[0][0];

        in[0][0] = 0.0L;
        in[0][1] = 0.0L;

        for (unsigned int l = 0; l < n; ++l) {
            in[0][0] = in[0][0] + sin_ml[0][l] * pml[0][l];
        }

        for (unsigned int m = 1; m < n; ++m) {
            in[m][0] = 0.0L;
            in[m][1] = 0.0L;
            for (unsigned int l = m; l < n; ++l) {
                in[m][0] = in[m][0] + 2.0L * sin_ml[m][l] * pml[m][l];
            }
        }

        fftwl_execute(p);

        for (unsigned int i = 0; i < npix; ++i) {
            map[i][j] = map[i][j] + out[i][1];
        }
        map[npix][j] = map[npix][j] + out[0][1];
    }

    n_matrix_destroyer(pml, n);

    fftwl_destroy_plan(p);
    fftwl_free(in);
    fftwl_free(out);
}


long double fft_point_yy_forward(long double phi, long double theta, long double** cos_ml, long double** sin_ml) {

    long double f = 0.0L;
    long double long_m;

    long double** pml = n_matrix_generator(nmod, nmod);

    pml_yy_gen(theta, pml);

    for (unsigned int m = 1; m < nmod; ++m) {

        long_m = static_cast<long double>(m);

        for (unsigned int l = m; l < nmod; ++l) {
            f = f + 2.0L * (cos_ml[m][l] * pml[m][l] * cosl(long_m * phi) -
                            sin_ml[m][l] * pml[m][l] * sinl(long_m * phi));
        }
    }

    for (unsigned int l = 0; l < nmod; ++l)
        f = f + cos_ml[0][l] * pml[0][l];

    n_matrix_destroyer(pml, nmod);

    return f;
}

long double fft_point_yy_forward(long double phi, long double theta, long double** cos_ml, long double** sin_ml,
                                 unsigned int n) {

    long double f = 0.0L;
    long double long_m;

    long double** pml = n_matrix_generator(n, n);

    pml_yy_gen(theta, pml);

    for (unsigned int m = 1; m < n; ++m) {

        long_m = static_cast<long double>(m);

        for (unsigned int l = m; l < n; ++l) {
            f = f + 2.0L * (cos_ml[m][l] * pml[m][l] * cosl(long_m * phi) -
                            sin_ml[m][l] * pml[m][l] * sinl(long_m * phi));
        }
    }

    for (unsigned int l = 0; l < n; ++l)
        f = f + cos_ml[0][l] * pml[0][l];

    n_matrix_destroyer(pml, n);

    return f;
}


void fft_map_xy_forward(long double** map, long double** cos_ml, long double** sin_ml)
{

    fftwl_complex* in = static_cast<fftwl_complex*>(fftwl_malloc(sizeof(*in) * npix));
    fftwl_complex* out = static_cast<fftwl_complex*>(fftwl_malloc(sizeof(*out) * npix));
    fftwl_plan p = fftwl_plan_dft_1d(int(npix), in, out, FFTW_FORWARD, FFTW_ESTIMATE);

    for (unsigned int i = 0; i < npix; ++i) {
        in[i][0] = 0.0L;
        in[i][1] = 0.0L;
        out[i][0] = 0.0L;
        out[i][1] = 0.0L;
    }

    for (unsigned int i = 0; i <= npix; ++i) {
        for (unsigned int j = 0; j <= npix / 2; ++j) {
            map[i][j] = 0.0L;
        }
    }

    long double** pml = n_matrix_generator(nmod, nmod);

    for (unsigned int j = 1; j < npix / 2; ++j) {

        pml_xy_gen(j, pml);

        in[0][0] = 0.0L;
        in[0][1] = 0.0L;

        for (unsigned int l = 0; l < nmod; ++l) {
            in[0][0] = in[0][0] + cos_ml[0][l] * pml[0][l];
        }

        for (unsigned int m = 1; m < nmod; ++m) {
            in[m][0] = 0.0L;
            in[m][1] = 0.0L;
            for (unsigned int l = m; l < nmod; ++l) {
                in[m][0] = in[m][0] + 2.0L * cos_ml[m][l] * pml[m][l];
            }
        }

        fftwl_execute(p);

        for (unsigned int i = 0; i < npix; ++i) {
            map[i][j] = out[i][1];
        }
        map[npix][j] = out[0][1];

        in[0][0] = 0.0L;
        in[0][1] = 0.0L;

        for (unsigned int l = 0; l < nmod; ++l) {
            in[0][0] = in[0][0] + sin_ml[0][l] * pml[0][l];
        }

        for (unsigned int m = 1; m < nmod; ++m) {
            in[m][0] = 0.0L;
            in[m][1] = 0.0L;
            for (unsigned int l = m; l < nmod; ++l) {
                in[m][0] = in[m][0] + 2.0L * sin_ml[m][l] * pml[m][l];
            }
        }

        fftwl_execute(p);

        for (unsigned int i = 0; i < npix; ++i) {
            map[i][j] = map[i][j] - out[i][0];
        }
        map[npix][j] = map[npix][j] - out[0][0];
    }

    n_matrix_destroyer(pml, nmod);

    fftwl_destroy_plan(p);
    fftwl_free(in);
    fftwl_free(out);
}

void fft_map_xy_forward(long double** map, long double** cos_ml, long double** sin_ml, unsigned int n)
{

    fftwl_complex* in = static_cast<fftwl_complex*>(fftwl_malloc(sizeof(*in) * npix));
    fftwl_complex* out = static_cast<fftwl_complex*>(fftwl_malloc(sizeof(*out) * npix));
    fftwl_plan p = fftwl_plan_dft_1d(int(npix), in, out, FFTW_FORWARD, FFTW_ESTIMATE);

    for (unsigned int i = 0; i < npix; ++i) {
        in[i][0] = 0.0L;
        in[i][1] = 0.0L;
        out[i][0] = 0.0L;
        out[i][1] = 0.0L;
    }

    for (unsigned int i = 0; i <= npix; ++i) {
        for (unsigned int j = 0; j <= npix / 2; ++j) {
            map[i][j] = 0.0L;
        }
    }

    long double** pml = n_matrix_generator(n, n);

    for (unsigned int j = 1; j < npix / 2; ++j) {

        pml_xy_gen(j, pml);

        in[0][0] = 0.0L;
        in[0][1] = 0.0L;

        for (unsigned int l = 0; l < n; ++l) {
            in[0][0] = in[0][0] + cos_ml[0][l] * pml[0][l];
        }

        for (unsigned int m = 1; m < n; ++m) {
            in[m][0] = 0.0L;
            in[m][1] = 0.0L;
            for (unsigned int l = m; l < n; ++l) {
                in[m][0] = in[m][0] + 2.0L * cos_ml[m][l] * pml[m][l];
            }
        }

        fftwl_execute(p);

        for (unsigned int i = 0; i < npix; ++i) {
            map[i][j] = out[i][1];
        }
        map[npix][j] = out[0][1];

        in[0][0] = 0.0L;
        in[0][1] = 0.0L;

        for (unsigned int l = 0; l < n; ++l) {
            in[0][0] = in[0][0] + sin_ml[0][l] * pml[0][l];
        }

        for (unsigned int m = 1; m < n; ++m) {
            in[m][0] = 0.0L;
            in[m][1] = 0.0L;
            for (unsigned int l = m; l < n; ++l) {
                in[m][0] = in[m][0] + 2.0L * sin_ml[m][l] * pml[m][l];
            }
        }

        fftwl_execute(p);

        for (unsigned int i = 0; i < npix; ++i) {
            map[i][j] = map[i][j] - out[i][0];
        }
        map[npix][j] = map[npix][j] - out[0][0];
    }

    n_matrix_destroyer(pml, n);

    fftwl_destroy_plan(p);
    fftwl_free(in);
    fftwl_free(out);
}


long double fft_point_xy_forward(long double phi, long double theta, long double** cos_ml, long double** sin_ml) {

    long double f = 0.0L;
    long double long_m;

    long double** pml = n_matrix_generator(nmod, nmod);

    pml_xy_gen(theta, pml);

    for (unsigned int m = 1; m < nmod; ++m) {

        long_m = static_cast<long double>(m);

        for (unsigned int l = m; l < nmod; ++l) {
            f = f + 2.0L * (- cos_ml[m][l] * pml[m][l] * sinl(long_m * phi) -
                            sin_ml[m][l] * pml[m][l] * cosl(long_m * phi));
        }
    }

    for (unsigned int l = 0; l < nmod; ++l)
        f = f + cos_ml[0][l] * pml[0][l];

    n_matrix_destroyer(pml, nmod);

    return f;
}

long double fft_point_xy_forward(long double phi, long double theta, long double** cos_ml, long double** sin_ml,
                                 unsigned int n) {

    long double f = 0.0L;
    long double long_m;

    long double** pml = n_matrix_generator(n, n);

    pml_xy_gen(theta, pml);

    for (unsigned int m = 1; m < n; ++m) {

        long_m = static_cast<long double>(m);

        for (unsigned int l = m; l < n; ++l) {
            f = f + 2.0L * (- cos_ml[m][l] * pml[m][l] * sinl(long_m * phi) -
                            sin_ml[m][l] * pml[m][l] * cosl(long_m * phi));
        }
    }

    for (unsigned int l = 0; l < n; ++l)
        f = f + cos_ml[0][l] * pml[0][l];

    n_matrix_destroyer(pml, n);

    return f;
}
