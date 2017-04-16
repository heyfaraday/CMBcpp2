#include <fftw3.h>
#include <math.h>

#include "fft.hpp"

#include "parameters.hpp"
#include "constants.hpp"
#include "pml.hpp"
#include "utils.hpp"

void fft_map_forward(long double** map, long double** cos_ml, long double** sin_ml)
{

    fftwl_complex* in = (fftwl_complex*) fftwl_malloc(sizeof(*in) * int(npix));
    fftwl_complex* out = (fftwl_complex*) fftwl_malloc(sizeof(*out) * int(npix));
    fftwl_plan p = fftwl_plan_dft_1d(int(npix), in, out, FFTW_FORWARD, FFTW_ESTIMATE);

    for (unsigned int i = 0; i < npix; ++i) {
        in[i][0] = 0.0;
        in[i][1] = 0.0;
        out[i][0] = 0.0;
        out[i][1] = 0.0;
    }

    for (unsigned int i = 0; i <= npix; ++i) {
        for (unsigned int j = 0; j <= npix / 2; ++j) {
            map[i][j] = 0.0;
        }
    }

    long double** pml = n_matrix_generator(nmod, nmod);

    for (unsigned int j = 1; j < npix/2; ++j) {

        pml_gen(j, pml);

        for (unsigned int m = 0; m < nmod; ++m) {
            in[m][0] = 0.0;
            in[m][1] = 0.0;
            for (unsigned int l = m; l < nmod; ++l) {
                in[m][0] += cos_ml[m][l] * pml[m][l];
            }
        }

        fftwl_execute(p);

        for (unsigned int i = 0; i < npix; ++i) {
            map[i][j] = out[i][0];
        }

        map[npix][j] = out[0][0];

        for (unsigned int m = 0; m < nmod; ++m) {
            in[m][0] = 0.0;
            in[m][1] = 0.0;
            for (unsigned int l = m; l < nmod; ++l) {
                in[m][0] += sin_ml[m][l] * pml[m][l];
            }
        }

        fftwl_execute(p);

        for (unsigned int i = 0; i < npix; ++i) {
            map[i][j] += out[i][1];
        }
        map[npix][j] += out[0][1];
    }

    n_matrix_destroyer(pml, nmod);

    fftwl_destroy_plan(p);
    fftwl_free(in);
    fftwl_free(out);
}

long double fft_point_forward(unsigned int i, unsigned int j, long double** cos_ml, long double** sin_ml) {

    double f = 0.0;
    double phi = i * 2.0 * PI / npix;

    long double** pml = n_matrix_generator(nmod, nmod);

    pml_gen(j, pml);

    for (unsigned int m = 0; m < nmod; ++m) {
        for (unsigned int l = m; l < nmod; ++l) {
            f += cos_ml[m][l] * pml[m][l] * cos(m * phi) +
                 sin_ml[m][l] * pml[m][l] * sin(m * phi);
        }
    }

    n_matrix_destroyer(pml, nmod);

    return f;
}

void fft_map_backward(long double** map, long double** cos_ml, long double** sin_ml) {

    fftwl_complex* in = (fftwl_complex*) fftwl_malloc(sizeof(*in) * int(npix));
    fftwl_complex* out = (fftwl_complex*) fftwl_malloc(sizeof(*out) * int(npix));
    fftwl_plan p = fftwl_plan_dft_1d(int(npix), in, out, FFTW_BACKWARD, FFTW_ESTIMATE);

    double theta;
    double norm = 0.0;

    for (unsigned int i = 0; i < npix; ++i) {
        in[i][0] = 0.0;
        in[i][1] = 0.0;
        out[i][0] = 0.0;
        out[i][1] = 0.0;
    }

    for (unsigned int m = 0; m < nmod; ++m) {
        for (unsigned int l = 0; l < nmod; ++l) {
            cos_ml[m][l] = 0.0;
            sin_ml[m][l] = 0.0;
        }
    }

    long double** pml = n_matrix_generator(nmod, nmod);

    for (unsigned int j = 0; j <= npix/2; ++j) {

        theta = 2.0 * j * PI / npix;
        norm += sin(theta);

        pml_gen(j, pml);

        for (unsigned int i = 0; i < npix; ++i) {
            in[i][0] = map[i][j];
        }

        fftwl_execute(p);

        for (unsigned int i = 1; i < npix / 2; ++i) {
            out[i][0] *= 2.0;
        }

        for (unsigned int i = 1; i < npix / 2 - 1; ++i) {
            out[i][1] *= 2.0;
        }

        for (unsigned int m = 0; m < nmod; ++m) {
            for (unsigned int l = 0; l < nmod; ++l) {
                cos_ml[m][l] += out[m][0] * pml[m][l] * sin(theta) * 4.0 * PI;
                sin_ml[m][l] += out[m][1] * pml[m][l] * sin(theta) * 4.0 * PI;
            }
        }
    }

    for (unsigned int m = 1; m < nmod; ++m) {
        for (unsigned int l = 0; l < nmod; ++l) {
            cos_ml[m][l] /= (norm * npix * 2.0);
            sin_ml[m][l] /= (norm * npix * 2.0);
        }
    }

    for (unsigned int l = 0; l < nmod; ++l) {
        cos_ml[0][l] /= (norm * npix);
        sin_ml[0][l] /= (norm * npix);
    }

    n_matrix_destroyer(pml, nmod);

    fftwl_destroy_plan(p);
    fftwl_free(in);
    fftwl_free(out);
}
