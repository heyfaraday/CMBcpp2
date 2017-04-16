#include <fftw3.h>
#include <math.h>

#include "fft.hpp"

#include "parameters.hpp"
#include "constants.hpp"
#include "pml.hpp"
#include "utils.hpp"

void fft_map_forward(long double** map, long double** cos_ml, long double** sin_ml)
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

        pml_gen(j, pml);

        for (unsigned int m = 0; m < nmod; ++m) {
            in[m][0] = 0.0L;
            in[m][1] = 0.0L;
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
            in[m][0] = 0.0L;
            in[m][1] = 0.0L;
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

    long double long_i = static_cast<long double>(i);

    long double f = 0.0L;
    long double phi = long_i * 2.0L * PI / long_npix;

    long double** pml = n_matrix_generator(nmod, nmod);

    pml_gen(j, pml);

    for (unsigned int m = 0; m < nmod; ++m) {

        long double long_m = static_cast<long double>(m);

        for (unsigned int l = m; l < nmod; ++l) {
            f += cos_ml[m][l] * pml[m][l] * cosl(long_m * phi) +
                 sin_ml[m][l] * pml[m][l] * sinl(long_m * phi);
        }
    }

    n_matrix_destroyer(pml, nmod);

    return f;
}

void fft_map_backward(long double** map, long double** cos_ml, long double** sin_ml) {

    fftwl_complex* in = static_cast<fftwl_complex*>(fftwl_malloc(sizeof(*in) * npix));
    fftwl_complex* out = static_cast<fftwl_complex*>(fftwl_malloc(sizeof(*out) * npix));
    fftwl_plan p = fftwl_plan_dft_1d(int(npix), in, out, FFTW_FORWARD, FFTW_ESTIMATE);

    long double theta;
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

    for (unsigned int j = 0; j <= npix/2; ++j) {

        long double long_j = static_cast<long double>(j);

        theta = 2.0L * long_j * PI / long_npix;
        norm += sinl(theta);

        pml_gen(j, pml);

        for (unsigned int i = 0; i < npix; ++i) {
            in[i][0] = map[i][j];
        }

        fftwl_execute(p);

        for (unsigned int i = 1; i < npix / 2; ++i) {
            out[i][0] *= 2.0L;
        }

        for (unsigned int i = 1; i < npix / 2 - 1; ++i) {
            out[i][1] *= 2.0L;
        }

        for (unsigned int m = 0; m < nmod; ++m) {
            for (unsigned int l = 0; l < nmod; ++l) {
                cos_ml[m][l] += out[m][0] * pml[m][l] * sinl(theta) * 4.0L * PI;
                sin_ml[m][l] += out[m][1] * pml[m][l] * sinl(theta) * 4.0L * PI;
            }
        }
    }

    for (unsigned int m = 1; m < nmod; ++m) {
        for (unsigned int l = 0; l < nmod; ++l) {
            cos_ml[m][l] /= (norm * long_npix * 2.0L);
            sin_ml[m][l] /= (norm * long_npix * 2.0L);
        }
    }

    for (unsigned int l = 0; l < nmod; ++l) {
        cos_ml[0][l] /= (norm * long_npix);
        sin_ml[0][l] /= (norm * long_npix);
    }

    n_matrix_destroyer(pml, nmod);

    fftwl_destroy_plan(p);
    fftwl_free(in);
    fftwl_free(out);
}
