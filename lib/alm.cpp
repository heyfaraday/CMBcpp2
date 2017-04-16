#include "alm.hpp"

void alm_gasdev(double** cos_lm, double** sin_lm, int nmod, int seed, double mean, double std)
{

    for (int i = 0; i < nmod; ++i) {
        for (int j = 0; j < nmod; ++j) {
            cos_lm[i][j] = 0.0;
            sin_lm[i][j] = 0.0;
        }
    }

    cos_lm[0][1] = 1.0;
}
