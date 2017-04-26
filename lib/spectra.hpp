#pragma once

long double sigma_0_map(long double** mapx);

long double sigma_0_cl(long double* cl);

long double sigma_1_cl(long double* cl);

long double sigma_2_cl(long double* cl);

long double sigma_0_aml(long double** cos_ml, long double** sin_ml);

long double sigma_1_aml(long double** cos_ml, long double** sin_ml);

long double sigma_2_aml(long double** cos_ml, long double** sin_ml);

void aml_to_healpix_cl(long double** cos_ml, long double** sin_ml, long double* cl);

void aml_to_cl(long double** cos_ml, long double** sin_ml, long double* cl);
