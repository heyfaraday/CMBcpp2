#pragma once

void aml_gasdev(std::default_random_engine generator,
                long double** cos_ml, long double** sin_ml, long double mean, long double std);

void aml_from_cl(std::default_random_engine generator,
                 long double** cos_ml, long double** sin_ml, long double* cl);

void healpix_transform(long double** aml);

void add_monopole(long double** cos_ml, long double** sin_ml, long double aml_0_0);

void remove_monopole(long double** cos_ml, long double** sin_ml);

void add_dipole(long double** cos_ml, long double** sin_ml,
                long double aml_cos_0_1, long double aml_cos_1_1, long double aml_sin_1_1);

void remove_dipole(long double** cos_ml, long double** sin_ml);
