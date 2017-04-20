#pragma once

void pml_gen(unsigned int j, long double **array);

void pml_gen(long double theta, long double **array);

void pml_x_gen(unsigned int j, long double** array);

void pml_x_gen(long double theta, long double** array);

void pml_y_gen(unsigned int j, long double** array);

void pml_y_gen(long double theta, long double** array);

void pml_xx_gen(unsigned int j, long double** array);

void pml_xx_gen(long double theta, long double** array);

void pml_yy_gen(unsigned int j, long double** array);

void pml_yy_gen(long double theta, long double** array);

void pml_xy_gen(unsigned int j, long double** array);

void pml_xy_gen(long double theta, long double** array);

void x_1_ml_gen(unsigned int j, long double** array);

void x_2_ml_gen(unsigned int j, long double** array);
