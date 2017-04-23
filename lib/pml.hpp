#pragma once

void pml_gen(unsigned int j, long double **array);

void pml_gen(long double theta, long double **array);

void pml_gen(unsigned int j, long double** array, unsigned int n);

void pml_gen(long double theta, long double** array, unsigned int n);



long double coef1(unsigned int l, unsigned int m);



void pml_x_gen(unsigned int j, long double** array);

void pml_x_gen(long double theta, long double** array);

void pml_x_gen(unsigned int j, long double** array, unsigned int n);

void pml_x_gen(long double theta, long double** array, unsigned int n);



void pml_y_gen(unsigned int j, long double** array);

void pml_y_gen(long double theta, long double** array);

void pml_y_gen(unsigned int j, long double** array, unsigned int n);

void pml_y_gen(long double theta, long double** array, unsigned int n);



void pml_xx_gen(unsigned int j, long double** array);

void pml_xx_gen(long double theta, long double** array);

void pml_xx_gen(unsigned int j, long double** array, unsigned int n);

void pml_xx_gen(long double theta, long double** array, unsigned int n);


void pml_yy_gen(unsigned int j, long double** array);

void pml_yy_gen(long double theta, long double** array);

void pml_yy_gen(unsigned int j, long double** array, unsigned int n);

void pml_yy_gen(long double theta, long double** array, unsigned int n);



void pml_xy_gen(unsigned int j, long double** array);

void pml_xy_gen(long double theta, long double** array);

void pml_xy_gen(unsigned int j, long double** array, unsigned int n);

void pml_xy_gen(long double theta, long double** array, unsigned int n);


void y_plus_gen(unsigned int j, long double** array);

void y_plus_gen(long double theta, long double** array);

void y_plus_gen(unsigned int j, long double** array, unsigned int n);

void y_plus_gen(long double theta, long double** array, unsigned int n);



void y_minus_gen(unsigned int j, long double** array);

void y_minus_gen(long double theta, long double** array);

void y_minus_gen(unsigned int j, long double** array, unsigned int n);

void y_minus_gen(long double theta, long double** array, unsigned int n);
