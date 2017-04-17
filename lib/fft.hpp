#pragma once

void fft_map_forward(long double** map, long double** cos_ml, long double** sin_ml);

long double fft_point_forward(unsigned int i, unsigned int j, long double** cos_ml, long double** sin_ml);

void fft_map_backward(long double** map, long double** cos_ml, long double** sin_ml);

void fft_map_x_forward(long double** map, long double** cos_ml, long double** sin_ml);

long double fft_point_x_forward(unsigned int i, unsigned int j, long double** cos_ml, long double** sin_ml);

void fft_map_y_forward(long double** map, long double** cos_ml, long double** sin_ml);

long double fft_point_y_forward(unsigned int i, unsigned int j, long double** cos_ml, long double** sin_ml);

void fft_map_xx_forward(long double** map, long double** cos_ml, long double** sin_ml);

long double fft_point_xx_forward(unsigned int i, unsigned int j, long double** cos_ml, long double** sin_ml);

void fft_map_yy_forward(long double** map, long double** cos_ml, long double** sin_ml);

long double fft_point_yy_forward(unsigned int i, unsigned int j, long double** cos_ml, long double** sin_ml);

void fft_map_xy_forward(long double** map, long double** cos_ml, long double** sin_ml);

long double fft_point_xy_forward(unsigned int i, unsigned int j, long double** cos_ml, long double** sin_ml);

void fft_map_spin1_backward(long double** map, long double** cos_ml, long double** sin_ml);

void fft_map_spin2_backward(long double** map, long double** cos_ml, long double** sin_ml);
