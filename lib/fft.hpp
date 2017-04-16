#pragma once

void fft_map_forward(long double** map, long double** cos_ml, long double** sin_ml);

long double fft_point_forward(unsigned int i, unsigned int j, long double** cos_ml, long double** sin_ml);

void fft_map_backward(long double** map, long double** cos_ml, long double** sin_ml);
