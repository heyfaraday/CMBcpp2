#pragma once

void circles_area(long double** map, long double level, std::string name);

void circles_length(long double** map, long double level, std::string name);

void circles_points(long double** map_x, long double** map_y, long double** cos_lm, long double** sin_lm,
                    std::string name);

void circles_singular_points(long double** q, long double** u, long double** q_cos_lm, long double** q_sin_lm,
                             long double** u_cos_lm, long double** u_sin_lm, std::string name);

void level_area(long double** map, std::string name);

void level_length(long double** map, std::string name);

void level_points(long double** map_x, long double** map_y, long double** cos_lm, long double** sin_lm,
                  std::string name);

void level_singular_points(long double** q, long double** u, long double** q_cos_lm, long double** q_sin_lm,
                           long double** u_cos_lm, long double** u_sin_lm, std::string name);
