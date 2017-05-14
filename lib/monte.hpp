#pragma once

void circles_area(long double** map, long double level, std::string name);

void circles_area_sigma(long double** map, std::string name, long double sigma);

void circles_area_map(long double** map, long double level, long double** out_map);

void circles_length(long double** map, long double level, std::string name);

void circles_length_sigma(long double** map, std::string name, long double sigma);

void circles_length_map(long double** map, long double level, long double** out_map);

void circles_points(long double** map_x, long double** map_y, long double** cos_ml, long double** sin_ml,
                    std::string name);

void circles_singular_points(long double** q, long double** u, long double** q_cos_ml, long double** q_sin_ml,
                             long double** u_cos_ml, long double** u_sin_ml, std::string name);

void level_area(long double** map, std::string name);

void level_length(long double** map, std::string name);

void level_points(long double** map_x, long double** map_y, long double** cos_ml, long double** sin_ml,
                  std::string name);

void level_singular_points(long double** q, long double** u, long double** q_cos_ml, long double** q_sin_ml,
                           long double** u_cos_ml, long double** u_sin_ml, std::string name);

void minkovscki_on_sphere(long double** map, long double** map_out, long double level);
