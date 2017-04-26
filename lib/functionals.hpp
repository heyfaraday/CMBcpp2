#pragma once

long double area(long double** map, long double level);

long double area(long double** map, long double level, unsigned int l_1, unsigned int l_2);

long double length(long double** map, long double level);

long double length(long double** map, long double level, unsigned int l_1, unsigned int l_2);

int condition_1(long double xx, long double yy, long double xy);

int condition_2(long double qx, long double qy, long double ux, long double uy);

void points_classifier(long double** map_x, long double** map_y, long double** cos_ml, long double** sin_ml,
                       std::string name);

void points_classifier(long double** map_x, long double** map_y, long double** cos_ml, long double** sin_ml,
                       std::string name, unsigned int l_1, unsigned int l_2);

void level_points_classifier(long double** map_x, long double** map_y, long double** cos_ml, long double** sin_ml,
                             long double level, std::string name);

void singular_points_classifier(long double** q, long double** u, long double** q_cos_ml, long double** q_sin_ml,
                                long double** u_cos_ml, long double** u_sin_ml, std::string name);

void level_singular_points_classifier(long double** q, long double** u, long double** q_cos_ml, long double** q_sin_ml,
                                      long double** u_cos_ml, long double** u_sin_ml, long double level, std::string name);

void singular_points_classifier(long double** q, long double** u, long double** q_cos_ml, long double** q_sin_ml,
                                long double** u_cos_ml, long double** u_sin_ml, std::string name,
                                unsigned int l_1, unsigned int l_2);
