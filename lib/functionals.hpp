#pragma once

long double area(long double** map, long double level);

long double length(long double** map, long double level);

int condition_1(long double xx, long double yy, long double xy);

int condition_2(long double qx, long double qy, long double ux, long double uy);

void points_classifier(long double** map_x, long double** map_y, long double** cos_lm, long double** sin_lm,
                       std::string name);

void singular_points_classifier(long double** q, long double** u, long double** q_cos_lm, long double** q_sin_lm,
                                long double** u_cos_lm, long double** u_sin_lm, std::string name);
