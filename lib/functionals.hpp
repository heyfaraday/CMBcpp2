#pragma once

long double area(long double** map, long double level);

long double area(long double** map, long double level, unsigned int l_1, unsigned int l_2);

long double area(long double** map, long double level, unsigned int l_1, unsigned int l_2,
                 unsigned int phi_1, unsigned int phi_2);

long double length(long double** map, long double level);

long double length(long double** map, long double level, unsigned int l_1, unsigned int l_2);

long double length(long double** map, long double level, unsigned int l_1, unsigned int l_2,
                   unsigned int phi_1, unsigned int phi_2);

int condition_1(long double xx, long double yy, long double xy);

int condition_2(long double qx, long double qy, long double ux, long double uy);

void points_classifier(long double** map_x, long double** map_y, long double** cos_ml, long double** sin_ml,
                       std::string name);

void points_classifier(long double** map_x, long double** map_y, long double** cos_ml, long double** sin_ml,
                       std::string name, unsigned int l_1, unsigned int l_2);

void points_classifier(long double** map_x, long double** map_y, long double** cos_ml, long double** sin_ml,
                       std::string name, int** whitelist);

void level_points_classifier(long double** map_x, long double** map_y, long double** cos_ml, long double** sin_ml,
                             long double level, std::string name);

void singular_points_classifier(long double** q, long double** u, long double** q_cos_ml, long double** q_sin_ml,
                                long double** u_cos_ml, long double** u_sin_ml, std::string name);

void singular_points_classifier(long double** q, long double** u, long double** qx, long double** ux,
                                long double** qy, long double** uy, std::string name, long double** whitelist,
                                unsigned int n);

void points_classifier_p(long double** map, long double** map_x, long double** map_y, long double** map_xx,
                         long double** map_yy, long double** map_xy, std::string name, long double** whitelist,
                         unsigned int n, long double sigma_p);

void level_singular_points_classifier(long double** q, long double** u, long double** q_cos_ml, long double** q_sin_ml,
                                      long double** u_cos_ml, long double** u_sin_ml, long double level, std::string name);

void singular_points_classifier(long double** q, long double** u, long double** q_cos_ml, long double** q_sin_ml,
                                long double** u_cos_ml, long double** u_sin_ml, std::string name,
                                unsigned int l_1, unsigned int l_2);
