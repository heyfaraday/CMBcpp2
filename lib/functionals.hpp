#pragma once

long double area(long double** map, long double level);

long double length(long double** map, long double level);

int condition_1(long double xx, long double yy, long double xy);

int condition_2(long double qx, long double qy, long double ux, long double uy);

void points_classifier(long double** map, long double** map_x, long double** map_y,
                       long double** map_xx, long double** map_yy, long double** map_xy);

void singular_points_classifier(long double** map, long double** map_x, long double** map_y,
                                long double** map_xx, long double** map_yy, long double** map_xy);
