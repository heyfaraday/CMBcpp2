#pragma once


long double s2(long double phi1, long double phi2, long double theta1, long double theta2);

long double r2(long double x1, long double x2, long double y1, long double y2);

void cross(long double& phi_out, long double& theta_out,
           long double phi1a, long double theta1a, long double phi1b, long double theta1b,
           long double phi2a, long double theta2a, long double phi2b, long double theta2b);
