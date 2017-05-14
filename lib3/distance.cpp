#include <cmath>

#include "distance.hpp"

// function computes length from point (phi1, theta1) to (phi2, theta2) on sphere
// returns value from [0, pi]
// takes phi1 and phi2 from [0, 2 * pi]
// takes theta1 and theta2 from [0, pi]
long double s2(long double phi1, long double phi2, long double theta1, long double theta2) {
    return acosl(cosl(theta1) * cosl(theta2) + sinl(theta1) * sinl(theta2) * cosl(phi1 - phi2));
}

// function compute length from point (x1, y1) to (x2, y2) on euclidean space
// return value in R
// takes x1, x2, y1 and y2 in R
long double r2(long double x1, long double x2, long double y1, long double y2) {
    return sqrtl(x1 * x2 + y1 * y2);
}

// function compute coordinates (phi_out, theta_out) of the point of intersection of two great circles
// the first circle passes through (phi1a, theta1a) and (phi1b, theta1b)
// and the second circle through (phi2a, theta2a) and (phi2b, theta2b)
// return phi_out in [0, 2 * pi], theta_out in [0, pi]
void cross(long double& phi_out, long double& theta_out,
           long double phi1a, long double theta1a, long double phi1b, long double theta1b,
           long double phi2a, long double theta2a, long double phi2b, long double theta2b) {

    long double a1[3] = {sinl(theta1a) * cosl(phi1a), sinl(theta1a) * sinl(phi1a), - cosl(theta1a)};
    long double b1[3] = {sinl(theta1b) * cosl(phi1b), sinl(theta1b) * sinl(phi1b), - cosl(theta1b)};
    long double c1[3] = {a1[1] * b1[2] - a1[2] * b1[1], a1[2] * b1[0] - a1[0] * b1[2], a1[0] * b1[1] - a1[1] * b1[0]};

    long double a2[3] = {sinl(theta2a) * cosl(phi2a), sinl(theta2a) * sinl(phi2a), - cosl(theta2a)};
    long double b2[3] = {sinl(theta2b) * cosl(phi2b), sinl(theta2b) * sinl(phi2b), - cosl(theta2b)};
    long double c2[3] = {a2[1] * b2[2] - a2[2] * b2[1], a2[2] * b2[0] - a2[0] * b2[2], a2[0] * b2[1] - a2[1] * b2[0]};

    long double a[3] = {c1[1] * c2[2] - c1[2] * c2[1], c1[2] * c2[0] - c1[0] * c2[2], c1[0] * c2[1] - c1[1] * c2[0]};
    long double mod = sqrtl(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);

    phi_out = atanl(a[1] / a[0]);
    theta_out = acosl(- a[2] / mod);
}

