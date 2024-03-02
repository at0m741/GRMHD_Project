#include "include.h"

long double surface_density(long double r) {
    long double rho0 = 1.0; 
    long double p = -3.0 / 2.0; 
    long double rho = rho0 * pow(r, p);

    return rho;
}

void metric_tensor_kerr(long double r, long double a, long double gmunu[4][4]) {
    long double delta = r * r - 2 * r + a * a;
    long double sigma = r * r + a * a * cos(PI / 2.0);
    long double rho2 = r * r + a * a * cos(PI / 2.0);
    long double rho = sqrt(rho2);
    long double gtt = -1 + 2 * r / rho2;
    long double grr = rho2 / delta;
    long double gthth = rho2;
    long double gphph = (r * r + a * a + 2 * a * a * r * sin(PI / 2.0) / rho2) * sin(PI / 2.0) * sin(PI / 2.0);
    long double gph = -2 * a * r * sin(PI / 2.0) / rho2;
    gmunu[0][0] = gtt;
    gmunu[1][1] = grr;
    gmunu[2][2] = gthth;
    gmunu[3][3] = gphph;
    gmunu[3][0] = gmunu[0][3] = gph;
}

void christoffel_symbols(long double r, long double a, long double gmunu[4][4], long double gammamu[4][4][4]) {
    long double delta = r * r - 2 * r + a * a;
    long double sigma = r * r + a * a * cos(PI / 2.0);
    long double rho2 = r * r + a * a * cos(PI / 2.0);
    long double rho = sqrt(rho2);
    long double gtt = -1 + 2 * r / rho2;
    long double grr = rho2 / delta;
    long double gthth = rho2;
    long double gphph = (r * r + a * a + 2 * a * a * r * sin(PI / 2.0) / rho2) * sin(PI / 2.0) * sin(PI / 2.0);
    long double gph = -2 * a * r * sin(PI / 2.0) / rho2;
    long double dgtt_dr = -2 * (r - a * a) / (rho2 * rho2);
    long double dgtt_da = 2 * a * (r - a * a) * sin(PI / 2.0) / (rho2 * rho2);
    long double dgrr_dr = (2 * r * (r - a * a) * (r - a * a) - 2 * r * rho2 * rho2) / (delta * delta * rho2);
    long double dgrr_da = -2 * a * (r * r - a * a * a * a) * sin(PI / 2.0) / (delta * rho2 * rho2);
    long double dgphph_dr = 2 * a * a * (r * r - a * a * a * a) * sin(PI / 2.0) / (rho2 * rho2);
    long double dgphph_da = -2 * a * a * a * (r - a * a) * sin(PI / 2.0) / (rho2 * rho2);
    long double dgrt_dr = -r / (rho2 * rho2);
    long double dgrt_da = a * sin(PI / 2.0) / (rho2 * rho2);
    long double dgthph_dr = -a / (rho2 * rho2);
    long double dgthph_da = -r * a * sin(PI / 2.0) / (rho2 * rho2);
    long double dgrph_dr = -a / (rho2 * rho2);
    long double dgrph_da = -r * sin(PI / 2.0) / (rho2 * rho2);
    long double dgtht_dr = 0;
    long double dgtht_da = 0;

    #pragma omp parallel for
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            for (int k = 0; k < 4; k++) {
                gammamu[i][j][k] = 0.5 * (gmunu[i][j] * (dgtt_dr + dgtt_da) - gmunu[i][k] * (dgtt_dr + dgtt_da) + gmunu[j][k] * (dgtt_dr + dgtt_da));
            }
        }
    }
}

void runge_kutta4(long double *x, long double *y, long double *z, long double t, long double rS, long double a) {
    long double vx, vy, vz;
    long double proton_mass = 1.6726219e-27;
    long double electron_mass = 9.10938356e-31;
    long double proton_energy = 938.2720813e6;
    long double k1x, k1y, k1z;
    long double k2x, k2y, k2z;
    long double k3x, k3y, k3z;
    long double k4x, k4y, k4z;
    long double f1x, f1y, f1z;
    long double gmunu[4][4];
    long double gammamu[4][4][4];

    christoffel_symbols(sqrt(*x * *x + *y * *y + *z * *z), a, gmunu, gammamu);
    calculate_new_position(&(*x), &(*y), &(*z), t, rS, a, proton_mass, proton_energy);
    k1x = dt * f1x;
    k1y = dt * f1y;
    k1z = dt * f1z;

    calculate_new_position(&(*x), &(*y), &(*z), t, rS, a, proton_mass, proton_energy);
    k2x = dt * f1x;
    k2y = dt * f1y;
    k2z = dt * f1z;

    calculate_new_position(&(*x), &(*y), &(*z), t, rS, a, proton_mass, proton_energy);
    k3x = dt * f1x;
    k3y = dt * f1y;
    k3z = dt * f1z;

    calculate_new_position(&(*x), &(*y), &(*z), t, rS, a, proton_mass, proton_energy);
    k4x = dt * f1x;
    k4y = dt * f1y;
    k4z = dt * f1z;

    *x += (k1x + 2.0 * k2x + 2.0 * k3x + k4x) / 6.0;
    *y += (k1y + 2.0 * k2y + 2.0 * k3y + k4y) / 6.0;
    *z += (k1z + 2.0 * k2z + 2.0 * k3z + k4z) / 6.0;

}

void calculate_new_position(long double *x, long double *y, long double *z, long double t, long double rS, long double a, long double mass, long double energy)
{
    long double theta = nextafter(PI / 2.0L, 0.0L);
    long double sin_theta = sin(theta);
    if (fabs(sin_theta) < 1e-10L)
        return;
    long double cos_theta = cos(theta);
    long double rho2 = rS * rS + a * a * cos(PI / 2.0);
    long double delta = rS * rS * (1.0L + pow(a * sin_theta / rho2, 2.0L)) * 0.5f;
    long double a_rS_Ms_delta = a * (rS * (rS + 2.0L * Ms) / delta);
    long double Psi = -(a * *x * *y) / rho2 + t * a_rS_Ms_delta;
    long double sqrt_delta_sin_theta = sqrt(delta) * sin_theta;
    long double Theta = atan2(sqrt_delta_sin_theta, rS * cos_theta + a * sin_theta);
    long double Phi = atan2(*y, *x) + a_rS_Ms_delta ;
    long double vx = *x * Psi + (a * *y) / sin_theta * sin(Theta);
    long double vy = *y * Psi - (a * *x) / sin_theta * -sin(Theta);
    long double vz = cos(Theta) / sin_theta * (a_rS_Ms_delta - Psi) - a * *x * *y / rho2 * sin_theta;
    long double distance = sqrt(*x * *x + *y * *y + *z * *z);
    long double rotation_factor = a / distance * 0.01; 
    long double v = sqrt(G * Ms / distance);
    long double kerrDrag = 2 * a * G * Ms / (distance * pow(C, 2.0));
    vx += kerrDrag * (*y) * rotation_factor * 0.000005;
    vy -= kerrDrag * (*x) * rotation_factor;

    double r = distance;
    double rs_a = rS + a;
    double rs_a_square = rs_a * rs_a;
    double omega = (2 * G * Ms * a) / (r * r * r + a * a * r + a * a * a * rs_a);
    double kappa = (2 * G * Ms * r) / (r * r * r + a * a * r + a * a * a * rs_a);
    double v_phi = r * sin_theta * omega;
    double v_t = sqrt(fabs((omega - kappa) * (omega + kappa)));
    double v_r = sqrt(fabs((omega * omega * r * r + kappa * kappa - 1) / (1 - r * r * (omega * omega - 1))));
    double v_theta = sqrt(fabs((kappa * kappa - cos_theta * cos_theta) / (1 - a * a * cos_theta * cos_theta / (r * r))));
    
    // Update the velocity of the particle with relativistic correction
    vx = v_r * sin_theta * cos(Phi) + v_theta * cos_theta * cos(Phi) - v_phi * sin(Phi)* 0.000006 ;
    vy = v_r * sin_theta * sin(Phi) + v_theta * cos_theta * sin(Phi) + v_phi * cos(Phi) * 0.000005;
    vz = v_r * cos_theta - v_theta * sin_theta * 0.000005;
    long double rho = surface_density(r);
    mass += 4 * PI * G * ((rand() % 10)) * rho / (sqrt(3) * pow(C, 3));
    // Update the energy of the particle
    energy += dt + (vx * vx + vy * vy + vz * vz) / 2.0 - G * Ms / distance;
    *x = *x + dt * vx * (1.0 + energy / (mass * C * C)) / (1.0 - 2.0 * G * Ms / distance);
    *y = *y + dt * vy * (1.0 + energy / (mass * C * C)) / (1.0 - 2.0 * G * Ms / distance);
    *z = *z + dt * vz / (1.0 + energy / (mass * C * C));
    
}

/* void calculate_new_position(long double *x, long double *y, long double *z, long double t, long double rS, long double a, long double mass, long double energy)
{
    long double gmunu[4][4];
    long double gammamu[4][4][4];
    christoffel_symbols(sqrt(*x * *x + *y * *y + *z * *z), a, gmunu, gammamu);

    long double dvx = 0.0, dvy = 0.0, dvz = 0.0;
    #pragma omp parallel
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            dvx -= gammamu[0][i][j] * gmunu[i][j];
            dvy -= gammamu[1][i][j] * gmunu[i][j];
            dvz -= gammamu[2][i][j] * gmunu[i][j];
            printf("dvx = %Le, dvy = %Le, dvz = %Le\n", dvx, dvy, dvz);
        }
    }

    // Mettre à jour la vitesse
    *x += dt * dvx;
    *y += dt * dvy;
    *z += dt * dvz;

    // Mettre à jour la position
    *x += dt * *x;
    *y += dt * *y;
    *z += dt * *z;
} */

void calculate_gravitational_lens_effect(long double *x, long double *y, long double *z, long double rS, long double lens_strength) 
{
    long double distance = sqrt(*x * *x + *y * *y + *z * *z);
    if (distance < rS) {
        return; 
    }

    long double xi_E = sqrt((2 * G * Ms / (C * C)) * (distance - rS) / (distance * rS)) / distance / C / C / 2.0;

    long double deviation_factor = lens_strength * rS / distance;
    long double distance_factor = pow((rS / distance), 4.0);
    *x -= deviation_factor * (*x) * distance_factor;
    *y += deviation_factor * (*y) * distance_factor;
    *z += deviation_factor * (*z) * distance_factor / 2.0;

    *x += xi_E * (*x);
    *y += xi_E * (*y);
    *z += xi_E * (*z);
}
