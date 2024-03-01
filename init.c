#include "include.h"


void calculate_new_position(long double *x, long double *y, long double *z, long double t, long double rS, long double a, long double mass, long double energy)
{
    long double theta = *z;
    long double sin_theta = sin(theta);
    if (fabs(sin_theta) < 1e-10L)
        return;
    long double cos_theta = cos(theta);
    long double rho2 = *x * *x + *y * *y / (sin_theta * sin_theta);
    long double delta = rho2 - rS * rS * (1.0L + pow(a * sin_theta / rho2, 2.0L)) * 0.5f;
    long double a_rS_Ms_delta = a * (rS * (rS + 2.0L * Ms) / delta);
    long double Psi = -(a * *x * *y) / rho2 + t * a_rS_Ms_delta;
    long double sqrt_delta_sin_theta = sqrt(delta) * sin_theta;
    long double Theta = atan2(sqrt_delta_sin_theta, rS * cos_theta + a * sin_theta);
    long double Phi = atan2(*y, *x) + a_rS_Ms_delta;
    long double vx = *x * Psi + (a * *y) / sin_theta * sin(Theta);
    long double vy = *y * Psi - (a * *x) / sin_theta * -sin(Theta);
    long double vz = cos(Theta) / sin_theta * (a_rS_Ms_delta - Psi) - a * *x * *y / rho2 * sin_theta;
    long double distance = sqrt(*x * *x + *y * *y + *z * *z);
    long double rotation_factor = a / distance * 0.01; 
    long double v = sqrt(G * Ms / distance);
    long double kerrDrag = 2 * a * G * Ms / (distance * pow(C, 2.0));
    vx += kerrDrag * (*y) * rotation_factor;
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
    vx = v_r * sin_theta * cos(Phi) + v_theta * cos_theta * cos(Phi) - v_phi * sin(Phi) ;
    vy = v_r * sin_theta * sin(Phi) + v_theta * cos_theta * sin(Phi) + v_phi * cos(Phi);
    vz = v_r * cos_theta - v_theta * sin_theta;

    // Update the energy of the particle
    energy += dt + (vx * vx + vy * vy + vz * vz) / 2.0 - G * Ms / distance;
    *x = *x + dt * vx * (1.0 + energy / (mass * C * C)) / (1.0 - 2.0 * G * Ms / distance);
    *y = *y + dt * vy * (1.0 + energy / (mass * C * C)) / (1.0 - 2.0 * G * Ms / distance);
    *z = *z + dt * vz / (1.0 + energy / (mass * C * C));
}

void calculate_gravitational_lens_effect(long double *x, long double *y, long double *z, long double rS, long double lens_strength) {
    long double distance = sqrt(*x * *x + *y * *y + *z * *z);
    if (distance < rS) {
        return; 
    }
    long double deviation_factor = lens_strength * rS / distance;
    long double distance_factor = pow((rS / distance), 2.0);
    *x += deviation_factor * (*x) * distance_factor;
    *y += deviation_factor * (*y) * distance_factor;
    *z += deviation_factor * (*z) * distance_factor / 2.0;
}