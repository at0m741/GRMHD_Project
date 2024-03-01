
void calculate_new_position(long double *x, long double *y, long double *z, long double t, long double rS, long double a) {
    long double x2 = *x * *x;
    long double y2 = *y * *y;
    long double z2 = *z * *z;
    long double theta = *z;
    long double sin_theta = sin(theta);
    if (fabs(sin_theta) < 1e-10L) {
        return;
    }
    long double cos_theta = cos(theta);
    long double sin_theta2 = sin_theta * sin_theta;
    long double rho2 = x2 + y2 / sin_theta2;
    long double a_sin_theta = a * sin_theta;
    long double a_sin_theta_div_rho2 = a_sin_theta / rho2;
    long double delta = rho2 - rS * rS * (1.0L + a_sin_theta_div_rho2 * a_sin_theta_div_rho2) / 2;
    long double a_rS_Ms_delta = a * (rS * (rS + 2.0L * Ms) / delta);
    long double Psi = -(a * *x * *y) / rho2 + t * a_rS_Ms_delta;
    long double sqrt_delta_sin_theta = sqrt(delta) * sin_theta;
    long double Theta = atan2(sqrt_delta_sin_theta, rS * cos_theta + a_sin_theta);
    long double Phi = atan2(*y, *x) + a_rS_Ms_delta;
    double vx = *x * Psi + a * *y * sin(Theta) / sin_theta;
    double vy = *y * Psi - a * *x * sin(Theta) / sin_theta;
    double vz = cos(Theta) / sin_theta * (a_rS_Ms_delta - Psi) - a * *x * *y / rho2;

    long double distance = sqrt(x2 + y2 + z2);
    long double rotation_factor = a / distance * 0.01; 
    double v = sqrt(G * Ms / distance);
    double kerrDrag = 2 * a * G * Ms / (distance * pow(C, 2.0));
    vx += kerrDrag * (*y) * rotation_factor;
    vy -= kerrDrag * (*x) * rotation_factor;

    *x += dt * vx;
    *y += dt * vy;
    *z += dt * vz;
}