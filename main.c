#include "include.h"

void runge_kutta4(long double *x, long double *y, long double *z, long double t, long double rS, long double a) {
    long double vx, vy, vz;

    long double k1x, k1y, k1z;
    long double k2x, k2y, k2z;
    long double k3x, k3y, k3z;
    long double k4x, k4y, k4z;

    long double f1x, f1y, f1z;
    calculate_new_position(&(*x), &(*y), &(*z), t, rS, a);
    k1x = dt * f1x;
    k1y = dt * f1y;
    k1z = dt * f1z;

    calculate_new_position(&(*x), &(*y), &(*z), t, rS, a);
    k2x = dt * f1x;
    k2y = dt * f1y;
    k2z = dt * f1z;

    calculate_new_position(&(*x), &(*y), &(*z), t, rS, a);
    k3x = dt * f1x;
    k3y = dt * f1y;
    k3z = dt * f1z;

    calculate_new_position(&(*x), &(*y), &(*z), t, rS, a);
    k4x = dt * f1x;
    k4y = dt * f1y;
    k4z = dt * f1z;

    *x += (k1x + 2.0 * k2x + 2.0 * k3x + k4x) / 6.0;
    *y += (k1y + 2.0 * k2y + 2.0 * k3y + k4y) / 6.0;
    *z += (k1z + 2.0 * k2z + 2.0 * k3z + k4z) / 6.0;
}


void init() 
{
    glClearColor(0.0, 0.0, 0.0, 1.0);
}

void cleanup() 
{
    printf("Cleaning up...\n");
    glutDestroyWindow(glutGetWindow());
    exit(0);    
}

int main(int argc, char** argv) {
    glutInit(&argc, argv);
    init();
    printf("Rs = %e\n", 2.0 * G * Ms / (C * C));
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB); 
    glEnable(GL_DEPTH_TEST);
    glutInitWindowSize(920, 920);
    glutInitWindowPosition(2000, 2000);
    glutCreateWindow("Black Hole Simulation");
    glEnable(GL_BLEND);
    glEnable(GL_POINT_SMOOTH);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glutTimerFunc(1000 / 10, display, 0);
    glutDisplayFunc(display);
    glutIdleFunc(display);
    glutMainLoop();
    return 0;
} 
