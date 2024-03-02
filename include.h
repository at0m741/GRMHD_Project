#ifndef INCLUDE_H
#define INCLUDE_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <GL/glut.h>
#include <sys/time.h> 
#include <stdbool.h>

#define C 299792458.0  
#define G 6.67430e-11  
#define PI 3.14159265358979323846
#define Ms 42.99e22
#define dt 0.00000000000155
#define NUM_PARTICLES 100000
#define MAX_PHI_VARIATION 0.0001

typedef struct {
    long double r;
    long double theta;
    long double phi;
} State;

//long double particles[NUM_PARTICLES][6];
void init();
void runge_kutta4(long double *x, long double *y, long double *z, long double t, long double rS, long double a);
void display();
void calculate_new_position(long double *x, long double *y, long double *z, long double t, long double rS, long double a, long double mass, long double energy);
void calculate_gravitational_lens_effect(long double *x, long double *y, long double *z, long double rS, long double lens_strength);
void cleanup();
void runge_kutta_fehlberg(long double *x, long double *y, long double *z, long double t, long double rS, long double a);
//void metric_tensor(double r, double m, double gmunu[4][4]);
//void christoffel_symbols(double r, double m, double gmunu[4][4], double gammamu[4][4][4]);
void euler(long double *x, long double *y, long double *z, long double t, long double rS, long double a);
void calculateBp(double r, double theta, double phi, double Bp[3]) ;
#endif