#ifndef INCLUDE_H
#define INCLUDE_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <GL/glut.h>
#include <sys/time.h> 

#define C 299792458.0  
#define G 6.67430e-11  
#define PI 3.14159265358979323846
#define Ms 6e15
#define dt 0.0002     
#define NUM_PARTICLES 200000

//long double particles[NUM_PARTICLES][6];
void init();
void runge_kutta4(long double *x, long double *y, long double *z, long double t, long double rS, long double a);
void display();
void calculate_new_position(long double *x, long double *y, long double *z, long double t, long double rS, long double a);
void calculate_gravitational_lens_effect(long double *x, long double *y, long double *z, long double rS, long double lens_strength);
void cleanup();
void metric_tensor(double r, double m, double gmunu[4][4]);
void christoffel_symbols(double r, double m, double gmunu[4][4], double gammamu[4][4][4]);

#endif