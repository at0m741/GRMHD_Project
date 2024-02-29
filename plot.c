#include "include.h"

void display() 
{
    static long double x[NUM_PARTICLES][3];
    static long double rS = 2.0 * G * Ms / (C * C);
    static long double a = 0.8;
    static long double disk_radius = 2.0 * G * Ms / (C * C) ;
    static int t = 1;
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(200.0, 1.6, 20.0, 1000.0);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    double zoom_factor = 0.5; 
    gluLookAt(40.0 * 2 * zoom_factor, 1.0, 600.0 * zoom_factor,   
            5.0, 40.0, 0.0,             
            0.0, 0.2, -100.0);    

    if (t == 1) {
        srand(time(NULL));
        long double r = (long double)rand() / RAND_MAX * 100.0 + disk_radius + 600.0; 
        #pragma omp parallel for

        for (int i = 0; i < NUM_PARTICLES; i++) {
            long double r = (long double)rand() / RAND_MAX * 100.0 + disk_radius + 600.0; 
            long double theta = 2.0 * PI * ((long double)rand() / RAND_MAX); 
            long double phi = acos(1 - 2 * ((long double)rand() / RAND_MAX)); 
            x[i][0] = r * sin(phi) * cos(theta) + 500; // x
            x[i][1] = r * sin(phi) * sin(theta); // y
            x[i][2] = r * cos(phi) * 2; // z
            long double v = sqrt(G * Ms / r);
            long double angular_velocity = sqrt(G * Ms / pow(r, 3)) * 0.0000000001; 
            long double vx = -angular_velocity * x[i][1];
            long double vy = angular_velocity * x[i][0];
            long double vz = angular_velocity * x[i][2];
            vx += (long double)(rand() % 1000) / 1000.0 * v * 0.001;
            vy += (long double)(rand() % 1000) / 1000.0 * v * 0.001;
            long double rotation_velocity = angular_velocity * 0.0000000001;
            long double vx_rotation = -rotation_velocity * x[i][2];
            long double vy_rotation = 0.0;
            long double vz_rotation = rotation_velocity * x[i][0];

            /* x[i][3] = vx + vx_rotation;
            x[i][4] = vy + vy_rotation;
            x[i][5] = vz + vz_rotation; */

        }
    }
    glClear(GL_COLOR_BUFFER_BIT);
    glPushMatrix(); 
    glTranslatef(0.0, 0.0, 0.0); 
    glColor3f(1.0, 1.0, 1.0); 
    glutSolidSphere(disk_radius * 100, 50, 50);
    glPopMatrix();    
    #pragma omp parallel for
    for (int i = 0; i < NUM_PARTICLES; i++) {
        long double old_x = x[i][0];
        long double old_y = x[i][1];
        long double old_z = x[i][2];
        long double distance = sqrt(old_x * old_x + old_y * old_y + old_z * old_z);

        //calculate_new_position(&x[i][0], &x[i][1], &x[i][2], t * dt, rS, a);
        runge_kutta4(&x[i][0], &x[i][1], &x[i][2], t * dt, rS, a);
        calculate_gravitational_lens_effect(&x[i][0], &x[i][1], &x[i][2], rS , 1000);

        if (distance <= rS * 100) continue;

        glColor3f(1.0, 0.0, 0.0);   
        float vx = (float)(x[i][0] - old_x) * (rand() % 1000);
        float vy = (float)(x[i][1] - old_y)  * (rand() % 1000); 
        float vz = (float)(x[i][2] - old_z) * (rand() % 1000);
        float speed = sqrt(vx * vx + vy * vy + vz * vz);

        float redshift = 0.01f - distance;
        if (redshift < 0.0f) redshift = 0.0f;
        if (redshift > 1.0f) redshift = 1.0f;

        float red = 2.0f + redshift * 0.5f; 
        float green = 0.5f + 0.5f * redshift; 
        float blue = 0.0f;

        glColor3f(red, green, blue);
        if (redshift < 0.0f) redshift = 0.0f;
        if (redshift > 1.0f) redshift = 1.0f;
        int num_blur_frames = 5;
        #pragma omp parallel for
        for (int j = 0; j < num_blur_frames; j++) {
            float t = (float)j / (num_blur_frames - 1);
            float fx = (float)old_x + t * vx;
            float fy = (float)old_y + t * vy;
            float fz = (float)old_z + t * vz;
            float distance_to_center = sqrt(fx * fx + fy * fy + fz * fz);

            if (distance_to_center < distance_to_center) continue;

            if (distance <= disk_radius)
            {
               continue;    
            }
            float red = 2.0f + redshift; 
            float green = 0.5f + 0.5f * redshift; 
            float blue = 0.0f;

            glColor3f(red, green, blue);
            glBegin(GL_POINTS);
            glVertex3f(fx, fy, fz);
            glEnd();
        }
    }

    t++;
    glFlush();
    glutPostRedisplay();
    glutSwapBuffers();
}