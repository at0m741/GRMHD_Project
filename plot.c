#include "include.h"

#include "include.h"

void display(void) {
    static long double x[NUM_PARTICLES][6];
    static long double rS = 2.0 * G * Ms / (C * C);
    static long double a = 0.2;
    static long double disk_radius = 2.0 * G * Ms / (C * C);
    static int t = 1;
    static bool particles_initialized = false;

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(200.0, 1.0, 20.0, 1600.0);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    double zoom_factor = 0.5; 
    glRotatef(95.0, 0.0, 0.0, -45.0);
    gluLookAt(10.0 * 2 * zoom_factor, 1.0, 1000.0,   
            4.0, 40.0, 0.0,             
            0.0, 0.2, -100.0);    

    if (!particles_initialized) {
        srand(time(NULL));
        #pragma omp parallel for num_threads(8)
        for (int i = 0; i < NUM_PARTICLES; i++) {
            long double r = (long double)rand() / RAND_MAX * 100.0 + disk_radius + 600.0;
			long double theta = 2.0 * PI * ((long double)rand() / RAND_MAX); 
            long double phi = acos(1 - 2 * ((long double)rand() / RAND_MAX)); 
			long double v = sqrt(G * Ms / r);
            x[i][0] = r * sin(phi) * cos(theta) + 570; // x
            x[i][1] = r * sin(phi) * sin(theta); // y
            x[i][2] = r * cos(phi) * 2; // z
            long double angular_velocity = sqrt(G * Ms / pow(r, 3)) * 0.001; 
            long double vx = 0.0;
            long double vy = 0.0;
            long double vz = 0.0;
            vx += r * angular_velocity;
            vy += -r * angular_velocity;
            long double rotation_velocity = angular_velocity;
            long double vx_rotation = -rotation_velocity * x[i][2];
            long double vy_rotation = 1.0 * x[i][2];
            long double vz_rotation = -rotation_velocity * x[i][0];

            x[i][3] = vx;
            x[i][4] = vy;
            x[i][5] = vz;
        }
        particles_initialized = true;
    }
/*     #pragma omp parallel for
    for (int i = 0; i < 100; i++) {
        #pragma omp parallel
        for (int j = i + 1; j < 100; j++) {
            long double dx = x[i][0] - x[j][0];
            long double dy = x[i][1] - x[j][1];
            long double dz = x[i][2] - x[j][2];
            long double distance = sqrt(dx * dx + dy * dy + dz * dz);
            float COLLISION_THRESHOLD = 0.1;
            if (distance < COLLISION_THRESHOLD) {
                long double temp = x[i][3];
                x[i][3] = x[j][3];
                x[j][3] = temp;

                temp = x[i][4];
                x[i][4] = x[j][4];
                x[j][4] = temp;

                temp = x[i][5];
                x[i][5] = x[j][5];
                x[j][5] = temp;
            }
        }
    } */
    glClear(GL_COLOR_BUFFER_BIT);
    glPushMatrix();  
    for (int i = 0; i < NUM_PARTICLES; i++) {

		long double old_x = x[i][0];
		long double old_y = x[i][1];
		long double old_z = x[i][2];
		long double distance = sqrt(old_x * old_x + old_y * old_y + old_z * old_z);
		long double v = sqrt(G * Ms / distance) * 0.0000000000001;
		long double v_factor = 1.0 / (1.0 + distance) * 0.0000001;  // Fonction décroissante de la distance
		x[i][3] = (x[i][0] - old_x) * v_factor;
		x[i][4] = (x[i][1] - old_y) * v_factor;
		x[i][5] = (x[i][2] - old_z) * v_factor;

		calculate_new_position(&x[i][0], &x[i][1], &x[i][2], t * dt, rS, a, 1.6726219e-27, 938.2720813e6);
        //runge_kutta4(&x[i][0], &x[i][1], &x[i][2], t * dt, rS, a);
        calculate_gravitational_lens_effect(&x[i][0], &x[i][1], &x[i][2], rS, 0.0005);
        if (distance <= rS * 100) 
			continue;
        glColor3f(0.0, 0.0, 0.0);   
        long double vx = (long double)(x[i][0] - old_x) * (rand() % 1000) * 0.001;
        long double vy = (long double)(x[i][1] - old_y)  * (rand() % 1000) * 0.001; 
        long double vz = (long double)(x[i][2] - old_z) * (rand() % 1000) * 0.001;
        long double speed = sqrt(vx * vx + vy * vy + vz * vz) * 0.000000000001;

        long double redshift = sqrt(1.0f - v * v / (C * C));
        if (redshift < 0.0f) redshift = 0.0f;
        if (redshift > 1.0f) redshift = 1.0f;

        long double red = 2.0f + redshift * 0.5f; 
		long double green = 0.5f - speed * redshift; 
		long double blue = 0.0f;

		glColor3f(red, green, blue);
        if (redshift < 0.0f) redshift = 0.1f;
        if (redshift > 1.0f) redshift = 1.0f;
        int num_blur_frames = 5;
        //#pragma omp parallel for
        for (int j = 0; j < num_blur_frames; j++) 
        {
            long double t = (long double)j * (num_blur_frames - 1);
            long double fx = (long double)old_x + t * vx * 0.001;
            long double fy = (long double)old_y + t * vy * 0.001;
            long double fz = (long double)old_z + t * vz * 0.001;
            long double distance_to_center = sqrt(fx * fx + fy * fy + fz * fz);
            if (distance_to_center < rS)
            {
                fx = 0.0;
                fy = 0.0;
                fz = 0.0;
                vx = 0.0;
                vy = 0.0;
                vz = 0.0;
                Ms * 1.6726219e-27;
                rS *= 2.0;
            }
            if (distance > rS * 10)
            {
                long double slowdown_factor = 1e-12;
                vx *= slowdown_factor;
                vy *= slowdown_factor;
                vz *= slowdown_factor;
            }
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