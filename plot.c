#include "include.h"

void display(void) {
    static long double x[NUM_PARTICLES][6];
    static long double rS = 2.0 * G * Ms / (C * C);
    static long double a = 0.3;
    static long double disk_radius = 2.0 * G * Ms / (C * C);
    static int t = 1;
    static bool particles_initialized = false;

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(190.0, 1.0, 1.0, 380.0);
    
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    double zoom_factor = 0.4; 
    glRotatef(95.0, -35.0, 1.0, -45.0);
/*     gluLookAt(10.0 * 2 * zoom_factor, 1.0, 100.0,   
            4.0, 40.0, 0.0,             
            0.0, 0.2, -100.0); */    

    if (!particles_initialized) {
        srand(time(NULL));
        #pragma omp parallel for num_threads(8)
        for (int i = 0; i < NUM_PARTICLES; i++) {
            long double r = (long double)rand() / RAND_MAX * 100.0 + disk_radius + 800.0;
			long double theta = 2.0 * PI * ((long double)rand() / RAND_MAX); 
            long double phi = acos(1 - 2 * ((long double)rand() / RAND_MAX)); 
			long double v = sqrt(G * Ms / r);
            x[i][0] = sqrt(r * r + a * a) * sin(theta) * cos(phi); // x
            x[i][1] = sqrt(r * r + a * a) * sin(theta) * sin(phi); // y
            x[i][2] = r * cos(phi) * 2.0;
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

            x[i][3] = r * sin(theta) * cos(phi) * v + vx_rotation;
            x[i][4] = r * sin(theta) * sin(phi) * v + vy_rotation;
            x[i][5] = r * cos(theta) * v + vz_rotation;
        }
        particles_initialized = true;
    }
    glClear(GL_COLOR_BUFFER_BIT);
    //glPushMatrix();  
    for (int i = 0; i < NUM_PARTICLES; i++) {

		long double old_x = x[i][0];
		long double old_y = x[i][1];
		long double old_z = x[i][2];
		long double distance = sqrt(old_x * old_x + old_y * old_y + old_z * old_z);
		long double v = sqrt(G * Ms / distance);
		long double v_factor = 1.0 / (1.0 + distance);
		x[i][3] = (x[i][0] - old_x) * v_factor;
		x[i][4] = (x[i][1] - old_y) * v_factor;
		x[i][5] = (x[i][2] - old_z) * v_factor;

		calculate_new_position(&x[i][0], &x[i][1], &x[i][2], t * dt, rS, a, 1.6726219e-27, 938.2720813e6);
        //runge_kutta(&x[i][0], &x[i][1], &x[i][2], t * dt, rS, a);
        //calculate_gravitational_lens_effect(&x[i][0], &x[i][1], &x[i][2], rS, 500);
        //runge_kutta_fehlberg(&x[i][0], &x[i][1], &x[i][2], t * dt, rS, a); 
        long double vx = (long double)(x[i][0] - old_x) * (rand() % 100) * 0.0001;
        long double vy = (long double)(x[i][1] - old_y)  * (rand() % 100) * 0.0001; 
        long double vz = (long double)(x[i][2] - old_z) * (rand() % 100) * 0.0001;
        long double speed = sqrt(vx * vx + vy * vy + vz * vz) * 0.00000000001;

        long double redshift = sqrt(1.0f - v * v / (C * C));
        if (redshift < 0.0f) redshift = 0.0f;
        if (redshift > 1.0f) redshift = 1.0f;

        long double red = 2.0f + redshift * 0.5f; 
		long double green = 0.5f - speed * redshift; 
		long double blue = 0.0f;

		glColor3f(redshift, green, blue);
        if (redshift < 0.0f) redshift = 0.1f;
        if (redshift > 1.0f) redshift = 1.0f;
        int num_blur_frames = 5;
        //#pragma omp parallel for
        if (distance < rS) 
        {
            //afficher une sphere 
            glutSolidSphere(rS, 100, 100);
            glutCopyColormap(0);    

            //glutSolidSphere(0.1, 100, 100);
        }
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