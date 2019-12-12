//
// Yusuf Bham, 2 December 2019
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>
#include <signal.h>

#include <mpi.h>
#include <time.h>

#define G 6.674e-11 // (m^3)(kg^-1)(s^-2)
#define M 5.972e+24 // kg
#define R 6.371e+6  // m

#define DT 1     // s

#define M_MOON 7.349e+22 // kg
#define R_MOON 1.7374e+6 // m
#define V_MOON 1023.157  // m/s

double* create_double_array(int n)
{
    return (double*) malloc(sizeof(double) * n);
}

double simulate(
    int,
    int,
    double,
    double*,
    double*,
    double*,
    double*,
    double*,
    double*,
    double*,
    double*,
    double*,
    double*,
    bool
);

void master(int size)
{
    const int n = (int)( 0.5 + ( 160 * 60 * 60 ) / DT );

    struct timespec start, end;
    clock_gettime(CLOCK_MONOTONIC, &start);

    MPI_Status status;

    double max_res[3] = { 0, 0, 0 };

    // Theta, v0, vf.
    double res[3];

    printf("Manager started...\n");

    // Recieving dummy values so we can ensure that all nodes get started...
    int dummy;
    
    for (int i = 0; i < (size - 1); i++)
    {
        MPI_Recv(&dummy, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
        printf("Node %i is online\n", status.MPI_SOURCE);
    }

    for (int i = 0; i < (size - 1); i++)
    {
        MPI_Recv(&res, 3, MPI_DOUBLE, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);

        if (max_res[2] < res[2])
        {
            max_res[0] = res[0];
            max_res[1] = res[1];
            max_res[2] = res[2];
        }
    }

    // TODO: How are we going to send the graph, when?
    fprintf(stderr, "Max: [%f deg, %f v0, %f vmax]", max_res[0], max_res[1], max_res[2]);

    double*  t = create_double_array(n);

    double*  x = create_double_array(n);
    double*  y = create_double_array(n);
    double* vx = create_double_array(n);
    double* vy = create_double_array(n);
    double* di = create_double_array(n);

    double* x_moon = create_double_array(n);
    double* y_moon = create_double_array(n);
    double* vx_moon = create_double_array(n);
    double* vy_moon = create_double_array(n);

    simulate(max_res[1], n, max_res[0], t, x, y, vx, vy, di, x_moon, y_moon, vx_moon, vy_moon, false);

    FILE *f = fopen("orbit.txt", "w");

    for(int i = 0; i < n; i += 30)
    {
        fprintf(f, "%d %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f\n", i, t[i], x[i], y[i], di[i], vx[i], vy[i], x_moon[i], y_moon[i]);
    }

    fclose(f);

    clock_gettime(CLOCK_MONOTONIC, &end);

    double delta_s = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) * pow(10, -9);

    printf("Time elapsed: %f\n", delta_s);
}

void slave(int size, int rank)
{
    const int n = (int)( 0.5 + ( 160 * 60 * 60 ) / DT );

    const int max_speed = 7e3;
    const int start_speed = 1e3;

    // Send dummy value to indicate that we're actually connected.
    {
        int dummy = 0;

        MPI_Send(&dummy, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
    }
    // end dummy

    double res[3] = {0, 0, 0};

    double*  t = create_double_array(n);

    double*  x = create_double_array(n);
    double*  y = create_double_array(n);
    double* vx = create_double_array(n);
    double* vy = create_double_array(n);
    double* di = create_double_array(n);

    double* x_moon = create_double_array(n);
    double* y_moon = create_double_array(n);
    double* vx_moon = create_double_array(n);
    double* vy_moon = create_double_array(n);

    const int iter_amount = (max_speed - start_speed) / (size - 1);

    const int worker_start = 1000 + iter_amount * (rank - 1);
    const int worker_end = worker_start + iter_amount;

    for (double theta = 180; theta < 269; theta += 0.05)
    {
        // Ranks start at 1 and the 1st is master, so we subtract 2.
        for (int speed = worker_start; speed < worker_end; speed += 10)
        {
            double vmag = simulate(speed, n, theta, t, x, y, vx, vy, di, x_moon, y_moon, vx_moon, vy_moon, true);

            if (vmag > res[2])
            {
                res[0] = theta;
                res[1] = speed;
                res[2] = vmag;
            }
        }
    }

    MPI_Send(&res, 3, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);

    free(x);
    free(t);
    free(y);
    free(di);
    free(vx);
    free(vy);
    free(x_moon);
    free(y_moon);
}

void sig_int(int sig) 
{
    printf("Caught Ctrl-C. Aborting.\n");
    MPI_Finalize();
}

int main(int argc, char* argv[])
{
    long long rseed = 1454734;

    int size;
    int rank;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    srand(rseed + rank);

    signal(SIGINT, sig_int);

    if (rank == 0)
        master(size);
    else
        slave(size, rank);

    MPI_Finalize();

    return 0;
}

double simulate(
    int speed,
    int n,
    double theta,
    double* t,
    double* x,
    double* y,
    double* vx,
    double* vy,
    double* di,
    double* x_moon,
    double* y_moon,
    double* vx_moon,
    double* vy_moon,
    bool print
)
{
    theta *= M_PI / 180;

    double v_mag = (double) speed;

    t[0]  = 0;
    x[0]  =                3.844e8 - 1.9e6;
    y[0]  =                              0;
    di[0] =    sqrt(x[0]*x[0] + y[0]*y[0]);
    vx[0] =             v_mag * cos(theta);
    vy[0] =             v_mag * sin(theta);

    x_moon[0] =     3.844e8;
    y_moon[0] =         0.0;
    vx_moon[0] =        0.0;
    vy_moon[0] =     V_MOON;

    bool to_earth = false;
    bool to_moon = false;
    bool behind_moon = false;
    bool velocity = false;

    int j = 1;

    for(; j < n; j++)
    {
        t[j] = t[j-1] + DT;

        // EARTH => MOON
        {
            x_moon[j] = x_moon[j-1] + DT * vx_moon[j-1];
            y_moon[j] = y_moon[j-1] + DT * vy_moon[j-1];

            double distSq = x_moon[j]*x_moon[j] + y_moon[j]*y_moon[j];

            double r = sqrt(distSq);

            double a = -(G*M)/distSq;

            vx_moon[j] = vx_moon[j - 1] + a * DT * x_moon[j] / r;
            vy_moon[j] = vy_moon[j - 1] + a * DT * y_moon[j] / r;
        } // END EARTH => MOON

        // EARTH, MOON => SHIP
        {
            x[j] = x[j-1] + DT * vx[j-1];
            y[j] = y[j-1] + DT * vy[j-1];

            double distSq = x[j]*x[j] + y[j]*y[j];

            double r = sqrt(distSq);

            if (r < R)
            {
                // Skip this speed.
                return 0;
            }

            di[j] = r;

            double a = -(G*M)/distSq;

            double distSq_moon = pow(x[j] - x_moon[j], 2) + pow(y[j] - y_moon[j], 2);

            double r_moon = sqrt(distSq_moon);

            if (r_moon < 1.7374e+6)
            {
                // Skip this speed.
                return 0;
            }

            double a_moon = -(G*M_MOON)/distSq_moon;

            vx[j] = vx[j-1] + (a * DT * x[j] / r) + (a_moon * DT * (x[j] - x_moon[j]) / r_moon);
            vy[j] = vy[j-1] + (a * DT * y[j] / r) + (a_moon * DT * (y[j] - y_moon[j]) / r_moon);

            if (!to_earth && r < 1e8)
            {
                to_earth = true;
            }

            if (to_earth && !to_moon && /* t[j] > 300 && x[j] < 3e8 && */ r_moon < 3e7) // && x_moon[j] - 1e5 > x[j])
            {
                if (x_moon[j] < x[j])
                {
                    // Skip.
                    return 0;
                }

                to_moon = true;
            }

            if (to_moon && r_moon < 9e6 && x_moon[j] > x[j])
            {
                behind_moon = true;
            }

            double vmag = sqrt(vx[j]*vx[j] + vy[j]*vy[j]);

            if (behind_moon && !velocity && vmag > sqrt(2.0 * G * M / r))
            {
                velocity = true;
            }
        } // END EARTH, MOON => SHIP

    } // END TIME

    if (!velocity) return 0;

    if (print)
    {
        fprintf(stderr, "No crash! (%i, %f) Ended with velocity: (%f, %f)\n", speed, theta, vx[j-1], vy[j-1]);
    }

    return sqrt(vx[j-1] * vx[j-1] + vy[j-1] * vy[j-1]);
}
//
// end of file
//
