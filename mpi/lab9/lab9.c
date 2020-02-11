//
// Yusuf Bham, 17 October 2019
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>

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

int main(int argc, char** argv)
{
    // time 
    int n = (int)( 0.5 + ( 160 * 60 * 60 ) / DT );

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

    FILE*  fout;

    // International Space Station
    //
    // https://www.nasa.gov/sites/default/files/atoms/files/np-2015-05-022-jsc-iss-guide-2015-update-111015-508c.pdf
    //
    // Page 54 - altitude : 370 km to 460 km
    // Page 54 - speed    : 28,000 km per hour

    bool success = false;

    for (double theta_deg = 190; theta_deg <= 205; theta_deg += .5)
    {
        double theta = theta_deg * (M_PI / 180);

        for (int speed = 1000; speed < 5e3; speed += 150)
        {
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

            bool skip = false;

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
                        skip = true;
                        break;
                    }

                    di[j] = r;

                    double a = -(G*M)/distSq;

                    double distSq_moon = pow(x[j] - x_moon[j], 2) + pow(y[j] - y_moon[j], 2);

                    double r_moon = sqrt(distSq_moon);

                    if (r_moon < 1.7374e+6)
                    {
                        // Skip this speed.
                        skip = true;
                        break;
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
                            skip = true;
                            break;
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

            if (skip) continue;

            if (velocity)
            {
                fprintf(stderr, "No crash! (%i, %f) Ended with velocity: (%f, %f)\n", speed, theta_deg, vx[j], vy[j]);
                success = true;
                goto write;
            }

        } // END SPEED

    } // END THETA

    if (!success)
    {
        fprintf(stderr, "Unable to find appropriate theta and vmag :(\n");
        goto free_mem;
    }

write:

    fout = fopen( "orbit.txt" , "w" );

    fprintf(stderr, "writing!");

    for(int i = 0; i < n; i += 30)
    {
        fprintf(fout, "%d %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f\n", i, t[i], x[i], y[i], di[i], vx[i], vy[i], x_moon[i], y_moon[i]);
    }

    fclose( fout );

free_mem:

    free(x);
    free(t);
    free(y);
    free(di);
    free(vx);
    free(vy);
    free(x_moon);
    free(y_moon);

    return success ? 0 : -1;
}

//
// end of file
//
