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

#define DT .5     // s

#define M_MOON 7.349e+22 // kg
#define R_MOON 1.7374e+6 // m
#define V_MOON 1023.157  // m/s

double* create_double_array(int n)
{
    return (double*) malloc(sizeof(double) * n);
}

int main(int argc, char** argv)
{
    if (argc == 1 || argc > 2)
    {
        fprintf(stderr, "Argument mismatch! Usage: [%s] THETA\n", argv[0]);
        return 1;
    }

    //
    // time intervals - duration is 90 minutes
    //
    int n = (int)( 0.5 + ( 300 * 60 * 60 ) / DT );

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

    t[0]  =         96302.0;

    double theta = 25.9 * (M_PI / 180);
    double dist_mag = R + 202751774.4;
    double v_mag = 1527.048;

    double scalar = atof(argv[1]);

    x[0]  =  dist_mag * cos(theta);
    y[0]  =  dist_mag * sin(theta);
    di[0] =               dist_mag;
    vx[0] =     v_mag * cos(theta);
    vy[0] =     v_mag * sin(theta);

    x_moon[0] =     3.844e8;
    y_moon[0] =         0.0;
    vx_moon[0] =        0.0;
    vy_moon[0] =     V_MOON;

    bool turned = false;

    for(int j = 1; j < n; j++)
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
        }

        // EARTH, MOON => SHIP
        {
            x[j] = x[j-1] + DT * vx[j-1];
            y[j] = y[j-1] + DT * vy[j-1];
    
            double distSq = x[j]*x[j] + y[j]*y[j];

            double r = sqrt(distSq);

            di[j] = r;
    
            double a = -(G*M)/distSq;

            double distSq_moon = pow(x[j] - x_moon[j], 2) + pow(y[j] - y_moon[j], 2);

            double r_moon = sqrt(distSq_moon);

            if (r_moon < 1.7374e+6)
            {
                fprintf(stderr, "Crashed into moon at t = %f!\n", t[j]);
                return -1;
            }

            double a_moon = -(G*M_MOON)/distSq_moon;

            vx[j] = vx[j-1] + (a * DT * x[j] / r) + (a_moon * DT * (x[j] - x_moon[j]) / r_moon);
            vy[j] = vy[j-1] + (a * DT * y[j] / r) + (a_moon * DT * (y[j] - y_moon[j]) / r_moon);

            if (!turned && di[j] < di[j - 1])
            {
                turned = true;

                vx[j] *= scalar;
                vy[j] *= scalar;
            }
        }
    }

    fout = fopen( "orbit.txt" , "w" );

    for(int j = 0; j < n; j++)
    {
        fprintf(fout, "%d %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f\n", j, t[j], x[j], y[j], di[j], vx[j], vy[j], x_moon[j], y_moon[j]);
    }

    fclose( fout );
    return 0;
}

//
// end of file
//
