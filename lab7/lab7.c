//
// Yusuf Bham, 17 October 2019
//

#include <stdio.h>
#include <math.h>

#define G 6.674e-11 // (m^3)(kg^-1)(s^-2)
#define M 5.972e+24 // kg
#define R 6.371e+6  // m

#define DT 0.25     // s

int main()
{
    //
    // time intervals - duration is 90 minutes
    //
    int n = (int)( 0.5 + ( 5 * 60 * 60 ) / DT );

    double  t[n];
    double  x[n];
    double  y[n];
    double vx[n];
    double vy[n];
    double di[n];

    FILE*  fout;

    // International Space Station
    //
    // https://www.nasa.gov/sites/default/files/atoms/files/np-2015-05-022-jsc-iss-guide-2015-update-111015-508c.pdf
    //
    // Page 54 - altitude : 370 km to 460 km
    // Page 54 - speed    : 28,000 km per hour

    t[0]  =          0.0;
    x[0]  =          0.0;
    y[0]  = R + 400000.0;
    y[0]  *=         1.5;
    vx[0] =       7670.1;
    vx[0] *=         1.5;
    vy[0] =          0.0;

    for(int j = 1; j < n; j++)
    {
        t[j] = t[j-1] + DT;

        x[j] = x[j-1] + DT * vx[j-1];
        y[j] = y[j-1] + DT * vy[j-1];

        double distSq = x[j]*x[j] + y[j]*y[j];

        double r = sqrt(distSq);

        di[j] = r;

        double a = -(G*M)/distSq;

        vx[j] = vx[j-1] + a * DT * x[j] / r;
        vy[j] = vy[j-1] + a * DT * y[j] / r;
    }

    fout = fopen( "orbit.txt" , "w" );

    for(int j = 0; j < n; j++)
    {
        fprintf(fout, "%d %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f\n", j, t[j], x[j], y[j], di[j], vx[j], vy[j]);
    }

    fclose( fout );
    return 0;
}

//
// end of file
//
