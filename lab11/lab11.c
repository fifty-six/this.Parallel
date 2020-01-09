//
// Yusuf Bham, 6 January 2020
//
#include <stdio.h>
#include <math.h>

#define M 4800
#define N 3600
#define ITER 1000

#define R_MAX 2.0
#define R_MIN -R_MAX
#define I_MAX 1.5
#define I_MIN -1.5

double lerp(double v0, double v1, double t) {
  return (1 - t) * v0 + t * v1;
}

struct Color
{
    int r;
    int g;
    int b;
};

struct Color mandelbrot(int x, int y, struct Color* palette)
{
    double re = R_MIN + x * ((R_MAX - R_MIN) / (M * 1.0));
    double im = I_MIN + y * ((I_MAX - I_MIN) / (N * 1.0));

    double a = 0;
    double b = 0;

    double i = 0;

    for (;i < ITER; i++)
    {
        double mag = a*a + b*b;

        if (mag >= 20)
            break;

        double a_new = a*a - b*b + re;

        b = 2*a*b + im;

        a = a_new;
    }

    if (i < ITER)
    {
        // sqrt of inner term removed using log simplification rules.
        double log_zn = log(a*a + b*b) / 2.0;
        double nu = log(log_zn / log(2.0)) / log(2.0);
        // Rearranging the potential function.
        // Dividing log_zn by log(2) instead of log(N = 1<<8)
        // because we want the entire palette to range from the
        // center to radius 2, NOT our bailout radius.
        i = i + 1 - nu;
    }

    struct Color c1 = palette[(int) i];
    struct Color c2 = palette[( (int) i + 1 ) > ITER ? (int) i : (int) i + 1];

    double frac = i - ((int) i);

    struct Color ret = (struct Color) {
        .r = (int) lerp(c1.r, c2.r, frac),
        .g = (int) lerp(c1.g, c2.g, frac),
        .b = (int) lerp(c1.b, c2.b, frac)
    };

    return ret;
}

int rgb[N][M][3]; // red-green-blue for each pixel

int main(void)
{
   struct Color palette[ITER + 1];

   for (size_t i = 0; i < ITER + 1; i++)
   {
       if (i >= ITER)
       {
           palette[i] = (struct Color) {
               .r = 0,
               .g = 0,
               .b = 0
           };

           continue;
       }

       double c = 3.0 * (i == 0 ? 0 : log(i)) / log(ITER - 1.0);

       if (c < 1)
       {
           palette[i] = (struct Color) { 
               .r = 0,
               .g = 0,
               .b = 255 * c
           };
       }
       else if (c < 2) 
       {
           palette[i] = (struct Color) { 
               .r = 0,
               .g = 255 * (c - 1),
               .b = 255
           };
       }
       else
       {
           palette[i] = (struct Color) { 
               .r = 255 * (c - 2),
               .g = 255,
               .b = 255
           };
       }
   }

   FILE* fout;

   for(size_t y = 0; y < N; y++ )
   {
      for(size_t x = 0; x < M; x++)
      {
         struct Color c = mandelbrot(x, y, palette);

         rgb[y][x][0] = c.r;
         rgb[y][x][1] = c.g;
         rgb[y][x][2] = c.b;
      }
   }

   fout = fopen("out.ppm", "w");

   fprintf(fout, "P3\n");
   fprintf(fout, "%d %d\n", M, N);
   fprintf(fout, "255\n");

   for (size_t y = 0; y < N; y++ )
   {
       for (size_t x = 0; x < M; x++)
       {
           fprintf(fout, "%d %d %d\n", rgb[y][x][0], rgb[y][x][1], rgb[y][x][2]);
       }
   }

   fclose(fout);

   return 0;
}
//
// end of file
//
