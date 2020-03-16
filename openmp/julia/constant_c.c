//
// Yusuf Bham, 6 January 2020
//
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAX_X 1440
#define MAX_Y 1080

#define ITER 3000

#define R_MAX +2.0
#define R_MIN -2.0
#define I_MAX +1.5
#define I_MIN -1.5

// 4 or 20 is what we were using fwiw
#define MAX_MAG 16

#define MAX_FRAMES 20

double lerp(double v0, double v1, double t) {
  return (1 - t) * v0 + t * v1;
}

struct Color
{
    int r;
    int g;
    int b;
};

struct View
{
    double r_min;
    double r_max;

    double i_min;
    double i_max;
};

struct Color mandelbrot(int x, int y, struct Color* palette, double re_c, double im_c)
{
    // z
    double a = R_MIN + x * ((R_MAX - R_MIN) / (MAX_X * 1.0));
    double b = I_MIN + y * ((I_MAX - I_MIN) / (MAX_Y * 1.0));

    double i = 0;

    for (;i < ITER; i++)
    {
        double mag = a*a + b*b;

        if (mag >= MAX_MAG)
            break;

        double a_new = a*a - b*b + re_c;

        b = 2*a*b + im_c;

        a = a_new;
    }

    if (i < ITER)
    {
        // sqrt of inner term removed using log simplification rules.
        double log_zn = log(a*a + b*b) / 2.0;
        double nu = log2(log_zn / M_LN2);
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

void gen_palette(struct Color* palette)
{
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

}

void write_file(struct Color** grid, char* fname, size_t max_y, size_t max_x)
{
    FILE *fout = fopen(fname, "w");

    fprintf(fout, "P3\n");
    fprintf(fout, "%d %d\n", (int) max_x, (int) max_y);
    fprintf(fout, "255\n");

    for (size_t y = 0; y < max_y; y++ )
    {
        for (size_t x = 0; x < max_x; x++)
        {
            struct Color c = grid[y][x];

            fprintf(fout, "%d %d %d\n", c.r, c.g, c.b);
        }
    }

    fclose(fout);
}

int main(void)
{
   struct Color palette[ITER + 1];

   gen_palette(palette);

   struct Color** grid = malloc(sizeof(struct Color*) * MAX_Y);

   for (size_t i = 0; i < MAX_Y; i++)
       grid[i] = malloc(sizeof(struct Color) * MAX_X);


   char fname[30];

   for (size_t frame = 0; frame < MAX_FRAMES; frame++)
   {
       double theta = (2 * M_PI / (double) MAX_FRAMES) * frame;

       // e^ix = rcis(theta) = r * ( cos(theta) + i sin(theta) )

       double re_z0 = cos(theta);
       double im_z0 = sin(theta);

       printf("Starting frame %03d (x = %f, y = %f).\n", (int) frame, re_z0, im_z0);

       #pragma omp parallel for
       for (size_t y = 0; y < MAX_Y; y++)
       {
           for (size_t x = 0; x < MAX_X; x++)
           {
               grid[y][x] = mandelbrot(x, y, palette, re_z0, im_z0);
           }
       }

       printf("Finished frame %03d.\n", (int) frame);

       sprintf(fname, "%03d.ppm", (int) frame);

       write_file(grid, fname, MAX_Y, MAX_X);

       printf("Wrote frame %03d.\n", (int) frame);
   }

   for (size_t i = 0; i < MAX_Y; i++)
       free(grid[i]);

   free(grid);

   return 0;
}
//
// end of file
//
