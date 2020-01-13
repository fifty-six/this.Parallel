//
// Yusuf Bham, 6 January 2020
//
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>
#include <signal.h>

#define M 28800
#define N 21600
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

void master(int size, struct Color* palette)
{
    int*** rgb = malloc(sizeof(int**) * N);

    for (int i = 0; i < N; i++)
    {
        rgb[i] = malloc(sizeof(int*) * M);

        for (int j = 0; j < M; j++)
        {
            rgb[i][j] = malloc(sizeof(int) * 3);
        }
    }

    struct timespec start, end;
    clock_gettime(CLOCK_MONOTONIC, &start);

    MPI_Status status;

    printf("Manager started...\n");

    // Recieving dummy values so we can ensure that all nodes get started...
    int dummy;

    for (int i = 0; i < (size - 1); i++)
    {
        MPI_Recv(&dummy, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
        printf("Node %i is online\n", status.MPI_SOURCE);
    }

    size_t amount = N / (size - 1);
    fprintf(stderr, "amount: %zu, N: %i\n", amount, N);
    struct Color* res = malloc(sizeof(struct Color) * amount * M);

    // Recieve actual values...
    for (int i = 0; i < (size - 1); i++)
    {
        MPI_Recv(res, sizeof(struct Color) * amount * M, MPI_CHAR, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, &status);

        fprintf(stderr, "Recieving value...\n");

        int recv = status.MPI_SOURCE - 1;

        for (size_t i = 0; i < amount; i++)
        {
            for (size_t j = 0; j < M; j++)
            {
                struct Color c = res[j * amount + i];

                rgb[recv * amount + i][j][0] = c.r;
                rgb[recv * amount + i][j][1] = c.g;
                rgb[recv * amount + i][j][2] = c.b;
            }
        }
    }

    FILE *fout = fopen("out.ppm", "w");

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


    clock_gettime(CLOCK_MONOTONIC, &end);

    double delta_s = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) * pow(10, -9);

    printf("Time elapsed: %f\n", delta_s);
}

void slave(int size, int rank, struct Color* palette)
{
    // Send dummy value to indicate that we're actually connected.
    {
        int dummy = 0;

        MPI_Send(&dummy, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
    }

    size_t amount = N / (size - 1);

    struct Color* res = malloc(sizeof(struct Color) * amount * M);

    for (size_t i = 0; i < amount; i++)
    {
        for (size_t j = 0; j < M; j++)
        {
            res[j * amount + i] = mandelbrot(j, (rank - 1) * amount + i, palette);
        }
    }

    MPI_Send(res, sizeof(struct Color) * amount * M, MPI_CHAR, 0, 1, MPI_COMM_WORLD);

    free(res);
}

void sig_int(int sig) 
{
    printf("Caught Ctrl-C. Aborting.\n");
    MPI_Finalize();
}

int main(int argc, char* argv[])
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

    long long rseed = 1454734;

    int size;
    int rank;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    srand(rseed + rank);

    signal(SIGINT, sig_int);

    if (rank == 0)
        master(size, palette);
    else
        slave(size, rank, palette);

    MPI_Finalize();

    return 0;
}

//
// end of file
//
