//
// Yusuf Bham, January 2020
//

#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>

#define MAX_X 1920
#define MAX_Y 1080
#define EPSILON 0.001

typedef struct
{
   double x ;
   double y ;
   double z ;

} Vector3;

typedef struct
{
    int r;
    int g;
    int b;

} Color;

typedef struct
{
    // Radius
    double r;

    // Center
    Vector3 c;

    // Color
    Color h;

} Sphere;

#define SPHERE_NUM 4

Sphere spheres[SPHERE_NUM];

Vector3 eye = { .x = 0.50, .y = 0.50, .z = -1.00 };
Vector3 light = { .x = 0.00, .y = 1.25, .z = -0.50 }; 

double dot_vec( Vector3 t , Vector3 u )
{
   return t.x * u.x + t.y * u.y + t.z * u.z ;
}

inline double min(double a, double b)
{
    return a < b ? a : b;
}

inline double max(double a, double b)
{
    return a > b ? a : b;
}

void init_objects()
{
   // Floor
   spheres[0].c = (Vector3) { 
       .x = 0.50,
       .y = -20000.00,
       .z = 0.50 
   };
   spheres[0].r = 20000.25;
   spheres[0].h = (Color) { 
       .r = 205,
       .g = 133,
       .b = 63 
   };

   // Blue sphere
   spheres[1].c = (Vector3) {
       .x = 0.50,
       .y = 0.50,
       .z = 0.50
   };
   spheres[1].r  = 0.25;
   spheres[1].h = (Color) {
       .r = 0,
       .g = 0,
       .b = 255 
   };

   // Green sphere 
   spheres[2].c = (Vector3) {
       .x = 1.00,
       .y = 0.50,
       .z = 1.00 
   };
   spheres[2].r = 0.25;
   spheres[2].h = (Color) {
       .r = 0,
       .g = 255,
       .b = 0 
   };

   // Red sphere
   spheres[3].c = (Vector3) {
       .x = 0.00,
       .y = 0.75,
       .z = 1.25 
   };
   spheres[3].r = 0.50;
   spheres[3].h = (Color) { 
       .r = 255,
       .g = 0,
       .b = 0 
   };
}

inline Vector3 add_vec(Vector3 a, Vector3 b)
{
    return (Vector3) {
        .x = a.x + b.x,
        .y = a.y + b.y,
        .z = a.z + b.z
    };
}

inline Vector3 mul_vec(Vector3 v, double scalar)
{
    return (Vector3) {
        .x = v.x * scalar,
        .y = v.y * scalar,
        .z = v.z * scalar
    };
}

inline Vector3 sub_vec(Vector3 a, Vector3 b)
{
    return add_vec(a, mul_vec(b, -1));
}

inline double sq(double d)
{
    return d * d;
}

void normalize(Vector3* ray_dir)
{
    double magnitude = sqrt(sq(ray_dir -> x) + sq(ray_dir -> y) + sq(ray_dir -> z));

    ray_dir -> x /= magnitude;
    ray_dir -> y /= magnitude;
    ray_dir -> z /= magnitude;
}

Vector3 create_ray(Vector3 origin, Vector3 end)
{
    Vector3 ray_dir = (Vector3)
    {
        .x = end.x - origin.x,
        .y = end.y - origin.y,
        .z = end.z - origin.z
    };

    normalize(&ray_dir);

    return ray_dir;
}

bool cast(Vector3 origin, Vector3 ray_dir, Sphere s, double* t)
{
    Vector3 ray_diff = sub_vec(origin, s.c);

    double b = 2 * dot_vec(ray_dir, ray_diff);
    double c = dot_vec(ray_diff, ray_diff) - sq(s.r);

    double discrim = sq(b) - 4 * c;

    if (discrim < 0) 
        return false;

    discrim = sqrt(discrim);

    double sub = (-b - discrim) / 2.0;
    double add = (-b + discrim) / 2.0;

    // If both are pos, we want the min
    // Otherwise if one is negative we want the max because negatives are smaller
    *t = (sub > 0 && add > 0)
        ? min(sub, add)
        : max(sub, add);

    // If both are negative we didn't hit
    return *t > 0;
}

void write_file(Color** grid)
{
    FILE *fout = fopen("out.ppm", "w");

    fprintf(fout, "P3\n");
    fprintf(fout, "%d %d\n", MAX_X, MAX_Y);
    fprintf(fout, "255\n");

    for (size_t y = 0; y < MAX_Y; y++ )
    {
        for (size_t x = 0; x < MAX_X; x++)
        {
            Color c = grid[y][x];

            fprintf(fout, "%d %d %d\n", c.r, c.g, c.b);
        }
    }

    fclose(fout);
}

int main(void)
{
    Color** grid = malloc(sizeof(Color*) * MAX_Y);

    for (size_t i = 0; i < MAX_Y; i++)
       grid[i] = malloc(sizeof(Color) * MAX_X);

    // TODO make this not global
    init_objects();

    double aspect_ratio = (MAX_X * 1.0) / MAX_Y;

    #pragma omp parallel for
    for (size_t y = 0; y < MAX_Y; y++)
    {
        for (size_t x = 0; x < MAX_X; x++)
        {
            double px_scaled = (x + 0.5) / (1.0 * MAX_X);
            double py_scaled = ((MAX_Y - y) + 0.5) / (1.0 * MAX_Y);

            px_scaled *= aspect_ratio;

            Vector3 ray_dir = create_ray(eye, (Vector3) { .x = px_scaled, .y = py_scaled });

            Color c = { .r = 255, .g = 255, .b = 255 };

            double min_t = INFINITY;
            double t = 0;

            // foreach (Sphere a[i] in a)
            for (size_t i = 0; i < SPHERE_NUM; i++)
            {
                if (cast(eye, ray_dir, spheres[i], &t))
                {
                   if (t >= min_t) continue;

                   min_t = t;

                   c = spheres[i].h;
                }
            }

            // If we haven't hit something, end
            if (min_t == INFINITY)
            {
                grid[y][x] = c;
                continue;
            }

            // Calculate the intersection point with i = origin + t * dir
            // t is subtracted by a bit so we aren't inside the sphere
            Vector3 intersection = add_vec(eye, mul_vec(ray_dir, min_t - EPSILON));

            ray_dir = create_ray(intersection, light);

            // And then check if intersection -> light hits anything
            for (size_t i = 0; i < SPHERE_NUM; i++)
            {
                if (!cast(intersection, ray_dir, spheres[i], &t)) continue;

                c.r /= 2;
                c.g /= 2;
                c.b /= 2;

                break;
            }

            grid[y][x] = c;
        }
    }

    write_file(grid);

    for (size_t i = 0; i < MAX_Y; i++)
        free(grid[i]);

    free(grid);

    return 0;
}

//
// end of file
//

