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
#define SHADOW 0.8
#define REFLECT 0.6
#define WIDTH 0.1

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

#define SPHERE_NUM 5

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

   // Light
   spheres[4].c = (Vector3) { .x = 0.00, .y = 1.25, .z = -0.50 }; 
   spheres[4].r = 0.2;
   spheres[4].h = (Color) {
       .r = 255,
       .g = 255,
       .b = 255
   };
}

int mercator_w, mercator_h;
Color** mercator;

void write_file(Color** grid, char* fname, size_t max_y, size_t max_x);

void init_globe()
{
    FILE* f = fopen("mercator82.ppm", "r");

    if (!f)
    {
        fprintf(stderr, "Meracator82 doesn't exist!");

        exit(-1);

        return;
    }

    fscanf(f, "P3\n");
    fscanf(f, "%d %d\n", &mercator_w, &mercator_h);

    mercator = malloc((mercator_h + 1)* sizeof(Color*));

    for (int i = 0; i < mercator_h; i++)
    {
        mercator[i] = malloc((mercator_w + 1) * sizeof(Color));
    }

    // We know the file is 255 color
    fscanf(f, "255\n");

    int r, g, b;

    for (int i = 0; i < mercator_h; i++)
    {
        for (int j = 0; j < mercator_w; j++)
        {
            fscanf(f, "%d %d %d\n", &r, &g, &b);

            mercator[i][j] = (Color) {
                .r = r,
                .g = g,
                .b = b
            };
        }
    }

    printf("mercator_h: %d, mercator_w: %d\n", mercator_h, mercator_w);

    write_file(mercator, "clone.ppm", mercator_h, mercator_w);
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

void write_file(Color** grid, char* fname, size_t max_y, size_t max_x)
{
    FILE *fout = fopen(fname, "w");

    fprintf(fout, "P3\n");
    fprintf(fout, "%d %d\n", max_x, max_y);
    fprintf(fout, "255\n");

    for (size_t y = 0; y < max_y; y++ )
    {
        for (size_t x = 0; x < max_x; x++)
        {
            Color c = grid[y][x];

            fprintf(fout, "%d %d %d\n", c.r, c.g, c.b);
        }
    }

    fclose(fout);
}

bool shadow(Vector3 intersection, Vector3 light_dir, Color* c)
{
    double t = 0.0;

    for (size_t i = 0; i < SPHERE_NUM - 1; i++)
    {
        if (!cast(intersection, light_dir, spheres[i], &t)) continue;

        c -> r *= 1 - SHADOW;
        c -> g *= 1 - SHADOW;
        c -> b *= 1 - SHADOW;

        return true;
    }

    return false;
}

Color get_color(Vector3 origin, Vector3 ray_dir, int depth, int max_depth)
{

    Color c = { .r = 0, .g = 0, .b = 0 };

    double min_t = INFINITY;
    double t = 0;

    int sphere_ind = 0;

    // foreach (Sphere a[i] in a)
    for (size_t i = 0; i < SPHERE_NUM; i++)
    {
        if (cast(origin, ray_dir, spheres[i], &t))
        {
            if (t >= min_t) continue;

            min_t = t;

            sphere_ind = i;

            c = spheres[i].h;
        }
    }

    // If we haven't hit something, end
    if (min_t == INFINITY || sphere_ind == SPHERE_NUM - 1)
    {
        return c;
    }

    Sphere sphere = spheres[sphere_ind];

    // Calculate the intersection point with i = origin + t * dir
    // t is subtracted by a bit so we aren't inside the sphere
    Vector3 intersection = add_vec(origin, mul_vec(ray_dir, min_t - EPSILON));

    if (sphere_ind == 0)
    {
        if ( ((int) (round(intersection.x / WIDTH) + round(intersection.z / WIDTH))) % 2 == 0)
            c = (Color) { .r = 255, .g = 255, .b = 255 };
    }

    Vector3 light_dir = create_ray(intersection, light);

    // Surface normal
    Vector3 gradient = sub_vec(intersection, sphere.c);

    normalize(&gradient);

    // Blue sphere becomes the earth
    if (sphere_ind == 1)
    {
        double lat_ang = acos(gradient.y);
        // -pi to pi -> 0 to 2pi
        double long_ang = atan2(gradient.z, gradient.x) + M_PI;

        int latitude  = (int) (lat_ang * mercator_h / M_PI);
        int longitude = (int) (long_ang * mercator_w / (2.0 * M_PI));

        c = mercator[latitude][longitude % mercator_w];
    }

    Vector3 new_ray = sub_vec(ray_dir, mul_vec(gradient, 2 * dot_vec(gradient, ray_dir)));

    Color reflect = get_color(intersection, new_ray, depth + 1, max_depth);

    c.r *= 1 - REFLECT;
    c.g *= 1 - REFLECT;
    c.b *= 1 - REFLECT;

    c.r += REFLECT * reflect.r;
    c.g += REFLECT * reflect.g;
    c.b += REFLECT * reflect.b;

    // And then check if intersection -> light hits anything
    if (!shadow(intersection, light_dir, &c))
    {
        // |grad| = 1, |ray_dir| = 1, 1 * 1 * cos(\theta) = cos(\theta).
        double cos_theta = dot_vec(gradient, light_dir);

        if (cos_theta < 0)
            cos_theta = 0;

        c.r *= (1 - SHADOW) + SHADOW * cos_theta;
        c.g *= (1 - SHADOW) + SHADOW * cos_theta;
        c.b *= (1 - SHADOW) + SHADOW * cos_theta;
    }

    return c;
}

int main(void)
{
    Color** grid = malloc(sizeof(Color*) * MAX_Y);

    for (size_t i = 0; i < MAX_Y; i++)
       grid[i] = malloc(sizeof(Color) * MAX_X);

    // TODO make this not global
    init_objects();
    init_globe();

    double aspect_ratio = (MAX_X * 1.0) / MAX_Y;

    #pragma omp parallel for
    for (size_t y = 0; y < MAX_Y; y++)
    {
        for (size_t x = 0; x < MAX_X; x++)
        {
            double px_scaled = -.25 + (x + 0.5) / (1.0 * MAX_X);
            double py_scaled = ((MAX_Y - y) + 0.5) / (1.0 * MAX_Y);

            px_scaled *= aspect_ratio;

            Vector3 ray_dir = create_ray(eye, (Vector3) { .x = px_scaled, .y = py_scaled });

            grid[y][x] = get_color(eye, ray_dir, 0, 6);
        }
    }

    write_file(grid, "out.ppm", MAX_Y, MAX_X);

    for (size_t i = 0; i < MAX_Y; i++)
        free(grid[i]);

    free(grid);

    return 0;
}

//
// end of file
//

