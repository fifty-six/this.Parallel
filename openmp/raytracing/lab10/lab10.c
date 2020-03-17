//
// Yusuf Bham, January 2020
//

#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#define MAX_X 1920
#define MAX_Y 1080
#define EPSILON 0.001
#define SHADOW 0.7
#define REFLECT 0.7
#define WIDTH 0.1
#define MAX_FRAMES 600

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

    // Angle;
    double theta;

} Sphere;

typedef struct
{
    Vector3 normal;

    Vector3 vertices[3];

    Color h;

} Triangle;

#define SPHERE_NUM 5

Sphere spheres[SPHERE_NUM];

Vector3 eye = { .x = 0.50, .y = 0.50, .z = -1.00 };
Vector3 light = { .x = 0.00, .y = 1.25, .z = -0.50 }; 

double dot_vec( Vector3 t , Vector3 u )
{
   return t.x * u.x + t.y * u.y + t.z * u.z ;
}

static inline double min(double a, double b)
{
    return a < b ? a : b;
}

static inline double max(double a, double b)
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
       .r = 234,
       .g = 21,
       .b = 81 
   };
   spheres[0].theta = 0;

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
   spheres[1].theta = 0;

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
   spheres[2].theta = 0;

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
   spheres[4].theta = 0;
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
}

Triangle* triangles;
size_t tri_count;

void init_triangles()
{
    FILE* f_vertices = fopen("vertices.txt", "r");

    if (!f_vertices)
    {
        fprintf(stderr, "vertices.txt not found");
        exit(-1);
    }

    int vert_count;

    fscanf(f_vertices, "%i\n", &vert_count);

    Vector3* verts = malloc(vert_count * sizeof(Vector3));

    for (int i = 0; i < vert_count; i++)
    {
        double x, y, z;

        fscanf(f_vertices, "%lf %lf %lf\n", &x, &y, &z);

        verts[i] = (Vector3) { .x = x, .y = y, .z = z };
    }

    fclose(f_vertices);

    FILE* f_faces = fopen("faces.txt", "r");
    
    if (!f_faces)
    {
        fprintf(stderr, "faces.txt not found");
        exit(-1);
    }

    // Number of faces == Number of triangles
    fscanf(f_faces, "%lu\n", &tri_count);

    // 3 vertices for each face
    size_t* faces = malloc(sizeof(size_t) * tri_count * 3);

    for (size_t i = 0; i < tri_count * 3; i += 3)
    {
        fscanf(f_faces, "%lu %lu %lu\n", &faces[i], &faces[i+1], &faces[i+2]);
    }

    FILE* f_normals = fopen("normals.txt", "r");

    if (!f_normals)
    {
        fprintf(stderr, "normals.txt not found");
        exit(-1);
    }

    Vector3* normals = malloc(sizeof(Vector3) * tri_count);

    size_t normal_count;

    fscanf(f_normals, "%lu\n", &normal_count);

    assert(normal_count == tri_count);

    for (size_t i = 0; i < tri_count; i++)
    {
        double x, y, z;

        fscanf(f_normals, "%lf %lf %lf\n", &x, &y, &z);

        normals[i] = (Vector3) { .x = x, .y = y, .z = z };
    }

    triangles = malloc(sizeof(Triangle) * tri_count);

    for (size_t i = 0; i < tri_count; i++)
    {
        triangles[i] = (Triangle) {
            .h      = (Color) { .r = 0, .g = 0, .b = 255 },
            .normal = normals[i]
        };

        for (size_t j = 0; j < 3; j++)
            triangles[i].vertices[j] = verts[faces[3 * i + j]];
    }
}

static inline Vector3 add_vec(Vector3 a, Vector3 b)
{
    return (Vector3) {
        .x = a.x + b.x,
        .y = a.y + b.y,
        .z = a.z + b.z
    };
}

static inline Vector3 mul_vec(Vector3 v, double scalar)
{
    return (Vector3) {
        .x = v.x * scalar,
        .y = v.y * scalar,
        .z = v.z * scalar
    };
}

static inline Vector3 sub_vec(Vector3 a, Vector3 b)
{
    return add_vec(a, mul_vec(b, -1));
}

static inline double sq(double d)
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

bool tri_cast(Vector3 origin, Vector3 ray_dir, Triangle tri, double* t, double min_t)
{
    *t = dot_vec(tri.normal, sub_vec(tri.vertices[0], origin)) / dot_vec(tri.normal, ray_dir);

    if (*t > min_t || *t < 0)
        return false;

    Vector3 p = add_vec(origin, mul_vec(ray_dir, *t));

    Vector3 w = sub_vec(p,               tri.vertices[0]);
    Vector3 u = sub_vec(tri.vertices[1], tri.vertices[0]);
    Vector3 v = sub_vec(tri.vertices[2], tri.vertices[0]);

    double wu = dot_vec(w, u);
    double vv = dot_vec(v, v);
    double uv = dot_vec(u, v);
    double wv = dot_vec(w, v);
    double uu = dot_vec(u, u);

    double denom = (uu*vv - uv*uv);

    double alpha = (wu*vv - uv*wv)/denom;

    // fprintf(stderr, "Alpha: %f\n", alpha);

    if (alpha < 0)//  || alpha > 1)
        return false;

    double beta  = (wv*uu - uv*wu)/denom;

//    fprintf(stderr, "Beta: %f\n", beta);

    if (beta < 0)// || beta > 1)
        return false;

 //   fprintf(stderr, "alpha+beta: %f\n", alpha+beta);

    return alpha + beta <= 1;
}

void write_file(Color** grid, char* fname, size_t max_y, size_t max_x)
{
    FILE *fout = fopen(fname, "w");

    fprintf(fout, "P3\n");
    fprintf(fout, "%d %d\n", (int) max_x, (int) max_y);
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

    for (size_t i = 0; i < tri_count; i++)
    {
        if (!tri_cast(intersection, light_dir, triangles[i], &t, 0)) continue;

        c -> r *= 1 - SHADOW;
        c -> g *= 1 - SHADOW;
        c -> b *= 1 - SHADOW;

        return true;
    }

    return false;
}

void rotate_sphere(Vector3 gradient, Sphere sphere, Color* c)
{
    Vector3 grad_rotate = gradient;

    double c_theta = cos(sphere.theta);
    double s_theta = sin(sphere.theta);

    grad_rotate.x = c_theta * gradient.x - s_theta * gradient.y;
    grad_rotate.y = s_theta * gradient.x + c_theta * gradient.y;

    double lat_ang = acos(grad_rotate.y);
    // -pi to pi -> 0 to 2pi
    double long_ang = atan2(grad_rotate.z, grad_rotate.x) + M_PI;

    int latitude  = (int) (lat_ang * mercator_h / M_PI);
    int longitude = (int) (long_ang * mercator_w / (2.0 * M_PI));

    *c = mercator[abs((latitude - 20) % mercator_h)][abs((longitude + 35) % mercator_w)];
}


Color get_color(Vector3 origin, Vector3 ray_dir, int depth, int max_depth)
{
    Color c = { .r = 0, .g = 0, .b = 0 };

    double min_t = INFINITY;
    double t = 0;

    int sphere_ind = 0;
    int triangle_ind = 0;

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

    for (size_t i = 0; i < tri_count; i++)
    {
        if (tri_cast(origin, ray_dir, triangles[i], &t, min_t))
        {
            min_t = t;

            triangle_ind = i;
            sphere_ind = -1;

            c = triangles[i].h;
        }
    }

    // If we haven't hit something, end
    if (min_t == INFINITY || sphere_ind == SPHERE_NUM - 1)
    {
        return c;
    }

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
    Vector3 gradient = sphere_ind != -1 
        ? sub_vec(intersection, spheres[sphere_ind].c) 
        : triangles[triangle_ind].normal;

    normalize(&gradient);

    // Blue sphere becomes the earth
    if (sphere_ind == 1)
        rotate_sphere(gradient, spheres[sphere_ind], &c);

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
    init_triangles();

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

            grid[y][x] = get_color(eye, ray_dir, 0, 10);
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

