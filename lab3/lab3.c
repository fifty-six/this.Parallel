//
// Bham, 3 September 2019
//

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define ON 'T'
#define OFF ' '
#define FIRE '*'
// #define SHOW

struct Point
{
    size_t x;
    size_t y;
};

struct Node
{
    struct Point pos;
    struct Node* next;
};

struct Queue
{
    struct Node* first;
    struct Node* last;
};

struct Point rem(struct Queue* q)
{
    if (q -> first == NULL)
    {
        fprintf(stderr, "Tried to pop from an empty queue!\n");

        __builtin_trap();
    }

    struct Node* first = q -> first;

    q -> first = first -> next;

    // Removing the only node in the Queue.
    if (q -> first == NULL)
        q -> last = NULL;

    struct Point p = first -> pos;

    free(first);

    return p;
}

void add(struct Queue* q, size_t x, size_t y)
{
    struct Point p = {.x = x, .y = y};

    struct Node* node = (struct Node*) malloc(sizeof(struct Node));
    node -> next = NULL;
    node -> pos = p;

    if (q -> first == NULL)
    {
        q -> first = q -> last = node;

        return;
    }

    q -> last -> next = node;
    q -> last = node;
}

double randDouble()
{
    return rand() * 1.0 / RAND_MAX;
}

void fill(size_t r, size_t c, char grid[][c], double p)
{
    for (size_t j = 0; j < r; j++)
    {
        for (size_t k = 0; k < c; k++)
        {
            grid[j][k] = randDouble() < p ? ON : OFF;
        }
    }
}

void msleep(int ms)  
{
    struct timespec ts = 
    {
        ms / 1000,
        (ms % 1000) * 1000000L
    };

    nanosleep(&ts, NULL);
}

void show(size_t r, size_t c, char grid[][c], int step)
{
#ifdef SHOW
    printf("\e[1;1H\e[2J");

    for (size_t j = 0; j < r; j++)
    {
        for (size_t k = 0; k < c; k++)
        {
            printf("%c", grid[j][k]);
        }
        printf("\n");
    }

    printf("%i\n", step);

    msleep(150);
#endif
}


void ignite(size_t r, size_t c, char grid[][c], struct Queue* q)
{
    for (size_t i = 0; i < r; i++)
    {
        if (grid[i][0] != ON) continue;

        grid[i][0] = FIRE;
        add(q, i, 0);
    } 
}

void spread(int* step, size_t r, size_t c, char grid[][c], struct Queue* fire)
{
    struct Queue newFire = {.first = NULL, .last = NULL};

    while (fire -> first != NULL)
    {
        struct Point p = rem(fire);

        int x = p.x;
        int y = p.y;

        if (y + 1 < c && grid[x][y + 1] == ON)
        {
            grid[x][y + 1] = FIRE;
            add(&newFire, x, y + 1);
        }
        if (y - 1 >= 0 && grid[x][y - 1] == ON)
        {
            grid[x][y - 1] = FIRE;
            add(&newFire, x, y - 1);
        }
        if (x + 1 < r && grid[x+1][y] == ON)
        {
            grid[x+1][y] = FIRE;
            add(&newFire, x + 1, y);
        }
        if (x - 1 >= 0 && grid[x-1][y] == ON)
        {
            grid[x-1][y] = FIRE;
            add(&newFire, x - 1, y);
        }

        grid[x][y] = OFF;
    }

    fire -> first = newFire.first;
    fire -> last = newFire.last;

    (*step)++;
}

double sim(double p, int r, int c)
{
    char grid[r][c];

    fill(r, c, grid, p);

    struct Queue fires = { .first = NULL, .last = NULL };

    int step = -1;

    show(r, c, grid, step);

    step = 0;

    ignite(r, c, grid, &fires);

    while(fires.first != NULL)
    {
        show(r, c, grid, step);
        spread(&step, r, c, grid, &fires);
    } 

    show(r, c, grid, step);

    return (step * 1.0)/c;
}

double avg(double p, int r, int c)
{
    const double count = 10000;

    double total = 0;

    for (int i = 0; i < count; i++)
        total += sim(p, r, c);

    return total / count;

}

int main()
{
    long long rseed = 1454734;

    srand(rseed);

    const int r = 30;
    const int c = 40;

    printf("P, 30x40, 60x80, 120x160\n");
    for (double p = 0.01; p <= 1; p += 0.005)
    {
        printf("%f, %f, %f, %f\n", p, avg(p, r, c), avg(p, r*2, c*2), avg(p, r*4, c*4));
        fprintf(stderr, "%f\n", p);
    }
    fflush(stdout);

    return 0;
}
