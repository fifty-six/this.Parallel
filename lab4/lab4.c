//
// Bham, 3 September 2019
//

#include <stdio.h>
#include <signal.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>

#define ON 'T'
#define OFF ' '
#define FIRE '*'
// Equivalent to (1/STEPS) probability = .005 = .5%
#define STEPS 201 

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

    step = 0;

    ignite(r, c, grid, &fires);

    while(fires.first != NULL)
    {
        spread(&step, r, c, grid, &fires);
    } 

    return (step * 1.0)/c;
}

int sum(double p, int r, int c, int count)
{
    double total = 0;

    for (int i = 0; i < count; i++)
        total += sim(p, r, c);

    return total;
}

const int totalTrials = 10000;

void master(int size)
{
    MPI_Status status;

    FILE *f = fopen("out.csv", "w");

    fprintf(f, "P, 30x40");

    int sums[STEPS];

    memset(&sums, 0, sizeof(int) * STEPS);

    int res[STEPS];

    for (int i = 0; i < (size - 1); i++)
    {
        MPI_Recv(&res, STEPS, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);

        printf("Recieved value from %i\n", status.MPI_SOURCE);

        for (int j = 0; j < STEPS; j++)
        {
            sums[j] += res[j];
        }
    }

    double dP = 1.0 / (STEPS - 1);
    for (int i = 0; i < STEPS; i++)
    {
        fprintf(f, "%f, %f\n", (dP * i), sums[i]*1.0 / totalTrials);
    }

    fclose(f);
}

void slave(int size)
{
    const int r = 30;
    const int c = 40;

    const int trials = (int) (totalTrials / (size - 1));

    printf("Size is %i so %i trials\n", size, trials);

    int res[STEPS];

    memset(&res, 0, sizeof(int) * STEPS);

    for (int p = 0; p < STEPS; p++)
    {
        res[p] = sum(p*1.0 / STEPS, r, c, trials);
    }

    MPI_Send(&res, STEPS, MPI_INT, 0, 0, MPI_COMM_WORLD);
}

void sigInt(int sig) 
{
    printf("Caught Ctrl-C. Aborting.\n");
    MPI_Finalize();
}

int main(int argc, char* argv[])
{
    long long rseed = 1454734;

    int size;
    int rank;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    srand(rseed + rank);

    signal(SIGINT, sigInt);

    if (rank == 0)
        master(size);
    else
        slave(size);

    MPI_Finalize();

    return 0;
}
