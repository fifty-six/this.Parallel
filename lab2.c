//
// Bham, 3 September 2019
//

#include <stdio.h>
#include <unistd.h>
#include <stdbool.h>
#include <stdlib.h>
#include <time.h>

#define ON 'T'
#define OFF ' '
#define FIRE '*'
#define PREFIRE '^'

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

void show(size_t r, size_t c, char grid[][c], int step)
{
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

    // sleep(1);
}


void ignite(size_t r, size_t c, char grid[][c])
{
    for (size_t i = 0; i < r; i++)
    {
        grid[i][0] = grid[i][0] == ON ? FIRE : grid[i][0];
    } 
}

void spread(int* step, size_t r, size_t c, char grid[][c])
{
    for (size_t i = 0; i < r; i++)
    {
        for (size_t j = 0; j < c; j++)
        {
            if (grid[i][j] != FIRE) continue;

            if (j < c && grid[i][j + 1] == ON)
                grid[i][j + 1] = PREFIRE;
            if (j > 0 && grid[i][j - 1] == ON)
                grid[i][j - 1] = PREFIRE;
            if (i < r && grid[i+1][j] == ON)
                grid[i+1][j] = PREFIRE;
            if (i > 0 && grid[i-1][j] == ON)
                grid[i-1][j] = PREFIRE;

            grid[i][j] = OFF;
        }
    }

    for (size_t i = 0; i < r; i++)
    {
        for (size_t j = 0; j < c; j++)
        {
            if (grid[i][j] == PREFIRE)
                grid[i][j] = FIRE;
        }
    }

    (*step)++;
}

bool doneSpread(size_t r, size_t c, char grid[][c])
{
    for (size_t i = 0; i < r; i++)
        for (size_t j = 0; j < c; j++)
            if (grid[i][j] == FIRE)
                return false;
    return true;
}

int main()
{
    long long rseed = 1454734;

    const size_t r = 30;
    const size_t c = 40;

    srand(rseed);

    char grid[r][c];

    fill(r, c, grid, 0.60);

    int step = -1;

    show(r, c, grid, step);

    step = 0;
    
    ignite(r, c, grid);

    while(!doneSpread(r, c, grid))
    {
        show(r, c, grid, step);
        spread(&step, r, c, grid);
    }

    show(r, c, grid, step);

    printf("Steps: %i/%lu, %f\n", step, c, (step * 1.0/c));

    return 0;
}
