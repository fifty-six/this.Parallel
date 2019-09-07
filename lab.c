//
// Bham, 3 September 2019
//

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

double randDouble()
{
    return rand() * 1.0 / RAND_MAX;
}

void fill(int r, int c, char grid[][c], double p)
{
    for (int j = 0; j < r; j++)
    {
        for (int k = 0; k < c; k++)
        {
            grid[j][k] = randDouble() < p ? 'T' : ' ';
        }
    }
}

void show(int r, int c, char grid[][c])
{
    for (int j = 0; j < r; j++)
    {
        for (int k = 0; k < c; k++)
        {
            printf("%c", grid[j][k]);
        }
        printf("\n");
    }
}

int count(int r, int c, char grid[][c])
{
    int count = 0;
    for (int i = 0; i < r; i++)
    {
        for (int j = 0; j < c; j++)
        {
            count += grid[i][j] == 'T';
        }
    }
    return count;
}

int colCount(int r, int c, char grid[][c], int col)
{
    int count = 0;
    for (int i = 0; i < r; i++)
    {
        count += grid[i][col] == 'T';
    }
    return count;
} 

int adjacentCol(int r, int c, char grid[][c], int col)
{
    int count = 0;
    for (int i = 0; i < r; i++)
    {
        count += grid[i][col] == 'T' && grid[i][col - 1] == 'T';
    }
    return count;
}

int main()
{
    long long int rseed = 1454734;

    const int r = 30;
    const int c = 40;

    srand(rseed);

    char grid[r][c];

    fill(r, c, grid, 0.60);
    show(r, c, grid);

    int countVal = count(r, c, grid);

    printf("Count: %i\n", countVal);
    printf("Percent: %f\n", (countVal*1.0/(r*c)) * 100);
    printf("Seed: %lli\n", rseed); 
    printf("%i/30 in the left column\n", colCount(r, c, grid, 0));
    printf("%i/30 in the next column\n", colCount(r, c, grid, 1));
    printf("%i of these are adjacent to activated cells in the first column.\n", adjacentCol(r, c, grid, 1));

    int totalCells = 0;
    int totalActive = 0;
    for (int i = 0; i < 100; i++)
    {
        fill(r, c, grid, 0.60);

        totalActive += count(r, c, grid);
        totalCells += r * c;
    }

    printf("Percent for 100 grids: %f\n", ((totalActive * 1.0) / totalCells) * 100);

    return 0;
}
