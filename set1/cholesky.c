#include <stdio.h>
#include <stdlib.h>
#include <errno.h>

/* Function Prototypes */
int    alloc_arrays(int ***a1, int ***a2, int **b1, int **b2, int n);
void   free_arrays(int **a1, int **a2, int *b1, int *b2, int n);
int   *alloc_1d_array(int n);
int  **alloc_2d_array(int n);
void   free_2d_array(int **arr, int n);
void   write_2d_array(char * filename, int **arr, int n);


int main(void)
{
	int i, j, n = 100;
	int **a1 = NULL,
	    **a2 = NULL,
	     *b1 = NULL,
	     *b2 = NULL;

	if (alloc_arrays(&a1, &a2, &b1, &b2, n) != 0)
	{
		return EXIT_FAILURE;
	}

	/* Init b1 */
	b1[0] = b1[n-1]   = 3;
	b1[1] = b1[n-2] = -1;


	/* Init b2 */
	b2[0] = b2[n-1] = 4;
	for (i = 2; i < n-2; i++)
		b2[i] = 1;

	/* Init a1, a2 */
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < n; j++)
		{
			if (i == j)
			{
				a1[i][j] = 6;
				a2[i][j] = 7;
			}
			else if (i+j == (i<<1)-1 || i+j == (j<<1)-1)
				a1[i][j] = a2[i][j] = -4;
			else if (i+j == (i<<1)-2 || i+j == (j<<1)-2)
				a1[i][j] = a2[i][j] = 1;
		}
	}

	write_2d_array("a1.txt", a1, n);
	write_2d_array("a2.txt", a2, n);

	free_arrays(a1, a2, b1, b2, n);

	return EXIT_SUCCESS;
}


int alloc_arrays(int ***a1, int ***a2, int **b1, int **b2, int n)
{
	*a1 = alloc_2d_array(n);
	if (!(**a1))
		return 1;

	*a2 = alloc_2d_array(n);
	if (!(*a2))
	{
		free_2d_array(*a1, n);
		return 1;
	}

	*b1 = alloc_1d_array(n);
	if (!(*b1))
	{
		free_2d_array(*a1, n);
		free_2d_array(*a2, n);
		return 1;
	}

	*b2 = alloc_1d_array(n);
	if (!(*b2))
	{
		free_2d_array(*a1, n);
		free_2d_array(*a2, n);
		free(b1);
		return 1;
	}

	return 0;
}


int **alloc_2d_array(int n)
{
	int i, j;
	int **arr = (int **) calloc(n, sizeof(int *));

	if (!arr)
	{
		perror("calloc");
		return NULL;
	}

	for (i = 0; i < n; i++)
	{
		arr[i] = (int *) calloc(n, sizeof(int));
		if (!arr[i])
		{
			perror("calloc");
			for (j = 0; j < i; j++)
				free(arr[j]);
			free(arr);
			return NULL;
		}
	}
	return arr;
}


int *alloc_1d_array(int n)
{
	int *arr = (int *) calloc(n, sizeof(int));

	if (!arr)
	{
		perror("calloc");
		return NULL;
	}

	return arr;
}


void free_2d_array(int **arr, int n)
{
	int i;
	for (i = 0; i < n; i++)
		free(arr[i]);
	free(arr);
}


void free_arrays(int **a1, int **a2, int *b1, int *b2, int n)
{
	free_2d_array(a1, n);
	free_2d_array(a2, n);
	free(b1);
	free(b2);
}


void write_2d_array(char *filename, int **arr, int n)
{
	int i, j;

	FILE *outfile = fopen(filename, "w");
	if (!outfile)
	{
		perror("fopen");
		return;
	}

	for (i = 0; i < n; i++)
	{
		for (j = 0; j < n; j++)
			fprintf(outfile, "%2d ", arr[i][j]);
		fprintf(outfile, "\n");
	}

	fclose(outfile);
}

