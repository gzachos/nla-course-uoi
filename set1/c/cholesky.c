#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <math.h>

/* Function Prototypes */
int      alloc_matrices(float ***a1, float ***a2, float **b1, float **b2, int n);
void     free_matrices(float **a1, float **a2, float *b1, float *b2, int n);
float   *alloc_1d_matrix(int n);
float  **alloc_2d_matrix(int n);
void     free_2d_matrix(float **arr, int n);
void     write_2d_matrix(char *filename, float **arr, int n);
void     write_1d_matrix(char *filename, float *arr, int n);
float  **cholesky_decomposition(float **a, int n);

int main(void)
{
	int i, j, n = 10;
	float **a1 = NULL,
	      **a2 = NULL,
	       *b1 = NULL,
	       *b2 = NULL,
	      **l1 = NULL;

	if (alloc_matrices(&a1, &a2, &b1, &b2, n) != 0)
	{
		return EXIT_FAILURE;
	}

	/* Init b1 */
	b1[0] = b1[n-1] = 3;
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

	write_2d_matrix("a1.txt", a1, n);
	write_2d_matrix("a2.txt", a2, n);
	write_1d_matrix("b1.txt", b1, n);
	write_1d_matrix("b2.txt", b2, n);

	l1 = cholesky_decomposition(a1, n);
#if 1
	printf("### L\n");
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < n; j++)
			printf("%10f ", l1[i][j]);
		printf("\n");
	}
	printf("\n");
#endif
	free_matrices(a1, a2, b1, b2, n);

	return EXIT_SUCCESS;
}


int alloc_matrices(float ***a1, float ***a2, float **b1, float **b2, int n)
{
	*a1 = alloc_2d_matrix(n);
	if (!(**a1))
		return 1;

	*a2 = alloc_2d_matrix(n);
	if (!(*a2))
	{
		free_2d_matrix(*a1, n);
		return 1;
	}

	*b1 = alloc_1d_matrix(n);
	if (!(*b1))
	{
		free_2d_matrix(*a1, n);
		free_2d_matrix(*a2, n);
		return 1;
	}

	*b2 = alloc_1d_matrix(n);
	if (!(*b2))
	{
		free_2d_matrix(*a1, n);
		free_2d_matrix(*a2, n);
		free(b1);
		return 1;
	}

	return 0;
}


float **alloc_2d_matrix(int n)
{
	int i, j;
	float **arr = (float **) malloc(n * sizeof(float *));

	if (!arr)
	{
		perror("malloc");
		return NULL;
	}

	for (i = 0; i < n; i++)
	{
		arr[i] = (float *) calloc(n, sizeof(float));
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


float  **alloc_lower_triangular_matrix(int n)
{
	int i, j;
	float **arr = (float **) malloc(n * sizeof(float *));

	if (!arr)
	{
		perror("malloc");
		return NULL;
	}

	for (i = 0; i < n; i++)
	{
		arr[i] = (float *) calloc(i+1, sizeof(float));
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


float *alloc_1d_matrix(int n)
{
	float *arr = (float *) calloc(n, sizeof(float));

	if (!arr)
	{
		perror("calloc");
		return NULL;
	}

	return arr;
}


void free_2d_matrix(float **arr, int n)
{
	int i;
	for (i = 0; i < n; i++)
		free(arr[i]);
	free(arr);
}


void free_matrices(float **a1, float **a2, float *b1, float *b2, int n)
{
	free_2d_matrix(a1, n);
	free_2d_matrix(a2, n);
	free(b1);
	free(b2);
}


void write_2d_matrix(char *filename, float **arr, int n)
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
			fprintf(outfile, "%10f  ", arr[i][j]);
		fprintf(outfile, "\n");
	}

	fclose(outfile);
}


void write_1d_matrix(char *filename, float *arr, int n)
{
	int i;

	FILE *outfile = fopen(filename, "w");
	if (!outfile)
	{
		perror("fopen");
		return;
	}

	for (i = 0; i < n; i++)
	{
		fprintf(outfile, "%f\n", arr[i]);
	}

	fclose(outfile);
}

#if 1

float  **cholesky_decomposition(float **a, int n)
{
	int i, j, k;
	float sum, **l = alloc_2d_matrix(n);

	for (i = 1; i <= n; i++)
	{
		for (j = 1; j <= i-1; j++)
		{
			sum = a[i-1][j-1];
			for (k = 1; k <= j-1; k++)
			{
				sum -= l[i-1][k-1] * l[j-1][k-1];
			}
			l[i-1][j-1] = ((float)sum) / l[j-1][j-1];
		}
		sum = a[i-1][i-1];
		for (j = 1; j <= i-1; j++)
		{
			sum -= l[i-1][j-1] * l[i-1][j-1];
		}
		l[i-1][i-1] = sqrt(sum);
	}
	return l;
}

#else

float  **cholesky_decomposition(float **a, int n)
{
	int i, j, k;
	float sum, **l = alloc_2d_matrix(n);

	for (i = 0; i < n; i++)
	{
		for (j = 0; j <= i; j++)
		{
			sum = 0;
			if (j == i)
			{
				for (k = 0; k < j; k++)
					sum += l[i][k] * l[i][k];
				l[i][i] = sqrt(a[i][i] - sum);
			}
			else
			{
				for (k = 0; k < j; k++)
					sum += l[i][k] * l[j][k];
				l[i][j] = (a[i][j] - sum) / l[j][j];
			}
		}
	}
	return l;
}

#endif
