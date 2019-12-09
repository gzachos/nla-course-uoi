#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <math.h>

/* Function Prototypes */
int       alloc_matrices(double ***a1, double ***a2, double **b1, double **b2, int n);
void      free_matrices(double **a1, double **a2, double *b1, double *b2, int n);
double   *alloc_1d_matrix(int n);
double  **alloc_2d_matrix(int n);
void      free_2d_matrix(double **arr, int n);
void      write_2d_matrix(char *filename, double **arr, int n);
void      write_1d_matrix(char *filename, double *arr, int n);
double  **cholesky_decomposition(double **a, int n);
double   *forward_substitution(double **l, double *b, int n);
double   *back_substitution(double **l, double *b, int n);
double  **transpose(double **arr, int n);


int main(void)
{
	int i, j, k, n = 10000;
	double **a1 = NULL,
	       **a2 = NULL,
	        *b1 = NULL,
	        *b2 = NULL,
	       **l1 = NULL,
	       **a  = NULL,
	         sum;

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

	/* Write a1, a2, b1 and b2 to files */
	write_2d_matrix("a1.txt", a1, n);
	write_2d_matrix("a2.txt", a2, n);
	write_1d_matrix("b1.txt", b1, n);
	write_1d_matrix("b2.txt", b2, n);

	l1 = cholesky_decomposition(a1, n);
#if 1
	printf("\nL =\n\n");
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < n; j++)
			printf("%10f ", l1[i][j]);
		printf("\n");
	}
	printf("\n");
#endif
#if 0
	a = alloc_2d_matrix(n);

	/* Calculate L * L-transpose */
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < n; j++)
		{
			sum = 0.0;
			for (k = 0; k < n; k++) {
				sum += l1[i][k] * l1[j][k];
			}
			a[i][j] = sum;
		}
	}

	printf("\nA = L * L-transpose\n\n");
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < n; j++)
			printf("%10f ", a[i][j]);
		printf("\n");
	}
	printf("\n");
#endif

	/* L*y=b (forward substitution) */
	double *y1 = forward_substitution(l1, b1, n);
#if 1
	printf("\nY =\n\n");
	for (i = 0; i < n; i++)
		printf("%10f\n", y1[i]);
#endif
	double **lt = transpose(l1, n);
#if 0
	printf("\nL-transpose\n\n");
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < n; j++)
			printf("%10f ", lt[i][j]);
		printf("\n");
	}
	printf("\n");
#endif
	/* LT*x=y (back substitution) */
	double *x1 = back_substitution(lt, y1, n);
#if 1
	printf("\nX =\n\n");
	for (i = 0; i < n; i++)
		printf("%10f\n", x1[i]);

	free_matrices(a1, a2, b1, b2, n);
#endif
	return EXIT_SUCCESS;
}

double **transpose(double **arr, int n)
{
	double **trn = alloc_2d_matrix(n);
	int i, j;

	for (i = 0; i < n; i++)
		for (j = 0; j < n; j++)
			trn[i][j] = arr[j][i];
	return trn;
}

double *forward_substitution(double **l, double *b, int n)
{
	double *y = alloc_1d_matrix(n), sum;
	int i, k;

	for (i = 0; i < n; i++)
	{
		sum = b[i];
		for (k = 0; k < i; k++)
			sum -= l[i][k] * y[k];
		y[i] = sum / l[i][i];
	}

	return y;
}


double *back_substitution(double **l, double *b, int n)
{
	double *y = alloc_1d_matrix(n), sum;
	int i, k;

	for (i = n; i >= 1; i--)
	{
		sum = b[i-1];
		for (k = n; k >= i+1; k--)
			sum -= l[i-1][k-1] * y[k-1];
		y[i-1] = sum / l[i-1][i-1];
	}

	return y;
}


int alloc_matrices(double ***a1, double ***a2, double **b1, double **b2, int n)
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


double **alloc_2d_matrix(int n)
{
	int i, j;
	double **arr = (double **) malloc(n * sizeof(double *));

	if (!arr)
	{
		perror("malloc");
		return NULL;
	}

	for (i = 0; i < n; i++)
	{
		arr[i] = (double *) calloc(n, sizeof(double));
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


double  **alloc_lower_triangular_matrix(int n)
{
	int i, j;
	double **arr = (double **) malloc(n * sizeof(double *));

	if (!arr)
	{
		perror("malloc");
		return NULL;
	}

	for (i = 0; i < n; i++)
	{
		arr[i] = (double *) calloc(i+1, sizeof(double));
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


double *alloc_1d_matrix(int n)
{
	double *arr = (double *) calloc(n, sizeof(double));

	if (!arr)
	{
		perror("calloc");
		return NULL;
	}

	return arr;
}


void free_2d_matrix(double **arr, int n)
{
	int i;
	for (i = 0; i < n; i++)
		free(arr[i]);
	free(arr);
}


void free_matrices(double **a1, double **a2, double *b1, double *b2, int n)
{
	free_2d_matrix(a1, n);
	free_2d_matrix(a2, n);
	free(b1);
	free(b2);
}


void write_2d_matrix(char *filename, double **arr, int n)
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


void write_1d_matrix(char *filename, double *arr, int n)
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

double  **cholesky_decomposition(double **a, int n)
{
	int i, j, k;
	double sum, **l = alloc_2d_matrix(n);

	for (i = 1; i <= n; i++)
	{
		for (j = 1; j <= i-1; j++)
		{
			sum = a[i-1][j-1];
			for (k = 1; k <= j-1; k++)
			{
				sum -= l[i-1][k-1] * l[j-1][k-1];
			}
			l[i-1][j-1] = sum / l[j-1][j-1];
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

double  **cholesky_decomposition(double **a, int n)
{
	int i, j, k;
	double sum, **l = alloc_2d_matrix(n);

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
