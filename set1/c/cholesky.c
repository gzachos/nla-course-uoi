#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <math.h>

#define BUFF_SIZE	64

#if 1
	typedef double fptype;
#else
	typedef float fptype;
#endif

/* Function Prototypes */
int       alloc_matrices(fptype ***a1, fptype ***a2, fptype **b1, fptype **b2, int n);
void      free_matrices(fptype **a1, fptype **a2, fptype *b1, fptype *b2, int n);
fptype   *alloc_1d_matrix(int n);
fptype  **alloc_2d_matrix(int n);
void      free_2d_matrix(fptype **arr, int n);
void      write_2d_matrix(char *filename, fptype **arr, int n);
void      write_1d_matrix(char *filename, fptype *arr, int n);
fptype  **cholesky_decomposition(fptype **a, int n);
fptype   *forward_substitution(fptype **l, fptype *b, int n);
fptype   *back_substitution(fptype **l, fptype *b, int n);
fptype  **transpose(fptype **arr, int n);
void      init_matrices(fptype **a1, fptype **a2, fptype *b1, fptype *b2, int n);


int main(int argc, char **argv)
{
	int i, j, k, n;
	fptype **a1  = NULL, **a2  = NULL,
	        *b1  = NULL,  *b2  = NULL,
	        *x1  = NULL,  *x2  = NULL,
	        *y1  = NULL,  *y2  = NULL,
	       **l1  = NULL, **l2  = NULL,
	       **lt1 = NULL, **lt2 = NULL,
	       **ra1 = NULL, **ra2 = NULL,
	         sum;
	char filename[64];

	if (argc != 2)
	{
		fprintf(stderr, "Usage: %s N\t(N > 0)\n", argv[0]);
		return EXIT_SUCCESS;
	}

	if ((n = atoi(argv[1])) <= 0)
	{
		fprintf(stderr, "Matrix size (N) should be positive!\n");
		return EXIT_SUCCESS;
	}
	else
		printf("N = %d\n", n);

	if (alloc_matrices(&a1, &a2, &b1, &b2, n) != 0)
	{
		return EXIT_FAILURE;
	}

	/* Initialize matrices */
	init_matrices(a1, a2, b1, b2, n);

	/* Write a1, a2, b1 and b2 to files */
#if 1
	snprintf(filename, BUFF_SIZE, "a1_%d.txt", n);
	write_2d_matrix(filename, a1, n);
	snprintf(filename, BUFF_SIZE, "a2_%d.txt", n);
	write_2d_matrix(filename, a2, n);
	snprintf(filename, BUFF_SIZE, "b1_%d.txt", n);
	write_1d_matrix(filename, b1, n);
	snprintf(filename, BUFF_SIZE, "b2_%d.txt", n);
	write_1d_matrix(filename, b2, n);
#endif
	l1 = cholesky_decomposition(a1, n);
#if 1
	ra1 = alloc_2d_matrix(n);

	/* Calculate L * L-transpose */
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < n; j++)
		{
			for (k = 0, sum = 0.0; k < n; k++)
				sum += l1[i][k] * l1[j][k];
			ra1[i][j] = sum;
		}
	}

	snprintf(filename, BUFF_SIZE, "ra1_%d.txt", n);
	write_2d_matrix(filename, ra1, n);
#endif

	/* L*y=b (forward substitution) */
	y1  = forward_substitution(l1, b1, n);
	lt1 = transpose(l1, n);

	/* LT*x=y (back substitution) */
	x1 = back_substitution(lt1, y1, n);

#if 0
	printf("\nY =\n\n");
	for (i = 0; i < n; i++)
		printf("%10f\n", y1[i]);
	printf("\nX =\n\n");
	for (i = 0; i < n; i++)
		printf("%10f\n", x1[i]);
#endif

	snprintf(filename, BUFF_SIZE, "l1_%d.txt", n);
	write_2d_matrix(filename, l1, n);
	snprintf(filename, BUFF_SIZE, "lt1_%d.txt", n);
	write_2d_matrix(filename, lt1, n);
	snprintf(filename, BUFF_SIZE, "x1_%d.txt", n);
	write_1d_matrix(filename, x1, n);
	snprintf(filename, BUFF_SIZE, "y1_%d.txt", n);
	write_1d_matrix(filename, y1, n);

	free_2d_matrix(l1, n);  free_2d_matrix(l2, n);
	free_2d_matrix(lt1, n); free_2d_matrix(lt2, n);
	free_2d_matrix(ra1, n); free_2d_matrix(ra2, n);
	free(y1); free(y2);
	free(x1); free(x2);
	free_matrices(a1, a2, b1, b2, n);

	return EXIT_SUCCESS;
}


void init_matrices(fptype **a1, fptype **a2, fptype *b1, fptype *b2, int n)
{
	int i, j;

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
}


fptype **transpose(fptype **arr, int n)
{
	fptype **trn = alloc_2d_matrix(n);
	int i, j;

	for (i = 0; i < n; i++)
		for (j = 0; j < n; j++)
			trn[i][j] = arr[j][i];
	return trn;
}


fptype *forward_substitution(fptype **l, fptype *b, int n)
{
	fptype *y = alloc_1d_matrix(n), sum;
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


fptype *back_substitution(fptype **l, fptype *b, int n)
{
	fptype *y = alloc_1d_matrix(n), sum;
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


int alloc_matrices(fptype ***a1, fptype ***a2, fptype **b1, fptype **b2, int n)
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


fptype **alloc_2d_matrix(int n)
{
	int i, j;
	fptype **arr = (fptype **) malloc(n * sizeof(fptype *));

	if (!arr)
	{
		perror("malloc");
		return NULL;
	}

	for (i = 0; i < n; i++)
	{
		arr[i] = (fptype *) calloc(n, sizeof(fptype));
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


fptype *alloc_1d_matrix(int n)
{
	fptype *arr = (fptype *) calloc(n, sizeof(fptype));

	if (!arr)
	{
		perror("calloc");
		return NULL;
	}

	return arr;
}


void free_2d_matrix(fptype **arr, int n)
{
	int i;

	if (!arr)
		return;

	for (i = 0; i < n; i++)
		free(arr[i]);
	free(arr);
}


void free_matrices(fptype **a1, fptype **a2, fptype *b1, fptype *b2, int n)
{
	free_2d_matrix(a1, n);
	free_2d_matrix(a2, n);
	free(b1);
	free(b2);
}


void write_2d_matrix(char *filename, fptype **arr, int n)
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


void write_1d_matrix(char *filename, fptype *arr, int n)
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
		fprintf(outfile, "%10f\n", arr[i]);
	}

	fclose(outfile);
}

#if 1

fptype  **cholesky_decomposition(fptype **a, int n)
{
	int i, j, k;
	fptype sum, **l = alloc_2d_matrix(n);

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

fptype  **cholesky_decomposition(fptype **a, int n)
{
	int i, j, k;
	fptype sum, **l = alloc_2d_matrix(n);

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
