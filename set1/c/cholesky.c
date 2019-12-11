#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <math.h>

#define BUFF_SIZE	16

#define PRINT_INPUT_MATRICES
#define VERIFY_CHOLESKY_DECOMP
#define PRINT_INTERMEDIATE_RESULTS
#define PRINT_RESULTS
#define PRINT_TOFILE

#undef PRINT_INPUT_MATRICED
// #undef VERIFY_CHOLESKY_DECOMP
#undef PRINT_INTERMEDIARE_RESULTS
// #undef PRINT_RESULTS
// #undef PRINT_TOFILE

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
void      solve_system(fptype **a, fptype *b, int n, int sysn);
void      verify_cholesky_decomposition(fptype **l, int n, int sysn);


int main(int argc, char **argv)
{
	int n;
	fptype **a1  = NULL, **a2  = NULL,
	        *b1  = NULL,  *b2  = NULL;
	char filename[BUFF_SIZE];

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
#ifdef PRINT_INPUT_MATRICES
	snprintf(filename, BUFF_SIZE, "a1_%d.txt", n);
	write_2d_matrix(filename, a1, n);
	snprintf(filename, BUFF_SIZE, "a2_%d.txt", n);
	write_2d_matrix(filename, a2, n);
	snprintf(filename, BUFF_SIZE, "b1_%d.txt", n);
	write_1d_matrix(filename, b1, n);
	snprintf(filename, BUFF_SIZE, "b2_%d.txt", n);
	write_1d_matrix(filename, b2, n);
#endif

	solve_system(a1, b1, n, 1);
	solve_system(a2, b2, n, 2);

	free_matrices(a1, a2, b1, b2, n);

	return EXIT_SUCCESS;
}


#ifdef VERIFY_CHOLESKY_DECOMP
void verify_cholesky_decomposition(fptype **l, int n, int sysn)
{
	int i, j, k;
	fptype **ra, sum;
	char filename[BUFF_SIZE];

	ra = alloc_2d_matrix(n);

	/* Calculate L * L-transpose */
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < n; j++)
		{
			for (k = 0, sum = 0.0; k < n; k++)
				sum += l[i][k] * l[j][k];
			ra[i][j] = sum;
		}
	}

	snprintf(filename, BUFF_SIZE, "ra%1d_%d.txt", sysn , n);
	write_2d_matrix(filename, ra, n);
	free_2d_matrix(ra, n);
}
#endif


void solve_system(fptype **a, fptype *b, int n, int sysn)
{
	fptype **l, **l_trn, *y, *x;
	char filename[BUFF_SIZE];

	l = cholesky_decomposition(a, n);

#ifdef VERIFY_CHOLESKY_DECOMP
	verify_cholesky_decomposition(l, n, sysn);
#endif

	/* L * y = b, solve for y */
	y = forward_substitution(l, b, n);

	/* Calculate transpose of matrix l */
	l_trn = transpose(l, n);

	/* L^T * x = y, solve for x */
	x = back_substitution(l_trn, y, n);

#ifdef PRINT_INTERMEDIATE_RESULTS
	snprintf(filename, BUFF_SIZE, "l%1d_%d.txt", sysn, n);
	write_2d_matrix(filename, l, n);
	snprintf(filename, BUFF_SIZE, "lt%1d_%d.txt", sysn, n);
	write_2d_matrix(filename, l_trn, n);
	snprintf(filename, BUFF_SIZE, "x%1d_%d.txt", sysn, n);
	write_1d_matrix(filename, x, n);
#endif

#ifdef PRINT_RESULTS
	snprintf(filename, BUFF_SIZE, "y%1d_%d.txt", sysn, n);
	write_1d_matrix(filename, y, n);
#endif

	free_2d_matrix(l, n);
	free_2d_matrix(l_trn, n);
	free(y);
	free(x);
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
#ifdef OPTIMIZED
	for (i = 0; i < n; i++)
		for (j = (i-2 >= 0) ? i-2 : 0; j < n; j++)
			trn[i][j] = arr[j][i];
#else
	for (i = 0; i < n; i++)
		for (j = 0; j < n; j++)
			trn[i][j] = arr[j][i];
#endif
	return trn;
}


fptype *forward_substitution(fptype **l, fptype *b, int n)
{
	fptype *y = alloc_1d_matrix(n), sum;
	int i, k;

	for (i = 0; i < n; i++)
	{
		sum = b[i];
#ifdef OPTIMIZED
		for (k = (i-2 >= 0) ? i-2 : 0; k < i; k++)
#else
		for (k = 0; k < i; k++)
#endif
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
#ifdef OPTIMIZED
		for (k = (i+2 <= n) ? i+2 : n; k >= i+1; k--)
#else
		for (k = n; k >= i+1; k--)
#endif
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
	FILE *outfile;

#ifdef PRINT_TOFILE
	outfile = fopen(filename, "w");
	if (!outfile)
	{
		perror("fopen");
		return;
	}
#else
	outfile = stdout;
#endif

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
	FILE *outfile;

#ifdef PRINT_TOFILE
	outfile = fopen(filename, "w");
	if (!outfile)
	{
		perror("fopen");
		return;
	}
#else
	outfile = stdout;
#endif

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
#ifdef OPTIMIZED
		for (j = (i-2 >= 1) ? i-2 : 1; j <= i-1; j++)
#else
		for (j = 1; j <= i-1; j++)
#endif
		{
			sum = a[i-1][j-1];
#ifdef OPTIMIZED
			for (k = (j-2 >= 1) ? j-2 : 1; k <= j-1; k++)
#else
			for (k = 1; k <= j-1; k++)
#endif
			{
				sum -= l[i-1][k-1] * l[j-1][k-1];
			}
			l[i-1][j-1] = sum / l[j-1][j-1];
		}
		sum = a[i-1][i-1];
#ifdef OPTIMIZED
		for (j = (i-2 >= 1) ? i-2 : 1; j <= i-1; j++)
#else
		for (j = 1; j <= i-1; j++)
#endif
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
