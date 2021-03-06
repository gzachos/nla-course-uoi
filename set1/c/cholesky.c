/*
 * +-----------------------------------------------------------------------+
 * |               Copyright (C) 2020 George Z. Zachos                     |
 * +-----------------------------------------------------------------------+
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * Contact Information:
 * Name: George Z. Zachos
 * Email: gzzachos <at> gmail.com
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <errno.h>
#include <math.h>

#define  BUFF_SIZE	16

// #define  PRINT_INPUT_MATRICES
// #define  VERIFY_CHOLESKY_DECOMP
// #define  PRINT_INTERMEDIATE_RESULTS
#define  PRINT_RESULTS
// #define  PRINT_TOFILE

#if 1
	#define FPTYPE_DOUBLE
	typedef double fptype;
#else
	#undef FPTYPE_DOUBLE
	typedef float fptype;
#endif

typedef enum {S1=1, S2} sys_id;

/* Function Prototypes */
int       alloc_2d_matrices(int n, int num_args, ...);
int       alloc_1d_matrices(int n, int num_args, ...);
fptype  **alloc_2d_matrix(int n);
fptype   *alloc_1d_matrix(int n);
void      init_matrices(fptype **a1, fptype **a2, fptype *b1, fptype *b2, int n);
void      write_2d_matrix(char *filename, fptype **mat, int n);
void      write_1d_matrix(char *filename, fptype *mat, int n);
void      print_input_matrices(void);
void      solve_system(fptype **a, fptype *b, int n, sys_id sid);
fptype  **cholesky_decomposition(fptype **a, int n);
void      verify_cholesky_decomposition(fptype **l, int n, sys_id sid);
fptype  **transpose(fptype **mat, int n);
fptype   *forward_substitution(fptype **l, fptype *b, int n);
fptype   *back_substitution(fptype **l, fptype *b, int n);
void      free_1d_matrices(int num_args, ...);
void      free_2d_matrices(int n, int num_args, ...);
void      free_2d_matrix(fptype **mat, int n);


int main(int argc, char **argv)
{
	int n;
	fptype **a1  = NULL, **a2  = NULL,
	        *b1  = NULL,  *b2  = NULL;
#ifdef PRINT_INPUT_MATRICES
	char filename[BUFF_SIZE];
#endif

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

	if (alloc_1d_matrices(n, 2, &b1, &b2) != 0 ||
			alloc_2d_matrices(n, 2, &a1, &a2))
		return EXIT_FAILURE;

	/* Initialize matrices */
	init_matrices(a1, a2, b1, b2, n);

	/* Print a1, a2, b1 and b2 to files or stdout */
#ifdef PRINT_INPUT_MATRICES
	print_input_matrices();
#endif

	solve_system(a1, b1, n, S1);
	solve_system(a2, b2, n, S2);

	free_2d_matrices(n, 2, a1, a2);
	free_1d_matrices(2, b1, b2);

	return EXIT_SUCCESS;
}


int alloc_1d_matrices(int n, int num_args, ...)
{
	va_list args;
	fptype ***ptr;
	int i, j;

	if (!(ptr = (fptype ***) malloc(num_args * sizeof(fptype **))))
	{
		perror("malloc");
		return EXIT_FAILURE;
	}

	va_start(args, num_args);
	for (i = 0; i < num_args; i++)
	{
		ptr[i] = va_arg(args, fptype **);
		if (!(*ptr[i] = alloc_1d_matrix(n)))
		{
			for (j = 0; j < i; j++)
				free(*ptr[j]);
			free(ptr);
			va_end(args);
			return EXIT_FAILURE;
		}
	}
	free(ptr);
	va_end(args);
	return EXIT_SUCCESS;
}


int alloc_2d_matrices(int n, int num_args, ...)
{
	va_list args;
	fptype ****ptr;
	int i, j;

	if (!(ptr = (fptype ****) malloc(num_args * sizeof(fptype ***))))
	{
		perror("malloc");
		return EXIT_FAILURE;
	}

	va_start(args, num_args);
	for (i = 0; i < num_args; i++)
	{
		ptr[i] = va_arg(args, fptype ***);
		if (!(*ptr[i] = alloc_2d_matrix(n)))
		{
			for (j = 0; j < i; j++)
				free_2d_matrix(*ptr[j], n);
			free(ptr);
			va_end(args);
			return EXIT_FAILURE;
		}
	}
	free(ptr);
	va_end(args);
	return EXIT_SUCCESS;
}


fptype **alloc_2d_matrix(int n)
{
	int i, j;
	fptype **mat = (fptype **) malloc(n * sizeof(fptype *));

	if (!mat)
	{
		perror("malloc");
		return NULL;
	}

	for (i = 0; i < n; i++)
	{
		mat[i] = (fptype *) calloc(n, sizeof(fptype));
		if (!mat[i])
		{
			perror("calloc");
			for (j = 0; j < i; j++)
				free(mat[j]);
			free(mat);
			return NULL;
		}
	}
	return mat;
}


fptype *alloc_1d_matrix(int n)
{
	fptype *mat = (fptype *) calloc(n, sizeof(fptype));

	if (!mat)
	{
		perror("calloc");
		return NULL;
	}

	return mat;
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
			else if (abs(i-j) == 1)
				a1[i][j] = a2[i][j] = -4;
			else if (abs(i-j) == 2)
				a1[i][j] = a2[i][j] = 1;
		}
	}
}


void write_2d_matrix(char *filename, fptype **mat, int n)
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
			fprintf(outfile, "%10f  ", mat[i][j]);
		fprintf(outfile, "\n");
	}

#ifdef PRINT_TOFILE
	fclose(outfile);
#endif
}


void write_1d_matrix(char *filename, fptype *mat, int n)
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
		fprintf(outfile, "%10f\n", mat[i]);

#ifdef PRINT_TOFILE
	fclose(outfile);
#endif
}


#ifdef PRINT_INPUT_MATRICES
void print_input_matrices(void)
{
	#ifndef PRINT_TOFILE
	printf("\nWriting A1...\n");
	#endif
	snprintf(filename, BUFF_SIZE, "a1_%d.txt", n);
	write_2d_matrix(filename, a1, n);

	#ifndef PRINT_TOFILE
	printf("\nWriting A2...\n");
	#endif
	snprintf(filename, BUFF_SIZE, "a2_%d.txt", n);
	write_2d_matrix(filename, a2, n);

	#ifndef PRINT_TOFILE
	printf("\nWriting B1...\n");
	#endif
	snprintf(filename, BUFF_SIZE, "b1_%d.txt", n);
	write_1d_matrix(filename, b1, n);

	#ifndef PRINT_TOFILE
	printf("\nWriting B2...\n");
	#endif
	snprintf(filename, BUFF_SIZE, "b2_%d.txt", n);
	write_1d_matrix(filename, b2, n);
}
#endif


void solve_system(fptype **a, fptype *b, int n, sys_id sid)
{
	fptype **l, **l_trn, *y, *x;
#if defined(PRINT_INTERMEDIATE_RESULTS) || defined(PRINT_RESULTS)
	char filename[BUFF_SIZE];
#endif

	if (!(l = cholesky_decomposition(a, n)))
		return;

#ifdef VERIFY_CHOLESKY_DECOMP
	verify_cholesky_decomposition(l, n, sid);
#endif

	/* L * y = b, solve for y */
	if (!(y = forward_substitution(l, b, n)))
	{
		free_2d_matrix(l, n);
		return;
	}

	/* Calculate transpose of matrix l */
	if (!(l_trn = transpose(l, n)))
	{
		free_2d_matrix(l, n);
		free(y);
		return;
	}

	/* L^T * x = y, solve for x */
	if (!(x = back_substitution(l_trn, y, n)))
	{
		free_2d_matrices(n, 2, l, l_trn);
		free(y);
		return;
	}

#ifdef PRINT_INTERMEDIATE_RESULTS
	#ifndef PRINT_TOFILE
	printf("\nWriting L%d...\n", sid);
	#endif
	snprintf(filename, BUFF_SIZE, "l%1d_%d.txt", sid, n);
	write_2d_matrix(filename, l, n);

	#ifndef PRINT_TOFILE
	printf("\nWriting L%d^T...\n", sid);
	#endif
	snprintf(filename, BUFF_SIZE, "lt%1d_%d.txt", sid, n);
	write_2d_matrix(filename, l_trn, n);

	#ifndef PRINT_TOFILE
	printf("\nWriting Y%d...\n", sid);
	#endif
	snprintf(filename, BUFF_SIZE, "y%1d_%d.txt", sid, n);
	write_1d_matrix(filename, y, n);
#endif

#ifdef PRINT_RESULTS
	#ifndef PRINT_TOFILE
	printf("\nWriting X%d...\n", sid);
	#endif
	snprintf(filename, BUFF_SIZE, "x%1d_%d.txt", sid, n);
	write_1d_matrix(filename, x, n);
#endif

	free_2d_matrices(n, 2, l, l_trn);
	free_1d_matrices(2, y, x);
}


fptype  **cholesky_decomposition(fptype **a, int n)
{
	int i, j, k;
	fptype sum, **l;

	if (!(l = alloc_2d_matrix(n)))
		return NULL;

	for (i = 1; i <= n; i++)
	{
#ifdef OPTIMIZED
		for (j = (i-2 >= 1) ? i-2 : 1; j <= i-1; j++)
#else
		for (j = 1; j <= i-1; j++)
#endif
		{
			sum = 0.0;
#ifdef OPTIMIZED
			for (k = (j-2 >= 1) ? j-2 : 1; k <= j-1; k++)
#else
			for (k = 1; k <= j-1; k++)
#endif
			{
				sum += l[i-1][k-1] * l[j-1][k-1];
			}
			l[i-1][j-1] = (a[i-1][j-1] - sum) / l[j-1][j-1];
		}
		sum = 0.0;
#ifdef OPTIMIZED
		for (j = (i-2 >= 1) ? i-2 : 1; j <= i-1; j++)
#else
		for (j = 1; j <= i-1; j++)
#endif
		{
			sum += l[i-1][j-1] * l[i-1][j-1];
		}
#ifdef FPTYPE_DOUBLE
		l[i-1][i-1] = sqrt(a[i-1][i-1] - sum);
#else
		l[i-1][i-1] = sqrtf(a[i-1][i-1] - sum);
#endif
	}
	return l;
}


#ifdef VERIFY_CHOLESKY_DECOMP
void verify_cholesky_decomposition(fptype **l, int n, sys_id sid)
{
	int i, j, k;
	fptype **ra, sum;
	char filename[BUFF_SIZE];

	if (!(ra = alloc_2d_matrix(n)))
	{
		printf("[WARNING]: Couldn't verify cholesky decomposition\n");
		return;
	}

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

	#ifndef PRINT_TOFILE
	printf("\nWriting A%d = L%d * L%d^T...\n", sid, sid, sid);
	#endif
	snprintf(filename, BUFF_SIZE, "ra%1d_%d.txt", sid , n);
	write_2d_matrix(filename, ra, n);
	free_2d_matrix(ra, n);
}
#endif


fptype **transpose(fptype **mat, int n)
{
	int i, j;
	fptype **trn;

	if (!(trn = alloc_2d_matrix(n)))
		return NULL;

#ifdef OPTIMIZED
	for (i = 0; i < n; i++)
		for (j = (i-2 >= 0) ? i-2 : 0; j < n; j++)
			trn[i][j] = mat[j][i];
#else
	for (i = 0; i < n; i++)
		for (j = 0; j < n; j++)
			trn[i][j] = mat[j][i];
#endif
	return trn;
}


fptype *forward_substitution(fptype **l, fptype *b, int n)
{
	int i, k;
	fptype *y, sum;

	if (!(y = alloc_1d_matrix(n)))
		return NULL;

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
	int i, k;
	fptype *y, sum;

	if (!(y = alloc_1d_matrix(n)))
		return NULL;

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


void free_1d_matrices(int num_args, ...)
{
	va_list args;
	fptype *ptr;
	int i;

	va_start(args, num_args);
	for (i = 0; i < num_args; i++)
	{
		ptr = va_arg(args, fptype *);
		free(ptr);
	}
	va_end(args);
}


void free_2d_matrices(int n, int num_args, ...)
{
	va_list args;
	fptype **ptr;
	int i;

	va_start(args, num_args);
	for (i = 0; i < num_args; i++)
	{
		ptr = va_arg(args, fptype **);
		free_2d_matrix(ptr, n);
	}
	va_end(args);
}


void free_2d_matrix(fptype **mat, int n)
{
	int i;

	if (!mat)
		return;

	for (i = 0; i < n; i++)
		free(mat[i]);
	free(mat);
}

