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
#include <string.h>
#include <errno.h>
#include <math.h>

#define  BUFF_SIZE	16

// #define  PRINT_INPUT_MATRICES
#define  PRINT_RESULTS
#define  PRINT_TOFILE

#if 1
	#define FPTYPE_DOUBLE
	typedef double fptype;
#else
	#undef FPTYPE_DOUBLE
	typedef float fptype;
#endif

typedef enum {A1, A2, B1, ALL} dealloc_level;
typedef enum {S1=1, S2} sys_id;

/* Function Prototypes */
int       alloc_matrices(fptype ***a1, fptype ***a2, fptype **b1, fptype **b2,
		int n);
int       alloc_sd_matrices(fptype **x, fptype **r, fptype **Ar, fptype **Ax,
		int n);
fptype  **alloc_2d_matrix(int n);
fptype   *alloc_1d_matrix(int n);
void      init_matrices(fptype **a1, fptype **a2, fptype *b1, fptype *b2, int n);
void      write_2d_matrix(char *filename, fptype **mat, int n);
void      write_1d_matrix(char *filename, fptype *mat, int n);
void      print_input_matrices(void);
void      solve_system(fptype **a, fptype *b, int n, sys_id sid);
fptype   *steepest_descent(fptype **A, fptype *b, fptype max_error, int n);
fptype    euclidean_norm(fptype *v, int n);
fptype    dot_product(fptype *v1, fptype *v2, int n);
fptype   *matrix_vector_multiplication(fptype *res, fptype **mat, fptype *v,
		int n);
fptype   *scalar_vector_multiplication(fptype s, fptype *v, int n);
fptype   *add_vectors(fptype *v1, fptype *v2, int n);
fptype   *subtract_vectors(fptype *res, fptype *v1, fptype *v2, int n);
void      free_matrices(fptype **a1, fptype **a2, fptype *b1, fptype *b2,
		int n, dealloc_level level);
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

	if (alloc_matrices(&a1, &a2, &b1, &b2, n) != 0)
		return EXIT_FAILURE;

	/* Initialize matrices */
	init_matrices(a1, a2, b1, b2, n);

	/* Print a1, a2, b1 and b2 to files or stdout */
#ifdef PRINT_INPUT_MATRICES
	print_input_matrices();
#endif

	solve_system(a1, b1, n, S1);
	solve_system(a2, b2, n, S2);

	free_matrices(a1, a2, b1, b2, n, ALL);

	return EXIT_SUCCESS;
}


int alloc_matrices(fptype ***a1, fptype ***a2, fptype **b1, fptype **b2, int n)
{
	if (!(*a1 = alloc_2d_matrix(n)))
		return EXIT_FAILURE;

	if (!(*a2 = alloc_2d_matrix(n)))
	{
		free_matrices(*a1, *a2, *b1, *b2, n, A1);
		return EXIT_FAILURE;
	}

	if (!(*b1 = alloc_1d_matrix(n)))
	{
		free_matrices(*a1, *a2, *b1, *b2, n, A2);
		return EXIT_FAILURE;
	}

	if (!(*b2 = alloc_1d_matrix(n)))
	{
		free_matrices(*a1, *a2, *b1, *b2, n, B1);
		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}


int alloc_sd_matrices(fptype **x, fptype **r, fptype **Ar, fptype **Ax, int n)
{
	if (!(*x = alloc_1d_matrix(n)))
		return EXIT_FAILURE;

	if (!(*r = alloc_1d_matrix(n)))
	{
		free(*x);
		return EXIT_FAILURE;
	}

	if (!(*Ar = alloc_1d_matrix(n)))
	{
		free(*x);
		free(*r);
		return EXIT_FAILURE;
	}

	if (!(*Ax = alloc_1d_matrix(n)))
	{
		free(*x);
		free(*x);
		free(*Ar);
		return EXIT_FAILURE;
	}

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
			fprintf(outfile, "%8.4f  ", mat[i][j]);
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
		fprintf(outfile, "%8.4f\n", mat[i]);

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
	fptype *x;
#ifdef PRINT_RESULTS
	char filename[BUFF_SIZE];
#endif

	if (!(x = steepest_descent(a, b, 0.0005, n)))
			return;

#ifdef PRINT_RESULTS
	#ifndef PRINT_TOFILE
	printf("\nWriting X%d...\n", sid);
	#endif
	snprintf(filename, BUFF_SIZE, "x%1d_%d.txt", sid, n);
	write_1d_matrix(filename, x, n);
#endif

	free(x);
}


fptype *steepest_descent(fptype **A, fptype *b, fptype max_error, int n)
{
	int k;
	fptype *x, *r, *Ar, *Ax, a;

	if (alloc_sd_matrices(&x, &r, &Ar, &Ax, n) != 0)
		return NULL;

	// Vector x^(0) is already the zero vector
	memcpy(r, b, n*sizeof(fptype)); // r^(0) = b;
	k = 0;
	while (euclidean_norm(r, n) > max_error)
	{
		k++;
		Ar = matrix_vector_multiplication(Ar, A, r, n);
		a  = dot_product(r, r, n) / dot_product(Ar, r, n);
		x  = add_vectors(x, scalar_vector_multiplication(a, r, n), n);
		Ax = matrix_vector_multiplication(Ax, A, x, n);
		r  = subtract_vectors(r, b, Ax, n);
	}
	
	printf("\nk = %d\n", k);

	free(r);
	free(Ar);
	free(Ax);
	return x;
}


fptype euclidean_norm(fptype *v, int n)
{
	int i;
	fptype sum = 0.0;

	if (!v)
	{
		fprintf(stderr, "euclidean_norm: argument is NULL!\n");
		return -1;
	}

	for (i = 0; i < n; i++)
		sum += v[i] * v[i];
#ifdef FPTYPE_DOUBLE
	return sqrt(sum);
#else
	return sqrtf(sum);
#endif
}


fptype dot_product(fptype *v1, fptype *v2, int n)
{
	int i;
	fptype prod = 0.0;

	if (!v1 || !v2)
	{
		fprintf(stderr, "euclidean_norm: argument is NULL!\n");
		// TODO longjmp
		return -1;
	}

	for (i = 0; i < n; i++)
		prod += v1[i] * v2[i];

	return prod;
}


fptype *matrix_vector_multiplication(fptype *res, fptype **mat, fptype *v, int n)
{
	int i, j;
	fptype tmp_sum;

	if (!res || !mat || !v)
		return NULL;

	for (i = 0; i < n; i++)
	{
		for (j = 0, tmp_sum = 0.0; j < n; j++)
			tmp_sum += mat[i][j] * v[j];
		res[i] = tmp_sum;
	}
	return res;
}


/* Multiply vector v by scalar  s. Doesn't alloc a new vector! */
fptype *scalar_vector_multiplication(fptype s, fptype *v, int n)
{
	int i;

	if (!v)
		return NULL;

	for (i = 0; i < n; i++)
		v[i] *= s;

	return v;
}


/* Add vector v2 to v1. Doesn't alloc a new vector! */
fptype *add_vectors(fptype *v1, fptype *v2, int n)
{
	int i;

	if (!v1 || !v2)
		return NULL;

	for (i = 0; i < n; i++)
		v1[i] += v2[i];

	return v1;
}


fptype *subtract_vectors(fptype *res, fptype *v1, fptype *v2, int n)
{
	int i;

	if (!res || !v1 || !v2)
		return NULL;

	for (i = 0; i < n; i++)
		res[i] = v1[i] - v2[i];

	return res;
}


void free_matrices(fptype **a1, fptype **a2, fptype *b1, fptype *b2,
		int n, dealloc_level level)
{
	switch (level)
	{
		case ALL:
			free(b2);
		case B1:
			free(b1);
		case A2:
			free_2d_matrix(a2, n);
		case A1:
			free_2d_matrix(a1, n);
	}
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

