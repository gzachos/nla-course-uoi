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

#define BUFF_SIZE      16
#define NUM_METHODS    2
#define MAX_ERROR      0.00005

// #define  PRINT_INPUT_MATRICES
#define  PRINT_RESULTS
// #define  PRINT_TOFILE

#if 1
	#define FPTYPE_DOUBLE
	typedef double fptype;
#else
	#undef FPTYPE_DOUBLE
	typedef float fptype;
#endif

typedef enum {A1, A2, B1, ALL} dealloc_level;
typedef enum {S1=1, S2} sys_id;
typedef enum {STEEPEST_DESCENT, CONJUGATE_GRADIENT} method_id;

/* Function Prototypes */
int       alloc_matrices(fptype ***a1, fptype ***a2, fptype **b1, fptype **b2,
		int n);
int       alloc_sd_matrices(fptype **x, fptype **r, fptype **Ar, fptype **Ax,
		fptype **tmp, int n);
int       alloc_cg_matrices(fptype **x, fptype **r, fptype **p, fptype **Ap,
		fptype **Ax, fptype **tmp, int n);
fptype  **alloc_2d_matrix(int n);
fptype   *alloc_1d_matrix(int n);
void      init_matrices(fptype **a1, fptype **a2, fptype *b1, fptype *b2, int n);
void      write_2d_matrix(char *filename, fptype **mat, int n);
void      write_1d_matrix(char *filename, fptype *mat, int n);
void      print_input_matrices(void);
void      solve_system(fptype **a, fptype *b, int n, sys_id sid);
fptype   *steepest_descent(fptype **A, fptype *b, fptype max_error, int n);
fptype   *conjugate_gradient(fptype **A, fptype *b, fptype max_error, int n);
fptype    euclidean_norm(fptype *v, int n);
fptype    dot_product(fptype *v1, fptype *v2, int n);
fptype   *matrix_vector_multiplication(fptype *res, fptype **mat, fptype *v,
		int n);
fptype   *scalar_vector_multiplication(fptype *res, fptype s, fptype *v, int n);
fptype   *add_vectors(fptype *res, fptype *v1, fptype *v2, int n);
fptype   *subtract_vectors(fptype *res, fptype *v1, fptype *v2, int n);
void      free_matrices(fptype **a1, fptype **a2, fptype *b1, fptype *b2,
		int n, dealloc_level level);
void      free_sd_matrices(fptype *r, fptype *Ar, fptype *Ax, fptype *tmp);
void      free_cg_matrices(fptype **r, fptype *p, fptype *Ap, fptype *Ax,
		fptype *tmp);
void      free_2d_matrix(fptype **mat, int n);

/* Global data */
fptype *(*methods[NUM_METHODS])(fptype **A, fptype *b, fptype max_error, int n) = {
	steepest_descent,
	conjugate_gradient
};

char *method_names[NUM_METHODS] = {"Steepest Descent", "Conjugate Gradient"};
char *method_initials[NUM_METHODS] = {"sd", "cg"};


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

int alloc_sd_matrices(fptype **x, fptype **r, fptype **Ar, fptype **Ax, fptype **tmp, int n)
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

	if (!(*tmp = alloc_1d_matrix(n)))
	{
		free(*x);
		free(*x);
		free(*Ar);
		free(*Ax);
		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}


int alloc_cg_matrices(fptype **x, fptype **r, fptype **p, fptype **Ap,
		fptype **Ax, fptype **tmp, int n)
{
	int i, j;

	if (alloc_sd_matrices(x, p, Ap, Ax, tmp, n) != 0)
		return EXIT_FAILURE;

	for (i = 0; i < 3; i++)
	{
		if (!(r[i] = alloc_1d_matrix(n)))
		{
			free(*x);
			free(*p);
			free(*Ap);
			free(*Ax);
			free(*tmp);
			for (j = 0; j < i; j++)
				free(r[j]);
			return EXIT_FAILURE;
		}
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
	int i;
	fptype *x;
#ifdef PRINT_RESULTS
	char filename[BUFF_SIZE];
#endif
	printf("\n######################\n");
	printf("# System: No.%1d       #\n", sid);
	printf("######################\n");
	for (i = 0; i < NUM_METHODS; i++)
	{
		printf("\n# Method: %s\n", method_names[i]);
		if (!(x = (methods[i])(a, b, MAX_ERROR, n)))
			return;
#ifdef PRINT_RESULTS
	#ifndef PRINT_TOFILE
		printf("\nWriting X%d...\n", sid);
	#endif
		snprintf(filename, BUFF_SIZE, "x%1d_%d_%2s.txt", sid, n, method_initials[i]);
		write_1d_matrix(filename, x, n);
#endif
		free(x);
	}
}


fptype *steepest_descent(fptype **A, fptype *b, fptype max_error, int n)
{
	int k;
	fptype *x, *r, *Ar, *Ax, *tmp, a;

	if (alloc_sd_matrices(&x, &r, &Ar, &Ax, &tmp, n) != 0)
		return NULL;

	// Vector x^(0) is already the zero vector
	memcpy(r, b, n*sizeof(fptype)); // r^(0) = b;
	k = 0;
	while (euclidean_norm(r, n) > max_error)
	{
		k++;
		// Ar = A * r^(k-1)
		Ar = matrix_vector_multiplication(Ar, A, r, n);
		// a_k = (r^(k-1), r^(k-1)) / (Ar, r^(k-1))
		a  = dot_product(r, r, n) / dot_product(Ar, r, n);
		// x^(k) = x^(k-1) + a_k * r^(k-1)
		x  = add_vectors(x, x, scalar_vector_multiplication(tmp, a, r, n), n);
		// Ax = A * x^(k)
		Ax = matrix_vector_multiplication(Ax, A, x, n);
#ifndef OPTIMIZED
		// r^(k) = b - Ax
		r  = subtract_vectors(r, b, Ax, n);
#else
		// r^(k) = r^(k-1) - a_k * Ar
		r  = subtract_vectors(r, r, scalar_vector_multiplication(Ar, a, Ar, n), n);
#endif
	}

	printf("\nk = %d\n", k);

	free_sd_matrices(r, Ar, Ax, tmp);

	return x;
}


fptype *conjugate_gradient(fptype **A, fptype *b, fptype max_error, int n)
{
	int k;
	fptype *x, *r[3], *p, a_k, b_k, *Ap, *Ax, *tmp, *tmp_ptr;

	if (alloc_cg_matrices(&x, r, &p, &Ap, &Ax, &tmp, n) != 0)
		return NULL;

	// Vector x^(0) is already the zero vector
	memcpy(r[0], b, n*sizeof(fptype)); // r^(0) = b;
	memcpy(p, r[0], n*sizeof(fptype)); // p^(0) = r^(0);
	// Ap = A * p^(1)
	Ap = matrix_vector_multiplication(Ap, A, p, n);
	// a_1 = (r^(0), r^(0)) / (Ap, p^(1))
	a_k = dot_product(r[0], r[0], n) / dot_product(Ap, p, n);
	// x^(1) = x^(0) + a_1 * p^(1)
	x = add_vectors(x, x, scalar_vector_multiplication(tmp, a_k, p, n), n);
	// Ax = A * x^(1)
	Ax = matrix_vector_multiplication(Ax, A, x, n);
	// r^(1) = b - Ax
	r[1] = subtract_vectors(r[1], b, Ax, n);
	k = 1;
	while (euclidean_norm(r[1], n) > max_error && k < n)
	{
		k++;
		// b_k = (r^(k-1), r^(k-1)) / (r^(k-2), r^(k-2))
		b_k = dot_product(r[1], r[1], n) / dot_product(r[0], r[0], n);
		// p^(k) = r^(k-1) + b_k * p^(k-1)
		p = add_vectors(p, r[1], scalar_vector_multiplication(tmp, b_k, p, n), n);
		// Ap = A * p^(k)
		Ap = matrix_vector_multiplication(Ap, A, p, n);
		// a_k = (r^(k-1), r^(k-1)) / (Ap, p^(k))
		a_k = dot_product(r[1], r[1], n) / dot_product(Ap, p, n);
		// x^(k) = x^(k-1) + a_k * p^(k)
		x = add_vectors(x, x, scalar_vector_multiplication(tmp, a_k, p, n), n);
		// Ax = A * x^(k)
		Ax = matrix_vector_multiplication(Ax, A, x, n);
#ifndef OPTIMIZED
		// r^(k) = b - Ax
		r[2] = subtract_vectors(r[2], b, Ax, n);
#else
		// r^(k) = r^(k-1) - a_k * Ap
		r[2] = subtract_vectors(r[2], r[1], scalar_vector_multiplication(Ap, a_k, Ap, n), n);
#endif
		// Shift left r[i] in a round-robin manner
		tmp_ptr = r[0];
		r[0] = r[1];
		r[1] = r[2];
		r[2] = tmp_ptr;
	}

	printf("\nk = %d\n", k);

	free_cg_matrices(r, p, Ap, Ax, tmp);

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
	int i, j, lb, ub;

	if (!res || !mat || !v)
		return NULL;

	for (i = 0; i < n; i++)
	{
#ifdef OPTIMIZED
		lb = (i >= 2)   ? i-2 : 0;
		ub = (i <= n-3) ? i+2 : n-1;
#else
		lb = 0;
		ub = n-1;
#endif
		for (j = lb, res[i] = 0.0; j <= ub; j++)
			res[i] += mat[i][j] * v[j];
	}
	return res;
}


fptype *scalar_vector_multiplication(fptype *res, fptype s, fptype *v, int n)
{
	int i;

	if (!res || !v)
		return NULL;

	for (i = 0; i < n; i++)
		res[i] = v[i] * s;

	return res;
}


fptype *add_vectors(fptype *res, fptype *v1, fptype *v2, int n)
{
	int i;

	if (!res || !v1 || !v2)
		return NULL;

	for (i = 0; i < n; i++)
		res[i] = v1[i] + v2[i];

	return res;
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


void free_sd_matrices(fptype *r, fptype *Ar, fptype *Ax, fptype *tmp)
{
	free(r);
	free(Ar);
	free(Ax);
	free(tmp);
}


void free_cg_matrices(fptype **r, fptype *p, fptype *Ap, fptype *Ax, fptype *tmp)
{
	free(r[0]);
	free(r[1]);
	free(r[2]);
	free_sd_matrices(p, Ap, Ax, tmp);
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

