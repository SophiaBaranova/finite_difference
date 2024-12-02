#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>
#include "error.h"
#include "color.h"

double px(double x);
double qx(double x);
double fx(double x);
double* finite_difference(float A[], float B[], int n, double x0, double xn, double *x);

int main()
{
	double x0, xn; //start and end points of the interval
	float A[3], B[3]; //coefficients in boundary condition equations
	double h, e = 1.; //length of subintervals and accuracy
    double r2 = 10., r1, delta; //calculation errors
	int i, n, e_out, iter = 0; //counters
	double *x1 = NULL, *y1 = NULL, *x2 = NULL, *y2 = NULL;
	FILE *p = NULL;
	int ans, stop = 0;
	printf("FINITE DIFFERENCE METHOD FOR THE SECOND ORDER ORDINARY DIFFERENTIAL EQUATION WITH A LINEAR BOUNDARY CONDITION\n");
	yellow();
	printf("y'' - y'/4 + 2y/x = x/2\n");
	reset();
	printf("boundary condition:\n1 - own\n2 - by default\n-> ");
	scanf("%d", &ans);
	if (ans == 1)
	{
		//enter boundary condition
		printf("x0: ");
		scanf("%lf", &x0);
		printf("xn: ");
		scanf("%lf", &xn);
		printf("alpha0: ");
		scanf("%f", &A[0]);
		printf("alpha1: ");
		scanf("%f", &A[1]);
		printf("A: ");
		scanf("%f", &A[2]);
		printf("beta0: ");
		scanf("%f", &B[0]);
		printf("beta1: ");
		scanf("%f", &B[1]);
		printf("B: ");
		scanf("%f", &B[2]);
	}
	else if (ans == 2)
	{
		//read boundary condition from input file
		p = fopen("input.txt", "r");
		if (!p)
		{
			error(1);
		}
		fscanf(p, "%f%f%f", &A[0], &A[1], &A[2]);
		fscanf(p, "%f%f%f", &B[0], &B[1], &B[2]);
		fscanf(p, "%lf%lf", &x0, &xn);
		fclose(p);
	}
    else
    {
		printf("invalid answer\n");
        return -1;
    }
	//print boundary condition
	printf("\nboundary condition:\n");
	yellow();
	printf("%.1fy(%.1lf) ", A[0], x0);
	if (A[1] > 0)
	{
		printf("+ %.1fy''(%.1lf) ", A[1], x0);
	}
	else if (A[1] < 0)
	{
		printf("- %.1fy''(%.1lf) ", fabs(A[1]), x0);
	}
	printf("= %.1f\n", A[2]);
	printf("%.1fy(%.1lf) ", B[0], xn);
	if (B[1] > 0)
	{
		printf("+ %.1fy''(%.1lf) ", B[1], xn);
	}
	else if (B[1] < 0)
	{
		printf("- %.1fy''(%.1lf) ", fabs(B[1]), xn);
	}
	printf("= %.1f\n", B[2]);
	reset();
	//enter initial length of subintervals and desired accuracy
	do
	{
		printf("\ninitial length of subintervals: ");
		scanf("%lf", &h);
	} while (h <= 0 || h > (xn - x0));
	printf("desired accuracy (number of digits after the decimal point): ");
	scanf("%d", &e_out);
	//calculate accuracy
	for (i = 0; i < e_out; i++)
	{
		e /= 10;
	}
	//calculate number of subintervals
	n = (xn - x0) / h;
	x2 = (double*)malloc((n + 1) * sizeof(double));
	if (!x2)
	{
		error(2);
	}
	//calculate x
	for (i = 0; i <= n; i++)
	{
		x2[i] = x0 + i * h;
	}
	//calculate y using finite difference method
	y2 = finite_difference(A, B, n, x0, xn, x2);
	//while the desired accuracy isn't reached
	while (1)
	{
        iter++;
		x1 = x2;
		y1 = y2;
		r1 = r2;
		//double number of subintervals
		n *= 2;
		//calculate length of subintervals
		h = (xn - x0) / n;
		x2 = (double*)malloc((n + 1) * sizeof(double));
		if (!x2)
		{
			error(2);
		}
		//calculate x
		for (i = 0; i <= n; i++)
		{
			x2[i] = x0 + i * h;
		}
		//calculate y using finite difference method
		y2 = finite_difference(A, B, n, x0, xn, x2);
		//calculate error using runge's rule
		r2 = fabs(y1[0] - y2[0]);
		for (i = 2; i <= n; i += 2)
		{
			delta = fabs(y1[i / 2] - y2[i]);
			//memorize max error for this iteration
			if (r2 < delta)
			{
				r2 = delta;
			}
		}
		//if calculation error stopped decreasing
		if (r2 >= r1)
		{
		    stop = 1;
		 	break;
		}
		//if calculation error is greater than desired accuracy
        if (r2 > e)
        {
            free(x1);
            free(y1);
            continue;
        }
		//if the desired accuracy is reached
        else
        {
            break;
        }
	}
	//if the desired accuracy wasn't reached
	if (stop)
	{
		red();
		printf("\nthe desired accuracy wasn't reached\n");
		reset();
	}
	//print results
	yellow();
	printf("\nsolution:\n");
	reset();
	printf("%10s|%10s\n", "x", "y");
	for (i = 0; i <= n / 2; i++)
	{
		printf("%10.*lf|%10.*lf\n", e_out, x1[i], e_out, y1[i]);
	}
	yellow();
	printf("\ncalculation error (using Runge's rule): %.4le\n", r2);
	printf("number of iterations: %d\n", iter);
    printf("number of subintervals: %d\n", n / 2);
	reset();
	//free the memory
	free(x1);
	free(y1);
	free(x2);
	free(y2);
	return 0;
}

double px(double x)
{
	return -1. / 4.;
}

double qx(double x)
{
	if (x == 0)
	{
		error(3);
	}
	return 2. / x;
}

double fx(double x)
{
	return x / 2.;
}

double* finite_difference(float A[], float B[], int n, double x0, double xn, double *x)
{
	double *y = NULL;
	double *a = NULL, *b = NULL, *c = NULL, *d = NULL;
	double *s = NULL, *t = NULL;
	double h = (xn - x0) / n;
	int i;
	a = (double*)malloc((n + 1) * sizeof(double));
	b = (double*)malloc((n + 1) * sizeof(double));
	c = (double*)malloc((n + 1) * sizeof(double));
	d = (double*)malloc((n + 1) * sizeof(double));
	s = (double*)malloc(n * sizeof(double));
	t = (double*)malloc((n + 1) * sizeof(double));
	y = (double*)malloc((n + 1) * sizeof(double));
	if (!a || !b || !c || !d || !s || !t || !y)
	{
		error(2);
	}
	//calculate coefficients of the tri-diagonal system of linear equations
	b[0] = h * A[0] - A[1];
	c[0] = A[1];
	d[0] = h * A[2];
	for (i = 1; i < n; i++)
	{
		a[i] = 1 - h * px(x[i]) / 2;
		b[i] = pow(h, 2) * qx(x[i]) - 2;
		c[i] = 1 + h * px(x[i]) / 2;
		d[i] = pow(h, 2) * fx(x[i]);
	}
	a[n] = -B[1];
	b[n] = h * B[0] + B[1];
	d[n] = h * B[2];
	//calculate coefficients for double sweep method
	s[0] = -c[0] / b[0];
	t[0] = d[0] / b[0];
	for (i = 1; i < n; i++)
	{
		s[i] = -c[i] / (b[i] + a[i] * s[i - 1]);
		t[i] = (d[i] - a[i] * t[i - 1]) / (b[i] + a[i] * s[i - 1]);
	}
	t[n] = (d[n] - a[n] * t[n - 1]) / (b[n] + a[n] * s[n - 1]);
	//find unknown y in the reverse order
	y[n] = t[n];
	for (i = n - 1; i >= 0; i--)
	{
		y[i] = s[i] * y[i + 1] + t[i];
	}
	//free the memory
	free(a);
	free(b);
	free(c);
	free(d);
	free(s);
	free(t);
	return y;
}
