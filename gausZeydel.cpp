/* 

Пример решения краевой задачи с помощью конечно-разностного метода:
(k(x)*u'(x))' - q(x)*u = - f(x),
u(0) = m1,
u(l) = m2,

где k(x) = (1+alpha*x^2); q(x) = exp(alpha*x); m1=0; m2=0; f(x) = 10*x*(l-x)*exp(alpha*x) + 20*(1+alpha*x^2) - 2*alpha*x*(10*l-20*x); 0<x<1

Для решения СЛАУ используется метод Гаусса-Зейделя
Описание метода и пример реализации: https://ru.wikipedia.org/wiki/%D0%9C%D0%B5%D1%82%D0%BE%D0%B4_%D0%93%D0%B0%D1%83%D1%81%D1%81%D0%B0_%E2%80%94_%D0%97%D0%B5%D0%B9%D0%B4%D0%B5%D0%BB%D1%8F_%D1%80%D0%B5%D1%88%D0%B5%D0%BD%D0%B8%D1%8F_%D1%81%D0%B8%D1%81%D1%82%D0%B5%D0%BC%D1%8B_%D0%BB%D0%B8%D0%BD%D0%B5%D0%B9%D0%BD%D1%8B%D1%85_%D1%83%D1%80%D0%B0%D0%B2%D0%BD%D0%B5%D0%BD%D0%B8%D0%B9#.D0.9F.D1.80.D0.B8.D0.BC.D0.B5.D1.80_.D1.80.D0.B5.D0.B0.D0.BB.D0.B8.D0.B7.D0.B0.D1.86.D0.B8.D0.B8_.D0.BD.D0.B0_C.2B.2B

*/

#include "stdafx.h"
#include <math.h>
#include <iostream>
using namespace std;
#define eps 0.026

// Условие окончания
bool converge(double *xk, double *xkp, int n)
{
	double norm = 0;
	for (int i = 0; i < n; i++)
	{
		norm += (xk[i] - xkp[i])*(xk[i] - xkp[i]);
	}
	if (sqrt(norm) >= eps)
		return false;
	return true;
}

//Известное решение задачи u(x) - используется для подсчета невязки
double funcU(double x) {
	double u = 10 * x*(1 - x);
	return u;
}

//Функция k(x)
double funcK(double x, double alp) {
	return (1 + alp*x*x);
}

/*
Ход метода, где:
a[n][n] - Матрица коэффициентов
x[n], p[n] - Текущее и предыдущее решения
b[n] - Столбец правых частей
Все перечисленные массивы вещественные и
должны быть определены в основной программе,
также в массив x[n] следует поместить начальное
приближение столбца решений (например, все нули)
*/
void gausZey(double *x, int n, double **a, double *b) {
	double *p = new double[n];
	int k = 1;
	x[0] = 0.0;
	x[n - 1] = 0.0;
	do
	{
		for (int i = 0; i < n; i++)
			p[i] = x[i];

		for (int i = 1; i < n-1; i++)
		{
			double var = 0;
			for (int j = 0; j < i; j++)
				var += (a[i][j] * x[j]);
			for (int j = i + 1; j < n; j++)
				var += (a[i][j] * p[j]);
			x[i] = (b[i] - var) / a[i][i];
		}
		double norm = 0;
		for (int i = 0; i < n; i++)
		{
			norm += (x[i] - p[i])*(x[i] - p[i]);
		}
		printf("%.3lf\n", sqrt(norm));
	} while (!converge(x, p, n));
}

// Задача записывается в виде СЛАУ
void init(int N, double *a, double *b, double *c, double *f, double *u, double alp, double h, double l)
{
	double xi;
	for (int i = 0; i < N; i++) {
		xi = 0.0 + i*h;
		if (i != 0) {
			a[i] = funcK(xi - 0.5*h, alp) / (h*h);
		}
		else {
			a[i] = 0.0;
		}

		if (i != (N - 1)) {
			b[i] = funcK(xi + 0.5*h, alp) / (h*h);
		}
		else {
			b[i] = 0.0;
		}

		c[i] = -(a[i] + b[i] + exp((double)(alp*xi)));
	}
	for (int i = 0; i < N; i++) {
		xi = 0.0 + i*h;
		f[i] = -(10.0 * xi * (l - xi) * exp((double)(alp * xi)) + 20.0 * (1.0 + alp * xi * xi) - 2.0 * alp * xi * (10.0 * l - 20.0 * xi));
		u[i] = funcU(xi);
	}
}

int main()
{
	double *a, *b, *c;
	double l = 1.0;
	double h = 0.05;
	double alp = 9;
	double *y, *f, *u;
	double xi;
	int N = (1.0 - 0.0) / h + 1;
	a = new double[N];
	b = new double[N];
	c = new double[N];
	y = new double[N];
	f = new double[N];
	u = new double[N];
	
	init(N, a, b, c, f, u, alp, h, l);
	//Заполняем трехдиагональную матрицу
	double **A = new double *[N];
	for (int i = 0; i < N; i++) {
		A[i] = new double[N];
		y[i] = 0;
	}

	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			A[i][j] = 0.0;
		}
	}
	for (int j = 0; j < N; j++) {
		if (j != 0) {
			A[j - 1][j] = b[j - 1];
		}
		A[j][j] = c[j];
		if (j != (N - 1)) {
			A[j + 1][j] = a[j + 1];
		}
	}
	gausZey(y, N, A, f);

	double norm = 0;
	for (int i = 0; i < N; i++)
	{
		norm += (y[i] - u[i])*(y[i] - u[i]);
	}
	cout << sqrt(norm) << endl;
	system("pause");
	return 0;
}

