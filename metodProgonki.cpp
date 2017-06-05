/* 

1) Пример решения краевой задачи с помощью конечно-разностного метода:
(k(x)*u'(x))' - q(x)*u(x) = - f(x),
u(0) = m1,
u(l) = m2,

где k(x) = (1+alpha*x^2); q(x) = exp(alpha*x); m1=0; m2=0; f(x) = 10*x*(l-x)*exp(alpha*x) + 20*(1+alpha*x^2) - 2*alpha*x*(10*l-20*x); 0<x<1

2) Пример решения начально-краевой задачи с помощью конечно-разностного метода:
du(x,t)/dt = d/dx(k(x)*du(x,t)/dx) - q(x)*u + f(x,t)
u(0,t) = m1(t),
u(l,t) = m2(t),
u(x,0) = q(x)

Для решения СЛАУ используется метод прогонки
Описание метода в книге: Самарский - Теория разностных схем

*/

#include "stdafx.h"
#include <math.h>
#include <iostream>
using namespace std;

/**
* n - число уравнений (строк матрицы)
* b - диагональ, лежащая над главной (нумеруется: [0;n-2])
* c - главная диагональ матрицы A (нумеруется: [0;n-1])
* a - диагональ, лежащая под главной (нумеруется: [1;n-1])
* f - правая часть (столбец)
* x - решение, массив x будет содержать ответ
*/
void solveMatrix(int n, double *a, double *c, double *b, double *f, double *x)
{
	double m;
	for (int i = 1; i < n; i++)
	{
		m = a[i] / c[i - 1];
		c[i] = c[i] - m*b[i - 1];
		f[i] = f[i] - m*f[i - 1];
	}
	x[0] = 0.0;
	x[n - 1] = 0.0;

	for (int i = n - 2; i >= 1; i--) {
		x[i] = (f[i] - b[i] * x[i + 1]) / c[i];
	}

}

//Функция k(x)
double funcK(double x, double alp) {
	return (1.0 + alp*x*x);
}

//Известное решение задачи u(x) - используется для подсчета невязки
double funcU(double x) {
	double u = 10.0 * x*(1.0 - x);
	return u;
}

// Задача записывается в виде СЛАУ
void init(int N, double *a, double *b, double *c, double *f, double *u, double alp, double h, double l)
{
	double xi;
	for (int i = 0; i < N; i++) {
		xi = 0.0 + i*h;
		if (i != (N - 1)) {
			a[i] = funcK(xi - 0.5*h, alp) / (h*h);
		}
		else {
			a[i] = 0.0;
		}

		if (i != 0) {
			b[i] = funcK(xi + 0.5*h, alp) / (h*h);
		}
		else {
			b[i] = 0.0;
		}
		if (i == 0) {
			c[i] = 1;
		}
		else if (i == (N - 1)) {
			c[i] = 1;
		}
		else {
			c[i] = -(a[i] + b[i] + exp((double)(alp*xi)));
		}
	}
	for (int i = 0; i < N; i++) {
		xi = 0.0 + i*h;
		if (i == 0 || i == (N - 1)) {
			f[i] = 0;
		}
		else {
			f[i] = -(10.0 * xi * (l - xi) * exp((double)(alp * xi)) + 20.0 * (1.0 + alp * xi * xi) - 2.0 * alp * xi * (10.0 * l - 20.0 * xi));
		}
		u[i] = funcU(xi);
	}
	cout << endl;
}

//Функция k(x) для уравнения теплопроводности
double funcK_Tao(double x, double alp) {
	return (cos(x) + 2 * alp);
}

//Функция решения u(x,t) для уравнения теплопроводности
double funcU_Tao(double x, double t, double alp, double l) {
	double u = (l-x)*x*(t*t + alp);
	return u;
}

//Функция f(x,t) для уравнения теплопроводности
double funcF_Tao(double x, double tao, double alp, double l, double t) 
{
	return ((t*t+alp)*((l-2*x)*sin(x) + 2*cos(x) + 4*alp) + x*(l-x)*(2*t+x*(t*t+alp)*exp((double)(alp*x))));
}

// Задача записывается в виде СЛАУ в случае уравнения теплопроводности
void init_Tao(int N, double *a, double *b, double *c, double *f, double *u, double alp, double h, double l, double tao, double *y, double t)
{
	double xi, fi;
	for (int i = 0; i < N; i++) {
		xi = 0.0 + i*h;
		if (i != (N - 1)) {
			a[i] = funcK_Tao(xi - 0.5*h, alp) / (h*h);
		}
		else {
			a[i] = 0.0;
		}

		if (i != 0) {
			b[i] = funcK_Tao(xi + 0.5*h, alp) / (h*h);
		}
		else {
			b[i] = 0.0;
		}
		if (i == 0) {
			c[i] = 1;
		}
		else if (i == (N - 1)) {
			c[i] = 1;
		}
		else {
			c[i] = -(a[i] + b[i] + xi*exp((double)(alp*xi)) + 1/(tao));
		}
	}
	for (int i = 0; i < N; i++) {
		xi = 0.0 + i*h;
		if (i == 0 || i == (N - 1)) {
			f[i] = 0;
		}
		else {
			f[i] = -(funcF_Tao(xi, tao, alp, l, t) + y[i]/tao);
		}
		u[i] = funcU_Tao(xi, t, alp, l);
	}
	cout << endl;
}

int main()
{
	double *a, *b, *c;
	double l = 1.0;
	double h = 0.05;
	double alp = 9.0;
	double tao = 0.05;//1.0/2*h*h;
	double *y, *f, *u;
	double xi;
	int N = (1.0 - 0.0) / h + 1.0;
	int K = (1.0 - 0.0) / tao + 1.0;
	cout << N;
	a = new double[N];
	b = new double[N];
	c = new double[N];
	y = new double[N];
	f = new double[N];
	u = new double[N];
	if (tao == 0) {
		init(N, a, b, c, f, u, alp, h, l);
		solveMatrix(N, a, c, b, f, y);
	}
	else { //решение для задачи с уравнением теплопроводности
		for (int j = 0; j < K; j++) {
			double t = j*tao;
			if (j == 0) {
				for (int j = 0; j < N; j++) {
					xi = 0.0 + j*h;
					y[j] = alp*xi*(l - xi);
				}
			}
			init_Tao(N, a, b, c, f, u, alp, h, l, tao, y, t);
			solveMatrix(N, a, c, b, f, y);
			double norm = 0;
			for (int i = 0; i < N; i++)
			{
				printf("%.3lf\n",u[i]);
				norm += (y[i] - u[i])*(y[i] - u[i]);
			}
			cout << sqrt(norm) << endl;
			cout <<j<< "---------------" << endl;
		}
	}
	double norm = 0;
	for (int i = 0; i < N; i++)
	{
		norm += (y[i] - u[i])*(y[i] - u[i]);
	}
	
	system("pause");
    return 0;
}

