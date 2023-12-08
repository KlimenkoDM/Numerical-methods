#pragma once
#include<iostream>
#include<iomanip>
using namespace std;
const int SIZE = 3;

void print(double(&matrix)[SIZE][SIZE + 1]);
void min_elem(double(&matrix)[SIZE][SIZE + 1], int a, int k);
void swap_rows(double(&matrix)[SIZE][SIZE + 1], int row1, int row2);
void gauss(double(&matrix)[SIZE][SIZE + 1], double(&result)[SIZE]);
void vektor_nevyazki(double(&matrix)[SIZE][SIZE + 1], double(&result)[SIZE]);
void rezult(double(&matrix)[SIZE][SIZE + 1], double(&result)[SIZE]);
void CalculationOfAx(double(&matrix)[SIZE][SIZE + 1], double* Ax, double** matrixs, double* x);
void CalculationOfdelta(int SIZE, double* x, double* xs, double& delta);

void print(double(&matrix)[SIZE][SIZE + 1])
{
	for (int i = 0; i < SIZE; i++)
	{
		for (int j = 0; j < SIZE + 1; j++)
		{
			cout << matrix[i][j] << " ";
		}
		cout << endl;
	}
}

void min_elem(double(&matrix)[SIZE][SIZE + 1], int a, int k)
{
	for (int i = 0; i < 3; i++)
	{
		cout << "Result = " << matrix[0][i] << endl;
		if (abs(a) < matrix[0][i])
		{
			a = matrix[0][i];
			k = i;
		}
	}
}

void swap_rows(double(&matrix)[SIZE][SIZE + 1], int row1, int row2)   //функци€ мен€юща€ строки
{
	for (int i = 0; i <= SIZE; ++i)
	{
		double temp = matrix[row1][i];
		matrix[row1][i] = matrix[row2][i];
		matrix[row2][i] = temp;
	}
}

void gauss(double(&matrix)[SIZE][SIZE + 1], double(&result)[SIZE])
{
	for (int i = 0; i < SIZE; ++i) {
		int max_row = i;
		for (int j = i + 1; j < SIZE; ++j)
		{
			if (matrix[j][i] > matrix[max_row][i])
			{
				max_row = j;
			}
		}
		swap_rows(matrix, i, max_row);

		for (int j = i + 1; j < SIZE; ++j) {
			double factor = matrix[j][i] / matrix[i][i];
			for (int k = i; k <= SIZE; ++k)
			{
				matrix[j][k] -= factor * matrix[i][k];
			}
		}
	}

	for (int i = SIZE - 1; i >= 0; --i)
	{
		result[i] = matrix[i][SIZE];
		for (int j = i + 1; j < SIZE; ++j)
		{
			result[i] -= matrix[i][j] * result[j];
		}
		result[i] /= matrix[i][i];
	}
}

void vektor_nevyazki(double(&matrix)[SIZE][SIZE + 1], double(&result)[SIZE])   //подсчЄт вектора нев€зки
{
	for (int i = 0; i < SIZE; ++i)
	{
		double resid = 0;
		for (int j = 0; j < SIZE; ++j)
		{
			resid += matrix[i][j] * result[j];
		}
		resid -= matrix[i][SIZE];
		cout << "Nevyazka " << i + 1 << ": " << setprecision(25) << resid << endl;
	}
}

void rezult(double(&matrix)[SIZE][SIZE + 1], double(&result)[SIZE])   //вывод результата
{
	for (int i = 0; i < SIZE; ++i)
	{
		cout << "x" << i + 1 << " = " << result[i] << endl;
	}
	int c = 0;
	for (int i = 0; i < SIZE; i++)
	{
		c += matrix[0][i];
	}
}

void CalculationOfAx(double(&matrix)[SIZE][SIZE + 1], double* Ax, double** matrixs, double* x) {
	for (int i = 0; i < SIZE; i++) {
		for (int j = 0; j < SIZE; j++) {
			Ax[i] = Ax[i] + (matrixs[i][j] * x[j]);
		}
	}
	for (int i = 0; i < SIZE; i++) {

		for (int j = 0; j < SIZE; j++) {

			matrix[i][j] = matrixs[i][j];
		}
	}
}

void CalculationOfdelta(int SIZE, double* x, double* xs, double& delta) {
	double max = 0;
	double max1 = 0;
	for (int i = 0; i < SIZE; i++) {
		if (abs(x[i]) > max) {
			max = abs(x[i]);
		}
		if (abs(xs[i] - x[i]) > max1) {
			max1 = abs(xs[i] - x[i]);
		}
		if (i == SIZE - 1) {
			delta = abs(max1) / abs(max);
		}
	}
}