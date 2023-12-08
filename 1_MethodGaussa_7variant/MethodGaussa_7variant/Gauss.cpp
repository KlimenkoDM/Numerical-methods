#include<iostream>
#include<iomanip>
#include"Gauss.h"
using namespace std;

int main()
{
	double delta = 0;
	double Delta = 0;
	double* Ax = new double[SIZE];
	double* xs = new double[SIZE];
	double* x = new double[SIZE];
	double* matrixs[SIZE];
	double matrix[SIZE][SIZE + 1] = {
		{2.6, -4.5, -2, 19.07},
		{3, 3, 4.3, 3.21},
		{-6, 3.5, 3, -18.25}
	};
	double result[SIZE];
	int l = 0, k = 0;
	int a = matrix[0][0];

	cout << "Your matrix: " << endl;
	print(matrix);

	cout << endl;
	min_elem(matrix, a, k);
	cout << endl;
	cout << "Leading element " << a << endl;
	cout << endl;

	swap_rows(matrix, 0, 0);
	print(matrix);

	cout << endl;

	gauss(matrix, result);

	cout << "\nResult:" << endl;
	rezult(matrix, result);

	cout << endl << "Vectora nevyazki: " << endl;
	vektor_nevyazki(matrix, result);

	CalculationOfAx(matrix, Ax, matrixs, x);
	CalculationOfdelta(SIZE, x, xs, delta);
	cout << "delta:" << endl;
	cout << delta << endl;
	system("pause");
	return 0;
}
