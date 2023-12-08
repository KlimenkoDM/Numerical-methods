#include<iostream>
#include<cmath>
using namespace std;

double eps1 = 1e-04;
double eps2 = 1e-05;

double func1(double x)
{
	return sqrt(1+2*x*x*x);
}

double func2(double x, double y)
{
	return x * x / (1 + y * y);
}

double Trapec(double a, double b)
{
	double h = (b - a) / 2;
	double sum1 = 0;
	double sum2 = 0;
	int k = 0;
	bool fl = false;


	while (true)
	{

		sum2 = sum1;

		double n = (b - a) / h;

		sum1 = func1(a);
		sum1 += func1(b);

		for (int i = 1; i < n - 1; i++)
		{
			sum1 += 2 * func1(a + h * i);
		}

		sum1 /= 2;
		sum1 *= h;


		h /= 2;
		k++;

		cout << sum2 << " " << sum1 << endl;

		if (abs(sum1 - sum2) <= 3 * eps1 && fl == false)
		{
			cout << "Trapec: " << sum1 << " K: " << k << " Del: " << sum1 - sum2 << " e: " << eps1 << endl;
			fl = true;
		}

		if (abs(sum1 - sum2) <= 3 * eps1)
			break;

	}

	cout << "Trapec: " << sum1 << " K: " << k << " Del: " << sum1 - sum2 << " e: " << eps2;

	return sum1;
}

double Simpson(double a, double b)
{

	double h = (b - a) / 2;
	double sum1 = 0;
	double sum2 = 0;
	int k = 0;
	bool fl = false;


	while (true)
	{

		sum2 = sum1;

		double n = (b - a) / h;

		sum1 = func1(a);
		sum1 += func1(b);

		for (int i = 1; i < n - 1; i += 2)
		{
			sum1 += 2 * func1(a + h * i);
		}

		for (int i = 1; i < n; i += 2)
		{
			sum1 += 4 * func1(a + h * (i - 1));
		}

		sum1 /= 3;
		sum1 *= h;


		h /= 2;
		k++;

		if (abs(sum1 - sum2) <= 15 * eps1 && fl == false)
		{
			cout << "Simpson: " << sum1 << " K: " << k << " Del: " << sum1 - sum2 << " e: " << eps1 << endl;
			fl = true;
		}

		if (abs(sum1 - sum2) <= 15 * eps2)
			break;

	}

	cout << "Simpson: " << sum1 << " K: " << k << " Del: " << sum1 - sum2 << " e: " << eps2;

	return sum1;
}


double CubeSimpson(double a, double b, double c, double d)
{
	double hx = (b - a) / 2;
	double hy = (d - c) / 2;
	double sum1 = 0;
	double sum2 = 0;
	int k = 0;
	bool fl = false;


	while (true)
	{

		sum2 = sum1;
		sum1 = 0;

		double n = (b - a) / hx;
		double m = (d - c) / hy;

		for (int i = 0; i < n - 1; i += 2)
		{
			for (int j = 0; j < m - 1; j += 2)
			{
				sum1 += func2(a + hx * i, c + hy * j);
				sum1 += 4 * func2(a + hx * (i + 1), c + hy * j);
				sum1 += func2(a + hx * (i + 2), c + hy * j);
				sum1 += 4 * func2(a + hx * i, c + hy * (j + 1));
				sum1 += 16 * func2(a + hx * (i + 1), c + hy * (j + 1));
				sum1 += 4 * func2(a + hx * (i + 2), c + hy * (j + 1));
				sum1 += func2(a + hx * i, c + hy * (j + 2));
				sum1 += 4 * func2(a + hx * (i + 1), c + hy * (j + 2));
				sum1 += func2(a + hx * (i + 2), c + hy * (j + 2));
			}
		}

		sum1 *= hx * hy;
		sum1 /= 9;

		hx /= 2;
		hy /= 2;

		k++;

		if (abs(sum1 - sum2) <= 15 * eps1 && fl == false)
		{
			cout << "SimpsonCube: " << sum1 << " K: " << k << " Del: " << sum1 - sum2 << " e: " << eps1 << endl;
			fl = true;
		}

		if (abs(sum1 - sum2) <= 15 * eps2)
			break;

	}

	cout << "SimpsonCube: " << sum1 << " K: " << k << " Del: " << sum1 - sum2 << " e: " << eps2;

	return sum1;
}
int main()
{
	double a = 1.2;
	double b = 2.471;
	double c = 0;       //c И d для 29 варианта
	double d = 4.0;

	cout << 1 << endl;
	Trapec(a, b);
	cout << endl << endl;

	cout << 1 << endl;
	Simpson(a, b);
	cout << endl << endl;
	CubeSimpson(a, b, c, d);


	return 0;
}