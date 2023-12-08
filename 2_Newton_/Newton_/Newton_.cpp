#include <iostream>
#include <cmath>
using namespace std;

double f1(double x1, double x2)
{
   return cos(0.4 * x2 + x1 * x1) + x2 * x2 + x1 * x1 - 1.6;
}

double f2(double x1, double x2)
{
    return 1.5 * x1 * x1 - (x2 * x2) / 0.36 - 1;
}

double df1_dx1(double x1, double x2)
{
    return -2 * x1 * sin(0.4 * x2 + x1 * x1) + 2 * x1;
}

double df1_dx2(double x1, double x2)
{
    return -0.4 * sin(0.4 * x2 + x1 * x1) + 2 * x2;
}

double df2_dx1(double x1, double x2)
{
    return 3 * x1;
}

double df2_dx2(double x1, double x2)
{
    return -(2 / 0.36) * x2;
}
void newtonMethod(double& x1, double& x2, double epsilon)
{
    double x1_prev, x2_prev;
    int iterations = 0;

    cout << "Номер итерации\tАбсолютная погрешность x1\tАбсолютная погрешность x2" << endl;

    do
    {
        x1_prev = x1;
        x2_prev = x2;

        double J11 = df1_dx1(x1_prev, x2_prev);
        double J12 = df1_dx2(x1_prev, x2_prev);
        double J21 = df2_dx1(x1_prev, x2_prev);
        double J22 = df2_dx2(x1_prev, x2_prev);
        double detJ = J11 * J22 - J12 * J21;

        double invJ11 = J22 / detJ;
        double invJ12 = -J12 / detJ;
        double invJ21 = -J21 / detJ;
        double invJ22 = J11 / detJ;

        double f1_val = f1(x1_prev, x2_prev);
        double f2_val = f2(x1_prev, x2_prev);

        x1 = x1_prev - invJ11 * f1_val - invJ12 * f2_val;
        x2 = x2_prev - invJ21 * f1_val - invJ22 * f2_val;

        iterations++;

        double abs_error_x1 = abs(x1 - x1_prev);
        double abs_error_x2 = abs(x2 - x2_prev);
        cout << iterations << "\t\t" << abs_error_x1 << "\t\t\t" << abs_error_x2 << endl;

    } while (abs(x1 - x1_prev) > epsilon || abs(x2 - x2_prev) > epsilon);

    cout << "Решение: x1 = " << x1 << ", x2 = " << x2 << endl;
    cout << "Количество итераций: " << iterations << endl;

    // Вычисление вектора невязки
    double residual1 = f1(x1, x2);
    double residual2 = f2(x1, x2);

    cout << "Вектор невязки: " << residual1 << ", " << residual2 << endl;
}
int main()
{
    setlocale(LC_ALL, "ru");
    double x1 = 1;
    double x2 = -1;
    double epsilon = 1e-9;
    newtonMethod(x1, x2, epsilon);
    double e1 = 1e-9;
    double e2 = 1e-9;

    double x11 = -1;
    double x22 = 1;
    newtonMethod(x11, x22, epsilon);
    return 0;
}



