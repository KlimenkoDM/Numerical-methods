#include <iostream>
#include <cmath>
#include <iomanip>
using namespace std;


double** create(int n, int m)
{
    double** arr = new double* [n];
    for (int i = 0; i < n; i++)
        arr[i] = new double[m];
    return arr;
}

double func1(double x1, double x2)
{
    return cos(0.4 * x2 + x1 * x1) + x2 * x2 + x1 * x1 - 1.6;
}

double func2(double x1, double x2)
{
    return 1.5 * x1 * x1 - (x2 * x2) / 0.36 - 1;
}

void Jacob(double** J, double x1, double x2)
{
    J[0][0] = -2 * x1 * sin(0.4 * x2 + x1 * x1) + 2 * x1;
    J[1][0] = 3 * x1;
    J[1][0] = -0.4 * sin(0.4 * x2 + x1 * x1) + 2 * x2;
    J[1][1] = -(2 / 0.36) * x2;
}

void Jacob2(double** J, double x1, double x2, double M) {
    J[0][0] = ((func1(x1 + M * x1, x2) - func1(x1, x2)) / (M * x1));
    J[0][1] = ((func1(x1, x2 + M * x2) - func1(x1, x2)) / (M * x2));
    J[1][0] = ((func2(x1 + M * x1, x2) - func2(x1, x2)) / (M * x1));
    J[1][1] = ((func2(x1, x2 + M * x2) - func2(x1, x2)) / (M * x2));
}

void Nevyaz(double* arr, double x1, double x2)
{
    arr[0] = -func1(x1, x2);
    arr[1] = -func2(x1, x2);
}

void SwapLines(int line1, int line2, double** arr)
{
    for (int j = 0; j < 4; j++)
        swap(arr[line1][j], arr[line2][j]);
}

void NevPlJ(double** a, double** J, double* Nev, int n)
{
    for (int i = 0; i < n; i++)
    {
        a[i][n] = Nev[i];

        for (int j = 0; j < n; j++)
            a[i][j] = J[i][j];

    }
}

double* Gauss(double** arr, int n, int m)
{

    double maxi = arr[n - 1][m - 1];


    for (int i = 0; i < n; i++)
    {

        for (int l = i; l < n; l++)
        {
            if (arr[i][i] < arr[l][i])
                SwapLines(i, l, arr);
        }
        if (arr[i][i] == 0)
        {
            throw "Деление на ноль";
            exit(0);
        }

        maxi = arr[i][i];

        for (int j = i; j < m; j++)
        {
            arr[i][j] /= maxi;
        }

        for (int k = i + 1; k < n; k++)
        {
            double dp = arr[k][i];

            for (int j = i; j < m; j++)
                arr[k][j] -= arr[i][j] * dp;
        }
    }

    double* ans = new double[n];
    ans[n - 1] = arr[n - 1][m - 1];

    for (int i = n - 2; i >= 0; i--)
    {
        double sum = 0;
        for (int k = i + 1; k < n; k++)
        {
            sum += arr[i][k] * ans[k];

        }
        ans[i] = arr[i][n] - sum;
    }

    return ans;
}

double* Newton(int n, double x1, double x2) {
    int NIT = 100;
    double e = 1e-10;
    double** J;
    double** matrix;
    double* nevyaz = new double[n];
    double* del = new double[n];
    double* ans = new double[n + 1];
    double M1 = 0.001;
    double M2 = 0.5;
    double M3 = 0.1;
    ans[0] = x1;
    ans[1] = x2;

    double* ansG = new double[n];

    double del1 = x1, del2 = x2;
    int k = 0;
    J = create(n, n);
    matrix = create(n, n + 1);

    cout << "x1" << setw(10) << "x2" << setw(20) << "delta1" << setw(20) << "delta2" << setw(20) << "k" << endl;

    while (del1 > e || del2 > e)
    {
        if (k > NIT)
        {
            return nullptr;
        }

        Nevyaz(nevyaz, x1, x2);

        Jacob2(J, x1, x2, M1);

        NevPlJ(matrix, J, nevyaz, n);

        try {
            ansG = Gauss(matrix, n, n + 1);
        }
        catch (const char* mess)
        {
            cout << mess << " ";
            exit(0);
        }

        double maxdel1 = 0.0;
        double maxdel2 = 0.0;

        for (int i = 0; i < n; i++)
            ans[i] += ansG[i];

        maxdel1 = max(maxdel1, fabs(func1(ans[0], ans[1])));
        maxdel1 = max(maxdel1, fabs(func2(ans[0], ans[1])));

        for (int i = 0; i < n; i++)
        {

            if (fabs(ansG[i]) < 1)
                maxdel2 = max(maxdel2, fabs(ans[i] - (ans[i] - ansG[i])));

            if (fabs(ansG[i]) >= 1)
                maxdel2 = max(maxdel2, fabs((ans[i] - (ans[i] - ansG[i])) / ans[i]));
        }

        del1 = maxdel1;
        del2 = maxdel2;

        x1 = ans[0];
        x2 = ans[1];

        cout << setprecision(3) << x1 << setw(10) << x2 << setw(20) << del1 << setw(20) << del2 << setw(20) << k << endl << endl;

        k++;
    }
    return ans;
}

int main()
{

    int n = 2;

    double* answer = new double[n];

    answer = Newton(n, 1, -1);
    answer = Newton(n, -1, 1);

    if (answer == nullptr)
        return 0;

    cout << answer[0] << " " << answer[1];

    return 0;
}