// APPROXIMATION OF FUNCTIONS BY THE METHOD OF LEAST SQUARES.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include <iostream>
#include <cmath>
#include <iomanip>

using namespace std;

void fillTandR(double* tarr, double* rarr)
{
    double dp1[21] = { 0.0, 5.0, 10.0, 15.0, 20.0, 25.0,30.0, 35.0, 40.0, 45.0, 50.0, 55.0, 60.0, 65.0, 70.0, 75.0, 80.0, 85.0, 90.0, 95.0, 100.0};
    double dp2[21] = { 1.00762, 1.00392, 1.00153, 1.00000, 0.99907, 0.99852, 0.99826, 0.99818, 0.99828, 0.99849, 0.99878, 0.99919, 0.99967, 1.00024, 1.0091, 1.00167, 1.00253, 1.00351, 1.00461, 1.00586, 1.00721 };

    for (int i = 0; i < 21; i++)
    {
        tarr[i] = dp1[i];
        rarr[i] = dp2[i];
    }
}

void CalcXSumm(double* POWERX, double* tarr, int m, int N)
{
    for (int i = 0; i < 2 * m; i++)
        POWERX[i] = 0;

    for (int j = 0; j < 2 * m; j++)
    {
        for (int i = 0; i < N; i++)
        {
            POWERX[j] += pow(tarr[i], j + 1);
        }
    }

}

double** CalcSumX(double* POWERX, int m, int N)
{
    double** SUMX = new double* [m + 1];

    for (int i = 0; i < m + 1; i++)
        SUMX[i] = new double[m + 1];



    for (int i = 0; i < m + 1; i++)
    {
        for (int j = 0; j < m + 1; j++)
            SUMX[i][j] = POWERX[i + j - 1];
    }

    SUMX[0][0] = N;

    return SUMX;
}

double* FormPRAW(double* tarr, double* rarr, int m, int N)
{
    double* PRAW = new double[m + 1];

    for (int i = 0; i < m + 1; i++)
    {
        PRAW[i] = 0;
    }

    for (int i = 0; i < m + 1; i++)
    {
        for (int j = 0; j < N; j++)
        {
            PRAW[i] += rarr[j] * pow(tarr[j], i);
        }
    }

    return PRAW;
}

double** FormMatrix(double* PRAW, double** SUMX, int m)
{
    double** Matrix = new double* [m + 1];

    for (int i = 0; i < m + 1; i++)
        Matrix[i] = new double[m + 2];

    for (int i = 0; i < m + 1; i++)
    {
        for (int j = 0; j < m + 1; j++)
            Matrix[i][j] = SUMX[i][j];
    }

    for (int i = 0; i < m + 1; i++)
    {
        Matrix[i][m + 1] = PRAW[i];
    }

    return Matrix;
}

void SwapLines(int x, int x2, double** a)
{
    for (int j = 0; j < 4; j++)
        swap(a[x][j], a[x2][j]);
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

void CoutAns(double* ans, int m)
{
    for (int i = 0; i < m + 1; i++)
        cout << ans[i] << " ";

    cout << endl;
}

double Dispersion(double* ans, double* tarr, double* rarr, int m, int N)
{
    double sum = 0;
    double calc = 0;

    for (int i = 0; i < N; i++)
    {
        sum = 0;

        for (int j = 0; j < m + 1; j++)
        {
            sum += ans[j] * pow(tarr[i], j);
        }

        calc += pow(rarr[i] - sum, 2);
    }

    return calc / (N - m - 1);
}

int main()
{
    setlocale(0, "");

    int m = 1;
    int N = 21;
    double disp;
    double omega;

    double* tarr = new double[21];
    double* rarr = new double[21];
    double* POWERX = new double[2 * m];
    double* ans = new double[m + 1];
    double** SUMX;
    double* PRAW;
    double** Matrix;


    fillTandR(tarr, rarr);

    CalcXSumm(POWERX, tarr, m, N);

    SUMX = CalcSumX(POWERX, m, N);

    PRAW = FormPRAW(tarr, rarr, m, N);

    Matrix = FormMatrix(PRAW, SUMX, m);

    ans = Gauss(Matrix, m + 1, m + 2);

    disp = Dispersion(ans, tarr, rarr, m, N);


    omega = sqrt(disp);


    cout << "Отклонение: " << omega << endl;
    CoutAns(ans, m);

    cout << ans[0]<<" + "<<ans[1]<<"x + "<<ans[2]<<"x^2 "<<ans[3]<<"x^3";


    return 0;
}