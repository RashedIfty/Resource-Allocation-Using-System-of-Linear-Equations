#include <bits/stdc++.h>
using namespace std;

const int N = 3;
double a[N][N];
double b[N];
double expectedResults[N];

// Function for Matrix Inversion Method
void matrixInversionMethod(double x[N])
{
    double A[N][N];
    double det_A = a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1]) -
                   a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0]) +
                   a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]);

    if (det_A == 0)
    {
        cout << "Matrix Inversion cannot be applied. Determinant is zero." << endl;
        return;
    }

    A[0][0] = (a[1][1] * a[2][2] - a[1][2] * a[2][1]) / det_A;
    A[0][1] = (a[0][2] * a[2][1] - a[0][1] * a[2][2]) / det_A;
    A[0][2] = (a[0][1] * a[1][2] - a[0][2] * a[1][1]) / det_A;
    A[1][0] = (a[1][2] * a[2][0] - a[1][0] * a[2][2]) / det_A;
    A[1][1] = (a[0][0] * a[2][2] - a[0][2] * a[2][0]) / det_A;
    A[1][2] = (a[0][2] * a[1][0] - a[0][0] * a[1][2]) / det_A;
    A[2][0] = (a[1][0] * a[2][1] - a[1][1] * a[2][0]) / det_A;
    A[2][1] = (a[0][1] * a[2][0] - a[0][0] * a[2][1]) / det_A;
    A[2][2] = (a[0][0] * a[1][1] - a[0][1] * a[1][0]) / det_A;

    for (int i = 0; i < N; i++)
    {
        x[i] = 0;
        for (int j = 0; j < N; j++)
        {
            x[i] += A[i][j] * b[j];
        }
    }

    cout << "Matrix Inversion Method: ";
    for (int i = 0; i < N; i++)
        cout << x[i] << " ";
    cout << endl;
}

// Function for Cramer's Rule
void cramersRule(double x[N])
{
    double tmp[N];
    double det = a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1]) -
                 a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0]) +
                 a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]);

    if (det == 0)
    {
        cout << "Cramer's Rule cannot be applied. Determinant is zero." << endl;
        return;
    }

    auto cd = [&](int row)
    {
        for (int col = 0; col < N; col++)
        {
            tmp[col] = a[col][row];
            a[col][row] = b[col];
        }

        double result = a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1]) -
                        a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0]) +
                        a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]);

        for (int col = 0; col < N; col++)
        {
            a[col][row] = tmp[col];
        }

        return result;
    };

    x[0] = cd(0) / det;
    x[1] = cd(1) / det;
    x[2] = cd(2) / det;

    cout << "Cramer's Rule: ";
    for (int i = 0; i < N; i++)
        cout << x[i] << " ";
    cout << endl;
}

// Function for Jacobi's Method
void jacobiMethod(double x[N])
{
    double preX, preY, preZ;
    x[0] = x[1] = x[2] = 0;
    do
    {
        preX = x[0];
        preY = x[1];
        preZ = x[2];
        x[0] = (b[0] - a[0][2] * preZ - a[0][1] * preY) / a[0][0];
        x[1] = (b[1] - a[1][2] * preZ - a[1][0] * preX) / a[1][1];
        x[2] = (b[2] - a[2][0] * preX - a[2][1] * preY) / a[2][2];
    }
    while (abs(preX - x[0]) > 0.001 || abs(preY - x[1]) > 0.001 || abs(preZ - x[2]) > 0.001);

    cout << "Jacobi's Method: ";
    for (int i = 0; i < N; i++)
        cout << x[i] << " ";
    cout << endl;
}

// Function for Gauss-Seidel Method
void gaussSeidelMethod(double x[N])
{
    x[0] = x[1] = x[2] = 0;
    int iterations = 10;
    while (iterations--)
    {
        x[0] = (b[0] - a[0][2] * x[2] - a[0][1] * x[1]) / a[0][0];
        x[1] = (b[1] - a[1][2] * x[2] - a[1][0] * x[0]) / a[1][1];
        x[2] = (b[2] - a[2][0] * x[0] - a[2][1] * x[1]) / a[2][2];
    }
    cout << "Gauss-Seidel Method: ";
    for (int i = 0; i < N; i++)
        cout << x[i] << " ";
    cout << endl;
}

// Function to calculate accuracy
void calculateAccuracy(double calculated[N], double expected[N], double accuracy[N])
{
    for (int i = 0; i < N; i++)
    {
        accuracy[i] = (abs(calculated[i] - expected[i]) / abs(expected[i])) * 100;
    }
}

// Function to get input
void getInput()
{
    cout << "Enter the coefficients of the matrix (" << N << "x" << N << "):" << endl;
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            cin >> a[i][j];
        }
    }
    cout << "Enter the constants (vector):" << endl;
    for (int i = 0; i < N; i++)
    {
        cin >> b[i];
    }
}

// Function to get expected results for accuracy check
void getExpectedResults()
{
    cout << "Enter the expected results (solution vector):" << endl;
    for (int i = 0; i < N; i++)
    {
        cin >> expectedResults[i];
    }
}

int main()
{
    getInput();

    int x;
    double results[N];
menu:
    cout << "1. Matrix Inversion" << endl;
    cout << "2. Cramer's Rule" << endl;
    cout << "3. Jacobi Method" << endl;
    cout << "4. Gauss-Seidel Method" << endl;
    cout << "5. Accuracy Check" << endl;
    cout << "6. Exit" << endl;

    cin >> x;

    switch (x)
    {
    case 1:
        matrixInversionMethod(results);
        break;
    case 2:
        cramersRule(results);
        break;
    case 3:
        jacobiMethod(results);
        break;
    case 4:
        gaussSeidelMethod(results);
        break;
    case 5:
    {
        getExpectedResults();
        double miResults[N], crResults[N], jmResults[N], gsResults[N];
        double miAccuracy[N], crAccuracy[N], jmAccuracy[N], gsAccuracy[N];

        matrixInversionMethod(miResults);
        cramersRule(crResults);
        jacobiMethod(jmResults);
        gaussSeidelMethod(gsResults);

        calculateAccuracy(miResults, expectedResults, miAccuracy);
        calculateAccuracy(crResults, expectedResults, crAccuracy);
        calculateAccuracy(jmResults, expectedResults, jmAccuracy);
        calculateAccuracy(gsResults, expectedResults, gsAccuracy);

        cout << "Matrix Inversion Method: ";
        for (int i = 0; i < N; i++)
        {
            cout << 100 - miAccuracy[i] << "% ";
        }
        cout << "Overall: " << 100 - (miAccuracy[0] + miAccuracy[1] + miAccuracy[2]) / 3 << "%" << endl;

        cout << "Cramer's Rule: ";
        for (int i = 0; i < N; i++)
        {
            cout << 100 - crAccuracy[i] << "% ";
        }
        cout << "Overall: " << 100 - (crAccuracy[0] + crAccuracy[1] + crAccuracy[2]) / 3 << "%" << endl;

        cout << "Jacobi Method: ";
        for (int i = 0; i < N; i++)
        {
            cout << 100 - jmAccuracy[i] << "% ";
        }
        cout << "Overall: " << 100 - (jmAccuracy[0] + jmAccuracy[1] + jmAccuracy[2]) / 3 << "%" << endl;

        cout << "Gauss-Seidel Method: ";
        for (int i = 0; i < N; i++)
        {
            cout << 100 - gsAccuracy[i] << "% ";
        }
        cout << "Overall: " << 100 - (gsAccuracy[0] + gsAccuracy[1] + gsAccuracy[2]) / 3 << "%" << endl;

        break;
    }
    case 6:
        cout << "Exiting program." << endl;
        return 0;
    default:
        cout << "Invalid option. Please try again." << endl;
    }

    goto menu;

    return 0;
}
