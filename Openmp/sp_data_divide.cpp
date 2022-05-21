#include <iostream>
#include <pthread.h>
#include<omp.h>
#include <sys/time.h>
#include <iomanip>
using namespace std;
const int X = 1024;
float A[X][X];
int NUM_THREADS = 7;
//  Initial array
void A_reset(float A[][X])
{
    for (int i = 0; i < X; i++)
    {
        for (int j = 0; j < i; j++)
            A[i][j] = 0;
        A[i][i] = 1.0;
        for (int j = i + 1; j < X; j++)
            A[i][j] = rand();
    }
    for (int k = 0; k < X; k++)
        for (int i = k + 1; i < X; i++)
            for (int j = 0; j < X; j++)
                A[i][j] += A[k][j];
}

//并行的默认划分
void normal(int n, float A[][X])
{
    #pragma omp parallel num_threads(NUM_THREADS), shared(A)
    for (int k = 0; k < n; k++)
    {
        for (int j = k + 1; j < n; j++)
        {
            A[k][j] /= A[k][k];
        }
        A[k][k] = 1.0;
        #pragma omp for 
        for (int i = k + 1; i < n; i++)
        {
            for (int j = k + 1; j < n; j++)
            {
                A[i][j] -= A[i][k] * A[k][j];
            }
            A[i][k] = 0;
        }
    }
}
//并行动态循环调度循环划分
void sp(int n, float A[][X])
{
    #pragma omp parallel num_threads(NUM_THREADS), shared(A)
    for (int k = 0; k < n; k++)
    {
        for (int j = k + 1; j < n; j++)
        {
            A[k][j] /= A[k][k];
        }
        A[k][k] = 1.0;
        #pragma omp for schedule(dynamic)
        for (int i = k + 1; i < n; i++)
        {
            for (int j = k + 1; j < n; j++)
            {
                A[i][j] -= A[i][k] * A[k][j];
            }
            A[i][k] = 0;
        }
    }
}
int main()
{
    int step = 64;
    double counter_1;
    double counter_2;
    struct timeval start_1;
    struct timeval end_1;
    struct timeval start_2;
    struct timeval end_2;
    for (int i = step; i <= X; i += step)
    {
        //串行
        counter_1 = 0;
        A_reset(A);
        gettimeofday(&start_1, NULL);
        gettimeofday(&end_1, NULL);
        while ((end_1.tv_sec - start_1.tv_sec) < 1)
        {
            counter_1++;
            normal(i, A);
            gettimeofday(&end_1, NULL);
        }
        //并行
        counter_2 = 0;
        A_reset(A);
        gettimeofday(&start_2, NULL);
        gettimeofday(&end_2, NULL);
        while ((end_2.tv_sec - start_2.tv_sec) < 1)
        {
            counter_2++;
            sp(i,A);
            gettimeofday(&end_2, NULL);
        }

        double time_1 = (end_1.tv_sec - start_1.tv_sec) + double((end_1.tv_usec - start_1.tv_usec)) / 1000000; //单位s;
        double time_2 = (end_2.tv_sec - start_2.tv_sec) + double((end_2.tv_usec - start_2.tv_usec)) / 1000000; //单位s;
        cout << fixed << setprecision(6);
        cout<<time_1/counter_1<<"    "<<time_2/counter_2<<endl;
    }
    cout << "hello" << endl;
    return 0;
}