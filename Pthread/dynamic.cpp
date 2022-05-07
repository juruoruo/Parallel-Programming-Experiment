#include <iostream>
#include <pthread.h>
//#include<arm_neon.h>
#include <sys/time.h>
#include <iomanip>
using namespace std;
const int X = 1024;
float A[X][X];

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

// serial
void normal(int n, float A[][X])
{
    for (int k = 0; k < n; k++)
    {
        for (int j = k + 1; j < n; j++)
        {
            A[k][j] /= A[k][k];
        }
        A[k][k] = 1.0;
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
// parallel
typedef struct
{

    int k; //消去的轮次

    int t_id; // 线程 id
    int n;    // the size of array
    int size;

} threadParam_t;
void *threadFunc(void *param)
{

    threadParam_t *p = (threadParam_t *)param;

    int k = p->k; //消去的轮次

    int t_id = p->t_id; //线程编号
    int N = p->n;

    for (int i = k + 1 + t_id; i < N; i += 7)
    {
        //消去
        for (int j = k + 1; j < N; ++j)
            A[i][j] = A[i][j] - A[i][k] * A[k][j];

        A[i][k] = 0.0;
    }
    pthread_exit(NULL);
}
void dynamic_pthread(int N)
{

    for (int k = 0; k < N; ++k)
    {

        //主线程做除法操作

        for (int j = k + 1; j < N; j++)
        {

            A[k][j] = A[k][j] / A[k][k];
        }

        A[k][k] = 1.0;

        //创建工作线程，进行消去操作

        int worker_count = 7; //工作线程数量

        pthread_t *handles = new pthread_t[worker_count]; // 创建对应的 Handle

        threadParam_t *param = new threadParam_t[worker_count]; // 创建对应的线程数据结构

        //分配任务

        for (int t_id = 0; t_id < worker_count; t_id++)
        {

            param[t_id].k = k;
            param[t_id].t_id = t_id;
            param[t_id].n = N; // the size of array;
        }

        //创建线程

        for (int t_id = 0; t_id < worker_count; t_id++)
            pthread_create(&handles[t_id], NULL, threadFunc, (void *)&param[t_id]);

        //主线程挂起等待所有的工作线程完成此轮消去工作
        for (int t_id = 0; t_id < worker_count; t_id++)
            pthread_join(handles[t_id], NULL);
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
            dynamic_pthread(i);
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