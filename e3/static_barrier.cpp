#include <iostream>
#include <pthread.h>
//#include<arm_neon.h>
#include <sys/time.h>
#include <iomanip>
using namespace std;
const int X = 1024;
float A[X][X];
const int NUM_THREADS = 7;
//初始化数组
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

//串行算法
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
//并行算法
//线程数据结构定义
typedef struct
{
    int t_id; //线程 id
    int n; //数组规模
} threadParam_t;
// barrier 定义
pthread_barrier_t barrier_Divsion;
pthread_barrier_t barrier_Elimination;
//线程函数定义
void *threadFunc(void *param)
{
    threadParam_t *p = (threadParam_t *)param;
    int t_id = p->t_id;
    int n = p->n;
    for (int k = 0; k < n; ++k)
    {
        // t_id 为 0 的线程做除法操作，其它工作线程先等待
        // 这里只采用了一个工作线程负责除法操作，可以尝试采用多个工作线程完成除法操作
        if (t_id == 0)
        {
            for (int j = k + 1; j < n; j++)
                A[k][j] = A[k][j] / A[k][k];
            A[k][k] = 1.0;
        }
        //第一个同步点
        pthread_barrier_wait(&barrier_Divsion);
        //循环划分任务（可以尝试多种任务划分方式）
        for (int i = k + 1 + t_id; i < n; i += NUM_THREADS)
        {
            //消去
            for (int j = k + 1; j < n; ++j)
                A[i][j] = A[i][j] - A[i][k] * A[k][j];
            A[i][k] = 0.0;
        }
        // 第二个同步点
        pthread_barrier_wait(&barrier_Elimination);
    }
    pthread_exit(NULL);
}
bool static_barrier(int n)
{
    //初始化 barrier
    pthread_barrier_init(&barrier_Divsion, NULL, NUM_THREADS);
    pthread_barrier_init(&barrier_Elimination, NULL, NUM_THREADS);
    //创建线程
    pthread_t handles[NUM_THREADS];   // 创建对应的 Handle
    threadParam_t param[NUM_THREADS]; // 创建对应的线程数据结构
    for(int t_id = 0; t_id < NUM_THREADS; t_id++)
    {
        param[t_id].t_id = t_id;
        param[t_id].n = n;
        pthread_create(&handles[t_id], NULL, threadFunc, (void *)&param[t_id]);
    }
    for(int t_id = 0; t_id < NUM_THREADS; t_id++) 
    pthread_join(handles[t_id], NULL);
    //销毁所有的 barrier
    pthread_barrier_destroy(&barrier_Divsion);
    pthread_barrier_destroy(&barrier_Elimination);
    return true;
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
            static_barrier(i);
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