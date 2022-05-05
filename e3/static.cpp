#include <iostream>
#include <pthread.h>
#include <sys/time.h>
#include <iomanip>
#include <semaphore.h>
using namespace std;
const int X = 1024;
float A[X][X];
// float B[N][N];
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

const int NUM_THREADS = 7;
// parallel
typedef struct
{

    int t_id; // 线程 id
    int n;    // the size of array

} threadParam_t;
//信号量定义

sem_t sem_main;

sem_t sem_workerstart[NUM_THREADS]; // 每个线程有自己专属的信号量

sem_t sem_workerend[NUM_THREADS];
//线程函数定义

void *threadFunc(void *param)
{

    threadParam_t *p = (threadParam_t *)param;
    int t_id = p->t_id;
    int n = p->n;
    for (int k = 0; k < n; ++k)
    {
        sem_wait(&sem_workerstart[t_id]); // 阻塞，等待主线完成除法操作（操作自己专属的信号量）

        //循环划分任务
        for (int i = k + 1 + t_id; i < n; i += NUM_THREADS)
        {
            //消去
            for (int j = k + 1; j < n; ++j)
                A[i][j] = A[i][j] - A[i][k] * A[k][j];

            A[i][k] = 0.0;
        }
        sem_post(&sem_main);            // 唤醒主线程
        sem_wait(&sem_workerend[t_id]); //阻塞，等待主线程唤醒进入下一轮
    }

    pthread_exit(NULL);
}

void static_pthread(int n)
{

    //初始化信号量

    sem_init(&sem_main, 0, 0);
    for (int i = 0; i < NUM_THREADS; ++i)
    {
        sem_init(&sem_workerstart[i], 0, 0);
        sem_init(&sem_workerend[i], 0, 0);
    }

    //创建线程

    pthread_t handles[NUM_THREADS]; // 创建对应的 Handle

    threadParam_t param[NUM_THREADS]; // 创建对应的线程数据结构

    for (int t_id = 0; t_id < NUM_THREADS; t_id++)
    {
        param[t_id].t_id = t_id;
        param[t_id].n = n;
        pthread_create(&handles[t_id],NULL,threadFunc,(void *)&param[t_id]);
    }

    for (int k = 0; k < n; ++k)
    {
        //主线程做除法操作
        for (int j = k + 1; j < n; j++)
            A[k][j] = A[k][j] / A[k][k];
        A[k][k] = 1.0;
        //开始唤醒工作线程
        for (int t_id = 0; t_id < NUM_THREADS; ++t_id)
            sem_post(&sem_workerstart[t_id]);
        //主线程睡眠（等待所有的工作线程完成此轮消去任务）
        for(int t_id = 0;t_id < NUM_THREADS; ++t_id)
            sem_wait(&sem_main);

        // 主线程再次唤醒工作线程进入下一轮次的消去任务
        for (int t_id = 0; t_id < NUM_THREADS; ++t_id)
            sem_post(&sem_workerend[t_id]);
    }

    for (int t_id = 0; t_id < NUM_THREADS; t_id++)
        pthread_join(handles[t_id],NULL);

    //销毁所有信号量
    sem_destroy(&sem_main);
    sem_destroy(sem_workerend);
    sem_destroy(sem_workerstart);


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
            static_pthread(i);
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