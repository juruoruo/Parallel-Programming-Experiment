#include<iostream>
#include<iomanip>
#include<pthread.h>
#include <sys/time.h>
#include<xmmintrin.h> //SSE
#include<emmintrin.h> //SSE2
#include<pmmintrin.h> //SSE3
#include<tmmintrin.h> //SSSE3
#include<smmintrin.h> //SSE4.1
#include<nmmintrin.h> //SSSE4.2
#include<immintrin.h> //AVX、AVX2
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


//串行算法:
void serialSolution(float** A, int N) {
	for (int k = 0; k < N; k++) {
		for (int j = k + 1; j < N; j++) {
			A[k][j] /= A[k][k];
		}
		A[k][k] = 1.0;
		for (int i = k + 1; i < N; i++) {
			for (int j = k + 1; j < N; j++) {
				A[i][j] -= A[i][k] * A[k][j];
			}
			A[i][k] = 0;
		}
	}
}

//并行算法:
typedef struct
{
	int N;//数组的规模
	int k;//消去的轮次
	int t_id;//线程id
	int worker_count;

} threadParam_t;

//线程函数:
void* threadFunc(void* param)
{
	//获取参数:
	threadParam_t* p = (threadParam_t*)param;
	int N = p->N;
	int k = p->k; //消去的轮次
	int t_id = p->t_id; //线程编号
	int worker_count = p->worker_count;

	for (int i = k + 1 + t_id; i < N; i += worker_count) {
		__m128 vaik = _mm_loadu_ps(&A[i][k]);
		for (int j = k + 1; j + 4 <= N; j += 4) {
			__m128 vakj = _mm_loadu_ps(&A[k][j]);
			__m128 vaij = _mm_loadu_ps(&A[i][j]);
			__m128 vx = _mm_mul_ps(vakj, vaik);
			vaij = _mm_sub_ps(vaij, vx);
			_mm_storeu_ps(&A[i][j], vaij);
			if (j + 8 > N) {//处理末尾
				while (j < N) {
					A[i][j] -= A[i][k] * A[k][j];
					j++;
				}
				break;
			}
		}
		A[i][k] = 0;
	}
	pthread_exit(NULL);
}

void sp(int N){
	for (int k = 0; k < N; ++k) {
		//主线程做除法操作,使用SIMD优化;
		__m128 vt = _mm_set1_ps(A[k][k]);
		for (int j = k + 1; j + 4 <= N; j += 4) {
			__m128 va = _mm_loadu_ps(&A[k][j]);
			va = _mm_div_ps(va, vt);
			_mm_storeu_ps(&A[k][j], va);
			if (j + 8 > N) {//处理末尾
				while (j < N) {
					A[k][j] /= A[k][k];
					j++;
				}
				break;
			}
		}
		A[k][k] = 1.0;

		//创建工作线程，进行消去操作
		int worker_count = 7; //工作线程数量
		pthread_t* handles = new pthread_t[worker_count]; // 创建对应的 Handle
		threadParam_t* param = new threadParam_t[worker_count]; // 创建对应的线程数据结构

		//分配任务
		for (int t_id = 0; t_id < worker_count; t_id++) {
			param[t_id].N = N;
			param[t_id].k = k;
			param[t_id].t_id = t_id;
			param[t_id].worker_count = worker_count;
		}

		//创建线程
		for (int t_id = 0; t_id < worker_count; t_id++) {
			pthread_create(&handles[t_id], NULL, threadFunc, (void*)&param[t_id]);
		}

		//主线程挂起等待所有的工作线程完成此轮消去工作
		for (int t_id = 0; t_id < worker_count; t_id++) {
			pthread_join(handles[t_id], NULL);
		}
	}
}
void *threadFunc_1(void *param)
{

    threadParam_t *p = (threadParam_t *)param;

    int k = p->k; //消去的轮次

    int t_id = p->t_id; //线程编号
    int N = p->N;

    for (int i = k + 1 + t_id; i < N; i += 7)
    {
        //消去
        for (int j = k + 1; j < N; ++j)
            A[i][j] = A[i][j] - A[i][k] * A[k][j];

        A[i][k] = 0.0;
    }
    pthread_exit(NULL);
}
void normal(int N)
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
            param[t_id].N = N; // the size of array;
        }

        //创建线程

        for (int t_id = 0; t_id < worker_count; t_id++)
            pthread_create(&handles[t_id], NULL, threadFunc_1, (void *)&param[t_id]);

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
            normal(i);
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
        	sp(i);
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