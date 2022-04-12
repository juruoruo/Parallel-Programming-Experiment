#include<iostream>
#include<xmmintrin.h> //SSE
#include<emmintrin.h> //SSE2
#include<pmmintrin.h> //SSE3
#include<tmmintrin.h> //SSSE3
#include<smmintrin.h> //SSE4.1
#include<nmmintrin.h> //SSSE4.2
#include<immintrin.h> //AVX、AVX2
#include<windows.h>
#include<ctime>
using namespace std;
const int N = 1024;
float A[N][N];


void A_reset(float A[][N])
{
    for(int i=0;i<N;i++)
    {
        for(int j=0;j<i;j++)
            A[i][j]=0;
        A[i][i]=1.0;
    for(int j=i+1;j<N;j++)
        A[i][j]=rand();
    }
    for(int k=0;k<N;k++)
        for(int i=k+1;i<N;i++)
            for(int j=0;j<N;j++)
                A[i][j]+=A[k][j];
}

void sp(int n, float A[][N])
{
    for(int k=0; k<n; k++)
    {
        __m128 vt = _mm_set1_ps(A[k][k]);
        int j=k+1;
        for(; j+4<=n;j+=4)
        {
            __m128 va = _mm_loadu_ps(&A[k][j]);
            va = _mm_div_ps(va,vt);
            _mm_storeu_ps(&A[k][j],va);
        }
        for(; j<n; j++)
            A[k][j]=A[k][j]/A[k][k];
        A[k][k] = 1.0;
        for(int i=k+1; i<n; i++)
        {
            __m128 vaik = _mm_set1_ps(A[i][k]);
            int j=k+1;
            for(; j+4<=n; j+=4)
            {
                __m128 vakj = _mm_loadu_ps(&A[k][j]);
                __m128 vaij = _mm_loadu_ps(&A[i][j]);
                __m128 vx = _mm_mul_ps(vakj,vaik);
                vaij = _mm_sub_ps(vaij,vx);
                _mm_storeu_ps(&A[i][j],vaij);
            }
            for(; j<n; j++)
                A[i][j]=A[i][j] - A[k][j] * A[i][k];
            A[i][k] =0;
        }
    }
}
//串行
void normal(int n,float A[][N]) 
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




int main()
{

    A_reset(A);
    int step = 64;
    float counter ;
    float counter_1;
    long long head, tail , freq;
    long long head_1, tail_1;
    QueryPerformanceFrequency((LARGE_INTEGER *)&freq);
    long long limit =   freq*0.1;
    for(int i=step; i<=N; i+=step)
    {
        //串行
        A_reset(A);
        counter_1 = 0;
        QueryPerformanceCounter((LARGE_INTEGER *)&head_1);
        QueryPerformanceCounter((LARGE_INTEGER *)&tail_1);
        while(tail_1-head_1<limit)
        {
            counter_1++;
            normal(i,A);
            QueryPerformanceCounter((LARGE_INTEGER *)&tail_1);
        }
        A_reset(A);

        //并行
        counter = 0;
        QueryPerformanceCounter((LARGE_INTEGER *)&head);
        QueryPerformanceCounter((LARGE_INTEGER *)&tail);
        while(tail - head < limit)
        {
            counter++;
            sp(i, A);
            QueryPerformanceCounter((LARGE_INTEGER *)&tail);
        }
        float seconds = (tail - head)*1000.0 / freq;
        float seconds_1 = (tail_1 - head_1)*1000.0 / freq;
        cout <<" "<<i<<"        ";
        cout <<counter<<"  "<<seconds<<"  "<<seconds / counter<<"        ";
        cout <<counter_1<<"  "<<seconds_1<<"  "<<seconds_1 / counter_1<< endl;
    }
    cout<<"hello world"<<endl;
    return 0;
}
