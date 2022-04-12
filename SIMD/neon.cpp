#include<iostream>
#include<arm_neon.h>
#include<sys/time.h>
#include<iomanip>
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

void sp(int n,float A[][N])
{
    for(int k=0; k<n; k++)
    {  
        float32x4_t vt = vdupq_n_f32(A[k][k]);
        int j=k+1;
        for(; j+4<=n;j+=4)
        {
            float32x4_t va = vld1q_f32(&A[k][j]);
            va = vdivq_f32(va,vt);
            vst1q_f32(&A[k][j],va);
        }
        for(; j<n; j++)
            A[k][j]=A[k][j]/A[k][k];
        A[k][k] = 1.0;
        for(int i=k+1; i<n; i++)
        {
            float32x4_t vaik = vdupq_n_f32(A[i][k]);
            int j=k+1;
            for(; j+4<=n; j+=4)
            {
                float32x4_t vakj = vld1q_f32(&A[k][j]);
                float32x4_t vaij = vld1q_f32(&A[i][j]);
                float32x4_t vx = vmulq_f32(vakj,vaik);
                vaij = vsubq_f32(vaij,vx);
                vst1q_f32(&A[i][j],vaij);
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
	double counter_1;
	double counter_2;
	struct timeval start_1;
	struct timeval end_1;
	struct timeval start_2;
	struct timeval end_2;
	for (int i = step; i <= N; i += step)
	{
		//串行
		counter_1 = 0;
		gettimeofday(&start_1, NULL);
		gettimeofday(&end_1, NULL);
		while ((end_1.tv_sec - start_1.tv_sec) < 1)
        {
			counter_1++;
			normal(i,A);
			gettimeofday(&end_1, NULL);
		}
		//并行
		counter_2 = 0;
		gettimeofday(&start_2, NULL);
		gettimeofday(&end_2, NULL);
		while ((end_2.tv_sec - start_2.tv_sec) < 1) 
        {
			counter_2++;
			sp(i,A);
			gettimeofday(&end_2, NULL);
		}

		double time_1 = (end_1.tv_sec - start_1.tv_sec)+double((end_1.tv_usec - start_1.tv_usec))/1000000;//单位s;
		double time_2 = (end_2.tv_sec - start_2.tv_sec)+double((end_2.tv_usec - start_2.tv_usec))/1000000;//单位s;
        cout << fixed << setprecision(6);
		cout << "数组规模" << i <<endl;
		cout << "   串行重复次数" << counter_1 << "  " << "串行总用时(s)" << time_1 << "  " << "串行平均用时(s)" << time_1 / counter_1 << endl;
		cout << "   并行重复次数" << counter_2 << "  " << "并行总用时(s)" << time_2 << "  " << "并行平均用时(s)" << time_2 / counter_2 << endl;
	}
    cout<<"hello"<<endl;
    return 0;
}