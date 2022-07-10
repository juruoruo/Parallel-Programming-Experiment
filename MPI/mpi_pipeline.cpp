#include <omp.h>
#include <sys/time.h>
#include <iomanip>
#include "mpi.h"
#include<iostream>
using namespace std;
int myid, numprocs;
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
void initMpi()
{
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
}
void sp(int N)
{
    
    MPI_Request request;
    MPI_Status status;
    int num = N / numprocs;
	int remainder = N % numprocs;
	if (myid < remainder) ++num;
	
        //进程成环
	int up = myid - 1; int down = myid + 1;
	if (myid == 0)up = numprocs - 1;
	if (myid == numprocs - 1)down = 0;
    int p = 0;
    for(int k=0; k<N; k++)
    {
        int source = k%numprocs;
        if(myid == source)
        {
            for(int j =k+1; j<N; j++)
                A[k][j] = A[k][j]/A[k][k];
            A[k][k] = 1;
            MPI_Isend(&A[k][0], N, MPI_FLOAT, down, k, MPI_COMM_WORLD, &request);
        }
        else
		{
			MPI_Irecv(&A[k][0], N, MPI_FLOAT, up, MPI_ANY_TAG, MPI_COMM_WORLD, &request);
            MPI_Wait(&request, &status);
            int picked = status.MPI_TAG;
			if (down != source)
				MPI_Isend(&A[k][0], N, MPI_FLOAT, down, picked, MPI_COMM_WORLD, &request);
        }
        p = myid - k%numprocs;
        if(p<=0) p+=numprocs;
        p += k;
        for (int i = p; i < N; i+=numprocs)
		{
            for(int j=k+1; j<N; j++)
                A[i][j] = A[i][j]-A[k][j]*A[i][k];
            A[i][k] = 0;
		}
    }
    
}
int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);

    
    initMpi();


    int step = 64;
    struct timeval start_1;
    struct timeval end_1;
    struct timeval start_2;
    struct timeval end_2;
    
    for(int i=step; i<=X; i+=step)
    {
        //串行
        A_reset(A);
        gettimeofday(&start_1, NULL);
        normal(i,A);
        MPI_Barrier(MPI_COMM_WORLD);
        gettimeofday(&end_1, NULL);
        //并行
        A_reset(A);
        gettimeofday(&start_2, NULL);
        sp(i);
        MPI_Barrier(MPI_COMM_WORLD);
        gettimeofday(&end_2, NULL);
        double time_1 = (end_1.tv_sec - start_1.tv_sec) + double((end_1.tv_usec - start_1.tv_usec)) / 1000000; //单位s;
        double time_2 = (end_2.tv_sec - start_2.tv_sec) + double((end_2.tv_usec - start_2.tv_usec)) / 1000000; //单位s;
        cout << fixed << setprecision(6);
        if(myid == 0)
        cout<<time_1<<"    "<<time_2<<endl;
    }
  
    MPI_Finalize();
    return 0;
}
