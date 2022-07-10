#include <omp.h>
#include <sys/time.h>
#include <iomanip>
#include "mpi.h"
#include<iostream>
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
void sp(int N)
{
    int myid, numprocs;
    //char processor_name[MPI_MAX_PROCESSOR_NAME];
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    //MPI_Get_processor_name(processor_name, &namelen);
    //MPI_Barrier(MPI_COMM_WORLD);
    

    int r1, r2;
    int num = N/numprocs;
    r1 = myid * num;
    if(myid==numprocs-1) r2 = (myid+1) * N/numprocs ;
    else r2 = (myid+1) * num;
    //cout<<r1<<" "<<r2<<endl;
    // double start , end;
    // start = MPI_Wtime();
    for(int k=0; k<N; k++)
    {
        if(k>=r1 && k<r2)
        {
            for(int j =k+1; j<N; j++)
                A[k][j] = A[k][j]/A[k][k];
            A[k][k] = 1;
            for(int j=0; j<numprocs; j++)
            {   
                if(j!=myid)
                    MPI_Send(&A[k][0], N, MPI_FLOAT, j, k, MPI_COMM_WORLD);
            }
        }
        else    MPI_Recv(&A[k][0], N, MPI_FLOAT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD,MPI_STATUSES_IGNORE);
        for(int i=r1; i<r2; i++)
        {
            if(i<=k) continue; //做完除法的行不需要再进行消去操作。
            for(int j=k+1; j<N; j++)
                A[i][j] = A[i][j]-A[k][j]*A[i][k];
            A[i][k] = 0;
        }
    }
    //end = MPI_Wtime();
    MPI_Barrier(MPI_COMM_WORLD); /* IMPORTANT */
    // if (myid == 0)/* use time on master node */
    // cout<<end-start<<endl;
}
int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);

    int myid;
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);


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
        gettimeofday(&end_1, NULL);
        MPI_Barrier(MPI_COMM_WORLD);
        //并行
        A_reset(A);
        gettimeofday(&start_2, NULL);
        sp(i);
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
