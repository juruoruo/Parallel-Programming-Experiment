#include<iostream>
#include<ctime>
using namespace std;

const int N = 2048;
float sum[N];
float b[N][N];
float a[N];

void random_init(int n, float a[], float b[][N])
{
    for(int i=0; i<n; i++)
    {
        a[i] = rand()/double(RAND_MAX);
        for(int j=0; j<n; j++)
        {
            b[i][j] = j + i;
        }
    }
    
    
}

void normal(int n, float sum[], float a[], float b[][N])
{
    for(int i=0; i<n; i++)
    {
        sum[i] = 0.0;
        for (int j=0; j<n; j++)
        {
            
            sum[i] += b[j][i] * a[j];
        }

    }
}
void sp(int n, float sum[], float a[], float b[][N])
{
    for (int i = 0; i < n; i++)
        sum[i] = 0.0;
    for (int j = 0; j < n; j++)
    {
        for (int i = 0; i < n; i++)
            sum[i] += b[j][i] * a[j];
    }
}
int main()
{
    normal(N,sum,a,b);
    sp(N,sum,a,b);
    cout<<"hello"<<endl;
    system("pause");
    return 0;
}