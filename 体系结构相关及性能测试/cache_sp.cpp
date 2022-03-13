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
    int step = 64;
    int counter;
    float seconds;
    clock_t start, finish;
    for(int i=step; i<=N; i+=step)
    {
        start = clock();
        counter = 0;
        while(clock() - start <10)
        {
            counter++;
            normal(i,sum,a,b);
        }
        finish = clock();
        seconds = (finish - start)/float(CLOCKS_PER_SEC);
        cout << i << " " << counter << " " << seconds
        << " " << seconds / counter << endl ;
    }

    
   
   
    cout<<"hello"<<endl;
    system("pause");
    return 0;
}