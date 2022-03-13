#include<iostream>
#include<windows.h>
#include<ctime>
using namespace std;

const int N = 4194304;
float a[N];
float sum;

void random_init(int n, float a[])
{
    for(int i=0; i<n; i++)
    {
        a[i] = rand()/double(RAND_MAX);
        a[i] = i;
    } 
}

void normal(int n, float a[])
{
    sum= 0;
    for(int i=0; i<n; i++)
    {
        sum += a[i];
    }
}
int main()
{
    int step = 262144;
    int counter ;
    
    long long head, tail , freq;
    QueryPerformanceFrequency((LARGE_INTEGER *)&freq);
    long long limit =   freq*0.1;
    for(int i=step; i<=N; i+=step)
    {
        counter = 0;
        QueryPerformanceCounter((LARGE_INTEGER *)&head);
        QueryPerformanceCounter((LARGE_INTEGER *)&tail);
        while(tail - head < limit)
        {
            counter++;
            normal(i, a);
            QueryPerformanceCounter((LARGE_INTEGER *)&tail);
        }
        float seconds = (tail - head)*1000 / freq;
        cout <<i<<"  "<<counter<<"  "<<seconds<<"  "<<seconds / counter<< endl ;
    }
    cout<<"hello"<<endl;
    system("pause");
    return 0;
}