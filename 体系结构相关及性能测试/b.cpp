#include<iostream>
#include<ctime>
using namespace std;

const int N = 4194304;
float a[N];
float sum,sum1,sum2;

void random_init(int n, float a[])
{
    for(int i=0; i<n; i++)
    {
        a[i] = rand()/double(RAND_MAX);
        a[i] = i;
    } 
}

void sp(int n, float a[])
{
    sum = sum1 = sum2 = 0;
    for(int i=0; i<n; i+=2)
    {
        sum1 += a[i];
        sum2 += a[i+1];
    }
    sum = sum1 + sum2;
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
    normal(N,a);
    sp(N,a);
    cout << "hello" << endl;
    system("pause");
    return 0;
}