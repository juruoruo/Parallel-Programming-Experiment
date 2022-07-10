#include <iostream>
#include <fstream>
#include <string>
#include <sys/time.h>
#include <iomanip>
#include <xmmintrin.h> //SSE
#include <emmintrin.h> //SSE2
#include <pmmintrin.h> //SSE3
#include <tmmintrin.h> //SSSE3
#include <smmintrin.h> //SSE4.1
#include <nmmintrin.h> //SSSE4.2
#include <immintrin.h> //AVX、AVX2
// x86
using namespace std;
int Column = 130;
int Elimination = 22;
int Result = 8;
int N = 1;
int quantity = 129;
int NUM_THREADS = 8;
bool normal(int x)
{
    // selection 决定读取哪个文件
    int selection = x;
    string Folders[] = {"1_130_22_8", "2_254_106_53", "3_562_170_53", "4_1011_539_263", "5_2362_1226_453",
                        "6_3799_2759_1953", "7_8399_6375_4535", "11_85401_5724_756"};
    struct Size
    {
        int a, b, c; //分别为矩阵列数，消元子个数和被消元行个数
    } Quantities[] = {{130, 22, 8}, {254, 106, 53}, {562, 170, 53}, {1011, 539, 262}, {2362, 1226, 453}, {3799, 2759, 1953}, {8399, 6375, 4535}, {85401, 5724, 756}};
    ifstream EliminationFile("/home/juruo/Desktop/e3/data/" + Folders[selection] + "/1.txt", std::ios::binary), ResultFile;
    ResultFile.open("/home/juruo/Desktop/e3/data/" + Folders[selection] + "/2.txt", std::ios::binary);
    // ofstream resFile("/home/juruo/Desktop/e4/res.txt", ios::trunc);
    Column = Quantities[selection].a, Elimination = Quantities[selection].b, Result = Quantities[selection].c;
    N = (Column + 31) >> 5, quantity = Column - 1;
    // if(!EliminationFile.is_open())
    // {
    // 	cout<<"no1 失败读取"<<endl;
    //     return false;
    // }else cout<<"yes1 成功读取"<<endl;
    // if(!ResultFile.is_open())
    // {
    // 	cout<<"no2 失败读取"<<endl;
    //     return false;
    // }else cout<<"yes2 成功读取"<<endl;
    int **row = new int *[Elimination], **out = new int *[Result], *gotEliminator = new int[Column]{0};

    for (int i = 0; i < Elimination; i++)
    {
        row[i] = new int[N]{0};
        int column;
        char ch = ' ';
        EliminationFile >> column;
        int r = quantity - column;
        row[i][r >> 5] = 1 << (31 - (r & 31));
        EliminationFile.get(ch);
        gotEliminator[column] = i + 1;
        while (EliminationFile.peek() != '\r')
        {
            EliminationFile >> column;
            int diff = quantity - column;
            row[i][diff >> 5] += 1 << (31 - (diff & 31));
            EliminationFile.get(ch);
        }
    }

    //首先读取所有的被消元行
    for (int i = 0; i < Result; i++)
    {
        out[i] = new int[N]{0};
        int column;
        char ch = ' ';
        while (ResultFile.peek() != '\r')
        {
            ResultFile >> column;
            int diff = quantity - column;
            out[i][diff >> 5] += 1 << (31 - (diff & 31));
            ResultFile.get(ch);
        }
        ResultFile.get(ch);
    }

    for (int i = 0; i < Result; i += 1)
    {
        int first = 0, index = 0;
        while (true)
        {
            while (first < N && out[i][first] == 0)
                first++, index = 0;
            if (first >= N)
                break;
            int tmp = out[i][first] << index;
            while (tmp >= 0)
                index++, tmp <<= 1;
            int &state = gotEliminator[quantity - (first << 5) - index];
            if (!state == 0)
            {
                int *ptr = state > 0 ? row[state - 1] : out[~state];
                //修改位sse四路地址对齐并行
                int j = 0;
                unsigned long addr = (unsigned long)(&out[i][j]); //获取第一个元素的地址;
                //计算偏移量，使之对齐
                int offset = 4 - (addr % 16) / 4; //计算偏移量;
                //不对齐部分串行处理
                for (int p = 0; p < offset && (i + 1 + offset) < N; p++)
                {
                    out[i][j] ^= ptr[j];
                    j++;
                }
                for (; j + 4 <= N; j += 4)
                {
                    __m128 v1 = _mm_loadu_ps(&out[i][j]);
                    __m128 v2 = _mm_loadu_ps(&ptr[j]);
                    __m128 v1 = _mm_sub_ps(v1, v2);
                    _mm_storeu_ps(&out[i][j], v1);
                }
                //最后的部分进行串行处理
                for (; j < N; j++)
                {
                    out[i][j] ^= ptr[j];
                }
            }
            else
            {
                state = ~i;
                break;
            }
        }
    }
    for (int i = 0; i < Elimination; i++)
        delete[] row[i];
    delete[] row;
    for (int i = 0; i < Result; i++)
        delete[] out[i];
    delete[] out;
    delete[] gotEliminator;
    return true;
}
bool sp(int x)
{
    // selection 决定读取哪个文件
    int selection = x;
    string Folders[] = {"1_130_22_8", "2_254_106_53", "3_562_170_53", "4_1011_539_263", "5_2362_1226_453",
                        "6_3799_2759_1953", "7_8399_6375_4535", "11_85401_5724_756", "8_23045_18748_14325", "9_37960_29304_14921",
                        "11_85401_5724_756"};
    struct Size
    {
        int a, b, c; //分别为矩阵列数，消元子个数和被消元行个数
    } Quantities[] = {{130, 22, 8}, {254, 106, 53}, {562, 170, 53}, {1011, 539, 262}, {2362, 1226, 453}, {3799, 2759, 1953}, {8399, 6375, 4535}, {85401, 5724, 756}, {23045, 18748, 14325}, {37960, 29304, 14921}, {85401, 5724, 756}};
    ifstream EliminationFile("/home/juruo/Desktop/e3/data/" + Folders[selection] + "/1.txt", std::ios::binary), ResultFile;
    ResultFile.open("/home/juruo/Desktop/e3/data/" + Folders[selection] + "/2.txt", std::ios::binary);
    // ofstream resFile("/home/juruo/Desktop/e4/res.txt", ios::trunc);
    Column = Quantities[selection].a, Elimination = Quantities[selection].b, Result = Quantities[selection].c;
    N = (Column + 31) >> 5, quantity = Column - 1;
    // if(!EliminationFile.is_open())
    // {
    // 	cout<<"no1 失败读取"<<endl;
    //     return false;
    // }else cout<<"yes1 成功读取"<<endl;
    // if(!ResultFile.is_open())
    // {
    // 	cout<<"no2 失败读取"<<endl;
    //     return false;
    // }else cout<<"yes2 成功读取"<<endl;
    int **row = new int *[Elimination], **out = new int *[Result], *gotEliminator = new int[Column]{0};

    for (int i = 0; i < Elimination; i++)
    {
        row[i] = new int[N]{0};
        int column;
        char ch = ' ';
        EliminationFile >> column;
        int r = quantity - column;
        row[i][r >> 5] = 1 << (31 - (r & 31));
        EliminationFile.get(ch);
        gotEliminator[column] = i + 1;
        while (EliminationFile.peek() != '\r')
        {
            EliminationFile >> column;
            int diff = quantity - column;
            row[i][diff >> 5] += 1 << (31 - (diff & 31));
            EliminationFile.get(ch);
        }
    }

    //首先读取所有的被消元行
    for (int i = 0; i < Result; i++)
    {
        out[i] = new int[N]{0};
        int column;
        char ch = ' ';
        while (ResultFile.peek() != '\r')
        {
            ResultFile >> column;
            int diff = quantity - column;
            out[i][diff >> 5] += 1 << (31 - (diff & 31));
            ResultFile.get(ch);
        }
        ResultFile.get(ch);
    }

#pragma omp parallel for num_threads(NUM_THREADS)
    for (int i = 0; i < Result; i += 1)
    {
        int first = 0, index = 0;
        while (true)
        {
            while (first < N && out[i][first] == 0)
                first++, index = 0;
            if (first >= N)
                break;
            int tmp = out[i][first] << index;
            while (tmp >= 0)
                index++, tmp <<= 1;
            int &state = gotEliminator[quantity - (first << 5) - index];
            if (!state == 0)
            {
                int *ptr = state > 0 ? row[state - 1] : out[~state];
                for (int j = 0; j < N; j++)
                {
                    out[i][j] ^= ptr[j];
                }
            }
            else
            {
                state = ~i;
                break;
            }
        }
    }
    for (int i = 0; i < Elimination; i++)
        delete[] row[i];
    delete[] row;
    for (int i = 0; i < Result; i++)
        delete[] out[i];
    delete[] out;
    delete[] gotEliminator;
    return true;
}

int main()
{
    double counter_1;
    double counter_2;
    struct timeval start_1;
    struct timeval end_1;
    struct timeval start_2;
    struct timeval end_2;
    for (int i = 0; i < 8; i += 1)
    {
        //串行
        gettimeofday(&start_1, NULL);
        normal(i);
        gettimeofday(&end_1, NULL);

        //并行
        gettimeofday(&start_2, NULL);
        sp(i);
        gettimeofday(&end_2, NULL);
        double time_1 = (end_1.tv_sec - start_1.tv_sec) + double((end_1.tv_usec - start_1.tv_usec)) / 1000000; //单位s;
        double time_2 = (end_2.tv_sec - start_2.tv_sec) + double((end_2.tv_usec - start_2.tv_usec)) / 1000000; //单位s;
        cout << fixed << setprecision(6);
        cout << time_1 << "    " << time_2 << endl;
    }
    cout << "hello" << endl;
    return 0;
}