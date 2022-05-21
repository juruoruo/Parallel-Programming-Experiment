#include <iostream>
#include <fstream>
#include <string>
#include <sys/time.h>
// x86
using namespace std;
int Column = 130;
int Elimination = 22;
int Result = 8;
int N = 1;
int quantity = 129;
bool spGauss(int x)
{
	//selection 决定读取哪个文件
    int selection = x;
    string Folders[] = {"1_130_22_8", "2_254_106_53", "3_562_170_53", "4_1011_539_263", "5_2362_1226_453",
    "6_3799_2759_1953","7_8399_6375_4535", "11_85401_5724_756"};
    struct Size
    {
        int a, b, c; //分别为矩阵列数，消元子个数和被消元行个数
    } Quantities[] = {{130, 22, 8}, {254, 106, 53}, {562, 170, 53}, {1011, 539, 262}, {2362, 1226, 453},
    {3799, 2759, 1953},{8399, 6375, 4535},{85401,5724,756}};
    cout<<Folders[selection]<<endl;
	ifstream EliminationFile("/home/juruo/Desktop/e3/data/" + Folders[selection] + "/1.txt", std::ios::binary), ResultFile;
	ResultFile.open("/home/juruo/Desktop/e3/data/" + Folders[selection] + "/2.txt", std::ios::binary);
	//ofstream resFile("/home/juruo/Desktop/e4/res.txt", ios::trunc);
	Column = Quantities[selection].a, Elimination = Quantities[selection].b, Result = Quantities[selection].c;
	N = (Column + 31) >> 5, quantity = Column - 1;
	if(!EliminationFile.is_open())
    {
		cout<<"no1 失败读取"<<endl;
        return false;
    }else cout<<"yes1 成功读取"<<endl;
	if(!ResultFile.is_open())
    {
		cout<<"no2 失败读取"<<endl;
        return false;
    }else cout<<"yes2 成功读取"<<endl;
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
		int first = 0, index = 0;
		while (true)
		{
			while (first < N && out[i][first] == 0)
            {
				first++;
                index = 0;
            }
			if (first >= N)
				break;
			int tmp = out[i][first] << index;
			while (tmp >= 0)
            {
				index++;
                tmp <<= 1;
            }
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
	    //将结果写回到文件中
		// int count = quantity;
		// for (int j = 0; j < N; j++)
		// {
		// 	int dense = out[i][j];
		// 	for (int l = 0; l < 32; l++)
		// 	{
		// 		if (dense == 0)
		// 		{
		// 			break;
		// 		}
		// 		else if (dense < 0)
		// 		{
		// 			resFile << count - l << ' ';
		// 		}
		// 		dense <<= 1;
		// 	}
		// 	count -= 32;
		// }
		// resFile << '\n';
	}
	for(int i=0; i<Elimination; i++)
        delete []row[i];
    delete []row;
    for(int i=0; i<Result; i++)
        delete []out[i];
    delete []out;
    delete []gotEliminator;
    return true;
}

int main()
{
	struct timeval start;
    struct timeval end;
    // gettimeofday(&start, NULL);
    // spGauss_pthread(7);
    // gettimeofday(&end, NULL);
    // double time = (end.tv_sec - start.tv_sec) + double((end.tv_usec - start.tv_usec)) / 1000000; //单位s;
    // cout << time << endl;
    for(int i=0; i<8; i++)
    {
        gettimeofday(&start, NULL);
        spGauss(i);
        gettimeofday(&end, NULL);
        double time = (end.tv_sec - start.tv_sec) + double((end.tv_usec - start.tv_usec)) / 1000000; //单位s;
        cout << time << endl;
    }

    return 0;
}