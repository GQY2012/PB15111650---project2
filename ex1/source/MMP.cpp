#include "stdio.h"
#include "stdlib.h"
#include "time.h" 
#include "windows.h"

#define INF 65536

LARGE_INTEGER nFreq;//时钟频率 
LARGE_INTEGER t1;//开始 
LARGE_INTEGER t2;//结束 
double dt[4]; //时间差 

int cost[32][32];
int s[32][32];

void MATRIX_CHAIN_ORDER(int p[],int n)
{
	int i,j,k,l,q;
	memset(s,0,sizeof(s));
	for(l = 2;l <= n;l++)
	{
		for(i = 1;i <= n - l + 1;i++)
		{
			j = i + l - 1;
			cost[i][j] = INF;
			for(k = i;k <= j - 1;k++)
			{
				q = cost[i][k] + cost[k + 1][j] + p[i - 1] * p[k] * p[j];
				if(q < cost[i][j])
				{
					cost[i][j] = q;
				    s[i][j] = k;
				}
			}		
		}
	}
}

void PRINT_OPTIMAL_PARENS(int i,int j,FILE *fp)
{
		if(i == j){
			fprintf(fp,"A%d",i);
//			printf("A%d",i);
		}
			else
			{
				fprintf(fp,"(");
//				printf("(");
				PRINT_OPTIMAL_PARENS(i, s[i][j], fp);
				PRINT_OPTIMAL_PARENS(s[i][j] + 1, j, fp);
				fprintf(fp,")");
//				printf(")");			
			}
}

void Random_Num(void)
{
	int i,m = 40;
	int num[m];
	srand((unsigned)time(NULL));//以系统时间作为随机数种子 
	for(i = 0;i < m;i++)
	{
		num[i] =1 + rand()%20;
	}
		FILE *fp;
		if((fp = fopen("../input/input.txt","w"))== NULL)
		{
			printf("error");
		}
		else
		{
			for(i = 0;i < m;i++)
			{
				fprintf(fp,"%d\n",num[i]);
			}
			printf("successful!\n");
			fclose(fp);
		}
}


int main(void)
{
	int i,j,count,m = 31,n = 4;
	int num[m];
	memset(num,0,sizeof(num));
	int arr[n] = {6,11,21,31};
//	Random_Num();
	QueryPerformanceFrequency(&nFreq);
	FILE *fp;
	for(count = 0;count < n;count++)
	{
		m = arr[count];
		if((fp = fopen("../input/input.txt","r"))== NULL)
			printf("error\n");
			else	
			{
				for(i = 0;i < m;i++)
				{
					fscanf(fp,"%d",&num[i]);
				}
				printf("read successful!\n");
				fclose(fp);
			}
		QueryPerformanceCounter(&t1);
		MATRIX_CHAIN_ORDER(num,m - 1);
		QueryPerformanceCounter(&t2);
		dt[count] = (t2.QuadPart - t1.QuadPart )/ (double)nFreq.QuadPart *1000000;
		if((fp = fopen("../output/result.txt","a"))== NULL)
		{
			printf("error");
		}
			else
			{
			PRINT_OPTIMAL_PARENS(1,m - 1,fp);
			fprintf(fp,"\n");
			printf("resuilt written successful!\n");
				}		
		fclose(fp);
		}
		if((fp = fopen("../output/time.txt","a"))== NULL)
			printf("error\n");
			else
			{
				for(count = 0;count < n;count++)
					fprintf(fp,"%f\n",dt[count]);
				printf("time written successful!\n");	
				}
		fclose(fp);	
	return 0;
}
