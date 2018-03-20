#include "stdio.h"
#include "string.h"
#include "math.h"
#include "time.h" 
#include "windows.h"

#define MAX_SIZE 65536
#define PI	(acos((double)-1))

LARGE_INTEGER nFreq;//时钟频率 
LARGE_INTEGER t1;//开始 
LARGE_INTEGER t2;//结束 
double dt[4]; //时间差

struct Complex
{
	double real;	//实部 
	double imag;	//虚部 
};	//complex 

Complex A[MAX_SIZE],B[MAX_SIZE],C[MAX_SIZE],w[MAX_SIZE];

//以下为重载定义复数的四则运算和虚部取反 
Complex operator+(Complex a,Complex b)
{
	Complex r;
	r.real = a.real + b.real;
	r.imag = a.imag + b.imag;
	return r;
}	//plus

Complex operator-(Complex a,Complex b)
{
	Complex r;
	r.real = a.real - b.real;
	r.imag = a.imag - b.imag;
	return r;
}	//minus

Complex operator*(Complex a,Complex b)
{
	Complex r;
	r.real = a.real * b.real - a.imag * b.imag;
	r.imag = a.real * b.imag + a.imag * b.real;
	return r;
}	//mul

Complex operator/(Complex a,double b)
{
	Complex r;
	r.real = a.real / b;
	r.imag = a.imag / b;
	return r;
}	//div

Complex operator~(Complex a)
{
	Complex r;
	r.real = a.real;
	r.imag = - a.imag;
	return r;
}	//neg

void Reverse(int* id,int size,int m)
{
    for(int i = 0;i < size;i++)
	{
        for(int j = 0;j < (m + 1) / 2;j++)
		{
            int v1 = (1 << (j) & i) << (m - 2 * j - 1);
            int v2 = (1 << (m - j - 1) & i) >> (m - 2 * j - 1);
            id[i] |= (v1 | v2);
        }
    }
}; 	// 重新排列输入数组的元素下标 

void Compute_W(Complex w[],int size)
{
	for(int i = 0;i < size/2;i++)
	{
		w[i].real = cos(2 * PI * i/size);
		w[i].imag = sin(2 * PI * i/size);
		w[i + size/2].real = -w[i].real;
		w[i + size/2].imag = -w[i].imag;
	}
};	//预先计算FFT中需要的w值 

void FFT(Complex in[],int size)
{
	int* id = new int[size];
	memset(id,0,sizeof(int)*size);
	int m = log2((double) size);
	Reverse(id,size,m);	
	Complex *resort = new Complex[size];
	memset(resort,0,sizeof(Complex)*size);
	int i,j,k,s;
	for(i = 0;i < size;i++)
		resort[i] = in[id[i]];
	for(i = 1;i <= m;i++){
		s = (int) pow((double)2,(double)i);
		for(j = 0;j < size/s;j++)
		{
			for(k = j * s;k < j * s + s/2;k++)
			{		
				Complex k1 = resort[k] + w[size/s * (k - j * s)] * resort[k + s/2];
				resort[k + s/2] = resort[k] - w[size/s * (k - j * s)] * resort[k + s/2];
				resort[k] = k1;
			}
		}
	}
	for(i = 0;i < size;i++)
		in[i] = resort[i];
	delete[] id;
	delete[] resort;
};	//FFT 

void IFFT(Complex in[],int size)
{
	int* id = new int[size];
	memset(id, 0, sizeof(int)*size);
	int m = log2((double) size);
	Reverse(id,size,m);
	Complex *resort = new Complex[size];
	memset( resort, 0, sizeof(Complex)*size);
	int i,j,k,s;
	for(i = 0;i < size;i++)
		resort[i] = in[id[i]];
	for(i = 1;i <= m;i++){
		s = (int) pow((double)2,(double)i);
		for(j = 0;j < size/s;j++)
		{
			for(k = j*s;k < j * s + s/2;k++)
			{
				Complex k1 = (resort[k] + (~w[size/s * (k - j * s)]) * resort[k + s/2]);
				resort[k + s/2] = (resort[k] - (~w[size/s * (k - j * s)]) * resort[k + s/2]);
				resort[k] = k1;
			}
		}
	}
	for(i = 0;i < size;i++)
		in[i] = resort[i]/size;
	delete[] id;
	delete[] resort;
};	//IFFT 

void Random_Num(void)
{
	int i,m = 130;
	double num[m];
	srand((unsigned)time(NULL)); 
	for(i = 0;i < m;i++)
	{
		num[i] =rand() /(double)(RAND_MAX/20);
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
				fprintf(fp,"%.2lf\n",num[i]);
			}
			printf("successful!\n");
			fclose(fp);
		}
}

int main(void)
{	
	int i,j,count,m = 130,n = 4;
	double num[m] = {0};
	int arr[n] = {4,16,32,60};
//	Random_Num();
	QueryPerformanceFrequency(&nFreq);
	FILE *fp;
	for(count = 0;count < n;count++)
	{
		m = arr[count];
 		int m1 = pow(2,ceil(log2(2*m)));  //处理非2的整数次幂 
		if((fp = fopen("../input/input.txt","r"))== NULL)
			printf("error\n");
			else	
			{
				for(i = 0;i < m1;i++)
				{
					fscanf(fp,"%lf",&num[i]);
				}
				printf("read successful!\n");
				fclose(fp);
			} 
		memset(A,0,sizeof(A));
		memset(B,0,sizeof(B));
		memset(w,0,sizeof(w));
		memset(C,0,sizeof(C));
		for(i = 0;i < m;i++)
			A[i].real = num[i];
		for(i = 0;i < m;i++)
			B[i].real = num[i + m];
	
		QueryPerformanceCounter(&t1);
		Compute_W(w,m1);
		FFT(A,m1);
		FFT(B,m1);
		for(i = 0;i < m1;i++)
			C[i] = A[i] * B[i];
		IFFT(C,m1);
		QueryPerformanceCounter(&t2);
		
		dt[count] = (t2.QuadPart - t1.QuadPart )/ (double)nFreq.QuadPart *1000000;
		
		if((fp = fopen("../output/result.txt","a"))== NULL)
			printf("error\n");
			else
			{
				for(i = 0;i < 2*m-1;i++)
					fprintf(fp,"%.2lf ",C[i].real);
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
