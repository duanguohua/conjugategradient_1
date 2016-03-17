/*
*function:conjugategradient_1.c
*purpose:学习共轭梯度法,来实现《陈正隆 (2007). 分子模拟的理论与实践, 化学工业出版社.》相应习题的练习
*version:v1.0
*date:15:13 2016/2/6
*author:dgh
*note:虽然这个文件名为，file.c，但它是原来 “conjugategradient_1.c”将两个文件组合在一起来。主要是出于学习understand这个软件方便这个目的

*/

#include <stdio.h>
#include <math.h>
/*----------------prototype---------------------*/
/*--------主干函数------------------------*/
double fun3d(double p[]);						//目标函数表达式
void gradient(double p[], double g[]);		//目标函数的梯度
double conjugateGradient(double p[], double pmin[]);	// 共轭梯度法
double steepDescend(double p[], double h[], double pmin[]);	//最速下降法实现
void goldenResearch3d(double a[], double b[], double xmin[], double *zmin); 	//某一维方向上的黄金分割搜索算法

/*--------枝叶函数------------------------*/
double getGamma(double gp0[], double gp1[]);	//采用P-R算法, 计算gamma
double getMod(double gp[]);					//计算当前自变量向量的梯度的模


void test(void);			//测试函数

int main(void)
{
	test();
	// printf("%lf\n", sqrt(-1));		//这时会报错：-1.#IND00

	
	return 0;
}


/*----------------函数实现---------------------*/
//测试函数
void test(void)
{
	double p[2] = {9, 9};
	double pmin[2];
	double zmin;
	double gp[2];
	// printf("%lf \t%lf \t%lf\n", p[0], p[1], fun(p));
	gradient(p, gp);
	// printf("%lf \t%lf\n", gp[0], gp[1]);
	// printf("%lf\n", getMod(gp));
	zmin = conjugateGradient(p, pmin);
	printf("x=(%lf, %lf), z = %lf\n", pmin[0], pmin[1], zmin);

	
}


// 共轭梯度法
/*
*参数说明：
double p[], 初始搜索位置（自变量向量）
double pmin[], 最小值所对应的自变量向量
返回值：最小值
*/
double conjugateGradient(double p[], double pmin[])
{
	double gp[2];
	gradient(p, gp);				//得到自变量向量p处的负梯度
	
	double p0[2], p1[2]; 			//当前, 下一个自变量向量
	double gp0[2], gp1[2];			//当前, 下一个自变量向量的负梯度
	double h0[2], h1[2];			//当前, 下一个的搜索方向	
	//当eps的值非常小时，如1e-5时，会导致-1.#1ND00错误，这个错误主要是发生在steepDescend()函数中，非负参数 lamda会变为负值，这时使用sqrt()函数就会出现问题，
	double eps = 1e-3;				//收敛标准
	
	int i, j;
	double modgp = 0.0;				//负梯度的模
	double gamma = 0.0;				//共轭梯度法中的一个系数
	double minz = 0.0;				//某一搜索方向上的最小值
	
	//初始值
	for (i=0; i<2; i++)
	{
		p0[i] = p[i];
		gp0[i] = gp[i];
		h0[i] = gp[i];
	}
	modgp = getMod(gp0);				//得到自变量向量p0处的负梯度的模
	j = 0;
	//当梯度的模大于某一个值时，继续搜索
	while (modgp >eps)
	{
		//改变p1, gp1, h1的值
		minz = steepDescend(p0, h0, p1);		//得到h0方向上最小值对应的自变量向量
		gradient(p1, gp1);				//得到自变量向量p1处的梯度
		printf("j=%d\t p1=(%lf, %lf), minz = %lf\n", j, p1[0], p1[1], minz);
		//采用P-R算法, 计算gamma
		gamma = getGamma(gp0, gp1);
		//计算下一个搜索方向
		for (i=0; i<2; i++)
		{
			h1[i] = gp1[i] + gamma * h0[i];
		}
		
		// 更新h0, gp0, p0
		for (i=0; i<2; i++)
		{
			h0[i] = h1[i];
			gp0[i] = gp1[i];
			p0[i] = p1[i];			
		}
		//更新循环判断条件
		modgp = getMod(gp0);				//得到自变量向量p0处的负梯度的模
		j++;
	}
	//返回最小值 对应的自变量向量
	for (i=0; i<2; i++)
	{
		pmin[i] = p1[i];
	}
	
	return minz;
	
}
//最速下降法实现
/*
*参数说明：
double p[], 初始自变量向量 
double h[], 初始搜索方向
double pmin[], 在这个方向上，找到的最小值相应的自变量向量
返回值：这个方向上，找到的最小值
*/
double steepDescend(double p[], double h[], double pmin[])
{
	int i, j;
	double xmin[2], zmin;		//最小值 对应的自变量和最小值 
	double ztemp, z0, z1;		//上一步，当前步，下一步的应变量
	double ptemp[2], p0[2], p1[2];			//上一步, 当前步，下一步的自变量
	//这个搜索步长如何确定呢？这要结合搜索的方向是单位向量还是其它，
	//如果是单位向量，可以取为1，这样可能会增加搜索的次数
	//如果非单位向量，可以取一个方向向量模的倒数。或者模倒数的几分之几。
	double lamda;				//一维搜索步长
	for (i=0; i<2; i++)
	{
		lamda += h[i] * h[i];
	}
	// printf("h = (%lf, %lf), lamda = %lf\n", h[0], h[1], lamda);
	//避免作用下面的计算表达，尤其当lamda为很小的数值时，
	// lamda = 1.0/sqrt(lamda);	//如此，后面的搜索步长相当于沿负梯度单位向量进行搜索
	lamda = sqrt(lamda)/lamda;
	// printf("lamda = %lf\n", lamda);
	
	//初始化：自变量和应变量
	j=0;
	z1 = ztemp = z0 = fun3d(p);
	for (i=0; i<2; i++)
	{
		p0[i] = p[i];
		p1[i] = p[i];
		ptemp[i] = p[i];
	}
	//沿某一维方向进行搜索，至到f(x_k)<f(x_k+1),
	while (z1<=z0)
	{		
		for (i=0; i<2; i++)		//记录上一步，当前步的自变量
		{
			ptemp[i] = p0[i];
			p0[i] = p1[i];
		}
		ztemp = z0;				//记录上一步，当前步的应变量
		z0 = z1;				
		j++;					//增加搜索步长
		for (i=0; i<2; i++)		//下一个步的自变量
		{
			//note: p1是在最初始位置基础上查找 的，而非上一步。
			//所以计算p1[i]时，要用到p[i]。
			//这也是之前没能注意到的，而出问题
			p1[i] = p[i] + j * lamda * h[i];
		}
		z1 = fun3d(p1);			//下一个步的应变量
	}

	//调用黄金分割法
	if (j == 1)		//表示第一步结束时，就出现z1>z0的情况
	{
		goldenResearch3d(p0, p1, xmin, &zmin);
	}
	else
	{
		goldenResearch3d(ptemp, p1, xmin, &zmin);
	}
	
	//返回最小值对应的自变量
	for (i=0; i<2; i++)
	{
		pmin[i] = xmin[i];
	}
	//返回最小值 
	return zmin;
}

//某一维方向上的黄金分割搜索算法
/*
double a, 为自变量向量，区间的左边
double b, 为自变量向量，区间的右边
double xmin[], 最小值 对应的x自变量向量
double *zmin  最小值 对应的y值
*/
void goldenResearch3d(double a[], double b[], double xmin[], double *zmin)
{
	int i;
	double tol[2];			//区间长度
	for (i=0; i<2; i++)
	{
		tol[i] = b[i] - a[i];
	}
	
	double eps = 1e-3;			//区间最小长度
	double golden = 0.618;		//黄金分割比
	int k = 0;
	
	// 插入的两个自变量向量与对应的z值
	double x1[2], x2[2], z1, z2;
	for (i=0; i<2; i++)
	{
		x1[i] = b[i] - golden * tol[i];
		x2[i] = a[i] + golden * tol[i];
	}
	z1 = fun3d(x1);
	z2 = fun3d(x2);

	while (fabs(tol[0]) > eps && fabs(tol[1]) > eps)
	{
		if (z1 >= z2)
		{
			for (i=0; i<2; i++)
			{
				a[i] = x1[i];
				x1[i] = x2[i];
				z1 = z2;
				tol[i] = b[i] - a[i];
				x2[i] = a[i] + golden * tol[i];
			}
			z2 = fun3d(x2);
		}
		else
		{
			for (i=0; i<2; i++)
			{
				b[i] = x2[i];
				x2[i] = x1[i];
				z2 = z1;
				tol[i] = b[i] - a[i];
				x1[i] = b[i] - golden * tol[i];			
			}
			z1 = fun3d(x1);
		}
		k++;
	}
	// printf("tol = %lf, k=%d\n", tol, k);
	//输出结果
	for (i=0; i<2; i++)
	{
		xmin[i] = (a[i] + b[i])/2;
	}	
	(*zmin) = fun3d(xmin);
}



//采用P-R算法, 计算gamma
/*
*参数说明：
double gp0, 自变量向量p0位置处的负梯度 
double gp1,	自变量向量p1位置处的负梯度
返回一个系数
*/
double getGamma(double gp0[], double gp1[])
{
	//采用P-R算法, 计算gamma
	double up, down;
	up = down = 0.0;
	double gamma = 0.0 ;
	int i;
	
	for (i=0; i<2; i++)
	{
		up += (gp1[i] - gp0[i]) * gp1[i];
		down += gp0[i] * gp0[i];
	}
	gamma = up / down;	
	return gamma;
}

//计算当前自变量向量的负梯度的模
double getMod(double g[])
{
	int i;
	double sumsqr = 0.0;
	
	for (i=0; i<2; i++)
	{
		sumsqr += g[i] * g[i];
	}
	sumsqr = sqrt(sumsqr);
	return sumsqr;
}



//目标函数表达式
//注意修改目标函数时，也要同时对其梯度函数进行修改
/*
*参数说明：
double p[], 为一个自变量向量
返回值为函数的应变量。
*/
double fun3d(double p[])
{
	double x, y, z;
	x = p[0];
	y = p[1];
	z = x * x + 2 * y * y;
	return z;
}

//目标函数的负梯度
//注意修改目标函数时，也要同时对其梯度函数进行修改
/*
*参数说明；
double p[], 为一个自变量向量 
double g[], 利用指南针来返回值， 返回在自变量向量处，函数的负梯度向量。
*/
void gradient(double p[], double g[])
{
	double x, y;
	double dfx, dfy;
	x = p[0];
	y = p[1];
	dfx = 2 * x;
	dfy = 4 * y;
	g[0] = - dfx;
	g[1] = - dfy;
}



/*
*answer:
j=0      p1=(3.999869, -1.000261), minz = 18.000000
j=1      p1=(-0.000869, -0.000100), minz = 0.000001
j=2      p1=(-0.000551, 0.000536), minz = 0.000001
j=3      p1=(-0.000033, 0.000213), minz = 0.000000
x=(-0.000033, 0.000213), z = 0.000000
请按任意键继续. . .

*/
