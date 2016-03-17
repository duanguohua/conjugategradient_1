/*
*function:conjugategradient_1.c
*purpose:ѧϰ�����ݶȷ�,��ʵ�֡�����¡ (2007). ����ģ���������ʵ��, ��ѧ��ҵ������.����Ӧϰ�����ϰ
*version:v1.0
*date:15:13 2016/2/6
*author:dgh
*note:��Ȼ����ļ���Ϊ��file.c��������ԭ�� ��conjugategradient_1.c���������ļ������һ��������Ҫ�ǳ���ѧϰunderstand�������������Ŀ��

*/

#include <stdio.h>
#include <math.h>
/*----------------prototype---------------------*/
/*--------���ɺ���------------------------*/
double fun3d(double p[]);						//Ŀ�꺯�����ʽ
void gradient(double p[], double g[]);		//Ŀ�꺯�����ݶ�
double conjugateGradient(double p[], double pmin[]);	// �����ݶȷ�
double steepDescend(double p[], double h[], double pmin[]);	//�����½���ʵ��
void goldenResearch3d(double a[], double b[], double xmin[], double *zmin); 	//ĳһά�����ϵĻƽ�ָ������㷨

/*--------֦Ҷ����------------------------*/
double getGamma(double gp0[], double gp1[]);	//����P-R�㷨, ����gamma
double getMod(double gp[]);					//���㵱ǰ�Ա����������ݶȵ�ģ


void test(void);			//���Ժ���

int main(void)
{
	test();
	// printf("%lf\n", sqrt(-1));		//��ʱ�ᱨ��-1.#IND00

	
	return 0;
}


/*----------------����ʵ��---------------------*/
//���Ժ���
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


// �����ݶȷ�
/*
*����˵����
double p[], ��ʼ����λ�ã��Ա���������
double pmin[], ��Сֵ����Ӧ���Ա�������
����ֵ����Сֵ
*/
double conjugateGradient(double p[], double pmin[])
{
	double gp[2];
	gradient(p, gp);				//�õ��Ա�������p���ĸ��ݶ�
	
	double p0[2], p1[2]; 			//��ǰ, ��һ���Ա�������
	double gp0[2], gp1[2];			//��ǰ, ��һ���Ա��������ĸ��ݶ�
	double h0[2], h1[2];			//��ǰ, ��һ������������	
	//��eps��ֵ�ǳ�Сʱ����1e-5ʱ���ᵼ��-1.#1ND00�������������Ҫ�Ƿ�����steepDescend()�����У��Ǹ����� lamda���Ϊ��ֵ����ʱʹ��sqrt()�����ͻ�������⣬
	double eps = 1e-3;				//������׼
	
	int i, j;
	double modgp = 0.0;				//���ݶȵ�ģ
	double gamma = 0.0;				//�����ݶȷ��е�һ��ϵ��
	double minz = 0.0;				//ĳһ���������ϵ���Сֵ
	
	//��ʼֵ
	for (i=0; i<2; i++)
	{
		p0[i] = p[i];
		gp0[i] = gp[i];
		h0[i] = gp[i];
	}
	modgp = getMod(gp0);				//�õ��Ա�������p0���ĸ��ݶȵ�ģ
	j = 0;
	//���ݶȵ�ģ����ĳһ��ֵʱ����������
	while (modgp >eps)
	{
		//�ı�p1, gp1, h1��ֵ
		minz = steepDescend(p0, h0, p1);		//�õ�h0��������Сֵ��Ӧ���Ա�������
		gradient(p1, gp1);				//�õ��Ա�������p1�����ݶ�
		printf("j=%d\t p1=(%lf, %lf), minz = %lf\n", j, p1[0], p1[1], minz);
		//����P-R�㷨, ����gamma
		gamma = getGamma(gp0, gp1);
		//������һ����������
		for (i=0; i<2; i++)
		{
			h1[i] = gp1[i] + gamma * h0[i];
		}
		
		// ����h0, gp0, p0
		for (i=0; i<2; i++)
		{
			h0[i] = h1[i];
			gp0[i] = gp1[i];
			p0[i] = p1[i];			
		}
		//����ѭ���ж�����
		modgp = getMod(gp0);				//�õ��Ա�������p0���ĸ��ݶȵ�ģ
		j++;
	}
	//������Сֵ ��Ӧ���Ա�������
	for (i=0; i<2; i++)
	{
		pmin[i] = p1[i];
	}
	
	return minz;
	
}
//�����½���ʵ��
/*
*����˵����
double p[], ��ʼ�Ա������� 
double h[], ��ʼ��������
double pmin[], ����������ϣ��ҵ�����Сֵ��Ӧ���Ա�������
����ֵ����������ϣ��ҵ�����Сֵ
*/
double steepDescend(double p[], double h[], double pmin[])
{
	int i, j;
	double xmin[2], zmin;		//��Сֵ ��Ӧ���Ա�������Сֵ 
	double ztemp, z0, z1;		//��һ������ǰ������һ����Ӧ����
	double ptemp[2], p0[2], p1[2];			//��һ��, ��ǰ������һ�����Ա���
	//��������������ȷ���أ���Ҫ��������ķ����ǵ�λ��������������
	//����ǵ�λ����������ȡΪ1���������ܻ����������Ĵ���
	//����ǵ�λ����������ȡһ����������ģ�ĵ���������ģ�����ļ���֮����
	double lamda;				//һά��������
	for (i=0; i<2; i++)
	{
		lamda += h[i] * h[i];
	}
	// printf("h = (%lf, %lf), lamda = %lf\n", h[0], h[1], lamda);
	//������������ļ�������䵱lamdaΪ��С����ֵʱ��
	// lamda = 1.0/sqrt(lamda);	//��ˣ���������������൱���ظ��ݶȵ�λ������������
	lamda = sqrt(lamda)/lamda;
	// printf("lamda = %lf\n", lamda);
	
	//��ʼ�����Ա�����Ӧ����
	j=0;
	z1 = ztemp = z0 = fun3d(p);
	for (i=0; i<2; i++)
	{
		p0[i] = p[i];
		p1[i] = p[i];
		ptemp[i] = p[i];
	}
	//��ĳһά�����������������f(x_k)<f(x_k+1),
	while (z1<=z0)
	{		
		for (i=0; i<2; i++)		//��¼��һ������ǰ�����Ա���
		{
			ptemp[i] = p0[i];
			p0[i] = p1[i];
		}
		ztemp = z0;				//��¼��һ������ǰ����Ӧ����
		z0 = z1;				
		j++;					//������������
		for (i=0; i<2; i++)		//��һ�������Ա���
		{
			//note: p1�������ʼλ�û����ϲ��� �ģ�������һ����
			//���Լ���p1[i]ʱ��Ҫ�õ�p[i]��
			//��Ҳ��֮ǰû��ע�⵽�ģ���������
			p1[i] = p[i] + j * lamda * h[i];
		}
		z1 = fun3d(p1);			//��һ������Ӧ����
	}

	//���ûƽ�ָ
	if (j == 1)		//��ʾ��һ������ʱ���ͳ���z1>z0�����
	{
		goldenResearch3d(p0, p1, xmin, &zmin);
	}
	else
	{
		goldenResearch3d(ptemp, p1, xmin, &zmin);
	}
	
	//������Сֵ��Ӧ���Ա���
	for (i=0; i<2; i++)
	{
		pmin[i] = xmin[i];
	}
	//������Сֵ 
	return zmin;
}

//ĳһά�����ϵĻƽ�ָ������㷨
/*
double a, Ϊ�Ա�����������������
double b, Ϊ�Ա���������������ұ�
double xmin[], ��Сֵ ��Ӧ��x�Ա�������
double *zmin  ��Сֵ ��Ӧ��yֵ
*/
void goldenResearch3d(double a[], double b[], double xmin[], double *zmin)
{
	int i;
	double tol[2];			//���䳤��
	for (i=0; i<2; i++)
	{
		tol[i] = b[i] - a[i];
	}
	
	double eps = 1e-3;			//������С����
	double golden = 0.618;		//�ƽ�ָ��
	int k = 0;
	
	// ����������Ա����������Ӧ��zֵ
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
	//������
	for (i=0; i<2; i++)
	{
		xmin[i] = (a[i] + b[i])/2;
	}	
	(*zmin) = fun3d(xmin);
}



//����P-R�㷨, ����gamma
/*
*����˵����
double gp0, �Ա�������p0λ�ô��ĸ��ݶ� 
double gp1,	�Ա�������p1λ�ô��ĸ��ݶ�
����һ��ϵ��
*/
double getGamma(double gp0[], double gp1[])
{
	//����P-R�㷨, ����gamma
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

//���㵱ǰ�Ա��������ĸ��ݶȵ�ģ
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



//Ŀ�꺯�����ʽ
//ע���޸�Ŀ�꺯��ʱ��ҲҪͬʱ�����ݶȺ��������޸�
/*
*����˵����
double p[], Ϊһ���Ա�������
����ֵΪ������Ӧ������
*/
double fun3d(double p[])
{
	double x, y, z;
	x = p[0];
	y = p[1];
	z = x * x + 2 * y * y;
	return z;
}

//Ŀ�꺯���ĸ��ݶ�
//ע���޸�Ŀ�꺯��ʱ��ҲҪͬʱ�����ݶȺ��������޸�
/*
*����˵����
double p[], Ϊһ���Ա������� 
double g[], ����ָ����������ֵ�� �������Ա����������������ĸ��ݶ�������
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
�밴���������. . .

*/
