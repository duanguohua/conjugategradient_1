/*
*function:onjugategradient_1.c
*purpose:ѧϰ�����ݶȷ�,��ʵ�֡�����¡ (2007). ����ģ���������ʵ��, ��ѧ��ҵ������.����Ӧϰ�����ϰ
*version:v1.0
*date:15:13 2016/2/6
*author:dgh
*/

#include <stdio.h>
#include "function.h"


void test(void);			//���Ժ���

int main(void)
{
	test();
	// printf("%lf\n", sqrt(-1));		//��ʱ�ᱨ��-1.#IND00

	
	return 0;
}
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




/*
*answer:
j=0      p1=(3.999869, -1.000261), minz = 18.000000
j=1      p1=(-0.000869, -0.000100), minz = 0.000001
j=2      p1=(-0.000551, 0.000536), minz = 0.000001
j=3      p1=(-0.000033, 0.000213), minz = 0.000000
x=(-0.000033, 0.000213), z = 0.000000
�밴���������. . .

*/