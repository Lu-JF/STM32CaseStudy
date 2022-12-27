#include "filtfilt.h"
#include "stdio.h"
/* x     */
void iir_filt(double *x,int xlen,double *a,double *b,int nfilt,double *y)
{
	int i,j;
	double c,d;
//	double Head[8];
//	for(i=0;i<4;i++)
//	{
//		Head[i]=x[0]-(4-i)*(x[0]-x[1]);
//	}
//	for(i=4;i<8;i++)
//	{
//		Head[i]=x[i-4];
//	}
//	for(i=0;i<5;i++)
//	{
//		y[i]=0.0;
//		c=0;
//		d=0;
//		for(j=0;j<5;j++)
//		{
//			d=b[j]*Head[4-j+i];
//			if(i-j<0)
//			{
//				d=0;
//			}
//			else 
//			{
//				c=-a[j+1]*y[i-j];
//			}
//			y[i]=c+d;
//		}
//		y[i]=y[i]+b[0]*x[i];
//	}
	y[0]=b[0]*x[0];
	for(i=1;i<xlen;i++)
	{
		y[i]=0.0;
		c=0;
		d=0;
		for(j=1;j<nfilt;j++)
		{
			if(i-j<0)
			{
				c=0;
				d=0;
			}
			else 
			{
				c=-a[j]*y[i-j];
				d=b[j]*x[i-j];
			}
			y[i]=c+d+y[i];
		}
		y[i]=y[i]+b[0]*x[i];
	}
	for(i=0;i<xlen;i++)
	{
		x[i]=y[xlen-i-1];
	}
}
//����λ�˲��㷨��Ϊ����iir�˲���
//x�������źţ�����ź�Ҳ��x
//xlen Ϊ�źų���
//a,bΪ������˹�˲�������
//nfiltΪa,b��������
//yΪ�м�����������洢����
void filtfilt(double *x,int xlen,double *a,double *b,int nfilt,double *y)
{
	iir_filt(x,xlen,a,b,nfilt,y);
	iir_filt(x,xlen,a,b,nfilt,y);
}



















