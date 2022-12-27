#include "hilbert.h"
#include "math.h"
void abs_hilbert(int n,double x[])
{
	double m[4096];
	int i;
	for(i=0;i<n;i++)
	{
		m[i]=x[i];
	}
	hilbert(n,x);
	for(i=0;i<n;i++)
	{
		x[i]=sqrt((x[i]*x[i]+m[i]*m[i]));
	}
}
//希尔伯特变换，n为长度,x为信号
void hilbert(int n,double x[])
{
	int i,n1,n2;
	double t;
	n1=n/2;
	n2=n1+1;
	fht(n,x);
	for (i=1;i<n1;i++)
	{
		t=x[i];
		x[i]=x[n-i];
		x[n-i]=t;
	}
	for(i=n2;i<n;i++)
		x[i]=-x[i];
	x[0]=0.0;
	x[n1]=0.0;
	fht(n,x);
	t=1.0/n;
	for(i=0;i<n;i++)
		x[i]*=t;
}
void fht(int n,double x[])
{
	int i,j,k,m,l1,l2,l3,l4,n1,n2,n4;
	double a,e,c,s,t,t1,t2;
	for (j=1,i=1;i<16;i++)
	{
		m=i;
		j=2*j;
		if(j==n) break;
	}
	n1=n-1;
	for (j=0,i=0;i<n1;i++)
	{
		if (i<j)
		{
			t=x[j];
			x[j]=x[i];
			x[i]=t;
		}
		k=n/2;
		while (k<(j+1))
		{
			j=j-k;
			k=k/2;			
		}
		j=j+k;
	}
	for (i=0;i<n;i+=2)
	{
		t=x[i];
		x[i]=t+x[i+1];
		x[i+1]=t-x[i+1];
	}
	n2=1;
	for (k=2;k<=m;k++)
	{
		n4=n2;
		n2=n4+n4;
		n1=n2+n2;
		e=6.28318530719586/n1;
		for(j=0;j<n;j+=n1)
		{
			l2=j+n2;
			l3=j+n4;
			l4=l2+n4;
			t=x[j];
			x[j]=t+x[l2];
			x[l2]=t-x[l2];
			t=x[l3];
			x[l3]=t+x[l4];
			x[l4]=t-x[l4];
			a=e;
			for (i=1;i<n4;i++)
			{
				l1=j+i;
				l2=j-i+n2;
				l3=l1+n2;
				l4=l2+n2;
				c=cos(a);
				s=sin(a);
				t1=x[l3]*c+x[l4]*s;
				t2=x[l3]*s-x[l4]*c;
				a=(i+1)*e;
				t=x[l1];
				x[l1]=t+t1;
				x[l3]=t-t1;
				t=x[l2];
				x[l2]=t+t2;
				x[l4]=t-t2;
			}
		}
			
	}
}
/*
		计算希尔伯特变换后的角度并存入x
		n为向量长度，x存储虚部


*/
void get_hilbert_angle(int n,double x[])
{
	double y[4096];
	int i;
	for(i=0;i<n;i++)
	{
		y[i]=x[i];
	}
	hilbert(n,x);
	for(i=0;i<n;i++)
	{
		if(y[i]>0)  
		{
			if(x[i]!=0)
			{
				x[i]=atan(x[i]/y[i]);  					//第一、四象限
			}
			else {x[i]=0;}										//正实轴
		}
		else if(y[i]<0)
		{
			if(x[i]>0)
			{
				x[i]=atan(x[i]/y[i])+3.1415926535;  //第二象限
			}
			else if(x[i]<0)
			{
				x[i]=atan(x[i]/y[i])-3.1415926535;  //第三象限
			}
			else {x[i]=3.1415926535;}              //负实轴
		}
		else if(y[i]==0)
		{
			if(x>0) {x[i]=3.1415926535/2;}					//正虚轴
			else if(x<0) {x[i]=-3.1415926535/2;}    //负虚轴
			else {x[i]=0;}
		}
	}
}
void unwrap(int n,double x[])
{
	int i;
	i=0;
	for(i=1;i<n;i++)
	{
		while(x[i-1]-x[i]>3.1415926535 || x[i-1]-x[i]<-3.1415926535)
		{
			if(x[i-1]-x[i]>3.1415926535)
			{
				x[i]=x[i]+6.283185307;
			}
			else 
			{
				x[i]=x[i]-6.283185307;
			}
		}
	}
}

//      x信号 y 角度  time中间变量用来记录原始时间点
void resample(int n,double x[],double y[],double time[])
{
	double m,l;
	int i=0,j=0;
	for(i=0;i<5120;i++)
	{
		time[i]=0;
	}
	for(i=0;i<4096;i++)
	{
		x[i]=x[i+512];
	}
	l=y[n-1]-y[0];
	l=l/n;
	m=y[0];
	i=0;
	//将角度信号等角度切割，并计算每个均分的角度点在角度信号上的时间(两坐标点线性插值)
	while(1)
	{
		if(y[j]<=m&&m<y[j+1])
		{
			time[i]=j/10000.0+0.0001*(m-y[j])/(y[j+1]-y[j]);//10000.0为采样频率
			m+=l;
			i++;
		}
		else
		{
			j=j+1;
			if(j==4095)
				break;
		}
	}
	i=0;
	j=0;
	m=time[0];
	//通过等角度时间点计算该时间点上的振动信号点(两坐标点线性插值)
	while(1)
	{
		if((j/10000.0)<=m&&m<((j+1)/10000.0))
		{
			y[i]=x[j]+(m-j/10000.0)*(x[j+1]-x[j])/0.0001;
			i++;
			m=time[i];
		}
		else
		{
			j=j+1;
			if(j==4095)
				break;
		}
	}
}
void dfft(int n,int k,double fr[],double fi[],double pr[],double pi[])
{
	int it,m,is,i,j,nv,l0;
	double p,q,s,vr,vi,poddr,poddi;
	for (it=0; it<=n-1; it++)  //将pr[0]和pi[0]循环赋值给fr[]和fi[]
	{ 
		m=it; 
		is=0;
		for(i=0; i<=k-1; i++)
		{ 
			j=m/2; 
			is=2*is+(m-2*j); 
			m=j;
		}
		fr[it]=pr[is]; 
		fi[it]=pi[is];
	}
	pr[0]=1.0; 
	pi[0]=0.0;
	p=6.283185306/(1.0*n);
	pr[1]=cos(p); //将w=e^-j2pi/n用欧拉公式表示
	pi[1]=-sin(p);
	for (i=2; i<=n-1; i++)  //计算pr[]
	{ 
		p=pr[i-1]*pr[1]; 
		q=pi[i-1]*pi[1];
		s=(pr[i-1]+pi[i-1])*(pr[1]+pi[1]);
		pr[i]=p-q; pi[i]=s-p-q;
	}
	for (it=0; it<=n-2; it=it+2)  
	{ 
		vr=fr[it]; 
		vi=fi[it];
		fr[it]=vr+fr[it+1]; 
		fi[it]=vi+fi[it+1];
		fr[it+1]=vr-fr[it+1]; 
		fi[it+1]=vi-fi[it+1];
	}
	m=n/2; 
	nv=2;
	for (l0=k-2; l0>=0; l0--) //蝴蝶操作
	{ 
	m=m/2; 
	nv=2*nv;
	for (it=0; it<=(m-1)*nv; it=it+nv)
		for (j=0; j<=(nv/2)-1; j++)
			{ 
				p=pr[m*j]*fr[it+j+nv/2];
				q=pi[m*j]*fi[it+j+nv/2];
				s=pr[m*j]+pi[m*j];
				s=s*(fr[it+j+nv/2]+fi[it+j+nv/2]);
				poddr=p-q; 
				poddi=s-p-q;
				fr[it+j+nv/2]=fr[it+j]-poddr;
				fi[it+j+nv/2]=fi[it+j]-poddi;
				fr[it+j]=fr[it+j]+poddr;
				fi[it+j]=fi[it+j]+poddi;
			}
	}
}








