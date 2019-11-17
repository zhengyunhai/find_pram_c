#include<NTL/ZZ.h>
#include<iostream>
#include<fstream>
#include<NTL/matrix.h> 
#include<NTL/mat_ZZ.h>
#include<NTL/vec_ZZ.h>
#include<NTL/vector.h>
#include<NTL/LLL.h>
#include<NTL/tools.h>
#include<NTL/RR.h>
#include <NTL/xdouble.h>
#include <NTL/quad_float.h>
#include<NTL/mat_RR.h>
#include<NTL/vec_RR.h>
#include<cstdlib>
#define LEN 10000
NTL_CLIENT

ZZ round_ZZ(ZZ a,ZZ b)//输出最接近a/b的整数 //没用到 
{
	ZZ d,r1,r2,res;
	d=operator/(a,b);
	r1=operator-(a,operator*(d,b));
	r2=operator-(operator*(d+1,b),a);
	if(operator>=(r1,r2)) res=d+1;
	else res=d;
	
	return res;
 } 

vec_ZZ close_vec_in2(mat_ZZ mat_c,vec_ZZ b)//m个向量中距离b最近的向量 
{
	vec_ZZ min,b_c;
	int i,l,m,n;
	m=mat_c.NumRows();//行数
	n=mat_c.NumCols();//列数,mat_c的行列数不相等 
	min.SetLength(n);
	b_c.SetLength(n); 
	
	ZZ mim,d2_i;
	mim=ZZ(0);//记录长度 
	//初始化 
	sub(b_c,b,mat_c(1));//向量减法 
	for(l=1;l<=n;l++)
	{
		mim=operator+(mim,sqr(b_c[l-1])); 
		
	}
		
	VectorCopy(min,mat_c(1),n); 
	
	for(i=2;i<=m;i++)
	{
		clear(b_c);
		d2_i=ZZ(0);
		sub(b_c,b,mat_c(i));//向量减法 
		for(l=1;l<=n;l++)
		{
			d2_i=operator+(d2_i,sqr(b_c[l-1])); 
		}	
	
		if(operator<(d2_i,mim)) 
		{
		VectorCopy(min,mat_c(i),n);
		mim=d2_i;
		}
	}
	return min;
}

int intcount(ZZ lb,ZZ ub)//范围内整数计数 
{
	
	int i=0;
	while(operator<=(operator+(lb,i),ub))
		{
			i++;
		}
		return i;
}

void close_vec(vec_ZZ& c,mat_ZZ& B,int k)
{
	RR mu1,mu2;
	vec_RR u;//存储正交化系数 
	vec_RR b1_,b2_,b3_;
	mat_RR BR;
	int m;
	m=B.NumRows();
	BR.SetDims(m,m);
	
	u.SetLength(4);
	b1_.SetLength(m); 
	b2_.SetLength(m); 
	b3_.SetLength(m); 
	
	int x,y;
	for(x=1;x<=m;x++)//每次都要把B转RR 
	{
		for(y=1;y<=m;y++)
		{
			BR(x,y)=MakeRR(B(x,y),0);
		}
	}
	
		//计算正交化，用向量序列存储之
		VectorCopy(b1_,BR(1),m);//b1*
		
		InnerProduct(mu1,BR(2),b1_);//<b2,b1*>
		InnerProduct(mu2,b1_,b1_);//<b1*,b1*>
		div(u[1],mu1,mu2);//u[1]=<b2,b1*>/<b1*,b1*>=u21
		
		b2_=operator-(BR(2),operator*(u[1],b1_));//b2*

		
		
		
		
		
	if(k==2)
	{
	//	cout<<"uk2="<<u<<endl;
		
		mul(c,RoundToZZ(u[1]),B(1));
	 }
	 
	if(k==3)
	{
		
		InnerProduct(mu1,BR(3),b1_);//<b3,b1*>
		InnerProduct(mu2,b1_,b1_);//<b1*,b1*> 
		 
		div(u[2],mu1,mu2);//u[2]=<b3,b1>/<b1,b1>=u31
		
		InnerProduct(mu1,BR(3),b2_);//<b3,b2*>
		InnerProduct(mu2,b2_,b2_);//<b2*,b2*> 
				
		div(u[3],mu1,mu2);//u[3]=<b3,b2*>/<b2*,b2*>=u32
		
	
		//u[2]=u[2]-u[3]*u[1]
		
		b3_=operator-(BR(3),operator-(operator*(u[2],b1_),operator*(u[3],b2_)));
		
			u[2]=operator-(u[2],operator*(u[3],u[1]));
		
		
	vec_ZZ c1,c2,c3,c4;

	vec_ZZ c_,temp,z;
	vec_RR tempr;
	RR d,normb2_,ct;
	ZZ lb,ub,tempd; 
	int num,j,i;
	temp.SetLength(m);
	tempr.SetLength(m);	
	c_.SetLength(m);
		
	c1.SetLength(m);
	c2.SetLength(m);
	c3.SetLength(m);
	c4.SetLength(m); 
	
	mat_ZZ mat_c,mat_c_;
	mat_c.SetDims(4,m);
	
	//cout<<"1="<<mat_c<<endl;
	
		c1=operator+(operator*(FloorToZZ(u[2])+1,B(1)),operator*(FloorToZZ(u[3])+1,B(2)));
		c2=operator+(operator*(FloorToZZ(u[2])+1,B(1)),operator*(FloorToZZ(u[3]),B(2)));
		c3=operator+(operator*(FloorToZZ(u[2]),B(1)),operator*(FloorToZZ(u[3])+1,B(2)));
		c4=operator+(operator*(FloorToZZ(u[2]),B(1)),operator*(FloorToZZ(u[3]),B(2)));
		
		VectorCopy(mat_c(1),c1,m);
		VectorCopy(mat_c(2),c2,m);
		VectorCopy(mat_c(3),c3,m);
		VectorCopy(mat_c(4),c4,m);
	
		c=close_vec_in2(mat_c,B(3));//返回ci四个向量中与b3最近的一个 
/*		
		temp=operator-(B(3),c);
		for(i=0;i<m;i++)
		{
			tempd=operator+(tempd,sqr(temp[i]));
		}
		
		d=MakeRR(tempd,0);
		d=sqrt(d);//d=||b3-c||
		for(i=0;i<m;i++)
		{
			normb2_=operator+(normb2_,sqr(b2_[i]));
		}
		
		normb2_=sqrt(normb2_);
		lb=CeilToZZ(operator-(u[3],operator/(d,normb2_)));
		ub=FloorToZZ(operator+(u[3],operator/(d,normb2_)));
	//	cout<<"lb...="<<lb<<" "<<ub<<endl;
		num=intcount(lb,ub);
	//	cout<<"num="<<num<<endl;
		mat_c_.SetDims(num+1,m);
		z.SetLength(num);
		
		i=0;	
		while(operator<=(operator+(lb,i),ub))
		{
			z[i]=operator+(lb,i);
			i++;
		}
		VectorCopy(mat_c_(1),c,m);
		
		
			for(j=0;j<num;j++)
		{
			temp=operator-(B(3),operator*(z[j],B(2)));

			//计算t3参与的正交化系数，先将t3转为RR型 
			for(i=0;i<m;i++)
			{
				tempr[i]=MakeRR(temp[i],0);
			}
			InnerProduct(mu1,tempr,b1_);//<t3i,b1*>
			InnerProduct(mu2,b1_,b1_);//<b2*,b2*>
			div(ct,mu1,mu2);//ct1=<t3i,b1*>/<b2*,b2*>
			
			mat_c_(j+2)=operator+(operator*(RoundToZZ(ct),B(1)),operator*(z[j],B(2)));
			
		}
		//	cout<<"mat_c_=..."<<endl;
			clear(c);
			c=close_vec_in2(mat_c_,B(3));//t3i
*/			

	 } 
	
	if(k==4)//m=4不用改了？ 
	{
	//	cout<<BR(3)<<endl<<b1_<<endl;
		InnerProduct(mu1,BR(3),b1_);//<b3,b1>
		InnerProduct(mu2,BR(1),b1_);//<b1,b1>
		div(u[2],mu1,mu2);//u[2]=<b3,b1*>/<b1*,b1*>=u31
	//	cout<<mu1<<" "<<mu2<<endl;
		InnerProduct(mu1,BR(3),b2_);//<b3,b2>
		InnerProduct(mu2,BR(2),b2_);//<b2,b2>
		div(u[3],mu1,mu2);//u[1]=<b3,b2*>/<b2*,b2*>=u32
	//	cout<<mu1<<" "<<mu2<<endl;
		b3_=operator-(BR(3),operator-(operator*(u[2],b1_),operator*(u[3],b2_)));
		
		InnerProduct(mu1,BR(4),b3_);//<b4,b3*>
		InnerProduct(mu2,b3_,b3_);//<b3*,b3*> 
	//	cout<<mu1<<" "<<mu2<<endl;
		div(u[0],mu1,mu2);//u[0]=<b4,b3*>/<b3*,b3*>=u43=c3
	//	cout<<"uk4="<<u<<endl;
		
		ZZ z3,tempd;
		int num;
		RR d,normb3_,ct1,ct2;
		ZZ lb,ub,ui;
		int i,j;
		vec_ZZ temp,s1,s2,s3,s4,t3z,temp_c,v3,z;
		mat_ZZ mat_s,mat_v;
		vec_RR t3,ut3;
		mat_s.SetDims(4,4); 

		temp.SetLength(4); 
		t3.SetLength(4);
		t3z.SetLength(4);
		v3.SetLength(4);
	
		ut3.SetLength(2);//存储正交化系数 
		z3=RoundToZZ(u[0]);
		
		temp=operator-(B(4),operator*(z3,B(3)));
		VectorCopy(t3z,temp,4);//t3 in ZZ
		//计算t3参与的正交化系数，先将t3转为RR型 
		for(i=0;i<4;i++)
		{
			t3[i]=MakeRR(temp[i],0);//t3 in RR
		}
		
		InnerProduct(mu1,t3,b1_);//<t3,b1*>
		InnerProduct(mu2,b1_,b1_);//<b1*,b1*> 
		
		div(ut3[0],mu1,mu2);//ut3[0]=<t3,b1*>/<b1*,b1*>=ut31
		
		InnerProduct(mu1,t3,b2_);//<t3,b2*>
		InnerProduct(mu2,b2_,b2_);//<b2*,b2*> 
		
		div(ut3[1],mu1,mu2);//ut3[1]=<t3,b2*>/<b2*,b2*>=ut32
		
		ut3[0]=operator-(ut3[0],operator*(ut3[1],u[1]));
		
		s1=operator+(operator*(FloorToZZ(ut3[0])+1,B(1)),operator*(FloorToZZ(ut3[1])+1,B(2)));
		s2=operator+(operator*(FloorToZZ(ut3[0])+1,B(1)),operator*(FloorToZZ(ut3[1]),B(2)));
		s3=operator+(operator*(FloorToZZ(ut3[0]),B(1)),operator*(FloorToZZ(ut3[1])+1,B(2)));
		s4=operator+(operator*(FloorToZZ(ut3[0]),B(1)),operator*(FloorToZZ(ut3[1]),B(2)));
		VectorCopy(mat_s(1),s1,4);
		VectorCopy(mat_s(2),s2,4);
		VectorCopy(mat_s(3),s3,4);
		VectorCopy(mat_s(4),s4,4);
		
		temp_c=close_vec_in2(mat_s,t3z);//返回mat_s四个向量中与t3最近的一个 
		
		v3=operator+(operator*(z3,B(3)),temp_c);
		
		temp=operator-(B(4),v3);
		tempd=operator+(operator+(sqr(temp[0]),sqr(temp[1])),operator+(temp[2],sqr(temp[3])));
		d=MakeRR(tempd,0);
		d=sqrt(d);//d=||b4-v3||
		normb3_=operator+(operator+(sqr(b3_[0]),sqr(b3_[1])),operator+(b3_[2],sqr(b3_[3])));//||b3*||
		normb3_=sqrt(normb3_);
		lb=CeilToZZ(operator-(u[0],operator/(d,normb3_)));
		ub=FloorToZZ(operator+(u[0],operator/(d,normb3_)));
	//	cout<<"lb...="<<lb<<" "<<ub<<endl;
		num=intcount(lb,ub);//计数
	//	cout<<"num="<<num<<endl;
	//	if(num>LEN+1) {cout<<"数组范围不足，应该动态分配？"<<endl;exit(0); }
		mat_v.SetDims(num+1,4);
		z.SetLength(num);
		
		i=0;	
		while(operator<=(operator+(lb,i),ub))
		{
			z[i]=operator+(lb,i);
			i++;
		}
		//最后共有num+1个v3 
		VectorCopy(mat_v(1),v3,4);
		
		for(j=0;j<num;j++)
		{
			temp=operator-(B(4),operator*(z[j],B(3)));
			VectorCopy(t3z,temp,4);
			//计算t3参与的正交化系数，先将t3转为RR型 
			for(i=0;i<4;i++)
			{
				t3[i]=MakeRR(temp[i],0);
			}
			InnerProduct(mu1,t3,b1_);//<t3i,b1*>
			InnerProduct(mu2,b1_,b1_);//<b2*,b2*>
			div(ct1,mu1,mu2);//ct1=<t3i,b1*>/<b2*,b2*>
			
			InnerProduct(mu1,t3,b2_);//<t3i,b2*>
			InnerProduct(mu2,b2_,b2_);//<b2*,b2*>
			div(ct2,mu1,mu2);//ct2=<t3i,b2*>/<b2*,b2*>
			
			ct1=operator-(ct1,operator*(ct2,u[1]));
			
			s1=operator+(operator*(CeilToZZ(ct1),B(1)),operator*(CeilToZZ(ct2),B(2)));
			s2=operator+(operator*(CeilToZZ(ct1),B(1)),operator*(FloorToZZ(ct2),B(2)));
			s3=operator+(operator*(FloorToZZ(ct1),B(1)),operator*(CeilToZZ(ct2),B(2)));
			s4=operator+(operator*(FloorToZZ(ct1),B(1)),operator*(FloorToZZ(ct2),B(2)));
			VectorCopy(mat_s(1),s1,4);
			VectorCopy(mat_s(2),s2,4);
			VectorCopy(mat_s(3),s3,4);
			VectorCopy(mat_s(4),s4,4);
		//	cout<<"mat_sk4="<<endl<<mat_s<<endl;
			
			mat_v(j+2)=operator+(operator*(z[j],B(3)),close_vec_in2(mat_s,t3z));//t3i
			
		}
		

		
		//从num+1个向量即（num+1,4）的矩阵中找出距离 b4最近的向量
		 
		c=close_vec_in2(mat_v,B(4));//返回num+1个向量中距离b4最近的 
		
		
	
	 } 
}

int order(mat_ZZ& B,int k)//将b1..bk按照长度排序(b1..bk-1已经有序),输出排序后的B和原先bk排序后的位置k_.
//默认B为4*4 ,输入的矩阵B的b1...bk应该是约化过的?（norm递增? 
{
	vec_ZZ norm,tempv;
	norm.SetLength(k);
	tempv.SetLength(4);
	ZZ tempn;
	int i,l;
	int j=1; 
	int m;
	m=B.NumRows();
	for(i=1;i<=k;i++)//将b1...bk的范数平方存入向量norm中 
	{
		norm[i-1]=ZZ(0);
		for(l=1;l<=m;l++)//计算bi范数 
		{
			norm[i-1]=operator+(norm[i-1],sqr(B(i,l)));
		 } 
	
	}
	
	//排序
	while(j<k)
	{
		if(operator>(norm[k-j-1],norm[k-j]))
		{
			
		tempn=norm[k-j-1];//数组内交换位置 
		norm[k-j-1]=norm[k-j];
		norm[k-j]=tempn;
		

		VectorCopy(tempv,B(k-j),m);//向量交换位置 
		VectorCopy(B(k-j),B(k-j+1),m);
		VectorCopy(B(k-j+1),tempv,m);
		
		j++;
	 	} 
	 	else break;
	}

	clear(norm);
	clear(tempv);
	return k-j+1;
	
}

int order_all(mat_ZZ& B)//默认B为4*4 方阵,将B按照范数递增排序 
{
	vec_ZZ norm,tempv;
	ZZ tempn;
	int i,m;
	m=B.NumRows();

	for(i=2;i<=m;i++)//将b1...bk的范数平方存入向量norm中 
	{
		order(B,i);
	}

	
}


void ig_reduce(mat_ZZ& B)//默认L为4*4 ,B与BR为同一个矩阵，数据类型不同，我不会转化！他妈的 
//BR仅用于计算正交化系数？ 
{
	vec_ZZ c,b1,b2,b3,b4,temp;
	int k,k_,d;
	ZZ normb_k_2,normb_k1_2;//范数用平方比较 
	//cout<<normb_k_2<<normb_k1_2<<endl;
	int i,l,m;
	m=B.NumRows();
	if((B.NumRows()!=B.NumCols()))
	{
	cout<<"Can't invoke Greedy_Reduce,Input should be Square Matrix!"<<endl;}
	else{	
	
	order_all(B); //先将B排序为范数递增
//	cout<<"orderallB="<<endl<<B<<endl; 
	k=2;
	d=B.NumRows();
	//d=4;

	c.SetLength(d);

	
	
	while(k<=d)
	{
		normb_k_2=ZZ(0);
		normb_k1_2=ZZ(0);//初始化为0 
		
		//cout<<"orderB="<<endl<<B<<endl;
		
		//cout<<"k="<<k<<endl;
		//cout<<"接近bk的向量："<<endl; 
		
		close_vec(c,B,k);//输出的c是b1..bk-1中与bk最近的向量
		
	//	cout<<"bk="<<B(k)<<endl<<"c="<<c<<endl; 
 		sub(B(k),B(k),c);//更新bk
	//	cout<<"bk="<<B(k)<<endl;
	
		 //计算bk bk-1的范数（平方 
	
			for(l=1;l<=m;l++)
			{
				normb_k_2=operator+(normb_k_2,sqr(B(k,l)));
			 } 
		
	
		
			for(l=1;l<=m;l++)
			{
				normb_k1_2=operator+(normb_k1_2,sqr(B(k-1,l)));
			 } 
		
 		//cout<<"比较："<<normb_k_2<<" "<<normb_k1_2<<endl;
 		
		 if(operator<(normb_k_2,normb_k1_2)) 
 		{
 			k_=order(B,k);//将b1..bk按照长度排序(b1..bk-1已经有序),输出排序后的B和原先bk排序后的位置k_.
			 //测试猜想
			 //if(k_==1)k=2;
			 //else k=k_;
			 
			 k=k_+1; 
	//		 cout<<"k_="<<k_<<endl;			 
		 }
		else k++;
		
	}
//	cout<<"ig_reduced="<<endl<<B<<endl;
}
	
}
