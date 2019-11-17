#include<NTL/ZZ.h>
#include<NTL/tools.h>
#include<time.h>
#include<iostream>
#include<fstream>
#include<NTL/matrix.h> 
#include<NTL/mat_ZZ.h>
#include<NTL/LLL.h>
#include<NTL/vec_ZZ.h>
#include<NTL/vector.h>
#include<NTL/RR.h>
#include <NTL/xdouble.h>
#include <NTL/quad_float.h>
#include<NTL/mat_RR.h>
#include<NTL/vec_RR.h>
#include<cstdlib>
using namespace std;
NTL_CLIENT
int order_all(mat_ZZ& B);
void ig_reduce(mat_ZZ& B);
int main(void)
{
	ZZ p,q,x1,x2,y1,y2,f1,f2,f3,f4,f1_,f2_,f3_,f4_,X,Y,Z,z,det,det2;
	ZZ seed;
	ZZ temp,temp2;
	RR Zr;
	float c; 
	int c8,temp1,zlen;
	int i,count;//循环轮数 
	int tag,tagq;//超出上界计数 
	ZZ upperbound_test_c,norml3_2,norml3_4,norml3_8;
	ZZ boundq2;
	int cc;//测试每个c   count=1:6
	 

	mat_ZZ L;
	L.SetDims(4,4);
	ifstream inf;
	ofstream outf; 	
	//读入pq
	inf.open("p 1024 q 160 g.txt");
	inf>>p>>q;
	inf.close();	
	boundq2=sqr(q/2);

	//c 
//	cout<<"输入参数c(float):";//1/4,3/4,7/8,1,3/2,2
//	cin>>c;
	for(cc=1;cc<=6;cc++)
	{
		if(cc==1) c=0.25;
		if(cc==2) c=0.75;
		if(cc==3) c=0.875;
		if(cc==4) c=1;
		if(cc==5) c=1.5;
		if(cc==6) c=2;
	//c=2;
	cout<<"c="<<c<<endl;
	 c8=ceil(8*c);
	 cout<<"c8="<<c8<<endl;
	//确定XYZ
	X=ZZ(11);
	zlen=ceil(74.5406-2*c);//Z的长度 
	RandomLen(Z,zlen);
	Y=operator*(X,Z);
	//初始化计数标志
	tag=0;	
	tagq=0;
	  cout<<"Z="<<Z<<endl;
	seed=time(0);
	SetSeed(seed);
//	cout<<"输入循环轮数：";
//	cin>>count;
count=10000;
	 //进入循环 
	 i=0;
	 cout<<"计数中。。。"<<endl;
	 while(i<count)
	 {
	 //	cout<<i<<endl; 
	 	//生成若干xi yi
	 	RandomBnd(x1,q);
	 	RandomBnd(x2,q);
	 	RandomBnd(y1,q);
	 	RandomBnd(y2,q);
	 //	cout<<"xi="<<x1<<" "<<x2<<" "<<y1<<" "<<y2<<endl;
	 	//计算多项式系数
		f1=-x2;
		f2=x1;
		f3=y1;
		f4=-y2;
		XGCD(z,f1_,temp,f1,q); 
	//	cout<<"gcd="<<z<<endl;
		f2_=operator*(f1_,f2)%q;
		f3_=operator*(f1_,f3)%q;
		f4_=operator*(f1_,f4)%q;
	//	cout<<f2_<<" "<<f3_<<" "<<f4_<<endl;
		
		//构造格基
		 L(1,1)= ZZ(1);L(1,2)=operator*(f2_,X);L(1,3)=operator*(f3_,Y);L(1,4)=operator*(f4_,Z);
		 L(2,1)=ZZ(0);L(2,2)=operator*(q,X);L(3,1)=ZZ(0);L(4,1)=ZZ(0);
		 L(3,1)=ZZ(0);L(3,2)=ZZ(0);L(3,3)=operator*(q,Y);L(3,4)=ZZ(0);
		 L(4,1)=ZZ(0);L(4,2)=ZZ(0);L(4,3)=ZZ(0);L(4,4)=operator*(q,Z);
	//	 cout<<L<<endl;
		 
		//IG约化
	//	ig_reduce(L); 
		
		//LLL约化
		LLL(det2,L,0); 
		order_all(L);
	//	cout<<L<<endl;
	//	order_all(L);
	//	cout<<L<<endl;
		//计算b3范数
		norml3_2=operator+(operator+(sqr(L(3,1)),sqr(L(3,2))),operator+(sqr(L(3,3)),sqr(L(3,4))));//b3范数平方
		norml3_4=sqr(norml3_2);
		norml3_8=sqr(norml3_4);
		//计算取c时的上界
		
		temp1=pow(2,c8);
		det=operator*(operator*(power(q,3),X),operator*(Y,Z));
		upperbound_test_c=operator*(temp1,sqr(det));
	//upperbound_test_c=operator*(4,sqr(det));
		//比较并计数  
		if(operator>(norml3_2,boundq2)) tagq++;
	 	if(operator>(norml3_8,upperbound_test_c)) tag++;
	 	
	 //	cout<<norml3_8<<endl<<"<?"<<endl<<upperbound_test_c<<endl;
	 	i++;
	  } 
	  	cout<<"c="<<c<<"构造"<<count<<"次不同的格基进行约化，有"
	  	  <<tag<<"次超出了2^c*detL^1/4的上界,有"<<tagq<<"次超出q/2"<<endl;
			 
	outf.open("compareresult_exp.txt",ios::app);
	outf<<"c="<<c<<"构造"<<count<<"次不同的格基进行约化，有"
	  	<<tag<<"次超出了2^c*detL^1/4的上界,有"<<tagq<<"次超出q/2"<<endl;
	outf.close();

	}
	 
 } 
