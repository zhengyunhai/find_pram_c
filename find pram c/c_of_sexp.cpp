#include<NTL/ZZ.h>
#include<NTL/tools.h>
#include<time.h>
#include<iostream>
#include<fstream>
#include<NTL/matrix.h> 
#include<NTL/mat_ZZ.h>
#include<NTL/LLL.h>
using namespace std;
NTL_CLIENT
int order_all(mat_ZZ& B);
void ig_reduce(mat_ZZ& B);
int main(void)
{
	ZZ p,q,x1,x1_,x2,x2_,y1,y2,r,f1,f2,f3,f1_,f2_,f3_,X,Y,Z,z,det,det2;
	ZZ seed;
	ZZ temp,temp1,temp2;
//	RR Zr;
	float c; 
	int c12,temp3,zlen;
	int i,count;//循环轮数 
	int tag,tagq;//超出上界计数 
	ZZ upperbound_test_c,norml2_2,norml2_4,norml2_8,norml2_12,norml1_2;
	ZZ boundq2;
	int cc;
	
	mat_ZZ L;
	L.SetDims(3,3);
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
		if(cc==1) c=1.00/6;
		if(cc==2) c=1.00/3;
		if(cc==3) c=2.00/3;
		if(cc==4) c=1.00;
		if(cc==5) c=4.00/3;
		if(cc==6) c=2.00;
	
	cout<<"c="<<c<<endl;
	//XY
	zlen=ceil(78.8113-1.5*c);//Z的长度 
	RandomLen(X,zlen);
	Y=X;
	cout<<"X="<<X<<" "<<Y<<endl;
	c12=ceil(12*c);
	 cout<<"c12="<<c12<<endl;
	tag=0;	
	tagq=0;
	  
	seed=time(0);
	SetSeed(seed);
//	cout<<"输入循环轮数：";
//	cin>>count;	
	count=10000;
	i=0;
	while(i<count)
	{
		//生成x y
		RandomBnd(x1,q);
	 	RandomBnd(x2,q);
	 	RandomBnd(x1_,q);
	 	RandomBnd(x2_,q);
	 	RandomBnd(y1,q);
	 	RandomBnd(y2,q); 
	 	
	 	//r
	 	temp=operator-(x2,x1);//x2-x1
	 	XGCD(z,temp1,temp2,temp,q);//temp1=(x2-x1)^-1
	 	r=operator*(temp1,operator-(x2_,x1_));
	// 	cout<<"r="<<r<<endl;
	 	//fi
		 f1=operator-(operator*(r,x1),x1_);
		 f2=operator*(r,y1);
		 f3=-y2;
		 //fi_
		 XGCD(z,f1_,temp,f1,q);
		 f2_=operator*(f1_,f2)%q; 
		f3_=operator*(f1_,f3)%q;
	//	cout<<f1_<<" "<<f2_<<" "<<f3_<<endl;
		//格基 
		L(1,1)=ZZ(1);L(1,2)=operator*(f2_,X);L(1,3)=operator*(f3_,Y);
		L(2,1)=ZZ(0);L(2,2)=operator*(q,X);L(2,3)=ZZ(0);
		L(3,1)=ZZ(0);L(3,2)=ZZ(0);L(3,3)=operator*(q,Y);
	//	norml2_2=operator+(operator+(sqr(L(2,1)),sqr(L(2,2))),sqr(L(2,3)));//b2范数平方
	//	norml1_2=operator+(operator+(sqr(L(1,1)),sqr(L(1,2))),sqr(L(1,3)));//b1范数平方
	//	cout<<norml1_2<<endl<<norml2_2<<endl;
		order_all(L);
	//	cout<<"L orignal="<<endl<<L<<endl;
		//约化
	//	mat_ZZ L_;
	//	L_.SetDims(3,3);
	//	for(i=1;i<=3;i++)
	//	{
	//		VectorCopy(L_(i),L(i),3);
	//	}
		
	
	//	LLL(det,L_,0);
	//	cout<<"LLL="<<endl<<L_<<endl;
		//LLL(det2,L,0);
		//order_all(L);
		ig_reduce(L);
	//	cout<<"ig="<<endl<<L<<endl;
		
		//b范数
		norml2_2=operator+(operator+(sqr(L(2,1)),sqr(L(2,2))),sqr(L(2,3)));//b2范数平方
		norml1_2=operator+(operator+(sqr(L(1,1)),sqr(L(1,2))),sqr(L(1,3)));//b1范数平方
		norml2_4=sqr(norml2_2);
		norml2_8=sqr(norml2_4);
		norml2_12=power(norml2_2,6);
		//norml2_12=operator*(norml2_4,norml2_8);
		//计算上界
		temp3=pow(2,c12);
		det=operator*(sqr(q),sqr(X));//det=q^2*X^2
		temp1=sqr(sqr(det));
		upperbound_test_c=operator*(temp3,temp1);//2^c*detL^1/3
		
	//	cout<<"范数与上界："<<endl<<norml2_12<<endl<<upperbound_test_c<<endl;
		//比较并计数  
		//if(operator>(norml2_2,boundq2)) tagq++;
	 	if(operator>(norml2_12,upperbound_test_c)) tag++;
	 	i++;
	 	
	}
	cout<<"c="<<c<<"构造"<<count<<"次不同的格基进行约化，有"
	  	  <<tag<<"次超出了2^c*detL^1/3的上界,有"<<tagq<<"次超出q/2"<<endl;
			 
	outf.open("compareresult_sexp.txt",ios::app);
	outf<<"c="<<c<<"构造"<<count<<"次不同的格基进行约化，有"
	  	<<tag<<"次超出了2^c*detL^1/3的上界,有"<<tagq<<"次超出q/2"<<endl;
	outf.close();
	
}
	
}
