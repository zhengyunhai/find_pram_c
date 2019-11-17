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
	int i,count;//ѭ������ 
	int tag,tagq;//�����Ͻ���� 
	ZZ upperbound_test_c,norml3_2,norml3_4,norml3_8;
	ZZ boundq2;
	int cc;//����ÿ��c   count=1:6
	 

	mat_ZZ L;
	L.SetDims(4,4);
	ifstream inf;
	ofstream outf; 	
	//����pq
	inf.open("p 1024 q 160 g.txt");
	inf>>p>>q;
	inf.close();	
	boundq2=sqr(q/2);

	//c 
//	cout<<"�������c(float):";//1/4,3/4,7/8,1,3/2,2
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
	//ȷ��XYZ
	X=ZZ(11);
	zlen=ceil(74.5406-2*c);//Z�ĳ��� 
	RandomLen(Z,zlen);
	Y=operator*(X,Z);
	//��ʼ��������־
	tag=0;	
	tagq=0;
	  cout<<"Z="<<Z<<endl;
	seed=time(0);
	SetSeed(seed);
//	cout<<"����ѭ��������";
//	cin>>count;
count=10000;
	 //����ѭ�� 
	 i=0;
	 cout<<"�����С�����"<<endl;
	 while(i<count)
	 {
	 //	cout<<i<<endl; 
	 	//��������xi yi
	 	RandomBnd(x1,q);
	 	RandomBnd(x2,q);
	 	RandomBnd(y1,q);
	 	RandomBnd(y2,q);
	 //	cout<<"xi="<<x1<<" "<<x2<<" "<<y1<<" "<<y2<<endl;
	 	//�������ʽϵ��
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
		
		//������
		 L(1,1)= ZZ(1);L(1,2)=operator*(f2_,X);L(1,3)=operator*(f3_,Y);L(1,4)=operator*(f4_,Z);
		 L(2,1)=ZZ(0);L(2,2)=operator*(q,X);L(3,1)=ZZ(0);L(4,1)=ZZ(0);
		 L(3,1)=ZZ(0);L(3,2)=ZZ(0);L(3,3)=operator*(q,Y);L(3,4)=ZZ(0);
		 L(4,1)=ZZ(0);L(4,2)=ZZ(0);L(4,3)=ZZ(0);L(4,4)=operator*(q,Z);
	//	 cout<<L<<endl;
		 
		//IGԼ��
	//	ig_reduce(L); 
		
		//LLLԼ��
		LLL(det2,L,0); 
		order_all(L);
	//	cout<<L<<endl;
	//	order_all(L);
	//	cout<<L<<endl;
		//����b3����
		norml3_2=operator+(operator+(sqr(L(3,1)),sqr(L(3,2))),operator+(sqr(L(3,3)),sqr(L(3,4))));//b3����ƽ��
		norml3_4=sqr(norml3_2);
		norml3_8=sqr(norml3_4);
		//����ȡcʱ���Ͻ�
		
		temp1=pow(2,c8);
		det=operator*(operator*(power(q,3),X),operator*(Y,Z));
		upperbound_test_c=operator*(temp1,sqr(det));
	//upperbound_test_c=operator*(4,sqr(det));
		//�Ƚϲ�����  
		if(operator>(norml3_2,boundq2)) tagq++;
	 	if(operator>(norml3_8,upperbound_test_c)) tag++;
	 	
	 //	cout<<norml3_8<<endl<<"<?"<<endl<<upperbound_test_c<<endl;
	 	i++;
	  } 
	  	cout<<"c="<<c<<"����"<<count<<"�β�ͬ�ĸ������Լ������"
	  	  <<tag<<"�γ�����2^c*detL^1/4���Ͻ�,��"<<tagq<<"�γ���q/2"<<endl;
			 
	outf.open("compareresult_exp.txt",ios::app);
	outf<<"c="<<c<<"����"<<count<<"�β�ͬ�ĸ������Լ������"
	  	<<tag<<"�γ�����2^c*detL^1/4���Ͻ�,��"<<tagq<<"�γ���q/2"<<endl;
	outf.close();

	}
	 
 } 
