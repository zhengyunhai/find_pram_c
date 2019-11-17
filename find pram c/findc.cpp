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
/*������� ki,ti,r,���� u,a�������β�ֵ� xi,yi,�õ�����ʽ Ax+By+Cz+D ��ϵ�� ABCD
ϵ����D��õ� A'�ȣ��̶�����õ����Ԫ��ֵ 
*/

//����ֻ�޸�Z��ȡֵ��ʹ֮������Ȼ���ٿ������������ĳ��� 
int main(void)
{
	ZZ k1,k2,k3,k4,k5,k6,p,q,t1,t2,r,a,u,y1,y2,seed,x1,x2,gk1,gk2,gk3,gk4,gk5,gk6,det2,det,h,z;
	ZZ k1t1,k5a,temp,temp1,temp2,temp3,temp4,temp5,k2t2,k6ra;
	ZZ A,B,C,D,X,Y,Z,A_,B_,C_,D_;
	ZZ upperbound_test1_4,upperbound_test3_8,upperbound_test1_2,upperbound_test5_8,upperbound_test3_4,upperbound_test7_8,upperbound_test1,norml3_2,norml3_4,norml3_8,norml2_2,norml2_4,norml2_8;
	ZZ boundq2;
	int tag3_1_4,tag3_3_8,tag3_1_2,tag3_5_8,tag3_3_4,tag3_7_8,tagq_3,tag2_1_4,tag2_3_8,tag2_1_2,tag2_5_8,tag2_3_4,tag2_7_8,tagq_2,tag3_1,tag2_1,tag3_3_2,tag2_3_2,tag3_2;//�����Ͻ���� 
	int i,count;
	int len_Z,len_q;
	mat_ZZ L;
	
	L.SetDims(4,4);
	
	ifstream inf;
	ofstream outf;
	
	inf.open("p and q.txt");//p,q,,ä����ȡֵ�̶�����r,ti,u,a,��ȡֵ�仯 
	inf>>p>>q;
	inf.close();
	inf.open("pairs.txt");
	inf>>k1>>gk1>>k2>>gk2>>k3>>gk3>>k4>>gk4>>k5>>gk5>>k6>>gk6;
	inf.close();
//����r,ti
	inf.open("ti and r.txt");
	inf>>t1>>t2>>r;
	inf.close();
//����u,a 
	inf.open("u and a.txt");
	inf>>u>>a;
	inf.close();
//����xi yi
	inf.open("xi and yi.txt");
	inf>>x1>>x2>>y1>>y2;
	inf.close();

 
	boundq2=sqr(q/2);

	seed=time(0);
	SetSeed(seed);

	len_q=NumBits(q);
	
	/*�ⲿ�ּ���xi yi*/	
	X=ZZ(11);
	
	cout<<"p="<<p<<endl<<"q="<<q<<endl<<"ti and r="<<t1<<endl<<t2<<endl<<r<<endl<<"u="<<u<<endl<<"a="<<a<<endl<<"xi="<<x1<<endl<<x2<<"yi="<<y1<<y2<<endl;
	cout<<"boundq2="<<boundq2<<endl;
	
	/*================ѭ����ʼ=============*/
	
//	len_Z=NumBits(t2)+1;
	len_Z=64;
	while(len_Z<len_q)
	{			
	i=0;
	count=10000;
cout<<"len_Z="<<len_Z<<endl;
//getchar();
cout<<"finding......"<<endl; 

	tag3_1_4=0;//������ʼ�� b3
	tag3_3_8=0;
	tag3_1_2=0;
	tag3_5_8=0;
	tag3_3_4=0;
	tag3_7_8=0;
	tag3_1=0;
	tag3_3_2=0;
	tag3_2=0;

	
	tagq_3=0;

	
	while(i<count)
	{
		i++;
	
	//cout<<i<<endl;
	
	 /*����Z��ֵ*/ 
	//do{
		RandomLen(Z,len_Z);
		Y=operator*(X,Z);
//	}while(operator>(t2,Z)||operator>(operator*(r,t1),Y));
	/*over*/ 
	

	/*-----�������ʽϵ��ABCD,�¶���ʽϵ��A'B'C'D'----- */
	
	A=x1;
	B=y1;
	C=-y2;
	D=-x2;
	
	

	/*ϵ���������*/
	
	
	B_=operator*(InvMod(A,q),B);//��A���� ww
	C_=operator*(InvMod(A,q),C);
	D_=operator*(InvMod(A,q),D);
		
	//cout<<"B_="<<endl<<B_<<endl<<C_<<endl<<D_<<endl;
	
	/*����������L*/ 

	L(1,1)=X;     L(1,2)=operator*(B_,Y); L(1,3)=operator*(C_,Z); L(1,4)=D_;
	L(2,1)=ZZ(0); L(2,2)=operator*(q,Y);  L(2,3)=ZZ(0);           L(2,4)=ZZ(0);
	L(3,1)=ZZ(0); L(3,2)=ZZ(0);           L(3,3)=operator*(q,Z);  L(3,4)=ZZ(0);
	L(4,1)=ZZ(0); L(4,2)=ZZ(0);           L(4,3)=ZZ(0);           L(4,4)=q;	
	
	//cout<<L<<endl;
	/*Լ��L*/
	LLL(det2,L,0);
	
	/*����b3������ƽ����²���Ͻ�*/
	
	det=operator*(operator*(power(q,3),X),operator*(Y,Z));//����ʽ 
	norml2_2=operator+(operator+(sqr(L(2,1)),sqr(L(2,2))),operator+(sqr(L(2,3)),sqr(L(2,4))));//b2����ƽ�� 
	norml3_2=operator+(operator+(sqr(L(3,1)),sqr(L(3,2))),operator+(sqr(L(3,3)),sqr(L(3,4))));//b3����ƽ�� 
	 norml3_4=sqr(norml3_2);
	 norml3_8=sqr(norml3_4);
	 norml2_4=sqr(norml2_2);
	 norml2_8=sqr(norml2_4);
	 	
	 upperbound_test1_4=operator*(2,det);//�µ��Ͻ�c=1/4 	 
	 upperbound_test3_8=operator*(8,sqr(det));//�µ��Ͻ�c=3/8 
	 upperbound_test1_2=operator*(4,det);//�µ��Ͻ�c=1/2
	 upperbound_test5_8=operator*(32,sqr(det));//�µ��Ͻ�c=5/8
	 upperbound_test3_4=operator*(8,det);//�µ��Ͻ�c=3/4 
	 upperbound_test7_8=operator*(128,sqr(det));//�µ��Ͻ�c=7/8 
	 upperbound_test1=operator*(16,det);//�µ��Ͻ�c=1
	 
	
	 if(operator>(norml3_2,boundq2)) {
	 tagq_3++;//cout<<"norml3_2="<<norml3_2<<endl<<"boundq2="<<boundq2<<endl;
}
	 if(operator>(norml2_2,boundq2)) tagq_2++;

	 
	 /*�Ĵη��ıȽ�*/
	 if(operator>(norml3_4,upperbound_test1_4)) tag3_1_4++; //b3
	 if(operator>(norml3_8,upperbound_test3_8)) tag3_3_8++; 
	 if(operator>(norml3_4,upperbound_test1_2)) tag3_1_2++; 
	 if(operator>(norml3_8,upperbound_test5_8)) tag3_5_8++; 
	 if(operator>(norml3_4,upperbound_test3_4)) tag3_3_4++; 
	 if(operator>(norml3_8,upperbound_test7_8)) tag3_7_8++;
	 if(operator>(norml3_4,upperbound_test1)) tag3_1++;
	 if(operator>(norml3_4,operator*(64,det))) tag3_3_2++;
	 if(operator>(norml3_4,operator*(256,det))) tag3_2++;
	 
	 if(operator>(norml2_4,upperbound_test1_4)) tag2_1_4++; //b2
	 if(operator>(norml2_8,upperbound_test3_8)) tag2_3_8++; 
	 if(operator>(norml2_4,upperbound_test1_2)) tag2_1_2++; 
	 if(operator>(norml2_8,upperbound_test5_8)) tag2_5_8++; 
	 if(operator>(norml2_4,upperbound_test3_4)) tag2_3_4++; 
	 if(operator>(norml2_8,upperbound_test7_8)) tag2_7_8++;
	 if(operator>(norml2_4,upperbound_test1)) tag2_1++;
	}
//	cout<<"tag="<<tag<<endl; 
	outf.open("result.txt",ios::app);
	outf<<"Zȡ"<<len_Z<<"bit������"<<count<<"��Լ����u,a,r,ti,Z�������ȡֵ���ظ���������Խ�С��������У�b3������"<<endl<<"��"<<tag3_1_4<<"�γ�����c=1/4���Ͻ�,"
	<<"��"<<tag3_3_8<<"�γ�����c=3/8���Ͻ�,"
	<<"��"<<tag3_1_2<<"�γ�����c=1/2���Ͻ�,"
	<<"��"<<tag3_5_8<<"�γ�����c=5/8���Ͻ�."
	<<"��"<<tag3_3_4<<"�γ�����c=3/4���Ͻ�."
	<<"��"<<tag3_7_8<<"�γ�����c=7/8���Ͻ�."
	<<"��"<<tag3_1<<"�γ�����c=1���Ͻ�."
	<<"��"<<tag3_3_2<<"�γ�����c=3/2���Ͻ�."
	<<"��"<<tag3_2<<"�γ�����c=2���Ͻ�."
	<<"��"<<tagq_3<<"�γ�����q/2."
	<<"b2������"<<endl<<"��"<<tag2_1_4<<"�γ�����c=1/4���Ͻ�,"
	<<"��"<<tag2_3_8<<"�γ�����c=3/8���Ͻ�,"
	<<"��"<<tag2_1_2<<"�γ�����c=1/2���Ͻ�,"
	<<"��"<<tag2_5_8<<"�γ�����c=5/8���Ͻ�."
	<<"��"<<tag2_3_4<<"�γ�����c=3/4���Ͻ�."
	<<"��"<<tag2_7_8<<"�γ�����c=7/8���Ͻ�."
	<<"��"<<tag2_1<<"�γ�����c=1���Ͻ�."
	<<"��"<<tagq_2<<"�γ�����q/2."
	<<endl;
	outf.close();
	len_Z++;
	}
}
