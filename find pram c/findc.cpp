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
/*读入参数 ki,ti,r,生成 u,a计算两次拆分的 xi,yi,得到多项式 Ax+By+Cz+D 的系数 ABCD
系数乘D逆得到 A'等，继而计算得到格基元素值 
*/

//这里只修改Z的取值，使之递增，然后再看第三个向量的长度 
int main(void)
{
	ZZ k1,k2,k3,k4,k5,k6,p,q,t1,t2,r,a,u,y1,y2,seed,x1,x2,gk1,gk2,gk3,gk4,gk5,gk6,det2,det,h,z;
	ZZ k1t1,k5a,temp,temp1,temp2,temp3,temp4,temp5,k2t2,k6ra;
	ZZ A,B,C,D,X,Y,Z,A_,B_,C_,D_;
	ZZ upperbound_test1_4,upperbound_test3_8,upperbound_test1_2,upperbound_test5_8,upperbound_test3_4,upperbound_test7_8,upperbound_test1,norml3_2,norml3_4,norml3_8,norml2_2,norml2_4,norml2_8;
	ZZ boundq2;
	int tag3_1_4,tag3_3_8,tag3_1_2,tag3_5_8,tag3_3_4,tag3_7_8,tagq_3,tag2_1_4,tag2_3_8,tag2_1_2,tag2_5_8,tag2_3_4,tag2_7_8,tagq_2,tag3_1,tag2_1,tag3_3_2,tag2_3_2,tag3_2;//超出上界计数 
	int i,count;
	int len_Z,len_q;
	mat_ZZ L;
	
	L.SetDims(4,4);
	
	ifstream inf;
	ofstream outf;
	
	inf.open("p and q.txt");//p,q,,盲化对取值固定，让r,ti,u,a,的取值变化 
	inf>>p>>q;
	inf.close();
	inf.open("pairs.txt");
	inf>>k1>>gk1>>k2>>gk2>>k3>>gk3>>k4>>gk4>>k5>>gk5>>k6>>gk6;
	inf.close();
//读入r,ti
	inf.open("ti and r.txt");
	inf>>t1>>t2>>r;
	inf.close();
//读入u,a 
	inf.open("u and a.txt");
	inf>>u>>a;
	inf.close();
//读入xi yi
	inf.open("xi and yi.txt");
	inf>>x1>>x2>>y1>>y2;
	inf.close();

 
	boundq2=sqr(q/2);

	seed=time(0);
	SetSeed(seed);

	len_q=NumBits(q);
	
	/*这部分计算xi yi*/	
	X=ZZ(11);
	
	cout<<"p="<<p<<endl<<"q="<<q<<endl<<"ti and r="<<t1<<endl<<t2<<endl<<r<<endl<<"u="<<u<<endl<<"a="<<a<<endl<<"xi="<<x1<<endl<<x2<<"yi="<<y1<<y2<<endl;
	cout<<"boundq2="<<boundq2<<endl;
	
	/*================循环开始=============*/
	
//	len_Z=NumBits(t2)+1;
	len_Z=64;
	while(len_Z<len_q)
	{			
	i=0;
	count=10000;
cout<<"len_Z="<<len_Z<<endl;
//getchar();
cout<<"finding......"<<endl; 

	tag3_1_4=0;//计数初始化 b3
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
	
	 /*设置Z的值*/ 
	//do{
		RandomLen(Z,len_Z);
		Y=operator*(X,Z);
//	}while(operator>(t2,Z)||operator>(operator*(r,t1),Y));
	/*over*/ 
	

	/*-----计算多项式系数ABCD,新多项式系数A'B'C'D'----- */
	
	A=x1;
	B=y1;
	C=-y2;
	D=-x2;
	
	

	/*系数计算完成*/
	
	
	B_=operator*(InvMod(A,q),B);//乘A的逆 ww
	C_=operator*(InvMod(A,q),C);
	D_=operator*(InvMod(A,q),D);
		
	//cout<<"B_="<<endl<<B_<<endl<<C_<<endl<<D_<<endl;
	
	/*构造格基矩阵L*/ 

	L(1,1)=X;     L(1,2)=operator*(B_,Y); L(1,3)=operator*(C_,Z); L(1,4)=D_;
	L(2,1)=ZZ(0); L(2,2)=operator*(q,Y);  L(2,3)=ZZ(0);           L(2,4)=ZZ(0);
	L(3,1)=ZZ(0); L(3,2)=ZZ(0);           L(3,3)=operator*(q,Z);  L(3,4)=ZZ(0);
	L(4,1)=ZZ(0); L(4,2)=ZZ(0);           L(4,3)=ZZ(0);           L(4,4)=q;	
	
	//cout<<L<<endl;
	/*约化L*/
	LLL(det2,L,0);
	
	/*计算b3范数的平方与猜测的上界*/
	
	det=operator*(operator*(power(q,3),X),operator*(Y,Z));//行列式 
	norml2_2=operator+(operator+(sqr(L(2,1)),sqr(L(2,2))),operator+(sqr(L(2,3)),sqr(L(2,4))));//b2范数平方 
	norml3_2=operator+(operator+(sqr(L(3,1)),sqr(L(3,2))),operator+(sqr(L(3,3)),sqr(L(3,4))));//b3范数平方 
	 norml3_4=sqr(norml3_2);
	 norml3_8=sqr(norml3_4);
	 norml2_4=sqr(norml2_2);
	 norml2_8=sqr(norml2_4);
	 	
	 upperbound_test1_4=operator*(2,det);//猜的上界c=1/4 	 
	 upperbound_test3_8=operator*(8,sqr(det));//猜的上界c=3/8 
	 upperbound_test1_2=operator*(4,det);//猜的上界c=1/2
	 upperbound_test5_8=operator*(32,sqr(det));//猜的上界c=5/8
	 upperbound_test3_4=operator*(8,det);//猜的上界c=3/4 
	 upperbound_test7_8=operator*(128,sqr(det));//猜的上界c=7/8 
	 upperbound_test1=operator*(16,det);//猜的上界c=1
	 
	
	 if(operator>(norml3_2,boundq2)) {
	 tagq_3++;//cout<<"norml3_2="<<norml3_2<<endl<<"boundq2="<<boundq2<<endl;
}
	 if(operator>(norml2_2,boundq2)) tagq_2++;

	 
	 /*四次方的比较*/
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
	outf<<"Z取"<<len_Z<<"bit，进行"<<count<<"次约化（u,a,r,ti,Z，均随机取值，重复格基可能性较小？）结果中，b3范数："<<endl<<"有"<<tag3_1_4<<"次超出了c=1/4的上界,"
	<<"有"<<tag3_3_8<<"次超出了c=3/8的上界,"
	<<"有"<<tag3_1_2<<"次超出了c=1/2的上界,"
	<<"有"<<tag3_5_8<<"次超出了c=5/8的上界."
	<<"有"<<tag3_3_4<<"次超出了c=3/4的上界."
	<<"有"<<tag3_7_8<<"次超出了c=7/8的上界."
	<<"有"<<tag3_1<<"次超出了c=1的上界."
	<<"有"<<tag3_3_2<<"次超出了c=3/2的上界."
	<<"有"<<tag3_2<<"次超出了c=2的上界."
	<<"有"<<tagq_3<<"次超出了q/2."
	<<"b2范数："<<endl<<"有"<<tag2_1_4<<"次超出了c=1/4的上界,"
	<<"有"<<tag2_3_8<<"次超出了c=3/8的上界,"
	<<"有"<<tag2_1_2<<"次超出了c=1/2的上界,"
	<<"有"<<tag2_5_8<<"次超出了c=5/8的上界."
	<<"有"<<tag2_3_4<<"次超出了c=3/4的上界."
	<<"有"<<tag2_7_8<<"次超出了c=7/8的上界."
	<<"有"<<tag2_1<<"次超出了c=1的上界."
	<<"有"<<tagq_2<<"次超出了q/2."
	<<endl;
	outf.close();
	len_Z++;
	}
}
